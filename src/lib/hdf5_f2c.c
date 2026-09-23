#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils_sub_names.h"

#include "hdf5.h"

// for ZFP compression
#ifdef ENABLE_ZFP_COMPRESSION
#include "H5Zzfp_plugin.h"
#endif

// Size of temporary hsize_t arrays
#define HDF5_MAX_DIMS 10

// Prototypes for internal routines
hid_t  fh5_set_hdf5_dtype (int type_code);
 void fh5_convert_array (int n,int *i_array,hsize_t *h_array,char *dir);

//********************** Deferred (multi-dataset) writes **********************
//
// Why this exists: in a parallel run every H5Dwrite of a compressed (filtered)
// dataset is a separate collective operation. The ranks exchange chunk lists,
// ship pieces of shared chunks to one owner rank, compress, agree on the new
// chunk sizes, allocate file space together, and then write. That coordination
// has a roughly fixed cost per dataset, which dominated RAMS analysis writes
// (about 0.4 s per variable at 256 ranks, while compressing the data itself
// took about 1 s for the whole file).
//
// How it works: between fh5d_begin_deferred_writes and
// fh5d_flush_deferred_writes, fh5d_write still creates each dataset and its
// selections immediately, but instead of writing it copies the caller's buffer
// and queues (dataset, memory type, memory space, file space, buffer copy). The
// flush then writes the whole queue with ONE collective H5Dwrite_multi call,
// so HDF5 does the chunk redistribution, compression bookkeeping, and file
// space allocation for all datasets together. H5Dwrite_multi needs HDF5 1.14.
//
// Why the buffer is copied: the Fortran caller (anal_write) frees its
// temporary array right after each shdf5_orec call, so the data must be kept
// alive here until the flush. This raises peak memory by one copy of every
// queued variable's subdomain. RAMS_HDF5_MULTI_BATCH (environment variable,
// default unlimited) caps the number of queued datasets; when the queue reaches
// the cap it is flushed early, trading some speed for lower peak memory.
//
// Invariants: every rank must queue the same datasets in the same order,
// because both the dataset creation and the flush are collective. anal_write
// satisfies this since all ranks walk the same variable table. Only parallel
// files (iphdf5 == 1) and numeric types are queued; anything else is written
// immediately exactly as before. With an older HDF5, begin/flush are no-ops
// and every write stays immediate.

#define DEFERRED_WRITE_MAX_DATASETS 2048

static int    deferred_writes_active = 0;       // 1 between begin and flush
static size_t deferred_batch_limit = DEFERRED_WRITE_MAX_DATASETS;
static size_t deferred_queue_length = 0;
static hid_t  deferred_dataset_ids[DEFERRED_WRITE_MAX_DATASETS];
static hid_t  deferred_memory_types[DEFERRED_WRITE_MAX_DATASETS];
static hid_t  deferred_memory_spaces[DEFERRED_WRITE_MAX_DATASETS];
static hid_t  deferred_file_spaces[DEFERRED_WRITE_MAX_DATASETS];
static void  *deferred_buffer_copies[DEFERRED_WRITE_MAX_DATASETS];

#if defined (RAMS_MPI) && H5_VERSION_GE(1,14,0)
#define RAMS_HAVE_MULTI_DATASET_WRITE 1
#else
#define RAMS_HAVE_MULTI_DATASET_WRITE 0
#endif

#if RAMS_HAVE_MULTI_DATASET_WRITE
//*****************************************************************************
// Write everything queued so far with one collective H5Dwrite_multi, then
// release the queued datasets, dataspaces, and buffer copies. Returns the
// H5Dwrite_multi status (negative on failure).
static herr_t fh5_write_deferred_queue (void)
  {
  herr_t herr = 0;
  hid_t transfer_plist;
  size_t i;

  if (deferred_queue_length == 0) return 0;

  // Collective transfer, same as the immediate path in fh5d_write.
  transfer_plist = H5Pcreate (H5P_DATASET_XFER);
  H5Pset_dxpl_mpio (transfer_plist,H5FD_MPIO_COLLECTIVE);

  herr = H5Dwrite_multi (deferred_queue_length, deferred_dataset_ids
                        ,deferred_memory_types, deferred_memory_spaces
                        ,deferred_file_spaces, transfer_plist
                        ,(const void **)deferred_buffer_copies);
  if (herr < 0)
    printf("Error in H5Dwrite_multi over %zu datasets, herr: %d\n"
          ,deferred_queue_length, herr);

  // Release everything the queue held, whether or not the write succeeded.
  H5Pclose (transfer_plist);
  for (i = 0; i < deferred_queue_length; i++)
    {
    H5Sclose (deferred_memory_spaces[i]);
    H5Sclose (deferred_file_spaces[i]);
    H5Dclose (deferred_dataset_ids[i]);
    free (deferred_buffer_copies[i]);
    }
  deferred_queue_length = 0;

  return herr;
  }

//*****************************************************************************
// Add one already-created dataset and its selections to the queue, keeping a
// private copy of the caller's memory buffer. Flushes first if the queue is
// at its batch limit. Returns 0 on success, negative on failure.
static int fh5_queue_deferred_write (hid_t dsetid,hid_t memtype,hid_t mspcid
                                    ,hid_t fspcid,void *buf)
  {
  // The memory dataspace describes the whole caller buffer (the selection
  // is only a part of it), so copy all of it to keep the offsets valid.
  hssize_t n_memory_points = H5Sget_simple_extent_npoints (mspcid);
  size_t bytes_per_point = H5Tget_size (memtype);
  size_t buffer_bytes;
  void *buffer_copy;

  if (n_memory_points < 0 || bytes_per_point == 0) return -1;
  buffer_bytes = (size_t)n_memory_points * bytes_per_point;

  if (deferred_queue_length >= deferred_batch_limit)
    {
    if (fh5_write_deferred_queue () < 0) return -1;
    }

  buffer_copy = malloc (buffer_bytes > 0 ? buffer_bytes : 1);
  if (buffer_copy == NULL)
    {
    printf("Error: could not allocate %zu bytes for a deferred HDF5 write\n"
          ,buffer_bytes);
    return -1;
    }
  memcpy (buffer_copy, buf, buffer_bytes);

  deferred_dataset_ids[deferred_queue_length]   = dsetid;
  deferred_memory_types[deferred_queue_length]  = memtype;
  deferred_memory_spaces[deferred_queue_length] = mspcid;
  deferred_file_spaces[deferred_queue_length]   = fspcid;
  deferred_buffer_copies[deferred_queue_length] = buffer_copy;
  deferred_queue_length++;

  return 0;
  }
#endif

//*****************************************************************************
// Start queuing writes for the currently open parallel file. Called from
// Fortran (shdf5_begin_deferred_writes) right after the file is created.
// Does nothing for serial files or when HDF5 is older than 1.14.
 void fh5d_begin_deferred_writes (int *iphdf5,int *hdferr)
  {
  static int reported = 0;
  *hdferr = 0;
#if RAMS_HAVE_MULTI_DATASET_WRITE
  if (*iphdf5 == 1)
    {
    // Optional cap on datasets per H5Dwrite_multi, to bound peak memory.
    const char *batch_setting = getenv ("RAMS_HDF5_MULTI_BATCH");
    long requested_limit = batch_setting ? atol (batch_setting) : 0;
    deferred_batch_limit = DEFERRED_WRITE_MAX_DATASETS;
    if (requested_limit > 0 && requested_limit < DEFERRED_WRITE_MAX_DATASETS)
      deferred_batch_limit = (size_t)requested_limit;
    deferred_queue_length = 0;
    deferred_writes_active = 1;
    }
#endif
  // Say once, from rank 0, which write path this build is using.
  if (!reported)
    {
    int my_rank = 0;
#if defined (RAMS_MPI)
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
#endif
    if (my_rank == 0)
      {
      if (deferred_writes_active)
        printf(" HDF5 analysis writes: multi-dataset (H5Dwrite_multi), batch limit %zu\n"
              ,deferred_batch_limit);
      else
        printf(" HDF5 analysis writes: one H5Dwrite per dataset\n");
      }
    reported = 1;
    }
  return;
  }

//*****************************************************************************
// Write everything queued since fh5d_begin_deferred_writes and stop queuing.
// Called from Fortran (shdf5_flush_deferred_writes) before the file is closed.
// Collective: all ranks must call it.
 void fh5d_flush_deferred_writes (int *hdferr)
  {
  *hdferr = 0;
#if RAMS_HAVE_MULTI_DATASET_WRITE
  if (deferred_writes_active)
    {
    if (fh5_write_deferred_queue () < 0) *hdferr = -1;
    deferred_writes_active = 0;
    }
#endif
  return;
  }

//************************** File routines ***********************************

//****************************************************************************
 void fh5f_open (char *locfn,int *iaccess,int *iphdf5,int64_t *fileid,int *hdferr)
  {
  unsigned flags;
  hid_t access_id;
  hid_t fileid_h5 = (hid_t)(*fileid);
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    // open a parallel file collectively
    access_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (access_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    // Do metadata collectively rather than letting every rank fetch and write
    // the same metadata blocks itself. Without this, all ranks contend on the
    // same small reads, which gets worse as rank count rises. Affects only how
    // metadata moves; the data written is unchanged.
    //
    // MPI-IO hints were tried alongside this and are deliberately NOT here: a
    // bisect at 256 ranks showed the hints, not these two calls, caused a bus
    // error partway through a run. Do not reintroduce them without testing at
    // full rank count.
    H5Pset_all_coll_metadata_ops (access_id,1);
    H5Pset_coll_metadata_write (access_id,1);
#else
    access_id = H5P_DEFAULT;
#endif
    }
  else
    {
    access_id = H5P_DEFAULT;
    }
  
  if(*iaccess == 1) flags = H5F_ACC_RDONLY;
  if(*iaccess == 2) flags = H5F_ACC_RDWR;
  
  fileid_h5 = H5Fopen (locfn,flags,access_id);
  *fileid = (int64_t)fileid_h5;
  //printf("fh5f_open_ - fileid: %d\n",*fileid);
  if (*fileid < 0) { *hdferr = *fileid; return;}
  
  // Parallel mode created a property which needs to be closed
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    *hdferr = H5Pclose (access_id);
#else
    *hdferr = 0;
#endif
    }
  else
    {
    // no errors if we made it here
    *hdferr = 0;
    }
  
  return;
  }

//*****************************************************************************
 void fh5f_create (char *locfn,int *iaccess,int *iphdf5,int64_t *fileid,int *hdferr)
  {
  unsigned flags;
  hid_t access_id,create_id;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    // open a parallel file collectively
    access_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (access_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    // Do metadata collectively rather than letting every rank fetch and write
    // the same metadata blocks itself. Without this, all ranks contend on the
    // same small reads, which gets worse as rank count rises. Affects only how
    // metadata moves; the data written is unchanged.
    //
    // MPI-IO hints were tried alongside this and are deliberately NOT here: a
    // bisect at 256 ranks showed the hints, not these two calls, caused a bus
    // error partway through a run. Do not reintroduce them without testing at
    // full rank count.
    H5Pset_all_coll_metadata_ops (access_id,1);
    H5Pset_coll_metadata_write (access_id,1);
#else
    access_id = H5P_DEFAULT;
#endif
    }
  else
    {
    access_id = H5P_DEFAULT;
    }
  create_id = H5P_DEFAULT;
  
  if(*iaccess == 1) flags = H5F_ACC_TRUNC;
  if(*iaccess == 2) flags = H5F_ACC_EXCL ;
  
  fileid_h5 = H5Fcreate (locfn,flags,create_id,access_id);
  *fileid = (int64_t)fileid_h5;
  //printf("fh5f_create: fileid: %d\n",*fileid);
  if (*fileid < 0) { *hdferr = *fileid; return;}
  
  // Parallel mode created a property which needs to be closed
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    *hdferr = H5Pclose (access_id);
#else
    *hdferr = 0;
#endif
    }
  else
    {
    // no errors if we made it here
    *hdferr = 0;
    }
  
  return;
  }

//****************************************************************************
 void fh5f_close (int64_t *fileid,int *hdferr)
  {
  herr_t herr;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);
  
  // Don't use H5close() to close everything since there might be more than one
  // file open (due to REVU for example).
  
  herr = H5Fclose (fileid_h5);
  
  *hdferr = herr;
  
  return;
  }

//******************************* Dataset routines ****************************

//*****************************************************************************
 void fh5d_info (int64_t *fileid,char *dname,int *ndims,int *dims,int *hdferr)
  {
  hid_t dsetid, dspcid;
  herr_t herr;
  hsize_t dimsc[HDF5_MAX_DIMS], maxdimsc[HDF5_MAX_DIMS];
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);

  int i;
  
  // Note that the code is the same whether or not accessing a
  // parallel HDF5 file.
  
  // Open the dataset, HDF5 1.8 API
  dsetid = H5Dopen (fileid_h5,dname,H5P_DEFAULT);
  //printf("fh5d_info: %d\n",dsetid);
  if (dsetid < 0) { *hdferr = dsetid; return;}
  
  // Get the description of the dataset
  dspcid = H5Dget_space (dsetid);
  //printf("fh5d_info: %d\n",dspcid);
  if (dspcid < 0) { *hdferr = dspcid; return;}
  
  // Get the number of dimensions and the sizes of each dimension
  // Place array values in temp space and then copy to output arguments
  *ndims = H5Sget_simple_extent_dims (dspcid, dimsc, maxdimsc);
  //printf("fh5d_info: %d %d %d %d\n",*ndims,dimsc[0],dimsc[1],dimsc[2]);
  if (*ndims < 0) { *hdferr = *ndims; return;}
  
  // Copy hsize_t arrays into FORTRAN arrays (passed in as type int)
  // since HDF5 routines expect hsize_t arrays and the lengths of
  // int and hsize_t may differ.
  fh5_convert_array (*ndims, dims, dimsc, "h2i");
  //printf("fh5d_info: %d %d %d %d\n",*ndims,dimsc[0],dimsc[1],dimsc[2]);
  
  // Close the dataset
  herr = H5Sclose (dspcid);
  herr = H5Dclose (dsetid);
  //printf("fh5d_info: %d\n",herr);
  
  *hdferr = herr;
  
  return;
  }

//****************************** Reading routines *****************************

//*****************************************************************************
// This routine will read a dataset from an HDF5 file. It accommodates hyperslab
// selection on both the input file space and output memory space. The selection
// of memory is contained in the input arguments m_ndims and m_select; ditto for
// the file selection and f_ndims, f_select.
//
// Since passing data in from FORTRAN, the hyperslab selection data has been
// joined together in the linear arrays m_select and f_select. This prevents the
// issue with column-major versus row-major. The data are stored in these
// arrays in groups that are the same size and the corresponding ndims value
// (m_ndims and f_ndims), and go in the order: dims, block, count, offset, stride.
//
// Another complication of the caller being FORTRAN is that the HDF5 routines
// want hsize_t arrays (instead of int arrays). FORTRAN has no idea what hsize_t
// is, so integer arrays are used instead. The lengths of int and hsize_t are
// likely to be different so it is necessary to convert between the argument int
// arrays and internal hsize_t arrays.
//
 void fh5d_read (int64_t *fileid,char *dname,int *h5type,int *m_ndims,int *m_select
               ,int *f_ndims,int *f_select,void *buf,int *hdferr)
  {
  hid_t dsetid, mspcid, fspcid;
  herr_t herr;
  hid_t memtype;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);

  int i;
  
  hsize_t m_dims[HDF5_MAX_DIMS];   // memory selection description
  hsize_t m_block[HDF5_MAX_DIMS];
  hsize_t m_count[HDF5_MAX_DIMS];
  hsize_t m_offset[HDF5_MAX_DIMS];
  hsize_t m_stride[HDF5_MAX_DIMS];

  hsize_t f_dims[HDF5_MAX_DIMS];   // file selection description
  hsize_t f_block[HDF5_MAX_DIMS];
  hsize_t f_count[HDF5_MAX_DIMS];
  hsize_t f_offset[HDF5_MAX_DIMS];
  hsize_t f_stride[HDF5_MAX_DIMS];
  
  // Note that the code is the same whether or not accessing a
  // parallel HDF5 file.
  
  //******************* CONVERT INPUT ARGUMENTS *************************

  // Extract and convert the input int arrays into hsize_t arrays
  // Note: &(array[i]) gives the address of the ith element of array[] making
  // it appear to fh5_convert_array() that this is the zeroth element of
  // its argument array.
  fh5_convert_array (*m_ndims, &(m_select[0]),            m_dims,   "i2h");
  fh5_convert_array (*m_ndims, &(m_select[*m_ndims]),     m_block,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*2]), m_count,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*3]), m_offset, "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*4]), m_stride, "i2h");

  fh5_convert_array (*f_ndims, &(f_select[0]),            f_dims,   "i2h");
  fh5_convert_array (*f_ndims, &(f_select[*m_ndims]),     f_block,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*2]), f_count,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*3]), f_offset, "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*4]), f_stride, "i2h");

  // Set HDF5 data type
  memtype = fh5_set_hdf5_dtype(*h5type);

  //******************* SET UP SELECTION ******************************
  
  // Create descriptions of the file and memory data spaces

  dsetid = H5Dopen (fileid_h5,dname,H5P_DEFAULT); // file data space
  if (dsetid < 0) { *hdferr = dsetid; return;}
  fspcid = H5Dget_space (dsetid);
  if (fspcid < 0) { *hdferr = fspcid; return;}
  
  mspcid = H5Screate_simple (*m_ndims,m_dims,NULL); // memory data space
  if (mspcid < 0) { *hdferr = mspcid; return;}
  
  // Make the hyperslab selections

  // file
  herr = H5Sselect_hyperslab (fspcid,H5S_SELECT_SET,f_offset,f_stride
                             ,f_count,f_block);
  if (herr < 0) { *hdferr = herr; return;}

  // memory
  herr = H5Sselect_hyperslab (mspcid,H5S_SELECT_SET,m_offset,m_stride
                             ,m_count,m_block);
  if (herr < 0) { *hdferr = herr; return;}
  
  //******************* DO THE READ ******************************

  herr = H5Dread (dsetid,memtype,mspcid,fspcid,H5P_DEFAULT,buf);
  if (herr < 0) { *hdferr = herr; return;}
  
  // Close

  herr = H5Sclose (mspcid);
  herr = H5Sclose (fspcid);
  herr = H5Dclose (dsetid);
  
  *hdferr = herr;
  
  return;
  }

//**************************** Writing routines *******************************

//*****************************************************************************
// This routine will write a dataset into an HDF5 file. It accommodates hyperslab
// selection on both the input file space and output memory space. The selection
// of memory is contained in the input arguments m_ndims and m_select; ditto for
// the file selection and f_ndims, f_select.
//
// See the comments for fh5d_read() for info about the purpose of the hsize_t arrays.
//
 void fh5d_write (int64_t *fileid,char *dname,int *h5type,int *iphdf5,int *m_ndims
                ,int *m_select,int *f_ndims,int *f_select,int *f_csize, int *enable_truncation,
                float *truncation_accuracy, void *buf, int *hdferr)
  {
  hid_t dsetid, fspcid, mspcid, propid;
  int i;
  herr_t herr;
  hid_t memtype;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);


  hsize_t m_dims[HDF5_MAX_DIMS];   // memory selection description
  hsize_t m_block[HDF5_MAX_DIMS];
  hsize_t m_count[HDF5_MAX_DIMS];
  hsize_t m_offset[HDF5_MAX_DIMS];
  hsize_t m_stride[HDF5_MAX_DIMS];

  hsize_t f_dims[HDF5_MAX_DIMS];   // file selection description
  hsize_t f_block[HDF5_MAX_DIMS];
  hsize_t f_count[HDF5_MAX_DIMS];
  hsize_t f_offset[HDF5_MAX_DIMS];
  hsize_t f_stride[HDF5_MAX_DIMS];

  hsize_t chunk_size[HDF5_MAX_DIMS]; // for the file data space
  //bool enable_truncation = true;
  // ZFP Parameter array
  int cd_nelmts = 10;
  unsigned int cd_vals[cd_nelmts];
  //printf("ZFP Options: enable? %i accuracy? %f \n", *enable_truncation, *truncation_accuracy);
  //float accuracy = 0.1;
  
  //******************* CONVERT INPUT ARGUMENTS *************************

  // Extract and convert the input int arrays into hsize_t arrays
  // Note: &(array[i]) gives the address of the ith element of array[] making
  // it appear to fh5_convert_array() that this is the zeroth element of
  // its argument array.
  fh5_convert_array (*m_ndims, &(m_select[0]),            m_dims,   "i2h");
  fh5_convert_array (*m_ndims, &(m_select[*m_ndims]),     m_block,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*2]), m_count,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*3]), m_offset, "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*4]), m_stride, "i2h");

  fh5_convert_array (*f_ndims, &(f_select[0]),            f_dims,   "i2h");
  fh5_convert_array (*f_ndims, &(f_select[*f_ndims]),     f_block,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*2]), f_count,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*3]), f_offset, "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*4]), f_stride, "i2h");

  fh5_convert_array (*f_ndims,f_csize,chunk_size,"i2h");

  // Set HDF5 data type
  memtype = fh5_set_hdf5_dtype (*h5type);

  //******************* CREATE THE FILE DATASET ***********************

  // This must be done before the selection section

  // Create description of data space for the file dataset
  fspcid = H5Screate_simple (*f_ndims,f_dims,f_dims);
  if (fspcid < 0) { 
    *hdferr = fspcid; 
    printf("Error in H5Screate_simple, fspcid: %d", fspcid);
    return;}
  
  // Create properties for chunking, gzip compression. Can't do
  // compression for PHDF5 file.
  propid = H5Pcreate (H5P_DATASET_CREATE);
  //printf("fh5d_write - propid: %d\n",propid);
  if (propid < 0) { 
    *hdferr = propid; 
    printf("Error in H5Pcreate, propid: %d", propid);
    return;}

  // Enable compression options. 
  if(*enable_truncation == 1){
    #ifdef ENABLE_ZFP_COMPRESSION
    herr = H5Pset_chunk (propid,*f_ndims,chunk_size);
      if (herr < 0) { 
    *hdferr = herr; 
    printf("Error in chunking, herr: %d", herr);
    return;}

    //printf("Truncating variable %s to accuracy %f\n", dname, *truncation_accuracy);
    H5Pset_zfp_accuracy_cdata(*truncation_accuracy, cd_nelmts, cd_vals);
    herr = H5Pset_filter(propid, H5Z_FILTER_ZFP, H5Z_FLAG_MANDATORY, cd_nelmts, cd_vals);

      if (herr < 0) { 
    *hdferr = herr; 
    printf("Error in ZFP, herr: %d", herr);
    return;}

    #else
    // we want truncation on, but the ZFP library isn't included.
    *hdferr=-9999;
    printf("Error: ZFP compression requested, but ZFP library not included.\n");
    return;
    #endif
  }
  else{
  herr = H5Pset_chunk (propid,*f_ndims,chunk_size);
  #if H5_VERSION_GE(1,10,3)
  #ifdef ENABLE_PARALLEL_COMPRESSION 
    
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
  
  #else
  if (*iphdf5 != 1)
    {
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
    }
  #endif//enable parallel compression ifdef
  #else//H5_VERSION_GE
  if (*iphdf5 != 1)
    {
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
    }
    #endif //H5_VERSION_GE if
  }
  if (herr < 0) { 
    *hdferr = herr; 
    printf("Error in chunking, herr: %d", herr);
    return;}
  
  // Create the dataset
  dsetid = H5Dcreate (fileid_h5,dname,memtype,fspcid,H5P_DEFAULT,propid,H5P_DEFAULT);
  if (dsetid < 0) { 
    *hdferr = dsetid; 
    printf("Error in H5Dcreate, dsetid: %lld \n",dsetid);
    printf("Fileid_h5: %lld, dname: %s \n",fileid_h5,dname);
    printf("memtype: %d, fspcid: %d, propid: %d \n",memtype,fspcid, propid);
    for(int i=0; i<*f_ndims; i++){
      printf("Chunk size [%d], %lld \n",i,chunk_size[i]);
    }
    return;}

    

  // Close the items used for this section, but keep dsetid open for
  // subsequent use
  H5Pclose (propid);
  H5Sclose (fspcid);
  
  //******************* SET UP SELECTION ******************************
  
  // Create descriptions of the file and memory data spaces

  fspcid = H5Dget_space (dsetid);
  if (fspcid < 0) { *hdferr = fspcid; return;}
  
  mspcid = H5Screate_simple (*m_ndims,m_dims,NULL); // memory data space
  if (mspcid < 0) { *hdferr = mspcid; return;}
  
  // Make the hyperslab selections

  // file
  herr = H5Sselect_hyperslab (fspcid,H5S_SELECT_SET,f_offset,f_stride
                             ,f_count,f_block);
  if (herr < 0) { *hdferr = herr; return;}

  // memory
  herr = H5Sselect_hyperslab (mspcid,H5S_SELECT_SET,m_offset,m_stride
                             ,m_count,m_block); 
  if (herr < 0) { *hdferr = herr; return;}
  
  //******************* DO THE WRITE ******************************

#if RAMS_HAVE_MULTI_DATASET_WRITE
  // Inside a deferred-write section, queue this dataset instead of writing
  // it now; the flush writes all queued datasets in one collective call. The
  // dataset and both dataspaces stay open until then. Strings (type 3) are
  // never queued, since H5T_STRING is not a concrete memory type.
  if (deferred_writes_active && *iphdf5 == 1 && *h5type != 3)
    {
    if (fh5_queue_deferred_write (dsetid,memtype,mspcid,fspcid,buf) < 0)
      { *hdferr = -1; return; }
    *hdferr = 0;
    return;
    }
#endif

  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    propid = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (propid,H5FD_MPIO_COLLECTIVE);
#else
    propid = H5P_DEFAULT;
#endif
    }
  else
    {
    propid = H5P_DEFAULT;
    }
  
  herr = H5Dwrite (dsetid,memtype,mspcid,fspcid,propid,buf);
  if (herr < 0) { *hdferr = herr; return;}
  
  // Close
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    H5Pclose (propid);
#endif
    }
  herr = H5Sclose (mspcid);
  herr = H5Sclose (fspcid);
  herr = H5Dclose (dsetid);
  
  *hdferr = herr;
  
  return;
  }

//*************************** Utilities **************************************

//****************************************************************************
// This routine will translate the type codes from FORTAN to the built in
// HDF5 type codes.

hid_t fh5_set_hdf5_dtype (int type_code)
  {
  hid_t dtype;

  if (type_code == 1) dtype = H5T_NATIVE_INT32;
  if (type_code == 2) dtype = HDF5_FLOAT_TYPE;
  if (type_code == 3) dtype = H5T_STRING;
  if (type_code == 4) dtype = H5T_NATIVE_DOUBLE;
  if (type_code == 5) dtype = H5T_NATIVE_HBOOL;

  return dtype;
  }

//****************************************************************************
// This routine will copy the contents between an integer array and an hsize_t
// array. It is likely to encounter a system where hsize_t and int are not the same
// length. As an example, this would necessitate the need to convert back and
// forth when the integer array is used as an argument from FORTRAN for example.
//
// dir == 'i2h' --> convert int to hsize_t
// dir == 'h2i' --> convert hsize_t to int

 void fh5_convert_array (int n,int *i_array,hsize_t *h_array,char *dir)
  {
  int i;

  if (strcmp(dir, "i2h") == 0)
    {
    // copy int array into the hsize_t array
    for (i = 0; i < n; i++)
      {
      h_array[i] = i_array[i];
      }
    }
  else if (strcmp(dir, "h2i") == 0)
    {
    // copy hsize_t array into the int array
    for (i = 0; i < n; i++)
      {
      i_array[i] = h_array[i];
      }
    }
  else
    {
    printf ("ERROR: fh5_convert_array:");
    printf ("Do not recognize dir argument value: %d\n", dir);
    }

  return;
  }
