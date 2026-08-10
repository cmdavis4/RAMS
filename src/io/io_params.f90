!##############################################################################
Module io_params

use grid_dims

implicit none

character(len=32) :: lite_vars(maxlite)
! this must be fixed at a 32 bit float because that is what the C library expects
! Default 0. = lossless, so a RAMSIN that omits ACC_LT_VAR behaves exactly as
! it did before ZFP support existed. Without this initializer the array is
! garbage for any such RAMSIN and shdf5_orec either aborts or compresses at a
! random tolerance.
real(kind=4) :: lite_var_acc(maxlite) = 0.
character(len=strl1) :: hfilin,afilepref

integer :: ipast_sfc
!-------------------------------------------------------------------------------
integer :: ioutput,iclobber,nlite_vars
integer :: iuvwtend = 0
real    :: frqstate(maxgrds),avgtim,frqlite,frqmean,frqboth,frqst_keep  
!-------------------------------------------------------------------------------
integer, dimension(maxgrds) :: itoptflg,isstflg,ivegtflg,isoilflg  &
                              ,ndviflg,itopsflg,iz0flg
real                        :: z0fact
real, dimension(maxgrds)    :: z0max,toptenh,toptwvl
!-------------------------------------------------------------------------------
character(len=strl1), dimension(maxgrds) :: itoptfn,isstfn,ivegtfn,isoilfn  &
                                        ,ndvifn
!-------------------------------------------------------------------------------
character(len=strl1) :: sstfpfx,sfcfiles,topfiles
character(len=strl1), dimension(maxsstfiles,maxgrds) :: fnames_sst
character(len=14), dimension(maxsstfiles,maxgrds) :: itotdate_sst
!-------------------------------------------------------------------------------
integer                             :: iupdsst,isstcyclic,isstcycdata
integer,dimension(maxgrds)          :: nsstfiles,isstflp,isstflf
real,dimension(maxgrds)             :: ssttime1,ssttime2
!-------------------------------------------------------------------------------
character(len=strl1) :: ndvifpfx
character(len=strl1), dimension(maxndvifiles,maxgrds) :: fnames_ndvi
character(len=14), dimension(maxndvifiles,maxgrds) :: itotdate_ndvi
!-------------------------------------------------------------------------------
integer                             :: iupdndvi,indvicyclic,indvicycdata
integer,dimension(maxgrds)          :: nndvifiles,indviflp,indviflf
real,dimension(maxgrds)             :: ndvitime1,ndvitime2

END MODULE io_params
