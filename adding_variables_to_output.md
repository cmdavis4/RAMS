# Adding 3D State Variables to RAMS Model Output

This document explains how to add an existing three-dimensional state variable to the RAMS model output files. This assumes the variable is already allocated and calculated internally within the model.

## Overview of the Variable Registration System

RAMS uses a dynamic variable table system that allows flexible control over which variables are written to output files. The system involves:

1. **Variable Tables** ([src/memory/var_tables.f90](src/memory/var_tables.f90)) - A registration system that tracks all outputtable variables
2. **Memory Modules** ([src/memory/mem_*.f90](src/memory/)) - Fortran derived types that group related variables
3. **I/O Routines** ([src/io/anal_write.f90](src/io/anal_write.f90)) - Code that writes registered variables to HDF5 files
4. **Namelist Control** - Runtime selection of variables via RAMSIN namelist file

## The Variable Table Structure

Each registered variable has the following attributes (defined in [var_tables.f90:10-17](src/memory/var_tables.f90#L10-L17)):

```fortran
type var_tables_r
   real, pointer :: var_p        ! Pointer to variable data
   real, pointer :: var_m        ! Pointer to time-averaged data (if enabled)
   integer :: npts               ! Total number of points
   integer :: idim_type          ! Dimension type (2=2D, 3=3D, etc.)
   integer :: ianal              ! Output to standard analysis files (1=yes)
   integer :: imean              ! Include in time averaging (1=yes)
   integer :: ilite              ! Available for LITE_VARS selection (1=yes)
   integer :: impti, impt1       ! MPI communication flags
   integer :: irecycle_sfc       ! Surface recycling flag
   character(len=32) :: name     ! Variable name
end type
```

## Step-by-Step Process to Add a Variable to Output

Assuming you have a 3D variable that already exists in memory and is being calculated during model integration, follow these steps:

### Step 1: Locate the Memory Module

Identify which memory module contains (or should contain) your variable. Variables are organized by functionality:

- **[src/memory/mem_basic.f90](src/memory/mem_basic.f90)** - Core state variables (u, v, w, theta, pressure, etc.)
- **[src/memory/mem_micro.f90](src/memory/mem_micro.f90)** - Microphysics variables (cloud/precipitation mixing ratios, number concentrations)
- **[src/memory/mem_turb.f90](src/memory/mem_turb.f90)** - Turbulence variables (TKE, diffusion coefficients)
- **[src/memory/mem_radiate.f90](src/memory/mem_radiate.f90)** - Radiation variables
- **[src/memory/mem_cuparm.f90](src/memory/mem_cuparm.f90)** - Cumulus parameterization variables

For this example, let's assume we want to output a variable called `my_3d_var` that's already in the `basic_vars` derived type in [mem_basic.f90](src/memory/mem_basic.f90).

### Step 2: Verify Variable is in the Derived Type

Check that your variable is declared in the appropriate derived type. In [mem_basic.f90:6-23](src/memory/mem_basic.f90#L6-L23):

```fortran
Type basic_vars
   ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
      up,uc,vp,vc,wp,wc,pp,pc, &
      rv,theta,thp,rtp, &
      pi0,th0,rvt0,dn0,dn0u,dn0v, &
      wp_buoy_theta,wp_buoy_cond,wp_advdif, &
      my_3d_var  ! <-- Your variable should be here
End Type
```

### Step 3: Verify Variable is Allocated

Check the allocation routine to ensure your variable is being allocated. In [mem_basic.f90:30-73](src/memory/mem_basic.f90#L30-L73), look for the `alloc_basic` subroutine:

```fortran
Subroutine alloc_basic (basic,n1,n2,n3)
   type (basic_vars) :: basic
   integer, intent(in) :: n1,n2,n3

   allocate (basic%up(n1,n2,n3))
   allocate (basic%vp(n1,n2,n3))
   ! ... other allocations ...
   allocate (basic%my_3d_var(n1,n2,n3))  ! <-- Your variable allocation
```

### Step 4: Register Variable in the Variable Table

This is the **critical step** for enabling output. Locate the `filltab_*` subroutine for your memory module (e.g., `filltab_basic` in [mem_basic.f90:114-239](src/memory/mem_basic.f90#L114-L239)).

Add a vtables2 call to register your variable:

```fortran
Subroutine filltab_basic (basic,basicm,imean,n1,n2,n3,ng)
   use var_tables
   type (basic_vars) :: basic,basicm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

   npts=n1*n2*n3  ! Total number of points for 3D variables

   ! Existing variable registrations...
   if (allocated(basic%up))  &
      CALL vtables2 (basic%up(1,1,1),basicm%up(1,1,1), &
                     ng, npts, imean, 'UP :3:anal:mpti')

   ! Add your variable registration here:
   if (allocated(basic%my_3d_var))  &
      CALL vtables2 (basic%my_3d_var(1,1,1),basicm%my_3d_var(1,1,1), &
                     ng, npts, imean, 'MY_3D_VAR :3:anal:mpti')
```

### Step 5: Understanding the vtables2 Registration String

The registration string format is: `'VARNAME :dimtype:flag1:flag2:...'`

**Components:**

1. **VARNAME** - Variable name as it appears in output files (max 32 characters)
   - Use descriptive, uppercase names
   - This name will appear in HDF5 files and is used in LITE_VARS

2. **:dimtype** - Dimension type code:
   - `:2` - 2D horizontal field (nxp, nyp)
   - `:3` - 3D field (nzp, nxp, nyp) ← Most common for atmospheric variables
   - `:4` - 4D soil field (nxp, nyp, nzg, npatch)
   - `:5` - 4D snow field (nxp, nyp, nzs, npatch)
   - `:6` - 4D microphysics bins (nxp, nyp, nzp, nbins)
   - `:7` - 4D aerosol categories (nxp, nyp, nzp, naer)
   - `:10` - 3D ocean field (nxp, nyp, nkppz)

3. **Output Control Flags:**
   - `:anal` - **Include in standard analysis output files** (most important!)
     - Without this flag, the variable is NOT written to output
     - With this flag, variable is always written to analysis files
   - `:lite` - Make available for selective output via LITE_VARS namelist
     - Allows users to choose this variable at runtime without recompiling

4. **MPI Communication Flags** (for parallel runs):
   - `:mpti` - Variable needs MPI communication (almost always needed)
   - `:mpt1` - Additional MPI communication flag for tendency variables

5. **Special Flags:**
   - `:recycle_sfc` - For surface variables that recycle during restarts

**Common Registration Patterns:**

```fortran
! Standard 3D state variable (always output)
'THETA :3:anal:mpti'

! 3D variable available for selective output
'BUOYANCY :3:lite:mpti'

! 3D tendency variable with extra MPI flag
'THP :3:anal:mpti:mpt1'

! 2D surface variable
'TOPO :2:anal:mpti'

! 3D variable that is both always output AND selectable
'THETA :3:anal:lite:mpti'
```

### Step 6: Rebuild the Model

After modifying the source code:

```bash
cd bin.rams
make clean
make
```

The compilation will pick up your changes and rebuild the executable.

### Step 7: Run and Verify Output

**For variables with `:anal` flag:**
Simply run RAMS normally. The variable will appear in all analysis output files (A-*.h5).

**For variables with `:lite` flag only:**
Add the variable name to the RAMSIN namelist:

```fortran
$MODEL_FILE_INFO
   ! ... other namelist parameters ...

   NLITE_VARS = 3,           ! Number of LITE variables to output
   LITE_VARS = 'MY_3D_VAR', 'THETA', 'RV',  ! List of variable names

   FRQLITE = 600.,           ! Output frequency for LITE files (seconds)
$END
```

This creates separate "LITE" output files (L-*.h5) containing only the selected variables, which can significantly reduce output size and I/O time.

### Step 8: Verify Variable is in Output

Check the HDF5 output files:

```bash
# List contents of an analysis file
h5dump -H A-2025-01-01-000000-g1.h5 | grep MY_3D_VAR

# Or use Python:
python3 << EOF
import h5py
f = h5py.File('A-2025-01-01-000000-g1.h5', 'r')
print('MY_3D_VAR' in f.keys())
print(f['MY_3D_VAR'].shape)
EOF
```

## Complete Example: Adding wp_buoy_theta

Here's a real example from the codebase showing how `wp_buoy_theta` (vertical velocity budget term) was added:

**In [mem_basic.f90:9-13](src/memory/mem_basic.f90#L9-L13):**
```fortran
Type basic_vars
   real, allocatable, dimension(:,:,:) :: &
      up,uc,vp,vc,wp,wc,pp,pc, &
      ! ... other variables ...
      wp_buoy_theta,wp_buoy_cond,wp_advdif  ! Budget terms
End Type
```

**In [mem_basic.f90:66-70](src/memory/mem_basic.f90#L66-L70) (conditional allocation):**
```fortran
if(imbudget>=1) then
   allocate (basic%wp_buoy_theta(n1,n2,n3))
   allocate (basic%wp_buoy_cond(n1,n2,n3))
   allocate (basic%wp_advdif(n1,n2,n3))
endif
```

**In [mem_basic.f90:214-225](src/memory/mem_basic.f90#L214-L225) (registration):**
```fortran
if (allocated(basic%wp_buoy_theta))  &
   CALL vtables2 (basic%wp_buoy_theta(1,1,1),basicm%wp_buoy_theta(1,1,1), &
                  ng, npts, imean, 'WP_BUOY_THETA :3:anal:mpti')
if (allocated(basic%wp_buoy_cond))  &
   CALL vtables2 (basic%wp_buoy_cond(1,1,1),basicm%wp_buoy_cond(1,1,1), &
                  ng, npts, imean, 'WP_BUOY_COND :3:anal:mpti')
if (allocated(basic%wp_advdif))  &
   CALL vtables2 (basic%wp_advdif(1,1,1),basicm%wp_advdif(1,1,1), &
                  ng, npts, imean, 'WP_ADVDIF :3:anal:mpti')
```

Note the use of `if (allocated(...))` to handle conditional compilation based on namelist options.

## Additional Considerations

### Time-Averaged Output

The variable table system supports time-averaged output. When `AVGTIM` is set in the RAMSIN namelist:

- The `basicm` counterpart stores accumulated values
- Variables with `imean=1` are time-averaged
- Use `FRQMEAN` or `FRQBOTH` to control averaging output frequency

### Memory Considerations

Each registered 3D variable requires memory for:
- Current value: `(nzp × nxp × nyp) × 8 bytes` (single precision)
- Time average (if enabled): Same size again

For a 100×100×40 grid, each 3D variable uses ~3.2 MB.

### Output File Size

Each 3D variable adds to HDF5 file size. Strategies to manage this:

1. **Use LITE_VARS** - Only output essential variables
2. **Reduce output frequency** - Increase FRQSTATE or FRQLITE
3. **Use compression** - Enable HDF5 compression (if compiled with `-DENABLE_PARALLEL_COMPRESSION`)
4. **Selective grids** - Use grid-dependent FRQSTATE settings

### Variable Naming Conventions

Follow RAMS conventions:
- **State variables**: UP, VP, WP, THETA, RV, etc. (present/prognostic)
- **Base state**: TH0, PI0, DN0 (reference state with '0')
- **Tendencies**: THP, RTP (with 'P' suffix)
- **Diagnostics**: Descriptive names (WP_BUOY_THETA, CLOUD_FRAC)
- **Fluxes**: Often use 'FLX' or end in 'F'

### Deallocation

Don't forget to add deallocation in the `dealloc_*` subroutine if you added allocation:

```fortran
Subroutine dealloc_basic (basic)
   type (basic_vars) :: basic

   if (allocated(basic%my_3d_var)) deallocate (basic%my_3d_var)
   ! ... other deallocations ...
```

## Summary Checklist

For adding an existing 3D variable to output:

- [ ] Variable is in the derived type definition (`Type basic_vars`)
- [ ] Variable is allocated in `alloc_*` subroutine
- [ ] Variable is deallocated in `dealloc_*` subroutine
- [ ] **Variable is registered with `vtables2` in `filltab_*` subroutine**
- [ ] Registration string includes `:3:anal:mpti` for standard 3D output
- [ ] Optionally include `:lite` flag for LITE_VARS capability
- [ ] Model rebuilt with `make clean && make`
- [ ] Output verified in HDF5 files

## Next Steps: Creating New Calculated Variables

Once you understand how to add existing variables to output, the next step is to:

1. Declare a new variable in the appropriate memory module
2. Allocate it during model initialization
3. Calculate its values during the model integration (in src/core/, src/micro/, etc.)
4. Register it for output using the steps above

The key is that the variable must be populated with data during the time integration loop before the output routines are called.
