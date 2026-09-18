!##############################################################################
Subroutine simdata (fname,fln)

use mem_grid
use node_mod, only:mmzp,mzp
use io_params, only:frqstate

implicit none

integer :: i
character(len=*) :: fname,fln
character(len=64) :: fullname

fullname=trim(fln)//trim(fname)

open(50,file=trim(fullname),status='unknown')
write(50,*)iyear1
write(50,*)imonth1
write(50,*)idate1
write(50,*)itime1
write(50,*)
write(50,*)dtlt
write(50,*)time
write(50,*)ngrid
write(50,*)mzp
write(50,*)
write(50,*)ngrids
do i=1,ngrids
write(50,*)mmzp(i)
enddo
do i=1,ngrids
write(50,*)frqstate(i)
enddo
close(50)

return
END SUBROUTINE simdata

!##############################################################################
Subroutine iofil (fname,x,ix,wf,fln)

use mem_grid, only:time,dtlt

implicit none

character(len=*) :: fname,fln
character(len=64) :: fullname
integer :: ix                    !number vertical levels
integer :: wf                    !read or write flag
real, dimension(ix) :: x         !variable to read/write

integer k,l

fullname=trim(fln)//trim(fname)

!READ SCM DATA
if(wf == 1) then 
  stop 'No reading of SCM data for full RAMS model'
!WRITE SCM DATA
elseif (wf >= 2) then
  open(50,file=trim(fullname),status='unknown')
  if(wf==2)write(50,*)'Elapsed time',time
  if(wf==3)write(50,*)'Elapsed time',time+dtlt
  do k=1,ix
    write(50,'(i5,e26.18)') k,x(k)
  enddo
  close(50)
endif

return
END SUBROUTINE iofil

!##############################################################################
Subroutine readwrite_scm (wf,fln)

use mem_basic
use mem_micro
use mem_grid
use mem_radiate
use mem_leaf
use micphys
use rconstants
use node_mod, only:mi0,mj0,my_rams_num,my_scm_num,mzp

implicit none

character(len=*) :: fln
integer :: m1,k,i,j,wf,ng

if(iprntstmt>=1 .and. my_rams_num==my_scm_num) &
 print*,'Writing out SCM variables at time:',time,' on grid:',ngrid

!Set the grid we are on to local variable
ng = ngrid
m1 = mzp

!Get sub-grid point associated with absolute point
i=iscmx-mi0(ng)
j=iscmy-mj0(ng)

!Compute pressure from Exner function
do k=1,m1
 pitot(k) = basic_g(ng)%pi0(k,i,j) + basic_g(ng)%pp(k,i,j)
 press(k) = p00 * (pitot(k) * cpi) ** cpor
enddo
CALL iofil ('pressure.txt'   ,press(:),m1,wf,fln)

!Output all the potential variables in the chosen column
CALL iofil ('soil_water.txt',leaf_g(ng)%soil_water(nzg,i,j,:),npatch,wf,fln)
CALL iofil ('patch_area.txt',leaf_g(ng)%patch_area(i,j,:),npatch,wf,fln)
CALL iofil ('leaf_class.txt',leaf_g(ng)%leaf_class(i,j,:),npatch,wf,fln)
CALL iofil ('soil_text.txt' ,leaf_g(ng)%soil_text(nzg,i,j,:),npatch,wf,fln)

CALL iofil ('zm.txt',zm,m1,wf,fln)
CALL iofil ('zt.txt',zt,m1,wf,fln)
CALL iofil ('glat.txt',grid_g(ng)%glat(i,j),1,wf,fln)
CALL iofil ('glon.txt',grid_g(ng)%glon(i,j),1,wf,fln)
CALL iofil ('topt.txt',grid_g(ng)%topt(i,j),1,wf,fln)

CALL iofil ('up.txt'   ,basic_g(ng)%up   (:,i,j),m1,wf,fln)
CALL iofil ('vp.txt'   ,basic_g(ng)%vp   (:,i,j),m1,wf,fln)
CALL iofil ('wp.txt'   ,basic_g(ng)%wp   (:,i,j),m1,wf,fln)
CALL iofil ('thp.txt'  ,basic_g(ng)%thp  (:,i,j),m1,wf,fln)
CALL iofil ('theta.txt',basic_g(ng)%theta(:,i,j),m1,wf,fln)
CALL iofil ('pp.txt'   ,basic_g(ng)%pp   (:,i,j),m1,wf,fln)
CALL iofil ('rtp.txt'  ,basic_g(ng)%rtp  (:,i,j),m1,wf,fln)
CALL iofil ('rv.txt'   ,basic_g(ng)%rv   (:,i,j),m1,wf,fln)
CALL iofil ('dn0.txt'  ,basic_g(ng)%dn0  (:,i,j),m1,wf,fln)
CALL iofil ('pi0.txt'  ,basic_g(ng)%pi0  (:,i,j),m1,wf,fln)

if(iswrtyp == 3 .or. ilwrtyp == 3) then
 CALL iofil ('fthrd.txt'  ,radiate_g(ng)%fthrd  (:,i,j),m1,wf,fln)
 CALL iofil ('fthrdp.txt' ,radiate_g(ng)%fthrdp (:,i,j),m1,wf,fln)
 CALL iofil ('bext.txt'   ,radiate_g(ng)%bext   (:,i,j),m1,wf,fln)
 CALL iofil ('swdn.txt'   ,radiate_g(ng)%swdn   (:,i,j),m1,wf,fln)
 CALL iofil ('swup.txt'   ,radiate_g(ng)%swup   (:,i,j),m1,wf,fln)
 CALL iofil ('lwdn.txt'   ,radiate_g(ng)%lwdn   (:,i,j),m1,wf,fln)
 CALL iofil ('lwup.txt'   ,radiate_g(ng)%lwup   (:,i,j),m1,wf,fln)
 CALL iofil ('rshort.txt' ,radiate_g(ng)%rshort (i,j),1,wf,fln)
 CALL iofil ('rlong.txt'  ,radiate_g(ng)%rlong  (i,j),1,wf,fln)
 CALL iofil ('rlongup.txt',radiate_g(ng)%rlongup(i,j),1,wf,fln)
 CALL iofil ('albedt.txt' ,radiate_g(ng)%albedt (i,j),1,wf,fln)
 CALL iofil ('cosz.txt'   ,radiate_g(ng)%cosz   (i,j),1,wf,fln)
 CALL iofil ('aodt.txt'   ,radiate_g(ng)%aodt   (i,j),1,wf,fln)
endif

if(iaerosol > 0) then
 CALL iofil ('cn1np.txt',micro_g(ng)%cn1np(:,i,j),m1,wf,fln)
 CALL iofil ('cn1mp.txt',micro_g(ng)%cn1mp(:,i,j),m1,wf,fln)
 CALL iofil ('cn2np.txt',micro_g(ng)%cn2np(:,i,j),m1,wf,fln)
 CALL iofil ('cn2mp.txt',micro_g(ng)%cn2mp(:,i,j),m1,wf,fln)
endif

if(idust > 0) then
 CALL iofil ('md1np.txt',micro_g(ng)%md1np(:,i,j),m1,wf,fln)
 CALL iofil ('md2np.txt',micro_g(ng)%md2np(:,i,j),m1,wf,fln)
 CALL iofil ('md1mp.txt',micro_g(ng)%md1mp(:,i,j),m1,wf,fln)
 CALL iofil ('md2mp.txt',micro_g(ng)%md2mp(:,i,j),m1,wf,fln)
 if(idust == 2) CALL iofil ('dustfrac.txt',micro_g(ng)%dustfrac(i,j),1,wf,fln)
endif

if(iabcarb > 0) then
 CALL iofil ('abc1np.txt',micro_g(ng)%abc1np(:,i,j),m1,wf,fln)
 CALL iofil ('abc2np.txt',micro_g(ng)%abc2np(:,i,j),m1,wf,fln)
 CALL iofil ('abc1mp.txt',micro_g(ng)%abc1mp(:,i,j),m1,wf,fln)
 CALL iofil ('abc2mp.txt',micro_g(ng)%abc2mp(:,i,j),m1,wf,fln)
endif

if(isalt > 0) then
 CALL iofil ('salt_film_np.txt',micro_g(ng)%salt_film_np(:,i,j),m1,wf,fln)
 CALL iofil ('salt_jet_np.txt' ,micro_g(ng)%salt_jet_np (:,i,j),m1,wf,fln)
 CALL iofil ('salt_spum_np.txt',micro_g(ng)%salt_spum_np(:,i,j),m1,wf,fln)
 CALL iofil ('salt_film_mp.txt',micro_g(ng)%salt_film_mp(:,i,j),m1,wf,fln)
 CALL iofil ('salt_jet_mp.txt' ,micro_g(ng)%salt_jet_mp (:,i,j),m1,wf,fln)
 CALL iofil ('salt_spum_mp.txt',micro_g(ng)%salt_spum_mp(:,i,j),m1,wf,fln)
endif

if(level == 2 .or. level == 3) then
 CALL iofil ('rcp.txt',micro_g(ng)%rcp(:,i,j),m1,wf,fln)
endif

if(level == 3) then
 CALL iofil ('pcpg.txt' ,micro_g(ng)%pcpg (i,j),1,wf,fln)
 CALL iofil ('qpcpg.txt',micro_g(ng)%qpcpg(i,j),1,wf,fln)
 CALL iofil ('dpcpg.txt',micro_g(ng)%dpcpg(i,j),1,wf,fln)
 if(ipris>=5 .and. (iifn==1.or.iifn==2)) then
  CALL iofil ('cifnp.txt',micro_g(ng)%cifnp(:,i,j),m1,wf,fln)
 endif
 if(idriz >= 1) then
  CALL iofil ('rdp.txt'  ,micro_g(ng)%rdp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvd.txt',micro_g(ng)%pcpvd(:,i,j),m1,wf,fln)
  CALL iofil ('pcprd.txt',micro_g(ng)%pcprd(i,j),1,wf,fln)
  CALL iofil ('accpd.txt',micro_g(ng)%accpd(i,j),1,wf,fln)
 endif
 if(irain >= 1) then 
  CALL iofil ('rrp.txt'  ,micro_g(ng)%rrp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvr.txt',micro_g(ng)%pcpvr(:,i,j),m1,wf,fln)
  CALL iofil ('q2.txt'   ,micro_g(ng)%q2   (:,i,j),m1,wf,fln)
  CALL iofil ('pcprr.txt',micro_g(ng)%pcprr(i,j),1,wf,fln)
  CALL iofil ('accpr.txt',micro_g(ng)%accpr(i,j),1,wf,fln)
 endif
 if(ipris >= 1) then
  CALL iofil ('rpp.txt'  ,micro_g(ng)%rpp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvp.txt',micro_g(ng)%pcpvp(:,i,j),m1,wf,fln)
  CALL iofil ('pcprp.txt',micro_g(ng)%pcprp(i,j),1,wf,fln)
  CALL iofil ('accpp.txt',micro_g(ng)%accpp(i,j),1,wf,fln)
 endif
 if(isnow >= 1) then
  CALL iofil ('rsp.txt'  ,micro_g(ng)%rsp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvs.txt',micro_g(ng)%pcpvs(:,i,j),m1,wf,fln)
  CALL iofil ('pcprs.txt',micro_g(ng)%pcprs(i,j),1,wf,fln)
  CALL iofil ('accps.txt',micro_g(ng)%accps(i,j),1,wf,fln)
 endif
 if(iaggr >= 1) then
  CALL iofil ('rap.txt'  ,micro_g(ng)%rap  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpva.txt',micro_g(ng)%pcpva(:,i,j),m1,wf,fln)
  CALL iofil ('pcpra.txt',micro_g(ng)%pcpra(i,j),1,wf,fln)
  CALL iofil ('accpa.txt',micro_g(ng)%accpa(i,j),1,wf,fln)
 endif
 if(igraup >= 1) then
  CALL iofil ('rgp.txt'  ,micro_g(ng)%rgp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvg.txt',micro_g(ng)%pcpvg(:,i,j),m1,wf,fln)
  CALL iofil ('q6.txt'   ,micro_g(ng)%q6   (:,i,j),m1,wf,fln)
  CALL iofil ('pcprg.txt',micro_g(ng)%pcprg(i,j),1,wf,fln)
  CALL iofil ('accpg.txt',micro_g(ng)%accpg(i,j),1,wf,fln)
 endif
 if(ihail >=1) then
  CALL iofil ('rhp.txt'  ,micro_g(ng)%rhp  (:,i,j),m1,wf,fln)
  CALL iofil ('pcpvh.txt',micro_g(ng)%pcpvh(:,i,j),m1,wf,fln)
  CALL iofil ('q7.txt'   ,micro_g(ng)%q7   (:,i,j),m1,wf,fln)
  CALL iofil ('pcprh.txt',micro_g(ng)%pcprh(i,j),1,wf,fln)
  CALL iofil ('accph.txt',micro_g(ng)%accph(i,j),1,wf,fln)
 endif
 if(jnmb(1) >= 5) CALL iofil ('ccp.txt',micro_g(ng)%ccp(:,i,j),m1,wf,fln)
 if(jnmb(8) >= 5) CALL iofil ('cdp.txt',micro_g(ng)%cdp(:,i,j),m1,wf,fln)
 if(jnmb(2) >= 5) CALL iofil ('crp.txt',micro_g(ng)%crp(:,i,j),m1,wf,fln)
 if(jnmb(3) >= 5) CALL iofil ('cpp.txt',micro_g(ng)%cpp(:,i,j),m1,wf,fln)
 if(jnmb(4) >= 5) CALL iofil ('csp.txt',micro_g(ng)%csp(:,i,j),m1,wf,fln)
 if(jnmb(5) >= 5) CALL iofil ('cap.txt',micro_g(ng)%cap(:,i,j),m1,wf,fln)
 if(jnmb(6) >= 5) CALL iofil ('cgp.txt',micro_g(ng)%cgp(:,i,j),m1,wf,fln)
 if(jnmb(7) >= 5) CALL iofil ('chp.txt',micro_g(ng)%chp(:,i,j),m1,wf,fln)
 if(iccnlev >= 2) then
  CALL iofil ('regen_aero1_np.txt',micro_g(ng)%regen_aero1_np(:,i,j),m1,wf,fln)
  CALL iofil ('regen_aero1_mp.txt',micro_g(ng)%regen_aero1_mp(:,i,j),m1,wf,fln)
  CALL iofil ('regen_aero2_np.txt',micro_g(ng)%regen_aero2_np(:,i,j),m1,wf,fln)
  CALL iofil ('regen_aero2_mp.txt',micro_g(ng)%regen_aero2_mp(:,i,j),m1,wf,fln)
  if(jnmb(1)>=1) CALL iofil ('cnmcp.txt',micro_g(ng)%cnmcp(:,i,j),m1,wf,fln)
  if(jnmb(8)>=1) CALL iofil ('cnmdp.txt',micro_g(ng)%cnmdp(:,i,j),m1,wf,fln)
  if(jnmb(2)>=1) CALL iofil ('cnmrp.txt',micro_g(ng)%cnmrp(:,i,j),m1,wf,fln)
  if(jnmb(3)>=1) CALL iofil ('cnmpp.txt',micro_g(ng)%cnmpp(:,i,j),m1,wf,fln)
  if(jnmb(4)>=1) CALL iofil ('cnmsp.txt',micro_g(ng)%cnmsp(:,i,j),m1,wf,fln)
  if(jnmb(5)>=1) CALL iofil ('cnmap.txt',micro_g(ng)%cnmap(:,i,j),m1,wf,fln)
  if(jnmb(6)>=1) CALL iofil ('cnmgp.txt',micro_g(ng)%cnmgp(:,i,j),m1,wf,fln)
  if(jnmb(7)>=1) CALL iofil ('cnmhp.txt',micro_g(ng)%cnmhp(:,i,j),m1,wf,fln)
  CALL iofil ('accpaero.txt',micro_g(ng)%accpaero(i,j),1,wf,fln)
  CALL iofil ('pcpraero.txt',micro_g(ng)%pcpraero(i,j),1,wf,fln)
  if(itrkdust==1 .and. idust>0) then
   CALL iofil ('accpdust.txt',micro_g(ng)%accpdust(i,j),1,wf,fln)
   CALL iofil ('pcprdust.txt',micro_g(ng)%pcprdust(i,j),1,wf,fln)
   if(jnmb(1)>=1) CALL iofil ('dnmcp.txt',micro_g(ng)%dnmcp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) CALL iofil ('dnmdp.txt',micro_g(ng)%dnmdp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) CALL iofil ('dnmrp.txt',micro_g(ng)%dnmrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) CALL iofil ('dnmpp.txt',micro_g(ng)%dnmpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) CALL iofil ('dnmsp.txt',micro_g(ng)%dnmsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) CALL iofil ('dnmap.txt',micro_g(ng)%dnmap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) CALL iofil ('dnmgp.txt',micro_g(ng)%dnmgp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) CALL iofil ('dnmhp.txt',micro_g(ng)%dnmhp(:,i,j),m1,wf,fln)
  endif
  if(itrkdustifn==1 .and. idust>0) then
   if(jnmb(1)>=1) CALL iofil ('dincp.txt',micro_g(ng)%dincp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) CALL iofil ('dindp.txt',micro_g(ng)%dindp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) CALL iofil ('dinrp.txt',micro_g(ng)%dinrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) CALL iofil ('dinpp.txt',micro_g(ng)%dinpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) CALL iofil ('dinsp.txt',micro_g(ng)%dinsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) CALL iofil ('dinap.txt',micro_g(ng)%dinap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) CALL iofil ('dingp.txt',micro_g(ng)%dingp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) CALL iofil ('dinhp.txt',micro_g(ng)%dinhp(:,i,j),m1,wf,fln)
  endif
  if(itrkregen==1) then
   CALL iofil ('accpregen.txt',micro_g(ng)%accpregen(i,j),1,wf,fln)
   CALL iofil ('pcprregen.txt',micro_g(ng)%pcprregen(i,j),1,wf,fln)
   if(jnmb(1)>=1) CALL iofil ('rnmcp.txt',micro_g(ng)%rnmcp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) CALL iofil ('rnmdp.txt',micro_g(ng)%rnmdp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) CALL iofil ('rnmrp.txt',micro_g(ng)%rnmrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) CALL iofil ('rnmpp.txt',micro_g(ng)%rnmpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) CALL iofil ('rnmsp.txt',micro_g(ng)%rnmsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) CALL iofil ('rnmap.txt',micro_g(ng)%rnmap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) CALL iofil ('rnmgp.txt',micro_g(ng)%rnmgp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) CALL iofil ('rnmhp.txt',micro_g(ng)%rnmhp(:,i,j),m1,wf,fln)
   if(jnmb(1)>=1) CALL iofil ('rincp.txt',micro_g(ng)%rincp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) CALL iofil ('rindp.txt',micro_g(ng)%rindp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) CALL iofil ('rinrp.txt',micro_g(ng)%rinrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) CALL iofil ('rinpp.txt',micro_g(ng)%rinpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) CALL iofil ('rinsp.txt',micro_g(ng)%rinsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) CALL iofil ('rinap.txt',micro_g(ng)%rinap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) CALL iofil ('ringp.txt',micro_g(ng)%ringp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) CALL iofil ('rinhp.txt',micro_g(ng)%rinhp(:,i,j),m1,wf,fln)
  endif
  if(itrkepsilon==1) then
   CALL iofil ('resol_aero1_mp.txt',micro_g(ng)%resol_aero1_mp(:,i,j),m1,wf,fln)
   CALL iofil ('resol_aero2_mp.txt',micro_g(ng)%resol_aero2_mp(:,i,j),m1,wf,fln)
   if(jnmb(1)>=1) CALL iofil ('snmcp.txt',micro_g(ng)%snmcp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) CALL iofil ('snmdp.txt',micro_g(ng)%snmdp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) CALL iofil ('snmrp.txt',micro_g(ng)%snmrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) CALL iofil ('snmpp.txt',micro_g(ng)%snmpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) CALL iofil ('snmsp.txt',micro_g(ng)%snmsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) CALL iofil ('snmap.txt',micro_g(ng)%snmap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) CALL iofil ('snmgp.txt',micro_g(ng)%snmgp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) CALL iofil ('snmhp.txt',micro_g(ng)%snmhp(:,i,j),m1,wf,fln)
  endif
 endif
 if(iifn==3 .and. iccnlev>=1) then
  if(jnmb(1)>=5) CALL iofil ('ifnnucp.txt',micro_g(ng)%ifnnucp(:,i,j),m1,wf,fln)
  if(jnmb(1)>=5) CALL iofil ('immercp.txt',micro_g(ng)%immercp(:,i,j),m1,wf,fln)
  if(jnmb(8)>=5) CALL iofil ('immerdp.txt',micro_g(ng)%immerdp(:,i,j),m1,wf,fln)
  if(jnmb(2)>=5) CALL iofil ('immerrp.txt',micro_g(ng)%immerrp(:,i,j),m1,wf,fln)
 endif
 if(imbudget>=1) then !18
  CALL iofil ('latheatvapr.txt' ,micro_g(ng)%latheatvapr (:,i,j),m1,wf,fln)
  CALL iofil ('latheatfrzr.txt' ,micro_g(ng)%latheatfrzr (:,i,j),m1,wf,fln)
  CALL iofil ('accnuccldr.txt'   ,micro_g(ng)%accnuccldr   (:,i,j),m1,wf,fln)
  CALL iofil ('acccld2rain.txt'  ,micro_g(ng)%acccld2rain  (:,i,j),m1,wf,fln)
  CALL iofil ('accice2rain.txt'  ,micro_g(ng)%accice2rain  (:,i,j),m1,wf,fln)
  CALL iofil ('accnucicer.txt'   ,micro_g(ng)%accnucicer   (:,i,j),m1,wf,fln)
  CALL iofil ('accvapliq.txt'    ,micro_g(ng)%accvapliq    (:,i,j),m1,wf,fln)
  CALL iofil ('accvapice.txt'    ,micro_g(ng)%accvapice    (:,i,j),m1,wf,fln)
  CALL iofil ('accevapliq.txt'   ,micro_g(ng)%accevapliq   (:,i,j),m1,wf,fln)
  CALL iofil ('accevapice.txt'   ,micro_g(ng)%accevapice   (:,i,j),m1,wf,fln)
  CALL iofil ('accfreezing.txt'  ,micro_g(ng)%accfreezing  (:,i,j),m1,wf,fln)
  CALL iofil ('accmelting.txt'   ,micro_g(ng)%accmelting   (:,i,j),m1,wf,fln)
  CALL iofil ('accmeltice.txt'   ,micro_g(ng)%accmeltice   (:,i,j),m1,wf,fln)
  CALL iofil ('accrimecld.txt'   ,micro_g(ng)%accrimecld   (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2ice.txt'  ,micro_g(ng)%accrain2ice  (:,i,j),m1,wf,fln)
  CALL iofil ('accaggregate.txt' ,micro_g(ng)%accaggregate (:,i,j),m1,wf,fln)
  CALL iofil ('acclatheatvap.txt',micro_g(ng)%acclatheatvap(:,i,j),m1,wf,fln)
  CALL iofil ('acclatheatfrz.txt',micro_g(ng)%acclatheatfrz(:,i,j),m1,wf,fln)
 endif
 if(imbudget>=2) then !37
  CALL iofil ('accinuchomr.txt'    ,micro_g(ng)%accinuchomr    (:,i,j),m1,wf,fln)
  CALL iofil ('accinuccontr.txt'   ,micro_g(ng)%accinuccontr   (:,i,j),m1,wf,fln)
  CALL iofil ('accinucifnr.txt'    ,micro_g(ng)%accinucifnr    (:,i,j),m1,wf,fln)
  CALL iofil ('accinuchazr.txt'    ,micro_g(ng)%accinuchazr    (:,i,j),m1,wf,fln)
  CALL iofil ('accvapcld.txt'      ,micro_g(ng)%accvapcld      (:,i,j),m1,wf,fln)
  CALL iofil ('accvaprain.txt'     ,micro_g(ng)%accvaprain     (:,i,j),m1,wf,fln)
  CALL iofil ('accvappris.txt'     ,micro_g(ng)%accvappris     (:,i,j),m1,wf,fln)
  CALL iofil ('accvapsnow.txt'     ,micro_g(ng)%accvapsnow     (:,i,j),m1,wf,fln)
  CALL iofil ('accvapaggr.txt'     ,micro_g(ng)%accvapaggr     (:,i,j),m1,wf,fln)
  CALL iofil ('accvapgrau.txt'     ,micro_g(ng)%accvapgrau     (:,i,j),m1,wf,fln)
  CALL iofil ('accvaphail.txt'     ,micro_g(ng)%accvaphail     (:,i,j),m1,wf,fln)
  CALL iofil ('accvapdriz.txt'     ,micro_g(ng)%accvapdriz     (:,i,j),m1,wf,fln)
  CALL iofil ('accevapcld.txt'     ,micro_g(ng)%accevapcld     (:,i,j),m1,wf,fln)
  CALL iofil ('accevaprain.txt'    ,micro_g(ng)%accevaprain    (:,i,j),m1,wf,fln)
  CALL iofil ('accevappris.txt'    ,micro_g(ng)%accevappris    (:,i,j),m1,wf,fln)
  CALL iofil ('accevapsnow.txt'    ,micro_g(ng)%accevapsnow    (:,i,j),m1,wf,fln)
  CALL iofil ('accevapaggr.txt'    ,micro_g(ng)%accevapaggr    (:,i,j),m1,wf,fln)
  CALL iofil ('accevapgrau.txt'    ,micro_g(ng)%accevapgrau    (:,i,j),m1,wf,fln)
  CALL iofil ('accevaphail.txt'    ,micro_g(ng)%accevaphail    (:,i,j),m1,wf,fln)
  CALL iofil ('accevapdriz.txt'    ,micro_g(ng)%accevapdriz    (:,i,j),m1,wf,fln)
  CALL iofil ('accmeltpris.txt'    ,micro_g(ng)%accmeltpris    (:,i,j),m1,wf,fln)
  CALL iofil ('accmeltsnow.txt'    ,micro_g(ng)%accmeltsnow    (:,i,j),m1,wf,fln)
  CALL iofil ('accmeltaggr.txt'    ,micro_g(ng)%accmeltaggr    (:,i,j),m1,wf,fln)
  CALL iofil ('accmeltgrau.txt'    ,micro_g(ng)%accmeltgrau    (:,i,j),m1,wf,fln)
  CALL iofil ('accmelthail.txt'    ,micro_g(ng)%accmelthail    (:,i,j),m1,wf,fln)
  CALL iofil ('accrimecldsnow.txt' ,micro_g(ng)%accrimecldsnow (:,i,j),m1,wf,fln)
  CALL iofil ('accrimecldaggr.txt' ,micro_g(ng)%accrimecldaggr (:,i,j),m1,wf,fln)
  CALL iofil ('accrimecldgrau.txt' ,micro_g(ng)%accrimecldgrau (:,i,j),m1,wf,fln)
  CALL iofil ('accrimecldhail.txt' ,micro_g(ng)%accrimecldhail (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2pr.txt'     ,micro_g(ng)%accrain2pr     (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2sn.txt'     ,micro_g(ng)%accrain2sn     (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2ag.txt'     ,micro_g(ng)%accrain2ag     (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2gr.txt'     ,micro_g(ng)%accrain2gr     (:,i,j),m1,wf,fln)
  CALL iofil ('accrain2ha.txt'     ,micro_g(ng)%accrain2ha     (:,i,j),m1,wf,fln)
  CALL iofil ('accaggrselfpris.txt',micro_g(ng)%accaggrselfpris(:,i,j),m1,wf,fln)
  CALL iofil ('accaggrselfsnow.txt',micro_g(ng)%accaggrselfsnow(:,i,j),m1,wf,fln)
  CALL iofil ('accaggrprissnow.txt',micro_g(ng)%accaggrprissnow(:,i,j),m1,wf,fln)
 endif
 if(imbudget==3 .and. idust>0) then !4
  CALL iofil ('accdust1cldr.txt',micro_g(ng)%accdust1cldr(:,i,j),m1,wf,fln)
  CALL iofil ('accdust2cldr.txt',micro_g(ng)%accdust2cldr(:,i,j),m1,wf,fln)
  CALL iofil ('accdust1drzr.txt',micro_g(ng)%accdust1drzr(:,i,j),m1,wf,fln)
  CALL iofil ('accdust2drzr.txt',micro_g(ng)%accdust2drzr(:,i,j),m1,wf,fln)
 endif

endif !if level == 3

return
END SUBROUTINE readwrite_scm
