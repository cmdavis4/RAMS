!##############################################################################
!A. Igel - 5/2013 - Adapting for CLOUD model
!This module adapts the SBM from wrf, originally from A. Khain. 

!##############################################################################
Subroutine micro_bin ()

use node_mod, only:mzp,mxp,myp,ia,iz,ja,jz

implicit none
integer, save :: ncall = 0

!Things to do on the first call to micro_bin
if (ncall == 0) then
  CALL micro_init2 ()
  ncall = 1
endif

!Things to do every time
CALL thermo ()

!Precipitation
CALL precipitation ()

!update rtp, thp - theta, tempc, pressure do not change
CALL update_thermo ()

!Zero Out Budgets
CALL zero_budgets ()

!Main microphysics driver
CALL micro_proc ()

!Turn this step's budget amounts into rates (and add to run totals)
CALL finish_bin_budgets ()

return
END SUBROUTINE micro_bin

!##############################################################################
Subroutine micro_proc ()

use module_hujisbm
use node_mod
use mem_grid
use mem_basic
use mem_micro
use rconstants
use micro_prm
use micphys, only: imbudget,ipris,igraup,ihail,ccn1_max

IMPLICIT NONE
                     
REAL::qv,tabs,pcgs,rhocgs
REAL::theta,exner
! SBM VARIABLES
REAL,DIMENSION (nkr) :: FF1IN,FF3IN,FF4IN,FF5IN,&
                        FF1R,FF3R,FF4R,FF5R,FCCN,FIN,FCCNIN
REAL,DIMENSION (nkr,icemax) :: FF2IN,FF2R

REAL :: SUP2_OLD 
DOUBLE PRECISION :: ss_max=0.003d0
double precision del1in_limit
DOUBLE PRECISION DEL1NR,DEL2NR,DEL12R,DEL12RD,ES1N,ES2N,EW1N,EW1PN
DOUBLE PRECISION DELSUP1,DELSUP2,DELDIV1,DELDIV2
DOUBLE PRECISION TT,QQ,TTA,QQA,PP,QQt
DOUBLE PRECISION DIV1,DIV2,DIV3,DIV4,DEL1IN,DEL2IN,DEL1AD,DEL2AD
!For bulk nucleation
REAL FACTZ,CONCCCN_XZ,CONCDROP
!Flags
INTEGER ISYM1,ISYM2,ISYM3,ISYM4,ISYM5
INTEGER DIFFU
!!! For CCN regeneration
INTEGER aeroregen,fixice
real fccn0(nkr)
real :: ndrop, subtot, frzfract, rndrop, tot_reg  ! for diagnostic CCN
!functions
real :: sum_pris, sum_snow
INTEGER :: k,i,j
integer itimestep, kr, ikl

do kr=1,nkr
     fccn0(kr) =FCCNR0(KR)*1000.*xccn(kr)*ccn1_max
enddo


!-------------BEGIN Loops over Space------------
do k = 1,mzp
   do i = ia,iz
      do j = ja,jz
         exner=(basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) 
         theta=basic_g(ngrid)%theta(k,i,j)
         tabs=theta*exner*cpi
         pcgs = 10. * p00 * (basic_g(ngrid)%pi0(k,i,j) * cpi) ** cpor
         qv=basic_g(ngrid)%rv(k,i,j)
         rhocgs=basic_g(ngrid)%dn0(k,i,j)*0.001

         !Convert distributions from d(mixing ratio)/dln(r) to d(#/cm3)/d(mass)
         ! CCN
         IF (diagCCN) then
            ndrop=0.0
            DO KR=1,NKR
               ndrop=ndrop+COL*micro_g(ngrid)%ffcd(k,i,j,KR)/XL(KR)*rhocgs
            ENDDO
            FCCN(:) = FCCN0(:)
            if (ndrop >= FCCN0(NKR)*COL) then
               subtot=0.0
               do kr = nkr, 1, -1
                  subtot=subtot+FCCN0(kr)*COL
                  FCCN(kr)=0.0
                  if (subtot >= ndrop) then
                     FCCN(kr)= (subtot-ndrop)/COL
                     exit
                  endif
               enddo
            endif                            
         ELSE      ! if (diagCCN) not true, i.e., prognostic             
            DO KR=1,NKR
               FCCN(KR)=micro_g(ngrid)%fncn(k,i,j,KR)*rhocgs/xccn(kr)
            END DO
         ENDIF      ! if (diagCCN)

         !Convert distributions from d(mixing ratio)/dln(r) to
         !d(#/cm3)/d(mass)/mass. This is a meaningless quantity, but it is
         !treated this way through the remainder of the bin micro code.
         !So be it. Adele

         ! LIQUID
         DO KR=1,NKR
            FF1R(KR)=micro_g(ngrid)%ffcd(k,i,j,KR)*rhocgs/xl(kr)/xl(kr)/3.0
         END DO
        
         ! ICE
         IF (ICEPROCS.EQ.1)THEN
            DO KR=1,NKR
               !ICE NUCLEI
               IF (ICEFLAG .eq. 1) THEN
                  FIN(KR)=micro_g(ngrid)%ffin(k,i,j,KR)*rhocgs/xccn(kr)
               ELSE
                  FIN(KR)=0.
               ENDIF
               ! COLUMNS 
               IF (IPRIS == 1 .or. IPRIS >= 4) THEN
                  FF2R(KR,1)=micro_g(ngrid)%ffic(k,i,j,KR)*rhocgs/xi(kr,1)/xi(kr,1)/3.0
               ELSE
                  FF2R(KR,1) = 0.
               ENDIF

               ! PLATES   
               IF (IPRIS == 2 .or. IPRIS >= 4) THEN
                  FF2R(KR,2)=micro_g(ngrid)%ffip(k,i,j,KR)*rhocgs/xi(kr,2)/xi(kr,2)/3.0
               ELSE
                  FF2R(KR,2) = 0.
               ENDIF

               ! DENDRITES
               IF (IPRIS >= 3) THEN
                  FF2R(KR,3)=micro_g(ngrid)%ffid(k,i,j,KR)*rhocgs/xi(kr,3)/xi(kr,3)/3.0
               ELSE
                  FF2R(KR,3) = 0.
               ENDIF

               ! SNOW/RAMS AGGREGATES
               FF3R(KR)=micro_g(ngrid)%ffsn(k,i,j,KR)*rhocgs/xs(kr)/xs(kr)/3.0

               ! GRAUPEL
               IF (IGRAUP > 0) THEN
                  FF4R(KR)=micro_g(ngrid)%ffgl(k,i,j,KR)*rhocgs/xg(kr)/xg(kr)/3.0
               ELSE
                  FF4R(KR) = 0.
               ENDIF

               ! HAIL
               IF (IHAIL > 0) THEN
                  FF5R(KR)=micro_g(ngrid)%ffhl(k,i,j,KR)*rhocgs/xh(kr)/xh(kr)/3.0
               ELSE
                  FF5R(KR) = 0.
               ENDIF
             END DO
!--------------------------------------------------------
!Freezing (homogeneous nucleation)
             IF (ICEPROCS .eq. 1) THEN
                CALL lhf_budget (FF1R,EXNER,RHOCGS,K,I,J,'beg')
                CALL FREEZ (FF1R,FF2R,FF5R, &
                   tabs,dtlt,rhocgs, &
                   K,I,J)
                IF (ifastsbm==1) CALL fastsbm_reassign (FF2R,FF3R,FF4R,FF5R)
!--------------------------------------------------------
!Melting
                IF(TABS > 273.15) CALL MELT (FF1R,FF2R,FF3R,FF4R,FF5R, &
                   tabs,rhocgs,K,I,J)
            !No need to call fastsbm_reassign - nothing is created in MELT
                CALL lhf_budget (FF1R,EXNER,RHOCGS,K,I,J,'end')
            ENDIF
         ENDIF  ! if iceprocess == 1
!-----------------------------------------------
!Condensation and Nucleation Preparation
         IF (micro_g(ngrid)%t_old(k,i,j).GT.233)THEN   
            TT=micro_g(ngrid)%t_old(k,i,j)
            QQ=micro_g(ngrid)%RV_OLD(k,i,j)/(1+micro_g(ngrid)%RV_OLD(k,i,j))
            PP=pcgs
            TTA=tabs
            QQA=qv/(1+qv)
            !old values, 1 liquid, 2 ice
            ES1N=AA1_MY*DEXP(-BB1_MY/TT)
            ES2N=AA2_MY*DEXP(-BB2_MY/TT)
            EW1N=QQ*PP/(0.622+0.378*QQ)
            DIV1=EW1N/ES1N        !Saturation ratio wrt liquid
            DEL1IN=EW1N/ES1N-1.
            DIV2=EW1N/ES2N        !Saturation ratio wrt ice
            DEL2IN=EW1N/ES2N-1.
            !new values
            ES1N=AA1_MY*DEXP(-BB1_MY/TTA)
            ES2N=AA2_MY*DEXP(-BB2_MY/TTA)
            EW1N=QQA*PP/(0.622+0.378*QQA)
            DIV3=EW1N/ES1N
            DEL1AD=EW1N/ES1N-1.
            DIV4=EW1N/ES2N
            DEL2AD=EW1N/ES2N-1.
            SUP2_OLD=DEL2IN
            !local change
            DELSUP1=(DEL1AD-DEL1IN)/NCOND
            DELSUP2=(DEL2AD-DEL2IN)/NCOND
            DELDIV1=(DIV3-DIV1)/NCOND
            DELDIV2=(DIV4-DIV2)/NCOND
!---------------------------------------------------
            !Latent heating budget prep
            CALL lhv_budget (FF1R,FF2R,FF3R,FF4R,FF5R, &
                            EXNER,RHOCGS,K,I,J,'beg')

!Start nucleation/condensation loop
            DO IKL=1,NCOND

               DEL1IN=DEL1IN+DELSUP1
               DEL2IN=DEL2IN+DELSUP2
               DIV1=DIV1+DELDIV1
               DIV2=DIV2+DELDIV2
               DIFFU=1
               IF (DIV1.GT.DIV2.AND.TT.LE.265)THEN
                  DIV2=0.99999*DIV1
                  DEL2IN=0.99999*DEL2IN
                  DIFFU=0
               END IF
               DEL1NR=A1_MYN*(100.*DIV1)
               DEL2NR=A2_MYN*(100.*DIV2)
               if (DEL2NR.EQ.0) then
                  print*, 'DEL2NR',tt,qq, tta,qqa
                  PRINT*,'DEL2NR = 0',k,i,j
                  STOP
               endif
               DEL12R=DEL1NR/DEL2NR
               DEL12RD=DEL12R**DEL_BBR
               EW1PN=AA1_MY*100.*DIV1*DEL12RD/100.
               IF (DEL12R.EQ.0)PRINT*,'DEL12R = 0'
               IF (DEL12R.EQ.0)STOP
               TT=-DEL_BB/DLOG(DEL12R)
               QQ=0.622*EW1PN/(PP-0.378*EW1PN)
               FF1IN=FF1R
               FF2IN=FF2R
               FF3IN=FF3R
               FF4IN=FF4R
               FF5IN=FF5R
!-------------------------------------------------
!Nucleation
               IF (BULKNUC.eq.1)THEN
                  IF (DEL1IN.GT.0)THEN
                     IF (zt(k).LE.500.)THEN
                        FACTZ=0.
                     ELSE
                        FACTZ=1
                     END IF
                  CONCCCN_XZ=FACTZ*ACCN*(100.*DEL1IN)**BCCN
                  CONCDROP = sum(FF1IN*XL)*3.d0*col
            
                  IF(CONCCCN_XZ.GT.CONCDROP) &
                     FF1IN(1)=FF1IN(1)+(CONCCCN_XZ-CONCDROP)/(3.D0*COL*XL(1))
                  END IF
               !If supersat wrt liq or ice
               ELSEIF (DEL1IN.GT.0.OR.(ICEPROCS.EQ.1.and.DEL2IN.GT.0))THEN
                     !mo Limit SS for the droplet nucleation
                     del1in_limit = min(DEL1IN,ss_max)
                     !NUCLEATION
                     CALL JERNUCL01 (FF1IN,FF2IN,FCCN &
                          ,TT,RHOCGS &
                          ,DEL1IN_limit,DEL2IN &
                          ,SUP2_OLD &
                          ,DTLT,FIN &
                          ,K,I,J)
                     IF (ifastsbm==1) CALL fastsbm_reassign (FF2R,FF3R,FF4R,FF5R)
               END IF !If bulknuc

!-------------------------------------------------------------------
!Condensation/Deposition/Evaporation/Sublimation
               FF1R=FF1IN
               FF2R=FF2IN

               !FF3IN, FF4IN, FF5IN have not changed
               !Flags to determine if deposition/sublimation is necessary
               ISYM1=0;ISYM2=0;ISYM3=0;ISYM4=0;ISYM5=0
               if (any(FF1R.gt.1.e-6))ISYM1=1
               if (iceprocs.eq.1) then
                  if(any(FF2R.gt.1.e-6))ISYM2=1
                  if(any(FF3R.gt.1.e-6))ISYM3=1
                  if(any(FF4R.gt.1.e-6))ISYM4=1
                  if(any(FF5R.gt.1.e-6))ISYM5=1
               endif
    
               ccnreg=0.
               inreg=0.

               IF (DIFFU.NE.0 .and. (ISYM1.eq.1.or.ISYM2.eq.1.or.ISYM3.eq.1 &
                                 .or.ISYM4.eq.1.or.ISYM5.eq.1))THEN 
                  !DIFFU = 0 if too cold and RH>RHi (does this make sense?)
                  CALL vap_budget (FF1R,FF2R,FF3R,FF4R,FF5R,RHOCGS,K,I,J,'beg')

                  !If liquid
                  IF (ISYM1.EQ.1.AND.((TT-273.15).GT.-0.187.OR. &
                     (ISYM2.EQ.0.AND.ISYM3.EQ.0.AND.ISYM4.EQ.0.AND.ISYM5.EQ.0)))THEN
                     CALL ONECOND1 (TT,QQ,PP,rhocgs &
                          ,VR1,pcgs &
                          ,DEL1IN,DEL2IN,DIV1,DIV2 &
                          ,FF1R,FF1IN,XL,RLEC,RO1BL &
                          ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                          ,COL,DTCOND,ICEMAX,NKR &
                          ,micro_g(ngrid)%vapcldr(k,i,j),micro_g(ngrid)%vaprainr(k,i,j))
                  !If ice
                  ELSE IF(ISYM1.EQ.0.AND.(TT-273.15).LE.-0.187.AND. &
                     (ISYM2.EQ.1.OR.ISYM3.EQ.1.OR.ISYM4.EQ.1.OR.ISYM5.EQ.1))THEN
                     CALL ONECOND2 (TT,QQ,PP,rhocgs &
                          ,VR2,VR3,VR4,VR5,pcgs &
                          ,DEL1IN,DEL2IN,DIV1,DIV2 &
                          ,FF2R,FF2IN,XI,RIEC,RO2BL &
                          ,FF3R,FF3IN,XS,RSEC,RO3BL &
                          ,FF4R,FF4IN,XG,RGEC,RO4BL &
                          ,FF5R,FF5IN,XH,RHEC,RO5BL &
                          ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                          ,COL,DTCOND,ICEMAX,NKR &
                          ,ISYM2,ISYM3,ISYM4,ISYM5)
                  !If mixed-phase
                  ELSE IF(ISYM1.EQ.1.AND.(TT-273.15).LE.-0.187.AND. &
                     (ISYM2.EQ.1.OR.ISYM3.EQ.1.OR.ISYM4.EQ.1.OR.ISYM5.EQ.1))THEN
                     CALL ONECOND3 (TT,QQ,PP,rhocgs &
                          ,VR1,VR2,VR3,VR4,VR5,pcgs &
                          ,DEL1IN,DEL2IN,DIV1,DIV2 &
                          ,FF1R,FF1IN,XL,RLEC,RO1BL &
                          ,FF2R,FF2IN,XI,RIEC,RO2BL &
                          ,FF3R,FF3IN,XS,RSEC,RO3BL &
                          ,FF4R,FF4IN,XG,RGEC,RO4BL &
                          ,FF5R,FF5IN,XH,RHEC,RO5BL &
                          ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
                          ,COL,DTCOND,ICEMAX,NKR &
                          ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5 &
                          ,micro_g(ngrid)%vapcldr(k,i,j),micro_g(ngrid)%vaprainr(k,i,j))
                  END IF
                  IF (ifastsbm==1) CALL fastsbm_reassign (FF2R,FF3R,FF4R,FF5R)

                  CALL vap_budget (FF1R,FF2R,FF3R,FF4R,FF5R,RHOCGS,K,I,J,'end')
               END IF

!-----------------------------------------------------------------------
!CCN and IN regeneration
! 1. Saleeby(2020):This section of regeneration should be checked.
! 2. Should the line with "rndropr(k,i,j) = rndrop" be commented out???
! 3. "fixice" was not a standard flag here. I added it, but am not sure of its usage.
 aeroregen=0
 fixice=1
 if(aeroregen==1)then
       !For ccn regeneration from evaporation (J. Fan Oct 2007)
       if (fixice .ne.1 .and. iceflag == 1) then 
          ! For drop freezing from evaporation (J. Fan Oct 2007)
          if (iceform == 1) then
             if (ccnreg > 0.0 .and. TT < (273.15-5.0)) then
                frzfract = 0.8e-5
                CALL evapfrz (ccnreg,frzfract,NKR,ICEMAX,TT,dtlt,xi,ff2r,rndrop)
                !rndropr(k,i,j) = rndrop
             endif
             if (inreg > 0.) then
                fin(nkr)=fin(nkr)+inreg/col
             endif
          ! part of drop evaporating residuals back to inreg
          else if (iceform == 2) then
             inreg= inreg + ccnreg*0.5e-5
             fin(nkr)=fin(nkr)+inreg/col
          endif
       endif   !if (.not.fixice .and. iceflag == 1)
        
       !Put ccnreg back to aerosol (CCN regeneration) if diagCCN = false
       if (diagCCN .eqv. .false.) then
           tot_reg  = ccnreg
           if (tot_reg > 0.0) then
              fccnin(:) = fccn(:)  
              CALL ccn_reg (fccn0,fccnin,fccn,nkr,tot_reg,ff1r,xl,  &
                   ff2r,xi,ff3r,xs,ff4r,xg)
           endif
        endif

 endif !aeroregen

            END DO    ! end NCOND
            !Finish latent heating budget
            CALL lhv_budget (FF1R,FF2R,FF3R,FF4R,FF5R, &
                            EXNER,RHOCGS,K,I,J,'end')

!---------------------------------------------------------------------------
!Collisions (called every NDTCOLL*dt)
            IF (mod(istp,NDTCOLL).eq.0 .and. ((iceprocs.eq.0 .and. any(FF1R>0.)) &
                 .or. any(FF1R>0.).or.any(FF2R>0.).or.any(FF3R>0.)))  THEN
               CALL lhf_budget (FF1R,EXNER,RHOCGS,K,I,J,'beg')
               !Graupel and hail do not self-collect, so don't need to call 
               !collection if they are the only species present
               CALL COAL_BOTT_NEW (FF1R,FF2R,FF3R, &
                FF4R,FF5R,TT,QQ,PP,rhocgs,DTCOLL,TCRIT,TTCOAL,K,I,J)
               IF (ifastsbm==1) CALL fastsbm_reassign (FF2R,FF3R,FF4R,FF5R)
               CALL lhf_budget (FF1R,EXNER,RHOCGS,K,I,J,'end')
            END IF
!---------------------------------------------------------------------------
!Updates

            micro_g(ngrid)%t_old(k,i,j)= tt
            micro_g(ngrid)%rv_old(k,i,j)=qq/(1-qq)
         END IF  !If t_old > 233

         ! Convert back to f(diam)
         DO KR=1,NKR
            micro_g(ngrid)%fncn(k,i,j,KR)=FCCN(KR)/rhocgs*xccn(kr)
            micro_g(ngrid)%ffcd(k,i,j,KR)=FF1R(KR)/rhocgs*xl(kr)*xl(kr)*3.0
         END DO 
  
         IF (ICEPROCS.EQ.1)THEN
            DO KR=1,NKR
              if (iceflag .eq. 1) &
                 micro_g(ngrid)%ffin(k,i,j,KR)=FIN(KR)/rhocgs*xccn(kr)
              if (ipris == 1 .or. ipris >= 4) &
                 micro_g(ngrid)%ffic(k,i,j,KR)=FF2R(KR,1)/rhocgs*xi(kr,1)*xi(kr,1)*3.0
              if (ipris == 2 .or. ipris >= 4) &
                 micro_g(ngrid)%ffip(k,i,j,KR)=FF2R(KR,2)/rhocgs*xi(kr,2)*xi(kr,2)*3.0
              if (ipris >= 3) &
                 micro_g(ngrid)%ffid(k,i,j,KR)=FF2R(KR,3)/rhocgs*xi(kr,3)*xi(kr,3)*3.0
              micro_g(ngrid)%ffsn(k,i,j,KR)=FF3R(KR)/rhocgs*xs(kr)*xs(kr)*3.0
              if (igraup > 0) &
                 micro_g(ngrid)%ffgl(k,i,j,KR)=FF4R(KR)/rhocgs*xg(kr)*xg(kr)*3.0
              if (ihail > 0) &
                 micro_g(ngrid)%ffhl(k,i,j,KR)=FF5R(KR)/rhocgs*xh(kr)*xh(kr)*3.0
            END DO
         END IF
      END DO !End loops over space
   END DO
END DO


RETURN
END SUBROUTINE micro_proc

!##############################################################################
Subroutine micro_init (bininit)

use node_mod
use micro_prm
use mem_micro
use mem_grid
use micphys, only:ccn1_max,cin_max,iaero_chem,aero_medrad,iifn
use module_hujisbm

IMPLICIT NONE

character(len=strl1) :: hname
INTEGER k,i,j,KR,bininit
INTEGER,parameter :: hujisbm_unit1 = 22
REAL PI
data pi/3.141592654/
! dtime - timestep of integration (calculated in main program) :
! ax - coefficient used for masses calculation 
! ima(i,j) - k-category number, c(i,j) - courant number 
REAL X0DROP,DEG01,X0CCN

if (my_mpi_num .eq. 0) PRINT*, 'INITIALIZING HUCM'  
if (my_mpi_num .eq. 0) print *, ' ****** HUCM *******'

dlnr=dlog(2.d0)/(3.d0*scal)

!--- Read in various lookup tables

900   FORMAT(6E13.5)
! MASSES
hname=trim(hucmfile)//'/masses.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,900) XL,XI,XS,XG,XH          
CLOSE(hujisbm_unit1)

! BULKRADIUS
hname=trim(hucmfile)//'/bulkradii.asc_s_0_03_0_9'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,*) RADXXO
CLOSE(hujisbm_unit1)

! INITIALIZE AEROSOLS
DEG01=1./3.
X0DROP=XL(ICCN)             !Maximum ccn mass (minimum drop mass)
X0CCN =X0DROP/(2.**(NKR-1)) !Minimum ccn mass

DROPRADII = RADXXO(:,1)

!Choose aerosol type
IF (IAERO_CHEM(1) == 1) THEN
!Ammonium sulfate
   ROCCN0 = 1.769
   RO_SOLUTE = 1.769
   MWAERO = 132.0
   IONS = 3
   RP_G = AERO_MEDRAD(1)*100.
   SIG_G = 1.8
ELSEIF (IAERO_CHEM(1) == 2) THEN
!Sodium chloride
   ROCCN0 = 2.165
   RO_SOLUTE = 2.165
   MWAERO = 58.4
   IONS = 2
   RP_G = AERO_MEDRAD(1)*100.
   SIG_G = 1.8
ELSE
   PRINT*, 'INVALID AERO_CHEM TYPE. MUST EQUAL 1 OR 2', IAERO_CHEM(1)
   STOP
ENDIF

!Find CCN density, mass, and diameter for each bin
DO KR=1,NKR
   ROCCN(KR)=ROCCN0
   XCCN(KR)=X0CCN*2.**(KR-1)
   RCCN(KR)=(3.*XCCN(KR)/4./3.141593/ROCCN(KR))**DEG01    
ENDDO

! Lognormal functional Form
DO KR=1,NKR
  fccnr0(kr) = 1./(sqrt(2.*pi)*log(sig_g))    &
               *exp(-(log(rccn(kr)/rp_g))**2./2.0/(log(sig_g))**2.)
ENDDO
fccnr0 = fccnr0 / (sum(fccnr0) * col)

if (my_mpi_num .eq. 0) then
   PRINT *, '*********  DROP RADII *******'
   PRINT 200, DROPRADII
   PRINT *, '*********  CCN RADII *******'
   PRINT 200, RCCN
   PRINT *, '********* CCN MASSES *******'
   PRINT 200, XCCN
   PRINT *, '********* INITIAL NORMALIZED CCN DISTRIBUTION *******'
   PRINT 200, FCCNR0
endif
200   FORMAT(6E13.5)

! Initialize arrays for the initial run. Do not run this for history restart.
if(bininit==1) then
!Factor of 1000 is to convert to CGS units. ccn1_max is #/mg
   DO k=1,mzp
      DO i=1,mxp
         DO j=1,myp
            DO KR=1,NKR
               IF (k<=2) THEN
                  micro_g(ngrid)%fncn(k,i,j,KR)=FCCNR0(KR)*1000.*xccn(kr)*ccn1_max
               ELSE
                  micro_g(ngrid)%fncn(k,i,j,KR)=FCCNR0(KR)*1000.*xccn(kr)* &
                                                ccn1_max*exp(-zt(k)/7000.)
               ENDIF
            ENDDO
         ENDDO
      END DO
   END DO
   if (iceprocs == 1 .and. iifn == 2) then
      fracin = cin_max/ccn1_max
      micro_g(ngrid)%ffin=micro_g(ngrid)%fncn*fracin
   endif
endif

return
END SUBROUTINE micro_init

!##############################################################################
Subroutine micro_init2 ()

use node_mod, only: my_mpi_num
use mem_grid, only: dtlt,ngrid,initial,runtype
use rconstants, only: cpi
use micro_prm
use mem_micro
use mem_basic
use micphys, only:iifn,ipris,igraup,ihail
use module_hujisbm

IMPLICIT NONE

character(len=strl1) :: hname
INTEGER KR
INTEGER,parameter :: hujisbm_unit1 = 22
REAL PI,NDTCOLL_REAL
data PI/3.141592654/

!--- Read in various lookup tables
! CAPACITIES :
! Capacities are used for the condensation rates
hname=trim(hucmfile)//'/capacity.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
900   FORMAT(6E13.5)
READ(hujisbm_unit1,900) RLEC,RIEC,RSEC,RGEC,RHEC
CLOSE(hujisbm_unit1)

! MASSES :
hname=trim(hucmfile)//'/masses.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,900) XL,XI,XS,XG,XH          
CLOSE(hujisbm_unit1)


! TERMINAL VELOCITY :
hname=trim(hucmfile)//'/termvels.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,*) VR1,VR2,VR3,VR4,VR5     
CLOSE(hujisbm_unit1)

! CONSTANTS :   !Adele - none of these constants are used anywhere
!OPEN(UNIT=hujisbm_unit1,FILE="../etc/HUCM-SBM/constants.asc", &
!     FORM="FORMATTED",STATUS="OLD")
!READ(hujisbm_unit1,900) SLIC,TLIC,COEFIN,C2,C3,C4
!CLOSE(hujisbm_unit1)

! KERNELS DEPENDING ON PRESSURE :
hname=trim(hucmfile)//'/kernels_z.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,900) YWLL_1000MB,YWLL_750MB,YWLL_500MB
CLOSE(hujisbm_unit1)

! KERNELS NOT DEPENDING ON PRESSURE :
hname=trim(hucmfile)//'/kernels.asc_s_0_03_0_9'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,900) &
   YWLL,YWLI,YWLS,YWLG,YWLH, &
   YWIL,YWII,YWIS,YWIG,YWIH, &
   YWSL,YWSI,YWSS,YWSG,YWSH, &
   YWGL,YWGI,YWGS,YWGG,YWGH, &
   YWHL,YWHI,YWHS,YWHG,YWHH
close (hujisbm_unit1)

! BULKDENSITY :
hname=trim(hucmfile)//'/bulkdens.asc_s_0_03_0_9'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,900) RO1BL,RO2BL,RO3BL,RO4BL,RO5BL
CLOSE(hujisbm_unit1)

! BULKRADIUS
hname=trim(hucmfile)//'/bulkradii.asc_s_0_03_0_9'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
READ(hujisbm_unit1,*) RADXXO
CLOSE(hujisbm_unit1)

do kr=1,nkr
   xl_mg(kr)=xl(kr)*1.e3
   xs_mg(kr)=xs(kr)*1.e3
   xg_mg(kr)=xg(kr)*1.e3
   xh_mg(kr)=xh(kr)*1.e3
   xi1_mg(kr)=xi(kr,1)*1.e3
   xi2_mg(kr)=xi(kr,2)*1.e3
   xi3_mg(kr)=xi(kr,3)*1.e3
enddo
CALL courant_bott ()

!Initialize collision-coalescence parameters
CALL BREAKINIT ()
! Collisions are called every NDTCOLL*dt
DTCOLL=dtlt*real(NDTCOLL)
CALL kernals (DTCOLL)

!Initialize some other parameters
DTCOND=dtlt/REAL(NCOND)
ICEFLAG = IIFN-1
IFASTSBM = 0
if (IPRIS < 4 .or. IGRAUP == 0 .or. IHAIL == 0) IFASTSBM = 1

if (iceprocs.eq.1 .and. my_mpi_num.eq.0)print*,'ICE PROCESSES ACTIVE'
if (iceprocs.eq.0 .and. my_mpi_num.eq.0)print*,'LIQUID PROCESSES ONLY'

!Initialize t_old and rv_old if not history initialization or restart
if (initial .ne. 3) then
   micro_g(ngrid)%t_old=basic_g(ngrid)%theta * &
      (basic_g(ngrid)%pi0(:,:,:)+basic_g(ngrid)%pp(:,:,:)) * cpi
   micro_g(ngrid)%rv_old=basic_g(ngrid)%rv
endif

return
END SUBROUTINE micro_init2

!##############################################################################
Subroutine sum_bins (speciesbin,rx,m1,m2,m3,krs,kre,izero)

use micro_prm, only:nkr, col
use node_mod, only:ia,iz,ja,jz

implicit none

integer :: i,j,k,kr,krs,kre,m1,m2,m3,izero
real, dimension(m1,m2,m3,nkr) :: speciesbin
real, dimension(m1,m2,m3) :: rx

if (izero==0) then
   rx=0.
   izero=1
endif

do kr=krs,kre
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            rx(k,i,j)=rx(k,i,j)+speciesbin(k,i,j,kr)*col
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE sum_bins

!##############################################################################
Subroutine sum_bins_conc (speciesbin,cx,xsp,m1,m2,m3,krs,kre,izero)

use micro_prm, only:nkr, col
use node_mod, only:ia,iz,ja,jz

implicit none

integer :: i,j,k,kr,krs,kre,m1,m2,m3,izero
real, dimension(m1,m2,m3,nkr) :: speciesbin
real, dimension(m1,m2,m3) :: cx
real, dimension(nkr) :: xsp

if (izero==0) then
   cx=0.
   izero = 1
endif
do k = 1,m1
   do j = 1,m3
      do i = 1,m2
         do kr=krs,kre
            cx(k,i,j)=cx(k,i,j)+speciesbin(k,i,j,kr)*col/xsp(kr)*1000.
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE sum_bins_conc

!##############################################################################
real Function sum_pris (ff2r,rhocgs)

use micro_prm, only:nkr, icemax, col3, krpris
use module_hujisbm, only:xi,sum_mass
use micphys, only:ipris

implicit none
real :: rhocgs
real, dimension(nkr,icemax)::ff2r

IF (IPRIS >=4) THEN
    sum_pris = sum_mass(ff2r(:,1),XI(:,1),rhocgs,1,KRPRIS(1)-1) &
              + sum_mass(ff2r(:,2),XI(:,2),rhocgs,1,KRPRIS(2)-1) &
              + sum_mass(ff2r(:,3),XI(:,3),rhocgs,1,KRPRIS(3)-1)

ELSEIF (IPRIS==1 .or. IPRIS==2 .or. IPRIS==3) THEN
  sum_pris = sum_mass(ff2r(:,IPRIS),XI(:,IPRIS),rhocgs,1,KRPRIS(IPRIS)-1)
ENDIF

END FUNCTION sum_pris

!##############################################################################
real Function sum_snow (ff2r,rhocgs)

use micro_prm, only:nkr, icemax, col3, krpris
use module_hujisbm, only:xi,sum_mass
use micphys, only:ipris

implicit none
real:: rhocgs
real, dimension(nkr,icemax)::ff2r

IF (IPRIS >=4) THEN
    sum_snow = sum_mass(ff2r(:,1),XI(:,1),rhocgs,KRPRIS(1),NKR) &
             + sum_mass(ff2r(:,2),XI(:,2),rhocgs,KRPRIS(2),NKR) &
             + sum_mass(ff2r(:,3),XI(:,3),rhocgs,KRPRIS(3),NKR)

ELSEIF (IPRIS==1 .or. IPRIS==2 .or. IPRIS==3) THEN
  sum_snow = sum_mass(ff2r(:,IPRIS),XI(:,IPRIS),rhocgs,KRPRIS(IPRIS),NKR)
ENDIF

END FUNCTION sum_snow

!##############################################################################
Subroutine precipitation ()

use module_hujisbm, only: vr1,vr2,vr3,vr4,vr5
use mem_micro
use mem_grid 
use mem_basic
use node_mod, only:mmxp,mmyp
use micro_prm, only:iceprocs,nkr,krpris
use micphys, only:ipris,igraup,ihail
use rconstants, only:cpi

implicit none
integer :: icecat
real,dimension(mmxp(ngrid),mmyp(ngrid)) :: tempc

!Precipitation (fall speeds are in cm/s)

!Zero out arrays 
micro_g(ngrid)%pcpg = 0.
micro_g(ngrid)%qpcpg = 0.
micro_g(ngrid)%dpcpg = 0.

!Calculate surface temperature
tempc = basic_g(ngrid)%theta(2,:,:) * &
        (basic_g(ngrid)%pi0(2,:,:) + basic_g(ngrid)%pp(2,:,:)) * cpi &
        - 273.15

!Liquid precipitation
micro_g(ngrid)%pcprr = 0.
micro_g(ngrid)%pcpvr = 0.

CALL sedim_bin (vr1,1,micro_g(ngrid)%ffcd,micro_g(ngrid)%accpr, &
               micro_g(ngrid)%pcpvr,micro_g(ngrid)%pcprr, &
               tempc,1,nkr)

!Ice
if (iceprocs .eq. 1) then
   if (ipris >= 1) then
      micro_g(ngrid)%pcprp = 0.
      micro_g(ngrid)%pcpvp = 0.
      micro_g(ngrid)%pcprs = 0.
      micro_g(ngrid)%pcpvs = 0.
   endif

   if (ipris == 1 .or. ipris >= 4) then
      icecat = 1
      micro_g(ngrid)%accpic = micro_g(ngrid)%accpic - &
                              micro_g(ngrid)%accpp - micro_g(ngrid)%accps
      micro_g(ngrid)%pcpric = 0. - micro_g(ngrid)%pcprp - micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvic = 0. - micro_g(ngrid)%pcpvp - micro_g(ngrid)%pcpvs

      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffic,micro_g(ngrid)%accpp, &
                     micro_g(ngrid)%pcpvp,micro_g(ngrid)%pcprp, &
                     tempc,1,krpris(icecat)-1)
      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffic,micro_g(ngrid)%accps, &
                     micro_g(ngrid)%pcpvs,micro_g(ngrid)%pcprs, &
                     tempc,krpris(icecat),nkr)

      micro_g(ngrid)%accpic = micro_g(ngrid)%accpic + &
                              micro_g(ngrid)%accpp + micro_g(ngrid)%accps
      micro_g(ngrid)%pcpric = micro_g(ngrid)%pcpric + micro_g(ngrid)%pcprp + micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvic = micro_g(ngrid)%pcpvic + micro_g(ngrid)%pcpvp + micro_g(ngrid)%pcpvs
   endif

   if (ipris == 2 .or. ipris >= 4) then
      icecat = 2
      micro_g(ngrid)%accpip = micro_g(ngrid)%accpip - &
                              micro_g(ngrid)%accpp - micro_g(ngrid)%accps
      micro_g(ngrid)%pcprip = 0. - micro_g(ngrid)%pcprp - micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvip = 0. - micro_g(ngrid)%pcpvp - micro_g(ngrid)%pcpvs

      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffip,micro_g(ngrid)%accpp, &
                     micro_g(ngrid)%pcpvp,micro_g(ngrid)%pcprp, &
                     tempc,1,krpris(icecat)-1)
      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffip,micro_g(ngrid)%accps, &
                     micro_g(ngrid)%pcpvs,micro_g(ngrid)%pcprs, &
                     tempc,krpris(icecat),nkr)
      micro_g(ngrid)%accpip = micro_g(ngrid)%accpip + &
                              micro_g(ngrid)%accpp + micro_g(ngrid)%accps
      micro_g(ngrid)%pcprip = micro_g(ngrid)%pcprip + micro_g(ngrid)%pcprp + micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvip = micro_g(ngrid)%pcpvip + micro_g(ngrid)%pcpvp + micro_g(ngrid)%pcpvs
   endif

   if (ipris >= 3) then
      icecat = 3
      micro_g(ngrid)%accpid = micro_g(ngrid)%accpid - &
                              micro_g(ngrid)%accpp - micro_g(ngrid)%accps
      micro_g(ngrid)%pcprid = 0. - micro_g(ngrid)%pcprp - micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvid = 0. - micro_g(ngrid)%pcpvp - micro_g(ngrid)%pcpvs

      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffid,micro_g(ngrid)%accpp, &
                     micro_g(ngrid)%pcpvp,micro_g(ngrid)%pcprp, &
                     tempc,1,krpris(icecat)-1)
      CALL sedim_bin (vr2(:,icecat),icecat+1,micro_g(ngrid)%ffid,micro_g(ngrid)%accps, &
                     micro_g(ngrid)%pcpvs,micro_g(ngrid)%pcprs, &
                     tempc,krpris(icecat),nkr)
      micro_g(ngrid)%accpid = micro_g(ngrid)%accpid + &
                              micro_g(ngrid)%accpp + micro_g(ngrid)%accps
      micro_g(ngrid)%pcprid = micro_g(ngrid)%pcprid + micro_g(ngrid)%pcprp + micro_g(ngrid)%pcprs
      micro_g(ngrid)%pcpvid = micro_g(ngrid)%pcpvid + micro_g(ngrid)%pcpvp + micro_g(ngrid)%pcpvs
   endif

   micro_g(ngrid)%pcpra = 0.
   micro_g(ngrid)%pcpva = 0.
   CALL sedim_bin (vr3,5,micro_g(ngrid)%ffsn,micro_g(ngrid)%accpa, &
                  micro_g(ngrid)%pcpva,micro_g(ngrid)%pcpra, &
                  tempc,1,nkr)

   if (igraup > 0) then
      micro_g(ngrid)%pcprg = 0.
      micro_g(ngrid)%pcpvg = 0.
      CALL sedim_bin (vr4,6,micro_g(ngrid)%ffgl,micro_g(ngrid)%accpg, &
                     micro_g(ngrid)%pcpvg,micro_g(ngrid)%pcprg, &
                     tempc,1,nkr)
   endif

   if (ihail > 0) then
      micro_g(ngrid)%pcprh = 0.
      micro_g(ngrid)%pcpvh = 0.
      CALL sedim_bin (vr5,7,micro_g(ngrid)%ffhl,micro_g(ngrid)%accph, &
                     micro_g(ngrid)%pcpvh,micro_g(ngrid)%pcprh, &
                     tempc,1,nkr)
   endif

endif

return
END SUBROUTINE precipitation

!##############################################################################
Subroutine sedim_bin (fallvelcm,lhcat,rx,precip,pcpv,pcpr,tempc,krs,kre)

use micro_prm
use mem_micro
use mem_grid
use node_mod
use mem_basic
use rconstants, only: alli
use micphys, only: igraup,ihail
implicit none

integer :: i,j,k,kk,kkf,kr,krs,kre,lhcat,maxkfall = 4
integer :: idensrtgt, ndensrtgt = 40
integer, save :: ialloc = 0, ncall = 0
real :: disp, ztopnew, zbotnew, fallin, rolddn0, newprecip
real :: densrtgt
real, dimension(nkr) :: fallvel, fallvelcm
real, dimension(mmxp(ngrid),mmyp(ngrid)) :: precip, pcpr, tempc
real, dimension(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)) :: pcpv
real, dimension(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),nkr) :: rx
real, allocatable, dimension(:,:,:,:) :: rnew
real, allocatable, dimension(:,:,:,:,:), save :: pcpfill
real, allocatable, dimension(:,:,:,:), save :: sfcpcp,allpcp
real, dimension(7) :: dpcp0
!For use in the leaf3 scheme
data dpcp0 /0.001, 0.010, 0.010, 0.010, 0.010, 0.003, 0.001/

!allocate 4D local array for "rnew"
allocate(rnew(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),nkr))

!If first call, make table
if (ncall == 0 .and. krs == 1) then
   fallvel=fallvelcm*0.01 !Convert velocicty from cm/s to m/s
   if (ialloc == 0) then
      allocate (pcpfill(mmzp(ngrid),maxkfall,nkr,ndensrtgt,7))
      allocate (sfcpcp(mmzp(ngrid),nkr,ndensrtgt,7))
      allocate (allpcp(mmzp(ngrid),nkr,7,ndensrtgt))
      pcpfill = 0.
      sfcpcp = 0.
      allpcp = 0.
      ialloc = 1
   endif
   ! Loop over vertical grid index.
   do k = 2,mmzp(ngrid)-1
   ! Loop over bins
      do kr = 1,nkr
         do idensrtgt=1,ndensrtgt
            densrtgt = idensrtgt / 10.0

            disp = dtlt * fallvel(kr)*densrtgt
            ztopnew = zm(k) - disp
            zbotnew = zm(k-1) - disp

            ! Loop over grid cells that a parcel falls into.
            do kkf = 1,min(k-1,maxkfall)
               !kk = k to k-3 (or 2)
               kk = k + 1 - kkf
               if (zbotnew .gt. zm(kk)) go to 50  !i.e. if precip does not fall into current box
               if (ztopnew .le. zm(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzt(kk) *  &
                     (min(zm(kk),ztopnew) - max(zm(kk-1),zbotnew))
               endif
               !Amount from level k that is now in level kk
               pcpfill(k,kkf,kr,idensrtgt,lhcat) = fallin
               
               !Precipitation rate (will divide by the time step later)
               allpcp(k,kr,lhcat,idensrtgt) = allpcp(k,kr,lhcat,idensrtgt) + disp 
            enddo
50       continue
      
            ! Compute surface precipitation.
            if (zbotnew .lt. 0.) sfcpcp(k,kr,idensrtgt,lhcat) = min(0.,ztopnew) - zbotnew
         enddo
      enddo
   enddo
   if( (iceprocs==0 .and. lhcat==1) .or. &
       (iceprocs==1 .and. ihail==0 .and. lhcat==6) .or. &
       (iceprocs==1 .and. ihail==1 .and. lhcat==7) ) then 
     ncall = 1
   endif
endif

!Now let precipitate
!kg/m2 = mm
rnew = 0.
do k = 2,mmzp(ngrid)
   do i = ia,iz
      do j = ja,jz
         do kr = krs,kre
            rolddn0 = rx(k,i,j,kr) * basic_g(ngrid)%dn0(k,i,j)
            idensrtgt = nint(10.0*(0.7/basic_g(ngrid)%dn0(k,i,j))**0.5/grid_g(ngrid)%rtgt(i,j))
            idensrtgt=max(1,min(40,idensrtgt))
            do kkf = 1,min(maxkfall,k-1)
               kk = k+1-kkf
               rnew(kk,i,j,kr) = rnew(kk,i,j,kr) + rolddn0 / &
                  basic_g(ngrid)%dn0(kk,i,j) * pcpfill(k,kkf,kr,idensrtgt,lhcat)               
            enddo
            newprecip = rolddn0 * col * sfcpcp(k,kr,idensrtgt,lhcat)

            !Precip is accumulated from the beginning of the simulation and is tracked 
            !separately for each species
            precip(i,j) = precip(i,j) + newprecip  
            pcpr(i,j) = pcpr(i,j) + newprecip / dtlt
            pcpv(k,i,j) = pcpv(k,i,j) + allpcp(k,kr,lhcat,idensrtgt) / dtlt * rolddn0

            !Pcpg is reset every timestep and is a total for all precipitation types
            micro_g(ngrid)%pcpg(i,j) = micro_g(ngrid)%pcpg(i,j) &
                + newprecip
            !dpcp0 converts mm to m
            micro_g(ngrid)%dpcpg(i,j) = micro_g(ngrid)%dpcpg(i,j) &
                + newprecip*dpcp0(lhcat)
            if (lhcat.eq.1) then
             micro_g(ngrid)%qpcpg(i,j) = micro_g(ngrid)%qpcpg(i,j) &
                + newprecip*(tempc(i,j)*4186.+alli)
            elseif (lhcat.gt.1) then
             micro_g(ngrid)%qpcpg(i,j) = micro_g(ngrid)%qpcpg(i,j) &
                + newprecip*(tempc(i,j)*2093.)
            endif
         enddo
      enddo
   enddo
enddo

rx = rnew

deallocate(rnew)

return
END SUBROUTINE sedim_bin

!##############################################################################
Subroutine update_thermo ()

use rconstants
use mem_micro
use mem_basic
use micro_prm
use mem_grid, only: ngrid
use node_mod
use micphys, only:ipris,igraup,ihail

implicit none

real, allocatable, dimension(:,:,:,:) :: ice_bins
real, allocatable, dimension(:,:,:) :: rliq,rice,tair,qhyd

integer:: k,i,j,izero

allocate(ice_bins(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),nkr))
allocate(rliq(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
allocate(rice(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
allocate(tair(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
allocate(qhyd(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

izero=0
CALL sum_bins (micro_g(ngrid)%ffcd,rliq,mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),1,nkr,izero)

if(iceprocs.eq.1) then
   ice_bins = micro_g(ngrid)%ffsn
   if (ipris == 1 .or. ipris >= 4) ice_bins = ice_bins + micro_g(ngrid)%ffic
   if (ipris == 2 .or. ipris >= 4) ice_bins = ice_bins + micro_g(ngrid)%ffip
   if (ipris >= 3) ice_bins = ice_bins + micro_g(ngrid)%ffid
   if (igraup > 0) ice_bins = ice_bins + micro_g(ngrid)%ffgl
   if (ihail > 0) ice_bins = ice_bins + micro_g(ngrid)%ffhl
   izero=0
   CALL sum_bins (ice_bins,rice,mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),1,nkr,izero)
else
   rice=0.
endif

do k=1,mmzp(ngrid)
   do i=ia,iz
      do j=ja,jz
         basic_g(ngrid)%rtp(k,i,j) = basic_g(ngrid)%rv(k,i,j) + rliq(k,i,j) + rice(k,i,j)

         tair(k,i,j) = basic_g(ngrid)%theta(k,i,j) * (basic_g(ngrid)%pi0(k,i,j) &
                       + basic_g(ngrid)%pp(k,i,j)) * cpi
         qhyd(k,i,j) = (alvl * rliq(k,i,j)+alvi * rice(k,i,j)) * cpi / (max(tair(k,i,j), 253.))
         basic_g(ngrid)%thp(k,i,j)= basic_g(ngrid)%theta(k,i,j) / (1. + qhyd(k,i,j))
      enddo
   enddo
enddo

deallocate(ice_bins,rliq,rice,tair,qhyd)

return
END SUBROUTINE update_thermo

!##############################################################################
Subroutine zero_budgets ()

! Start a new model step for the bin-microphysics budgets.
!
! Every budget rate array (<name>r) collects the amount its process adds
! during THIS step, so it is zeroed here at the start of every step.
! finish_bin_budgets() then adds that amount to the run total (acc<name>,
! only when IMBUDGETACC=1) and divides it by dtlt to give the rate per
! second. No output frequency is involved, so changing FRQSTATE or FRQLITE
! cannot change the budget values.

use mem_micro
use mem_grid, only:ngrid

implicit none

! Zero this step's amount of every allocated budget process
if (allocated(micro_g(ngrid)%latheatvapr)) micro_g(ngrid)%latheatvapr = 0.
if (allocated(micro_g(ngrid)%latheatfrzr)) micro_g(ngrid)%latheatfrzr = 0.
if (allocated(micro_g(ngrid)%nuccldrr)) micro_g(ngrid)%nuccldrr = 0.
if (allocated(micro_g(ngrid)%cld2rainr)) micro_g(ngrid)%cld2rainr = 0.
if (allocated(micro_g(ngrid)%ice2rainr)) micro_g(ngrid)%ice2rainr = 0.
if (allocated(micro_g(ngrid)%nucicerr)) micro_g(ngrid)%nucicerr = 0.
if (allocated(micro_g(ngrid)%vapliqr)) micro_g(ngrid)%vapliqr = 0.
if (allocated(micro_g(ngrid)%vapicer)) micro_g(ngrid)%vapicer = 0.
if (allocated(micro_g(ngrid)%evapliqr)) micro_g(ngrid)%evapliqr = 0.
if (allocated(micro_g(ngrid)%evapicer)) micro_g(ngrid)%evapicer = 0.
if (allocated(micro_g(ngrid)%freezingr)) micro_g(ngrid)%freezingr = 0.
if (allocated(micro_g(ngrid)%meltingr)) micro_g(ngrid)%meltingr = 0.
if (allocated(micro_g(ngrid)%melticer)) micro_g(ngrid)%melticer = 0.
if (allocated(micro_g(ngrid)%rimecldr)) micro_g(ngrid)%rimecldr = 0.
if (allocated(micro_g(ngrid)%rain2icer)) micro_g(ngrid)%rain2icer = 0.
if (allocated(micro_g(ngrid)%aggregater)) micro_g(ngrid)%aggregater = 0.
if (allocated(micro_g(ngrid)%inuchomrr)) micro_g(ngrid)%inuchomrr = 0.
if (allocated(micro_g(ngrid)%inuccontrr)) micro_g(ngrid)%inuccontrr = 0.
if (allocated(micro_g(ngrid)%inucifnrr)) micro_g(ngrid)%inucifnrr = 0.
if (allocated(micro_g(ngrid)%inuchazrr)) micro_g(ngrid)%inuchazrr = 0.
if (allocated(micro_g(ngrid)%vapcldr)) micro_g(ngrid)%vapcldr = 0.
if (allocated(micro_g(ngrid)%vaprainr)) micro_g(ngrid)%vaprainr = 0.
if (allocated(micro_g(ngrid)%vapprisr)) micro_g(ngrid)%vapprisr = 0.
if (allocated(micro_g(ngrid)%vapsnowr)) micro_g(ngrid)%vapsnowr = 0.
if (allocated(micro_g(ngrid)%vapaggrr)) micro_g(ngrid)%vapaggrr = 0.
if (allocated(micro_g(ngrid)%vapgraur)) micro_g(ngrid)%vapgraur = 0.
if (allocated(micro_g(ngrid)%vaphailr)) micro_g(ngrid)%vaphailr = 0.
if (allocated(micro_g(ngrid)%vapdrizr)) micro_g(ngrid)%vapdrizr = 0.
if (allocated(micro_g(ngrid)%evapcldr)) micro_g(ngrid)%evapcldr = 0.
if (allocated(micro_g(ngrid)%evaprainr)) micro_g(ngrid)%evaprainr = 0.
if (allocated(micro_g(ngrid)%evapprisr)) micro_g(ngrid)%evapprisr = 0.
if (allocated(micro_g(ngrid)%evapsnowr)) micro_g(ngrid)%evapsnowr = 0.
if (allocated(micro_g(ngrid)%evapaggrr)) micro_g(ngrid)%evapaggrr = 0.
if (allocated(micro_g(ngrid)%evapgraur)) micro_g(ngrid)%evapgraur = 0.
if (allocated(micro_g(ngrid)%evaphailr)) micro_g(ngrid)%evaphailr = 0.
if (allocated(micro_g(ngrid)%evapdrizr)) micro_g(ngrid)%evapdrizr = 0.
if (allocated(micro_g(ngrid)%meltprisr)) micro_g(ngrid)%meltprisr = 0.
if (allocated(micro_g(ngrid)%meltsnowr)) micro_g(ngrid)%meltsnowr = 0.
if (allocated(micro_g(ngrid)%meltaggrr)) micro_g(ngrid)%meltaggrr = 0.
if (allocated(micro_g(ngrid)%meltgraur)) micro_g(ngrid)%meltgraur = 0.
if (allocated(micro_g(ngrid)%melthailr)) micro_g(ngrid)%melthailr = 0.
if (allocated(micro_g(ngrid)%rimecldsnowr)) micro_g(ngrid)%rimecldsnowr = 0.
if (allocated(micro_g(ngrid)%rimecldaggrr)) micro_g(ngrid)%rimecldaggrr = 0.
if (allocated(micro_g(ngrid)%rimecldgraur)) micro_g(ngrid)%rimecldgraur = 0.
if (allocated(micro_g(ngrid)%rimecldhailr)) micro_g(ngrid)%rimecldhailr = 0.
if (allocated(micro_g(ngrid)%rain2prr)) micro_g(ngrid)%rain2prr = 0.
if (allocated(micro_g(ngrid)%rain2snr)) micro_g(ngrid)%rain2snr = 0.
if (allocated(micro_g(ngrid)%rain2agr)) micro_g(ngrid)%rain2agr = 0.
if (allocated(micro_g(ngrid)%rain2grr)) micro_g(ngrid)%rain2grr = 0.
if (allocated(micro_g(ngrid)%rain2har)) micro_g(ngrid)%rain2har = 0.
if (allocated(micro_g(ngrid)%aggrselfprisr)) micro_g(ngrid)%aggrselfprisr = 0.
if (allocated(micro_g(ngrid)%aggrselfsnowr)) micro_g(ngrid)%aggrselfsnowr = 0.
if (allocated(micro_g(ngrid)%aggrprissnowr)) micro_g(ngrid)%aggrprissnowr = 0.
if (allocated(micro_g(ngrid)%dust1cldrr)) micro_g(ngrid)%dust1cldrr = 0.
if (allocated(micro_g(ngrid)%dust2cldrr)) micro_g(ngrid)%dust2cldrr = 0.
if (allocated(micro_g(ngrid)%dust1drzrr)) micro_g(ngrid)%dust1drzrr = 0.
if (allocated(micro_g(ngrid)%dust2drzrr)) micro_g(ngrid)%dust2drzrr = 0.
if (allocated(micro_g(ngrid)%regen1cldrr)) micro_g(ngrid)%regen1cldrr = 0.
if (allocated(micro_g(ngrid)%regen2cldrr)) micro_g(ngrid)%regen2cldrr = 0.
if (allocated(micro_g(ngrid)%regen1drzrr)) micro_g(ngrid)%regen1drzrr = 0.
if (allocated(micro_g(ngrid)%regen2drzrr)) micro_g(ngrid)%regen2drzrr = 0.
if (allocated(micro_g(ngrid)%nuccldcr)) micro_g(ngrid)%nuccldcr = 0.
if (allocated(micro_g(ngrid)%nucicecr)) micro_g(ngrid)%nucicecr = 0.
if (allocated(micro_g(ngrid)%inuchomcr)) micro_g(ngrid)%inuchomcr = 0.
if (allocated(micro_g(ngrid)%inucifncr)) micro_g(ngrid)%inucifncr = 0.

return
END SUBROUTINE zero_budgets

!##############################################################################
Subroutine finish_bin_budgets ()

! Finish the model step for the bin-microphysics budgets.
!
! Each budget rate array (<name>r) now holds the amount its process added
! during this step. Add that amount to the never-reset run total
! (acc<name>, allocated only when IMBUDGETACC=1), then divide it by the
! step length dtlt so the output holds the instantaneous rate per second.

use mem_micro
use mem_grid, only:ngrid,dtlt
use micphys, only:imbudgetacc

implicit none

if (allocated(micro_g(ngrid)%latheatvapr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%acclatheatvap)) micro_g(ngrid)%acclatheatvap = micro_g(ngrid)%acclatheatvap + micro_g(ngrid)%latheatvapr
   micro_g(ngrid)%latheatvapr = micro_g(ngrid)%latheatvapr / dtlt
endif
if (allocated(micro_g(ngrid)%latheatfrzr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%acclatheatfrz)) micro_g(ngrid)%acclatheatfrz = micro_g(ngrid)%acclatheatfrz + micro_g(ngrid)%latheatfrzr
   micro_g(ngrid)%latheatfrzr = micro_g(ngrid)%latheatfrzr / dtlt
endif
if (allocated(micro_g(ngrid)%nuccldrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accnuccldr)) micro_g(ngrid)%accnuccldr = micro_g(ngrid)%accnuccldr + micro_g(ngrid)%nuccldrr
   micro_g(ngrid)%nuccldrr = micro_g(ngrid)%nuccldrr / dtlt
endif
if (allocated(micro_g(ngrid)%cld2rainr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%acccld2rain)) micro_g(ngrid)%acccld2rain = micro_g(ngrid)%acccld2rain + micro_g(ngrid)%cld2rainr
   micro_g(ngrid)%cld2rainr = micro_g(ngrid)%cld2rainr / dtlt
endif
if (allocated(micro_g(ngrid)%ice2rainr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accice2rain)) micro_g(ngrid)%accice2rain = micro_g(ngrid)%accice2rain + micro_g(ngrid)%ice2rainr
   micro_g(ngrid)%ice2rainr = micro_g(ngrid)%ice2rainr / dtlt
endif
if (allocated(micro_g(ngrid)%nucicerr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accnucicer)) micro_g(ngrid)%accnucicer = micro_g(ngrid)%accnucicer + micro_g(ngrid)%nucicerr
   micro_g(ngrid)%nucicerr = micro_g(ngrid)%nucicerr / dtlt
endif
if (allocated(micro_g(ngrid)%vapliqr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapliq)) micro_g(ngrid)%accvapliq = micro_g(ngrid)%accvapliq + micro_g(ngrid)%vapliqr
   micro_g(ngrid)%vapliqr = micro_g(ngrid)%vapliqr / dtlt
endif
if (allocated(micro_g(ngrid)%vapicer)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapice)) micro_g(ngrid)%accvapice = micro_g(ngrid)%accvapice + micro_g(ngrid)%vapicer
   micro_g(ngrid)%vapicer = micro_g(ngrid)%vapicer / dtlt
endif
if (allocated(micro_g(ngrid)%evapliqr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapliq)) micro_g(ngrid)%accevapliq = micro_g(ngrid)%accevapliq + micro_g(ngrid)%evapliqr
   micro_g(ngrid)%evapliqr = micro_g(ngrid)%evapliqr / dtlt
endif
if (allocated(micro_g(ngrid)%evapicer)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapice)) micro_g(ngrid)%accevapice = micro_g(ngrid)%accevapice + micro_g(ngrid)%evapicer
   micro_g(ngrid)%evapicer = micro_g(ngrid)%evapicer / dtlt
endif
if (allocated(micro_g(ngrid)%freezingr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accfreezing)) micro_g(ngrid)%accfreezing = micro_g(ngrid)%accfreezing + micro_g(ngrid)%freezingr
   micro_g(ngrid)%freezingr = micro_g(ngrid)%freezingr / dtlt
endif
if (allocated(micro_g(ngrid)%meltingr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmelting)) micro_g(ngrid)%accmelting = micro_g(ngrid)%accmelting + micro_g(ngrid)%meltingr
   micro_g(ngrid)%meltingr = micro_g(ngrid)%meltingr / dtlt
endif
if (allocated(micro_g(ngrid)%melticer)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmeltice)) micro_g(ngrid)%accmeltice = micro_g(ngrid)%accmeltice + micro_g(ngrid)%melticer
   micro_g(ngrid)%melticer = micro_g(ngrid)%melticer / dtlt
endif
if (allocated(micro_g(ngrid)%rimecldr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrimecld)) micro_g(ngrid)%accrimecld = micro_g(ngrid)%accrimecld + micro_g(ngrid)%rimecldr
   micro_g(ngrid)%rimecldr = micro_g(ngrid)%rimecldr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2icer)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2ice)) micro_g(ngrid)%accrain2ice = micro_g(ngrid)%accrain2ice + micro_g(ngrid)%rain2icer
   micro_g(ngrid)%rain2icer = micro_g(ngrid)%rain2icer / dtlt
endif
if (allocated(micro_g(ngrid)%aggregater)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accaggregate)) micro_g(ngrid)%accaggregate = micro_g(ngrid)%accaggregate + micro_g(ngrid)%aggregater
   micro_g(ngrid)%aggregater = micro_g(ngrid)%aggregater / dtlt
endif
if (allocated(micro_g(ngrid)%inuchomrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinuchomr)) micro_g(ngrid)%accinuchomr = micro_g(ngrid)%accinuchomr + micro_g(ngrid)%inuchomrr
   micro_g(ngrid)%inuchomrr = micro_g(ngrid)%inuchomrr / dtlt
endif
if (allocated(micro_g(ngrid)%inuccontrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinuccontr)) micro_g(ngrid)%accinuccontr = micro_g(ngrid)%accinuccontr + micro_g(ngrid)%inuccontrr
   micro_g(ngrid)%inuccontrr = micro_g(ngrid)%inuccontrr / dtlt
endif
if (allocated(micro_g(ngrid)%inucifnrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinucifnr)) micro_g(ngrid)%accinucifnr = micro_g(ngrid)%accinucifnr + micro_g(ngrid)%inucifnrr
   micro_g(ngrid)%inucifnrr = micro_g(ngrid)%inucifnrr / dtlt
endif
if (allocated(micro_g(ngrid)%inuchazrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinuchazr)) micro_g(ngrid)%accinuchazr = micro_g(ngrid)%accinuchazr + micro_g(ngrid)%inuchazrr
   micro_g(ngrid)%inuchazrr = micro_g(ngrid)%inuchazrr / dtlt
endif
if (allocated(micro_g(ngrid)%vapcldr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapcld)) micro_g(ngrid)%accvapcld = micro_g(ngrid)%accvapcld + micro_g(ngrid)%vapcldr
   micro_g(ngrid)%vapcldr = micro_g(ngrid)%vapcldr / dtlt
endif
if (allocated(micro_g(ngrid)%vaprainr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvaprain)) micro_g(ngrid)%accvaprain = micro_g(ngrid)%accvaprain + micro_g(ngrid)%vaprainr
   micro_g(ngrid)%vaprainr = micro_g(ngrid)%vaprainr / dtlt
endif
if (allocated(micro_g(ngrid)%vapprisr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvappris)) micro_g(ngrid)%accvappris = micro_g(ngrid)%accvappris + micro_g(ngrid)%vapprisr
   micro_g(ngrid)%vapprisr = micro_g(ngrid)%vapprisr / dtlt
endif
if (allocated(micro_g(ngrid)%vapsnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapsnow)) micro_g(ngrid)%accvapsnow = micro_g(ngrid)%accvapsnow + micro_g(ngrid)%vapsnowr
   micro_g(ngrid)%vapsnowr = micro_g(ngrid)%vapsnowr / dtlt
endif
if (allocated(micro_g(ngrid)%vapaggrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapaggr)) micro_g(ngrid)%accvapaggr = micro_g(ngrid)%accvapaggr + micro_g(ngrid)%vapaggrr
   micro_g(ngrid)%vapaggrr = micro_g(ngrid)%vapaggrr / dtlt
endif
if (allocated(micro_g(ngrid)%vapgraur)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapgrau)) micro_g(ngrid)%accvapgrau = micro_g(ngrid)%accvapgrau + micro_g(ngrid)%vapgraur
   micro_g(ngrid)%vapgraur = micro_g(ngrid)%vapgraur / dtlt
endif
if (allocated(micro_g(ngrid)%vaphailr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvaphail)) micro_g(ngrid)%accvaphail = micro_g(ngrid)%accvaphail + micro_g(ngrid)%vaphailr
   micro_g(ngrid)%vaphailr = micro_g(ngrid)%vaphailr / dtlt
endif
if (allocated(micro_g(ngrid)%vapdrizr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accvapdriz)) micro_g(ngrid)%accvapdriz = micro_g(ngrid)%accvapdriz + micro_g(ngrid)%vapdrizr
   micro_g(ngrid)%vapdrizr = micro_g(ngrid)%vapdrizr / dtlt
endif
if (allocated(micro_g(ngrid)%evapcldr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapcld)) micro_g(ngrid)%accevapcld = micro_g(ngrid)%accevapcld + micro_g(ngrid)%evapcldr
   micro_g(ngrid)%evapcldr = micro_g(ngrid)%evapcldr / dtlt
endif
if (allocated(micro_g(ngrid)%evaprainr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevaprain)) micro_g(ngrid)%accevaprain = micro_g(ngrid)%accevaprain + micro_g(ngrid)%evaprainr
   micro_g(ngrid)%evaprainr = micro_g(ngrid)%evaprainr / dtlt
endif
if (allocated(micro_g(ngrid)%evapprisr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevappris)) micro_g(ngrid)%accevappris = micro_g(ngrid)%accevappris + micro_g(ngrid)%evapprisr
   micro_g(ngrid)%evapprisr = micro_g(ngrid)%evapprisr / dtlt
endif
if (allocated(micro_g(ngrid)%evapsnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapsnow)) micro_g(ngrid)%accevapsnow = micro_g(ngrid)%accevapsnow + micro_g(ngrid)%evapsnowr
   micro_g(ngrid)%evapsnowr = micro_g(ngrid)%evapsnowr / dtlt
endif
if (allocated(micro_g(ngrid)%evapaggrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapaggr)) micro_g(ngrid)%accevapaggr = micro_g(ngrid)%accevapaggr + micro_g(ngrid)%evapaggrr
   micro_g(ngrid)%evapaggrr = micro_g(ngrid)%evapaggrr / dtlt
endif
if (allocated(micro_g(ngrid)%evapgraur)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapgrau)) micro_g(ngrid)%accevapgrau = micro_g(ngrid)%accevapgrau + micro_g(ngrid)%evapgraur
   micro_g(ngrid)%evapgraur = micro_g(ngrid)%evapgraur / dtlt
endif
if (allocated(micro_g(ngrid)%evaphailr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevaphail)) micro_g(ngrid)%accevaphail = micro_g(ngrid)%accevaphail + micro_g(ngrid)%evaphailr
   micro_g(ngrid)%evaphailr = micro_g(ngrid)%evaphailr / dtlt
endif
if (allocated(micro_g(ngrid)%evapdrizr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accevapdriz)) micro_g(ngrid)%accevapdriz = micro_g(ngrid)%accevapdriz + micro_g(ngrid)%evapdrizr
   micro_g(ngrid)%evapdrizr = micro_g(ngrid)%evapdrizr / dtlt
endif
if (allocated(micro_g(ngrid)%meltprisr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmeltpris)) micro_g(ngrid)%accmeltpris = micro_g(ngrid)%accmeltpris + micro_g(ngrid)%meltprisr
   micro_g(ngrid)%meltprisr = micro_g(ngrid)%meltprisr / dtlt
endif
if (allocated(micro_g(ngrid)%meltsnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmeltsnow)) micro_g(ngrid)%accmeltsnow = micro_g(ngrid)%accmeltsnow + micro_g(ngrid)%meltsnowr
   micro_g(ngrid)%meltsnowr = micro_g(ngrid)%meltsnowr / dtlt
endif
if (allocated(micro_g(ngrid)%meltaggrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmeltaggr)) micro_g(ngrid)%accmeltaggr = micro_g(ngrid)%accmeltaggr + micro_g(ngrid)%meltaggrr
   micro_g(ngrid)%meltaggrr = micro_g(ngrid)%meltaggrr / dtlt
endif
if (allocated(micro_g(ngrid)%meltgraur)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmeltgrau)) micro_g(ngrid)%accmeltgrau = micro_g(ngrid)%accmeltgrau + micro_g(ngrid)%meltgraur
   micro_g(ngrid)%meltgraur = micro_g(ngrid)%meltgraur / dtlt
endif
if (allocated(micro_g(ngrid)%melthailr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accmelthail)) micro_g(ngrid)%accmelthail = micro_g(ngrid)%accmelthail + micro_g(ngrid)%melthailr
   micro_g(ngrid)%melthailr = micro_g(ngrid)%melthailr / dtlt
endif
if (allocated(micro_g(ngrid)%rimecldsnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrimecldsnow)) micro_g(ngrid)%accrimecldsnow = micro_g(ngrid)%accrimecldsnow + micro_g(ngrid)%rimecldsnowr
   micro_g(ngrid)%rimecldsnowr = micro_g(ngrid)%rimecldsnowr / dtlt
endif
if (allocated(micro_g(ngrid)%rimecldaggrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrimecldaggr)) micro_g(ngrid)%accrimecldaggr = micro_g(ngrid)%accrimecldaggr + micro_g(ngrid)%rimecldaggrr
   micro_g(ngrid)%rimecldaggrr = micro_g(ngrid)%rimecldaggrr / dtlt
endif
if (allocated(micro_g(ngrid)%rimecldgraur)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrimecldgrau)) micro_g(ngrid)%accrimecldgrau = micro_g(ngrid)%accrimecldgrau + micro_g(ngrid)%rimecldgraur
   micro_g(ngrid)%rimecldgraur = micro_g(ngrid)%rimecldgraur / dtlt
endif
if (allocated(micro_g(ngrid)%rimecldhailr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrimecldhail)) micro_g(ngrid)%accrimecldhail = micro_g(ngrid)%accrimecldhail + micro_g(ngrid)%rimecldhailr
   micro_g(ngrid)%rimecldhailr = micro_g(ngrid)%rimecldhailr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2prr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2pr)) micro_g(ngrid)%accrain2pr = micro_g(ngrid)%accrain2pr + micro_g(ngrid)%rain2prr
   micro_g(ngrid)%rain2prr = micro_g(ngrid)%rain2prr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2snr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2sn)) micro_g(ngrid)%accrain2sn = micro_g(ngrid)%accrain2sn + micro_g(ngrid)%rain2snr
   micro_g(ngrid)%rain2snr = micro_g(ngrid)%rain2snr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2agr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2ag)) micro_g(ngrid)%accrain2ag = micro_g(ngrid)%accrain2ag + micro_g(ngrid)%rain2agr
   micro_g(ngrid)%rain2agr = micro_g(ngrid)%rain2agr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2grr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2gr)) micro_g(ngrid)%accrain2gr = micro_g(ngrid)%accrain2gr + micro_g(ngrid)%rain2grr
   micro_g(ngrid)%rain2grr = micro_g(ngrid)%rain2grr / dtlt
endif
if (allocated(micro_g(ngrid)%rain2har)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accrain2ha)) micro_g(ngrid)%accrain2ha = micro_g(ngrid)%accrain2ha + micro_g(ngrid)%rain2har
   micro_g(ngrid)%rain2har = micro_g(ngrid)%rain2har / dtlt
endif
if (allocated(micro_g(ngrid)%aggrselfprisr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accaggrselfpris)) micro_g(ngrid)%accaggrselfpris = micro_g(ngrid)%accaggrselfpris + micro_g(ngrid)%aggrselfprisr
   micro_g(ngrid)%aggrselfprisr = micro_g(ngrid)%aggrselfprisr / dtlt
endif
if (allocated(micro_g(ngrid)%aggrselfsnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accaggrselfsnow)) micro_g(ngrid)%accaggrselfsnow = micro_g(ngrid)%accaggrselfsnow + micro_g(ngrid)%aggrselfsnowr
   micro_g(ngrid)%aggrselfsnowr = micro_g(ngrid)%aggrselfsnowr / dtlt
endif
if (allocated(micro_g(ngrid)%aggrprissnowr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accaggrprissnow)) micro_g(ngrid)%accaggrprissnow = micro_g(ngrid)%accaggrprissnow + micro_g(ngrid)%aggrprissnowr
   micro_g(ngrid)%aggrprissnowr = micro_g(ngrid)%aggrprissnowr / dtlt
endif
if (allocated(micro_g(ngrid)%dust1cldrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accdust1cldr)) micro_g(ngrid)%accdust1cldr = micro_g(ngrid)%accdust1cldr + micro_g(ngrid)%dust1cldrr
   micro_g(ngrid)%dust1cldrr = micro_g(ngrid)%dust1cldrr / dtlt
endif
if (allocated(micro_g(ngrid)%dust2cldrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accdust2cldr)) micro_g(ngrid)%accdust2cldr = micro_g(ngrid)%accdust2cldr + micro_g(ngrid)%dust2cldrr
   micro_g(ngrid)%dust2cldrr = micro_g(ngrid)%dust2cldrr / dtlt
endif
if (allocated(micro_g(ngrid)%dust1drzrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accdust1drzr)) micro_g(ngrid)%accdust1drzr = micro_g(ngrid)%accdust1drzr + micro_g(ngrid)%dust1drzrr
   micro_g(ngrid)%dust1drzrr = micro_g(ngrid)%dust1drzrr / dtlt
endif
if (allocated(micro_g(ngrid)%dust2drzrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accdust2drzr)) micro_g(ngrid)%accdust2drzr = micro_g(ngrid)%accdust2drzr + micro_g(ngrid)%dust2drzrr
   micro_g(ngrid)%dust2drzrr = micro_g(ngrid)%dust2drzrr / dtlt
endif
if (allocated(micro_g(ngrid)%regen1cldrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accregen1cldr)) micro_g(ngrid)%accregen1cldr = micro_g(ngrid)%accregen1cldr + micro_g(ngrid)%regen1cldrr
   micro_g(ngrid)%regen1cldrr = micro_g(ngrid)%regen1cldrr / dtlt
endif
if (allocated(micro_g(ngrid)%regen2cldrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accregen2cldr)) micro_g(ngrid)%accregen2cldr = micro_g(ngrid)%accregen2cldr + micro_g(ngrid)%regen2cldrr
   micro_g(ngrid)%regen2cldrr = micro_g(ngrid)%regen2cldrr / dtlt
endif
if (allocated(micro_g(ngrid)%regen1drzrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accregen1drzr)) micro_g(ngrid)%accregen1drzr = micro_g(ngrid)%accregen1drzr + micro_g(ngrid)%regen1drzrr
   micro_g(ngrid)%regen1drzrr = micro_g(ngrid)%regen1drzrr / dtlt
endif
if (allocated(micro_g(ngrid)%regen2drzrr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accregen2drzr)) micro_g(ngrid)%accregen2drzr = micro_g(ngrid)%accregen2drzr + micro_g(ngrid)%regen2drzrr
   micro_g(ngrid)%regen2drzrr = micro_g(ngrid)%regen2drzrr / dtlt
endif
if (allocated(micro_g(ngrid)%nuccldcr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accnuccldc)) micro_g(ngrid)%accnuccldc = micro_g(ngrid)%accnuccldc + micro_g(ngrid)%nuccldcr
   micro_g(ngrid)%nuccldcr = micro_g(ngrid)%nuccldcr / dtlt
endif
if (allocated(micro_g(ngrid)%nucicecr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accnucicec)) micro_g(ngrid)%accnucicec = micro_g(ngrid)%accnucicec + micro_g(ngrid)%nucicecr
   micro_g(ngrid)%nucicecr = micro_g(ngrid)%nucicecr / dtlt
endif
if (allocated(micro_g(ngrid)%inuchomcr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinuchomc)) micro_g(ngrid)%accinuchomc = micro_g(ngrid)%accinuchomc + micro_g(ngrid)%inuchomcr
   micro_g(ngrid)%inuchomcr = micro_g(ngrid)%inuchomcr / dtlt
endif
if (allocated(micro_g(ngrid)%inucifncr)) then
   if (imbudgetacc == 1 .and. allocated(micro_g(ngrid)%accinucifnc)) micro_g(ngrid)%accinucifnc = micro_g(ngrid)%accinucifnc + micro_g(ngrid)%inucifncr
   micro_g(ngrid)%inucifncr = micro_g(ngrid)%inucifncr / dtlt
endif

return
END SUBROUTINE finish_bin_budgets

!##############################################################################
Subroutine adj1_bin (m1,m2,m3,rtp,micro)

use mem_micro
use micro_prm, only:nkr,col,iceprocs,iceflag
use micphys, only:ipris,igraup,ihail

implicit none

type (micro_vars) :: micro

integer :: m1,m2,m3
integer :: i,j,k,kr
real :: frac, totbef
real, dimension(m1,m2,m3) :: rtp
real, dimension(m1) :: vctr9

do j=1,m3
   do i=1,m2
      do k=1,m1
         totbef=0.
         do kr = 1,nkr
            totbef = totbef + micro%ffcd(k,i,j,kr)
            if (iceprocs == 1) then  
               totbef = totbef + micro%ffsn(k,i,j,kr)
               if (ipris == 1 .or. ipris >= 4) totbef = totbef + micro%ffic(k,i,j,kr)
               if (ipris == 2 .or. ipris >= 4) totbef = totbef + micro%ffip(k,i,j,kr)
               if (ipris >= 3) totbef = totbef + micro%ffid(k,i,j,kr)
               if (igraup > 0) totbef = totbef + micro%ffgl(k,i,j,kr)
               if (ihail > 0) totbef = totbef + micro%ffhl(k,i,j,kr)
            endif
         enddo
         do kr = 1,nkr
            if (micro%ffcd(k,i,j,kr)*col < 1.e-12) micro%ffcd(k,i,j,kr) = 0.
            if (iceprocs == 1) then
               if (ipris == 1 .or. ipris >= 4) then
                  if (micro%ffic(k,i,j,kr)*col < 1.e-12) micro%ffic(k,i,j,kr) = 0.
               endif
               if (ipris == 2 .or. ipris >= 4) then
                  if (micro%ffip(k,i,j,kr)*col < 1.e-12) micro%ffip(k,i,j,kr) = 0.
               endif
               if (ipris >= 3) then
                  if (micro%ffid(k,i,j,kr)*col < 1.e-12) micro%ffid(k,i,j,kr) = 0.
               endif
               if (micro%ffsn(k,i,j,kr)*col < 1.e-12) micro%ffsn(k,i,j,kr) = 0.
               if (igraup > 0) then
                  if (micro%ffgl(k,i,j,kr)*col < 1.e-12) micro%ffgl(k,i,j,kr) = 0.
               endif
               if (ihail > 0) then
                  if (micro%ffhl(k,i,j,kr)*col < 1.e-12) micro%ffhl(k,i,j,kr) = 0.
               endif
               if (iceflag == 1) then
                  if (micro%ffin(k,i,j,kr)*col < 1.e-12) micro%ffin(k,i,j,kr) = 0.
               endif
             endif
            if (micro%fncn(k,i,j,kr)*col < 1.e-12) micro%fncn(k,i,j,kr) = 0.
         enddo

         if (totbef < 0.) then
            totbef = 0.0
         endif   
         rtp(k,i,j) = max(totbef * col * 1.01,rtp(k,i,j))
         vctr9(k)=0.
         do kr = 1,nkr
            vctr9(k) = vctr9(k) + micro%ffcd(k,i,j,kr) 
            if (iceprocs == 1) then  
               vctr9(k) = vctr9(k) + micro%ffsn(k,i,j,kr)
               if (ipris == 1 .or. ipris >= 4) vctr9(k) = vctr9(k) + micro%ffic(k,i,j,kr)
               if (ipris == 2 .or. ipris >= 4) vctr9(k) = vctr9(k) + micro%ffip(k,i,j,kr)
               if (ipris >= 3) vctr9(k) = vctr9(k) + micro%ffid(k,i,j,kr)
               if (igraup > 0) vctr9(k) = vctr9(k) + micro%ffgl(k,i,j,kr)
               if (ihail > 0) vctr9(k) = vctr9(k) + micro%ffhl(k,i,j,kr)
            endif
         enddo

         if (vctr9(k) > totbef .and. vctr9(k) > 0.) then
            frac = totbef / (vctr9(k))
            do kr = 1,nkr
              micro%ffcd(k,i,j,kr) = micro%ffcd(k,i,j,kr) * frac
              if (iceprocs == 1) then
               if (ipris==1 .or. ipris>=4) micro%ffic(k,i,j,kr) = micro%ffic(k,i,j,kr)*frac
               if (ipris==2 .or. ipris>=4) micro%ffip(k,i,j,kr) = micro%ffip(k,i,j,kr)*frac
               if (ipris>=3) micro%ffid(k,i,j,kr) = micro%ffid(k,i,j,kr) * frac
               micro%ffsn(k,i,j,kr) = micro%ffsn(k,i,j,kr) * frac
               if (igraup > 0) micro%ffgl(k,i,j,kr) = micro%ffgl(k,i,j,kr) * frac
               if (ihail > 0) micro%ffhl(k,i,j,kr) = micro%ffhl(k,i,j,kr) * frac
              endif
            enddo
         endif
      enddo
   enddo
enddo

return
END SUBROUTINE adj1_bin

!##############################################################################
Subroutine fastsbm_reassign (FF2R,FF3R,FF4R,FF5R)

use micphys, only:ipris,igraup,ihail
use micro_prm, only:nkr,icemax

implicit none

real, dimension(nkr,icemax) :: FF2R
real, dimension(nkr) :: FF3R,FF4R,FF5R
integer :: ikl

IF (IPRIS == 1 .or. IPRIS == 2 .or. IPRIS == 3) THEN
   FF2R(:,IPRIS) = FF2R(:,1)+FF2R(:,2)+FF2R(:,3)
   DO IKL=1,ICEMAX
      IF (IKL .ne. IPRIS) FF2R(:,IKL)=0.
   ENDDO
ELSEIF (IPRIS == 0) THEN
   FF3R(:) = FF3R(:)+FF2R(:,1)+FF2R(:,2)+FF2R(:,3)
   DO IKL=1,ICEMAX
      FF2R(:,IKL)=0.
   ENDDO
ENDIF
IF (IGRAUP == 0) THEN
   FF5R(:)=FF5R(:)+FF4R(:)
   FF4R(:)=0.
ENDIF
IF (IHAIL == 0) THEN
   FF4R(:)=FF4R(:)+FF5R(:)
   FF5R(:)=0.
ENDIF

return
END SUBROUTINE fastsbm_reassign

!##############################################################################
Subroutine vap_budget (ff1r,ff2r,ff3r,ff4r,ff5r,dens,k,i,j,str)

use mem_micro
use micro_prm, only:nkr,iceprocs,icemax,krdrop,krpris
use module_hujisbm, only:xl,xi,xs,xg,xh,sum_mass
use micphys, only:imbudget, ipris, igraup, ihail
use mem_grid, only:ngrid

implicit none

real, dimension(nkr) :: ff1r,ff3r,ff4r,ff5r
real, dimension(nkr,icemax) :: ff2r
real :: dens,plusminus,sum_pris,sum_snow
integer :: k,i,j
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (imbudget >= 1) then
   micro_g(ngrid)%vapliqr(k,i,j) = micro_g(ngrid)%vapliqr(k,i,j) + &
      plusminus * sum_mass(ff1r,xl(:),dens,1,nkr)

   if(iceprocs.eq.1) &
      micro_g(ngrid)%vapicer(k,i,j) = micro_g(ngrid)%vapicer(k,i,j) + plusminus * &
          sum_mass(ff2r(:,1)+ff2r(:,2)+ff2r(:,3)+ff3r+ff4r+ff5r,xs,dens,1,nkr)
endif
if (imbudget >= 2) then
   ! micro_g(ngrid)%vapcldr(k,i,j) = micro_g(ngrid)%vapcldr(k,i,j) +  &
   !              plusminus * sum_mass(ff1r,xl,dens,1,krdrop-1)
   ! micro_g(ngrid)%vaprainr(k,i,j) = micro_g(ngrid)%vaprainr(k,i,j) + &
   !              plusminus * sum_mass(ff1r,xl,dens,krdrop,nkr)
    micro_g(ngrid)%cld2rainr(k,i,j) = micro_g(ngrid)%cld2rainr(k,i,j) + &
                  plusminus * sum_mass(ff1r,xl,dens,krdrop,nkr) + &
                  plusminus * -1. * micro_g(ngrid)%vaprainr(k,i,j)

    if(iceprocs.eq.1) then
       if (ipris > 0) &
       micro_g(ngrid)%vapprisr(k,i,j) = micro_g(ngrid)%vapprisr(k,i,j) + &
                      plusminus * sum_pris(ff2r,dens)
       if (ipris > 0) &
       micro_g(ngrid)%vapsnowr(k,i,j) = micro_g(ngrid)%vapsnowr(k,i,j) + &
                      plusminus * sum_snow(ff2r,dens)
       micro_g(ngrid)%vapaggrr(k,i,j) = micro_g(ngrid)%vapaggrr(k,i,j) + &
                      plusminus * sum_mass(ff3r,xs,dens,1,nkr)
       if (igraup > 0) &
       micro_g(ngrid)%vapgraur(k,i,j) = micro_g(ngrid)%vapgraur(k,i,j) + &
                      plusminus * sum_mass(ff4r,xg,dens,1,nkr)
       if (ihail > 0) &
       micro_g(ngrid)%vaphailr(k,i,j) = micro_g(ngrid)%vaphailr(k,i,j) + &
                      plusminus * sum_mass(ff5r,xh,dens,1,nkr)
    endif
endif

return
END SUBROUTINE vap_budget

!##############################################################################
Subroutine lhv_budget (ff1r,ff2r,ff3r,ff4r,ff5r,pitot,dens,k,i,j,str)

use mem_micro
use micro_prm, only:nkr,iceprocs,icemax
use module_hujisbm, only:xl,xi,xs,xg,xh,sum_mass
use micphys, only:imbudget,lhrtheta
use mem_grid, only:ngrid
use rconstants, only:cpi,alvl,alvi

implicit none

real, dimension(nkr) :: ff1r,ff3r,ff4r,ff5r
real, dimension(nkr,icemax) :: ff2r
real :: pitot,dens,fac,plusminus
integer :: k,i,j
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (lhrtheta) then
   fac = 1./pitot
else
   fac = cpi
endif

if (imbudget >= 1) then
   micro_g(ngrid)%latheatvapr(k,i,j) = micro_g(ngrid)%latheatvapr(k,i,j) + &
      plusminus * alvl * fac * sum_mass(ff1r,xl(:),dens,1,nkr) 

   if(iceprocs.eq.1) &
      micro_g(ngrid)%latheatvapr(k,i,j) = micro_g(ngrid)%latheatvapr(k,i,j) + &
         plusminus * alvi * fac * &
         sum_mass(ff2r(:,1)+ff2r(:,2)+ff2r(:,3)+ff3r+ff4r+ff5r,xs,dens,1,nkr) 
endif


return
END SUBROUTINE lhv_budget

!##############################################################################
Subroutine lhf_budget (ff1r,pitot,dens,k,i,j,str)

use mem_micro
use micro_prm, only:nkr,iceprocs
use module_hujisbm, only:xl,sum_mass
use micphys, only:imbudget,lhrtheta
use mem_grid, only:ngrid
use rconstants, only:cpi,alli

implicit none

real, dimension(nkr) :: ff1r
real :: pitot,dens,fac,plusminus
integer :: k,i,j
character(len=3) :: str

!plusminus = -1 for budget prep, = +1 for budget finish
if (str .eq. 'beg') then
   plusminus = -1.
elseif (str .eq. 'end') then
   plusminus = 1.
else
   print*, 'Invalid string for budgets'
   stop
endif

if (lhrtheta) then
   fac = 1./pitot
else
   fac = cpi
endif

if (imbudget >= 1 .and. iceprocs == 1) then
   micro_g(ngrid)%latheatfrzr(k,i,j) = micro_g(ngrid)%latheatfrzr(k,i,j) + &
      plusminus * (-1.) * alli * fac * sum_mass(ff1r,xl(:),dens,1,nkr) 
endif


return
END SUBROUTINE lhf_budget
