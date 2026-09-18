!##############################################################################
Module mem_micro

implicit none

   Type micro_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          rcp,rdp,rrp,rpp,rsp,rap,rgp,rhp           &
                         ,ccp,cdp,crp,cpp,csp,cap,cgp,chp           &
                         ,cifnp,q2,q6,q7                            &
                         ,pcpvd,pcpvr,pcpvp,pcpvs,pcpva,pcpvg,pcpvh &
                         !Bin precip vars
                         ,pcpvic,pcpvip,pcpvid                      &
                         !Aerosol categories mass and number
                         ,cn1np,cn2np,cn1mp,cn2mp                   &
                         ,md1np,md2np,md1mp,md2mp                   &
                         ,salt_film_np,salt_jet_np,salt_spum_np     &
                         ,salt_film_mp,salt_jet_mp,salt_spum_mp     &
                         ,abc1np,abc2np,abc1mp,abc2mp               &
                         ,regen_aero1_np,regen_aero1_mp             &
                         ,regen_aero2_np,regen_aero2_mp             &
                         ,regen_aero1_nr,regen_aero1_mr             &
                         ,regen_aero2_nr,regen_aero2_mr             &
                         ,accregen_aero1_n,accregen_aero1_m         &
                         ,accregen_aero2_n,accregen_aero2_m         &
                         !Immersion freezing nuclei tracking
                         ,immercp,immerdp,immerrp,ifnnucp           &
                         !Aerosol tracking variables
                         ,cnmcp,cnmdp,cnmrp,cnmpp,cnmsp             &
                         ,cnmap,cnmgp,cnmhp                         &
                         ,dnmcp,dnmdp,dnmrp,dnmpp,dnmsp             &
                         ,dnmap,dnmgp,dnmhp                         &
                         ,dincp,dindp,dinrp,dinpp,dinsp             &
                         ,dinap,dingp,dinhp                         &
                         !Regenerated aerosol tracking variables
                         ,rnmcp,rnmdp,rnmrp,rnmpp,rnmsp             &
                         ,rnmap,rnmgp,rnmhp                         &
                         ,rincp,rindp,rinrp,rinpp,rinsp             &
                         ,rinap,ringp,rinhp                         &
                         ,snmcp,snmdp,snmrp,snmpp,snmsp             &
                         ,snmap,snmgp,snmhp                         &
                         ,resol_aero1_mp,resol_aero2_mp &
           ! MICRO BUDGET PROCESS RATES, per second (imbudget>=1..3)
           ,latheatvapr,latheatfrzr,nuccldrr,cld2rainr,ice2rainr,nucicerr &
           ,vapliqr,vapicer,evapliqr,evapicer,freezingr,meltingr,melticer &
           ,rimecldr,rain2icer,aggregater,inuchomrr,inuccontrr,inucifnrr,inuchazrr &
           ,vapcldr,vaprainr,vapprisr,vapsnowr,vapaggrr,vapgraur,vaphailr &
           ,vapdrizr,evapcldr,evaprainr,evapprisr,evapsnowr,evapaggrr,evapgraur &
           ,evaphailr,evapdrizr,meltprisr,meltsnowr,meltaggrr,meltgraur,melthailr &
           ,rimecldsnowr,rimecldaggrr,rimecldgraur,rimecldhailr,rain2prr,rain2snr &
           ,rain2agr,rain2grr,rain2har,aggrselfprisr,aggrselfsnowr,aggrprissnowr &
           ,dust1cldrr,dust2cldrr,dust1drzrr,dust2drzrr,regen1cldrr,regen2cldrr &
           ,regen1drzrr,regen2drzrr &
           ! MICRO BUDGET RUN TOTALS (imbudget>=1..3 and imbudgetacc==1)
           ,acclatheatvap,acclatheatfrz,accnuccldr,acccld2rain,accice2rain &
           ,accnucicer,accvapliq,accvapice,accevapliq,accevapice,accfreezing &
           ,accmelting,accmeltice,accrimecld,accrain2ice,accaggregate,accinuchomr &
           ,accinuccontr,accinucifnr,accinuchazr,accvapcld,accvaprain,accvappris &
           ,accvapsnow,accvapaggr,accvapgrau,accvaphail,accvapdriz,accevapcld &
           ,accevaprain,accevappris,accevapsnow,accevapaggr,accevapgrau,accevaphail &
           ,accevapdriz,accmeltpris,accmeltsnow,accmeltaggr,accmeltgrau,accmelthail &
           ,accrimecldsnow,accrimecldaggr,accrimecldgrau,accrimecldhail,accrain2pr &
           ,accrain2sn,accrain2ag,accrain2gr,accrain2ha,accaggrselfpris,accaggrselfsnow &
           ,accaggrprissnow,accdust1cldr,accdust2cldr,accdust1drzr,accdust2drzr &
           ,accregen1cldr,accregen2cldr,accregen1drzr,accregen2drzr &
           ! BIN MICROPHYSICS, and its extra number budgets
           ,t_old,rv_old &
           ,nuccldcr,nucicecr,inuchomcr,inucifncr &
           ,accnuccldc,accnucicec,accinuchomc,accinucifnc

   ! Variables to be dimenstioned by (nzp,nxp,nyp,nkr) for bin microphysics
   real, allocatable, dimension(:,:,:,:) :: &
            fncn,ffcd,ffic,ffip,ffid,ffsn,ffgl,ffhl,ffin

   ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          accpr,accpp,accps,accpa,accpg,accph,accpd &
                         ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh,pcprd &
                         ,pcpg,qpcpg,dpcpg                          &
                         !Accumulated aerosols and accumulation rate
                         ,accpdust,pcprdust,accpaero,pcpraero       &
                         ,accpregen,pcprregen                       &
                         !Bin precip vars
                         ,accpic,accpip,accpid,pcpric,pcprid,pcprip &
                         !Dust erodible fraction
                         ,dustfrac
                          
   End Type               
                          
   type (micro_vars), allocatable :: micro_g(:), microm_g(:)

   !Sedimentation table variables
   Type pcp_tab_type
   real, allocatable, dimension(:,:,:,:,:,:) :: pcpfillc,pcpfillr
   real, allocatable, dimension(:,:,:,:,:)   :: sfcpcp
   real, allocatable, dimension(:,:,:,:,:)   :: allpcp
   End Type
   type (pcp_tab_type), allocatable :: pcp_tab(:)

Contains                  

!##############################################################################
Subroutine alloc_sedim (pcp_tab,n1)

use micphys

implicit none

   type (pcp_tab_type) :: pcp_tab
   integer, intent(in) :: n1

! Micro Level=3 sedimentation tables
   if (level == 3) then
    allocate (pcp_tab%pcpfillc(n1,maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%pcpfillr(n1,maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%sfcpcp(maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%allpcp(n1,nembfall,nhcat,ndensrtgt,nband))
   endif

return
END SUBROUTINE alloc_sedim

!##############################################################################
Subroutine dealloc_sedim (pcp_tab)

implicit none

   type (pcp_tab_type) :: pcp_tab

   if (allocated(pcp_tab%pcpfillc)) deallocate(pcp_tab%pcpfillc)
   if (allocated(pcp_tab%pcpfillr)) deallocate(pcp_tab%pcpfillr)
   if (allocated(pcp_tab%sfcpcp))   deallocate(pcp_tab%sfcpcp)
   if (allocated(pcp_tab%allpcp))   deallocate(pcp_tab%allpcp)

return
END SUBROUTINE dealloc_sedim

!##############################################################################
Subroutine alloc_micro (micro,n1,n2,n3,n4)

use micphys
use micro_prm, only:iceprocs,iceflag

implicit none

   type (micro_vars) :: micro
   integer, intent(in) :: n1,n2,n3,n4

! Allocate arrays based on options (if necessary)
      if (level >= 0 .and. level .ne. 4) then
         if(iaerosol > 0) then
            allocate (micro%cn1np(n1,n2,n3))
            allocate (micro%cn1mp(n1,n2,n3))
            allocate (micro%cn2np(n1,n2,n3))
            allocate (micro%cn2mp(n1,n2,n3))
         endif
         if(idust > 0) then
            allocate (micro%md1np(n1,n2,n3))
            allocate (micro%md2np(n1,n2,n3))
            allocate (micro%md1mp(n1,n2,n3))
            allocate (micro%md2mp(n1,n2,n3))
            if(idust == 2) allocate (micro%dustfrac(n2,n3))
         endif
         if(isalt > 0) then
            allocate (micro%salt_film_np(n1,n2,n3))
            allocate (micro%salt_jet_np(n1,n2,n3))
            allocate (micro%salt_spum_np(n1,n2,n3))
            allocate (micro%salt_film_mp(n1,n2,n3))
            allocate (micro%salt_jet_mp(n1,n2,n3))
            allocate (micro%salt_spum_mp(n1,n2,n3))
         endif
         if(iabcarb > 0) then
            allocate (micro%abc1np(n1,n2,n3))
            allocate (micro%abc2np(n1,n2,n3))
            allocate (micro%abc1mp(n1,n2,n3))
            allocate (micro%abc2mp(n1,n2,n3))
         endif
      endif

      if (level >= 2 .and. level .ne. 4) then
         allocate (micro%rcp(n1,n2,n3))
      endif

      if (level == 3) then
         if(ipris>=5 .and. (iifn==1.or.iifn==2)) allocate (micro%cifnp(n1,n2,n3))
         if(idriz >= 1)  then
            allocate (micro%rdp(n1,n2,n3))
            allocate (micro%accpd(n2,n3))
            allocate (micro%pcprd(n2,n3))
            allocate (micro%pcpvd(n1,n2,n3))
         endif
         if(irain >= 1)  then
            allocate (micro%rrp(n1,n2,n3))
            allocate (micro%accpr(n2,n3))
            allocate (micro%pcprr(n2,n3))
            allocate (micro%pcpvr(n1,n2,n3))
            allocate (micro%q2(n1,n2,n3))
         endif
         if(ipris >= 1)  then
            allocate (micro%rpp(n1,n2,n3))
            allocate (micro%accpp(n2,n3))
            allocate (micro%pcprp(n2,n3))
            allocate (micro%pcpvp(n1,n2,n3))
         endif
         if(isnow >= 1)  then
            allocate (micro%rsp(n1,n2,n3))
            allocate (micro%accps(n2,n3))
            allocate (micro%pcprs(n2,n3))
            allocate (micro%pcpvs(n1,n2,n3))
         endif
         if(iaggr >= 1)  then
            allocate (micro%rap(n1,n2,n3))
            allocate (micro%accpa(n2,n3))
            allocate (micro%pcpra(n2,n3))
            allocate (micro%pcpva(n1,n2,n3))
         endif
         if(igraup >= 1) then
            allocate (micro%rgp(n1,n2,n3))
            allocate (micro%accpg(n2,n3))
            allocate (micro%pcprg(n2,n3))
            allocate (micro%pcpvg(n1,n2,n3))
            allocate (micro%q6(n1,n2,n3))
         endif
         if(ihail >= 1)  then
            allocate (micro%rhp(n1,n2,n3))
            allocate (micro%accph(n2,n3))
            allocate (micro%pcprh(n2,n3))
            allocate (micro%pcpvh(n1,n2,n3))
            allocate (micro%q7(n1,n2,n3))
         endif

         if(jnmb(1) >= 5) allocate (micro%ccp(n1,n2,n3))
         if(jnmb(8) >= 5) allocate (micro%cdp(n1,n2,n3))
         if(jnmb(2) >= 5) allocate (micro%crp(n1,n2,n3))
         if(jnmb(3) >= 5) allocate (micro%cpp(n1,n2,n3))
         if(jnmb(4) >= 5) allocate (micro%csp(n1,n2,n3))
         if(jnmb(5) >= 5) allocate (micro%cap(n1,n2,n3))
         if(jnmb(6) >= 5) allocate (micro%cgp(n1,n2,n3))
         if(jnmb(7) >= 5) allocate (micro%chp(n1,n2,n3))

         allocate (micro%pcpg(n2,n3))
         allocate (micro%qpcpg(n2,n3))
         allocate (micro%dpcpg(n2,n3))

         if(iccnlev>=2) then
           allocate (micro%regen_aero1_np(n1,n2,n3))
           allocate (micro%regen_aero1_mp(n1,n2,n3))
           allocate (micro%regen_aero2_np(n1,n2,n3))
           allocate (micro%regen_aero2_mp(n1,n2,n3))
           !Aerosol regeneration SOURCE, as an instantaneous rate per second
           !(regen_aero*_nr, _mr) and, when IMBUDGETACC=1, as a never-reset run
           !total (accregen_aero*_n, _m) -- the same convention as every other
           !IMBUDGET process (see regen_budget_step_start/finish, mic_adj.f90).
           if(imbudget>=1) then
             allocate (micro%regen_aero1_nr(n1,n2,n3)) ; micro%regen_aero1_nr = 0.
             allocate (micro%regen_aero1_mr(n1,n2,n3)) ; micro%regen_aero1_mr = 0.
             allocate (micro%regen_aero2_nr(n1,n2,n3)) ; micro%regen_aero2_nr = 0.
             allocate (micro%regen_aero2_mr(n1,n2,n3)) ; micro%regen_aero2_mr = 0.
             if(imbudgetacc==1) then
               allocate (micro%accregen_aero1_n(n1,n2,n3)) ; micro%accregen_aero1_n = 0.
               allocate (micro%accregen_aero1_m(n1,n2,n3)) ; micro%accregen_aero1_m = 0.
               allocate (micro%accregen_aero2_n(n1,n2,n3)) ; micro%accregen_aero2_n = 0.
               allocate (micro%accregen_aero2_m(n1,n2,n3)) ; micro%accregen_aero2_m = 0.
             endif
           endif
           if(jnmb(1) >= 1) allocate (micro%cnmcp(n1,n2,n3))
           if(jnmb(8) >= 1) allocate (micro%cnmdp(n1,n2,n3))
           if(jnmb(2) >= 1) allocate (micro%cnmrp(n1,n2,n3))
           if(jnmb(3) >= 1) allocate (micro%cnmpp(n1,n2,n3))
           if(jnmb(4) >= 1) allocate (micro%cnmsp(n1,n2,n3))
           if(jnmb(5) >= 1) allocate (micro%cnmap(n1,n2,n3))
           if(jnmb(6) >= 1) allocate (micro%cnmgp(n1,n2,n3))
           if(jnmb(7) >= 1) allocate (micro%cnmhp(n1,n2,n3))
           allocate (micro%accpaero(n2,n3))
           allocate (micro%pcpraero(n2,n3))
           if(itrkdust==1 .and. idust>0) then
            allocate (micro%accpdust(n2,n3))
            allocate (micro%pcprdust(n2,n3))
            if(jnmb(1) >= 1) allocate (micro%dnmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dnmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dnmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dnmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dnmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dnmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dnmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dnmhp(n1,n2,n3))
           endif
           if(itrkdustifn==1 .and. idust>0) then
            if(jnmb(1) >= 1) allocate (micro%dincp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dindp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dinrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dinpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dinsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dinap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dingp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dinhp(n1,n2,n3))
           endif
           !Regenerated aerosol mass in hydrometeors, and the part of it
           !that has served as ice nuclei. Mirrors ITRKDUST/ITRKDUSTIFN.
           if(itrkregen==1 .and. iccnlev>=2) then
            allocate (micro%accpregen(n2,n3))
            allocate (micro%pcprregen(n2,n3))
            if(jnmb(1) >= 1) allocate (micro%rnmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%rnmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%rnmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%rnmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%rnmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%rnmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%rnmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%rnmhp(n1,n2,n3))
            if(jnmb(1) >= 1) allocate (micro%rincp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%rindp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%rinrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%rinpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%rinsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%rinap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%ringp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%rinhp(n1,n2,n3))
           endif
           if(itrkepsilon==1) then
            allocate (micro%resol_aero1_mp(n1,n2,n3))
            allocate (micro%resol_aero2_mp(n1,n2,n3))
            if(jnmb(1) >= 1) allocate (micro%snmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%snmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%snmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%snmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%snmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%snmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%snmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%snmhp(n1,n2,n3))
           endif
         endif

         if(iifn==3 .and. iccnlev>=1) then
           if(jnmb(1) >= 5) allocate (micro%ifnnucp(n1,n2,n3))
           if(jnmb(1) >= 5) allocate (micro%immercp(n1,n2,n3))
           if(jnmb(8) >= 5) allocate (micro%immerdp(n1,n2,n3))
           if(jnmb(2) >= 5) allocate (micro%immerrp(n1,n2,n3))
         endif

         !MICRO BUDGET PROCESSES
         !Each process gets an instantaneous rate array (<name>r, per second),
         !rewritten every step, and -- only when IMBUDGETACC=1 -- a run-total
         !array (acc<name>) that accumulates for the whole run and is never
         !reset, like ACCPR. Both start at zero, which also keeps the lateral
         !boundary columns (never visited by the microphysics) at zero.
         if(imbudget>=1) then
           allocate (micro%latheatvapr(n1,n2,n3)) ; micro%latheatvapr = 0.
           allocate (micro%latheatfrzr(n1,n2,n3)) ; micro%latheatfrzr = 0.
           allocate (micro%nuccldrr(n1,n2,n3)) ; micro%nuccldrr = 0.
           allocate (micro%cld2rainr(n1,n2,n3)) ; micro%cld2rainr = 0.
           allocate (micro%ice2rainr(n1,n2,n3)) ; micro%ice2rainr = 0.
           allocate (micro%nucicerr(n1,n2,n3)) ; micro%nucicerr = 0.
           allocate (micro%vapliqr(n1,n2,n3)) ; micro%vapliqr = 0.
           allocate (micro%vapicer(n1,n2,n3)) ; micro%vapicer = 0.
           allocate (micro%evapliqr(n1,n2,n3)) ; micro%evapliqr = 0.
           allocate (micro%evapicer(n1,n2,n3)) ; micro%evapicer = 0.
           allocate (micro%freezingr(n1,n2,n3)) ; micro%freezingr = 0.
           allocate (micro%meltingr(n1,n2,n3)) ; micro%meltingr = 0.
           allocate (micro%melticer(n1,n2,n3)) ; micro%melticer = 0.
           allocate (micro%rimecldr(n1,n2,n3)) ; micro%rimecldr = 0.
           allocate (micro%rain2icer(n1,n2,n3)) ; micro%rain2icer = 0.
           allocate (micro%aggregater(n1,n2,n3)) ; micro%aggregater = 0.
           if(imbudgetacc==1) then
             allocate (micro%acclatheatvap(n1,n2,n3)) ; micro%acclatheatvap = 0.
             allocate (micro%acclatheatfrz(n1,n2,n3)) ; micro%acclatheatfrz = 0.
             allocate (micro%accnuccldr(n1,n2,n3)) ; micro%accnuccldr = 0.
             allocate (micro%acccld2rain(n1,n2,n3)) ; micro%acccld2rain = 0.
             allocate (micro%accice2rain(n1,n2,n3)) ; micro%accice2rain = 0.
             allocate (micro%accnucicer(n1,n2,n3)) ; micro%accnucicer = 0.
             allocate (micro%accvapliq(n1,n2,n3)) ; micro%accvapliq = 0.
             allocate (micro%accvapice(n1,n2,n3)) ; micro%accvapice = 0.
             allocate (micro%accevapliq(n1,n2,n3)) ; micro%accevapliq = 0.
             allocate (micro%accevapice(n1,n2,n3)) ; micro%accevapice = 0.
             allocate (micro%accfreezing(n1,n2,n3)) ; micro%accfreezing = 0.
             allocate (micro%accmelting(n1,n2,n3)) ; micro%accmelting = 0.
             allocate (micro%accmeltice(n1,n2,n3)) ; micro%accmeltice = 0.
             allocate (micro%accrimecld(n1,n2,n3)) ; micro%accrimecld = 0.
             allocate (micro%accrain2ice(n1,n2,n3)) ; micro%accrain2ice = 0.
             allocate (micro%accaggregate(n1,n2,n3)) ; micro%accaggregate = 0.
           endif
         endif
         if(imbudget>=2) then
           allocate (micro%inuchomrr(n1,n2,n3)) ; micro%inuchomrr = 0.
           allocate (micro%inuccontrr(n1,n2,n3)) ; micro%inuccontrr = 0.
           allocate (micro%inucifnrr(n1,n2,n3)) ; micro%inucifnrr = 0.
           allocate (micro%inuchazrr(n1,n2,n3)) ; micro%inuchazrr = 0.
           allocate (micro%vapcldr(n1,n2,n3)) ; micro%vapcldr = 0.
           allocate (micro%vaprainr(n1,n2,n3)) ; micro%vaprainr = 0.
           allocate (micro%vapprisr(n1,n2,n3)) ; micro%vapprisr = 0.
           allocate (micro%vapsnowr(n1,n2,n3)) ; micro%vapsnowr = 0.
           allocate (micro%vapaggrr(n1,n2,n3)) ; micro%vapaggrr = 0.
           allocate (micro%vapgraur(n1,n2,n3)) ; micro%vapgraur = 0.
           allocate (micro%vaphailr(n1,n2,n3)) ; micro%vaphailr = 0.
           allocate (micro%vapdrizr(n1,n2,n3)) ; micro%vapdrizr = 0.
           allocate (micro%evapcldr(n1,n2,n3)) ; micro%evapcldr = 0.
           allocate (micro%evaprainr(n1,n2,n3)) ; micro%evaprainr = 0.
           allocate (micro%evapprisr(n1,n2,n3)) ; micro%evapprisr = 0.
           allocate (micro%evapsnowr(n1,n2,n3)) ; micro%evapsnowr = 0.
           allocate (micro%evapaggrr(n1,n2,n3)) ; micro%evapaggrr = 0.
           allocate (micro%evapgraur(n1,n2,n3)) ; micro%evapgraur = 0.
           allocate (micro%evaphailr(n1,n2,n3)) ; micro%evaphailr = 0.
           allocate (micro%evapdrizr(n1,n2,n3)) ; micro%evapdrizr = 0.
           allocate (micro%meltprisr(n1,n2,n3)) ; micro%meltprisr = 0.
           allocate (micro%meltsnowr(n1,n2,n3)) ; micro%meltsnowr = 0.
           allocate (micro%meltaggrr(n1,n2,n3)) ; micro%meltaggrr = 0.
           allocate (micro%meltgraur(n1,n2,n3)) ; micro%meltgraur = 0.
           allocate (micro%melthailr(n1,n2,n3)) ; micro%melthailr = 0.
           allocate (micro%rimecldsnowr(n1,n2,n3)) ; micro%rimecldsnowr = 0.
           allocate (micro%rimecldaggrr(n1,n2,n3)) ; micro%rimecldaggrr = 0.
           allocate (micro%rimecldgraur(n1,n2,n3)) ; micro%rimecldgraur = 0.
           allocate (micro%rimecldhailr(n1,n2,n3)) ; micro%rimecldhailr = 0.
           allocate (micro%rain2prr(n1,n2,n3)) ; micro%rain2prr = 0.
           allocate (micro%rain2snr(n1,n2,n3)) ; micro%rain2snr = 0.
           allocate (micro%rain2agr(n1,n2,n3)) ; micro%rain2agr = 0.
           allocate (micro%rain2grr(n1,n2,n3)) ; micro%rain2grr = 0.
           allocate (micro%rain2har(n1,n2,n3)) ; micro%rain2har = 0.
           allocate (micro%aggrselfprisr(n1,n2,n3)) ; micro%aggrselfprisr = 0.
           allocate (micro%aggrselfsnowr(n1,n2,n3)) ; micro%aggrselfsnowr = 0.
           allocate (micro%aggrprissnowr(n1,n2,n3)) ; micro%aggrprissnowr = 0.
           if(imbudgetacc==1) then
             allocate (micro%accinuchomr(n1,n2,n3)) ; micro%accinuchomr = 0.
             allocate (micro%accinuccontr(n1,n2,n3)) ; micro%accinuccontr = 0.
             allocate (micro%accinucifnr(n1,n2,n3)) ; micro%accinucifnr = 0.
             allocate (micro%accinuchazr(n1,n2,n3)) ; micro%accinuchazr = 0.
             allocate (micro%accvapcld(n1,n2,n3)) ; micro%accvapcld = 0.
             allocate (micro%accvaprain(n1,n2,n3)) ; micro%accvaprain = 0.
             allocate (micro%accvappris(n1,n2,n3)) ; micro%accvappris = 0.
             allocate (micro%accvapsnow(n1,n2,n3)) ; micro%accvapsnow = 0.
             allocate (micro%accvapaggr(n1,n2,n3)) ; micro%accvapaggr = 0.
             allocate (micro%accvapgrau(n1,n2,n3)) ; micro%accvapgrau = 0.
             allocate (micro%accvaphail(n1,n2,n3)) ; micro%accvaphail = 0.
             allocate (micro%accvapdriz(n1,n2,n3)) ; micro%accvapdriz = 0.
             allocate (micro%accevapcld(n1,n2,n3)) ; micro%accevapcld = 0.
             allocate (micro%accevaprain(n1,n2,n3)) ; micro%accevaprain = 0.
             allocate (micro%accevappris(n1,n2,n3)) ; micro%accevappris = 0.
             allocate (micro%accevapsnow(n1,n2,n3)) ; micro%accevapsnow = 0.
             allocate (micro%accevapaggr(n1,n2,n3)) ; micro%accevapaggr = 0.
             allocate (micro%accevapgrau(n1,n2,n3)) ; micro%accevapgrau = 0.
             allocate (micro%accevaphail(n1,n2,n3)) ; micro%accevaphail = 0.
             allocate (micro%accevapdriz(n1,n2,n3)) ; micro%accevapdriz = 0.
             allocate (micro%accmeltpris(n1,n2,n3)) ; micro%accmeltpris = 0.
             allocate (micro%accmeltsnow(n1,n2,n3)) ; micro%accmeltsnow = 0.
             allocate (micro%accmeltaggr(n1,n2,n3)) ; micro%accmeltaggr = 0.
             allocate (micro%accmeltgrau(n1,n2,n3)) ; micro%accmeltgrau = 0.
             allocate (micro%accmelthail(n1,n2,n3)) ; micro%accmelthail = 0.
             allocate (micro%accrimecldsnow(n1,n2,n3)) ; micro%accrimecldsnow = 0.
             allocate (micro%accrimecldaggr(n1,n2,n3)) ; micro%accrimecldaggr = 0.
             allocate (micro%accrimecldgrau(n1,n2,n3)) ; micro%accrimecldgrau = 0.
             allocate (micro%accrimecldhail(n1,n2,n3)) ; micro%accrimecldhail = 0.
             allocate (micro%accrain2pr(n1,n2,n3)) ; micro%accrain2pr = 0.
             allocate (micro%accrain2sn(n1,n2,n3)) ; micro%accrain2sn = 0.
             allocate (micro%accrain2ag(n1,n2,n3)) ; micro%accrain2ag = 0.
             allocate (micro%accrain2gr(n1,n2,n3)) ; micro%accrain2gr = 0.
             allocate (micro%accrain2ha(n1,n2,n3)) ; micro%accrain2ha = 0.
             allocate (micro%accaggrselfpris(n1,n2,n3)) ; micro%accaggrselfpris = 0.
             allocate (micro%accaggrselfsnow(n1,n2,n3)) ; micro%accaggrselfsnow = 0.
             allocate (micro%accaggrprissnow(n1,n2,n3)) ; micro%accaggrprissnow = 0.
           endif
         endif
         if(imbudget==3 .and. idust>=1) then
           allocate (micro%dust1cldrr(n1,n2,n3)) ; micro%dust1cldrr = 0.
           allocate (micro%dust2cldrr(n1,n2,n3)) ; micro%dust2cldrr = 0.
           allocate (micro%dust1drzrr(n1,n2,n3)) ; micro%dust1drzrr = 0.
           allocate (micro%dust2drzrr(n1,n2,n3)) ; micro%dust2drzrr = 0.
           if(imbudgetacc==1) then
             allocate (micro%accdust1cldr(n1,n2,n3)) ; micro%accdust1cldr = 0.
             allocate (micro%accdust2cldr(n1,n2,n3)) ; micro%accdust2cldr = 0.
             allocate (micro%accdust1drzr(n1,n2,n3)) ; micro%accdust1drzr = 0.
             allocate (micro%accdust2drzr(n1,n2,n3)) ; micro%accdust2drzr = 0.
           endif
         endif
         if(imbudget==3 .and. iccnlev>=2) then
           allocate (micro%regen1cldrr(n1,n2,n3)) ; micro%regen1cldrr = 0.
           allocate (micro%regen2cldrr(n1,n2,n3)) ; micro%regen2cldrr = 0.
           allocate (micro%regen1drzrr(n1,n2,n3)) ; micro%regen1drzrr = 0.
           allocate (micro%regen2drzrr(n1,n2,n3)) ; micro%regen2drzrr = 0.
           if(imbudgetacc==1) then
             allocate (micro%accregen1cldr(n1,n2,n3)) ; micro%accregen1cldr = 0.
             allocate (micro%accregen2cldr(n1,n2,n3)) ; micro%accregen2cldr = 0.
             allocate (micro%accregen1drzr(n1,n2,n3)) ; micro%accregen1drzr = 0.
             allocate (micro%accregen2drzr(n1,n2,n3)) ; micro%accregen2drzr = 0.
           endif
         endif
      endif
      if (level == 4) then
         allocate (micro%fncn(n1,n2,n3,n4))
         allocate (micro%ffcd(n1,n2,n3,n4))
         allocate (micro%accpr(n2,n3))
         allocate (micro%pcprr(n2,n3))
         allocate (micro%pcpvr(n1,n2,n3))
         allocate (micro%pcpg(n2,n3))
         allocate (micro%qpcpg(n2,n3))
         allocate (micro%dpcpg(n2,n3))

         if (iceprocs == 1) then
            if (ipris == 1 .or. ipris >= 4) then
               allocate (micro%ffic(n1,n2,n3,n4))
               allocate (micro%accpic(n2,n3))
               allocate (micro%pcpric(n2,n3))
               allocate (micro%pcpvic(n1,n2,n3))
            endif

            if (ipris == 2 .or. ipris >= 4) then
               allocate (micro%ffip(n1,n2,n3,n4))
               allocate (micro%accpip(n2,n3))
               allocate (micro%pcprip(n2,n3))
               allocate (micro%pcpvip(n1,n2,n3))
            endif

            if (ipris >= 3) then 
               allocate (micro%ffid(n1,n2,n3,n4))
               allocate (micro%accpid(n2,n3))
               allocate (micro%pcprid(n2,n3))
               allocate (micro%pcpvid(n1,n2,n3))
            endif

            if (ipris >= 1) then
               allocate (micro%accpp(n2,n3))
               allocate (micro%pcprp(n2,n3))
               allocate (micro%pcpvp(n1,n2,n3))
               allocate (micro%accps(n2,n3))
               allocate (micro%pcprs(n2,n3))
               allocate (micro%pcpvs(n1,n2,n3))
            endif

            allocate (micro%ffsn(n1,n2,n3,n4))
            allocate (micro%accpa(n2,n3))
            allocate (micro%pcpra(n2,n3))
            allocate (micro%pcpva(n1,n2,n3))

            if (igraup > 0) then
               allocate (micro%ffgl(n1,n2,n3,n4))
               allocate (micro%accpg(n2,n3))
               allocate (micro%pcprg(n2,n3))
               allocate (micro%pcpvg(n1,n2,n3))
            endif

            if (ihail > 0) then
               allocate (micro%ffhl(n1,n2,n3,n4))
               allocate (micro%accph(n2,n3))
               allocate (micro%pcprh(n2,n3))
               allocate (micro%pcpvh(n1,n2,n3))
            endif

            if (iifn == 2) allocate (micro%ffin(n1,n2,n3,n4))
         endif

         allocate (micro%rv_old(n1,n2,n3))
         allocate (micro%t_old(n1,n2,n3))

         if(imbudget>=1) then
            allocate (micro%latheatvapr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%acclatheatvap(n1,n2,n3))
            allocate (micro%nuccldrr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%accnuccldr(n1,n2,n3))
            allocate (micro%nuccldcr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%accnuccldc(n1,n2,n3))
            allocate (micro%cld2rainr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%acccld2rain(n1,n2,n3))
            allocate (micro%vapliqr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%accvapliq(n1,n2,n3))

            if(iceprocs==1) then
               allocate (micro%latheatfrzr(n1,n2,n3))
              ! allocate (micro%accice2rain(n1,n2,n3))
               allocate (micro%nucicerr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accnucicer(n1,n2,n3))
               allocate (micro%nucicecr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accnucicec(n1,n2,n3))
               allocate (micro%vapicer(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accvapice(n1,n2,n3))
               allocate (micro%melticer(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accmeltice(n1,n2,n3))
               allocate (micro%rimecldr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accrimecld(n1,n2,n3))
               allocate (micro%rain2icer(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accrain2ice(n1,n2,n3))
               allocate (micro%aggregater(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accaggregate(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%acclatheatfrz(n1,n2,n3))
            endif
         endif
         if(imbudget>=2) then
            allocate (micro%vapcldr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%accvapcld(n1,n2,n3))
            allocate (micro%vaprainr(n1,n2,n3))
            if (imbudgetacc==1) allocate (micro%accvaprain(n1,n2,n3))
            if(iceprocs==1) then
              !Budget variables are the same as for the bulk
              !microphysics scheme. Those that are commented 
              !out have not been added to the bin microphysics.
               allocate (micro%inuchomrr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accinuchomr(n1,n2,n3))
               allocate (micro%inuchomcr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accinuchomc(n1,n2,n3))
              ! allocate (micro%accinuccontr(n1,n2,n3))
              ! allocate (micro%accinuccontc(n1,n2,n3))
               allocate (micro%inucifnrr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accinucifnr(n1,n2,n3))
               allocate (micro%inucifncr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accinucifnc(n1,n2,n3))
              ! allocate (micro%accinuchazr(n1,n2,n3))
              ! allocate (micro%accinuchazc(n1,n2,n3))
               if (ipris > 0) allocate (micro%vapprisr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accvappris(n1,n2,n3))
               if (ipris > 0) allocate (micro%vapsnowr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accvapsnow(n1,n2,n3))
               allocate (micro%vapaggrr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accvapaggr(n1,n2,n3))
               if (igraup > 0) allocate (micro%vapgraur(n1,n2,n3))
               if (igraup > 0 .and. imbudgetacc==1) allocate (micro%accvapgrau(n1,n2,n3))
               if (ihail > 0) allocate (micro%vaphailr(n1,n2,n3))
               if (ihail > 0 .and. imbudgetacc==1) allocate (micro%accvaphail(n1,n2,n3))
               if (ipris > 0) allocate (micro%meltprisr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accmeltpris(n1,n2,n3))
               if (ipris > 0) allocate (micro%meltsnowr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accmeltsnow(n1,n2,n3))
               allocate (micro%meltaggrr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accmeltaggr(n1,n2,n3))
               if (igraup > 0) allocate (micro%meltgraur(n1,n2,n3))
               if (igraup > 0 .and. imbudgetacc==1) allocate (micro%accmeltgrau(n1,n2,n3))
               if (ihail > 0) allocate (micro%melthailr(n1,n2,n3))
               if (ihail > 0 .and. imbudgetacc==1) allocate (micro%accmelthail(n1,n2,n3))
               if (ipris > 0) allocate (micro%rimecldsnowr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accrimecldsnow(n1,n2,n3))
               allocate (micro%rimecldaggrr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accrimecldaggr(n1,n2,n3))
               if (igraup > 0) allocate (micro%rimecldgraur(n1,n2,n3))
               if (igraup > 0 .and. imbudgetacc==1) allocate (micro%accrimecldgrau(n1,n2,n3))
               if (ihail > 0) allocate (micro%rimecldhailr(n1,n2,n3))
               if (ihail > 0 .and. imbudgetacc==1) allocate (micro%accrimecldhail(n1,n2,n3))
              ! allocate (micro%accrain2pr(n1,n2,n3))
               if (ipris > 0) allocate (micro%rain2snr(n1,n2,n3))
               if (ipris > 0 .and. imbudgetacc==1) allocate (micro%accrain2sn(n1,n2,n3))
               allocate (micro%rain2agr(n1,n2,n3))
               if (imbudgetacc==1) allocate (micro%accrain2ag(n1,n2,n3))
               if (igraup > 0) allocate (micro%rain2grr(n1,n2,n3))
               if (igraup > 0 .and. imbudgetacc==1) allocate (micro%accrain2gr(n1,n2,n3))
               if (ihail > 0) allocate (micro%rain2har(n1,n2,n3))
               if (ihail > 0 .and. imbudgetacc==1) allocate (micro%accrain2ha(n1,n2,n3))
              ! allocate (micro%accaggrselfpris(n1,n2,n3))
              ! allocate (micro%accaggrselfsnow(n1,n2,n3))
              ! allocate (micro%accaggrprissnow(n1,n2,n3))
            endif
         endif
      endif


   !Start every microphysics budget rate and run total at zero
   if (allocated(micro%latheatvapr)) micro%latheatvapr = 0.
   if (allocated(micro%acclatheatvap)) micro%acclatheatvap = 0.
   if (allocated(micro%latheatfrzr)) micro%latheatfrzr = 0.
   if (allocated(micro%acclatheatfrz)) micro%acclatheatfrz = 0.
   if (allocated(micro%nuccldrr)) micro%nuccldrr = 0.
   if (allocated(micro%accnuccldr)) micro%accnuccldr = 0.
   if (allocated(micro%cld2rainr)) micro%cld2rainr = 0.
   if (allocated(micro%acccld2rain)) micro%acccld2rain = 0.
   if (allocated(micro%ice2rainr)) micro%ice2rainr = 0.
   if (allocated(micro%accice2rain)) micro%accice2rain = 0.
   if (allocated(micro%nucicerr)) micro%nucicerr = 0.
   if (allocated(micro%accnucicer)) micro%accnucicer = 0.
   if (allocated(micro%vapliqr)) micro%vapliqr = 0.
   if (allocated(micro%accvapliq)) micro%accvapliq = 0.
   if (allocated(micro%vapicer)) micro%vapicer = 0.
   if (allocated(micro%accvapice)) micro%accvapice = 0.
   if (allocated(micro%evapliqr)) micro%evapliqr = 0.
   if (allocated(micro%accevapliq)) micro%accevapliq = 0.
   if (allocated(micro%evapicer)) micro%evapicer = 0.
   if (allocated(micro%accevapice)) micro%accevapice = 0.
   if (allocated(micro%freezingr)) micro%freezingr = 0.
   if (allocated(micro%accfreezing)) micro%accfreezing = 0.
   if (allocated(micro%meltingr)) micro%meltingr = 0.
   if (allocated(micro%accmelting)) micro%accmelting = 0.
   if (allocated(micro%melticer)) micro%melticer = 0.
   if (allocated(micro%accmeltice)) micro%accmeltice = 0.
   if (allocated(micro%rimecldr)) micro%rimecldr = 0.
   if (allocated(micro%accrimecld)) micro%accrimecld = 0.
   if (allocated(micro%rain2icer)) micro%rain2icer = 0.
   if (allocated(micro%accrain2ice)) micro%accrain2ice = 0.
   if (allocated(micro%aggregater)) micro%aggregater = 0.
   if (allocated(micro%accaggregate)) micro%accaggregate = 0.
   if (allocated(micro%inuchomrr)) micro%inuchomrr = 0.
   if (allocated(micro%accinuchomr)) micro%accinuchomr = 0.
   if (allocated(micro%inuccontrr)) micro%inuccontrr = 0.
   if (allocated(micro%accinuccontr)) micro%accinuccontr = 0.
   if (allocated(micro%inucifnrr)) micro%inucifnrr = 0.
   if (allocated(micro%accinucifnr)) micro%accinucifnr = 0.
   if (allocated(micro%inuchazrr)) micro%inuchazrr = 0.
   if (allocated(micro%accinuchazr)) micro%accinuchazr = 0.
   if (allocated(micro%vapcldr)) micro%vapcldr = 0.
   if (allocated(micro%accvapcld)) micro%accvapcld = 0.
   if (allocated(micro%vaprainr)) micro%vaprainr = 0.
   if (allocated(micro%accvaprain)) micro%accvaprain = 0.
   if (allocated(micro%vapprisr)) micro%vapprisr = 0.
   if (allocated(micro%accvappris)) micro%accvappris = 0.
   if (allocated(micro%vapsnowr)) micro%vapsnowr = 0.
   if (allocated(micro%accvapsnow)) micro%accvapsnow = 0.
   if (allocated(micro%vapaggrr)) micro%vapaggrr = 0.
   if (allocated(micro%accvapaggr)) micro%accvapaggr = 0.
   if (allocated(micro%vapgraur)) micro%vapgraur = 0.
   if (allocated(micro%accvapgrau)) micro%accvapgrau = 0.
   if (allocated(micro%vaphailr)) micro%vaphailr = 0.
   if (allocated(micro%accvaphail)) micro%accvaphail = 0.
   if (allocated(micro%vapdrizr)) micro%vapdrizr = 0.
   if (allocated(micro%accvapdriz)) micro%accvapdriz = 0.
   if (allocated(micro%evapcldr)) micro%evapcldr = 0.
   if (allocated(micro%accevapcld)) micro%accevapcld = 0.
   if (allocated(micro%evaprainr)) micro%evaprainr = 0.
   if (allocated(micro%accevaprain)) micro%accevaprain = 0.
   if (allocated(micro%evapprisr)) micro%evapprisr = 0.
   if (allocated(micro%accevappris)) micro%accevappris = 0.
   if (allocated(micro%evapsnowr)) micro%evapsnowr = 0.
   if (allocated(micro%accevapsnow)) micro%accevapsnow = 0.
   if (allocated(micro%evapaggrr)) micro%evapaggrr = 0.
   if (allocated(micro%accevapaggr)) micro%accevapaggr = 0.
   if (allocated(micro%evapgraur)) micro%evapgraur = 0.
   if (allocated(micro%accevapgrau)) micro%accevapgrau = 0.
   if (allocated(micro%evaphailr)) micro%evaphailr = 0.
   if (allocated(micro%accevaphail)) micro%accevaphail = 0.
   if (allocated(micro%evapdrizr)) micro%evapdrizr = 0.
   if (allocated(micro%accevapdriz)) micro%accevapdriz = 0.
   if (allocated(micro%meltprisr)) micro%meltprisr = 0.
   if (allocated(micro%accmeltpris)) micro%accmeltpris = 0.
   if (allocated(micro%meltsnowr)) micro%meltsnowr = 0.
   if (allocated(micro%accmeltsnow)) micro%accmeltsnow = 0.
   if (allocated(micro%meltaggrr)) micro%meltaggrr = 0.
   if (allocated(micro%accmeltaggr)) micro%accmeltaggr = 0.
   if (allocated(micro%meltgraur)) micro%meltgraur = 0.
   if (allocated(micro%accmeltgrau)) micro%accmeltgrau = 0.
   if (allocated(micro%melthailr)) micro%melthailr = 0.
   if (allocated(micro%accmelthail)) micro%accmelthail = 0.
   if (allocated(micro%rimecldsnowr)) micro%rimecldsnowr = 0.
   if (allocated(micro%accrimecldsnow)) micro%accrimecldsnow = 0.
   if (allocated(micro%rimecldaggrr)) micro%rimecldaggrr = 0.
   if (allocated(micro%accrimecldaggr)) micro%accrimecldaggr = 0.
   if (allocated(micro%rimecldgraur)) micro%rimecldgraur = 0.
   if (allocated(micro%accrimecldgrau)) micro%accrimecldgrau = 0.
   if (allocated(micro%rimecldhailr)) micro%rimecldhailr = 0.
   if (allocated(micro%accrimecldhail)) micro%accrimecldhail = 0.
   if (allocated(micro%rain2prr)) micro%rain2prr = 0.
   if (allocated(micro%accrain2pr)) micro%accrain2pr = 0.
   if (allocated(micro%rain2snr)) micro%rain2snr = 0.
   if (allocated(micro%accrain2sn)) micro%accrain2sn = 0.
   if (allocated(micro%rain2agr)) micro%rain2agr = 0.
   if (allocated(micro%accrain2ag)) micro%accrain2ag = 0.
   if (allocated(micro%rain2grr)) micro%rain2grr = 0.
   if (allocated(micro%accrain2gr)) micro%accrain2gr = 0.
   if (allocated(micro%rain2har)) micro%rain2har = 0.
   if (allocated(micro%accrain2ha)) micro%accrain2ha = 0.
   if (allocated(micro%aggrselfprisr)) micro%aggrselfprisr = 0.
   if (allocated(micro%accaggrselfpris)) micro%accaggrselfpris = 0.
   if (allocated(micro%aggrselfsnowr)) micro%aggrselfsnowr = 0.
   if (allocated(micro%accaggrselfsnow)) micro%accaggrselfsnow = 0.
   if (allocated(micro%aggrprissnowr)) micro%aggrprissnowr = 0.
   if (allocated(micro%accaggrprissnow)) micro%accaggrprissnow = 0.
   if (allocated(micro%dust1cldrr)) micro%dust1cldrr = 0.
   if (allocated(micro%accdust1cldr)) micro%accdust1cldr = 0.
   if (allocated(micro%dust2cldrr)) micro%dust2cldrr = 0.
   if (allocated(micro%accdust2cldr)) micro%accdust2cldr = 0.
   if (allocated(micro%dust1drzrr)) micro%dust1drzrr = 0.
   if (allocated(micro%accdust1drzr)) micro%accdust1drzr = 0.
   if (allocated(micro%dust2drzrr)) micro%dust2drzrr = 0.
   if (allocated(micro%accdust2drzr)) micro%accdust2drzr = 0.
   if (allocated(micro%regen1cldrr)) micro%regen1cldrr = 0.
   if (allocated(micro%accregen1cldr)) micro%accregen1cldr = 0.
   if (allocated(micro%regen2cldrr)) micro%regen2cldrr = 0.
   if (allocated(micro%accregen2cldr)) micro%accregen2cldr = 0.
   if (allocated(micro%regen1drzrr)) micro%regen1drzrr = 0.
   if (allocated(micro%accregen1drzr)) micro%accregen1drzr = 0.
   if (allocated(micro%regen2drzrr)) micro%regen2drzrr = 0.
   if (allocated(micro%accregen2drzr)) micro%accregen2drzr = 0.
   if (allocated(micro%nuccldcr)) micro%nuccldcr = 0.
   if (allocated(micro%accnuccldc)) micro%accnuccldc = 0.
   if (allocated(micro%nucicecr)) micro%nucicecr = 0.
   if (allocated(micro%accnucicec)) micro%accnucicec = 0.
   if (allocated(micro%inuchomcr)) micro%inuchomcr = 0.
   if (allocated(micro%accinuchomc)) micro%accinuchomc = 0.
   if (allocated(micro%inucifncr)) micro%inucifncr = 0.
   if (allocated(micro%accinucifnc)) micro%accinucifnc = 0.

return
END SUBROUTINE alloc_micro

!##############################################################################
Subroutine dealloc_micro (micro)

implicit none

   type (micro_vars) :: micro
   
   if (allocated(micro%rcp))     deallocate (micro%rcp)
   if (allocated(micro%rdp))     deallocate (micro%rdp)
   if (allocated(micro%rrp))     deallocate (micro%rrp)
   if (allocated(micro%rpp))     deallocate (micro%rpp)
   if (allocated(micro%rsp))     deallocate (micro%rsp)
   if (allocated(micro%rap))     deallocate (micro%rap)
   if (allocated(micro%rgp))     deallocate (micro%rgp)
   if (allocated(micro%rhp))     deallocate (micro%rhp)
   if (allocated(micro%ccp))     deallocate (micro%ccp)
   if (allocated(micro%cdp))     deallocate (micro%cdp)
   if (allocated(micro%crp))     deallocate (micro%crp)
   if (allocated(micro%cpp))     deallocate (micro%cpp)
   if (allocated(micro%csp))     deallocate (micro%csp)
   if (allocated(micro%cap))     deallocate (micro%cap)
   if (allocated(micro%cgp))     deallocate (micro%cgp)
   if (allocated(micro%chp))     deallocate (micro%chp)
   if (allocated(micro%cifnp))   deallocate (micro%cifnp)
   if (allocated(micro%q2))      deallocate (micro%q2)
   if (allocated(micro%q6))      deallocate (micro%q6)
   if (allocated(micro%q7))      deallocate (micro%q7)

   if (allocated(micro%cn1np))   deallocate (micro%cn1np)
   if (allocated(micro%cn2np))   deallocate (micro%cn2np)
   if (allocated(micro%cn1mp))   deallocate (micro%cn1mp)
   if (allocated(micro%cn2mp))   deallocate (micro%cn2mp)
   if (allocated(micro%md1np))   deallocate (micro%md1np)
   if (allocated(micro%md2np))   deallocate (micro%md2np)
   if (allocated(micro%md1mp))   deallocate (micro%md1mp)
   if (allocated(micro%md2mp))   deallocate (micro%md2mp)
   if (allocated(micro%dustfrac))deallocate (micro%dustfrac)
   if (allocated(micro%salt_film_np))  deallocate (micro%salt_film_np)
   if (allocated(micro%salt_jet_np))   deallocate (micro%salt_jet_np)
   if (allocated(micro%salt_spum_np))  deallocate (micro%salt_spum_np)
   if (allocated(micro%salt_film_mp))  deallocate (micro%salt_film_mp)
   if (allocated(micro%salt_jet_mp))   deallocate (micro%salt_jet_mp)
   if (allocated(micro%salt_spum_mp))  deallocate (micro%salt_spum_mp)
   if (allocated(micro%abc1np))   deallocate (micro%abc1np)
   if (allocated(micro%abc2np))   deallocate (micro%abc2np)
   if (allocated(micro%abc1mp))   deallocate (micro%abc1mp)
   if (allocated(micro%abc2mp))   deallocate (micro%abc2mp)
   if (allocated(micro%regen_aero1_np)) deallocate (micro%regen_aero1_np)
   if (allocated(micro%regen_aero1_nr)) deallocate (micro%regen_aero1_nr)
   if (allocated(micro%regen_aero1_mr)) deallocate (micro%regen_aero1_mr)
   if (allocated(micro%regen_aero2_nr)) deallocate (micro%regen_aero2_nr)
   if (allocated(micro%regen_aero2_mr)) deallocate (micro%regen_aero2_mr)
   if (allocated(micro%accregen_aero1_n)) deallocate (micro%accregen_aero1_n)
   if (allocated(micro%accregen_aero1_m)) deallocate (micro%accregen_aero1_m)
   if (allocated(micro%accregen_aero2_n)) deallocate (micro%accregen_aero2_n)
   if (allocated(micro%accregen_aero2_m)) deallocate (micro%accregen_aero2_m)
   if (allocated(micro%regen_aero1_mp)) deallocate (micro%regen_aero1_mp)
   if (allocated(micro%regen_aero2_np)) deallocate (micro%regen_aero2_np)
   if (allocated(micro%regen_aero2_mp)) deallocate (micro%regen_aero2_mp)

   if (allocated(micro%immercp)) deallocate (micro%immercp)
   if (allocated(micro%immerdp)) deallocate (micro%immerdp)
   if (allocated(micro%immerrp)) deallocate (micro%immerrp)
   if (allocated(micro%ifnnucp)) deallocate (micro%ifnnucp)

   if (allocated(micro%cnmcp))   deallocate (micro%cnmcp)
   if (allocated(micro%cnmdp))   deallocate (micro%cnmdp)
   if (allocated(micro%cnmrp))   deallocate (micro%cnmrp)
   if (allocated(micro%cnmpp))   deallocate (micro%cnmpp)
   if (allocated(micro%cnmsp))   deallocate (micro%cnmsp)
   if (allocated(micro%cnmap))   deallocate (micro%cnmap)
   if (allocated(micro%cnmgp))   deallocate (micro%cnmgp)
   if (allocated(micro%cnmhp))   deallocate (micro%cnmhp)
   if (allocated(micro%accpaero))deallocate (micro%accpaero)
   if (allocated(micro%pcpraero))deallocate (micro%pcpraero)

   if (allocated(micro%dnmcp))   deallocate (micro%dnmcp)
   if (allocated(micro%dnmdp))   deallocate (micro%dnmdp)
   if (allocated(micro%dnmrp))   deallocate (micro%dnmrp)
   if (allocated(micro%dnmpp))   deallocate (micro%dnmpp)
   if (allocated(micro%dnmsp))   deallocate (micro%dnmsp)
   if (allocated(micro%dnmap))   deallocate (micro%dnmap)
   if (allocated(micro%dnmgp))   deallocate (micro%dnmgp)
   if (allocated(micro%dnmhp))   deallocate (micro%dnmhp)
   if (allocated(micro%accpdust))deallocate (micro%accpdust)
   if (allocated(micro%pcprdust))deallocate (micro%pcprdust)

   if (allocated(micro%dincp))   deallocate (micro%dincp)
   if (allocated(micro%dindp))   deallocate (micro%dindp)
   if (allocated(micro%dinrp))   deallocate (micro%dinrp)
   if (allocated(micro%dinpp))   deallocate (micro%dinpp)
   if (allocated(micro%dinsp))   deallocate (micro%dinsp)
   if (allocated(micro%dinap))   deallocate (micro%dinap)
   if (allocated(micro%dingp))   deallocate (micro%dingp)
   if (allocated(micro%dinhp))   deallocate (micro%dinhp)

   if (allocated(micro%rnmcp))   deallocate (micro%rnmcp)
   if (allocated(micro%rnmdp))   deallocate (micro%rnmdp)
   if (allocated(micro%rnmrp))   deallocate (micro%rnmrp)
   if (allocated(micro%rnmpp))   deallocate (micro%rnmpp)
   if (allocated(micro%rnmsp))   deallocate (micro%rnmsp)
   if (allocated(micro%rnmap))   deallocate (micro%rnmap)
   if (allocated(micro%rnmgp))   deallocate (micro%rnmgp)
   if (allocated(micro%rnmhp))   deallocate (micro%rnmhp)
   if (allocated(micro%rincp))   deallocate (micro%rincp)
   if (allocated(micro%rindp))   deallocate (micro%rindp)
   if (allocated(micro%rinrp))   deallocate (micro%rinrp)
   if (allocated(micro%rinpp))   deallocate (micro%rinpp)
   if (allocated(micro%rinsp))   deallocate (micro%rinsp)
   if (allocated(micro%rinap))   deallocate (micro%rinap)
   if (allocated(micro%ringp))   deallocate (micro%ringp)
   if (allocated(micro%rinhp))   deallocate (micro%rinhp)
   if (allocated(micro%accpregen))deallocate (micro%accpregen)
   if (allocated(micro%pcprregen))deallocate (micro%pcprregen)

   if (allocated(micro%snmcp))   deallocate (micro%snmcp)
   if (allocated(micro%snmdp))   deallocate (micro%snmdp)
   if (allocated(micro%snmrp))   deallocate (micro%snmrp)
   if (allocated(micro%snmpp))   deallocate (micro%snmpp)
   if (allocated(micro%snmsp))   deallocate (micro%snmsp)
   if (allocated(micro%snmap))   deallocate (micro%snmap)
   if (allocated(micro%snmgp))   deallocate (micro%snmgp)
   if (allocated(micro%snmhp))   deallocate (micro%snmhp)
   if (allocated(micro%resol_aero1_mp)) deallocate (micro%resol_aero1_mp)
   if (allocated(micro%resol_aero2_mp)) deallocate (micro%resol_aero2_mp)

   if (allocated(micro%pcpvr))   deallocate (micro%pcpvr)
   if (allocated(micro%pcpvp))   deallocate (micro%pcpvp)
   if (allocated(micro%pcpvs))   deallocate (micro%pcpvs)
   if (allocated(micro%pcpva))   deallocate (micro%pcpva)
   if (allocated(micro%pcpvg))   deallocate (micro%pcpvg)
   if (allocated(micro%pcpvh))   deallocate (micro%pcpvh)
   if (allocated(micro%pcpvd))   deallocate (micro%pcpvd)

   if (allocated(micro%accpr))   deallocate (micro%accpr)
   if (allocated(micro%accpp))   deallocate (micro%accpp)
   if (allocated(micro%accps))   deallocate (micro%accps)
   if (allocated(micro%accpa))   deallocate (micro%accpa)
   if (allocated(micro%accpg))   deallocate (micro%accpg)
   if (allocated(micro%accph))   deallocate (micro%accph)
   if (allocated(micro%accpd))   deallocate (micro%accpd)
   if (allocated(micro%pcprr))   deallocate (micro%pcprr)
   if (allocated(micro%pcprp))   deallocate (micro%pcprp)
   if (allocated(micro%pcprs))   deallocate (micro%pcprs)
   if (allocated(micro%pcpra))   deallocate (micro%pcpra)
   if (allocated(micro%pcprg))   deallocate (micro%pcprg)
   if (allocated(micro%pcprh))   deallocate (micro%pcprh)
   if (allocated(micro%pcprd))   deallocate (micro%pcprd)
   if (allocated(micro%pcpg))    deallocate (micro%pcpg)
   if (allocated(micro%qpcpg))   deallocate (micro%qpcpg)
   if (allocated(micro%dpcpg))   deallocate (micro%dpcpg)

    !MICRO BUDGET PROCESS RATES AND RUN TOTALS
    if (allocated(micro%latheatvapr)) deallocate (micro%latheatvapr)
    if (allocated(micro%acclatheatvap)) deallocate (micro%acclatheatvap)
    if (allocated(micro%latheatfrzr)) deallocate (micro%latheatfrzr)
    if (allocated(micro%acclatheatfrz)) deallocate (micro%acclatheatfrz)
    if (allocated(micro%nuccldrr)) deallocate (micro%nuccldrr)
    if (allocated(micro%accnuccldr)) deallocate (micro%accnuccldr)
    if (allocated(micro%cld2rainr)) deallocate (micro%cld2rainr)
    if (allocated(micro%acccld2rain)) deallocate (micro%acccld2rain)
    if (allocated(micro%ice2rainr)) deallocate (micro%ice2rainr)
    if (allocated(micro%accice2rain)) deallocate (micro%accice2rain)
    if (allocated(micro%nucicerr)) deallocate (micro%nucicerr)
    if (allocated(micro%accnucicer)) deallocate (micro%accnucicer)
    if (allocated(micro%vapliqr)) deallocate (micro%vapliqr)
    if (allocated(micro%accvapliq)) deallocate (micro%accvapliq)
    if (allocated(micro%vapicer)) deallocate (micro%vapicer)
    if (allocated(micro%accvapice)) deallocate (micro%accvapice)
    if (allocated(micro%evapliqr)) deallocate (micro%evapliqr)
    if (allocated(micro%accevapliq)) deallocate (micro%accevapliq)
    if (allocated(micro%evapicer)) deallocate (micro%evapicer)
    if (allocated(micro%accevapice)) deallocate (micro%accevapice)
    if (allocated(micro%freezingr)) deallocate (micro%freezingr)
    if (allocated(micro%accfreezing)) deallocate (micro%accfreezing)
    if (allocated(micro%meltingr)) deallocate (micro%meltingr)
    if (allocated(micro%accmelting)) deallocate (micro%accmelting)
    if (allocated(micro%melticer)) deallocate (micro%melticer)
    if (allocated(micro%accmeltice)) deallocate (micro%accmeltice)
    if (allocated(micro%rimecldr)) deallocate (micro%rimecldr)
    if (allocated(micro%accrimecld)) deallocate (micro%accrimecld)
    if (allocated(micro%rain2icer)) deallocate (micro%rain2icer)
    if (allocated(micro%accrain2ice)) deallocate (micro%accrain2ice)
    if (allocated(micro%aggregater)) deallocate (micro%aggregater)
    if (allocated(micro%accaggregate)) deallocate (micro%accaggregate)
    if (allocated(micro%inuchomrr)) deallocate (micro%inuchomrr)
    if (allocated(micro%accinuchomr)) deallocate (micro%accinuchomr)
    if (allocated(micro%inuccontrr)) deallocate (micro%inuccontrr)
    if (allocated(micro%accinuccontr)) deallocate (micro%accinuccontr)
    if (allocated(micro%inucifnrr)) deallocate (micro%inucifnrr)
    if (allocated(micro%accinucifnr)) deallocate (micro%accinucifnr)
    if (allocated(micro%inuchazrr)) deallocate (micro%inuchazrr)
    if (allocated(micro%accinuchazr)) deallocate (micro%accinuchazr)
    if (allocated(micro%vapcldr)) deallocate (micro%vapcldr)
    if (allocated(micro%accvapcld)) deallocate (micro%accvapcld)
    if (allocated(micro%vaprainr)) deallocate (micro%vaprainr)
    if (allocated(micro%accvaprain)) deallocate (micro%accvaprain)
    if (allocated(micro%vapprisr)) deallocate (micro%vapprisr)
    if (allocated(micro%accvappris)) deallocate (micro%accvappris)
    if (allocated(micro%vapsnowr)) deallocate (micro%vapsnowr)
    if (allocated(micro%accvapsnow)) deallocate (micro%accvapsnow)
    if (allocated(micro%vapaggrr)) deallocate (micro%vapaggrr)
    if (allocated(micro%accvapaggr)) deallocate (micro%accvapaggr)
    if (allocated(micro%vapgraur)) deallocate (micro%vapgraur)
    if (allocated(micro%accvapgrau)) deallocate (micro%accvapgrau)
    if (allocated(micro%vaphailr)) deallocate (micro%vaphailr)
    if (allocated(micro%accvaphail)) deallocate (micro%accvaphail)
    if (allocated(micro%vapdrizr)) deallocate (micro%vapdrizr)
    if (allocated(micro%accvapdriz)) deallocate (micro%accvapdriz)
    if (allocated(micro%evapcldr)) deallocate (micro%evapcldr)
    if (allocated(micro%accevapcld)) deallocate (micro%accevapcld)
    if (allocated(micro%evaprainr)) deallocate (micro%evaprainr)
    if (allocated(micro%accevaprain)) deallocate (micro%accevaprain)
    if (allocated(micro%evapprisr)) deallocate (micro%evapprisr)
    if (allocated(micro%accevappris)) deallocate (micro%accevappris)
    if (allocated(micro%evapsnowr)) deallocate (micro%evapsnowr)
    if (allocated(micro%accevapsnow)) deallocate (micro%accevapsnow)
    if (allocated(micro%evapaggrr)) deallocate (micro%evapaggrr)
    if (allocated(micro%accevapaggr)) deallocate (micro%accevapaggr)
    if (allocated(micro%evapgraur)) deallocate (micro%evapgraur)
    if (allocated(micro%accevapgrau)) deallocate (micro%accevapgrau)
    if (allocated(micro%evaphailr)) deallocate (micro%evaphailr)
    if (allocated(micro%accevaphail)) deallocate (micro%accevaphail)
    if (allocated(micro%evapdrizr)) deallocate (micro%evapdrizr)
    if (allocated(micro%accevapdriz)) deallocate (micro%accevapdriz)
    if (allocated(micro%meltprisr)) deallocate (micro%meltprisr)
    if (allocated(micro%accmeltpris)) deallocate (micro%accmeltpris)
    if (allocated(micro%meltsnowr)) deallocate (micro%meltsnowr)
    if (allocated(micro%accmeltsnow)) deallocate (micro%accmeltsnow)
    if (allocated(micro%meltaggrr)) deallocate (micro%meltaggrr)
    if (allocated(micro%accmeltaggr)) deallocate (micro%accmeltaggr)
    if (allocated(micro%meltgraur)) deallocate (micro%meltgraur)
    if (allocated(micro%accmeltgrau)) deallocate (micro%accmeltgrau)
    if (allocated(micro%melthailr)) deallocate (micro%melthailr)
    if (allocated(micro%accmelthail)) deallocate (micro%accmelthail)
    if (allocated(micro%rimecldsnowr)) deallocate (micro%rimecldsnowr)
    if (allocated(micro%accrimecldsnow)) deallocate (micro%accrimecldsnow)
    if (allocated(micro%rimecldaggrr)) deallocate (micro%rimecldaggrr)
    if (allocated(micro%accrimecldaggr)) deallocate (micro%accrimecldaggr)
    if (allocated(micro%rimecldgraur)) deallocate (micro%rimecldgraur)
    if (allocated(micro%accrimecldgrau)) deallocate (micro%accrimecldgrau)
    if (allocated(micro%rimecldhailr)) deallocate (micro%rimecldhailr)
    if (allocated(micro%accrimecldhail)) deallocate (micro%accrimecldhail)
    if (allocated(micro%rain2prr)) deallocate (micro%rain2prr)
    if (allocated(micro%accrain2pr)) deallocate (micro%accrain2pr)
    if (allocated(micro%rain2snr)) deallocate (micro%rain2snr)
    if (allocated(micro%accrain2sn)) deallocate (micro%accrain2sn)
    if (allocated(micro%rain2agr)) deallocate (micro%rain2agr)
    if (allocated(micro%accrain2ag)) deallocate (micro%accrain2ag)
    if (allocated(micro%rain2grr)) deallocate (micro%rain2grr)
    if (allocated(micro%accrain2gr)) deallocate (micro%accrain2gr)
    if (allocated(micro%rain2har)) deallocate (micro%rain2har)
    if (allocated(micro%accrain2ha)) deallocate (micro%accrain2ha)
    if (allocated(micro%aggrselfprisr)) deallocate (micro%aggrselfprisr)
    if (allocated(micro%accaggrselfpris)) deallocate (micro%accaggrselfpris)
    if (allocated(micro%aggrselfsnowr)) deallocate (micro%aggrselfsnowr)
    if (allocated(micro%accaggrselfsnow)) deallocate (micro%accaggrselfsnow)
    if (allocated(micro%aggrprissnowr)) deallocate (micro%aggrprissnowr)
    if (allocated(micro%accaggrprissnow)) deallocate (micro%accaggrprissnow)
    if (allocated(micro%dust1cldrr)) deallocate (micro%dust1cldrr)
    if (allocated(micro%accdust1cldr)) deallocate (micro%accdust1cldr)
    if (allocated(micro%dust2cldrr)) deallocate (micro%dust2cldrr)
    if (allocated(micro%accdust2cldr)) deallocate (micro%accdust2cldr)
    if (allocated(micro%dust1drzrr)) deallocate (micro%dust1drzrr)
    if (allocated(micro%accdust1drzr)) deallocate (micro%accdust1drzr)
    if (allocated(micro%dust2drzrr)) deallocate (micro%dust2drzrr)
    if (allocated(micro%accdust2drzr)) deallocate (micro%accdust2drzr)
    if (allocated(micro%regen1cldrr)) deallocate (micro%regen1cldrr)
    if (allocated(micro%accregen1cldr)) deallocate (micro%accregen1cldr)
    if (allocated(micro%regen2cldrr)) deallocate (micro%regen2cldrr)
    if (allocated(micro%accregen2cldr)) deallocate (micro%accregen2cldr)
    if (allocated(micro%regen1drzrr)) deallocate (micro%regen1drzrr)
    if (allocated(micro%accregen1drzr)) deallocate (micro%accregen1drzr)
    if (allocated(micro%regen2drzrr)) deallocate (micro%regen2drzrr)
    if (allocated(micro%accregen2drzr)) deallocate (micro%accregen2drzr)
    if (allocated(micro%nuccldcr)) deallocate (micro%nuccldcr)
    if (allocated(micro%nucicecr)) deallocate (micro%nucicecr)
    if (allocated(micro%inuchomcr)) deallocate (micro%inuchomcr)
    if (allocated(micro%inucifncr)) deallocate (micro%inucifncr)

    !Bin microphysics variables
    if (allocated(micro%pcpvip))       deallocate (micro%pcpvip)
    if (allocated(micro%pcpvic))       deallocate (micro%pcpvic)
    if (allocated(micro%pcpvid))       deallocate (micro%pcpvid)
    if (allocated(micro%pcprip))       deallocate (micro%pcprip)
    if (allocated(micro%pcpric))       deallocate (micro%pcpric)
    if (allocated(micro%pcprid))       deallocate (micro%pcprid)
    if (allocated(micro%accnuccldc))     deallocate (micro%accnuccldc)
    if (allocated(micro%accnucicec))     deallocate (micro%accnucicec)
    if (allocated(micro%accinuchomc))    deallocate (micro%accinuchomc)
    if (allocated(micro%accinucifnc))    deallocate (micro%accinucifnc)
    if (allocated(micro%accpic))       deallocate (micro%accpic)
    if (allocated(micro%accpip))       deallocate (micro%accpip)
    if (allocated(micro%accpid))       deallocate (micro%accpid)
    if (allocated(micro%fncn))         deallocate (micro%fncn)
    if (allocated(micro%ffcd))         deallocate (micro%ffcd)
    if (allocated(micro%ffic))         deallocate (micro%ffic)
    if (allocated(micro%ffip))         deallocate (micro%ffip)
    if (allocated(micro%ffid))         deallocate (micro%ffid)
    if (allocated(micro%ffsn))         deallocate (micro%ffsn)
    if (allocated(micro%ffgl))         deallocate (micro%ffgl)
    if (allocated(micro%ffhl))         deallocate (micro%ffhl)
    if (allocated(micro%ffin))         deallocate (micro%ffin)
    if (allocated(micro%t_old))        deallocate (micro%t_old)
    if (allocated(micro%rv_old))       deallocate (micro%rv_old)

return
END SUBROUTINE dealloc_micro

!##############################################################################
Subroutine filltab_micro (micro,microm,imean,n1,n2,n3,ng)

use var_tables
use micro_prm, only:nkr

implicit none

   type (micro_vars) :: micro,microm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3
   if (allocated(micro%rcp))   &
      CALL vtables2 (micro%rcp(1,1,1),microm%rcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RCP :3:anal:mpti:mpt1')
   if (allocated(micro%rdp))   &
      CALL vtables2 (micro%rdp(1,1,1),microm%rdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RDP :3:anal:mpti:mpt1')
   if (allocated(micro%rrp))   &
      CALL vtables2 (micro%rrp(1,1,1),microm%rrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RRP :3:anal:mpti:mpt1')
   if (allocated(micro%rpp))   &
      CALL vtables2 (micro%rpp(1,1,1),microm%rpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RPP :3:anal:mpti:mpt1')
   if (allocated(micro%rsp))   &
      CALL vtables2 (micro%rsp(1,1,1),microm%rsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RSP :3:anal:mpti:mpt1')
   if (allocated(micro%rap))   &
      CALL vtables2 (micro%rap(1,1,1),microm%rap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAP :3:anal:mpti:mpt1')
   if (allocated(micro%rgp))   &
      CALL vtables2 (micro%rgp(1,1,1),microm%rgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RGP :3:anal:mpti:mpt1')
   if (allocated(micro%rhp))   &
      CALL vtables2 (micro%rhp(1,1,1),microm%rhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RHP :3:anal:mpti:mpt1')
   if (allocated(micro%ccp))   &
      CALL vtables2 (micro%ccp(1,1,1),microm%ccp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CCP :3:anal:mpti:mpt1')
   if (allocated(micro%cdp))   &
      CALL vtables2 (micro%cdp(1,1,1),microm%cdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CDP :3:anal:mpti:mpt1')
   if (allocated(micro%crp))   &
      CALL vtables2 (micro%crp(1,1,1),microm%crp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CRP :3:anal:mpti:mpt1')
   if (allocated(micro%cpp))   &
      CALL vtables2 (micro%cpp(1,1,1),microm%cpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CPP :3:anal:mpti:mpt1')
   if (allocated(micro%csp))   &
      CALL vtables2 (micro%csp(1,1,1),microm%csp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CSP :3:anal:mpti:mpt1')
   if (allocated(micro%cap))   &
      CALL vtables2 (micro%cap(1,1,1),microm%cap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CAP :3:anal:mpti:mpt1')
   if (allocated(micro%cgp))   &
      CALL vtables2 (micro%cgp(1,1,1),microm%cgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CGP :3:anal:mpti:mpt1')
   if (allocated(micro%chp))   &
      CALL vtables2 (micro%chp(1,1,1),microm%chp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CHP :3:anal:mpti:mpt1')
   if (allocated(micro%cifnp)) &
      CALL vtables2 (micro%cifnp(1,1,1),microm%cifnp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CIFNP :3:anal:mpti:mpt1')

   if (allocated(micro%q2))   &
      CALL vtables2 (micro%q2(1,1,1),microm%q2(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q2 :3:anal:mpti:mpt1')
   if (allocated(micro%q6)) &
      CALL vtables2 (micro%q6(1,1,1),microm%q6(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q6 :3:anal:mpti:mpt1')
   if (allocated(micro%q7)) &
      CALL vtables2 (micro%q7(1,1,1),microm%q7(1,1,1)  &
                 ,ng, npts, imean,  &
                 'Q7 :3:anal:mpti:mpt1')

!Aerosol categories mass and number
   if (allocated(micro%cn1np)) &
      CALL vtables2 (micro%cn1np(1,1,1),microm%cn1np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CN1NP :3:anal:mpti:mpt1')
   if (allocated(micro%cn2np)) &
      CALL vtables2 (micro%cn2np(1,1,1),microm%cn2np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CN2NP :3:anal:mpti:mpt1')
   if (allocated(micro%cn1mp)) &
      CALL vtables2 (micro%cn1mp(1,1,1),microm%cn1mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CN1MP :3:anal:mpti:mpt1')
   if (allocated(micro%cn2mp)) &
      CALL vtables2 (micro%cn2mp(1,1,1),microm%cn2mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CN2MP :3:anal:mpti:mpt1')
   if (allocated(micro%md1np)) &
      CALL vtables2 (micro%md1np(1,1,1),microm%md1np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD1NP :3:anal:mpti:mpt1')
   if (allocated(micro%md2np)) &
      CALL vtables2 (micro%md2np(1,1,1),microm%md2np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD2NP :3:anal:mpti:mpt1')
   if (allocated(micro%md1mp)) &
      CALL vtables2 (micro%md1mp(1,1,1),microm%md1mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD1MP :3:anal:mpti:mpt1')
   if (allocated(micro%md2mp)) &
      CALL vtables2 (micro%md2mp(1,1,1),microm%md2mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MD2MP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_film_np)) &
      CALL vtables2 (micro%salt_film_np(1,1,1),microm%salt_film_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_FILM_NP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_jet_np)) &
      CALL vtables2 (micro%salt_jet_np(1,1,1),microm%salt_jet_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_JET_NP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_spum_np)) &
      CALL vtables2 (micro%salt_spum_np(1,1,1),microm%salt_spum_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_SPUM_NP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_film_mp)) &
      CALL vtables2 (micro%salt_film_mp(1,1,1),microm%salt_film_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_FILM_MP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_jet_mp)) &
      CALL vtables2 (micro%salt_jet_mp(1,1,1),microm%salt_jet_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_JET_MP :3:anal:mpti:mpt1')
   if (allocated(micro%salt_spum_mp)) &
      CALL vtables2 (micro%salt_spum_mp(1,1,1),microm%salt_spum_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SALT_SPUM_MP :3:anal:mpti:mpt1')
   if (allocated(micro%abc1np)) &
      CALL vtables2 (micro%abc1np(1,1,1),microm%abc1np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ABC1NP :3:anal:mpti:mpt1')
   if (allocated(micro%abc2np)) &
      CALL vtables2 (micro%abc2np(1,1,1),microm%abc2np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ABC2NP :3:anal:mpti:mpt1')
   if (allocated(micro%abc1mp)) &
      CALL vtables2 (micro%abc1mp(1,1,1),microm%abc1mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ABC1MP :3:anal:mpti:mpt1')
   if (allocated(micro%abc2mp)) &
      CALL vtables2 (micro%abc2mp(1,1,1),microm%abc2mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ABC2MP :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero1_np)) &
      CALL vtables2 (micro%regen_aero1_np(1,1,1),microm%regen_aero1_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_NP :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero1_nr)) &
      CALL vtables2 (micro%regen_aero1_nr(1,1,1),microm%regen_aero1_nr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_NR :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero1_mr)) &
      CALL vtables2 (micro%regen_aero1_mr(1,1,1),microm%regen_aero1_mr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_MR :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero2_nr)) &
      CALL vtables2 (micro%regen_aero2_nr(1,1,1),microm%regen_aero2_nr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_NR :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero2_mr)) &
      CALL vtables2 (micro%regen_aero2_mr(1,1,1),microm%regen_aero2_mr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_MR :3:anal:mpti:mpt1')
   if (allocated(micro%accregen_aero1_n)) &
      CALL vtables2 (micro%accregen_aero1_n(1,1,1),microm%accregen_aero1_n(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN_AERO1_N :3:anal:mpti:mpt1')
   if (allocated(micro%accregen_aero1_m)) &
      CALL vtables2 (micro%accregen_aero1_m(1,1,1),microm%accregen_aero1_m(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN_AERO1_M :3:anal:mpti:mpt1')
   if (allocated(micro%accregen_aero2_n)) &
      CALL vtables2 (micro%accregen_aero2_n(1,1,1),microm%accregen_aero2_n(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN_AERO2_N :3:anal:mpti:mpt1')
   if (allocated(micro%accregen_aero2_m)) &
      CALL vtables2 (micro%accregen_aero2_m(1,1,1),microm%accregen_aero2_m(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN_AERO2_M :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero1_mp)) &
      CALL vtables2 (micro%regen_aero1_mp(1,1,1),microm%regen_aero1_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO1_MP :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero2_np)) &
      CALL vtables2 (micro%regen_aero2_np(1,1,1),microm%regen_aero2_np(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_NP :3:anal:mpti:mpt1')
   if (allocated(micro%regen_aero2_mp)) &
      CALL vtables2 (micro%regen_aero2_mp(1,1,1),microm%regen_aero2_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN_AERO2_MP :3:anal:mpti:mpt1')

!Immersion freezing nuclei number tracking variables
   if (allocated(micro%immercp)) &
      CALL vtables2 (micro%immercp(1,1,1),microm%immercp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERCP :3:anal:mpti:mpt1')
   if (allocated(micro%immerdp)) &
      CALL vtables2 (micro%immerdp(1,1,1),microm%immerdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERDP :3:anal:mpti:mpt1')
   if (allocated(micro%immerrp)) &
      CALL vtables2 (micro%immerrp(1,1,1),microm%immerrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IMMERRP :3:anal:mpti:mpt1')
   if (allocated(micro%ifnnucp)) &
      CALL vtables2 (micro%ifnnucp(1,1,1),microm%ifnnucp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'IFNNUCP :3:anal:mpti:mpt1')

!Total aerosol mass-in-hydrometeors tracking variables
   if (allocated(micro%cnmcp)) &
      CALL vtables2 (micro%cnmcp(1,1,1),microm%cnmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMCP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmdp)) &
      CALL vtables2 (micro%cnmdp(1,1,1),microm%cnmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMDP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmrp)) &
      CALL vtables2 (micro%cnmrp(1,1,1),microm%cnmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMRP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmpp)) &
      CALL vtables2 (micro%cnmpp(1,1,1),microm%cnmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMPP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmsp)) &
      CALL vtables2 (micro%cnmsp(1,1,1),microm%cnmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMSP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmap)) &
      CALL vtables2 (micro%cnmap(1,1,1),microm%cnmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMAP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmgp)) &
      CALL vtables2 (micro%cnmgp(1,1,1),microm%cnmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMGP :3:anal:mpti:mpt1')
   if (allocated(micro%cnmhp)) &
      CALL vtables2 (micro%cnmhp(1,1,1),microm%cnmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CNMHP :3:anal:mpti:mpt1')

!Dust aerosol mass-in-hydrometeors tracking variables
   if (allocated(micro%dnmcp)) &
      CALL vtables2 (micro%dnmcp(1,1,1),microm%dnmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMCP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmdp)) &
      CALL vtables2 (micro%dnmdp(1,1,1),microm%dnmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMDP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmrp)) &
      CALL vtables2 (micro%dnmrp(1,1,1),microm%dnmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMRP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmpp)) &
      CALL vtables2 (micro%dnmpp(1,1,1),microm%dnmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMPP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmsp)) &
      CALL vtables2 (micro%dnmsp(1,1,1),microm%dnmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMSP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmap)) &
      CALL vtables2 (micro%dnmap(1,1,1),microm%dnmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMAP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmgp)) &
      CALL vtables2 (micro%dnmgp(1,1,1),microm%dnmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMGP :3:anal:mpti:mpt1')
   if (allocated(micro%dnmhp)) &
      CALL vtables2 (micro%dnmhp(1,1,1),microm%dnmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DNMHP :3:anal:mpti:mpt1')

!Dust as Ice nuclei mass-in-hydrometeors tracking variables
   if (allocated(micro%dincp)) &
      CALL vtables2 (micro%dincp(1,1,1),microm%dincp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINCP :3:anal:mpti:mpt1')
   if (allocated(micro%dindp)) &
      CALL vtables2 (micro%dindp(1,1,1),microm%dindp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINDP :3:anal:mpti:mpt1')
   if (allocated(micro%dinrp)) &
      CALL vtables2 (micro%dinrp(1,1,1),microm%dinrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINRP :3:anal:mpti:mpt1')
   if (allocated(micro%dinpp)) &
      CALL vtables2 (micro%dinpp(1,1,1),microm%dinpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINPP :3:anal:mpti:mpt1')
   if (allocated(micro%dinsp)) &
      CALL vtables2 (micro%dinsp(1,1,1),microm%dinsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINSP :3:anal:mpti:mpt1')
   if (allocated(micro%dinap)) &
      CALL vtables2 (micro%dinap(1,1,1),microm%dinap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINAP :3:anal:mpti:mpt1')
   if (allocated(micro%dingp)) &
      CALL vtables2 (micro%dingp(1,1,1),microm%dingp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINGP :3:anal:mpti:mpt1')
   if (allocated(micro%dinhp)) &
      CALL vtables2 (micro%dinhp(1,1,1),microm%dinhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DINHP :3:anal:mpti:mpt1')

!Regenerated aerosol mass-in-hydrometeors tracking variables
   if (allocated(micro%rnmcp)) &
      CALL vtables2 (micro%rnmcp(1,1,1),microm%rnmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMCP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmdp)) &
      CALL vtables2 (micro%rnmdp(1,1,1),microm%rnmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMDP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmrp)) &
      CALL vtables2 (micro%rnmrp(1,1,1),microm%rnmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMRP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmpp)) &
      CALL vtables2 (micro%rnmpp(1,1,1),microm%rnmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMPP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmsp)) &
      CALL vtables2 (micro%rnmsp(1,1,1),microm%rnmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMSP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmap)) &
      CALL vtables2 (micro%rnmap(1,1,1),microm%rnmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMAP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmgp)) &
      CALL vtables2 (micro%rnmgp(1,1,1),microm%rnmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMGP :3:anal:mpti:mpt1')
   if (allocated(micro%rnmhp)) &
      CALL vtables2 (micro%rnmhp(1,1,1),microm%rnmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RNMHP :3:anal:mpti:mpt1')

!Regenerated aerosol as Ice nuclei mass-in-hydrometeors tracking variables
   if (allocated(micro%rincp)) &
      CALL vtables2 (micro%rincp(1,1,1),microm%rincp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINCP :3:anal:mpti:mpt1')
   if (allocated(micro%rindp)) &
      CALL vtables2 (micro%rindp(1,1,1),microm%rindp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINDP :3:anal:mpti:mpt1')
   if (allocated(micro%rinrp)) &
      CALL vtables2 (micro%rinrp(1,1,1),microm%rinrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINRP :3:anal:mpti:mpt1')
   if (allocated(micro%rinpp)) &
      CALL vtables2 (micro%rinpp(1,1,1),microm%rinpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINPP :3:anal:mpti:mpt1')
   if (allocated(micro%rinsp)) &
      CALL vtables2 (micro%rinsp(1,1,1),microm%rinsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINSP :3:anal:mpti:mpt1')
   if (allocated(micro%rinap)) &
      CALL vtables2 (micro%rinap(1,1,1),microm%rinap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINAP :3:anal:mpti:mpt1')
   if (allocated(micro%ringp)) &
      CALL vtables2 (micro%ringp(1,1,1),microm%ringp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINGP :3:anal:mpti:mpt1')
   if (allocated(micro%rinhp)) &
      CALL vtables2 (micro%rinhp(1,1,1),microm%rinhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RINHP :3:anal:mpti:mpt1')

!Total aerosol soluble mass-in-hydrometeors for aerosol tracking
   if (allocated(micro%snmcp)) &
      CALL vtables2 (micro%snmcp(1,1,1),microm%snmcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMCP :3:anal:mpti:mpt1')
   if (allocated(micro%snmdp)) &
      CALL vtables2 (micro%snmdp(1,1,1),microm%snmdp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMDP :3:anal:mpti:mpt1')
   if (allocated(micro%snmrp)) &
      CALL vtables2 (micro%snmrp(1,1,1),microm%snmrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMRP :3:anal:mpti:mpt1')
   if (allocated(micro%snmpp)) &
      CALL vtables2 (micro%snmpp(1,1,1),microm%snmpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMPP :3:anal:mpti:mpt1')
   if (allocated(micro%snmsp)) &
      CALL vtables2 (micro%snmsp(1,1,1),microm%snmsp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMSP :3:anal:mpti:mpt1')
   if (allocated(micro%snmap)) &
      CALL vtables2 (micro%snmap(1,1,1),microm%snmap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMAP :3:anal:mpti:mpt1')
   if (allocated(micro%snmgp)) &
      CALL vtables2 (micro%snmgp(1,1,1),microm%snmgp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMGP :3:anal:mpti:mpt1')
   if (allocated(micro%snmhp)) &
      CALL vtables2 (micro%snmhp(1,1,1),microm%snmhp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SNMHP :3:anal:mpti:mpt1')
   if (allocated(micro%resol_aero1_mp)) &
      CALL vtables2 (micro%resol_aero1_mp(1,1,1),microm%resol_aero1_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RESOL_AERO1_MP :3:anal:mpti:mpt1')
   if (allocated(micro%resol_aero2_mp)) &
      CALL vtables2 (micro%resol_aero2_mp(1,1,1),microm%resol_aero2_mp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RESOL_AERO2_MP :3:anal:mpti:mpt1')

!Vertical precipitation rates
   if (allocated(micro%pcpvr)) &
      CALL vtables2 (micro%pcpvr(1,1,1),microm%pcpvr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVR :3:anal:mpti')
   if (allocated(micro%pcpvp)) &
      CALL vtables2 (micro%pcpvp(1,1,1),microm%pcpvp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVP :3:anal:mpti')
   if (allocated(micro%pcpvs)) &
      CALL vtables2 (micro%pcpvs(1,1,1),microm%pcpvs(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVS :3:anal:mpti')
   if (allocated(micro%pcpva)) &
      CALL vtables2 (micro%pcpva(1,1,1),microm%pcpva(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVA :3:anal:mpti')
   if (allocated(micro%pcpvg)) &
      CALL vtables2 (micro%pcpvg(1,1,1),microm%pcpvg(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVG :3:anal:mpti')
   if (allocated(micro%pcpvh)) &
      CALL vtables2 (micro%pcpvh(1,1,1),microm%pcpvh(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVH :3:anal:mpti')
   if (allocated(micro%pcpvd)) &
      CALL vtables2 (micro%pcpvd(1,1,1),microm%pcpvd(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVD :3:anal:mpti')

   !MICRO BUDGET PROCESS RATES (<NAME>R, per second) AND RUN TOTALS (ACC<NAME>)
   if (allocated(micro%latheatvapr)) &
      CALL vtables2 (micro%latheatvapr(1,1,1),microm%latheatvapr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATVAPR :3:anal:mpti')
   if (allocated(micro%acclatheatvap)) &
      CALL vtables2 (micro%acclatheatvap(1,1,1),microm%acclatheatvap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCLATHEATVAP :3:anal:mpti')
   if (allocated(micro%latheatfrzr)) &
      CALL vtables2 (micro%latheatfrzr(1,1,1),microm%latheatfrzr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LATHEATFRZR :3:anal:mpti')
   if (allocated(micro%acclatheatfrz)) &
      CALL vtables2 (micro%acclatheatfrz(1,1,1),microm%acclatheatfrz(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCLATHEATFRZ :3:anal:mpti')
   if (allocated(micro%nuccldrr)) &
      CALL vtables2 (micro%nuccldrr(1,1,1),microm%nuccldrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCCLDRR :3:anal:mpti')
   if (allocated(micro%accnuccldr)) &
      CALL vtables2 (micro%accnuccldr(1,1,1),microm%accnuccldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCNUCCLDR :3:anal:mpti')
   if (allocated(micro%cld2rainr)) &
      CALL vtables2 (micro%cld2rainr(1,1,1),microm%cld2rainr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CLD2RAINR :3:anal:mpti')
   if (allocated(micro%acccld2rain)) &
      CALL vtables2 (micro%acccld2rain(1,1,1),microm%acccld2rain(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCCLD2RAIN :3:anal:mpti')
   if (allocated(micro%ice2rainr)) &
      CALL vtables2 (micro%ice2rainr(1,1,1),microm%ice2rainr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ICE2RAINR :3:anal:mpti')
   if (allocated(micro%accice2rain)) &
      CALL vtables2 (micro%accice2rain(1,1,1),microm%accice2rain(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCICE2RAIN :3:anal:mpti')
   if (allocated(micro%nucicerr)) &
      CALL vtables2 (micro%nucicerr(1,1,1),microm%nucicerr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCICERR :3:anal:mpti')
   if (allocated(micro%accnucicer)) &
      CALL vtables2 (micro%accnucicer(1,1,1),microm%accnucicer(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCNUCICER :3:anal:mpti')
   if (allocated(micro%vapliqr)) &
      CALL vtables2 (micro%vapliqr(1,1,1),microm%vapliqr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPLIQR :3:anal:mpti')
   if (allocated(micro%accvapliq)) &
      CALL vtables2 (micro%accvapliq(1,1,1),microm%accvapliq(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPLIQ :3:anal:mpti')
   if (allocated(micro%vapicer)) &
      CALL vtables2 (micro%vapicer(1,1,1),microm%vapicer(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPICER :3:anal:mpti')
   if (allocated(micro%accvapice)) &
      CALL vtables2 (micro%accvapice(1,1,1),microm%accvapice(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPICE :3:anal:mpti')
   if (allocated(micro%evapliqr)) &
      CALL vtables2 (micro%evapliqr(1,1,1),microm%evapliqr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPLIQR :3:anal:mpti')
   if (allocated(micro%accevapliq)) &
      CALL vtables2 (micro%accevapliq(1,1,1),microm%accevapliq(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPLIQ :3:anal:mpti')
   if (allocated(micro%evapicer)) &
      CALL vtables2 (micro%evapicer(1,1,1),microm%evapicer(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPICER :3:anal:mpti')
   if (allocated(micro%accevapice)) &
      CALL vtables2 (micro%accevapice(1,1,1),microm%accevapice(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPICE :3:anal:mpti')
   if (allocated(micro%freezingr)) &
      CALL vtables2 (micro%freezingr(1,1,1),microm%freezingr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'FREEZINGR :3:anal:mpti')
   if (allocated(micro%accfreezing)) &
      CALL vtables2 (micro%accfreezing(1,1,1),microm%accfreezing(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCFREEZING :3:anal:mpti')
   if (allocated(micro%meltingr)) &
      CALL vtables2 (micro%meltingr(1,1,1),microm%meltingr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTINGR :3:anal:mpti')
   if (allocated(micro%accmelting)) &
      CALL vtables2 (micro%accmelting(1,1,1),microm%accmelting(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTING :3:anal:mpti')
   if (allocated(micro%melticer)) &
      CALL vtables2 (micro%melticer(1,1,1),microm%melticer(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTICER :3:anal:mpti')
   if (allocated(micro%accmeltice)) &
      CALL vtables2 (micro%accmeltice(1,1,1),microm%accmeltice(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTICE :3:anal:mpti')
   if (allocated(micro%rimecldr)) &
      CALL vtables2 (micro%rimecldr(1,1,1),microm%rimecldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDR :3:anal:mpti')
   if (allocated(micro%accrimecld)) &
      CALL vtables2 (micro%accrimecld(1,1,1),microm%accrimecld(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRIMECLD :3:anal:mpti')
   if (allocated(micro%rain2icer)) &
      CALL vtables2 (micro%rain2icer(1,1,1),microm%rain2icer(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2ICER :3:anal:mpti')
   if (allocated(micro%accrain2ice)) &
      CALL vtables2 (micro%accrain2ice(1,1,1),microm%accrain2ice(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2ICE :3:anal:mpti')
   if (allocated(micro%aggregater)) &
      CALL vtables2 (micro%aggregater(1,1,1),microm%aggregater(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGREGATER :3:anal:mpti')
   if (allocated(micro%accaggregate)) &
      CALL vtables2 (micro%accaggregate(1,1,1),microm%accaggregate(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCAGGREGATE :3:anal:mpti')
   if (allocated(micro%inuchomrr)) &
      CALL vtables2 (micro%inuchomrr(1,1,1),microm%inuchomrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHOMRR :3:anal:mpti')
   if (allocated(micro%accinuchomr)) &
      CALL vtables2 (micro%accinuchomr(1,1,1),microm%accinuchomr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCHOMR :3:anal:mpti')
   if (allocated(micro%inuccontrr)) &
      CALL vtables2 (micro%inuccontrr(1,1,1),microm%inuccontrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCCONTRR :3:anal:mpti')
   if (allocated(micro%accinuccontr)) &
      CALL vtables2 (micro%accinuccontr(1,1,1),microm%accinuccontr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCCONTR :3:anal:mpti')
   if (allocated(micro%inucifnrr)) &
      CALL vtables2 (micro%inucifnrr(1,1,1),microm%inucifnrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCIFNRR :3:anal:mpti')
   if (allocated(micro%accinucifnr)) &
      CALL vtables2 (micro%accinucifnr(1,1,1),microm%accinucifnr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCIFNR :3:anal:mpti')
   if (allocated(micro%inuchazrr)) &
      CALL vtables2 (micro%inuchazrr(1,1,1),microm%inuchazrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHAZRR :3:anal:mpti')
   if (allocated(micro%accinuchazr)) &
      CALL vtables2 (micro%accinuchazr(1,1,1),microm%accinuchazr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCHAZR :3:anal:mpti')
   if (allocated(micro%vapcldr)) &
      CALL vtables2 (micro%vapcldr(1,1,1),microm%vapcldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPCLDR :3:anal:mpti')
   if (allocated(micro%accvapcld)) &
      CALL vtables2 (micro%accvapcld(1,1,1),microm%accvapcld(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPCLD :3:anal:mpti')
   if (allocated(micro%vaprainr)) &
      CALL vtables2 (micro%vaprainr(1,1,1),microm%vaprainr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPRAINR :3:anal:mpti')
   if (allocated(micro%accvaprain)) &
      CALL vtables2 (micro%accvaprain(1,1,1),microm%accvaprain(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPRAIN :3:anal:mpti')
   if (allocated(micro%vapprisr)) &
      CALL vtables2 (micro%vapprisr(1,1,1),microm%vapprisr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPPRISR :3:anal:mpti')
   if (allocated(micro%accvappris)) &
      CALL vtables2 (micro%accvappris(1,1,1),microm%accvappris(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPPRIS :3:anal:mpti')
   if (allocated(micro%vapsnowr)) &
      CALL vtables2 (micro%vapsnowr(1,1,1),microm%vapsnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPSNOWR :3:anal:mpti')
   if (allocated(micro%accvapsnow)) &
      CALL vtables2 (micro%accvapsnow(1,1,1),microm%accvapsnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPSNOW :3:anal:mpti')
   if (allocated(micro%vapaggrr)) &
      CALL vtables2 (micro%vapaggrr(1,1,1),microm%vapaggrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPAGGRR :3:anal:mpti')
   if (allocated(micro%accvapaggr)) &
      CALL vtables2 (micro%accvapaggr(1,1,1),microm%accvapaggr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPAGGR :3:anal:mpti')
   if (allocated(micro%vapgraur)) &
      CALL vtables2 (micro%vapgraur(1,1,1),microm%vapgraur(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPGRAUR :3:anal:mpti')
   if (allocated(micro%accvapgrau)) &
      CALL vtables2 (micro%accvapgrau(1,1,1),microm%accvapgrau(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPGRAU :3:anal:mpti')
   if (allocated(micro%vaphailr)) &
      CALL vtables2 (micro%vaphailr(1,1,1),microm%vaphailr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPHAILR :3:anal:mpti')
   if (allocated(micro%accvaphail)) &
      CALL vtables2 (micro%accvaphail(1,1,1),microm%accvaphail(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPHAIL :3:anal:mpti')
   if (allocated(micro%vapdrizr)) &
      CALL vtables2 (micro%vapdrizr(1,1,1),microm%vapdrizr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VAPDRIZR :3:anal:mpti')
   if (allocated(micro%accvapdriz)) &
      CALL vtables2 (micro%accvapdriz(1,1,1),microm%accvapdriz(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCVAPDRIZ :3:anal:mpti')
   if (allocated(micro%evapcldr)) &
      CALL vtables2 (micro%evapcldr(1,1,1),microm%evapcldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPCLDR :3:anal:mpti')
   if (allocated(micro%accevapcld)) &
      CALL vtables2 (micro%accevapcld(1,1,1),microm%accevapcld(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPCLD :3:anal:mpti')
   if (allocated(micro%evaprainr)) &
      CALL vtables2 (micro%evaprainr(1,1,1),microm%evaprainr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPRAINR :3:anal:mpti')
   if (allocated(micro%accevaprain)) &
      CALL vtables2 (micro%accevaprain(1,1,1),microm%accevaprain(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPRAIN :3:anal:mpti')
   if (allocated(micro%evapprisr)) &
      CALL vtables2 (micro%evapprisr(1,1,1),microm%evapprisr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPPRISR :3:anal:mpti')
   if (allocated(micro%accevappris)) &
      CALL vtables2 (micro%accevappris(1,1,1),microm%accevappris(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPPRIS :3:anal:mpti')
   if (allocated(micro%evapsnowr)) &
      CALL vtables2 (micro%evapsnowr(1,1,1),microm%evapsnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPSNOWR :3:anal:mpti')
   if (allocated(micro%accevapsnow)) &
      CALL vtables2 (micro%accevapsnow(1,1,1),microm%accevapsnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPSNOW :3:anal:mpti')
   if (allocated(micro%evapaggrr)) &
      CALL vtables2 (micro%evapaggrr(1,1,1),microm%evapaggrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPAGGRR :3:anal:mpti')
   if (allocated(micro%accevapaggr)) &
      CALL vtables2 (micro%accevapaggr(1,1,1),microm%accevapaggr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPAGGR :3:anal:mpti')
   if (allocated(micro%evapgraur)) &
      CALL vtables2 (micro%evapgraur(1,1,1),microm%evapgraur(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPGRAUR :3:anal:mpti')
   if (allocated(micro%accevapgrau)) &
      CALL vtables2 (micro%accevapgrau(1,1,1),microm%accevapgrau(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPGRAU :3:anal:mpti')
   if (allocated(micro%evaphailr)) &
      CALL vtables2 (micro%evaphailr(1,1,1),microm%evaphailr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPHAILR :3:anal:mpti')
   if (allocated(micro%accevaphail)) &
      CALL vtables2 (micro%accevaphail(1,1,1),microm%accevaphail(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPHAIL :3:anal:mpti')
   if (allocated(micro%evapdrizr)) &
      CALL vtables2 (micro%evapdrizr(1,1,1),microm%evapdrizr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'EVAPDRIZR :3:anal:mpti')
   if (allocated(micro%accevapdriz)) &
      CALL vtables2 (micro%accevapdriz(1,1,1),microm%accevapdriz(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCEVAPDRIZ :3:anal:mpti')
   if (allocated(micro%meltprisr)) &
      CALL vtables2 (micro%meltprisr(1,1,1),microm%meltprisr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTPRISR :3:anal:mpti')
   if (allocated(micro%accmeltpris)) &
      CALL vtables2 (micro%accmeltpris(1,1,1),microm%accmeltpris(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTPRIS :3:anal:mpti')
   if (allocated(micro%meltsnowr)) &
      CALL vtables2 (micro%meltsnowr(1,1,1),microm%meltsnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTSNOWR :3:anal:mpti')
   if (allocated(micro%accmeltsnow)) &
      CALL vtables2 (micro%accmeltsnow(1,1,1),microm%accmeltsnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTSNOW :3:anal:mpti')
   if (allocated(micro%meltaggrr)) &
      CALL vtables2 (micro%meltaggrr(1,1,1),microm%meltaggrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTAGGRR :3:anal:mpti')
   if (allocated(micro%accmeltaggr)) &
      CALL vtables2 (micro%accmeltaggr(1,1,1),microm%accmeltaggr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTAGGR :3:anal:mpti')
   if (allocated(micro%meltgraur)) &
      CALL vtables2 (micro%meltgraur(1,1,1),microm%meltgraur(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTGRAUR :3:anal:mpti')
   if (allocated(micro%accmeltgrau)) &
      CALL vtables2 (micro%accmeltgrau(1,1,1),microm%accmeltgrau(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTGRAU :3:anal:mpti')
   if (allocated(micro%melthailr)) &
      CALL vtables2 (micro%melthailr(1,1,1),microm%melthailr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'MELTHAILR :3:anal:mpti')
   if (allocated(micro%accmelthail)) &
      CALL vtables2 (micro%accmelthail(1,1,1),microm%accmelthail(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCMELTHAIL :3:anal:mpti')
   if (allocated(micro%rimecldsnowr)) &
      CALL vtables2 (micro%rimecldsnowr(1,1,1),microm%rimecldsnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDSNOWR :3:anal:mpti')
   if (allocated(micro%accrimecldsnow)) &
      CALL vtables2 (micro%accrimecldsnow(1,1,1),microm%accrimecldsnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRIMECLDSNOW :3:anal:mpti')
   if (allocated(micro%rimecldaggrr)) &
      CALL vtables2 (micro%rimecldaggrr(1,1,1),microm%rimecldaggrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDAGGRR :3:anal:mpti')
   if (allocated(micro%accrimecldaggr)) &
      CALL vtables2 (micro%accrimecldaggr(1,1,1),microm%accrimecldaggr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRIMECLDAGGR :3:anal:mpti')
   if (allocated(micro%rimecldgraur)) &
      CALL vtables2 (micro%rimecldgraur(1,1,1),microm%rimecldgraur(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDGRAUR :3:anal:mpti')
   if (allocated(micro%accrimecldgrau)) &
      CALL vtables2 (micro%accrimecldgrau(1,1,1),microm%accrimecldgrau(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRIMECLDGRAU :3:anal:mpti')
   if (allocated(micro%rimecldhailr)) &
      CALL vtables2 (micro%rimecldhailr(1,1,1),microm%rimecldhailr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RIMECLDHAILR :3:anal:mpti')
   if (allocated(micro%accrimecldhail)) &
      CALL vtables2 (micro%accrimecldhail(1,1,1),microm%accrimecldhail(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRIMECLDHAIL :3:anal:mpti')
   if (allocated(micro%rain2prr)) &
      CALL vtables2 (micro%rain2prr(1,1,1),microm%rain2prr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2PRR :3:anal:mpti')
   if (allocated(micro%accrain2pr)) &
      CALL vtables2 (micro%accrain2pr(1,1,1),microm%accrain2pr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2PR :3:anal:mpti')
   if (allocated(micro%rain2snr)) &
      CALL vtables2 (micro%rain2snr(1,1,1),microm%rain2snr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2SNR :3:anal:mpti')
   if (allocated(micro%accrain2sn)) &
      CALL vtables2 (micro%accrain2sn(1,1,1),microm%accrain2sn(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2SN :3:anal:mpti')
   if (allocated(micro%rain2agr)) &
      CALL vtables2 (micro%rain2agr(1,1,1),microm%rain2agr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2AGR :3:anal:mpti')
   if (allocated(micro%accrain2ag)) &
      CALL vtables2 (micro%accrain2ag(1,1,1),microm%accrain2ag(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2AG :3:anal:mpti')
   if (allocated(micro%rain2grr)) &
      CALL vtables2 (micro%rain2grr(1,1,1),microm%rain2grr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2GRR :3:anal:mpti')
   if (allocated(micro%accrain2gr)) &
      CALL vtables2 (micro%accrain2gr(1,1,1),microm%accrain2gr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2GR :3:anal:mpti')
   if (allocated(micro%rain2har)) &
      CALL vtables2 (micro%rain2har(1,1,1),microm%rain2har(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RAIN2HAR :3:anal:mpti')
   if (allocated(micro%accrain2ha)) &
      CALL vtables2 (micro%accrain2ha(1,1,1),microm%accrain2ha(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCRAIN2HA :3:anal:mpti')
   if (allocated(micro%aggrselfprisr)) &
      CALL vtables2 (micro%aggrselfprisr(1,1,1),microm%aggrselfprisr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRSELFPRISR :3:anal:mpti')
   if (allocated(micro%accaggrselfpris)) &
      CALL vtables2 (micro%accaggrselfpris(1,1,1),microm%accaggrselfpris(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCAGGRSELFPRIS :3:anal:mpti')
   if (allocated(micro%aggrselfsnowr)) &
      CALL vtables2 (micro%aggrselfsnowr(1,1,1),microm%aggrselfsnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRSELFSNOWR :3:anal:mpti')
   if (allocated(micro%accaggrselfsnow)) &
      CALL vtables2 (micro%accaggrselfsnow(1,1,1),microm%accaggrselfsnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCAGGRSELFSNOW :3:anal:mpti')
   if (allocated(micro%aggrprissnowr)) &
      CALL vtables2 (micro%aggrprissnowr(1,1,1),microm%aggrprissnowr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'AGGRPRISSNOWR :3:anal:mpti')
   if (allocated(micro%accaggrprissnow)) &
      CALL vtables2 (micro%accaggrprissnow(1,1,1),microm%accaggrprissnow(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCAGGRPRISSNOW :3:anal:mpti')
   if (allocated(micro%dust1cldrr)) &
      CALL vtables2 (micro%dust1cldrr(1,1,1),microm%dust1cldrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DUST1CLDRR :3:anal:mpti')
   if (allocated(micro%accdust1cldr)) &
      CALL vtables2 (micro%accdust1cldr(1,1,1),microm%accdust1cldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCDUST1CLDR :3:anal:mpti')
   if (allocated(micro%dust2cldrr)) &
      CALL vtables2 (micro%dust2cldrr(1,1,1),microm%dust2cldrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DUST2CLDRR :3:anal:mpti')
   if (allocated(micro%accdust2cldr)) &
      CALL vtables2 (micro%accdust2cldr(1,1,1),microm%accdust2cldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCDUST2CLDR :3:anal:mpti')
   if (allocated(micro%dust1drzrr)) &
      CALL vtables2 (micro%dust1drzrr(1,1,1),microm%dust1drzrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DUST1DRZRR :3:anal:mpti')
   if (allocated(micro%accdust1drzr)) &
      CALL vtables2 (micro%accdust1drzr(1,1,1),microm%accdust1drzr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCDUST1DRZR :3:anal:mpti')
   if (allocated(micro%dust2drzrr)) &
      CALL vtables2 (micro%dust2drzrr(1,1,1),microm%dust2drzrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DUST2DRZRR :3:anal:mpti')
   if (allocated(micro%accdust2drzr)) &
      CALL vtables2 (micro%accdust2drzr(1,1,1),microm%accdust2drzr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCDUST2DRZR :3:anal:mpti')
   if (allocated(micro%regen1cldrr)) &
      CALL vtables2 (micro%regen1cldrr(1,1,1),microm%regen1cldrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN1CLDRR :3:anal:mpti')
   if (allocated(micro%accregen1cldr)) &
      CALL vtables2 (micro%accregen1cldr(1,1,1),microm%accregen1cldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN1CLDR :3:anal:mpti')
   if (allocated(micro%regen2cldrr)) &
      CALL vtables2 (micro%regen2cldrr(1,1,1),microm%regen2cldrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN2CLDRR :3:anal:mpti')
   if (allocated(micro%accregen2cldr)) &
      CALL vtables2 (micro%accregen2cldr(1,1,1),microm%accregen2cldr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN2CLDR :3:anal:mpti')
   if (allocated(micro%regen1drzrr)) &
      CALL vtables2 (micro%regen1drzrr(1,1,1),microm%regen1drzrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN1DRZRR :3:anal:mpti')
   if (allocated(micro%accregen1drzr)) &
      CALL vtables2 (micro%accregen1drzr(1,1,1),microm%accregen1drzr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN1DRZR :3:anal:mpti')
   if (allocated(micro%regen2drzrr)) &
      CALL vtables2 (micro%regen2drzrr(1,1,1),microm%regen2drzrr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'REGEN2DRZRR :3:anal:mpti')
   if (allocated(micro%accregen2drzr)) &
      CALL vtables2 (micro%accregen2drzr(1,1,1),microm%accregen2drzr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCREGEN2DRZR :3:anal:mpti')

   npts=n2*n3
   if (allocated(micro%accpr)) &
      CALL vtables2 (micro%accpr(1,1),microm%accpr(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPR :2:anal:mpti')
   if (allocated(micro%accpp)) &
      CALL vtables2 (micro%accpp(1,1),microm%accpp(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPP :2:anal:mpti')
   if (allocated(micro%accps)) &
      CALL vtables2 (micro%accps(1,1),microm%accps(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPS :2:anal:mpti')
   if (allocated(micro%accpa)) &
      CALL vtables2 (micro%accpa(1,1),microm%accpa(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPA :2:anal:mpti')
   if (allocated(micro%accpg)) &
      CALL vtables2 (micro%accpg(1,1),microm%accpg(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPG :2:anal:mpti')
   if (allocated(micro%accph)) &
      CALL vtables2 (micro%accph(1,1),microm%accph(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPH :2:anal:mpti')
   if (allocated(micro%accpd)) &
      CALL vtables2 (micro%accpd(1,1),microm%accpd(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPD :2:anal:mpti')
   if (allocated(micro%pcprr)) &
      CALL vtables2 (micro%pcprr(1,1),microm%pcprr(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRR :2:anal')
   if (allocated(micro%pcprp)) &
      CALL vtables2 (micro%pcprp(1,1),microm%pcprp(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRP :2:anal')
   if (allocated(micro%pcprs)) &
      CALL vtables2 (micro%pcprs(1,1),microm%pcprs(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRS :2:anal')
   if (allocated(micro%pcpra)) &
      CALL vtables2 (micro%pcpra(1,1),microm%pcpra(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRA :2:anal')
   if (allocated(micro%pcprg)) &
      CALL vtables2 (micro%pcprg(1,1),microm%pcprg(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRG :2:anal')
   if (allocated(micro%pcprh)) &
      CALL vtables2 (micro%pcprh(1,1),microm%pcprh(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRH :2:anal')
   if (allocated(micro%pcprd)) &
      CALL vtables2 (micro%pcprd(1,1),microm%pcprd(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRD :2:anal')
   if (allocated(micro%pcpg)) &
      CALL vtables2 (micro%pcpg(1,1),microm%pcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPG :2:anal:mpti')
   if (allocated(micro%qpcpg)) &
      CALL vtables2 (micro%qpcpg(1,1),microm%qpcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'QPCPG :2:anal:mpti')
   if (allocated(micro%dpcpg)) &
      CALL vtables2 (micro%dpcpg(1,1),microm%dpcpg(1,1)  &
                 ,ng, npts, imean,  &
                 'DPCPG :2:anal:mpti')

!Total aerosol mass precipitation to the surface
   if (allocated(micro%accpaero)) &
      CALL vtables2 (micro%accpaero(1,1),microm%accpaero(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPAERO :2:anal:mpti')
   if (allocated(micro%pcpraero)) &
      CALL vtables2 (micro%pcpraero(1,1),microm%pcpraero(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRAERO :2:anal')

!Dust mass precipitation to the surface
   if (allocated(micro%accpdust)) &
      CALL vtables2 (micro%accpdust(1,1),microm%accpdust(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPDUST :2:anal:mpti')
   if (allocated(micro%pcprdust)) &
      CALL vtables2 (micro%pcprdust(1,1),microm%pcprdust(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRDUST :2:anal')

!Regenerated aerosol mass precipitation to the surface
   if (allocated(micro%accpregen)) &
      CALL vtables2 (micro%accpregen(1,1),microm%accpregen(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPREGEN :2:anal:mpti')
   if (allocated(micro%pcprregen)) &
      CALL vtables2 (micro%pcprregen(1,1),microm%pcprregen(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRREGEN :2:anal')

!Dust erodible fraction
   if (allocated(micro%dustfrac)) &
      CALL vtables2 (micro%dustfrac(1,1),microm%dustfrac(1,1)  &
                 ,ng, npts, imean,  &
                 'DUSTFRAC :2:anal:mpti')

!Bin microphysics variables
   npts=n2*n3
   if (allocated(micro%accpic)) &
      CALL vtables2 (micro%accpic(1,1),microm%accpic(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPIC :2:anal:mpti')
   if (allocated(micro%accpip)) &
      CALL vtables2 (micro%accpip(1,1),microm%accpip(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPIP :2:anal:mpti')
   if (allocated(micro%accpid)) &
      CALL vtables2 (micro%accpid(1,1),microm%accpid(1,1)  &
                 ,ng, npts, imean,  &
                 'ACCPID :2:anal:mpti')
   if (allocated(micro%pcprip)) &
      CALL vtables2 (micro%pcprip(1,1),microm%pcprip(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRIP :2:anal')
   if (allocated(micro%pcprip)) &
      CALL vtables2 (micro%pcpric(1,1),microm%pcpric(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRIC :2:anal')
   if (allocated(micro%pcprid)) &
      CALL vtables2 (micro%pcprid(1,1),microm%pcprid(1,1)  &
                 ,ng, npts, imean,  &
                 'PCPRID :2:anal')

   npts = n1*n2*n3
   if (allocated(micro%pcpvip)) &
      CALL vtables2 (micro%pcpvip(1,1,1),microm%pcpvip(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVIP :3:anal:mpti')
   if (allocated(micro%pcpvic)) &
      CALL vtables2 (micro%pcpvic(1,1,1),microm%pcpvic(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVIC :3:anal:mpti')
   if (allocated(micro%pcpvid)) &
      CALL vtables2 (micro%pcpvid(1,1,1),microm%pcpvid(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PCPVID :3:anal:mpti')
   if (allocated(micro%nuccldcr)) &
      CALL vtables2 (micro%nuccldcr(1,1,1),microm%nuccldcr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCCLDCR :3:anal:mpti')
   if (allocated(micro%accnuccldc)) &
      CALL vtables2 (micro%accnuccldc(1,1,1),microm%accnuccldc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCNUCCLDC :3:anal:mpti')
   if (allocated(micro%nucicecr)) &
      CALL vtables2 (micro%nucicecr(1,1,1),microm%nucicecr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'NUCICECR :3:anal:mpti')
   if (allocated(micro%accnucicec)) &
      CALL vtables2 (micro%accnucicec(1,1,1),microm%accnucicec(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCNUCICEC :3:anal:mpti')
   if (allocated(micro%inuchomcr)) &
      CALL vtables2 (micro%inuchomcr(1,1,1),microm%inuchomcr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCHOMCR :3:anal:mpti')
   if (allocated(micro%accinuchomc)) &
      CALL vtables2 (micro%accinuchomc(1,1,1),microm%accinuchomc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCHOMC :3:anal:mpti')
   if (allocated(micro%inucifncr)) &
      CALL vtables2 (micro%inucifncr(1,1,1),microm%inucifncr(1,1,1)  &
                 ,ng, npts, imean,  &
                 'INUCIFNCR :3:anal:mpti')
   if (allocated(micro%accinucifnc)) &
      CALL vtables2 (micro%accinucifnc(1,1,1),microm%accinucifnc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'ACCINUCIFNC :3:anal:mpti')

   npts=n1*n2*n3*nkr
   if (allocated(micro%fncn)) &
      CALL vtables2 (micro%fncn(1,1,1,1),microm%fncn(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FNCN :7:anal:mpti:mpt1')
   if (allocated(micro%ffcd)) &
      CALL vtables2 (micro%ffcd(1,1,1,1),microm%ffcd(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFCD :7:anal:mpti:mpt1')
   if (allocated(micro%ffic)) &
      CALL vtables2 (micro%ffic(1,1,1,1),microm%ffic(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFIC :7:anal:mpti:mpt1')
   if (allocated(micro%ffip)) &
      CALL vtables2 (micro%ffip(1,1,1,1),microm%ffip(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFIP :7:anal:mpti:mpt1')
   if (allocated(micro%ffid)) &
      CALL vtables2 (micro%ffid(1,1,1,1),microm%ffid(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFID :7:anal:mpti:mpt1')
   if (allocated(micro%ffsn)) &
      CALL vtables2 (micro%ffsn(1,1,1,1),microm%ffsn(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFSN :7:anal:mpti:mpt1')
   if (allocated(micro%ffgl)) &
      CALL vtables2 (micro%ffgl(1,1,1,1),microm%ffgl(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFGL :7:anal:mpti:mpt1')
   if (allocated(micro%ffhl)) &
      CALL vtables2 (micro%ffhl(1,1,1,1),microm%ffhl(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFHL :7:anal:mpti:mpt1')
   if (allocated(micro%ffin)) &
      CALL vtables2 (micro%ffin(1,1,1,1),microm%ffin(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'FFIN :7:anal:mpti:mpt1')
   npts=n1*n2*n3
   if (allocated(micro%t_old)) &
      CALL vtables2 (micro%t_old(1,1,1),microm%t_old(1,1,1)  &
                 ,ng, npts, imean,  &
                 'T_OLD :3:anal:mpti')
   if (allocated(micro%rv_old)) &
      CALL vtables2 (micro%rv_old(1,1,1),microm%rv_old(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RV_OLD :3:anal:mpti')
return
END SUBROUTINE filltab_micro

END MODULE mem_micro
