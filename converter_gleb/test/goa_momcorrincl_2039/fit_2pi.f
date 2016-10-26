
      SUBROUTINE fit_2pi
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      include "select.inc"

      integerj1,j2

        real raddeg, degrad
        data degrad/0.0174533/
        data raddeg/57.2958/

c+ local variables
c++ KINEFIT variables
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      INTEGER*4 ID_PART,IER_KIN,NST_KIN
      REAL*4 PM_KIN,TH_KIN,FI_KIN,
     +      EP_KIN,ET_KIN,EF_KIN,PM_GAM,EP_GAM,
     +      PP_INV,C4_KIN,CHI_KIN,CL_KIN

      logical flag_proton, flag_fit
c--      


c+ Starting KINEFIT for two pions PHOTO production analysis
        if(.not.beam_type) then
           flag_fit = .false.
           flag_proton = .false.
           nfit = 0
c+ part 1 = proton
c+ part 2 = pi+
c+ part 3 = pi-
c++ Input initialisations
        PM_GAM = Egamma*1000.
        EP_GAM = 5. ! Assuming 5 MeV photon E error    
c+++ Particle type
        ID_PART(1) = 1 
        ID_PART(2) = 2 
        ID_PART(3) = 2 
c+++ Proton variables filling
c++++ Unmesured proton
         EP_KIN(1)  = 0.
         ET_KIN(1)  = 0.      
         EF_KIN(1)  = 0.
c++++ Mesured proton
        if(nNU.eq.1.AND.q_NU(1).eq.1.) then 
         flag_proton = .true.
         PM_KIN(1) = P_NU(1)*1e3
         EP_KIN(1) = 0.01*P_NU(1)*1e3 ! Assuming 1% mom error 
         TH_KIN(1) = TH_NU(1)*degrad
         ET_KIN(1) = 0.01          ! Assuming 0.01 rad on theta error
         FI_KIN(1) = PH_NU(1)*degrad
         EF_KIN(1) = 0.2          ! Assuming 0.2 rad on phi error
        endif
c--- Ending proton init
c+++ Pion vars filling
c++++ Unmesured pions
         EP_KIN(2)  = 0. ! pion 1
         ET_KIN(2)  = 0.      
         EF_KIN(2)  = 0.
         EP_KIN(3)  = 0. ! pion 2
         ET_KIN(3)  = 0.      
         EF_KIN(3)  = 0.
c++++ pi+ .and. pi- mesured
        if(npi.eq.2.and.npi_plus.eq.1.AND.npi_minus.eq.1) then ! two pions
           flag_fit = .true.
c++++ pion 1 is always positive (2 is negative)
           if(q_pi(1).gt.0) then
              j1 = 1
              j2 = 2
           else
              j1 = 2
              j2 = 1
           endif   
c++++ pion 1
         PM_KIN(2) = P_pi(j1)*1e3
         EP_KIN(2) = 0.01*P_pi(j1)*1e3 ! Assuming 1% mom error 
         TH_KIN(2) = TH_pi(j1)*degrad
         ET_KIN(2) = 0.01          ! Assuming 0.01 rad on theta error
         FI_KIN(2) = PH_pi(j1)*degrad
         EF_KIN(2) = 0.2          ! Assuming 0.2 rad on phi error
c++++ pion 2
         PM_KIN(3) = P_pi(j2)*1e3
         EP_KIN(3) = 0.01*P_pi(j2)*1e3 ! Assuming 1% mom error 
         TH_KIN(3) = TH_pi(j2)*degrad
         ET_KIN(3) = 0.01          ! Assuming 0.01 rad on theta error
         FI_KIN(3) = PH_pi(j2)*degrad
         EF_KIN(3) = 0.2          ! Assuming 0.2 rad on phi error         
        endif
c++++ pi+ .or. pi- mesured
        if(npi.eq.1.and.nNU.eq.1.AND.q_NU(1).eq.1.) then ! pion and proton
           flag_fit = .true.
c++++ pion 1 is always positive (2 is negative)
           if(q_pi(1).gt.0) then
              j1 = 2
           else
              j1 = 3
           endif   
c++++ pion mesured
         PM_KIN(j1) = P_pi(1)*1e3
         EP_KIN(j1) = 0.01*P_pi(1)*1e3 ! Assuming 1% mom error 
         TH_KIN(j1) = TH_pi(1)*degrad
         ET_KIN(j1) = 0.01          ! Assuming 0.01 rad on theta error
         FI_KIN(j1) = PH_pi(1)*degrad
         EF_KIN(j1) = 0.2          ! Assuming 0.2 rad on phi error
        endif
c--- Ending pion section
c-- Ending variables initialisation

        if (flag_fit) then
         call kinfit(1) ! 1 = proton target
c++ getting results
        nfit = 3
c+++ photon
        fit_Eg(1)   = PM_GAM*1e-3
        fit_DEg(1)  = EP_GAM*1e-3
        fit_Eg(2)  = -1000. ! for the ntuple
        fit_Eg(3)  = -1000.
        fit_DEg(2) = -1000. 
        fit_DEg(3) = -1000. 
c+++ proton
        fit_p (1) = PM_KIN(1)*1e-3
        fit_dp(1) = EP_KIN(1)*1e-3
        fit_th(1) = TH_KIN(1)*raddeg
        fit_dt(1) = ET_KIN(1)*raddeg
        fit_fi(1) = FI_KIN(1)*raddeg
        fit_df(1) = EF_KIN(1)*raddeg
c+++ pi+
        fit_p (2) = PM_KIN(2)*1e-3
        fit_dp(2) = EP_KIN(2)*1e-3
        fit_th(2) = TH_KIN(2)*raddeg
        fit_dt(2) = ET_KIN(2)*raddeg
        fit_fi(2) = FI_KIN(2)*raddeg
        fit_df(2) = EF_KIN(2)*raddeg
c+++ pi-
        fit_p (3) = PM_KIN(3)*1e-3
        fit_dp(3) = EP_KIN(3)*1e-3
        fit_th(3) = TH_KIN(3)*raddeg
        fit_dt(3) = ET_KIN(3)*raddeg
        fit_fi(3) = FI_KIN(3)*raddeg
        fit_df(3) = EF_KIN(3)*raddeg
c+++ Invariant masses
        fit_W(1)  = PP_INV(1)*1e-3 ! 
        fit_W(2)  = PP_INV(2)*1e-3
        fit_W(3)  = PP_INV(3)*1e-3
c+++ Fit status
        fit_ch2(1) = CHI_KIN
        fit_err(1) = IER_KIN
        fit_ch2(2) = -1. ! for the ntuple
        fit_ch2(3) = -1.
        fit_err(2) = -1
        fit_err(3) = -1

c--
        endif ! fit performed 
        endif

      return
      end 


