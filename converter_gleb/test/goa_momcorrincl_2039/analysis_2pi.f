
      SUBROUTINE analysis_2pi(MC_ana)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"  
      include "select.inc"
       integer MC_ana,ntype,sector
       real En(3),Mo(3,3), Be
        real raddeg, degrad
        real acc_missing_part
        data degrad/0.0174533/
        data raddeg/57.2958/
c+ Variables for CM
c+ i index 
c+         0 = Energy of cm-rest particle
c+         1 = px of cm-rest   particle
c+         2 = py of  cm-rest  particle
c+         3 = pz of  cm-rest  particle
c+         4 = sqared mass of cm-rest   particle
c+         5 = momentum of cm-rest  particle
c+         6 = theta of cm-rest  particle
c+         7 = phi of cm-rest  particle
c+ j index runs over possibles combination (e.g. j=1 (p pi1) j=2 (pi1 pi2) j=3 (p pi2) j=4 pi1 for singol meson) 
      REAL*4 p_vec(0:7,3),p_mes(0:7,4),psi(3)
c+ local variables
c++ KINEFIT variables
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      INTEGER*4 ID_PART,IER_KIN,NST_KIN
      REAL*4 PM_KIN,TH_KIN,FI_KIN,
     +      EP_KIN,ET_KIN,EF_KIN,PM_GAM,EP_GAM,
     +      PP_INV,C4_KIN,CHI_KIN,CL_KIN
c--      


c+ Some general setting
c++ Choosing wich entry in the ntuple to fill
c++ Photon beam:
c++ MC only          : 1 <=> MC;  2 <=>fit on MC;
c++ DATA only        : 1 <=> DA;  2 <=>fit on DA; 
c++ MC.and.DATA      : 1 <=> MC;  2 <=>fit on MC; 3 <=> DA;  4 <=>fit on DA;
c++ Electron beam:
c++ MC only          : 1 <=> MC;  
c++ DATA only        : 1 <=> DA;  
c++ MC.and.DATA      : 1 <=> MC; 2 <=> DA 
            ntype = 1 ! By default it fills 1 variable (DATA.or.MC)
            if(.not.beam_type) then ! Photon beam
             if (bit_mc.eq.2) then ! If both are present: 1st = MC 3th = data
               if(MC_ana.eq.0) ntype = 3 
             endif 
            else                   ! Electron beam
             if (bit_mc.eq.2) then ! If both are present: 1st = MC 2nd = data
               if(MC_ana.eq.0) ntype = 2 
             endif 
            endif
c-- Ending Choosing wich entry in the ntuple to fill
c++ Filling general run variables
            indtype = ntype
            call sas_fill_genrun(ntype)
c-- Ending filling general run variables


c++ Filling beam related  variables
            call sas_fill_beamrel(ntype)
c-- Ending filling beam related variables    
c- ending general setting

c+ Starting cmh/rest frame meson and baryons calc
c++ e p -> e' p pi+ + X  Exactly this channel (proton pi+ detected) 
       if((bit2pi.eq.10.or.bit2pi.eq.20).and.W.gt.0) then	

c        call eff_pi_minus(1,point_pm)
            
            ncmrest=1 
c	     if((W.gt.0.93).AND.(W.lt.0.95).AND.(ECout_EL/P_EL.gt.0.25-1.39*ECin_EL/P_EL)
c     %	     .AND.(NpheCC_EL.gt.25)) then
c            write(*,*) 'LAB!!!!!!!!!!!!', 'E0=', E0, 'E_el=', E_el, 'th_el', th_el, 'ph_el=', ph_el,
c     %	                'W=', W, 'Q2=', Q2,
c     %                   'Pmiss_eNp(0,1)=', Pmiss_eNp(0,1), 'Pmiss_eNp(1,1)=', Pmiss_eNp(1,1),
c     %                   'Pmiss_eNp(2,1)=', Pmiss_eNp(2,1), 'Pmiss_eNp(3,1)=', Pmiss_eNp(3,1),
c     %                    'E_PI(1)=', E_PI(1), 'PX_PI(1)=', PX_PI(1), 'PY_PI(1)=', PY_PI(1),
c     %                    'PZ_PI(1)=', PZ_PI(1), 
c     %                   'E_NU(1)=', E_NU(1), 'PX_NU(1)=', PX_NU(1), 'PY_NU(1)=', PY_NU(1),
c     %                    'PZ_NU(1)=', PZ_NU(1),
c     %                   'p_vec=', p_vec, 'p_mes=', p_mes,
c     %                    'psi=',psi   
c                          endif       
            call lab2cms( 5
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eNp(0,1),Pmiss_eNp(1,1),Pmiss_eNp(2,1),Pmiss_eNp(3,1) ! pi- miss
     %                   ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)                          ! pi+ 
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)                          ! p
     %                   ,p_vec,p_mes
     %                   ,psi)
c            if((W.gt.0.93).AND.(W.lt.0.95).AND.(ECout_EL/P_EL.gt.0.25-1.39*ECin_EL/P_EL)
c     %	    .AND.(NpheCC_EL.gt.25).AND.(Pmiss_eNp(0,1).gt.0)) then
c            write(*,*) 'CMS!!!!!!', 'E0=', E0, 'E_el=', E_el, 'th_el', th_el, 'ph_el=', ph_el,
c     %	                'W=', W, 'Q2=', Q2,
c     %                   'Pmiss_eNp(0,1)=', Pmiss_eNp(0,1), 'Pmiss_eNp(1,1)=', Pmiss_eNp(1,1),
c     %                   'Pmiss_eNp(2,1)=', Pmiss_eNp(2,1), 'Pmiss_eNp(3,1)=', Pmiss_eNp(3,1),
c     %                    'E_PI(1)=', E_PI(1), 'PX_PI(1)=', PX_PI(1), 'PY_PI(1)=', PY_PI(1),
c     %                    'PZ_PI(1)=', PZ_PI(1), 
c     %                   'E_NU(1)=', E_NU(1), 'PX_NU(1)=', PX_NU(1), 'PY_NU(1)=', PY_NU(1),
c     %                    'PZ_NU(1)=', PZ_NU(1),
c     %                   'p_vec=', p_vec, 'p_mes=', p_mes,
c     %                    'psi=',psi
                          
c			  endif
               
	       	   
	   
C	   if((W.gt.0.93).AND.(W.lt.0.95).AND.(ECout_EL/P_EL.gt.0.25-1.39*ECin_EL/P_EL)
C     %	    .AND.(NpheCC_EL.gt.25)) then
C            write(*,*) 'CMS!!!!!!', 'E_EL=', E_EL,  'PX_EL=', PX_EL, 'PY_EL=', PY_EL,
C     %	                'PZ_EL', PZ_EL, 'Q2=', Q2,
C     %                   'PMiss_eNp(0)=', PMiss_eNp(0,1), 'PMiss_eNp(1)=', PMiss_eNp(1,1),
C     %                   'PMiss_eNp(2)=', PMiss_eNp(2,1), 'PMiss_eNp(3)=', PMiss_eNp(3,1),
C     %                   'PMiss_eNp(4)=', PMiss_eNp(4,1), 'PMiss_eNp(5)=', PMiss_eNp(5,1),
C     %                   'PMiss_eNp(6)=', PMiss_eNp(6,1), 'PMiss_eNp(7)=', PMiss_eNp(7,1),
C     %                    'E_PI=', E_PI(1), 'PX_PI=', PX_PI(1), 'PY_PI=', PY_PI(1),
C     %                    'PZ_PI=', PZ_PI(1), 
C     %                   'E_NU=', E_NU(1), 'PX_NU(1)=', PX_NU(1), 'PY_NU(1)=', PY_NU(1),
C     %                    'PZ_NU(1)=', PZ_NU(1)
C     %                   
                          
C			  endif
			  
	       
c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(1),P_PI(1),B_PI(1),M2_PI(1)
     &                           ,th_PI(1),ph_PI(1)
     &                           ,acc_PI(1)
     &                           ,x_PI(1),y_PI(1),z_PI(1)
     &                           ,ECtot_PI(1), ECin_PI(1), ECout_PI(1)
     &                           ,PdHit_PI(1) ) 
c+++ pi-
           call pseudo_spa2(-1,-211, Pmiss_eNp(5,1),Pmiss_eNp(6,1),
     &                     (mod(Pmiss_eNp(7,1)+30.,60.)-30.),                           
     &                     I_torus,1.,sector(Pmiss_eNp(7,1)),acc_missing_part)
            
            det_pim(ntype) = -1
            call sas_fill_pim_lab(ntype
     &          ,Pmiss_eNp(0,1),Pmiss_eNp(5,1),Pmiss_eNp(5,1)/Pmiss_eNp(0,1)
     &          ,Pmiss_eNp(4,1),Pmiss_eNp(6,1),Pmiss_eNp(7,1)
     &          ,acc_missing_part
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)

c+++ SC Trigger Threshold 
c           trig_SC(ntype) = max(thrSC_NU(1),thrSC_PI(1))
           id_SC_p(ntype) =  id_SC_NU(1)
           sec_SC_p(ntype) =  sec_SC_NU(1)
           thrSC_p(ntype) = thrSC_NU(1)
           id_SC_pip(ntype) = id_SC_PI(1)
           sec_SC_pip(ntype) = sec_SC_PI(1)
           thrSC_pip(ntype) = thrSC_PI(1)
           id_SC_pim(ntype) =  -1
           sec_SC_pim(ntype) = -1
           thrSC_pim(ntype) =  -1
c+++ eNpp
           call  sas_fill_enpp(ntype)
c+++Rho
            M2_rho(ntype) = p_vec(4,2)
            call sas_fill_rho_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_rho_hcm(ntype, p_vec(0,2) ,p_vec(5,2) ,p_vec(6,2), p_vec(7,2)) 
c+++ Delta++
            M2_dpp(ntype) = p_vec(4,3)
            call sas_fill_dpp_lab(ntype, W_Np(0,1), W_Np(5,1), W_Np(6,1), W_Np(7,1))
            call sas_fill_dpp_hcm(ntype,p_vec(0,3),p_vec(5,3),p_vec(6,3), p_vec(7,3)) 
C+++ Delta0
            M2_d0(ntype) = p_vec(4,1)
            call sas_fill_d0_lab(ntype,Pmiss_ep(0,1),Pmiss_ep(5,1),Pmiss_ep(6,1),Pmiss_ep(7,1))
            call sas_fill_d0_hcm(ntype,   p_vec(0,1) ,  p_vec(5,1) ,  p_vec(6,1),   p_vec(7,1)) 
C+++ Decaying pion
            call sas_fill_pip_rho(ntype,p_mes(0,2),p_mes(5,2),p_mes(6,2),p_mes(7,2),psi(2))
            call sas_fill_pip_dpp(ntype,p_mes(0,3),p_mes(5,3),p_mes(6,3),p_mes(7,3),psi(3))
            call sas_fill_pim_d0 (ntype,p_mes(0,1),p_mes(5,1),p_mes(6,1),p_mes(7,1),psi(1))

            E_lab(1)  = Pmiss_eN(0) ! Vector meson in lab
            P_lab(1)  = Pmiss_eN(5)
            th_lab(1) = Pmiss_eN(6)
            ph_lab(1) = Pmiss_eN(7)

            E_cm(1) = p_vec(0,2) ! Vector meson in cm
            P_cm(1) = p_vec(5,2)
           th_cm(1) = p_vec(6,2)
           ph_cm(1) = p_vec(7,2)

            E_rest(1) = p_mes(0,2)
            P_rest(1) = p_mes(5,2)
           th_rest(1) = p_mes(6,2)
           ph_rest(1) = p_mes(7,2)

	endif

c++ e p -> e' p pi- + X  Exactly this channel (proton pi- detected) 
	if((bit2pi.eq.11.or.bit2pi.eq.21).and.W.gt.0) then
	
c        call eff_pi_minus(1,point_pm)

           ncmrest=1
            call lab2cms( 5
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)                            ! pi-
     %                   ,Pmiss_eNp(0,1),Pmiss_eNp(1,1),Pmiss_eNp(2,1),Pmiss_eNp(3,1)   ! pi+ miss
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)                            ! p 
     %                   ,p_vec ,p_mes
     %                   ,psi)


c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ pi+
           call pseudo_spa2(1,211, Pmiss_eNp(5,1),Pmiss_eNp(6,1),
     &                     (mod(Pmiss_eNp(7,1)+30.,60.)-30.),
     &                     I_torus,1.,sector(Pmiss_eNp(7,1)),acc_missing_part)

            det_pip(ntype) = -1
            call sas_fill_pip_lab(ntype
     &          ,Pmiss_eNp(0,1),Pmiss_eNp(5,1),Pmiss_eNp(5,1)/Pmiss_eNp(0,1)
     &          ,Pmiss_eNp(4,1),Pmiss_eNp(6,1),Pmiss_eNp(7,1)
     &          ,acc_missing_part
     &          ,-1000.,-1000.,1000.
     &          ,-1000.,-1000.,1000.
     &          ,-1000)
c+++ pi-
            det_pim(ntype) = 1
            call sas_fill_pim_lab(ntype
     &          ,E_PI(1),P_PI(1),B_PI(1)
     &          ,M2_PI(1),th_PI(1),ph_PI(1)
     &          ,acc_PI(1)
     &          ,x_PI(1),y_PI(1),z_PI(1)
     &          ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &          ,PdHit_PI(1))
c+++ SC Trigger Threshold 
c           trig_SC(ntype) = max(thrSC_NU(1),thrSC_PI(1))
           id_SC_p(ntype) =  id_SC_NU(1)
           sec_SC_p(ntype) =  sec_SC_NU(1)
           thrSC_p(ntype) = thrSC_NU(1)
           id_SC_pip(ntype) = -1
           sec_SC_pip(ntype) = -1
           thrSC_pip(ntype) = -1
           id_SC_pim(ntype) = id_SC_PI(1)
           sec_SC_pim(ntype) = sec_SC_PI(1)
           thrSC_pim(ntype) = thrSC_PI(1)
c+++ eNpp
           call  sas_fill_enpp(ntype)
c+++Rho
            M2_rho(ntype) = p_vec(4,2)
            call sas_fill_rho_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_rho_hcm(ntype, p_vec(0,2) ,p_vec(5,2) ,p_vec(6,2), p_vec(7,2)) 
c+++ Delta++
            M2_dpp(ntype) = p_vec(4,3)
            call sas_fill_dpp_lab(ntype, W_Np(0,1), W_Np(5,1), W_Np(6,1), W_Np(7,1))
            call sas_fill_dpp_hcm(ntype,p_vec(0,3),p_vec(5,3),p_vec(6,3), p_vec(7,3)) 
C+++ Delta0
            M2_d0(ntype) = p_vec(4,1)
            call sas_fill_d0_lab(ntype,Pmiss_ep(0,1),Pmiss_ep(5,1),Pmiss_ep(6,1),Pmiss_ep(7,1))
            call sas_fill_d0_hcm(ntype,   p_vec(0,1) ,  p_vec(5,1) ,  p_vec(6,1),   p_vec(7,1)) 
C+++ Decaying pion
            call sas_fill_pip_rho(ntype,p_mes(0,2),p_mes(5,2),p_mes(6,2),p_mes(7,2),psi(2))
            call sas_fill_pip_dpp(ntype,p_mes(0,3),p_mes(5,3),p_mes(6,3),p_mes(7,3),psi(3))
            call sas_fill_pim_d0 (ntype,p_mes(0,1),p_mes(5,1),p_mes(6,1),p_mes(7,1),psi(1))


            E_lab(1)  = Pmiss_eN(0) ! Vector meson in lab
            P_lab(1)  = Pmiss_eN(5)
            th_lab(1) = Pmiss_eN(6)
            ph_lab(1) = Pmiss_eN(7)

            E_cm(1) = p_vec(0,2)
            P_cm(1) = p_vec(5,2)
           th_cm(1) = p_vec(6,2)
           ph_cm(1) = p_vec(7,2)

            E_rest(1) = p_mes(0,2)
            P_rest(1) = p_mes(5,2)
           th_rest(1) = p_mes(6,2)
           ph_rest(1) = p_mes(7,2)

	endif
c++ e p -> e' pi+ pi- + X  Exactly this channel (pi+ pi- detected) 
        if((bit2pi.eq.12.or.bit2pi.eq.22).and.W.gt.0) then
c        call eff_pi_minus(1,point_pm)

           ncmrest=1
            call lab2cms( 5
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)                            ! pi-
     %                   ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)                            ! pi+
     %                   ,Pmiss_epp(0,1),Pmiss_epp(1,1),Pmiss_epp(2,1),Pmiss_epp(3,1)   ! proton miss
     %                   ,p_vec ,p_mes
     %                   ,psi)
c+++ proton  
         call pseudo_spa2(1,2212, Pmiss_epp(5,1),Pmiss_epp(6,1),
     &                     (mod(Pmiss_epp(7,1)+30.,60.)-30.),
     &                     I_torus,1.,sector(Pmiss_epp(7,1)),acc_missing_part)
            det_p(ntype) = -1
            call sas_fill_p_lab(ntype
     &      ,Pmiss_epp(0,1),Pmiss_epp(5,1),Pmiss_epp(5,1)/Pmiss_epp(0,1)
     &      ,Pmiss_epp(4,1),Pmiss_epp(6,1),Pmiss_epp(7,1)
     &      ,acc_missing_part
     &      ,-1000.,-1000.,-1000.
     &      ,-1000.,-1000.
     &      ,-1000)
c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(2),P_PI(2),B_PI(2)
     &          ,M2_PI(2),th_PI(2),ph_PI(2)
     &          ,acc_PI(2)
     &          ,x_PI(2),y_PI(2),z_PI(2)
     &          ,ECtot_PI(2), ECin_PI(2), ECout_PI(2)
     &          ,PdHit_PI(2) ) 
c+++ pi-
            det_pim(ntype) = 1
            call sas_fill_pim_lab(ntype,E_PI(1),P_PI(1),B_PI(1)
     &          ,M2_PI(1),th_PI(1),ph_PI(1)
     &          ,acc_PI(1)
     &          ,x_PI(1),y_PI(1),z_PI(1)
     &          ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &          ,PdHit_PI(1))
c+++ SC Trigger Threshold 
c           trig_SC(ntype) = max(thrSC_PI(1),thrSC_PI(2))
           id_SC_p(ntype) =  -1
           sec_SC_p(ntype) =  -1
           thrSC_p(ntype) = -1
           id_SC_pip(ntype) = id_SC_PI(2)
           sec_SC_pip(ntype) = id_SC_PI(2)
           thrSC_pip(ntype) = thrSC_PI(2)
           id_SC_pim(ntype) = id_SC_PI(1)
           sec_SC_pim(ntype) = sec_SC_PI(1)
           thrSC_pim(ntype) = thrSC_PI(1)
c+++ eNpp
           call  sas_fill_enpp(ntype)
c+++Rho
            M2_rho(ntype) = p_vec(4,2)
            call sas_fill_rho_lab(ntype,W_pp(0,1), W_pp(5,1), W_pp(6,1), W_pp(7,1))
            call sas_fill_rho_hcm(ntype, p_vec(0,2) ,p_vec(5,2) ,p_vec(6,2), p_vec(7,2)) 
c+++ Delta++
            M2_dpp(ntype) = p_vec(4,3)
            call sas_fill_dpp_lab(ntype,Pmiss_ep(0,1),Pmiss_ep(5,1), Pmiss_ep(6,1),Pmiss_ep(7,1))
            call sas_fill_dpp_hcm(ntype,   p_vec(0,3),   p_vec(5,3),    p_vec(6,3),   p_vec(7,3)) 
C+++ Delta0
            M2_d0(ntype) = p_vec(4,1)
            call sas_fill_d0_lab(ntype,Pmiss_ep(0,2),Pmiss_ep(5,2), Pmiss_ep(6,2),Pmiss_ep(7,2)) 
            call sas_fill_d0_hcm(ntype,   p_vec(0,1),   p_vec(5,1),    p_vec(6,1),   p_vec(7,1)) 
C+++ Decaying pion
            call sas_fill_pip_rho(ntype,p_mes(0,2),p_mes(5,2),p_mes(6,2),p_mes(7,2),psi(2))
            call sas_fill_pip_dpp(ntype,p_mes(0,3),p_mes(5,3),p_mes(6,3),p_mes(7,3),psi(3))
            call sas_fill_pim_d0 (ntype,p_mes(0,1),p_mes(5,1),p_mes(6,1),p_mes(7,1),psi(1))

            E_lab(1)  = W_pp(0,1) ! Vector meson in lab
            P_lab(1)  = W_pp(5,1)
            th_lab(1) = W_pp(6,1)
            ph_lab(1) = W_pp(7,1)

            E_cm(1) = p_vec(0,2)
            P_cm(1) = p_vec(5,2)
           th_cm(1) = p_vec(6,2)
           ph_cm(1) = p_vec(7,2)

            E_rest(1) = p_mes(0,2)
            P_rest(1) = p_mes(5,2)
           th_rest(1) = p_mes(6,2)
           ph_rest(1) = p_mes(7,2)

	endif
c--
c++ e p -> e' p pi+ pi-  Exactly this channel (proton pi+ pi- detected) 
	if((bit2pi.eq.13.or.bit2pi.eq.23).and.W.gt.0) then	

c        call eff_pi_minus(2,point_pm)
           ncmrest=1
            call lab2cms( 5
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1) ! pi-
     %                   ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2) ! pi+
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1) ! p
     %                   ,p_vec ,p_mes
     %                   ,psi)
c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(2),P_PI(2),B_PI(2)
     &                           ,M2_PI(2),th_PI(2),ph_PI(2)
     &                           ,acc_PI(2)
     &                           ,x_PI(2),y_PI(2),z_PI(2)
     &                           ,ECtot_PI(2), ECin_PI(2), ECout_PI(2) 
     &                           ,PdHit_PI(2)) 
c+++ pi-
            det_pim(ntype) = 1
            call sas_fill_pim_lab(ntype,E_PI(1),P_PI(1),B_PI(1)
     &                           ,M2_PI(1),th_PI(1),ph_PI(1)
     &                           ,acc_PI(1)
     &                           ,x_PI(1),y_PI(1),z_PI(1)
     &                           ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &                           ,PdHit_PI(1))
c+++ SC Trigger Threshold 
c           trig_SC(ntype) = max(thrSC_NU(1),thrSC_PI(1),thrSC_PI(2))
           id_SC_p(ntype) =  id_SC_NU(1)
           sec_SC_p(ntype) =  sec_SC_NU(1)
           thrSC_p(ntype) = thrSC_PI(1)
           id_SC_pip(ntype) = id_SC_PI(2)
           sec_SC_pip(ntype) = sec_SC_PI(2)
           thrSC_pip(ntype) = thrSC_PI(2)
           id_SC_pim(ntype) = id_SC_PI(1)
           sec_SC_pim(ntype) = sec_SC_PI(1)
           thrSC_pim(ntype) = thrSC_PI(1)
c+++ eNpp
           call  sas_fill_enpp(ntype)
c+++Rho
            M2_rho(ntype) = p_vec(4,2)
            call sas_fill_rho_lab(ntype,W_pp(0,1), W_pp(5,1), W_pp(6,1), W_pp(7,1))
            call sas_fill_rho_hcm(ntype,p_vec(0,2),p_vec(5,2),p_vec(6,2), p_vec(7,2)) 
c+++ Delta++
            M2_dpp(ntype) = p_vec(4,3)
            call sas_fill_dpp_lab(ntype, W_Np(0,2), W_Np(5,2), W_Np(6,2), W_Np(7,2))
            call sas_fill_dpp_hcm(ntype,p_vec(0,3),p_vec(5,3),p_vec(6,3),p_vec(7,3)) 
C+++ Delta0
            M2_d0(ntype) = p_vec(4,1)
            call sas_fill_d0_lab(ntype, W_Np(0,1), W_Np(5,1), W_Np(6,1), W_Np(7,1))
            call sas_fill_d0_hcm(ntype,p_vec(0,1),p_vec(5,1),p_vec(6,1),p_vec(7,1)) 
C+++ Decaying pion
            call sas_fill_pip_rho(ntype,p_mes(0,2),p_mes(5,2),p_mes(6,2),p_mes(7,2),psi(2))
            call sas_fill_pip_dpp(ntype,p_mes(0,3),p_mes(5,3),p_mes(6,3),p_mes(7,3),psi(3))
            call sas_fill_pim_d0 (ntype,p_mes(0,1),p_mes(5,1),p_mes(6,1),p_mes(7,1),psi(1))

             E_lab(1) = W_pp(0,1) ! Vector meson in lab
             P_lab(1) = W_pp(5,1)
            th_lab(1) = W_pp(6,1)
            ph_lab(1) = W_pp(7,1)

            E_cm(1) = p_vec(0,2)
            P_cm(1) = p_vec(5,2)
           th_cm(1) = p_vec(6,2)
           ph_cm(1) = p_vec(7,2)

            E_rest(1) = p_mes(0,2)
            P_rest(1) = p_mes(5,2)
           th_rest(1) = p_mes(6,2)
           ph_rest(1) = p_mes(7,2)

	endif
c--
c- End cmh/rest frame meson and baryons calc
c+ Calling Kinematic fit. IT's always performed so that
c+ MC in 1 DATA in 2place 2on MC var
        
        if(.not.beam_type) then  
         call fit_2pi
         indtype = ntype+1
c+++ proton
            det_p(ntype+1) = 0
            En(1) = sqrt(PM_KIN(1)**2+938**2)*1e-3
            Mo(1,1)=PM_KIN(1)*1e-3*sin(TH_KIN(1))*cos(FI_KIN(1))
            Mo(1,2)=PM_KIN(1)*1e-3*sin(TH_KIN(1))*sin(FI_KIN(1))
            Mo(1,3)=PM_KIN(1)*1e-3*cos(TH_KIN(1))
            Be    = PM_KIN(1)*1e-3/En(1)
            call sas_fill_p_lab(ntype+1, En(1) ,PM_KIN(1)*1e-3, Be
     &           ,0.938 ,TH_KIN(1)*raddeg ,FI_KIN(1)*raddeg
     &           ,-1000.
     &           ,-1000.,-1000.,-1000.
     &           ,-1000.,-1000.
     &           ,-1000)
c+++ pi+
            det_pip(ntype+1) = 0
            En(2) = sqrt(PM_KIN(2)**2+140**2)*1e-3
            Mo(2,1)=PM_KIN(2)*1e-3*sin(TH_KIN(2))*cos(FI_KIN(2))
            Mo(2,2)=PM_KIN(2)*1e-3*sin(TH_KIN(2))*sin(FI_KIN(2))
            Mo(2,3)=PM_KIN(2)*1e-3*cos(TH_KIN(2))
            Be    = PM_KIN(2)*1e-3/En(2)
            call sas_fill_pip_lab(ntype+1,En(2) ,PM_KIN(2)*1e-3
     &          ,Be, .140 ,TH_KIN(2)*raddeg ,FI_KIN(2)*raddeg
     &          ,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)
c+++ pi-
            det_pim(ntype+1) = 0
            En(3) = sqrt(PM_KIN(3)**2+140**2)*1e-3
            Mo(3,1)=PM_KIN(3)*1e-3*sin(TH_KIN(3))*cos(FI_KIN(3))
            Mo(3,2)=PM_KIN(3)*1e-3*sin(TH_KIN(3))*sin(FI_KIN(3))
            Mo(3,3)=PM_KIN(3)*1e-3*cos(TH_KIN(3))
            Be    = PM_KIN(3)*1e-3/En(3)
            call sas_fill_pim_lab(ntype+1,En(3)
     &          ,PM_KIN(3)*1e-3, Be
     &          ,0.140 ,TH_KIN(3)*raddeg ,FI_KIN(3)*raddeg
     &          ,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)
c+++ eNpp
           call  sas_fill_enpp(ntype+1)


            call lab2cms( 5
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,En(3),Mo(3,1),Mo(3,2),Mo(3,3)      ! pi-
     %                   ,En(2),Mo(2,1),Mo(2,2),Mo(2,3)      ! pi+
     %                   ,En(1),Mo(1,1),Mo(1,2),Mo(1,3)      ! p
     %                   ,p_vec ,p_mes
     %                   ,psi)


c+++Rho
            M2_rho(ntype+1) = p_vec(4,2)
            call sas_fill_rho_lab(ntype+1,-1000., -1000., -1000., -1000.)
            call sas_fill_rho_hcm(ntype+1,p_vec(0,2),p_vec(5,2),p_vec(6,2), p_vec(7,2)) 
c+++ Delta++
            M2_dpp(ntype+1) = p_vec(4,3)
            call sas_fill_dpp_lab(ntype+1,-1000., -1000., -1000., -1000.)
            call sas_fill_dpp_hcm(ntype+1,p_vec(0,3),p_vec(5,3),p_vec(6,3),p_vec(7,3)) 
C+++ Delta0
            M2_d0(ntype+1) = p_vec(4,1)
            call sas_fill_d0_lab(ntype+1,-1000., -1000., -1000., -1000.)
            call sas_fill_d0_hcm(ntype+1,p_vec(0,1),p_vec(5,1),p_vec(6,1),p_vec(7,1)) 
C+++ Decaying pion
            call sas_fill_pip_rho(ntype+1,p_mes(0,2),p_mes(5,2),p_mes(6,2),p_mes(7,2),psi(2))
            call sas_fill_pip_dpp(ntype+1,p_mes(0,3),p_mes(5,3),p_mes(6,3),p_mes(7,3),psi(3))
            call sas_fill_pim_d0 (ntype+1,p_mes(0,1),p_mes(5,1),p_mes(6,1),p_mes(7,1),psi(1))
C+++ Specif fit variables (in MeV and rad)
            call sas_fill_fit2pi (ntype+1,PM_GAM, EP_GAM, EP_KIN, ET_KIN, EF_KIN, CHI_KIN, IER_KIN)
           endif
c-

      return
      end 


