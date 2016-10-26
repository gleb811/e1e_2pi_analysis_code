
      SUBROUTINE analysis_1pi_b(MC_ana)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      include "select.inc"
      integer MC_ana,ntype
        real raddeg, degrad
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

c+ Some general setting
c++ Choosing wich entry in the ntuple to fill
c++ Photon beam:
c++ MC only          : 1 <=> MC;  
c++ DATA only        : 1 <=> DA;  
c++ MC.and.DATA      : 1 <=> MC;  2 <=> DA 
c++ Electron beam:
c++ MC only          : 1 <=> MC;  
c++ DATA only        : 1 <=> DA;  
c++ MC.and.DATA      : 1 <=> MC; 2 <=> DA
            ntype = 1 ! By default it fills 1 variable (DATA.or.MC)
             if (bit_mc.eq.2) then ! If both are present: 1st = MC 2nd = data
               if(MC_ana.eq.0) ntype = 2 
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

c+++ Starting e p -> n pi+
	if((bit1pi.eq.211.or.bit1pi.eq.221).and.W.gt.0)  then	
            ncmrest=1
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eN(0),Pmiss_eN(1),Pmiss_eN(2),Pmiss_eN(3)
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)
     %                   ,p_vec,p_mes,psi)
c+++ neutron
            det_n(ntype) = 1
            call sas_fill_n_lab(ntype,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &                         ,th_NU(1),ph_NU(1)
     &                         ,acc_NU(1)
     &                         ,x_NU(1),y_NU(1),z_NU(1)) 
c+++ pip
            det_pip(ntype) = -1
            call sas_fill_pip_lab(ntype
     &          ,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(5)/Pmiss_eN(0)
     &          ,Pmiss_eN(4),Pmiss_eN(6),Pmiss_eN(7)
     &          ,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)

c+++ eNp
           call  sas_fill_enp_b(ntype)
C+++ Delta+
            M2_dp(ntype) = Pmiss_e(4)
            call sas_fill_dp_lab(ntype, Pmiss_e(0), Pmiss_e(5),Pmiss_e(6),Pmiss_e(7))
C+++ Decaying pion
            call sas_fill_pip_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

            E_lab(1)  = Pmiss_e(0)   ! N* in LAB frame
            P_lab(1)  = Pmiss_e(5) 
            th_lab(1) = Pmiss_e(6)
            ph_lab(1) = Pmiss_e(7) 

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) !  Pion in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif
	if((bit1pi.eq.212.or.bit1pi.eq.222).and.W.gt.0)  then	
            ncmrest=1
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_pi(1),Px_pi(1),Py_pi(1),Pz_pi(1)                       !pi+ det 
     %                   ,0.,0.,0.,0.
     %                   ,PMiss_ep(0,1),PMiss_ep(1,1),PMiss_ep(2,1),PMiss_ep(3,1) !n miss
     %                   ,p_vec ,p_mes, psi)
c+++ neutron
            det_n(ntype) = -1
            call sas_fill_n_lab(ntype
     &          ,Pmiss_ep(0,1),Pmiss_ep(5,1),Pmiss_ep(5,1)/Pmiss_ep(0,1)
     &          ,Pmiss_ep(4,1),Pmiss_ep(6,1),Pmiss_ep(7,1)
     &          ,-1000.
     &          ,-1000.,-1000.,-1000.)

c+++ pip
            det_pip(ntype) = 1
            call sas_fill_pip_lab (ntype
     &          ,E_PI(1),P_PI(1),B_PI(1)
     &          ,M2_PI(1),th_PI(1),ph_PI(1)
     &          ,acc_PI(1)
     &          ,x_PI(1),y_PI(1),z_PI(1)
     &          ,ECtot_PI(1), ECin_PI(1), ECout_PI(1) 
     &          ,PdHit_PI(1)) 

c+++ eNp
           call  sas_fill_enp_b(ntype)
C+++ Delta+
            M2_dp(ntype) = Pmiss_e(4)
            call sas_fill_dp_lab(ntype, Pmiss_e(0), Pmiss_e(5),Pmiss_e(6),Pmiss_e(7))
C+++ Decaying pion
            call sas_fill_pip_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

            E_lab(1)  = Pmiss_e(0)   ! N* in LAB frame
            P_lab(1)  = Pmiss_e(5) 
            th_lab(1) = Pmiss_e(6)
            ph_lab(1) = Pmiss_e(7) 

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) !  Pion in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif
	if((bit1pi.eq.213.or.bit1pi.eq.223).and.W.gt.0)  then	
            ncmrest=1
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_pi(1),Px_pi(1),Py_pi(1),Pz_pi(1)                       !pi+ det 
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)            !n det
     %                   ,p_vec,p_mes,psi)
c+++ neutron
            det_n(ntype) = 1
            call sas_fill_n_lab(ntype,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &                         ,th_NU(1),ph_NU(1)
     &                         ,acc_NU(1)
     &                         ,x_NU(1),y_NU(1),z_NU(1)) 
c+++ pip
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(1),P_PI(1),B_PI(1)
     &                           ,M2_PI(1),th_PI(1),ph_PI(1)
     &                           ,acc_PI(1)
     &                           ,x_PI(1),y_PI(1),z_PI(1)
     &                           ,ECtot_PI(1), ECin_PI(1), ECout_PI(1) 
     &                           ,PdHit_PI(1)) 

c+++ eNp
           call  sas_fill_enp_b(ntype)
C+++ Delta+
            M2_dp(ntype) = Pmiss_e(4)
            call sas_fill_dp_lab(ntype, Pmiss_e(0), Pmiss_e(5),Pmiss_e(6),Pmiss_e(7))
C+++ Decaying pion
            call sas_fill_pip_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))


            E_lab(1)  = Pmiss_e(0)   ! N* in LAB frame
            P_lab(1)  = Pmiss_e(5) 
            th_lab(1) = Pmiss_e(6)
            ph_lab(1) = Pmiss_e(7) 

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) !  Pion in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif
c+
c--- Ending e p -> p pi+


      return
      end 


