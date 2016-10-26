
      SUBROUTINE analysis_omega(MC_ana)
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
c+++ Starting e p -> p omega analysis
	if( ( (bitomega.eq.400.or.bitomega.eq.430).and.W.gt.0).or.
     %      ( (bitomega.eq.401.or.bitomega.eq.431).and.W.gt.0).or.
     %      ( (bitomega.eq.402.or.bitomega.eq.432).and.W.gt.0) )     then	
            ncmrest=1            
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eN(0),Pmiss_eN(1),Pmiss_eN(2),Pmiss_eN(3)  !miss omega
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)               ! p detected
     %                   ,p_vec ,p_mes, psi)

c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype
     &          ,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ gamma1&2
            if(bitomega.eq.400.or.bitomega.eq.430) then
             det_g1(ntype) = -1
             call sas_fill_g1_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.401.or.bitomega.eq.431) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.402.or.bitomega.eq.432) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = 1
             call sas_fill_g2_lab(ntype,E_g(2),th_g(2),ph_g(2)
     &                           ,acc_G(2)
     &                           ,x_G(2),y_G(2),z_G(2))

            endif

c+++ pi0
            if(det_g1(ntype).eq.1.and.det_g2(ntype).eq.1) then
             det_pi0(ntype) = 1
             call sas_fill_pi0_lab(ntype,W_GG(0,1),W_GG(5,1),W_GG(5,1)/W_GG(0,1)
     @                                 ,W_GG(4,1),W_GG(6,1),W_GG(7,1))
            else
             det_pi0(ntype) = -1
             call sas_fill_pi0_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.)
            endif
c+++ pi+
            det_pip(ntype) = -1
            call sas_fill_pip_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)
c+++ pi-
            det_pim(ntype) = -1
            call sas_fill_pim_lab(ntype
     &          ,-1000.0, -1000.0, -1000.0
     &          ,-1000.0, -1000.0, -1000.0
     &          ,-1000.0,
     &          ,-1000.0, -1000.0, -1000.0
     &          ,-1000.0, -1000.0, -1000.0
     &          ,-1000)
c+++ eNp
           call  sas_fill_enp_z(ntype)
c+++ eNpp
           call  sas_fill_enpp_z(ntype)

c+++Omega
            M2_omega(ntype) = p_mes(4,4)
            call sas_fill_omega_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_omega_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

C+++ Decaying pion
            call sas_fill_pip_omega(ntype,-1000.,-1000.,-1000.,-1000.,-1000.)

            E_lab(1)  = 0.  
            P_lab(1)  = 0.   
            th_lab(1) = 0.  
            ph_lab(1) = 0.  

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) ! Omega in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif

	if( ( (bitomega.eq.405.or.bitomega.eq.435).and.W.gt.0).or.
     %      ( (bitomega.eq.406.or.bitomega.eq.436).and.W.gt.0).or.
     %      ( (bitomega.eq.407.or.bitomega.eq.437).and.W.gt.0) )     then	
            ncmrest=1            
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eN(0),Pmiss_eN(1),Pmiss_eN(2),Pmiss_eN(3)  !miss omega
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)               ! p detected
     %                   ,p_vec ,p_mes, psi)

c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype
     &          ,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ gamma1&2
            if(bitomega.eq.405.or.bitomega.eq.435) then
             det_g1(ntype) = -1
             call sas_fill_g1_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.406.or.bitomega.eq.436) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.407.or.bitomega.eq.437) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = 1
             call sas_fill_g2_lab(ntype,E_g(2),th_g(2),ph_g(2)
     &                           ,acc_G(2)
     &                           ,x_G(2),y_G(2),z_G(2))

            endif

c+++ pi0
            if(det_g1(ntype).eq.1.and.det_g2(ntype).eq.1) then
             det_pi0(ntype) = 1
             call sas_fill_pi0_lab(ntype,W_GG(0,1),W_GG(5,1),W_GG(5,1)/W_GG(0,1)
     @                                 ,W_GG(4,1),W_GG(6,1),W_GG(7,1))
            else
             det_pi0(ntype) = -1
             call sas_fill_pi0_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.)
            endif
c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,
     &                         ,E_PI(1),P_PI(1),B_PI(1)
     &                         ,M2_PI(1),th_PI(1),ph_PI(1)
     &                         ,acc_PI(1)
     &                         ,x_PI(1),y_PI(1),z_PI(1)
     &                         ,ECtot_PI(1), ECin_PI(1), ECout_PI(1) 
     &                         ,PdHit_PI(1)) 
c+++ pi-
            det_pim(ntype) = -1
            call sas_fill_pim_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)
c+++ eNp
           call  sas_fill_enp_c(ntype)
c+++ eNpp
           call  sas_fill_enpp_z(ntype)


c+++Omega
            M2_omega(ntype) = p_mes(4,4)
            call sas_fill_omega_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_omega_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

C+++ Decaying pion
            call sas_fill_pip_omega(ntype,-1000.,-1000.,-1000.,-1000.,-1000.)

            E_lab(1)  = 0.  
            P_lab(1)  = 0.   
            th_lab(1) = 0.  
            ph_lab(1) = 0.  

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) ! Omega in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif

	if( ( (bitomega.eq.410.or.bitomega.eq.440).and.W.gt.0).or.
     %      ( (bitomega.eq.411.or.bitomega.eq.441).and.W.gt.0).or.
     %      ( (bitomega.eq.412.or.bitomega.eq.442).and.W.gt.0) )     then	
            ncmrest=1            
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eN(0),Pmiss_eN(1),Pmiss_eN(2),Pmiss_eN(3)  !miss omega
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)               ! p detected
     %                   ,p_vec ,p_mes, psi)

c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype
     &          ,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ gamma1&2
            if(bitomega.eq.410.or.bitomega.eq.440) then
             det_g1(ntype) = -1
             call sas_fill_g1_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.411.or.bitomega.eq.441) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.412.or.bitomega.eq.442) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = 1
             call sas_fill_g2_lab(ntype,E_g(2),th_g(2),ph_g(2)
     &                           ,acc_G(2)
     &                           ,x_G(2),y_G(2),z_G(2))

            endif

c+++ pi0
            if(det_g1(ntype).eq.1.and.det_g2(ntype).eq.1) then
             det_pi0(ntype) = 1
             call sas_fill_pi0_lab(ntype,W_GG(0,1),W_GG(5,1),W_GG(5,1)/W_GG(0,1)
     @                                 ,W_GG(4,1),W_GG(6,1),W_GG(7,1))
            else
             det_pi0(ntype) = -1
             call sas_fill_pi0_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.)
            endif
c+++ pi+
            det_pip(ntype) = -1
            call sas_fill_pip_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.
     &          ,-1000)
c+++ pi-
            det_pim(ntype) = +1
            call sas_fill_pim_lab(ntype
     &           ,E_PI(1), P_PI(1), B_PI(1)
     &           ,M2_PI(1),th_PI(1),ph_PI(1)
     &           ,acc_PI(1)
     &           ,x_PI(1),    y_PI(1),   z_PI(1)
     &           ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &           ,PdHit_PI(1))
c+++ eNp
           call  sas_fill_enp_d(ntype)
c+++ eNpp
           call  sas_fill_enpp_z(ntype)

c+++Omega
            M2_omega(ntype) = p_mes(4,4)
            call sas_fill_omega_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_omega_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

C+++ Decaying pion
            call sas_fill_pip_omega(ntype,-1000.,-1000.,-1000.,-1000.,-1000.)

            E_lab(1)  = 0.  
            P_lab(1)  = 0.   
            th_lab(1) = 0.  
            ph_lab(1) = 0.  

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) ! Omega in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif
	if( ( (bitomega.eq.415.or.bitomega.eq.445).and.W.gt.0).or.
     %      ( (bitomega.eq.416.or.bitomega.eq.446).and.W.gt.0).or.
     %      ( (bitomega.eq.417.or.bitomega.eq.447).and.W.gt.0) )     then	
            ncmrest=1            
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_eN(0),Pmiss_eN(1),Pmiss_eN(2),Pmiss_eN(3)  !miss omega
     %                   ,0.,0.,0.,0.
     %                   ,E_NU(1),PX_NU(1),PY_NU(1),PZ_NU(1)               ! p detected
     %                   ,p_vec ,p_mes, psi)

c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype
     &          ,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ gamma1&2
            if(bitomega.eq.415.or.bitomega.eq.445) then
             det_g1(ntype) = -1
             call sas_fill_g1_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.416.or.bitomega.eq.446) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = -1
             call  sas_fill_g2_lab(ntype,-1000.,-1000.,-1000.
     &                          ,-1000.
     &                          ,-1000.,-1000.,-1000.)
            else if (bitomega.eq.417.or.bitomega.eq.447) then
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = 1
             call sas_fill_g2_lab(ntype,E_g(2),th_g(2),ph_g(2)
     &                           ,acc_G(2)
     &                           ,x_G(2),y_G(2),z_G(2))

            endif

c+++ pi0
            if(det_g1(ntype).eq.1.and.det_g2(ntype).eq.1) then
             det_pi0(ntype) = 1
             call sas_fill_pi0_lab(ntype,W_GG(0,1),W_GG(5,1),W_GG(5,1)/W_GG(0,1)
     @                                 ,W_GG(4,1),W_GG(6,1),W_GG(7,1))
            else
             det_pi0(ntype) = -1
             call sas_fill_pi0_lab(ntype
     &          ,-1000.,-1000.,-1000.
     &          ,-1000.,-1000.,-1000.)
            endif
c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(1),P_PI(1),B_PI(1),M2_PI(1)
     &                           ,th_PI(1),ph_PI(1)
     &                           ,acc_PI(1)
     &                           ,x_PI(1),y_PI(1),z_PI(1)
     &                         ,ECtot_PI(1), ECin_PI(1), ECout_PI(1)
     &                         ,PdHit_PI(1) ) 
c+++ pi-
            det_pim(ntype) = +1
            call sas_fill_pim_lab(ntype
     &          ,E_PI(1),P_PI(1),B_PI(1)
     &          ,M2_PI(1),th_PI(1),ph_PI(1)
     &          ,acc_PI(1)
     &          ,x_PI(1),y_PI(1),z_PI(1)
     &          ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &          ,PdHit_PI(1))
c+++ eNp
           call  sas_fill_enp_z(ntype)
c+++ eNpp
           call  sas_fill_enpp(ntype)


c+++Omega
            M2_omega(ntype) = p_mes(4,4)
            call sas_fill_omega_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_omega_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

C+++ Decaying pion
            call sas_fill_pip_omega(ntype,-1000.,-1000.,-1000.,-1000.,-1000.)

            E_lab(1)  = 0.  
            P_lab(1)  = 0.   
            th_lab(1) = 0.  
            ph_lab(1) = 0.  

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) ! Omega in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif

	if((bitomega.eq.419.or.bitomega.eq.449).and.W.gt.0)     then	
            ncmrest=1            
            call lab2cms( 3
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,Pmiss_epp0(0,1),Pmiss_epp0(1,1),Pmiss_epp0(2,1),Pmiss_epp0(3,1)  !miss proton
     %                   ,0.,0.,0.,0.
     %                   ,W_pp0(0,1),W_pp0(1,1),W_pp0(2,1),W_pp0(3,1)       ! omega detected
     %                   ,p_vec ,p_mes, psi)

c+++ proton
            det_p(ntype) = 1
            call sas_fill_p_lab(ntype
     &          ,E_NU(1),P_NU(1),B_NU(1),M2_NU(1)
     &          ,th_NU(1),ph_NU(1)
     &          ,acc_NU(1)
     &          ,x_NU(1),y_NU(1),z_NU(1)
     &          ,th1_NU(1),th2_NU(1)
     &          ,PdHit_NU(1))
c+++ gamma1&2
             det_g1(ntype) = 1
             call sas_fill_g1_lab(ntype,E_g(1),th_g(1),ph_g(1)
     &                           ,acc_G(1)
     &                           ,x_G(1),y_G(1),z_G(1))
             det_g2(ntype) = 1
             call sas_fill_g2_lab(ntype,E_g(2),th_g(2),ph_g(2)
     &                           ,acc_G(2)
     &                           ,x_G(2),y_G(2),z_G(2))


c+++ pi0
             det_pi0(ntype) = 1
             call sas_fill_pi0_lab(ntype,W_GG(0,1),W_GG(5,1),W_GG(5,1)/W_GG(0,1)
     @                                 ,W_GG(4,1),W_GG(6,1),W_GG(7,1))

c+++ pi+
            det_pip(ntype) = 1
            call sas_fill_pip_lab(ntype,E_PI(1),P_PI(1),B_PI(1),M2_PI(1)
     &                           ,th_PI(1),ph_PI(1)
     &                           ,acc_PI(1)
     &                           ,x_PI(1),y_PI(1),z_PI(1)
     &                         ,ECtot_PI(1), ECin_PI(1), ECout_PI(1) 
     &                         ,PdHit_PI(1)) 
c+++ pi-
            det_pim(ntype) = +1
            call sas_fill_pim_lab(ntype
     &          ,E_PI(1),P_PI(1),B_PI(1)
     &          ,M2_PI(1),th_PI(1),ph_PI(1)
     &          ,acc_PI(1)
     &          ,x_PI(1),y_PI(1),z_PI(1)
     &          ,ECtot_PI(1),ECin_PI(1),ECout_PI(1)
     &          ,PdHit_PI(1))
c+++ eNp
           call  sas_fill_enp_z(ntype)
c+++ eNpp
           call  sas_fill_enpp(ntype)


c+++Omega
            M2_omega(ntype) = p_mes(4,4)
            call sas_fill_omega_lab(ntype,Pmiss_eN(0),Pmiss_eN(5),Pmiss_eN(6),Pmiss_eN(7))
            call sas_fill_omega_hcm(ntype,p_mes(0,4),p_mes(5,4),p_mes(6,4),p_mes(7,4))

C+++ Decaying pion
            call sas_fill_pip_omega(ntype,-1000.,-1000.,-1000.,-1000.,-1000.)

            E_lab(1)  = 0.  
            P_lab(1)  = 0.   
            th_lab(1) = 0.  
            ph_lab(1) = 0.  

            E_cm(1) = 0. 
            P_cm(1) = 0.
           th_cm(1) = 0.
           ph_cm(1) = 0.

            E_rest(1) = p_mes(0,4) ! Omega in hadronic cm frame
            P_rest(1) = p_mes(5,4)
           th_rest(1) = p_mes(6,4)
           ph_rest(1) = p_mes(7,4)

	endif

      return
      end 


