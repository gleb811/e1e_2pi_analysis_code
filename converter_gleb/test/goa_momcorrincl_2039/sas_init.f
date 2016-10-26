c------------------------------------------------------------------------- 
      SUBROUTINE sas_init(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE 
      include "ntpl_goa.inc"
      integer itype

      n_ev_run_a(itype) = -1000.
      trig_type_a(itype)= -1000.
      trig_clas_a(itype)= -1000.
      n_ev_nt_a(itype)  = -1000.
      n_part_a(itype)   = -1000.
      ntot_a(itype)     = -1000.
      bit_a(itype)      = -1000.
      bit1pi_a(itype)   = -1000.
      bit2pi_a(itype)   = -1000.
      bitomega_a(itype)   = -1000.

cc NOT      sec_coin_SC_a(itype)   = -1000.
cc NOT      sec_coin_ST_a(itype)   = -1000.
cc NOT      sec_coin_SC_thr_a(itype)   = -1000.
cc NOT      trig_SC(itype)             = -1000.

                 if(beam_type) then !skipping for photons
       Qgated_total   =  -1000.
       Qgated_plus    =  -1000.
       Omega_a(itype) =  -1000.
          Q2_a(itype) = -1000.
        W_el_a(itype) =  -1000.
         eps_a(itype)  =-1000.
        eps_l_a(itype) = -1000.
      Flux_vp_a(itype) = -1000.
         E_EL_a(itype) =  -1000.
         P_EL_a(itype) = -1000.
         B_EL_a(itype) = -1000.
        th_EL_a(itype) =  -1000.
        ph_EL_a(itype) = -1000. 
        acc_EL_a(itype)= -1000. 
        x_EL_a(itype) =  -1000.
        y_EL_a(itype) =  -1000. 
        z_EL_a(itype) =   -1000.
        ECtot_EL_a(itype) = -1000.
        ECin_EL_a(itype)  = -1000.
        ECout_EL_a(itype) = -1000.
        NpheCC_EL_a(itype) = -1000.
        EffCC_EL_a(itype) = -1000.
        PdHit_EL_a(itype) = -1000.
      else
          Egamma_a(itype)  = -1000.
            W_ph_a(itype)  =  -1000.   
        Tag_time_a(itype)  = -1000.
         ST_time_a(itype)  =-1000. 
        RF1_time_a(itype)  =-1000.
        RF2_time_a(itype)  =-1000. 
      hit_in_time_a(itype) =-1000.   
             chan_a(itype)=-1000. 
             true_a(itype)= -1000.
      endif
 
       det_P(itype)   =-1000.
       E_P_lab(itype) =-1000.
       P_P_lab(itype) =-1000.
       B_P_lab(itype) =-1000.
      M2_P_lab(itype) =-1000.
      th_P_lab(itype) =-1000.    
      ph_P_lab(itype) =-1000.
      acc_P_lab(itype)=-1000.
      x_P_lab(itype)  =-1000.
      y_P_lab(itype)  =-1000.
      z_P_lab(itype)  =-1000.
      th1_P_lab(itype)=-1000.
      th2_P_lab(itype)=-1000.
      id_SC_p(itype)      =-1000.
      sec_SC_p(itype)     =-1000.
      thrSC_p(itype)      =-1000.
      PdHit_P_lab(indtype)=-1000.

       det_n(itype)   =-1000.
       E_n_lab(itype) =-1000.
       P_n_lab(itype) =-1000.
       B_n_lab(itype) =-1000.
      M2_n_lab(itype) =-1000.
      th_n_lab(itype) =-1000.    
      ph_n_lab(itype) =-1000. 
      acc_n_lab(itype)=-1000. 
      x_n_lab(itype)  =-1000. 
      y_n_lab(itype)  =-1000. 
      z_n_lab(itype)  =-1000.
 
       det_Pip(itype)   =-1000. 
       E_PIp_lab(itype) =-1000.
       P_PIp_lab(itype) =-1000.
       B_PIp_lab(itype) =-1000.
      M2_PIp_lab(itype) =-1000.
      th_PIp_lab(itype) =-1000. 
      ph_PIp_lab(itype) =-1000.
      acc_PIp_lab(itype)=-1000.
      x_PIp_lab(itype)  =-1000.
      y_PIp_lab(itype)  =-1000.
      z_PIp_lab(itype)  =-1000.
      ECtot_PIp_lab(itype) =-1000.
      ECin_PIp_lab(itype)  =-1000.
      ECout_PIp_lab(itype) =-1000.
      PdHit_PIp_lab(indtype)=-1000.
      id_SC_pip(indtype)    =-1000.
      sec_SC_pip(indtype)   =-1000.
      thrSC_pip(indtype)    =-1000.
       E_PIp_hcm(itype) = -1000.
       P_PIp_hcm(itype) = -1000.
      th_PIp_hcm(itype) = -1000.    
      ph_PIp_hcm(itype) = -1000.

       det_PIm(itype)    =-1000.
       E_PIm_lab(itype) =-1000.
       P_PIm_lab(itype) =-1000.
       B_PIm_lab(itype) =-1000.
      M2_PIm_lab(itype) =-1000.
      th_PIm_lab(itype) =-1000.
      ph_PIm_lab(itype) =-1000.
      acc_PIm_lab(itype)=-1000.
      x_PIm_lab(itype)  =-1000.
      y_PIm_lab(itype)  =-1000.
      z_PIm_lab(itype)  =-1000.
      ECtot_PIm_lab(itype) =-1000.
      ECin_PIm_lab(itype)  =-1000.
      ECout_PIm_lab(itype) =-1000.
      id_SC_pim(indtype)   =-1000.
      sec_SC_pim(indtype)  =-1000.
      thrSC_pim(indtype)   =-1000.

      det_PI0(itype)    = -1000.
       E_PI0_lab(itype) = -1000.
       P_PI0_lab(itype) = -1000.
       B_PI0_lab(itype) = -1000.
      M2_PI0_lab(itype) = -1000.
      th_PI0_lab(itype) = -1000.      
      ph_PI0_lab(itype) = -1000.
       E_PI0_hcm(itype) = -1000.
       P_PI0_hcm(itype) = -1000.
      th_PI0_hcm(itype) = -1000.
      ph_PI0_hcm(itype) = -1000.


      det_G1(itype)    = -1000.
       E_G1_lab(itype) = -1000.
      th_G1_lab(itype) = -1000.    
      ph_G1_lab(itype) = -1000.
      acc_G1_lab(itype)= -1000.
      x_G1_lab(itype)  = -1000.
      y_G1_lab(itype)  = -1000.
      z_G1_lab(itype)  = -1000.

      det_G2(itype)    = -1000.
       E_G2_lab(itype) = -1000.
      th_G2_lab(itype) = -1000.   
      ph_G2_lab(itype) = -1000.
      acc_G2_lab(itype)= -1000.
      x_G2_lab(itype)  = -1000.
      y_G2_lab(itype)  = -1000.
      z_G2_lab(itype)  = -1000.

      M2_eN_a2(itype)   = -1000.
      M2_epip_a2(itype) = -1000.
      M2_epim_a2(itype) = -1000.
      M2_eNpip_a2(itype)= -1000.
      M2_eNpim_a2(itype)= -1000.

      E_eNpp_a(itype)  = -1000.
      M2_eNpp_a(itype) = -1000.
      W2_Npp_a(itype)  = -1000.
      M2_eNp_a(itype)  = -1000.
      W2_Np_a(itype)   = -1000.
       M2_eN_a(itype)  = -1000.
       E_eN_a(itype)   = -1000.
       P_eN_a(itype)   = -1000.
      th_eN_a(itype)   = -1000.    
      ph_eN_a(itype)   = -1000.
 
      M2_rho(itype)     = -1000.
       E_rho_lab(itype) = -1000.
       P_rho_lab(itype) = -1000.
      th_rho_lab(itype) = -1000.      
      ph_rho_lab(itype) = -1000.
       E_rho_hcm(itype) = -1000.
       P_rho_hcm(itype) = -1000. 
      th_rho_hcm(itype) = -1000.      
      ph_rho_hcm(itype) = -1000. 
 
      M2_Dpp(itype)     = -1000.
       E_dpp_lab(itype) = -1000.
       P_dpp_lab(itype) = -1000.
      th_dpp_lab(itype) = -1000.     
      ph_dpp_lab(itype) = -1000.
       E_dpp_hcm(itype) = -1000. 
       P_dpp_hcm(itype) = -1000.
      th_dpp_hcm(itype) = -1000.     
      ph_dpp_hcm(itype) = -1000.

      M2_dp(itype)     =  -1000.
       E_dp_lab(itype) =  -1000.
       P_dp_lab(itype) =  -1000.
      th_dp_lab(itype) =  -1000.     
      ph_dp_lab(itype) =  -1000.

       M2_d0(itype)    = -1000.
       E_d0_lab(itype) = -1000.
       P_d0_lab(itype) = -1000.
      th_d0_lab(itype) = -1000.     
      ph_d0_lab(itype) = -1000.
       E_d0_hcm(itype) = -1000.
       P_d0_hcm(itype) = -1000.
      th_d0_hcm(itype) = -1000.     
      ph_d0_hcm(itype) = -1000.

       E_PIp_rho(itype) = -1000.
       P_PIp_rho(itype) = -1000.
      th_PIp_rho(itype) = -1000.     
      ph_PIp_rho(itype) = -1000.
      psi_rho(itype)    = -1000.

       E_PIp_dpp(itype) = -1000.
       P_PIp_dpp(itype) = -1000.
      th_PIp_dpp(itype) = -1000.     
      ph_PIp_dpp(itype) = -1000.
      psi_dpp(itype)    = -1000.

      E_PIm_d0(itype)  = -1000.
       P_PIm_d0(itype) = -1000.
      th_PIm_d0(itype) = -1000.
      ph_PIm_d0(itype) = -1000.
      psi_d0(itype)    = -1000.


      if(.not.beam_type) then

        Eg_fit(itype) =  -1000.
       DEg_fit(itype) =  -1000.
         DP_P_fit(itype) =  -1000.
       DP_PIp_fit(itype) =  -1000.
       DP_PIm_fit(itype) =  -1000.
         DTH_P_fit(itype) =  -1000.
       DTH_PIp_fit(itype) =  -1000.
       DTH_PIm_fit(itype) =  -1000.
         DPH_P_fit(itype) =  -1000. 
       DPH_PIp_fit(itype) =  -1000.
       DPH_PIm_fit(itype) =  -1000.
        CHI_fit(itype) =  -1000.
        Err_fit(itype) =  -1000.

        M2_omega(itype) =  -1000.
        E_omega_lab(itype) =  -1000.
        P_omega_lab(itype) =  -1000.

        th_omega_lab(itype) =  -1000.
        ph_omega_lab(itype) =  -1000.
        E_omega_hcm(itype) =  -1000.
        P_omega_hcm(itype) =  -1000.
        th_omega_hcm(itype) =  -1000.
        ph_omega_hcm(itype) =  -1000.
        E_PIp_omega(itype) =  -1000.
        P_PIp_omega(itype) =  -1000.
        th_PIp_omega(itype) =  -1000.
        ph_PIp_omega(itype) =  -1000.
        psi_omega(itype) =  -1000.
        
        endif
      return
      end 
