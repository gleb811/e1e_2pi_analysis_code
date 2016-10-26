c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_genrun(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE 
      include "ntpl_goa.inc"
      integer itype
              
      n_ev_run_a(itype) = n_ev_run
      trig_type_a(itype)= tg_hevt !trig_type
      trig_clas_a(itype)= trig_clas
      n_ev_nt_a(itype)  = n_ev_nt
      n_part_a(itype)   = n_part
      ntot_a(itype)     = ntot
      bit_a(itype)      = bit
      bit1pi_a(itype)   = bit1pi
      bit2pi_a(itype)   = bit2pi
      bitomega_a(itype)   = bitomega
c      sec_coin_SC_a(itype)   = sec_coin_SC
      sec_coin_ST_a(itype)   = sec_coin_ST

      return
      end 


c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_beamrel(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
                 if(beam_type) then !skipping for photons
       Omega_a(itype) =  Omega
          Q2_a(itype) =   Q2
        W_el_a(itype) =   W
         eps_a(itype)  =  eps
        eps_l_a(itype) =  eps_l
      Flux_vp_a(itype) =  Flux_vp
         E_EL_a(itype) =  E_EL
         P_EL_a(itype) =  P_EL
         B_EL_a(itype) =  B_EL
        th_EL_a(itype) =  th_EL
        ph_EL_a(itype) =  ph_EL
        acc_EL_a(itype)=  acc_EL
        x_EL_a(itype) =   x_EL
        y_EL_a(itype) =   y_EL
        z_EL_a(itype) =   z_EL
        ECtot_EL_a(itype) = ECtot_EL
        ECin_EL_a(itype)  = ECin_EL
        ECout_EL_a(itype) = ECout_EL
        NpheCC_EL_a(itype)= NpheCC_EL 
        EffCC_EL_a(itype) = EffCC_EL
        PdHit_EL_a(itype) = PdHit_EL
      else
          Egamma_a(itype)  = Egamma
            W_ph_a(itype)  =  W      
        Tag_time_a(itype)  = Tag_time
         ST_time_a(itype)  = ST_time
        RF1_time_a(itype)  = RF1_time
        RF2_time_a(itype)  = RF2_time
      hit_in_time_a(itype) = hit_in_time  
             chan_a(itype)= chan
             true_a(itype)= true
             T_id_a(itype)= T_id 
             E_id_a(itype)= E_id
      endif

      M2_eN_a2(itype)      =-1000
      M2_epim_a2(itype)    =-1000
      M2_epip_a2(itype)    =-1000
      M2_eNpim_a2(itype)   =-1000
      M2_eNpip_a2(itype)   =-1000
      if(nnu.eq.1.and.neN.eq.1) then
        M2_eN_a2(itype) = M2_eN(1)
      endif
      if(npi.ge.1.and.nep.ge.1)then
        if(q_pi(1).lt.0) M2_epim_a2(itype)=M2_ep(1)
        if(q_pi(1).gt.0) M2_epip_a2(itype)=M2_ep(1)
      endif
      if(npi.ge.2.and.nep.ge.2)then
        if(q_pi(2).lt.0) M2_epim_a2(itype)=M2_ep(2)
        if(q_pi(2).gt.0) M2_epip_a2(itype)=M2_ep(2)
      endif
      if(npi.ge.1.and.neNp.ge.1)then
        if(q_pi(1).lt.0) M2_eNpim_a2(itype)=M2_eNp(1)
        if(q_pi(1).gt.0) M2_eNpip_a2(itype)=M2_eNp(1)
      endif
      if(npi.ge.2.and.neNp.ge.2)then
        if(q_pi(2).lt.0) M2_eNpim_a2(itype)=M2_eNp(2)
        if(q_pi(2).gt.0) M2_eNpip_a2(itype)=M2_eNp(2)
      endif

      return
      end

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_p_lab(itype,
     &           E,P,B,M2,th,ph,acc,x,y,z,th1,th2,PdH)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype,PdH
      real E,P,B,M2,th,ph,acc,x,y,z,th1,th2
       E_P_lab(itype) = E
       P_P_lab(itype) = P
       B_P_lab(itype) = B
      M2_P_lab(itype) = M2
      th_P_lab(itype) = th      
      ph_P_lab(itype) = ph
      acc_P_lab(itype)= acc
      x_P_lab(itype)  = x
      y_P_lab(itype)  = y
      z_P_lab(itype)  = z
      th1_P_lab(itype)= th1
      th2_P_lab(itype)= th2
      PdHit_P_lab(itype)=PdH
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_n_lab(itype,E,P,B,M2,th,ph,acc,x,y,z)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,B,M2,th,ph,acc,x,y,z
       E_n_lab(itype) = E
       P_n_lab(itype) = P
       B_n_lab(itype) = B
      M2_n_lab(itype) = M2
      th_n_lab(itype) = th      
      ph_n_lab(itype) = ph
      acc_n_lab(itype)= acc
      x_n_lab(itype)  = x
      y_n_lab(itype)  = y
      z_n_lab(itype)  = z
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pip_lab(itype,E,P,B,M2,th,ph,
     &                            acc,x,y,z,ECtot,ECin,ECout,PdH)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype,Pdh
      real E,P,B,M2,th,ph,acc,x,y,z,ECtot,ECin,ECout
       E_PIp_lab(itype) = E
       P_PIp_lab(itype) = P
       B_PIp_lab(itype) = B
      M2_PIp_lab(itype) = M2
      th_PIp_lab(itype) = th
      ph_PIp_lab(itype) = ph
      acc_PIp_lab(itype) = acc
      x_PIp_lab(itype) = x
      y_PIp_lab(itype) = y
      z_PIp_lab(itype) = z
      ECtot_PIp_lab(itype) = ECtot
      ECin_PIp_lab(itype)  = ECin
      ECout_PIp_lab(itype) = ECout
      PdHit_PIp_lab(itype) = PdH
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pim_lab(itype,E,P,B,M2,th,ph,
     &                            acc,x,y,z,ECtot,ECin,ECout,PdH)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype,PdH
      real E,P,B,M2,th,ph,acc,x,y,z,ECtot,ECin,ECout
       E_PIm_lab(itype) = E
       P_PIm_lab(itype) = P
       B_PIm_lab(itype) = B
      M2_PIm_lab(itype) = M2
      th_PIm_lab(itype) = th      
      ph_PIm_lab(itype) = ph
      acc_PIm_lab(itype)= acc
      x_PIm_lab(itype) = x
      y_PIm_lab(itype) = y
      z_PIm_lab(itype) = z
      ECtot_PIm_lab(itype) = ECtot
      ECin_PIm_lab(itype)  = ECin
      ECout_PIm_lab(itype) = ECout
      PdHit_PIm_lab(itype) = PdH
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pi0_lab(itype,E,P,B,M2,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,B,M2,th,ph
       E_PI0_lab(itype) = E
       P_PI0_lab(itype) = P
       B_PI0_lab(itype) = B
      M2_PI0_lab(itype) = M2
      th_PI0_lab(itype) = th      
      ph_PI0_lab(itype) = ph
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_g1_lab(itype,E,th,ph,acc,x,y,z)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,th,ph,x,y,z,acc
       E_G1_lab(itype) = E
      th_G1_lab(itype) = th      
      ph_G1_lab(itype) = ph
      acc_G1_lab(itype)  = x
      x_G1_lab(itype)  = x
      y_G1_lab(itype)  = y
      z_G1_lab(itype)  = z
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_g2_lab(itype,E,th,ph,acc,x,y,z)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,th,ph,x,y,z,acc
       E_G2_lab(itype) = E
      th_G2_lab(itype) = th      
      ph_G2_lab(itype) = ph
      acc_G2_lab(itype)  = x
      x_G2_lab(itype)  = x
      y_G2_lab(itype)  = y
      z_G2_lab(itype)  = z
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enpp(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real P_a(3,3),c,dummy(7)
      data c/0.01745/
      
c+ Calculating some variables using precedents variables
            P_a(1,1)=P_P_lab(itype)*sin(c*Th_P_lab(itype))*cos(c*Ph_P_lab(itype)) ! Proton in lab
            P_a(1,2)=P_P_lab(itype)*sin(c*Th_P_lab(itype))*sin(c*Ph_P_lab(itype))
            P_a(1,3)=P_P_lab(itype)*cos(c*Th_P_lab(itype))
            P_a(2,1)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*cos(c*Ph_Pip_lab(itype)) ! Pip in lab
            P_a(2,2)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*sin(c*Ph_Pip_lab(itype))
            P_a(2,3)=P_Pip_lab(itype)*cos(c*Th_Pip_lab(itype))
            P_a(3,1)=P_Pim_lab(itype)*sin(c*Th_Pim_lab(itype))*cos(c*Ph_Pim_lab(itype)) ! Pim in lab
            P_a(3,2)=P_Pim_lab(itype)*sin(c*Th_pim_lab(itype))*sin(c*Ph_Pim_lab(itype))
            P_a(3,3)=P_Pim_lab(itype)*cos(c*Th_Pim_lab(itype))
          call miss_mass_4part(0.9382,E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                        ,E_pim_lab(itype),P_a(3,1),P_a(3,2),P_a(3,3)
     %                        ,E_P_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                        ,E_eNpp_a(itype),dummy(2),dummy(3),dummy(4)
     %                        ,M2_eNpp_a(itype),dummy(5),dummy(6),dummy(7))
           call invariant_mass_3part
     %                     (E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                     ,E_pim_lab(itype),P_a(3,1),P_a(3,2),P_a(3,3)
     %                     ,E_P_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                     ,dummy(1)          ,dummy(2),dummy(3),dummy(4)
     %                     ,W2_Npp_a(itype),dummy(5),dummy(6),dummy(7))

c-
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enpp_z(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      
       E_eNpp_a(itype) = -1000.
       M2_eNpp_a(itype)= -1000.
       W2_Npp_a(itype) = -1000.
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enp_a(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real P_a(2,3),c,dummy(7)
      data c/0.01745/
      
c+ Calculating some variables using precedents variables
            P_a(1,1)=P_P_lab(itype)*sin(c*Th_P_lab(itype))*cos(c*Ph_P_lab(itype)) ! Proton in lab
            P_a(1,2)=P_P_lab(itype)*sin(c*Th_P_lab(itype))*sin(c*Ph_P_lab(itype))
            P_a(1,3)=P_P_lab(itype)*cos(c*Th_P_lab(itype))
            P_a(2,1)=P_Pi0_lab(itype)*sin(c*Th_Pi0_lab(itype))*cos(c*Ph_Pi0_lab(itype)) ! Pi0 in lab
            P_a(2,2)=P_Pi0_lab(itype)*sin(c*Th_Pi0_lab(itype))*sin(c*Ph_Pi0_lab(itype))
            P_a(2,3)=P_Pi0_lab(itype)*cos(c*Th_Pi0_lab(itype))
          call miss_mass_3part(0.9382,E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_pi0_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                        ,E_P_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                        ,dummy(1),dummy(2),dummy(3),dummy(4)
     %                        ,M2_eNp_a(itype),dummy(5),dummy(6),dummy(7))
           call invariant_mass_2part
     %                     (E_pi0_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                     ,E_P_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                     ,dummy(1)          ,dummy(2),dummy(3),dummy(4)
     %                     ,W2_Np_a(itype),dummy(5),dummy(6),dummy(7))


c-
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enp_b(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real P_a(2,3),c,dummy(7)
      data c/0.01745/
      
c+ Calculating some variables using precedents variables
            P_a(1,1)=P_n_lab(itype)*sin(c*Th_n_lab(itype))*cos(c*Ph_n_lab(itype)) ! neutron in lab
            P_a(1,2)=P_n_lab(itype)*sin(c*Th_n_lab(itype))*sin(c*Ph_n_lab(itype))
            P_a(1,3)=P_n_lab(itype)*cos(c*Th_n_lab(itype))
            P_a(2,1)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*cos(c*Ph_Pip_lab(itype)) ! Pip in lab
            P_a(2,2)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*sin(c*Ph_Pip_lab(itype))
            P_a(2,3)=P_Pip_lab(itype)*cos(c*Th_Pip_lab(itype))
          call miss_mass_3part(0.9382,E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                        ,E_n_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                        ,dummy(1),dummy(2),dummy(3),dummy(4)
     %                        ,M2_eNp_a(itype),dummy(5),dummy(6),dummy(7))
           call invariant_mass_2part
     %                     (E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                     ,E_n_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                     ,dummy(1)          ,dummy(2),dummy(3),dummy(4)
     %                     ,W2_Np_a(itype),dummy(5),dummy(6),dummy(7))

c-
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enp_c(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real P_a(2,3),c,dummy(7)
      data c/0.01745/
      
c+ Calculating some variables using precedents variables
            P_a(1,1)=P_p_lab(itype)*sin(c*Th_p_lab(itype))*cos(c*Ph_p_lab(itype)) ! proton in lab
            P_a(1,2)=P_p_lab(itype)*sin(c*Th_p_lab(itype))*sin(c*Ph_p_lab(itype))
            P_a(1,3)=P_p_lab(itype)*cos(c*Th_p_lab(itype))
            P_a(2,1)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*cos(c*Ph_Pip_lab(itype)) ! Pip in lab
            P_a(2,2)=P_Pip_lab(itype)*sin(c*Th_Pip_lab(itype))*sin(c*Ph_Pip_lab(itype))
            P_a(2,3)=P_Pip_lab(itype)*cos(c*Th_Pip_lab(itype))
          call miss_mass_3part(0.9382,E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                        ,E_p_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                        ,dummy(1),dummy(2),dummy(3),dummy(4)
     %                        ,M2_eNp_a(itype),dummy(5),dummy(6),dummy(7))
           call invariant_mass_2part
     %                     (E_pip_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                     ,E_p_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                     ,dummy(1)          ,dummy(2),dummy(3),dummy(4)
     %                     ,W2_Np_a(itype),dummy(5),dummy(6),dummy(7))

c-
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enp_d(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real P_a(2,3),c,dummy(7)
      data c/0.01745/
      
c+ Calculating some variables using precedents variables
            P_a(1,1)=P_p_lab(itype)*sin(c*Th_p_lab(itype))*cos(c*Ph_p_lab(itype)) ! proton in lab
            P_a(1,2)=P_p_lab(itype)*sin(c*Th_p_lab(itype))*sin(c*Ph_p_lab(itype))
            P_a(1,3)=P_p_lab(itype)*cos(c*Th_p_lab(itype))
            P_a(2,1)=P_Pim_lab(itype)*sin(c*Th_Pim_lab(itype))*cos(c*Ph_Pim_lab(itype)) ! Pim in lab
            P_a(2,2)=P_Pim_lab(itype)*sin(c*Th_Pim_lab(itype))*sin(c*Ph_Pim_lab(itype))
            P_a(2,3)=P_Pim_lab(itype)*cos(c*Th_Pim_lab(itype))
          call miss_mass_3part(0.9382,E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_pim_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                        ,E_p_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                        ,dummy(1),dummy(2),dummy(3),dummy(4)
     %                        ,M2_eNp_a(itype),dummy(5),dummy(6),dummy(7))
           call invariant_mass_2part
     %                     (E_pim_lab(itype),P_a(2,1),P_a(2,2),P_a(2,3)
     %                     ,E_p_lab(itype),  P_a(1,1),P_a(1,2),P_a(1,3)
     %                     ,dummy(1)          ,dummy(2),dummy(3),dummy(4)
     %                     ,W2_Np_a(itype),dummy(5),dummy(6),dummy(7))

c-
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_enp_z(itype)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
        M2_eNp_a(itype)=-1000
        W2_Np_a(itype)=-1000
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_eN_a(itype,M2,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real M2,E,P,th,ph
       M2_eN_a(itype) = M2
       E_eN_a(itype)  = E
       P_eN_a(itype)  = P
      th_eN_a(itype)  = th      
      ph_eN_a(itype)  = ph
      return
      end 



c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_rho_lab(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_rho_lab(itype) = E
       P_rho_lab(itype) = P
      th_rho_lab(itype) = th      
      ph_rho_lab(itype) = ph
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_omega_lab(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_OMEGA_lab(itype) = E
       P_OMEGA_lab(itype) = P
      th_OMEGA_lab(itype) = th      
      ph_OMEGA_lab(itype) = ph
      return
      end 


c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_dpp_lab(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_dpp_lab(itype) = E
       P_dpp_lab(itype) = P
      th_dpp_lab(itype) = th      
      ph_dpp_lab(itype) = ph
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_dp_lab(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_dp_lab(itype) = E
       P_dp_lab(itype) = P
      th_dp_lab(itype) = th      
      ph_dp_lab(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_d0_lab(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_d0_lab(itype) = E
       P_d0_lab(itype) = P
      th_d0_lab(itype) = th      
      ph_d0_lab(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_rho_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_rho_hcm(itype) = E
       P_rho_hcm(itype) = P
      th_rho_hcm(itype) = th      
      ph_rho_hcm(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_omega_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_omega_hcm(itype) = E
       P_omega_hcm(itype) = P
      th_omega_hcm(itype) = th      
      ph_omega_hcm(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_dpp_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_dpp_hcm(itype) = E
       P_dpp_hcm(itype) = P
      th_dpp_hcm(itype) = th      
      ph_dpp_hcm(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_d0_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_d0_hcm(itype) = E
       P_d0_hcm(itype) = P
      th_d0_hcm(itype) = th      
      ph_d0_hcm(itype) = ph
      return
      end 

c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pip_rho(itype,E,P,th,ph,ps)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph,ps
       E_PIp_rho(itype) = E
       P_PIp_rho(itype) = P
      th_PIp_rho(itype) = th      
      ph_PIp_rho(itype) = ph
      psi_rho(itype) = ps
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pip_omega(itype,E,P,th,ph,ps)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph,ps
       E_PIp_omega(itype) = E
       P_PIp_omega(itype) = P
      th_PIp_omega(itype) = th      
      ph_PIp_omega(itype) = ph
      psi_omega(itype) = ps
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pip_dpp(itype,E,P,th,ph,ps)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph,ps
       E_PIp_dpp(itype) = E
       P_PIp_dpp(itype) = P
      th_PIp_dpp(itype) = th      
      ph_PIp_dpp(itype) = ph
      psi_dpp(itype)    = ps
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pim_d0(itype,E,P,th,ph,ps)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph,ps
       E_PIm_d0(itype) = E
       P_PIm_d0(itype) = P
      th_PIm_d0(itype) = th      
      ph_PIm_d0(itype) = ph
      psi_d0(itype)    = ps
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pi0_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_PI0_hcm(itype) = E
       P_PI0_hcm(itype) = P
      th_PI0_hcm(itype) = th      
      ph_PI0_hcm(itype) = ph
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_pip_hcm(itype,E,P,th,ph)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real E,P,th,ph
       E_PIp_hcm(itype) = E
       P_PIp_hcm(itype) = P
      th_PIp_hcm(itype) = th      
      ph_PIp_hcm(itype) = ph
      return
      end 
c------------------------------------------------------------------------- 
      SUBROUTINE sas_fill_fit2pi(itype,Eg,de,dp,dt,df,chi,ier)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      integer itype
      real eg,de,dp(3),dt(3),df(3),chi
      integer  ier   
        Eg_fit(itype) = Eg*1e-3
       DEg_fit(itype) = de*1e-3
         DP_P_fit(itype) = dp(1)*1e-3 
       DP_PIp_fit(itype) = dp(2)*1e-3
       DP_PIm_fit(itype) = dp(3)*1e-3
         DTH_P_fit(itype) = dt(1)*57.2958 
       DTH_PIp_fit(itype) = dt(2)*57.2958 
       DTH_PIm_fit(itype) = dt(3)*57.2958 
         DPH_P_fit(itype) = df(1)*57.2958  
       DPH_PIp_fit(itype) = df(2)*57.2958 
       DPH_PIm_fit(itype) = df(3)*57.2958 
        CHI_fit(itype) = chi
        Err_fit(itype) = ier
      return
      end 
