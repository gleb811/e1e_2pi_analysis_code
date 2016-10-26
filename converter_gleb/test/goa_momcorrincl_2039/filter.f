c*************************************************
	SUBROUTINE filter(OK)
c*************************************************
	include "ntpl_goa.inc"
	logical OK, flag_g6a, flag_g6b		      	
	real C_Ebeam,C_phimin,C_Wmin ,C_Wmax ,C_Wgamma,
     %     C_M2_eN_min,C_M2_eN_1,C_M2_eN_2 ,C_Elast,C_Q2min 
	integer ind_MC_data
c	real cut1,cut2,cut3
	OK = .FALSE.

	ind_MC_data=1 ! Only real data (no filter on MC data)
	

c+ Filtering events
c+ Initialisation
	      
	     C_Ebeam     = 4.045
	     C_phimin    =   30
	     C_Wmin      = 1.04
	     C_Wmax      = 1.4
	     C_Wgamma    = 0.15
	     C_M2_eN_min = -.2
	     C_M2_eN_1   = 0.005
	     C_M2_eN_2   = 0.06
	     C_Elast     = -1.
	     C_Q2min     = 0.5     
c-	     
c+ Cuts according g6 analysis 2/11/99 
	flag_g6a=  
     % abs(tag_time_a(ind_MC_data)-st_time_a(ind_MC_data)).lt.4.and.
     % hit_in_time_a(ind_MC_data).eq.1.and.
     % Egamma_a(ind_MC_data).gt.3.197.and.Egamma_a(ind_MC_data).lt.3.917.and.w_ph_a(ind_MC_data).gt.2.6.and.
     % sqrt(abs((det_pip(ind_MC_data)-1.)/2.*M2_PIp_lab(ind_MC_data))).lt.0.4.and.
     % sqrt(abs((det_pim(ind_MC_data)-1.)/2.*M2_PIm_lab(ind_MC_data))).lt.0.4.and.
     % sqrt(abs((det_P(ind_MC_data)-1.)/2.*M2_P_lab(ind_MC_data))).lt.1.15.and.
     %  sqrt(abs(M2_eNpp_a(ind_MC_data))).lt.0.15
c     % .and.-2*0.938*(e_p_lab(ind_MC_data)-0.938).lt.-0.5  

	flag_g6b=  
     % abs(tag_time_a(ind_MC_data)-st_time_a(ind_MC_data)).lt.3.and.
     % hit_in_time_a(ind_MC_data).eq.1.and.
c     % Egamma_a(ind_MC_data).gt.3.197.and.Egamma_a(ind_MC_data).lt.3.917.and.w_ph_a(ind_MC_data).gt.2.6.and.
     % sqrt(abs((det_pip(ind_MC_data)-1.)/2.*M2_PIp_lab(ind_MC_data))).lt.0.4.and.
     % sqrt(abs((det_pim(ind_MC_data)-1.)/2.*M2_PIm_lab(ind_MC_data))).lt.0.4.and.
     % sqrt(abs((det_P(ind_MC_data)-1.)/2.*M2_P_lab(ind_MC_data))).lt.1.15.and.
     %  sqrt(abs(M2_eNpp_a(ind_MC_data))).lt.0.15
c     % .and.-2*0.938*(e_p_lab(ind_MC_data)-0.938).lt.-0.5  
      
          if(flag_g6b)    OK = .true.


	return
	end      	
	      	
	      	
	      	
	      	
