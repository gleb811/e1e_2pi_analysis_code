c+ Misssing mass matrix structure Mij
c+ i index 
c+         0 = Energy of mmising particle
c+         1 = px of mmising particle
c+         2 = py of mmising particle
c+         3 = pz of mmising particle
c+         4 = sqared mass of missing particle
c+         5 = momentum of mmising particle
c+         6 = theta of mmising particle
c+         7 = phi of mmising particle
c+ j index runs over possibles combination (e.g. j=1 (e pi1) j=2 (e pi2) j=3 (e pi3)) 
c+         
	REAL*4         PMiss_e, PMiss_eN, PMiss_ep, PMiss_ek
     &               ,PMiss_eNp, PMiss_epp, PMiss_eNk, PMiss_epk, PMiss_ekk
     &               ,PMiss_eNpp, PMiss_eppp, PMiss_eNkk, PMiss_eNpk
     &               ,PMiss_eNppp
       common/missing/   PMiss_e(0:7), PMiss_eN(0:7), PMiss_ep(0:7,3), PMiss_ek(0:7,3)
     &               ,PMiss_eNp(0:7,3), PMiss_epp(0:7,3), PMiss_eNk(0:7,3), PMiss_epk(0:7,3), PMiss_ekk(0:7,3)
     &               ,PMiss_eNpp(0:7,3), PMiss_eppp(0:7,3), PMiss_eNkk(0:7,3), PMiss_eNpk(0:7,3)
     &               ,PMiss_eNppp(0:7,3)
     
c+ 0 = (GG) combination
 	REAL*4        
     &               PMiss_eNG
     &               ,PMiss_e0, PMiss_eN0, PMiss_ep0, PMiss_ek0
     &               ,PMiss_eNp0, PMiss_epp0, PMiss_eNk0, PMiss_epk0, PMiss_ekk0
     &               ,PMiss_eNpp0, PMiss_eppp0, PMiss_eNkk0, PMiss_eNpk0
     &               ,PMiss_eNppp0       
       common/miss_neutral/   
     &               PMiss_eNG(0:7,4)
     &               ,PMiss_e0(0:7), PMiss_eN0(0:7), PMiss_ep0(0:7,3), PMiss_ek0(0:7,3)
     &               ,PMiss_eNp0(0:7,3), PMiss_epp0(0:7,3), PMiss_eNk0(0:7,3), PMiss_epk0(0:7,3), PMiss_ekk0(0:7,3)
     &               ,PMiss_eNpp0(0:7,3), PMiss_eppp0(0:7,3), PMiss_eNkk0(0:7,3), PMiss_eNpk0(0:7,3)
     &               ,PMiss_eNppp0(0:7,3)       
       
c-
c+ Invariant mass matrix structure Mij:
c+ i index 
c+         0 = Energy of invariant particle
c+         1 = px of invariant particle
c+         2 = py of invariant particle
c+         3 = pz of invariant particle
c+         4 = sqared mass of invariant particle
c+         5 = momentum of invariant particle
c+         6 = theta of invariant particle
c+         7 = phi of invariant particle
c+ j index runs over possibles combination (e.g. j=1 (N pi1) j=2 (N pi2) j=3 (N pi3)) 
c+         
	REAL*4      
     &               W_Np, W_pp, W_Nk, W_pk, W_kk
     &               ,W_Npp, W_ppp, W_Nkk, W_Npk
     &               ,W_Nppp
       common/invariant/   
     &               W_Np(0:7,3), W_pp(0:7,3), W_Nk(0:7,3), W_pk(0:7,3), W_kk(0:7,3)
     &               ,W_Npp(0:7,3), W_ppp(0:7,3), W_Nkk(0:7,3), W_Npk(0:7,3)
     &               ,W_Nppp(0:7,3)
     

        REAL*4 W_GG, W_GGG, W_GGGG
       common/inv_gamma/     W_GG(0:7,6), W_GGG(0:7,4), W_GGGG(0:7,1)   
c+ 0 = (GG) combination
	REAL*4  
     &               W_NG
     &               ,W_N0, W_p0, W_k0
     &               ,W_Np0, W_pp0, W_Nk0, W_pk0, W_kk0
     &               ,W_Npp0, W_ppp0, W_Nkk0, W_Npk0
     &               ,W_Nppp0
       common/inv_neutral/   
     &               W_NG(0:7,4)
     &               ,W_N0(0:7), W_p0(0:7,3), W_k0(0:7,3) 
     &               ,W_Np0(0:7,3), W_pp0(0:7,3), W_Nk0(0:7,3), W_pk0(0:7,3), W_kk0(0:7,3)
     &               ,W_Npp0(0:7,3), W_ppp0(0:7,3), W_Nkk0(0:7,3), W_Npk0(0:7,3)
     &               ,W_Nppp0(0:7,3)

c-
