c*******************************************************************************
	subroutine angles(iflag,v1,v2,v3,theta,phi)
c**********	
c+ INPUT	
c++	iflag = 0 calculating th and phi in degrees
c++	      = 1 calculating th and phi in rad
c++	v1,v2,v3 = input vector components
c+ OUTPUT
c++     theta, phi = calcalated angles      
c*******************************************************************************

	implicit none
	integer iflag
	real v1,v2,v3,v_m,theta,phi,php1

	v_m = sqrt(v1**2+v2**2+v3**2)

c-----------------------------------------
c	angles in the new frame (degrees)
c-----------------------------------------

	if(iflag.eq.0) then

	  if(v_m.gt.0.) then

	    theta = acos((v3/v_m))*180./3.14159265

	    if (v1.eq.0..and.v2.eq.0.) phi = 0

	    if (v1.gt.0..and.v2.ge.0.) then
 	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265
	    endif

	    if (v1.eq.0..and.v2.gt.0.) phi = 90.

	    if (v1.lt.0..and.v2.ge.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 180.
	    endif
	
	    if (v1.lt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 180.
	    endif

	    if (v1.eq.0..and.v2.lt.0.) phi = 270.

	    if (v1.gt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1)*180./3.14159265 + 360.
	    endif

	  else

	    theta = 0.
	    phi = 0.

	  endif


	endif

c-----------------------------------------
c	angles in the new frame (radians)
c-----------------------------------------

	if(iflag.eq.1) then
	

	  if(v_m.gt.0.) then

	    theta = acos(v3/v_m)

	    if (v1.eq.0..and.v2.eq.0.) phi = 0.

	    if (v1.gt.0..and.v2.ge.0.) then
	    php1 = v2/v1
	    phi = atan(php1)
	    endif

	    if (v1.eq.0..and.v2.gt.0.) phi = 1.5708

	    if (v1.lt.0..and.v2.ge.0.) then
	    php1 = v2/v1
 	    phi = atan(php1) + 3.14159265
	    endif
	
	    if (v1.lt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1) + 3.14159265
	    endif

	    if (v1.eq.0..and.v2.lt.0.) phi = 4.7124

	    if (v1.gt.0..and.v2.lt.0.) then
	    php1 = v2/v1
	    phi = atan(php1) + 6.2831852
	    endif

	  else

	    theta = 0.
	    phi = 0.

	  endif


	endif
c-------------------------------------------------------------------------------


	return

	end



c------------------------------------------------------------------------------
c       *************************
        function sector(angle)
c       ******c       Modified to be in accordance with PHI range from
c       -180 - +180;
c       
	real angle
        integer sector

        if (angle.le.30.or.angle.ge.330)   sector = 1
        if (angle.ge.30.and.angle.le.90)   sector = 2
        if (angle.ge.90.and.angle.le.150)  sector = 3
        if (angle.ge.150.and.angle.le.210) sector = 4
        if (angle.ge.210.and.angle.le.270) sector = 5
        if (angle.ge.270.and.angle.le.330) sector = 6

	return
	end


c ******************************************************************
c      sind,cosd,tand - functions with an argument being in radians
c      They should be intrinsic but in LINUX they are not 
c
c ******************************************************************
c      function sind(angle)
c      real angle,sind
c      sind = sin(angle*3.14159/180.)
c      return
c      end
c      function cosd()
c      real angle,cosd
c      cosd = cos(angle*3.14159/180.)
c      return
c      end
c      function tand()
c      real angle,tand
c      tand = tan(angle*3.14159/180.)
c      return
c      end

c ******************************************************************


	subroutine opposite_sec(mc_f)
	implicit none
        include "ntpl_goa.inc"
	integer j,jj,sec(6),sector,ST_all,ST_pair(3),mc_f,g6_max_TOF
	REAL*4   t1,t2,a1,a2,adc_thr_max_SC(50),ot(4,6),tau,mean_before,sigma
	REAL*4    open_opposite_trig, open_opposite_trig_1(4),open_opposite_trig_2(4),open_opposite_trig_3(4)
        INTEGER sec_coin_1,sec_coin_2,sec_coin_3
        INTEGER*4 j_SCPB,j_SC
        INTEGER*4 ll,sec_thr(4,6),sec_coin_1_thr(4),sec_coin_2_thr(4),sec_coin_3_thr(4)
	logical dead_g6a,dead_g6b
c	type *, sec_coin_ST, Egamma, mc_f

	If(mc_f.eq.1) then ! only generated MC data
c+ This part overwrite sec_coin_ST only for MC generated data
c+ assuming that ST eff = 1, it loops over EVNT particles (=MC generated)
c+ initialization
c	sec_coin = 0
	 do j=1,6
	   sec(j) = 0	
	 enddo
	 ST_all = 0
	 do j=1,3
	  ST_pair(j) = 0
	 enddo
c++
c+ Loop on charged particles (only hadrons)
	 do j = 2,nevnt
	  if(Charge_EVNT(j).ne.0) sec(sector(ph_EVNT(j))) = sec(sector(ph_EVNT(j))) + 1
	 enddo
c- Check for opposite sector coincidence
c	 if(sec(1).ne.0.and.sec(4).ne.0) sec_coin  = 1
c	 if(sec(2).ne.0.and.sec(5).ne.0) sec_coin  = 1
c	 if(sec(3).ne.0.and.sec(6).ne.0) sec_coin  = 1
c--

c++ Check for ST pairs hits (assuming ST_1=sector1-2; ST_2=sector3-4;ST_3=sector5-6;

	 if (sec(1).ne.0.or.sec(2).ne.0) ST_pair(1) =  1
	 if (sec(3).ne.0.or.sec(4).ne.0) ST_pair(2) =  1
	 if (sec(5).ne.0.or.sec(6).ne.0) ST_pair(3) =  1
	 ST_all = ST_pair(1) + ST_pair(2) + ST_pair(3)
c+ !!!! USING SEC_COIN TO HAVE ST PAIR INFO (OVERWRITING ON ORIGINAL SEC_COIN info)!!!!
	 sec_coin_ST = ST_all
c--	 

c++ JUST USING CEB SCR-based opposite sector coincidence
	
c	sec_coin = sec_coin_SCR


	else ! REC_MC/REAL_DATA
c Finding opposite sector coincidence
c use the same variable sec_coin_SCR as in SCR based coincidence
	sec_coin_SC = 0
        sec_coin_1 = 0
        sec_coin_2 = 0
        sec_coin_3 = 0
	do jj=1,6
	 sec(jj) = 0
	 do ll = 1,4
          sec_coin_SC_a(ll) = 0
	  sec_thr(ll,jj) = 0
	  ot(ll,jj) = 0
          sec_coin_SC_thr_a(ll) = 0
	  sec_coin_1_thr(ll) = 0
	  sec_coin_2_thr(ll) = 0
	  sec_coin_3_thr(ll) = 0
	  open_opposite_trig_1(ll)=0
	  open_opposite_trig_2(ll)=0
	  open_opposite_trig_3(ll)=0
	  trig_SC(ll)=0
	 enddo
	enddo
c+ Loop on SC hits 
	do jj = 1,nSC
c++ Clearing dead SC or dead pretrigger signals (according Buring note) 25/11/99
c	   if(
c     %    (sector_sc(jj).eq.1.and.
c     %      (id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.33.or.id_SC(jj).eq.34))
c     %.or.(sector_sc(jj).eq.2.and.
c     %      (id_SC(jj).eq.18.or.id_SC(jj).eq.24.or.id_SC(jj).eq.27.or.id_SC(jj).eq.29.or.id_SC(jj).eq.35))
c     %.or.(sector_sc(jj).eq.3.and.
c     %      (id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.36))
c     %.or.(sector_sc(jj).eq.4.and.
c     %      (id_SC(jj).eq.14.or.id_SC(jj).eq.21.or.id_SC(jj).eq.25))
c     %.or.(sector_sc(jj).eq.5.and.
c     %      (id_SC(jj).eq.15.or.id_SC(jj).eq.16.or.id_SC(jj).eq.29.or
c     %           .id_SC(jj).eq.31.or.id_SC(jj).eq.33.or.id_SC(jj).eq.37))
c     % .or.(sector_sc(jj).eq.6.and.
c     %      (id_SC(jj).eq.3.or.id_SC(jj).eq.25.or.id_SC(jj).eq.35))
c     %        ) then
c	    tdcr_SC(jj) = 0
c	    tdcl_SC(jj) = 0
c	    adcr_SC(jj) = 0
c	    adcl_SC(jj) = 0
c	   endif
c+++ Only dead channels
c	   if(
c     %    (sector_sc(jj).eq.1.and.
c     %      (id_SC(jj).eq.16))
c     %.or.(sector_sc(jj).eq.3.and.
c     %      (id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.36))
c     %.or.(sector_sc(jj).eq.4.and.
c     %      (id_SC(jj).eq.21))
c     %.or.(sector_sc(jj).eq.5.and.
c     %      (id_SC(jj).eq.16))
c     % .or.(sector_sc(jj).eq.6.and.
c     %      (id_SC(jj).eq.35))
c     %        ) then
c	      do ll=2,nEVNT
c	       j_SCPB = SCstat_EVNT(ll)
c               if(j_SCPB.gt.0) then
c                j_SC = SCstat_SCPB(j_SCPB)     ! Pointer SCPB->SC
c                if (j_SC.eq.jj) id_hadr(ll) = -971 ! killing events in dead channels
c	       endif	
c 	      enddo
c	      endif

c+++ Only dead channel
	   dead_g6a=
     %    (sector_sc(jj).eq.1.and.
     %      (id_SC(jj).eq.12.or.id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.23))
     %.or.(sector_sc(jj).eq.2.and.
     %      (id_SC(jj).eq.18.or.id_SC(jj).eq.22.or.id_SC(jj).eq.27.))
     %.or.(sector_sc(jj).eq.3.and.
     %      (id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.23.or.id_SC(jj).eq.36))
     %.or.(sector_sc(jj).eq.4.and.
     %      (id_SC(jj).eq.14.or.id_SC(jj).eq.16.or.id_SC(jj).eq.21.or.id_SC(jj).eq.25))
     %.or.(sector_sc(jj).eq.5.and.
     %      (id_SC(jj).eq.12.or.id_SC(jj).eq.15.or.id_SC(jj).eq.16.or.id_SC(jj).eq.29.or.id_SC(jj).eq.31))
     % .or.(sector_sc(jj).eq.6.and.
     %      (id_SC(jj).eq.13.or.id_SC(jj).eq.25.or.id_SC(jj).eq.35))
 
	   dead_g6b=
c     %    (sector_sc(jj).eq.1.and.
c     %      (id_SC(jj).eq.16.or.id_SC(jj).eq.34))
c     %.or.(sector_sc(jj).eq.2.and.
c     %      (id_SC(jj).eq.30))
c     %.or.
     %(sector_sc(jj).eq.3.and.
     %      (id_SC(jj).eq.18.or.id_SC(jj).eq.23.or.id_SC(jj).eq.37.or.id_SC(jj).eq.38
     %                      .or.id_SC(jj).eq.42.or.id_SC(jj).eq.43.or.id_SC(jj).eq.48))
     %.or.(sector_sc(jj).eq.4.and.
     %      (id_SC(jj).eq.21.or.id_SC(jj).eq.25))
     %.or.(sector_sc(jj).eq.5.and.
     %      (id_SC(jj).eq.2.or.id_SC(jj).eq.16.or.id_SC(jj).eq.47))
c     %.or.(sector_sc(jj).eq.6.and.
c     %      (id_SC(jj).eq.13.or.id_SC(jj).eq.25.or.id_SC(jj).eq.35))


        	   if( dead_g6b) then
	      do ll=2,nEVNT
	       j_SCPB = SCstat_EVNT(ll)
               if(j_SCPB.gt.0) then
                j_SC = SCstat_SCPB(j_SCPB)     ! Pointer SCPB->SC
                if (j_SC.eq.jj) id_hadr(ll) = -971 ! killing events in dead channels
	       endif	
 	      enddo
	      endif

c+++ According to the survay
c	   if(
c     %    (sector_sc(jj).eq.1.and.
c     %      (id_SC(jj).eq.5.or.id_SC(jj).eq.6.or.id_SC(jj).eq.14.or.id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.23))
c     %.or.(sector_sc(jj).eq.2.and.
c     %      (id_SC(jj).eq.18.or.id_SC(jj).eq.22.or.id_SC(jj).eq.24.or.id_SC(jj).eq.27.))
c     %.or.(sector_sc(jj).eq.3.and.
c     %      (id_SC(jj).eq.5.or.id_SC(jj).eq.7.or.id_SC(jj).eq.16.or.id_SC(jj).eq.18.or.id_SC(jj).eq.22.or
c     %          .id_SC(jj).eq.23.or.id_SC(jj).eq.36))
c     %.or.(sector_sc(jj).eq.4.and.
c     %      (id_SC(jj).eq.12.or.id_SC(jj).eq.14.or.id_SC(jj).eq.16.or.id_SC(jj).eq.21.or.id_SC(jj).eq.25.or.id_SC(jj).eq.8))
c     %.or.(sector_sc(jj).eq.5.and.
c     %      (id_SC(jj).eq.15.or.id_SC(jj).eq.16.or.id_SC(jj).eq.25.or.id_SC(jj).eq.27.or.id_SC(jj).eq.29.or.id_SC(jj).eq.31))
c     %.or.(sector_sc(jj).eq.6.and.
c     %      (id_SC(jj).eq.13.or.id_SC(jj).eq.24.or.id_SC(jj).eq.25.or
c     %           .id_SC(jj).eq.35))
c     %        ) then
c	    tdcr_SC(jj) = 0
c	    tdcl_SC(jj) = 0
c	    adcr_SC(jj) = 0
c	    adcl_SC(jj) = 0
c	   endif
	     
c--
c++ Finding threshold condition
        if (tdcl_SC(jj).gt.tdcr_SC(jj)) then
          t1 = tdcr_SC(jj)
          t2 = tdcl_SC(jj)
          a1 = adcr_SC(jj) - ped_SC_R((sector_sc(jj)-1)*48+id_sc(jj))
          a2 = adcl_SC(jj) - ped_SC_L((sector_sc(jj)-1)*48+id_sc(jj))
c	   write(*,*) sector_sc(jj),id_sc(jj),(sector_sc(jj)-1)*48+id_sc(jj),
c     %       ped_SC_R((sector_sc(jj)-1)*48+id_sc(jj)),ped_SC_L((sector_sc(jj)-1)*48+id_sc(jj)),a1,a2
	  
        else
          t1 = tdcl_SC(jj)
          t2 = tdcr_SC(jj)
          a1 = adcl_SC(jj) - ped_SC_L((sector_sc(jj)-1)*48+id_sc(jj))
          a2 = adcr_SC(jj) - ped_SC_R((sector_sc(jj)-1)*48+id_sc(jj))
        endif
c++ checking a1 and a2 after pedestal subtraction
	if (a1.lt.0) a1=0
	if (a2.lt.0) a2=0

c--
c-- 1/13ns = time constatnt corresponding to 30 ns @ 1/10*Amplitude and 0.05 ns/ch
         adc_thr_max_SC(jj) = max(a1,a1*exp(-1./(13./0.05)*(t2-t1))+a2)    
c++ Tyring a new signal shape
	 tau=(t2-t1)*0.05
         adc_thr_max_SC(jj) = max(a1,a1*((1-0.5875)*
     %      exp(-tau**2/2/5.326**2)+0.5875*exp(-tau/15.034))+a2) 
c+++ Changing MC threshold: size and sigma
	 if(bit_mc.eq.2) then
c	    mean_before = 1.45*a1
c	    sigma = 7*sqrt(mean_before)
	    mean_before = a1
	    sigma = 2.5*sqrt(mean_before)
            CALL GAUSS(sigma,mean_before,a1)
c	    write (*,*) mean_before, a1
	    mean_before = a2
	    sigma = 2.5*sqrt(mean_before)
            CALL GAUSS(sigma,mean_before,a2)
            adc_thr_max_SC(jj) = max(a1,a1*((1-0.5875)*
     %       exp(-tau**2/2/5.326**2)+0.5875*exp(-tau/15.034))+a2) 
	 endif
c---	    
c-- Ending  threshold condition

ccc        if(tdcl_sc(jj).gt.0.and.tdcr_sc(jj).gt.0) then
cc        if(adcl_sc(jj).gt.0.and.adcr_sc(jj).gt.0) then
ccc Var sec_SC contains the open_trigger condition with a fixed threshold (thesame as sec_thr)
c        if(adc_thr_max_SC(jj).gt.850) then
c          sec(sector_SC(jj)) = sec(sector_SC(jj)) + 1
c	  if(id_SC(jj).le.33) ot(sector_SC(jj)) = ot(sector_SC(jj)) +1
c	 endif
cccccccccccccccccccccccccccccccccccc
ccccccc     RENAMING :24/11/99
ccccccc sec_coin_SC_a(1-4)     = open trigger w fix thr(TOF<=33 in g6a; <=38in g6b)   -> SEC_SC
ccccccc sec_coin_SC_thr_a(1-4) = opposite sector (any TOF)                            -> SEC_SC_thr
ccccccc trig_SC(1-4)           = open.and.opposite                                    -> TRIG_SC 
cccccccccccccccccccccccccccccccccccc
c g6a values: 1 =850 2 =900 3=900 4=750
c g6a values:  g6_max_TOF= 33

	 g6_max_TOF= 38


         if(adc_thr_max_SC(jj).gt.850) then ! TO BE COMPATIBLE WITH g6a ANALYSIS
	  if(id_SC(jj).le.g6_max_TOF) then
	   sec_coin_SC_a(1) = 1
 	   ot(1,sector_SC(jj)) = ot(1,sector_SC(jj)) + 1
	  endif
          sec_thr(1,sector_SC(jj)) = sec_thr(1,sector_SC(jj)) + 1
	 endif

         if(adc_thr_max_SC(jj).gt.150) then ! thr = 
	  if(id_SC(jj).le.g6_max_TOF) then
       	   sec_coin_SC_a(2) = 1
 	   ot(2,sector_SC(jj)) = ot(2,sector_SC(jj)) + 1
	  endif
          sec_thr(2,sector_SC(jj)) = sec_thr(2,sector_SC(jj)) + 1
	 endif

         if(adc_thr_max_SC(jj).gt.10) then ! thr = 
	  if(id_SC(jj).le.g6_max_TOF)then
           sec_coin_SC_a(3) = 1
 	   ot(3,sector_SC(jj)) = ot(3,sector_SC(jj)) + 1
	  endif
          sec_thr(3,sector_SC(jj)) = sec_thr(3,sector_SC(jj)) + 1
	 endif

         if(adc_thr_max_SC(jj).gt.5) then ! thr = 
	  if(id_SC(jj).le.g6_max_TOF) then
	   sec_coin_SC_a(4) = 1
 	   ot(4,sector_SC(jj)) = ot(4,sector_SC(jj)) + 1
	  endif
          sec_thr(4,sector_SC(jj)) = sec_thr(4,sector_SC(jj)) + 1
	 endif
ccc        endif
c+++ Open_async trigger condition (no threshold):
c	 if(adc_thr_max_SC(jj).gt.1500.and.id_SC(jj).le.32) open_trig=1
        enddo ! closing loop on nSC
c-
c+ Check for opposite sector coincidence
c	 if(sec(1).ne.0.and.sec(4).ne.0) then 
c          sec_coin_1  = 1
c	  if(ot(1).ne.0.or.ot(4).ne.0) open_opposite_trig=1
c	 endif
c	 if(sec(2).ne.0.and.sec(5).ne.0) then 
c          sec_coin_2  = 1
c	  if(ot(2).ne.0.or.ot(5).ne.0) open_opposite_trig=1
c	 endif
c	 if(sec(3).ne.0.and.sec(6).ne.0) then
c          sec_coin_3  = 1
c	  if(ot(3).ne.0.or.ot(6).ne.0) open_opposite_trig=1
c	 endif
c	sec_coin_SC = sec_coin_1 + sec_coin_2 + sec_coin_3 ! replaced by sec_coin_SC_a(..)
	 do ll=1,4
	  if(sec_thr(ll,1).ne.0.and.sec_thr(ll,4).ne.0)then
           sec_coin_1_thr(ll)= 1
	   if(ot(ll,1).ne.0.or.ot(ll,4).ne.0) open_opposite_trig_1(ll)=1
	  endif
	  if(sec_thr(ll,2).ne.0.and.sec_thr(ll,5).ne.0)then
           sec_coin_2_thr(ll)= 1
	   if(ot(ll,2).ne.0.or.ot(ll,5).ne.0) open_opposite_trig_2(ll)=1
	  endif
	  if(sec_thr(ll,3).ne.0.and.sec_thr(ll,6).ne.0)then
           sec_coin_3_thr(ll)= 1
	   if(ot(ll,3).ne.0.or.ot(ll,6).ne.0) open_opposite_trig_3(ll)=1
	  endif 
c+ New variable in ntuple 20 containing different threshold for opposite sec coincid
c+ it makes sense only for MC data because the scales are differents
         sec_coin_SC_thr_a(ll) = sec_coin_1_thr(ll) + sec_coin_2_thr(ll) + sec_coin_3_thr(ll) 
	 trig_SC(ll) = open_opposite_trig_1(ll) + open_opposite_trig_2(ll) + open_opposite_trig_3(ll)

	 enddo
c-
c	 sec_coin_SC_thr_a(1) = open_trig


c+ Looping over EVNT bank to find  Thr_trigSC_EVNT
	do j = 2,nEVNT
         Thr_trigSC_EVNT(j) = -1000.
         j_SCPB = SCstat_EVNT(j)        ! Pointer EVNT->SCPB bank
          if(j_SCPB.gt.0) then
            j_SC = SCstat_SCPB(j_SCPB)     ! Pointer SCPB->SC
            if (j_SC.gt.0) then
             Thr_trigSC_EVNT(j) = adc_thr_max_SC(j_SC)
             id_SC_EVNT(j) = id_SC(j_SC)
             sec_SC_EVNT(j) =sector_SC(j_SC)
             endif
            endif
         ENDDO ! next particle

	endif

c	type *,'dopo', sec_coin_ST, Egamma, mc_f

	return
	end

C*********************************
	SUBROUTINE GAUSS(S,AM,V)
C********************************
	IMPLICIT NONE
	REAL RRAN
	REAL A,AM,V,S
	INTEGER*4 IZ,I
	iz=337699231
	A=0.
	DO 50 I=1,12
50	A=A+rran()
	V=(A-6.0)*S+AM
	END    

	FUNCTION RRAN()
	REAL RRAN
	CALL RANLUX(RRAN,1)
	RETURN
	END
