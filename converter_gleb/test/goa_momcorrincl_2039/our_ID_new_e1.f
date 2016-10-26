
	subroutine our_ID_new_e1
c----------------------------------------------------------
c	"homemade" particle ID!
c----------------------------------------------------------
      implicit NONE	
      include "ntpl_goa.inc"

      real*4 c,tcorr
      INTEGER*4 sc_el,sc_h
      REAL*4 tof_e,tdc_e,tdc_h,l_h 
        data c,tcorr/30.0,0.00/	
	integer j
	REAL*4 x_ID,y_ID,pions_1,pions_2,protons_1,protons_2
c-
c+ The First particle is always assumed to be electron although it's not true
           jopa=-1000.
           mass_our(1) = mass_EVNT(1)
           beta_our(1) = beta_EVNT(1)   
	   id_our(1)   = id_EVNT(1)
           DO j=2,nEVNT
	     id_our(j)  = -979 ! init mass_our
             mass_our(j)= -979 !  init id_our
	   ENDDO


c+ Trigger time

	   sc_el = SCstat_EVNT(1)
	   IF(sc_el.eq.0) then ! if no e- sc time take SEB mass and PID for NEUTRAL particles
	    DO j=2,nEVNT
	     mass_our(j) = mass_EVNT(j)
             beta_our(j) = beta_EVNT(j)
	     IF(charge_evnt(j).eq.0) id_our(j) = id_EVNT(j)! our=evnt for Neutral particles   
	    ENDDO
	   RETURN
	   ENDIF

           tof_e=Path_SCPB(sc_el)/c
           tdc_e=Time_SCPB (sc_el)

c+ Looping on hadrons           
	DO j=2,nEVNT
	  
	  IF(charge_evnt(j).ne.0.and.SCstat_EVNT(j).gt.0) then
	   sc_h=SCstat_EVNT (j)     
           tdc_h=Time_SCPB (sc_h)
           l_h =Path_SCPB(sc_h)
           beta_our(j)=l_h/(tof_e - tdc_e + tdc_h + tcorr)/c
           mass_our(j)=((pmom_EVNT(j)/beta_our(j))**2)*(1. -(beta_our(j))**2)
	  ELSEIF(charge_evnt(j).eq.0) then
	   beta_our(j) = beta_EVNT(j)
           mass_our(j) = mass_EVNT(j)
	   id_our(j)   = id_EVNT(j)! our=evnt for NEUTRAL particles
	  ENDIF
	
	x_ID = pmom_EVNT(j)
	y_ID = Beta_our(j)


c	pions_1 = 1.1

c	pions_2 = (1+5*1.07*(x_ID-0.07))/(1+5*(x_ID-0.07))
c	pions_2 = pions_2*(x_ID-0.07)/sqrt((x_ID-0.07)**2+0.138**2)
c	pions_2 = pions_2 - 0.1

	pions_1 = 1.4

	pions_2 = (1+5*1.4*(x_ID-0.07))/(1+5*(x_ID-0.07))
	pions_2 = pions_2*(x_ID-0.07)/sqrt((x_ID-0.07)**2+0.138**2)
	pions_2 = pions_2 - 0.4

	protons_1 = x_ID/sqrt(x_ID**2+0.938**2) + 0.03
	protons_1 = protons_1*(1.2+0.92*x_ID)/(1+x_ID)

	protons_2 = x_ID/sqrt(x_ID**2+0.938**2) - 0.075
	protons_2 = protons_2*(0.9+1.06*x_ID)/(1+x_ID)

c----------------------------------------------------------
	if(charge_EVNT(j).eq.1) then
c----------------------------------------------------------	

	  if(y_ID.gt.pions_2.AND.y_ID.lt.pions_1) then
	  ID_our(j) = 211
	  jopa=(Time_SCPB(SCstat_EVNT (j)) -
     &     (Path_SCPB(SCstat_EVNT(j))/((x_ID)/(sqrt(x_ID**2+0.139**2))))) - joppa
	   elseif(y_ID.gt.protons_2.AND.y_ID.lt.protons_1) then
	  ID_our(j) = 2212
	  else
	  ID_our(j) = 0
	  endif

c----------------------------------------------------------
	elseif(charge_EVNT(j).eq.-1) then
c----------------------------------------------------------

	  if(y_ID.gt.pions_2.AND.y_ID.lt.pions_1) then
	  ID_our(j) = -211
	  else
	  ID_our(j) = 0
	  endif

c----------------------------------------------------------
	endif	
c----------------------------------------------------------

	ENDDO
	  

	RETURN
	END
