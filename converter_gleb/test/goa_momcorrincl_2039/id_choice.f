c*************************************************
	SUBROUTINE id_choice
c*************************************************
c  Fill vector ID_HADR and MASS_HADR using:
c  choice = 1 --> _OUR
c  choice = 2 --> _SEB
c*************************************************
      IMPLICIT none
      include "ntpl_goa.inc"
      integer i
		      	
      do i=1,10
	ID_hadr(i)  = 0
	mass_hadr(i)= -1000
	beta_hadr(i)= -1000
      enddo

c  HADRon is to be understood as a charged particle 
c  not first one in EVNT (i.e. not the trigger electron)
c  with SCstat_EVNT(j).gt.0
c  otherwise OUR=EVNT
	      	
      if(choice) then  ! OUR id,m,b
	DO i = 2,nEVNT
	  id_hadr(i)  = id_OUR(i)
	  mass_hadr(i)= mass_OUR(i)  ! it cames already squared
	  beta_hadr(i)= beta_OUR(i)
        ENDDO
      else
	DO i = 2,nEVNT ! SEB id,m,b
	  ID_hadr(i)  = ID_EVNT(i)
	  mass_hadr(i)= mass_EVNT(i) ! it cames already squared
	  beta_hadr(i)= beta_EVNT(i)
	ENDDO
      endif
	      	
      RETURN
      END      	
	      	
	      	
	      	
	      	
