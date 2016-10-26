c------------------------------------------------------------------------ 
C	Kinematical processing from EVNT bank information
c------------------------------------------------------------------------ 

      SUBROUTINE kinemaEVNT
     
      implicit NONE
     
	include "ntpl_goa.inc"
	
c------------------------------------------------------
C	local variables
c------------------------------------------------------
      integer*4 i,j,k
      REAL*4    pmx,pmy,pmz,m_cx,m_cy,m_cz,Ema
      REAL*4    E1,p1,p1x,p1y,p1z,E2,p2,p2x,p2y,p2z,W2_EVNT
c------------------------------------------------------
      
	W_EVNT = -1000
	omega_EVNT = -1000
	Q2_EVNT = -1000
	Tp_el_EVNT = -1000
	thp_el_EVNT = -1000
	Mm2a_EVNT = -1000
	Mm2b_EVNT = -1000
	Mm2c_EVNT = -1000
	s1_EVNT = -1000
	s2_EVNT = -1000
	s3_EVNT = -1000
	p_miss_EVNT = -1000
	th_miss_EVNT = -1000
	ph_miss_EVNT = -1000
c	do i=1,10
c	 ID_hadr(i) = 0
c	 mass_hadr(i)= -1000
c	enddo
c------------------------------------------------------------------------ 
C	Choice of SEB or OUR particle ID for hadrons
c------------------------------------------------------------------------ 
c	if(choice) then
c	DO i = 2,nEVNT
c	ID_hadr(i) = ID_our(i)
c	mass_hadr(i) = mass_our(i)
c	ENDDO
c	else
c	DO i = 2,nEVNT
c	ID_hadr(i) = ID_EVNT(i)
c	mass_hadr(i) = mass_EVNT(i)
c	ENDDO
c	endif

c------------------------------------------------------------------------ 
C	Calculation of typical electron scattering kinematical quantities
c------------------------------------------------------------------------ 
c!!!! START !!!!  Marco Battaglieri:10/2/98
c+ Calculating proton variables in elastic hypotesis
	Tp_el_EVNT = E0**2*sin(th_evnt(1)*3.1416/2./180.)**2/(E0*sin(th_evnt(1)*3.1416/2./180.)**2+0.938/2.)
	thp_el_EVNT = (180./3.1416)*
     &	acos((1+0.938/E0)/sqrt(1.+2*0.938/E0+(0.938/E0)**2/sin(th_evnt(1)*3.1416/2./180.)**2))
c-
c!!!! END !!!!
	omega_EVNT = E0 - pmom_EVNT(1)
	Q2_EVNT = 4.*E0*pmom_EVNT(1)*(sin(th_EVNT(1)*3.1416/2./180.))**2
	W2_EVNT = 0.938**2+2.*0.938*omega_EVNT-Q2_EVNT
	if(W2_EVNT.ge.0) then
	W_EVNT = sqrt(W2_EVNT)
	else
	continue
	endif
c!!!! START !!!!  Marco Battaglieri:10/2/98
c+ Kinematics for 2 charged particles detected and 1 hadron identified
	IF (nEVNT.eq.2.and.id_hadr(2).ne.0) then

c++ Missing particle calculation
c+++ Energy
	Ema = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+mass_hadr(2)**2)
c+++ 3-Momentum
	pmx = -pmomx_EVNT(1)-pmomx_EVNT(2)
	pmy = -pmomy_EVNT(1)-pmomy_EVNT(2)
	pmz = E0-pmomz_EVNT(1)-pmomz_EVNT(2)
	p_miss_EVNT = sqrt(pmx**2+pmy**2+pmz**2)
c+++ Missing mass
	Mm2a_EVNT = Ema**2-p_miss_EVNT**2
		
	m_cx = pmx/p_miss_EVNT 
	m_cy = pmy/p_miss_EVNT
	m_cz = pmz/p_miss_EVNT

		if(m_CZ.gt.1.or.m_CZ.lt.-1) then
         	print *, 'EVNT theta angle undefined for missing particle'
         	continue
        	else
		th_miss_EVNT=(180./3.1416)*acos(m_CZ)
		endif

		if(m_CX.gt.1.or.m_CX.lt.-1.or.m_CY.gt.1.or.m_CY.lt.-1) then
         	print *, 'EVNT x or y component of missing momentum not well defined'
         	continue
        	else
		ph_miss_EVNT = (180./3.1416)*atan(m_CY/m_CX)
		if(m_CX.lt.0.and.m_CY.gt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.lt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.gt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 360.
		endif
c---
c--
c++ Invarant mass p pi0 or pi+ n

	  E1 = sqrt(pmom_EVNT(2)**2+mass_hadr(2)**2) + Ema
	  p1x = pmomx_EVNT(2) + pmx
	  p1y = pmomy_EVNT(2) + pmy
	  p1z = pmomz_EVNT(2) + pmz	
	  p1 = sqrt(p1x**2+p1y**2+p1z**2)
	  s1_EVNT = E1**2-p1**2
c-
c!!!! END !!!!
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 
	ELSE IF(nEVNT.eq.3) then
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

c------------------------------------------------------------------------ 
C	Calculation of missing quantities in the case (e-)(p pi+/-)
c------------------------------------------------------------------------ 
	if((charge_EVNT(2).eq.1.OR.charge_EVNT(2).eq.-1).AND.(ID_hadr(2).eq.2212.OR.ID_hadr(2).eq.211)) then
	 if((charge_EVNT(3).eq.1.OR.charge_EVNT(3).eq.-1).AND.(ID_hadr(3).eq.2212.OR.ID_hadr(3).eq.211)) then
c------------------------------------------------------------------------ 
	
	Ema = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+mass_hadr(2)**2)-
     &	      sqrt(pmom_EVNT(3)**2+mass_hadr(3)**2)

	pmx = -pmomx_EVNT(1)-pmomx_EVNT(2)-pmomx_EVNT(3)
	pmy = -pmomy_EVNT(1)-pmomy_EVNT(2)-pmomy_EVNT(3)
	pmz = E0-pmomz_EVNT(1)-pmomz_EVNT(2)-pmomz_EVNT(3)
	
	p_miss_EVNT = sqrt(pmx**2+pmy**2+pmz**2)
	Mm2a_EVNT = Ema**2-p_miss_EVNT**2
		
	m_cx = pmx/p_miss_EVNT
	m_cy = pmy/p_miss_EVNT
	m_cz = pmz/p_miss_EVNT

		if(m_CZ.gt.1.or.m_CZ.lt.-1) then
         	print *, 'EVNT theta angle undefined for missing particle'
         	continue
        	else
		th_miss_EVNT=(180./3.1416)*acos(m_CZ)
		endif

		if(m_CX.gt.1.or.m_CX.lt.-1.or.m_CY.gt.1.or.m_CY.lt.-1) then
         	print *, 'EVNT x or y component of missing momentum not well defined'
         	continue
        	else
		ph_miss_EVNT = (180./3.1416)*atan(m_CY/m_CX)
		if(m_CX.lt.0.and.m_CY.gt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.lt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.gt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 360.
		endif


c------------------------------------------------------------------------ 
C	Calculation of invariant masses in the case (e-)(p pi+/-)
C	i = proton, j = pi+, k = pi-
c------------------------------------------------------------------------ 

	  if(ID_hadr(2).eq.2212.AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.1) then
	  i = 2
	  j = 3	  
	  elseif(ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.1.AND.ID_hadr(3).eq.2212) then
	  i = 3
	  j = 2
	  else
	  goto 111
	  endif

	  E1 = sqrt(pmom_EVNT(i)**2+mass_hadr(i)**2) + sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2)
	  p1x = pmomx_EVNT(i) + pmomx_EVNT(j)
	  p1y = pmomy_EVNT(i) + pmomy_EVNT(j)
	  p1z = pmomz_EVNT(i) + pmomz_EVNT(j)	
	  p1 = sqrt(p1x**2+p1y**2+p1z**2)
	  s1_EVNT = E1**2-p1**2
	  E2 = sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2) + sqrt(p_miss_EVNT**2+0.138**2)
	  p2x = pmomx_EVNT(j) + pmx
	  p2y = pmomy_EVNT(j) + pmy
	  p2z = pmomz_EVNT(j) + pmy	
	  p2 = sqrt(p2x**2+p2y**2+p2z**2)
	  s2_EVNT = E2**2-p2**2

	  goto 999

111	  if(ID_hadr(2).eq.2212.AND.ID_hadr(3).eq.211.and.charge_EVNT(3).eq.-1) then
	  i = 2
	  k = 3
	  elseif(ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.-1.AND.ID_hadr(3).eq.2212) then
	  i = 3
	  k = 2
	  else
	  goto 222
	  endif

	  E1 = sqrt(pmom_EVNT(i)**2+mass_hadr(i)**2) + sqrt(p_miss_EVNT**2+0.138**2)
	  p1x = pmomx_EVNT(i) + pmx
	  p1y = pmomy_EVNT(i) + pmy
	  p1z = pmomz_EVNT(i) + pmy	
	  p1 = sqrt(p1x**2+p1y**2+p1z**2)
	  s1_EVNT = E1**2-p1**2
	  E2 = sqrt(pmom_EVNT(k)**2+mass_hadr(k)**2) + sqrt(p_miss_EVNT**2+0.138**2)
	  p2x = pmomx_EVNT(k) + pmx
	  p2y = pmomy_EVNT(k) + pmy
	  p2z = pmomz_EVNT(k) + pmy	
	  p2 = sqrt(p2x**2+p2y**2+p2z**2)
	  s2_EVNT = E2**2-p2**2
	  
	  goto 999

222	  if(ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.1.AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.-1) then
	  j = 2
	  k = 3
	  elseif(ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.-1.AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.1) then
	  j = 3
	  k = 2
	  else
	  goto 999
	  endif
	  
	  E1 = sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2) + sqrt(p_miss_EVNT**2+0.938**2)
	  p1x = pmomx_EVNT(j) + pmx
	  p1y = pmomy_EVNT(j) + pmy
	  p1z = pmomz_EVNT(j) + pmy	
	  p1 = sqrt(p1x**2+p1y**2+p1z**2)
	  s1_EVNT = E1**2-p1**2
	  E2 = sqrt(pmom_EVNT(k)**2+mass_hadr(k)**2) + sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2)
	  p2x = pmomx_EVNT(k) + pmomx_EVNT(j)
	  p2y = pmomy_EVNT(k) + pmomy_EVNT(j)
	  p2z = pmomz_EVNT(k) + pmomz_EVNT(j)	
	  p2 = sqrt(p2x**2+p2y**2+p2z**2)
	  s2_EVNT = E2**2-p2**2

c------------------------------------------------------------------------ 
	 endif
	endif		

c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 
	ELSEIF(nEVNT.eq.4) then
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

c------------------------------------------------------------------------ 
C	Calculation of missing quantities in the case (e-)(p pi+/- pi+/-)
c------------------------------------------------------------------------ 
	if((charge_EVNT(2).eq.1.OR.charge_EVNT(2).eq.-1).AND.(ID_hadr(2).eq.2212.OR.ID_hadr(2).eq.211)) then
	 if((charge_EVNT(3).eq.1.OR.charge_EVNT(3).eq.-1).AND.(ID_hadr(3).eq.2212.OR.ID_hadr(3).eq.211)) then
	  if((charge_EVNT(4).eq.1.OR.charge_EVNT(4).eq.-1).AND.(ID_hadr(4).eq.2212.OR.ID_hadr(4).eq.211)) then
c------------------------------------------------------------------------ 
	
	Ema = E0+0.938-pmom_EVNT(1) - sqrt(pmom_EVNT(2)**2+mass_hadr(2)**2) -
     &	      sqrt(pmom_EVNT(3)**2+mass_hadr(3)**2) - sqrt(pmom_EVNT(4)**2+mass_hadr(4)**2)

	pmx = -pmomx_EVNT(1)-pmomx_EVNT(2)-pmomx_EVNT(3)-pmomx_EVNT(4)
	pmy = -pmomy_EVNT(1)-pmomy_EVNT(2)-pmomy_EVNT(3)-pmomy_EVNT(4)
	pmz = E0-pmomz_EVNT(1)-pmomz_EVNT(2)-pmomz_EVNT(3)-pmomz_EVNT(4)
	
	p_miss_EVNT = sqrt(pmx**2+pmy**2+pmz**2)
	Mm2a_EVNT = Ema**2-p_miss_EVNT**2
		
	m_cx = pmx/p_miss_EVNT
	m_cy = pmy/p_miss_EVNT
	m_cz = pmz/p_miss_EVNT

		if(m_CZ.gt.1.or.m_CZ.lt.-1) then
         	print *, 'EVNT theta angle undefined for missing particle'
         	continue
        	else
		th_miss_EVNT=(180./3.1416)*acos(m_CZ)
		endif

		if(m_CX.gt.1.or.m_CX.lt.-1.or.m_CY.gt.1.or.m_CY.lt.-1) then
         	print *, 'EVNT x or y component of missing momentum not well defined'
         	continue
        	else
		ph_miss_EVNT = (180./3.1416)*atan(m_CY/m_CX)
		if(m_CX.lt.0.and.m_CY.gt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.lt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
		if(m_CX.gt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 360.
		endif


c------------------------------------------------------------------------ 
C	Calculation of invariant masses in the case (e-)(p pi+/- pi+/-)
c------------------------------------------------------------------------ 

	  if(ID_hadr(2).eq.2212.AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.1
     &	     .AND.ID_hadr(4).eq.211.AND.charge_EVNT(4).eq.-1) then
	  i = 2
	  j = 3
	  k = 4
	  elseif(ID_hadr(2).eq.2212.AND.ID_hadr(4).eq.211.AND.charge_EVNT(4).eq.1
     &	     .AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.-1) then
	  i = 2
	  j = 4
	  k = 3
	  elseif(ID_hadr(3).eq.2212.AND.ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.1
     &	     .AND.ID_hadr(4).eq.211.AND.charge_EVNT(4).eq.-1) then
	  i = 3
	  j = 2
	  k = 4
	  elseif(ID_hadr(4).eq.2212.AND.ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.1
     &	     .AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.-1) then
	  i = 4
	  j = 2
	  k = 3
	  elseif(ID_hadr(3).eq.2212.AND.ID_hadr(4).eq.211.AND.charge_EVNT(4).eq.1
     &	     .AND.ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.-1) then
	  i = 3
	  j = 4
	  k = 2
	  elseif(ID_hadr(4).eq.2212.AND.ID_hadr(3).eq.211.AND.charge_EVNT(3).eq.1
     &	     .AND.ID_hadr(2).eq.211.AND.charge_EVNT(2).eq.-1) then
	  i = 4
	  j = 3
	  k = 2
	  else
	  goto 999
	  endif

	  E1 = sqrt(pmom_EVNT(i)**2+mass_hadr(i)**2) + sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2)
	  p1x = pmomx_EVNT(i) + pmomx_EVNT(j)
	  p1y = pmomy_EVNT(i) + pmomy_EVNT(j)
	  p1z = pmomz_EVNT(i) + pmomz_EVNT(j)	
	  p1 = sqrt(p1x**2+p1y**2+p1z**2)
	  s1_EVNT = E1**2-p1**2
	  E2 = sqrt(pmom_EVNT(k)**2+mass_hadr(k)**2) + sqrt(pmom_EVNT(j)**2+mass_hadr(j)**2)
	  p2x = pmomx_EVNT(k) + pmomx_EVNT(j)
	  p2y = pmomy_EVNT(k) + pmomy_EVNT(j)
	  p2z = pmomz_EVNT(k) + pmomz_EVNT(j)	
	  p2 = sqrt(p2x**2+p2y**2+p2z**2)
	  s2_EVNT = E2**2-p2**2

c------------------------------------------------------------------------ 
	  endif
	 endif
	endif		

c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 
	ENDIF
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

c------------------------------------------------------------------------ 
c	IF(nEVNT.eq.3) then
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

c------------------------------------------------------------------------ 
C	Request that second particle in the bank has positive charge
c------------------------------------------------------------------------ 
c	if(charge_EVNT(2).ne.1) return
c------------------------------------------------------------------------ 
C	Request that third particle in the bank is charged
c------------------------------------------------------------------------ 
c	if(charge_EVNT(3).ne.-1.AND.charge_EVNT(3).ne.1) return
	
c------------------------------------------------------------------------ 

c	if(charge_EVNT(3).eq.-1) then
c------------------------------------------------------------------------ 
C	Calculation of missing quantities in the case (e-)(+)(-)
C	First hypothesis is that second particle is a proton
C	and third is a pi-
C	Second hypothesis is that second particle is a pi+
C	and third is a pi-
c------------------------------------------------------------------------ 
	
c	Ema = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+0.938**2)-
c     &	      sqrt(pmom_EVNT(3)**2+0.138**2)
c	Emb = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+0.138**2)-
c     &	      sqrt(pmom_EVNT(3)**2+0.138**2)

c	pmx = -pmom_EVNT(1)*CX_EVNT(1)-pmom_EVNT(2)*CX_EVNT(2)-pmom_EVNT(3)*CX_EVNT(3)
c	pmy = -pmom_EVNT(1)*CY_EVNT(1)-pmom_EVNT(2)*CY_EVNT(2)-pmom_EVNT(3)*CY_EVNT(3)
c	pmz = E0-pmom_EVNT(1)*CZ_EVNT(1)-pmom_EVNT(2)*CZ_EVNT(2)-pmom_EVNT(3)*CZ_EVNT(3)
	
c	p_miss_EVNT = sqrt(pmx**2+pmy**2+pmz**2)
c	Mm2a_EVNT = Ema**2-p_miss_EVNT**2
c	Mm2b_EVNT = Emb**2-p_miss_EVNT**2
		
c	m_cx = pmx/p_miss_EVNT
c	m_cy = pmy/p_miss_EVNT
c	m_cz = pmz/p_miss_EVNT

c		if(m_CZ.gt.1.or.m_CZ.lt.-1) then
c         	print *, 'EVNT theta angle undefined for missing particle'
c         	continue
c        	else
c		th_miss_EVNT=(180./3.1416)*acos(m_CZ)
c		endif

c		if(m_CX.gt.1.or.m_CX.lt.-1.or.m_CY.gt.1.or.m_CY.lt.-1) then
c         	print *, 'EVNT x or y component of missing momentum not well defined'
c         	continue
c        	else
c		ph_miss_EVNT = (180./3.1416)*atan(m_CY/m_CX)
c		if(m_CX.lt.0.and.m_CY.gt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
c		if(m_CX.lt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
c		if(m_CX.gt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 360.
c		endif
		
c------------------------------------------------------------------------ 
c	elseif(charge_EVNT(3).eq.1) then	
		
c------------------------------------------------------------------------ 
C	Calculation of missing quantities in the case (e-)(+)(+)
C	First hypothesis is that second particle is a proton
C	and third is a pi+
C	Second hypothesis is that second particle is a pi+
C	and third is a proton
c------------------------------------------------------------------------ 
	
c	Ema = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+0.938**2)-
c     &	      sqrt(pmom_EVNT(3)**2+0.138**2)
c	Emb = E0+0.938-pmom_EVNT(1)-sqrt(pmom_EVNT(2)**2+0.138**2)-
c     &	      sqrt(pmom_EVNT(3)**2+0.938**2)

c	pmx = -pmom_EVNT(1)*CX_EVNT(1)-pmom_EVNT(2)*CX_EVNT(2)-pmom_EVNT(3)*CX_EVNT(3)
c	pmy = -pmom_EVNT(1)*CY_EVNT(1)-pmom_EVNT(2)*CY_EVNT(2)-pmom_EVNT(3)*CY_EVNT(3)
c	pmz = E0-pmom_EVNT(1)*CZ_EVNT(1)-pmom_EVNT(2)*CZ_EVNT(2)-pmom_EVNT(3)*CZ_EVNT(3)
	
c	p_miss_EVNT = sqrt(pmx**2+pmy**2+pmz**2)
c	Mm2a_EVNT = Ema**2-p_miss_EVNT**2
c	Mm2b_EVNT = Emb**2-p_miss_EVNT**2

c	m_cx = pmx/p_miss_EVNT
c	m_cy = pmy/p_miss_EVNT
c	m_cz = pmz/p_miss_EVNT

c		if(m_CZ.gt.1.or.m_CZ.lt.-1) return
c		th_miss_EVNT=(180./3.1416)*acos(m_CZ)

c		if(m_CX.gt.1.or.m_CX.lt.-1) return
c		if(m_CY.gt.1.or.m_CY.lt.-1) return
		
c		ph_miss_EVNT = (180./3.1416)*atan(m_CY/m_CX)
c		if(m_CX.lt.0.and.m_CY.gt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
c		if(m_CX.lt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 180.
c		if(m_CX.gt.0.and.m_CY.lt.0) ph_miss_EVNT=ph_miss_EVNT + 360.

c	endif	

c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 
c	ENDIF
c------------------------------------------------------------------------ 
c------------------------------------------------------------------------ 

c------------------------------------------------------------------------ 

999	continue

	return
	end
