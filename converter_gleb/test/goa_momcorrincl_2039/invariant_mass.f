c************************************************************************
	SUBROUTINE invariant_mass_2part(E1,px1,py1,pz1
     %                                 ,E2,px2,py2,pz2
     %                                 ,Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw)
c*************************************************************************
c+ Ivariant Mass and others variables of resulting "particle"  from 2 particles
c+    INPUT:
c+          p1(0-3)    = particle 1  4-momentum (Energy,px,py,pz)
c+          p2(0-3)    = particle 2  4-momentum (Energy,px,py,pz)
c+    OUTPUT
c+          M2w         = squared invariant mass
c+          Pw          = invariant momentum 
c+          Thw,PHw     = theta and phi of invariant particle
c+          pxm,pym,pzm = invariant 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw
	real E1,px1,py1,pz1,E2,px2,py2,pz2
	Ew  =  E1 +  E2 
        pxw =  px1 + px2 
        pyw =  py1 + py2 
        pzw =  pz1 + pz2 
	Pw  =  sqrt(Pxw**2+Pyw**2+Pzw**2)
        M2w =  Ew**2 - Pw**2
	call angles(0,Pxw,Pyw,Pzw,thw,phw)
	return
	end

c************************************************************************
	SUBROUTINE invariant_mass_3part(E1,px1,py1,pz1
     %                                 ,E2,px2,py2,pz2
     %                                 ,E3,px3,py3,pz3
     %                                 ,Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw)
c*************************************************************************
c+ Ivariant Mass and others variables of resulting "particle"  from 2 particles
c+    INPUT:
c+          p1(0-3)    = particle 1  4-momentum (Energy,px,py,pz)
c+          p2(0-3)    = particle 2  4-momentum (Energy,px,py,pz)
c+          p3(0-3)    = particle 3  4-momentum (Energy,px,py,pz)
c+    OUTPUT
c+          M2w         = squared invariant mass
c+          Pw          = invariant momentum 
c+          Thw,PHw     = theta and phi of invariant particle
c+          pxm,pym,pzm = invariant 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw
	real E1,px1,py1,pz1,E2,px2,py2,pz2,E3,px3,py3,pz3
	Ew  =  E1 +  E2 + E3 
        pxw =  px1 + px2 + px3
        pyw =  py1 + py2 + py3
        pzw =  pz1 + pz2 + pz3
	Pw  =  sqrt(Pxw**2+Pyw**2+Pzw**2)
        M2w =  Ew**2 - Pw**2
	call angles(0,Pxw,Pyw,Pzw,thw,phw)
	return
	end

c************************************************************************
	SUBROUTINE invariant_mass_4part(E1,px1,py1,pz1
     %                                 ,E2,px2,py2,pz2
     %                                 ,E3,px3,py3,pz3
     $                                 ,E4,px4,py4,pz4
     %                                 ,Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw)
c*************************************************************************
c+ Ivariant Mass and others variables of resulting "particle"  from 2 particles
c+    INPUT:
c+          p1(0-3)    = particle 1  4-momentum (Energy,px,py,pz)
c+          p2(0-3)    = particle 2  4-momentum (Energy,px,py,pz)
c+          p3(0-3)    = particle 3  4-momentum (Energy,px,py,pz)
c+          p4(0-3)    = particle 4  4-momentum (Energy,px,py,pz)
c+    OUTPUT
c+          M2w         = squared invariant mass
c+          Pw          = invariant momentum 
c+          Thw,PHw     = theta and phi of invariant particle
c+          pxm,pym,pzm = invariant 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw
	real E1,px1,py1,pz1,E2,px2,py2,pz2,E3,px3,py3,pz3,E4,px4,py4,pz4
	Ew  =  E1  +  E2 + E3  + E4
        pxw =  px1 + px2 + px3 + px4
        pyw =  py1 + py2 + py3 + py4
        pzw =  pz1 + pz2 + pz3 + pz4
	Pw  =  sqrt(Pxw**2+Pyw**2+Pzw**2)
        M2w =  Ew**2 - Pw**2
	call angles(0,Pxw,Pyw,Pzw,thw,phw)
	return
	end

c************************************************************************
	SUBROUTINE invariant_mass_5part(E1,px1,py1,pz1
     %                                 ,E2,px2,py2,pz2
     %                                 ,E3,px3,py3,pz3
     $                                 ,E4,px4,py4,pz4
     $                                 ,E5,px5,py5,pz5
     %                                 ,Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw)
c*************************************************************************
c+ Ivariant Mass and others variables of resulting "particle"  from 2 particles
c+    INPUT:
c+          p1(0-3)    = particle 1  4-momentum (Energy,px,py,pz)
c+          p2(0-3)    = particle 2  4-momentum (Energy,px,py,pz)
c+          p3(0-3)    = particle 3  4-momentum (Energy,px,py,pz)
c+          p4(0-3)    = particle 4  4-momentum (Energy,px,py,pz)
c+          p5(0-3)    = particle 5  4-momentum (Energy,px,py,pz)
c+    OUTPUT
c+          M2w         = squared invariant mass
c+          Pw          = invariant momentum 
c+          Thw,PHw     = theta and phi of invariant particle
c+          pxm,pym,pzm = invariant 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Ew,pxw,pyw,pzw,M2w,Pw,THw,PHw
	real E1,px1,py1,pz1,E2,px2,py2,pz2,E3,px3,py3,pz3,E4,px4,py4,pz4
     %  ,E5,px5,py5,pz5     
	Ew  =  E1  +  E2 + E3  + E4  + E5
        pxw =  px1 + px2 + px3 + px4 + px5
        pyw =  py1 + py2 + py3 + py4 + py5
        pzw =  pz1 + pz2 + pz3 + pz4 + pz5
	Pw  =  sqrt(Pxw**2+Pyw**2+Pzw**2)
        M2w =  Ew**2 - Pw**2
	call angles(0,Pxw,Pyw,Pzw,thw,phw)
	return
	end
