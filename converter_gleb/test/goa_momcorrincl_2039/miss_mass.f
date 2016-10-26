c************************************************************************
	SUBROUTINE miss_mass_2part(Mtarget,E_beam
     &                            ,Eel,pxel,pyel,pzel
     &                            ,E1,px1,py1,pz1
     @                            ,Em,pxm,pym,pzm,M2m,Pm,THm,PHm)
c*************************************************************************
c+ 4-momentum and mass of missing particle calcutation.
c+ The "missing" particle is defined respect to electron and particle 1
c+    INPUT:
c+          Mtarget    = target mass
c+          Ebeam      = Beam energy 
c+          Eel,Pxel...= scattered electron 4-momentum (Energy,px,py,pz)
c+          E1,Px1...  = 4-momentum of particle 1
c+    OUTPUT
c+          M2m         = squared missing mass
c+          Pm          = missing momentum
c+          Thm,PHm     = theta and phi of missing particle
c+          pxm,pym,pzm = missing 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Mtarget,E_beam
	real Eel,pxel,pyel,pzel,M2m,Pm,THm,PHm,E1,px1,py1,pz1
	real Em,Pxm,Pym,Pzm
	Em  =  E_beam + Mtarget - Eel  - E1
        Pxm =                   - pxel - px1 
        Pym =                   - pyel - py1
        Pzm =  E_beam           - pzel - pz1
	Pm  = sqrt(Pxm**2+Pym**2+Pzm**2)
        M2m = Em**2 - Pm**2
	call angles(0,Pxm,Pym,Pzm,thm,phm)
	return
	end


c************************************************************************
	SUBROUTINE miss_mass_3part(Mtarget,Ebeam
     $                            ,Eel,pxel,pyel,pzel
     $                            ,E1,px1,py1,pz1
     $                            ,E2,px2,py2,pz2
     $                            ,Em,pxm,pym,pzm,M2m,Pm,THm,PHm)
c*************************************************************************
c+ 4-momentum and mass of missing particle calcutation.
c+ The "missing" particle is defined respect to electron particle 1 and particle 2
c+    INPUT:
c+          Mtarget    = target mass
c+          Ebeam      = Beam energy 
c+          Eel,Pxel...= scattered electron 4-momentum (Energy,px,py,pz)
c+          E1,Px1...  = 4-momentum of particle 1
c+          E2,Px2...  = 4-momentum of particle 2
c+    OUTPUT
c+          M2m         = squared missing mass
c+          Pm          = missing momentum
c+          Thm,PHm     = theta and phi of missing particle
c+          pxm,pym,pzm = missing 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Mtarget,Ebeam
	real Eel,pxel,pyel,pzel,M2m,Pm,THm,PHm,E1,px1,py1,pz1,E2,px2,py2,pz2
	real Em,Pxm,Pym,Pzm
	Em  =  Ebeam + Mtarget - Eel  - E1  - E2
        Pxm =                  - pxel - px1 - px2 
        Pym =                  - pyel - py1 - py2
        Pzm =  Ebeam           - pzel - pz1 - pz2
	Pm  = sqrt(Pxm**2+Pym**2+Pzm**2)
        M2m = Em**2 - Pm**2
	call angles(0,Pxm,Pym,Pzm,thm,phm)
	return
	end

c************************************************************************
	SUBROUTINE miss_mass_4part(Mtarget,Ebeam
     $                            ,Eel,pxel,pyel,pzel
     $                            ,E1,px1,py1,pz1
     $                            ,E2,px2,py2,pz2
     $                            ,E3,px3,py3,pz3
     $                            ,Em,pxm,pym,pzm,M2m,Pm,THm,PHm)
c*************************************************************************
c+ 4-momentum and mass of missing particle calcutation.
c+ The "missing" particle is defined respect to electron particle 1 and particle 2
c+    INPUT:
c+          Mtarget    = target mass
c+          Ebeam      = Beam energy 
c+          Eel,Pxel...= scattered electron 4-momentum (Energy,px,py,pz)
c+          E1,Px1...  = 4-momentum of particle 1
c+          E2,Px2...  = 4-momentum of particle 2
c+          E3,Px3...  = 4-momentum of particle 3
c+    OUTPUT
c+          M2m         = squared missing mass
c+          Pm          = missing momentum
c+          Thm,PHm     = theta and phi of missing particle
c+          pxm,pym,pzm = missing 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Mtarget,Ebeam
	real Eel,pxel,pyel,pzel,M2m,Pm,THm,PHm,E1,px1,py1,pz1,E2,px2,py2,pz2
     #    ,E3,px3,py3,pz3    
	real Em,Pxm,Pym,Pzm
	Em  =  Ebeam + Mtarget - Eel  - E1  - E2  - E3
        Pxm =                  - pxel - px1 - px2 - px3
        Pym =                  - pyel - py1 - py2 - py3
        Pzm =  Ebeam           - pzel - pz1 - pz2 - pz3
	Pm  = sqrt(Pxm**2+Pym**2+Pzm**2)
        M2m = Em**2 - Pm**2
	call angles(0,Pxm,Pym,Pzm,thm,phm)
	return
	end



c************************************************************************
	SUBROUTINE miss_mass_5part(Mtarget,Ebeam
     $                            ,Eel,pxel,pyel,pzel
     $                            ,E1,px1,py1,pz1
     $                            ,E2,px2,py2,pz2
     $                            ,E3,px3,py3,pz3
     $                            ,E4,px4,py4,pz4
     $                            ,Em,pxm,pym,pzm,M2m,Pm,THm,PHm)
c*************************************************************************
c+ 4-momentum and mass of missing particle calcutation.
c+ The "missing" particle is defined respect to electron particle 1 and particle 2
c+    INPUT:
c+          Mtarget    = target mass
c+          Ebeam      = Beam energy 
c+          Eel,Pxel...= scattered electron 4-momentum (Energy,px,py,pz)
c+          E1,Px1...  = 4-momentum of particle 1
c+          E2,Px2...  = 4-momentum of particle 2
c+          E3,Px3...  = 4-momentum of particle 3
c+          E4,Px4...  = 4-momentum of particle 4
c+    OUTPUT
c+          M2m         = squared missing mass
c+          Pm          = missing momentum
c+          Thm,PHm     = theta and phi of missing particle
c+          pxm,pym,pzm = missing 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Mtarget,Ebeam
	real Eel,pxel,pyel,pzel,M2m,Pm,THm,PHm,E1,px1,py1,pz1,E2,px2,py2,pz2
     #    ,E3,px3,py3,pz3,E4,px4,py4,pz4  
	real Em,Pxm,Pym,Pzm
	Em  =  Ebeam + Mtarget - Eel  - E1  - E2  - E3  - E4
        Pxm =                  - pxel - px1 - px2 - px3 - px4
        Pym =                  - pyel - py1 - py2 - py3 - py4
        Pzm =  Ebeam           - pzel - pz1 - pz2 - pz3 - pz4
	Pm  = sqrt(Pxm**2+Pym**2+Pzm**2)
        M2m = Em**2 - Pm**2
	call angles(0,Pxm,Pym,Pzm,thm,phm)
	return
	end


c************************************************************************
	SUBROUTINE miss_mass_6part(Mtarget,Ebeam
     $                            ,Eel,pxel,pyel,pzel
     $                            ,E1,px1,py1,pz1
     $                            ,E2,px2,py2,pz2
     $                            ,E3,px3,py3,pz3
     $                            ,E4,px4,py4,pz4
     $                            ,E5,px5,py5,pz5
     $                            ,Em,pxm,pym,pzm,M2m,Pm,THm,PHm)
c*************************************************************************
c+ 4-momentum and mass of missing particle calcutation.
c+ The "missing" particle is defined respect to electron particle 1 and particle 2
c+    INPUT:
c+          Mtarget    = target mass
c+          Ebeam      = Beam energy 
c+          Eel,Pxel...= scattered electron 4-momentum (Energy,px,py,pz)
c+          E1,Px1...  = 4-momentum of particle 1
c+          E2,Px2...  = 4-momentum of particle 2
c+          E3,Px3...  = 4-momentum of particle 3
c+          E4,Px4...  = 4-momentum of particle 4
c+          E5,Px5...  = 4-momentum of particle 5
c+    OUTPUT
c+          M2m         = squared missing mass
c+          Pm          = missing momentum
c+          Thm,PHm     = theta and phi of missing particle
c+          pxm,pym,pzm = missing 3-momentum
c*************************************************************************
        IMPLICIT NONE
	real Mtarget,Ebeam
	real Eel,pxel,pyel,pzel,M2m,Pm,THm,PHm,E1,px1,py1,pz1,E2,px2,py2,pz2
     #    ,E3,px3,py3,pz3, E4,px4,py4,pz4, E5,px5,py5,pz5  
	real Em,Pxm,Pym,Pzm
	Em  =  Ebeam + Mtarget - Eel  - E1  - E2  - E3  - E4  - E5
        Pxm =                  - pxel - px1 - px2 - px3 - px4 - px5
        Pym =                  - pyel - py1 - py2 - py3 - py4 - py5
        Pzm =  Ebeam           - pzel - pz1 - pz2 - pz3 - pz4 - pz5
	Pm  = sqrt(Pxm**2+Pym**2+Pzm**2)
        M2m = Em**2 - Pm**2
	call angles(0,Pxm,Pym,Pzm,thm,phm)
	return
	end


