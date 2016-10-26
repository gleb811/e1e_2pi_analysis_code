c
c Momentum correction routines for eg1, for e,(p,pi+,pi-) 
c ID   -- particle ID 
c         11   - electron 
c         2212 - proton
c         211  - pi+
c        -211  - pi-
c pmom -- measured momentum on input / corrected momentum on output 
c theta,phi -- particle theta(0-180deg) and phi(0-360deg) angles 
c Ebeam,Torus -- beam energy and torus current
c
      SUBROUTINE cormom_eg1(ID,pmom,theta,phi,Ebeam,Torus) 
      IMPLICIT none
      INTEGER ID
      REAL pmom,theta,phi,Ebeam,Torus
      REAL correct_pel_eg1_2565_1500
      REAL correct_pel_eg1_2565_2250
      REAL correct_pel_eg1_4278_m2250
    
c --- corrections for electrons ---
      IF(ID.eq.11) THEN
      if    (abs(Ebeam-2.5).lt.0.2 .and. abs(Torus-1500.).lt.200. )then
        pmom=correct_pel_eg1_2565_1500(pmom,theta,phi)
      elseif(abs(Ebeam-2.5).lt.0.2 .and. abs(Torus-2250.).lt.200. )then
        pmom=correct_pel_eg1_2565_2250(pmom,theta,phi)
      elseif(abs(Ebeam-4.2).lt.0.2 .and. abs(Torus+2250.).lt.200. )then
        pmom=correct_pel_eg1_4278_m2250(pmom,theta,phi)
      else
       print *,'cormom_eg1: Cannot correct for: EbeamTorus=',Ebeam,Torus
      endif
      RETURN

c --- corrections for protons ---
      ELSEIF(ID.eq.2212) THEN
      pmom = pmom
      RETURN

c --- corrections for pi+ ---
      ELSEIF(ID.eq.211) THEN
      pmom = pmom
      RETURN

c --- corrections for pi- ---
      ELSEIF(ID.eq.-211) THEN
      pmom = pmom
      RETURN

      ELSE
        print *,' cormom_eg1: Cannot correct for ID=',ID
        stop
      ENDIF
      RETURN
      END 




c----------------------------------------------------------------------
c         ELECTRON momentum corrections for EG1: E=2.565, Tor=1500
c----------------------------------------------------------------------
      FUNCTION correct_pel_eg1_2565_1500(pe,theta,phi) 
      REAL pe,theta,phi, phisec,fact
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      IF(phi.lt.  0..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99119E+00+(-0.54721E-03*x)+( 0.30698E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99119E+00+(-0.54721E-03*(x2-d))+( 0.30698E-04*(x2-d)**2)
       f2= 0.99041E+00+(-0.69220E-03*(x2+d))+( 0.22101E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99041E+00+(-0.69220E-03*x)+( 0.22101E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94624E+00+( 0.59791E-02*y)+(-0.17179E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94624E+00+( 0.59791E-02*(y2-d))+(-0.17179E-03*(y2-d)**2)
       f2= 0.10058E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10058E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98642E+00+( 0.27660E-03*x)+( 0.76416E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98642E+00+( 0.27660E-03*(x2-d))+( 0.76416E-04*(x2-d)**2)
       f2= 0.98518E+00+(-0.50972E-03*(x2+d))+( 0.38228E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98518E+00+(-0.50972E-03*x)+( 0.38228E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.86038E+00+( 0.17801E-01*y)+(-0.59855E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.86038E+00+( 0.17801E-01*(y2-d))+(-0.59855E-03*(y2-d)**2)
       f2= 0.99309E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99309E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99847E+00+( 0.50254E-03*x)+( 0.37138E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99847E+00+( 0.50254E-03*(x2-d))+( 0.37138E-04*(x2-d)**2)
       f2= 0.99853E+00+(-0.68549E-06*(x2+d))+( 0.40928E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99853E+00+(-0.68549E-06*x)+( 0.40928E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95459E+00+( 0.52528E-02*y)+(-0.14300E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95459E+00+( 0.52528E-02*(y2-d))+(-0.14300E-03*(y2-d)**2)
       f2= 0.99780E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99780E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99719E+00+( 0.14544E-02*x)+( 0.91894E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99719E+00+( 0.14544E-02*(x2-d))+( 0.91894E-04*(x2-d)**2)
       f2= 0.99477E+00+( 0.16881E-03*(x2+d))+( 0.15959E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99477E+00+( 0.16881E-03*x)+( 0.15959E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90398E+00+( 0.11581E-01*y)+(-0.35022E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90398E+00+( 0.11581E-01*(y2-d))+(-0.35022E-03*(y2-d)**2)
       f2= 0.99885E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99885E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99038E+00+(-0.75910E-04*x)+( 0.44001E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99038E+00+(-0.75910E-04*(x2-d))+( 0.44001E-04*(x2-d)**2)
       f2= 0.99204E+00+( 0.77386E-03*(x2+d))+(-0.31750E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99204E+00+( 0.77386E-03*x)+(-0.31750E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87785E+00+( 0.15148E-01*y)+(-0.46536E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87785E+00+( 0.15148E-01*(y2-d))+(-0.46536E-03*(y2-d)**2)
       f2= 0.99781E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99781E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99720E+00+( 0.49895E-03*x)+( 0.41449E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99720E+00+( 0.49895E-03*(x2-d))+( 0.41449E-04*(x2-d)**2)
       f2= 0.99547E+00+(-0.10284E-02*(x2+d))+( 0.17209E-03*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99547E+00+(-0.10284E-02*x)+( 0.17209E-03*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97400E+00+( 0.27244E-02*y)+(-0.74209E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97400E+00+( 0.27244E-02*(y2-d))+(-0.74209E-04*(y2-d)**2)
       f2= 0.99832E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99832E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_eg1_2565_1500=fact*pe
      RETURN
      END 


C---------------- OLD ----------------
      SUBROUTINE momcorr_eg1_2565_1500_OLD(ch1,ebeam,torcur,sector,pmom,theta,phi,q)
      IMPLICIT none
      CHARACTER*1 ch1
      INTEGER     sector
      REAL        ebeam,torcur,pmom,theta,phi,q
      REAL        ppp, correct_pel_eg1_2565_1500_OLD
      ppp=correct_pel_eg1_2565_1500_OLD(pmom,theta,phi) 
      pmom=ppp
      RETURN
      END
      FUNCTION correct_pel_eg1_2565_1500_OLD(pe,theta,phi) 
      REAL pe,theta,phi, phisec,fact
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      IF(phi.lt.  0..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10001E+01+(-0.76907E-04*x)+( 0.39462E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10001E+01+(-0.76907E-04*(x2-d))+( 0.39462E-05*(x2-d)**2)
       f2= 0.99999E+00+(-0.46639E-03*(x2+d))+( 0.30662E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99999E+00+(-0.46639E-03*x)+( 0.30662E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98406E+00+( 0.20401E-02*y)+(-0.60224E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98406E+00+( 0.20401E-02*(y2-d))+(-0.60224E-04*(y2-d)**2)
       f2= 0.99965E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99965E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99565E+00+(-0.71815E-05*x)+( 0.90285E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99565E+00+(-0.71815E-05*(x2-d))+( 0.90285E-05*(x2-d)**2)
       f2= 0.99538E+00+(-0.15996E-03*(x2+d))+( 0.17694E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99538E+00+(-0.15996E-03*x)+( 0.17694E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97602E+00+( 0.27537E-02*y)+(-0.91049E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97602E+00+( 0.27537E-02*(y2-d))+(-0.91049E-04*(y2-d)**2)
       f2= 0.99701E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99701E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10042E+01+(-0.31014E-03*x)+(-0.16097E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10042E+01+(-0.31014E-03*(x2-d))+(-0.16097E-04*(x2-d)**2)
       f2= 0.10051E+01+( 0.68582E-04*(x2+d))+( 0.56696E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10051E+01+( 0.68582E-04*x)+( 0.56696E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99643E+00+( 0.11154E-02*y)+(-0.33129E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99643E+00+( 0.11154E-02*(y2-d))+(-0.33129E-04*(y2-d)**2)
       f2= 0.10069E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10069E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10038E+01+( 0.18597E-04*x)+( 0.33168E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10038E+01+( 0.18597E-04*(x2-d))+( 0.33168E-05*(x2-d)**2)
       f2= 0.10036E+01+( 0.70895E-04*(x2+d))+(-0.73764E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10036E+01+( 0.70895E-04*x)+(-0.73764E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99309E+00+( 0.12611E-02*y)+(-0.34016E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99309E+00+( 0.12611E-02*(y2-d))+(-0.34016E-04*(y2-d)**2)
       f2= 0.10032E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10032E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10015E+01+(-0.21139E-04*x)+( 0.95377E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10015E+01+(-0.21139E-04*(x2-d))+( 0.95377E-05*(x2-d)**2)
       f2= 0.10018E+01+(-0.14763E-04*(x2+d))+( 0.11893E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10018E+01+(-0.14763E-04*x)+( 0.11893E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98549E+00+( 0.22347E-02*y)+(-0.72396E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98549E+00+( 0.22347E-02*(y2-d))+(-0.72396E-04*(y2-d)**2)
       f2= 0.10028E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10028E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10022E+01+(-0.19988E-03*x)+(-0.68159E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10022E+01+(-0.19988E-03*(x2-d))+(-0.68159E-05*(x2-d)**2)
       f2= 0.10024E+01+( 0.19204E-03*(x2+d))+(-0.11151E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10024E+01+( 0.19204E-03*x)+(-0.11151E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10077E+01+(-0.55799E-03*y)+( 0.15777E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10077E+01+(-0.55799E-03*(y2-d))+( 0.15777E-04*(y2-d)**2)
       f2= 0.10026E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10026E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_eg1_2565_1500_OLD=fact*pe
      RETURN
      END 




c----------------------------------------------------------------------
c         ELECTRON momentum corrections for EG1: E=2.565, Tor=2250
c----------------------------------------------------------------------
      SUBROUTINE momcorr_eg1_2565_2250(ch1,ebeam,torcur,sector,pmom,theta,phi,q)
      IMPLICIT none
      CHARACTER*1 ch1
      INTEGER     sector
      REAL        ebeam,torcur,pmom,theta,phi,q
      REAL        ppp, correct_pel_eg1_2565_2250
      ppp=correct_pel_eg1_2565_2250(pmom,theta,phi) 
      pmom=ppp
      RETURN
      END
      FUNCTION correct_pel_eg1_2565_2250(pe,theta,phi) 
      REAL pe,theta,phi, phisec,fact
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      IF(phi.lt.  0..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99329E+00+(-0.75085E-03*x)+( 0.53410E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99329E+00+(-0.75085E-03*(x2-d))+( 0.53410E-05*(x2-d)**2)
       f2= 0.99363E+00+(-0.24261E-03*(x2+d))+(-0.20681E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99363E+00+(-0.24261E-03*x)+(-0.20681E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94035E+00+( 0.71929E-02*y)+(-0.22219E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94035E+00+( 0.71929E-02*(y2-d))+(-0.22219E-03*(y2-d)**2)
       f2= 0.99493E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99493E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98911E+00+(-0.76298E-03*x)+( 0.17465E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98911E+00+(-0.76298E-03*(x2-d))+( 0.17465E-05*(x2-d)**2)
       f2= 0.98852E+00+(-0.41982E-03*(x2+d))+(-0.10590E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98852E+00+(-0.41982E-03*x)+(-0.10590E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93179E+00+( 0.74324E-02*y)+(-0.22344E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93179E+00+( 0.74324E-02*(y2-d))+(-0.22344E-03*(y2-d)**2)
       f2= 0.99501E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99501E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99775E+00+(-0.32500E-03*x)+(-0.13499E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99775E+00+(-0.32500E-03*(x2-d))+(-0.13499E-04*(x2-d)**2)
       f2= 0.99885E+00+(-0.18775E-03*(x2+d))+( 0.42259E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99885E+00+(-0.18775E-03*x)+( 0.42259E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99070E+00+( 0.43373E-03*y)+( 0.70680E-05*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99070E+00+( 0.43373E-03*(y2-d))+( 0.70680E-05*(y2-d)**2)
       f2= 0.10036E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10036E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99837E+00+( 0.87080E-03*x)+( 0.66595E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99837E+00+( 0.87080E-03*(x2-d))+( 0.66595E-04*(x2-d)**2)
       f2= 0.99617E+00+(-0.15055E-03*(x2+d))+( 0.23526E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99617E+00+(-0.15055E-03*x)+( 0.23526E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94059E+00+( 0.65654E-02*y)+(-0.18523E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94059E+00+( 0.65654E-02*(y2-d))+(-0.18523E-03*(y2-d)**2)
       f2= 0.99981E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99981E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99802E+00+( 0.43280E-04*x)+( 0.22276E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99802E+00+( 0.43280E-04*(x2-d))+( 0.22276E-04*(x2-d)**2)
       f2= 0.99725E+00+(-0.38210E-03*(x2+d))+( 0.33756E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99725E+00+(-0.38210E-03*x)+( 0.33756E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89859E+00+( 0.11633E-01*y)+(-0.32975E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89859E+00+( 0.11633E-01*(y2-d))+(-0.32975E-03*(y2-d)**2)
       f2= 0.10011E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10011E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99720E+00+( 0.11060E-03*x)+( 0.27089E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99720E+00+( 0.11060E-03*(x2-d))+( 0.27089E-04*(x2-d)**2)
       f2= 0.99744E+00+(-0.16652E-04*(x2+d))+( 0.13121E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99744E+00+(-0.16652E-04*x)+( 0.13121E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95136E+00+( 0.56091E-02*y)+(-0.16215E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95136E+00+( 0.56091E-02*(y2-d))+(-0.16215E-03*(y2-d)**2)
       f2= 0.10008E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10008E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_eg1_2565_2250=fact*pe
      RETURN
      END 



c----------------------------------------------------------------------
c         ELECTRON momentum corrections for EG1: E=4.278, Tor=-2250
c----------------------------------------------------------------------
      SUBROUTINE momcorr_eg1_4278_m2250(ch1,ebeam,torcur,sector,pmom,theta,phi,q)
      IMPLICIT none
      CHARACTER*1 ch1
      INTEGER     sector
      REAL        ebeam,torcur,pmom,theta,phi,q
      REAL        ppp, correct_pel_eg1_4278_m2250
      ppp=correct_pel_eg1_4278_m2250(pmom,theta,phi) 
      pmom=ppp
      RETURN
      END
      FUNCTION correct_pel_eg1_4278_m2250(pe,theta,phi) 
      REAL pe,theta,phi, phisec,fact
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      IF(phi.lt.  0..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99196E+00+(-0.57348E-04*x)+( 0.69586E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99196E+00+(-0.57348E-04*(x2-d))+( 0.69586E-04*(x2-d)**2)
       f2= 0.99109E+00+(-0.58745E-03*(x2+d))+( 0.54329E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99109E+00+(-0.58745E-03*x)+( 0.54329E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90213E+00+( 0.14537E-01*y)+(-0.51633E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90213E+00+( 0.14537E-01*(y2-d))+(-0.51633E-03*(y2-d)**2)
       f2= 0.10047E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10047E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99506E+00+( 0.58438E-03*x)+( 0.11770E-03*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99506E+00+( 0.58438E-03*(x2-d))+( 0.11770E-03*(x2-d)**2)
       f2= 0.99337E+00+(-0.30116E-03*(x2+d))+( 0.24332E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99337E+00+(-0.30116E-03*x)+( 0.24332E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91800E+00+( 0.11708E-01*y)+(-0.38019E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91800E+00+( 0.11708E-01*(y2-d))+(-0.38019E-03*(y2-d)**2)
       f2= 0.10111E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10111E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99156E+00+(-0.88216E-03*x)+( 0.10908E-03*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99156E+00+(-0.88216E-03*(x2-d))+( 0.10908E-03*(x2-d)**2)
       f2= 0.99020E+00+(-0.18587E-02*(x2+d))+( 0.53004E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99020E+00+(-0.18587E-02*x)+( 0.53004E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91272E+00+( 0.97536E-02*y)+(-0.23681E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91272E+00+( 0.97536E-02*(y2-d))+(-0.23681E-03*(y2-d)**2)
       f2= 0.10042E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10042E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99388E+00+(-0.18929E-02*x)+( 0.12524E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99388E+00+(-0.18929E-02*(x2-d))+( 0.12524E-04*(x2-d)**2)
       f2= 0.99582E+00+(-0.11192E-02*(x2+d))+( 0.23241E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99582E+00+(-0.11192E-02*x)+( 0.23241E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98314E+00+( 0.15021E-02*y)+(-0.30045E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98314E+00+( 0.15021E-02*(y2-d))+(-0.30045E-04*(y2-d)**2)
       f2= 0.10049E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10049E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99209E+00+(-0.77886E-03*x)+( 0.53973E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99209E+00+(-0.77886E-03*(x2-d))+( 0.53973E-04*(x2-d)**2)
       f2= 0.99152E+00+(-0.96963E-03*(x2+d))+( 0.58517E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99152E+00+(-0.96963E-03*x)+( 0.58517E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90665E+00+( 0.13830E-01*y)+(-0.49760E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90665E+00+( 0.13830E-01*(y2-d))+(-0.49760E-03*(y2-d)**2)
       f2= 0.10073E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10073E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 17.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99244E+00+(-0.65043E-03*x)+( 0.82857E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99244E+00+(-0.65043E-03*(x2-d))+( 0.82857E-04*(x2-d)**2)
       f2= 0.99243E+00+(-0.71359E-03*(x2+d))+( 0.19391E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99243E+00+(-0.71359E-03*x)+( 0.19391E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93370E+00+( 0.84975E-02*y)+(-0.23699E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93370E+00+( 0.84975E-02*(y2-d))+(-0.23699E-03*(y2-d)**2)
       f2= 0.99862E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99862E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_eg1_4278_m2250=fact*pe
      RETURN
      END 








c----------------------------------------------------------------------
c BURKERT's
c ELECTRON momentum corrections for EG1
c----------------------------------------------------------------------
      SUBROUTINE momcorr_eg1(ch1,ebeam,torcur,sector,pmom,theta,phi,q)
      IMPLICIT none
      CHARACTER*1 ch1
      INTEGER     sector,secte
      REAL        pi,ebeam,torcur,pmom,theta,phi,q
      REAL        pel,phiel,phield,thetael,thetaeld
      REAL        cx,cy,cz,pi180


      pi180=acos(-1.)/180.
        secte=sector
        pel=pmom
        cx=sin(theta*pi180)*cos(phi*pi180)
        cy=sin(theta*pi180)*sin(phi*pi180)
        cz=cos(theta*pi180)

        thetael=acos(cz)
        thetaeld=thetael*90./acos(0.)
        phiel=atan2(cy,cx)
        phield=atan2(cy,cx)*90./acos(0.)


       IF( ch1 .eq. 'E' ) THEN
*******************************************************************
* momentum  corrections electron start here
********************************************************************         
*
*eg1 run 2.565GeV
*
********************************************************************         
       if(ebeam.ge.2.5.and.ebeam.lt.2.7)then

      
* -50% field run 
*********************************************************************
* MA corrections
*********************************************************************  
         if(torcur.lt.-1450.and.torcur.gt.-2350.)then
     
         if(secte.eq.1.AND.thetaeld.gt.6)then 
            pel = 0.995*pel/(0.979+0.06/(thetaeld-5.5))
         endif
             
         if(secte.eq.2.AND.thetaeld.gt.6)then
            pel= pel*(1.0391-0.0017*phield+0.000015*(phield**2))*0.98
            pel= pel/(0.957+0.147/(thetaeld-4.26))
         endif
        
         if(secte.eq.3.AND.thetaeld.gt.6)then
            pel= pel*(1.5746-0.0093*phield+0.000035*(phield**2))*0.991
            pel= pel/(0.979+0.06/(thetaeld-5.5))
         endif
        
         if(secte.eq.4.AND.thetaeld.gt.6)then
         if(phield.gt.0) phield=phield
         if(phield.lt.0) phield=phield+360
            pel= pel*(1.6925-0.00711*phield+0.0000178*(phield**2))*0.9915
            pel= pel/(0.979+0.06/(thetaeld-5.5))
         endif
        
         if(secte.eq.5.AND.thetaeld.gt.6)then
            pel= pel*(2.316-0.0109*(phield+360)+0.0000225*
     &		((phield+360)**2))*0.9875
            pel= pel/(0.979+0.06/(thetaeld-5.5))
         endif
        
         if(secte.eq.6.AND.thetaeld.gt.6)then
            pel= pel*(4.663-0.0243*(phield+360)+0.00004*
     &          ((phield+360)**2))*0.991
            pel= pel/(0.979+0.06/(thetaeld-5.5))
         endif

         endif

*********************************************************************
* end of MA corrections
*********************************************************************        
           
*+60% B-field run
         if(torcur.gt.2200.and.torcur.lt.2300.)then

         if(secte.eq.1)then 
         if(thetaeld.lt.25.)pel=pel*(1. - phiel*0.04 + phiel**2*0.04)
         pel = pel*(1. - 0.014)
         endif
           
         if(secte.eq.2)then
         if(thetaeld.lt.25.)pel = pel*(1. -(phiel-pi/3.)*0.03 + 
     *           (phiel-pi/3.)**2*0.05)
         pel = pel*(1. - 0.013)
         endif
             
         if(secte.eq.3)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel-2.*pi/3.)*0.0+
     *     (phiel-2.*pi/3.)**2*0.04)  
         pel = pel*(1. + 0.008)
         endif

         if(secte.eq.4)then
         if(phiel.lt.pi.and.phiel.gt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel-pi)*0.0+
     *      (phiel-pi)**2*0.05)
         endif
         if(phiel.gt.-pi.and.phiel.lt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+pi)*0.015+
     *      (phiel+pi)**2*0.05)
         endif  
         pel = pel*(1. - 0.0095)
         endif

         if(secte.eq.5)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+2.*pi/3.)*0.015+
     *      (phiel+2.*pi/3.)**2*0.05)
         pel = pel*(1. - 0.008)
         endif

         if(secte.eq.6)then
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel+pi/3.)*0.015+
     *      (phiel+pi/3.)**2*0.04)
         pel = pel*(1. + 0.0035)
         endif  
         
         
         endif
            
       endif
       
********************************************************************         
*
*eg1 run 4.235GeV
*
********************************************************************         
       if(ebeam.ge.4.0.and.ebeam.lt.4.5)then
           
*+50% B-field run
         if(torcur.gt.1900.and.torcur.lt.2000.)then
         
         if(secte.eq.1)then 
         if(thetaeld.lt.25.)pel=pel*(1. - phiel*0.06 + phiel**2*0.025)
         pel = pel*(1. - 0.019)
         endif
           
         if(secte.eq.2)then
         if(thetaeld.lt.25.)pel = pel*(1. -(phiel-pi/3.)*0.05 + 
     *      (phiel-pi/3.)**2*0.05)
         if(ebeam.lt.4.1)pel = pel*(1. - 0.018)
         if(ebeam.gt.4.1)pel = pel*(1. - 0.0195)
         endif
             
         if(secte.eq.3)then
         if(thetaeld.lt.25.)pel = pel*(1.+(phiel-2.*pi/3.)*0.02+
     *      (phiel-2.*pi/3.)**2*0.04)  
         if(ebeam.gt.4.1)pel = pel*(1. + 0.020)
         if(ebeam.lt.4.1)pel = pel*(1. + 0.016)
         endif

         if(secte.eq.4)then
         if(phiel.lt.pi.and.phiel.gt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1. + (phiel-pi)*0.025+
     *      (phiel-pi)**2*0.05)
         endif
         if(phiel.gt.-pi.and.phiel.lt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1.+(phiel+pi)*0.025+
     *      (phiel+pi)**2*0.05)
         endif
         if(ebeam.gt.4.1)pel = pel*(1. - 0.0115)
         if(ebeam.lt.4.1)pel = pel*(1. - 0.0130)     
         endif

         if(secte.eq.5)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+2.*pi/3.)*0.015+
     *      (phiel+2.*pi/3.)**2*0.09)
         pel = pel*(1. - 0.016)
         endif

         if(secte.eq.6)then
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel+pi/3.)*0.005+
     *      (phiel+pi/3.)**2*0.04)
         if(ebeam.gt.4.1)pel = pel*(1. + 0.0145)
         if(ebeam.lt.4.1)pel = pel*(1. + 0.0090)
         endif 
          
         endif

           
*+60% B-field run
         if(torcur.gt.2200.and.torcur.lt.2300.)then
         
         if(secte.eq.1)then 
         if(thetaeld.lt.18.)pel = pel*(1 -(18.-thetaeld)**2.5*0.00009) 
         if(thetaeld.lt.25.)pel=pel*(1. - phiel*0.06 + phiel**2*0.035)
         pel = pel*(1. - 0.018)
         endif
           
         if(secte.eq.2)then
         if(thetaeld.lt.18.)pel =pel*(1-(18.-thetaeld)**2.5*0.00009) 
         if(thetaeld.lt.25.)pel =pel*(1. -(phiel-pi/3.)*0.04 + 
     *      (phiel-pi/3.)**2*0.05)
         pel = pel*(1. - 0.0171)
         endif
             
         if(secte.eq.3)then
         if(thetaeld.lt.18.)pel =pel*(1-(18.-thetaeld)**2.5*0.00010) 
         if(thetaeld.lt.25.)pel =pel*(1.+(phiel-2.*pi/3.)*0.015+
     *      (phiel-2.*pi/3.)**2*0.05)  
         pel = pel*(1. + 0.0105)
         endif

         if(secte.eq.4)then
         if(thetaeld.lt.18.)pel = pel*(1 -(18.-thetaeld)**2.5*0.00012) 
         if(phiel.lt.pi.and.phiel.gt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1. + (phiel-pi)*0.025+
     *      (phiel-pi)**2*0.05)
         endif
         if(phiel.gt.-pi.and.phiel.lt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1.+(phiel+pi)*(-0.02)+
     *      (phiel+pi)**2*0.05)
         endif   
         pel = pel*(1. - 0.0113)
         endif

         if(secte.eq.5)then
         if(thetaeld.lt.18.)pel =pel*(1 - (18.-thetaeld)**2.5*0.00012) 
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+2.*pi/3.)*0.015+
     *     (phiel+2.*pi/3.)**2*0.07)
         pel = pel*(1. - 0.0143)
         endif

         if(secte.eq.6)then
         if(thetaeld.lt.18.)pel=pel*(1 - (18.-thetaeld)**2.5*0.00010)            
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel+pi/3.)*0.005+
     *      (phiel+pi/3.)**2*0.04)
         pel = pel*(1. + 0.005)
         endif   
           
         endif

** -60% field run  
         if(torcur.lt.-2100.and.torcur.gt.-2350.)then
     
         if(secte.eq.1)then 
         if(thetaeld.lt.25.)pel=pel*(1. + phiel*0.013 + phiel**2*0.05)
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**2.*0.00030) 
         pel = pel*(1. + 0.016)
         endif
           
         if(secte.eq.2)then
         if(thetaeld.lt.25.)pel = pel*(1. +(phiel-pi/3.)*0.01 + 
     *      (phiel-pi/3.)**2*0.05)
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**1.3*0.0013)
         pel = pel*(1. + 0.006)
         endif
             
         if(secte.eq.3)then
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**1.5*0.0009)
         pel = pel*(1. - 0.0188)   
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel-2.*pi/3.)*0.075+
     *      (phiel-2.*pi/3.)**2*0.06)
         endif

         if(secte.eq.4)then           
         if(phiel.lt.pi.and.phiel.gt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel-pi)*0.06+
     *      (phiel-pi)**2*0.02)
         endif
         if(phiel.gt.-pi.and.phiel.lt.0.)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+pi)*0.020+
     *      (phiel+pi)**2*0.02)
         endif
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**2*0.0002)
         pel = pel*(1. + 0.00)
         endif

         if(secte.eq.5)then
         if(thetaeld.lt.25.)pel = pel*(1.-(phiel+2.*pi/3.)*0.0175+(phiel+2.*pi/3.)**2*0.06)
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**2.*0.0002)
         pel = pel*(1.+ 0.0035)
         endif

         if(secte.eq.6)then
         if(thetaeld.lt.25.)pel = pel*(1. - (phiel+pi/3.)*0.035+
     *      (phiel+pi/3.)**2*0.06)
         if(thetaeld.lt.18.)pel = pel*(1 - (18.-thetaeld)**2*0.0001)
         pel = pel*(1. - 0.02)
         endif  
         
         endif

       endif 

       pmom=pel

******************************************************************
*  corrections end here        
******************************************************************     
       ELSEIF(ch1 .eq. 'H') THEN   
              continue
**********************************************************************
*  no momentum correction hadron
**********************************************************************                    
       ELSE
              print *,' ch1 is bad =',ch1
              stop
       ENDIF


      RETURN
      END

