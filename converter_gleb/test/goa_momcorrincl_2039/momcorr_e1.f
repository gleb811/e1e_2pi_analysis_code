c
c Momentum correction routines for e1, for e,p,pi+,pi- 
c ID   -- particle ID 
c         11   - electron 
c         2212 - proton
c         211  - pi+
c        -211  - pi-
c pmom -> measured momentum on input / corrected momentum on output 
c theta,phi -> particle theta(0-180deg) and phi(0-360deg) angles 
c Ebeam,Torus -> beam energy and torus current
c
c E1a and E1B runs are distinguished with Ebeam,Torus
c they must be set exactly (4 digits)
c
      SUBROUTINE cormom_e1(ID,pmom,theta,phi,Ebeam,Torus) 
      IMPLICIT none
      INTEGER ID
      REAL pmom,theta,phi,Ebeam,Torus
      REAL correct_ppr_e1_feb98
      REAL correct_ppip_e1_feb98
      REAL correct_ppim_e1_feb98
      REAL correct_pel_e1_1645_1500_feb98,correct_ppr_e1_1645_1500_feb98
      REAL correct_pel_e1_1645_2250_feb98,correct_ppr_e1_1645_2250_feb98
      REAL correct_pel_e1_2445_1500_feb98,correct_ppr_e1_2445_1500_feb98
      REAL correct_pel_e1_2445_2250_feb98,correct_ppr_e1_2445_2250_feb98
      REAL correct_pel_e1_4045_2250_feb98,correct_ppr_e1_4045_2250_feb98
      REAL correct_pel_e1_4045_3375_feb98,correct_ppr_e1_4045_3375_feb98
      REAL correct_pel_e1_1515_0750_jan99,correct_ppr_e1_1515_0750_jan99
      REAL correct_pel_e1_1515_1500_jan99,correct_ppr_e1_1515_1500_jan99
      REAL correct_pel_e1_1515_2250_jan99,correct_ppr_e1_1515_2250_jan99
      REAL correct_pel_e1_2567_1500_jan99,correct_ppr_e1_2567_1500_jan99
      REAL correct_pel_e1_2567_2250_jan99,correct_ppr_e1_2567_2250_jan99
      REAL correct_pel_e1_4056_2250_jan99,correct_ppr_e1_4056_2250_jan99
      REAL correct_pel_e1_4247_2250_jan99,correct_ppr_e1_4247_2250_jan99
      REAL correct_pel_e1_4462_2250_jan99,correct_ppr_e1_4462_2250_jan99
      REAL correct_pel_e1_4462_3375_jan99,correct_ppr_e1_4462_3375_jan99
      REAL correct_pel_e1_4817_3375_feb00,correct_ppr_e1_4817_3375_feb00
      REAL correct_pel_e1_2039_2250_feb09,correct_thel_e1_2039_2250_feb09
      
      INTEGER n1
      DATA n1/0/

      n1=n1+1
      IF(ID.ne.11.and.ID.ne.2212.and.ID.ne.211.and.ID.ne.-211)THEN
        print *,' cormom_e1: Cannot correct for ID=',ID
        stop
      ENDIF
    



      IF    (abs(Ebeam-1.645).lt.0.002 .and. abs(Torus-1500.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 1645-1500 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_1645_1500_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_1645_1500_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)

      ELSEIF(abs(Ebeam-1.645).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 1645-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_1645_2250_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_1645_2250_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)

      ELSEIF(abs(Ebeam-2.445).lt.0.002 .and. abs(Torus-1500.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 2445-1500 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_2445_1500_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_2445_1500_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)

      ELSEIF(abs(Ebeam-2.445).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 2445-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_2445_2250_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_2445_2250_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)

      ELSEIF(abs(Ebeam-4.045).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 4045-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4045_2250_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_4045_2250_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)

      ELSEIF(abs(Ebeam-4.045).lt.0.002 .and. abs(Torus-3375.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1A 4045-3375 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4045_3375_feb98(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_feb98(pmom)
      if(ID.eq.2212) pmom=correct_ppr_e1_4045_3375_feb98(pmom,theta,phi)
      if(ID.eq.211)  pmom=correct_ppip_e1_feb98(pmom)
      if(ID.eq.-211) pmom=correct_ppim_e1_feb98(pmom)



      ELSEIF(abs(Ebeam-1.515).lt.0.002 .and. abs(Torus-0750.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 1515-0750 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_1515_0750_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_1515_0750_jan99(pmom,theta,phi)
      

      ELSEIF(abs(Ebeam-1.515).lt.0.002 .and. abs(Torus-1500.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 1515-1500 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_1515_1500_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_1515_1500_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=pmom+0.000425/pmom/pmom/pmom+(0.0095-0.01*
     & pmom)*cos(1.8*(theta-7.)*3.14/180./(10.*pmom-1.5)+0.95)+0.003 !Gleb
c     7.11.03


      ELSEIF(abs(Ebeam-2.039).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 2039-2250 to use'
      if(ID.eq.11)  then 
      theta=correct_thel_e1_2039_2250_feb09(pmom,theta,phi)
      pmom=correct_pel_e1_2039_2250_feb09(pmom,theta,phi)
      endif !Gleb     6.02.09
      
	
      ELSEIF(abs(Ebeam-1.515).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 1515-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_1515_2250_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_1515_2250_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-2.567).lt.0.002 .and. abs(Torus-1500.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 2567-1500 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_2567_1500_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_2567_1500_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-2.567).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 2567-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_2567_2250_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_2567_2250_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-4.056).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 4056-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4056_2250_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_4056_2250_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-4.247).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 4247-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4247_2250_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_4247_2250_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-4.462).lt.0.002 .and. abs(Torus-2250.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 4462-2250 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4462_2250_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_4462_2250_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-4.462).lt.0.002 .and. abs(Torus-3375.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1B 4462-3375 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4462_3375_jan99(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_4462_3375_jan99(pmom,theta,phi)

      ELSEIF(abs(Ebeam-4.817).lt.0.002 .and. abs(Torus-3375.).lt.2.)THEN
      if(n1.eq.1) print *,'INFO: mom.corr. for E1D 4817-3375 to use'
      if(ID.eq.11)   pmom=correct_pel_e1_4817_3375_feb00(pmom,theta,phi)
      if(ID.eq.2212) pmom=correct_ppr_e1_4817_3375_feb00(pmom,theta,phi)


      ELSE
       print *,'cormom_e1: Cannot correct for: EbeamTorus=',Ebeam,Torus
       stop
      ENDIF
      RETURN
      END 




      FUNCTION correct_ppr_e1_feb98(x)
      IMPLICIT none
      REAL correct_ppr_e1_feb98,x,fact      
      fact=  0.99273 + exp(  0.09722+( -8.21994)*x)
      correct_ppr_e1_feb98 = fact*x
      RETURN
      END
      FUNCTION correct_ppip_e1_feb98(x)
      IMPLICIT none
      REAL correct_ppip_e1_feb98,x,fact
      fact=  0.99823
      correct_ppip_e1_feb98 = fact*x
      RETURN
      END
      FUNCTION correct_ppim_e1_feb98(x)
      IMPLICIT none
      REAL correct_ppim_e1_feb98,x,fact
      fact=  0.99849
      correct_ppim_e1_feb98 = fact*x
      RETURN
      END




      FUNCTION correct_pel_e1_1645_1500_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_1645_1500_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 24.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99936E+00+(-0.64973E-03*x)+(-0.19003E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99936E+00+(-0.64973E-03*(x2-d))+(-0.19003E-04*(x2-d)**2)
       f2= 0.99995E+00+(-0.43507E-03*(x2+d))+( 0.34318E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99995E+00+(-0.43507E-03*x)+( 0.34318E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94580E+00+( 0.50879E-02*y)+(-0.11363E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94580E+00+( 0.50879E-02*(y2-d))+(-0.11363E-03*(y2-d)**2)
       f2= 0.10014E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10014E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99760E+00+(-0.42329E-03*x)+(-0.61233E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99760E+00+(-0.42329E-03*(x2-d))+(-0.61233E-05*(x2-d)**2)
       f2= 0.99820E+00+(-0.18622E-03*(x2+d))+( 0.21691E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99820E+00+(-0.18622E-03*x)+( 0.21691E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92701E+00+( 0.67240E-02*y)+(-0.15212E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92701E+00+( 0.67240E-02*(y2-d))+(-0.15212E-03*(y2-d)**2)
       f2= 0.99992E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99992E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 25.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99786E+00+( 0.49067E-04*x)+(-0.11273E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99786E+00+( 0.49067E-04*(x2-d))+(-0.11273E-05*(x2-d)**2)
       f2= 0.99815E+00+( 0.10976E-03*(x2+d))+( 0.14762E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99815E+00+( 0.10976E-03*x)+( 0.14762E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95745E+00+( 0.38197E-02*y)+(-0.86591E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95745E+00+( 0.38197E-02*(y2-d))+(-0.86591E-04*(y2-d)**2)
       f2= 0.99872E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99872E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99794E+00+(-0.29931E-03*x)+(-0.74919E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99794E+00+(-0.29931E-03*(x2-d))+(-0.74919E-05*(x2-d)**2)
       f2= 0.99865E+00+( 0.16602E-03*(x2+d))+( 0.24172E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99865E+00+( 0.16602E-03*x)+( 0.24172E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93540E+00+( 0.61294E-02*y)+(-0.14256E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93540E+00+( 0.61294E-02*(y2-d))+(-0.14256E-03*(y2-d)**2)
       f2= 0.10000E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99809E+00+(-0.28867E-03*x)+( 0.27985E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99809E+00+(-0.28867E-03*(x2-d))+( 0.27985E-05*(x2-d)**2)
       f2= 0.99859E+00+(-0.13738E-03*(x2+d))+( 0.15286E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99859E+00+(-0.13738E-03*x)+( 0.15286E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93006E+00+( 0.65380E-02*y)+(-0.15088E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93006E+00+( 0.65380E-02*(y2-d))+(-0.15088E-03*(y2-d)**2)
       f2= 0.10001E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10001E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99894E+00+(-0.11803E-03*x)+( 0.58102E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99894E+00+(-0.11803E-03*(x2-d))+( 0.58102E-05*(x2-d)**2)
       f2= 0.99970E+00+( 0.19080E-03*(x2+d))+(-0.23286E-06*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99970E+00+( 0.19080E-03*x)+(-0.23286E-06*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94488E+00+( 0.52121E-02*y)+(-0.12020E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94488E+00+( 0.52121E-02*(y2-d))+(-0.12020E-03*(y2-d)**2)
       f2= 0.10018E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10018E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_1645_1500_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_1645_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_1645_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 24.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10006E+01+(-0.29145E-03*x)+(-0.11545E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10006E+01+(-0.29145E-03*(x2-d))+(-0.11545E-04*(x2-d)**2)
       f2= 0.10007E+01+(-0.27430E-03*(x2+d))+( 0.12005E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+(-0.27430E-03*x)+( 0.12005E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94522E+00+( 0.49157E-02*y)+(-0.10739E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94522E+00+( 0.49157E-02*(y2-d))+(-0.10739E-03*(y2-d)**2)
       f2= 0.10007E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99751E+00+(-0.32995E-03*x)+(-0.90206E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99751E+00+(-0.32995E-03*(x2-d))+(-0.90206E-05*(x2-d)**2)
       f2= 0.99796E+00+(-0.43782E-04*(x2+d))+( 0.44459E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99796E+00+(-0.43782E-04*x)+( 0.44459E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90574E+00+( 0.84052E-02*y)+(-0.18942E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90574E+00+( 0.84052E-02*(y2-d))+(-0.18942E-03*(y2-d)**2)
       f2= 0.99863E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99863E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 25.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99784E+00+( 0.11602E-03*x)+( 0.55825E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99784E+00+( 0.11602E-03*(x2-d))+( 0.55825E-05*(x2-d)**2)
       f2= 0.99781E+00+( 0.12361E-03*(x2+d))+( 0.62069E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99781E+00+( 0.12361E-03*x)+( 0.62069E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95906E+00+( 0.33298E-02*y)+(-0.70187E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95906E+00+( 0.33298E-02*(y2-d))+(-0.70187E-04*(y2-d)**2)
       f2= 0.99820E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99820E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99940E+00+(-0.54918E-04*x)+(-0.16128E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99940E+00+(-0.54918E-04*(x2-d))+(-0.16128E-05*(x2-d)**2)
       f2= 0.99950E+00+( 0.57119E-04*(x2+d))+(-0.15313E-07*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99950E+00+( 0.57119E-04*x)+(-0.15313E-07*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92731E+00+( 0.66959E-02*y)+(-0.15321E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92731E+00+( 0.66959E-02*(y2-d))+(-0.15321E-03*(y2-d)**2)
       f2= 0.99962E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99962E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99899E+00+(-0.16902E-03*x)+( 0.31202E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99899E+00+(-0.16902E-03*(x2-d))+( 0.31202E-05*(x2-d)**2)
       f2= 0.99916E+00+(-0.11756E-03*(x2+d))+( 0.69546E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99916E+00+(-0.11756E-03*x)+( 0.69546E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90220E+00+( 0.89001E-02*y)+(-0.20241E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90220E+00+( 0.89001E-02*(y2-d))+(-0.20241E-03*(y2-d)**2)
       f2= 0.99974E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99974E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99950E+00+(-0.11772E-03*x)+(-0.34171E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99950E+00+(-0.11772E-03*(x2-d))+(-0.34171E-05*(x2-d)**2)
       f2= 0.99983E+00+( 0.53398E-04*(x2+d))+(-0.11587E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99983E+00+( 0.53398E-04*x)+(-0.11587E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92406E+00+( 0.72713E-02*y)+(-0.17426E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92406E+00+( 0.72713E-02*(y2-d))+(-0.17426E-03*(y2-d)**2)
       f2= 0.10004E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10004E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_1645_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_2445_1500_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_2445_1500_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.10002E+01+(-0.81667E-03*x)+(-0.32988E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10002E+01+(-0.81667E-03*(x2-d))+(-0.32988E-04*(x2-d)**2)
       f2= 0.10007E+01+(-0.61884E-03*(x2+d))+( 0.30477E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+(-0.61884E-03*x)+( 0.30477E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97217E+00+( 0.27974E-02*y)+(-0.66075E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97217E+00+( 0.27974E-02*(y2-d))+(-0.66075E-04*(y2-d)**2)
       f2= 0.10014E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10014E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99840E+00+(-0.30060E-03*x)+( 0.66853E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99840E+00+(-0.30060E-03*(x2-d))+( 0.66853E-05*(x2-d)**2)
       f2= 0.99878E+00+(-0.27046E-03*(x2+d))+( 0.18533E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99878E+00+(-0.27046E-03*x)+( 0.18533E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97522E+00+( 0.18720E-02*y)+(-0.27655E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97522E+00+( 0.18720E-02*(y2-d))+(-0.27655E-04*(y2-d)**2)
       f2= 0.10000E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99800E+00+( 0.14873E-03*x)+( 0.28774E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99800E+00+( 0.14873E-03*(x2-d))+( 0.28774E-05*(x2-d)**2)
       f2= 0.99811E+00+( 0.21780E-03*(x2+d))+( 0.11960E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99811E+00+( 0.21780E-03*x)+( 0.11960E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.96341E+00+( 0.34411E-02*y)+(-0.82161E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.96341E+00+( 0.34411E-02*(y2-d))+(-0.82161E-04*(y2-d)**2)
       f2= 0.99894E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99894E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99876E+00+(-0.25433E-04*x)+( 0.11394E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99876E+00+(-0.25433E-04*(x2-d))+( 0.11394E-05*(x2-d)**2)
       f2= 0.99905E+00+( 0.56483E-04*(x2+d))+( 0.13479E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99905E+00+( 0.56483E-04*x)+( 0.13479E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94567E+00+( 0.52899E-02*y)+(-0.12615E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94567E+00+( 0.52899E-02*(y2-d))+(-0.12615E-03*(y2-d)**2)
       f2= 0.99996E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99996E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99761E+00+(-0.47177E-03*x)+(-0.14013E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99761E+00+(-0.47177E-03*(x2-d))+(-0.14013E-04*(x2-d)**2)
       f2= 0.99858E+00+(-0.13128E-03*(x2+d))+( 0.10495E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99858E+00+(-0.13128E-03*x)+( 0.10495E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.96045E+00+( 0.34802E-02*y)+(-0.74744E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.96045E+00+( 0.34802E-02*(y2-d))+(-0.74744E-04*(y2-d)**2)
       f2= 0.99967E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99967E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10001E+01+( 0.12008E-03*x)+( 0.12930E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10001E+01+( 0.12008E-03*(x2-d))+( 0.12930E-04*(x2-d)**2)
       f2= 0.99991E+00+( 0.10552E-03*(x2+d))+( 0.93861E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99991E+00+( 0.10552E-03*x)+( 0.93861E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10032E+01+(-0.13094E-02*y)+( 0.63147E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10032E+01+(-0.13094E-02*(y2-d))+( 0.63147E-04*(y2-d)**2)
       f2= 0.10019E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10019E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_2445_1500_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_2445_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_2445_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.10000E+01+(-0.54681E-03*x)+(-0.25103E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10000E+01+(-0.54681E-03*(x2-d))+(-0.25103E-04*(x2-d)**2)
       f2= 0.10002E+01+(-0.42091E-03*(x2+d))+( 0.21770E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+(-0.42091E-03*x)+( 0.21770E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95579E+00+( 0.42112E-02*y)+(-0.97963E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95579E+00+( 0.42112E-02*(y2-d))+(-0.97963E-04*(y2-d)**2)
       f2= 0.10007E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99788E+00+(-0.30749E-03*x)+(-0.18367E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99788E+00+(-0.30749E-03*(x2-d))+(-0.18367E-05*(x2-d)**2)
       f2= 0.99848E+00+(-0.63557E-04*(x2+d))+( 0.55381E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99848E+00+(-0.63557E-04*x)+( 0.55381E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97129E+00+( 0.21313E-02*y)+(-0.35032E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97129E+00+( 0.21313E-02*(y2-d))+(-0.35032E-04*(y2-d)**2)
       f2= 0.99932E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99932E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99745E+00+( 0.30469E-04*x)+(-0.20229E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99745E+00+( 0.30469E-04*(x2-d))+(-0.20229E-05*(x2-d)**2)
       f2= 0.99777E+00+( 0.13921E-03*(x2+d))+( 0.94922E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99777E+00+( 0.13921E-03*x)+( 0.94922E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95458E+00+( 0.41239E-02*y)+(-0.96783E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95458E+00+( 0.41239E-02*(y2-d))+(-0.96783E-04*(y2-d)**2)
       f2= 0.99820E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99820E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99875E+00+(-0.31380E-04*x)+( 0.10903E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99875E+00+(-0.31380E-04*(x2-d))+( 0.10903E-05*(x2-d)**2)
       f2= 0.99893E+00+( 0.43280E-04*(x2+d))+( 0.72802E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99893E+00+( 0.43280E-04*x)+( 0.72802E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94928E+00+( 0.47659E-02*y)+(-0.11138E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94928E+00+( 0.47659E-02*(y2-d))+(-0.11138E-03*(y2-d)**2)
       f2= 0.99949E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99949E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99851E+00+(-0.16716E-03*x)+( 0.48566E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99851E+00+(-0.16716E-03*(x2-d))+( 0.48566E-05*(x2-d)**2)
       f2= 0.99872E+00+(-0.13243E-03*(x2+d))+( 0.11358E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99872E+00+(-0.13243E-03*x)+( 0.11358E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93254E+00+( 0.63525E-02*y)+(-0.14954E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93254E+00+( 0.63525E-02*(y2-d))+(-0.14954E-03*(y2-d)**2)
       f2= 0.99946E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99946E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99941E+00+(-0.63481E-04*x)+( 0.14471E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99941E+00+(-0.63481E-04*(x2-d))+( 0.14471E-05*(x2-d)**2)
       f2= 0.99970E+00+( 0.80771E-04*(x2+d))+( 0.22055E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99970E+00+( 0.80771E-04*x)+( 0.22055E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98771E+00+( 0.65659E-03*y)+(-0.10941E-05*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98771E+00+( 0.65659E-03*(y2-d))+(-0.10941E-05*(y2-d)**2)
       f2= 0.10006E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10006E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_2445_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_4045_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4045_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.99832E+00+(-0.69250E-03*x)+(-0.16079E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99832E+00+(-0.69250E-03*(x2-d))+(-0.16079E-04*(x2-d)**2)
       f2= 0.99856E+00+(-0.85013E-03*(x2+d))+( 0.54743E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99856E+00+(-0.85013E-03*x)+( 0.54743E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95606E+00+( 0.44187E-02*y)+(-0.10596E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95606E+00+( 0.44187E-02*(y2-d))+(-0.10596E-03*(y2-d)**2)
       f2= 0.10017E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10017E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99507E+00+(-0.54582E-03*x)+( 0.14303E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99507E+00+(-0.54582E-03*(x2-d))+( 0.14303E-04*(x2-d)**2)
       f2= 0.99605E+00+(-0.49060E-03*(x2+d))+( 0.52641E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99605E+00+(-0.49060E-03*x)+( 0.52641E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92843E+00+( 0.72016E-02*y)+(-0.17520E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92843E+00+( 0.72016E-02*(y2-d))+(-0.17520E-03*(y2-d)**2)
       f2= 0.10002E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99586E+00+(-0.22353E-03*x)+(-0.97580E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99586E+00+(-0.22353E-03*(x2-d))+(-0.97580E-05*(x2-d)**2)
       f2= 0.99644E+00+(-0.94174E-04*(x2+d))+( 0.35394E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99644E+00+(-0.94174E-04*x)+( 0.35394E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94681E+00+( 0.55336E-02*y)+(-0.14506E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94681E+00+( 0.55336E-02*(y2-d))+(-0.14506E-03*(y2-d)**2)
       f2= 0.99838E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99838E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99620E+00+( 0.24036E-03*x)+( 0.33061E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99620E+00+( 0.24036E-03*(x2-d))+( 0.33061E-04*(x2-d)**2)
       f2= 0.99637E+00+( 0.80363E-05*(x2+d))+( 0.34197E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99637E+00+( 0.80363E-05*x)+( 0.34197E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91854E+00+( 0.85857E-02*y)+(-0.22231E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91854E+00+( 0.85857E-02*(y2-d))+(-0.22231E-03*(y2-d)**2)
       f2= 0.10001E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10001E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99609E+00+( 0.10402E-03*x)+( 0.38572E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99609E+00+( 0.10402E-03*(x2-d))+( 0.38572E-04*(x2-d)**2)
       f2= 0.99617E+00+(-0.27580E-03*(x2+d))+( 0.31518E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99617E+00+(-0.27580E-03*x)+( 0.31518E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93887E+00+( 0.58215E-02*y)+(-0.13552E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93887E+00+( 0.58215E-02*(y2-d))+(-0.13552E-03*(y2-d)**2)
       f2= 0.10003E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10003E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99693E+00+( 0.11162E-03*x)+( 0.30892E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99693E+00+( 0.11162E-03*(x2-d))+( 0.30892E-04*(x2-d)**2)
       f2= 0.99734E+00+( 0.26211E-03*(x2+d))+( 0.86266E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99734E+00+( 0.26211E-03*x)+( 0.86266E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94791E+00+( 0.52384E-02*y)+(-0.12647E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94791E+00+( 0.52384E-02*(y2-d))+(-0.12647E-03*(y2-d)**2)
       f2= 0.10025E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10025E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4045_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_4045_3375_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4045_3375_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.99960E+00+(-0.50306E-03*x)+(-0.23015E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99960E+00+(-0.50306E-03*(x2-d))+(-0.23015E-04*(x2-d)**2)
       f2= 0.99951E+00+(-0.47673E-03*(x2+d))+( 0.25228E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99951E+00+(-0.47673E-03*x)+( 0.25228E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99490E+00+( 0.15784E-03*y)+( 0.53493E-05*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99490E+00+( 0.15784E-03*(y2-d))+( 0.53493E-05*(y2-d)**2)
       f2= 0.10004E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10004E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99633E+00+(-0.42823E-03*x)+( 0.12525E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99633E+00+(-0.42823E-03*(x2-d))+( 0.12525E-05*(x2-d)**2)
       f2= 0.99719E+00+(-0.18927E-03*(x2+d))+( 0.21606E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99719E+00+(-0.18927E-03*x)+( 0.21606E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94967E+00+( 0.47126E-02*y)+(-0.11035E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94967E+00+( 0.47126E-02*(y2-d))+(-0.11035E-03*(y2-d)**2)
       f2= 0.99925E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99925E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99692E+00+(-0.90842E-04*x)+(-0.82156E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99692E+00+(-0.90842E-04*(x2-d))+(-0.82156E-05*(x2-d)**2)
       f2= 0.99728E+00+(-0.32890E-05*(x2+d))+( 0.17760E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99728E+00+(-0.32890E-05*x)+( 0.17760E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.96419E+00+( 0.33557E-02*y)+(-0.82398E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.96419E+00+( 0.33557E-02*(y2-d))+(-0.82398E-04*(y2-d)**2)
       f2= 0.99807E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99807E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 21.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99703E+00+(-0.12782E-03*x)+( 0.21406E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99703E+00+(-0.12782E-03*(x2-d))+( 0.21406E-05*(x2-d)**2)
       f2= 0.99799E+00+( 0.23832E-03*(x2+d))+(-0.10893E-06*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99799E+00+( 0.23832E-03*x)+(-0.10893E-06*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94730E+00+( 0.52290E-02*y)+(-0.12956E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94730E+00+( 0.52290E-02*(y2-d))+(-0.12956E-03*(y2-d)**2)
       f2= 0.99922E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99922E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99736E+00+( 0.11615E-04*x)+( 0.26611E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99736E+00+( 0.11615E-04*(x2-d))+( 0.26611E-04*(x2-d)**2)
       f2= 0.99732E+00+(-0.16707E-03*(x2+d))+( 0.19552E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99732E+00+(-0.16707E-03*x)+( 0.19552E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94597E+00+( 0.51349E-02*y)+(-0.12272E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94597E+00+( 0.51349E-02*(y2-d))+(-0.12272E-03*(y2-d)**2)
       f2= 0.99930E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99930E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99813E+00+(-0.87278E-04*x)+( 0.24255E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99813E+00+(-0.87278E-04*(x2-d))+( 0.24255E-05*(x2-d)**2)
       f2= 0.99871E+00+( 0.17929E-03*(x2+d))+(-0.39571E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99871E+00+( 0.17929E-03*x)+(-0.39571E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95981E+00+( 0.41010E-02*y)+(-0.10509E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95981E+00+( 0.41010E-02*(y2-d))+(-0.10509E-03*(y2-d)**2)
       f2= 0.10003E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10003E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4045_3375_feb98=fact*pe
      RETURN
      END 

      FUNCTION correct_pel_e1_1515_0750_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_1515_0750_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.98789E+00+(-0.10049E-02*x)+( 0.24523E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98789E+00+(-0.10049E-02*(x2-d))+( 0.24523E-04*(x2-d)**2)
       f2= 0.98911E+00+(-0.41700E-03*(x2+d))+( 0.43321E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98911E+00+(-0.41700E-03*x)+( 0.43321E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.83176E+00+( 0.17589E-01*y)+(-0.46334E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.83176E+00+( 0.17589E-01*(y2-d))+(-0.46334E-03*(y2-d)**2)
       f2= 0.10197E+01+(-0.10342E-02*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10197E+01+(-0.10342E-02*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98344E+00+(-0.45686E-03*x)+( 0.43023E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98344E+00+(-0.45686E-03*(x2-d))+( 0.43023E-04*(x2-d)**2)
       f2= 0.98332E+00+(-0.89893E-03*(x2+d))+( 0.70038E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98332E+00+(-0.89893E-03*x)+( 0.70038E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.83388E+00+( 0.16002E-01*y)+(-0.40814E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.83388E+00+( 0.16002E-01*(y2-d))+(-0.40814E-03*(y2-d)**2)
       f2= 0.10007E+01+(-0.50397E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+(-0.50397E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99364E+00+( 0.14556E-03*x)+( 0.11373E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99364E+00+( 0.14556E-03*(x2-d))+( 0.11373E-04*(x2-d)**2)
       f2= 0.99459E+00+( 0.14509E-03*(x2+d))+( 0.45459E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99459E+00+( 0.14509E-03*x)+( 0.45459E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91028E+00+( 0.84242E-02*y)+(-0.19634E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91028E+00+( 0.84242E-02*(y2-d))+(-0.19634E-03*(y2-d)**2)
       f2= 0.99959E+00+( 0.74955E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99959E+00+( 0.74955E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98935E+00+(-0.23651E-05*x)+( 0.36501E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98935E+00+(-0.23651E-05*(x2-d))+( 0.36501E-04*(x2-d)**2)
       f2= 0.98968E+00+(-0.49565E-05*(x2+d))+( 0.50809E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98968E+00+(-0.49565E-05*x)+( 0.50809E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.83445E+00+( 0.16390E-01*y)+(-0.41157E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.83445E+00+( 0.16390E-01*(y2-d))+(-0.41157E-03*(y2-d)**2)
       f2= 0.10007E+01+(-0.13593E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+(-0.13593E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98879E+00+(-0.15100E-02*x)+(-0.25132E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98879E+00+(-0.15100E-02*(x2-d))+(-0.25132E-04*(x2-d)**2)
       f2= 0.99005E+00+(-0.67608E-03*(x2+d))+( 0.77020E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99005E+00+(-0.67608E-03*x)+( 0.77020E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.82114E+00+( 0.17554E-01*y)+(-0.43332E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.82114E+00+( 0.17554E-01*(y2-d))+(-0.43332E-03*(y2-d)**2)
       f2= 0.10006E+01+(-0.76020E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10006E+01+(-0.76020E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99304E+00+( 0.85386E-03*x)+( 0.82201E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99304E+00+( 0.85386E-03*(x2-d))+( 0.82201E-04*(x2-d)**2)
       f2= 0.99281E+00+( 0.56885E-03*(x2+d))+( 0.20139E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99281E+00+( 0.56885E-03*x)+( 0.20139E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.88040E+00+( 0.11128E-01*y)+(-0.26235E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.88040E+00+( 0.11128E-01*(y2-d))+(-0.26235E-03*(y2-d)**2)
       f2= 0.99233E+00+( 0.28047E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99233E+00+( 0.28047E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_1515_0750_jan99=fact*pe
      RETURN
      END 



      FUNCTION correct_thel_e1_2039_2250_feb09(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_thel_e1_2039_2250_feb09
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN

        if ((theta.gt.16.).and.(theta.lt.32.5)) then
	 x=0.2051420E-01-.5001128E-03*(theta)
	 y=-.2507647E-02+0.9725749E-04*(theta)
	 x2=0.1526173E-03+0.4581319E-05*(theta)+0.1142428E-03*SIN(449./theta+8.054586)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	
         elseif ((theta.gt.32.5).and.(theta.lt.44.)) then
	 x=6.252384-.4711988*(theta)+0.1178777E-01*((theta)**2)-.9834436E-04*((theta)**3)
	 y=-.2758633+0.2160061E-01*(theta)-.5579384E-03*((theta)**2)+0.4775051E-05*((theta)**3)
	 x2=-.3277342E-03+0.1174220E-04*(theta)-.2276303E-03*SIN(449./theta-8.186739)
	 fact=x+y*(phisec)+x2*((phisec)**2) 

	
	 
	
	 else
	 fact=0
	 endif   	
	    
      ENDIF
      
      IF(nsector.eq. 2)THEN
      
        if ((theta.gt.16.).and.(theta.lt.32.5)) then
	 x=-.1339374+0.1013706E-01*(theta)-.1723397E-03*((theta)**2)
	 y=-.1638717E-02+0.9508414E-04*(theta)
	 x2=0.2754756E-03+0.7588073E-07*(theta)-.1080132E-03*SIN(449./theta-1.257899)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	
         elseif ((theta.gt.32.5).and.(theta.lt.44.)) then
	 x=5.225882-.4008103*(theta)+0.1023979E-01*((theta)**2)-.8721992E-04*((theta)**3)
	 y=-.6548676+0.5129098E-01*(theta)-.1329902E-02*((theta)**2)+0.1144865E-04*((theta)**3)
	 x2=-.3514004E-01+0.1304112E-02*(theta)+0.2458228E-04*((theta)**2)-.7083146E-06*((theta)**3)
     &-.2314476E-07*((theta)**4)+0.4713921E-09*((theta)**5)
	 fact=x+y*(phisec)+x2*((phisec)**2) 

	
	 
	
	 else
	 fact=0
	 endif   		

      ENDIF
      
      
      IF(nsector.eq. 3)THEN
      

        if ((theta.gt.16.).and.(theta.lt.44.)) then
	 x=0.5706769-.9922977E-01*(theta)+0.5942143E-02*((theta)**2)-.1460487E-03*((theta)**3)+
     &0.1269023E-05*((theta)**4)
	 y=0.7726650E-02-.1776088E-03*(theta)
	 x2=0.5845903E-03-.1155695E-04*(theta)-.1278195E-03*SIN(449./theta-1.534029)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	

	
	 else
	 fact=0
	 endif   	
      
    
      ENDIF
      
      
      
      IF(nsector.eq. 4)THEN


        if ((theta.gt.16.).and.(theta.lt.44.)) then
	 x=1.767244-.2542930*(theta)+0.1246248E-01*((theta)**2)-.1430444E-03*((theta)**3)
     &-.6262982E-05*((theta)**4)+0.1960881E-06*((theta)**5)-.1566946E-08*((theta)**6)
	 y=-.1807822E-01+0.3758978E-02*(theta)-.2207941E-03*((theta)**2)+0.5200121E-05*((theta)**3)
     &-.4353816E-07*((theta)**4)
	 x2=0.4989405E-03-.8992112E-05*(theta)+0.1226031E-03*SIN(449./theta+1.448175)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	

	
	 else
	 fact=0
	 endif   	 
	
 
 
      ENDIF
      IF(nsector.eq. 5)THEN
      
      
        if ((theta.gt.16.).and.(theta.lt.32.5)) then
	 x=0.6389222E-01-.4801639E-03*(theta)
	 y=0.3421679E-02-.8440115E-04*(theta)
	 x2=0.3923497E-03-.4957983E-05*(theta)-.9631179E-04*SIN(449./theta-1.264236)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	
         elseif ((theta.gt.32.5).and.(theta.lt.44.)) then
	 x=24.81553-1.917555*(theta)+0.4922545E-01*((theta)**2)-.4198686E-03*((theta)**3)
	 y=17.04936-1.342118*(theta)+0.2436030E-01*((theta)**2)+0.5353051E-03*((theta)**3)
     &-.2185274E-04*((theta)**4)+0.1889016E-06*((theta)**5)
	 x2=0.2119622E-02-.5375929E-04*(theta)-.2181619E-03*SIN(449./theta-.1648105)
	 fact=x+y*(phisec)+x2*((phisec)**2) 

	
	 
	
	 else
	 fact=0
	 endif   	
      
 
       ENDIF
       
       
      IF(nsector.eq. 6)THEN
      
       if ((theta.gt.16.).and.(theta.lt.32.5)) then
	 x=0.2883536E-01-.4278146E-04*(theta)
	 y=0.4266386E-02-.7680542E-04*(theta)
	 x2=0.2714948E-03+0.6075697E-06*(theta)+0.7332848E-04*SIN(449./theta+8.132975)
	 fact=x+y*(phisec)+x2*((phisec)**2)
	
	
         elseif ((theta.gt.32.5).and.(theta.lt.44.)) then
	 x=8.539058-.6286885*(theta)+0.1536358E-01*((theta)**2)-.1245493E-03*((theta)**3)
	 y=-.1018162E-02+0.6171986E-04*(theta)+0.8155919E-06*((theta)**2)+0.1453145E-07*((theta)**3)
     &+0.2333810E-09*((theta)**4)-.2149097E-10*((theta)**5)
	 x2=0.1190329E-02-.2835207E-04*(theta)-.1359806E-03*SIN(449./theta-.8999598)
	 fact=x+y*(phisec)+x2*((phisec)**2) 

	
	 
	
	 else
	 fact=0
	 endif         
      

      ENDIF
       
      correct_thel_e1_2039_2250_feb09=theta-fact
      RETURN
      END 





      FUNCTION correct_pel_e1_2039_2250_feb09(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_2039_2250_feb09
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN

	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=1.008854-.3422894E-03*(theta)-.6954598E-06*((theta)**2)+0.8156167E-07*((theta)**3)
	 y=-.1025022E-02+0.8506941E-05*(theta)+0.1281147E-05*((theta)**2)-.2162752E-07*((theta)**3)
	 x2=-.9028474E-04+0.8689741E-05*(theta)-.2580860E-06*((theta)**2)+0.2403976E-08*((theta)**3)
	 fact=(x+y*(phisec)+x2*((phisec)**2))
	
	
	 
	
	 else
	 fact=1.
	 endif  
	    
      ENDIF
      
      IF(nsector.eq. 2)THEN
      
	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=0.9984451+0.5805102E-04*(theta)-.2909986E-05*((theta)**2)+0.4793942E-07*((theta)**3)
	 y=-.1578927E-02+0.9188605E-04*(theta)-.1745132E-05*((theta)**2)+0.8876872E-08*((theta)**3)
	 x2=-.4820524E-04+0.5224319E-05*(theta)-.1622990E-06*((theta)**2)+0.1537213E-08*((theta)**3)
	 fact=(x+y*(phisec)+x2*((phisec)**2))
	
	
	 
	
	 else
	 fact=1.
	 endif        
      

      ENDIF
      
      
      IF(nsector.eq. 3)THEN
      
	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=1.017301-.2568450E-02*(theta)+0.1060422E-03*((theta)**2)-.1080502E-05*((theta)**3)
     &-.1366504E-07*((theta)**4)+0.2249102E-09*((theta)**5)
	 y=-.2755583E-02+0.2811991E-03*(theta)-.6515745E-05*((theta)**2)-.5426127E-07*((theta)**3)
     &+0.3184020E-08*((theta)**4)-.2669610E-10*((theta)**5)
	 x2=-.8553622E-04+0.8231474E-05*(theta)-.2452097E-06*((theta)**2)+0.2279487E-08*((theta)**3)
	 fact=(x+y*(phisec)+x2*((phisec)**2))
	
	
	 
	
	 else
	 fact=1.
	 endif        
            
      
    
      ENDIF
      
      
      
      IF(nsector.eq. 4)THEN
 
	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=1.016233-.7399115E-03*(theta)+0.8491787E-05*((theta)**2)
	 y=-.3036835E-02+0.3991334E-03*(theta)-.1802434E-04*((theta)**2)+0.3370880E-06*((theta)**3)
     &-.2288362E-08*((theta)**4)
	 x2=-.4231546E-04+0.4222708E-05*(theta)-.1275686E-06*((theta)**2)+0.1193804E-08*((theta)**3)
	 fact=(x+y*(phisec)+x2*((phisec)**2))
	
	
	 
	
	 else
	 fact=1.
	 endif   
 
 
 
      ENDIF
      IF(nsector.eq. 5)THEN
      
	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=0.9879246+0.1914278E-02*(theta)-.7166770E-04*((theta)**2)+0.8048290E-06*((theta)**3)
	 y=-.1823859E-02+0.1565340E-03*(theta)-.4768085E-05*((theta)**2)+0.4759397E-07*((theta)**3)
	 x2=-.5335819E-04+0.5686952E-05*(theta)-.1885548E-06*((theta)**2)+0.1940291E-08*((theta)**3)
	 y2= -.2186460E-07-.6860175E-09*(theta)
	 d=0.5041720E-08+0.1139555E-09*(theta)
	 fact=(x+y*(phisec)+x2*((phisec)**2)+y2*((phisec)**3)+d*((phisec)**4))
	
	
	 
	
	 else
	 fact=1.
	 endif        
      
 
       ENDIF
       
       
      IF(nsector.eq. 6)THEN
      
      
	 if ((theta.gt.16.).and.(theta.lt.52.)) then
	 
	 x=1.007478-.2812020E-03*(theta)+0.3874791E-05*((theta)**2)
	 y=-.7835627E-03+0.9022003E-04*(theta)-.2670348E-05*((theta)**2)+0.2166390E-07*((theta)**3)
	 x2=0.3509918E-05-.1477360E-06*(theta)
	 fact=(x+y*(phisec)+x2*((phisec)**2))
	
	
	 
	
	 else
	 fact=1.
	 endif         
      

      ENDIF
       
      correct_pel_e1_2039_2250_feb09=fact*pe
      RETURN
      END 
      
            FUNCTION correct_pel_e1_1515_1500_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_1515_1500_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.99358E+00+(-0.70804E-03*x)+(-0.30910E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99358E+00+(-0.70804E-03*(x2-d))+(-0.30910E-05*(x2-d)**2)
       f2= 0.99390E+00+(-0.50093E-03*(x2+d))+( 0.36917E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99390E+00+(-0.50093E-03*x)+( 0.36917E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.88308E+00+( 0.10843E-01*y)+(-0.25609E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.88308E+00+( 0.10843E-01*(y2-d))+(-0.25609E-03*(y2-d)**2)
       f2= 0.10051E+01+(-0.32993E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10051E+01+(-0.32993E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98904E+00+(-0.48208E-03*x)+( 0.54538E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98904E+00+(-0.48208E-03*(x2-d))+( 0.54538E-05*(x2-d)**2)
       f2= 0.98966E+00+(-0.29416E-03*(x2+d))+( 0.18120E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98966E+00+(-0.29416E-03*x)+( 0.18120E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.85975E+00+( 0.13069E-01*y)+(-0.32241E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.85975E+00+( 0.13069E-01*(y2-d))+(-0.32241E-03*(y2-d)**2)
       f2= 0.99546E+00+(-0.13865E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99546E+00+(-0.13865E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99416E+00+(-0.42043E-03*x)+(-0.21998E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99416E+00+(-0.42043E-03*(x2-d))+(-0.21998E-04*(x2-d)**2)
       f2= 0.99550E+00+( 0.18387E-03*(x2+d))+( 0.15942E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99550E+00+( 0.18387E-03*x)+( 0.15942E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95092E+00+( 0.36291E-02*y)+(-0.65581E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95092E+00+( 0.36291E-02*(y2-d))+(-0.65581E-04*(y2-d)**2)
       f2= 0.99638E+00+( 0.76828E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99638E+00+( 0.76828E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99433E+00+( 0.81421E-04*x)+( 0.22693E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99433E+00+( 0.81421E-04*(x2-d))+( 0.22693E-04*(x2-d)**2)
       f2= 0.99453E+00+( 0.10881E-03*(x2+d))+( 0.85190E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99453E+00+( 0.10881E-03*x)+( 0.85190E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92787E+00+( 0.59262E-02*y)+(-0.12430E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92787E+00+( 0.59262E-02*(y2-d))+(-0.12430E-03*(y2-d)**2)
       f2= 0.99900E+00+(-0.83302E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99900E+00+(-0.83302E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99398E+00+(-0.32676E-03*x)+( 0.16235E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99398E+00+(-0.32676E-03*(x2-d))+( 0.16235E-04*(x2-d)**2)
       f2= 0.99454E+00+(-0.95435E-04*(x2+d))+( 0.16601E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99454E+00+(-0.95435E-04*x)+( 0.16601E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89510E+00+( 0.93581E-02*y)+(-0.21302E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89510E+00+( 0.93581E-02*(y2-d))+(-0.21302E-03*(y2-d)**2)
       f2= 0.99786E+00+( 0.11810E-05*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99786E+00+( 0.11810E-05*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99461E+00+(-0.11352E-03*x)+( 0.93415E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99461E+00+(-0.11352E-03*(x2-d))+( 0.93415E-05*(x2-d)**2)
       f2= 0.99526E+00+( 0.25471E-03*(x2+d))+( 0.10319E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99526E+00+( 0.25471E-03*x)+( 0.10319E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89672E+00+( 0.96887E-02*y)+(-0.23380E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89672E+00+( 0.96887E-02*(y2-d))+(-0.23380E-03*(y2-d)**2)
       f2= 0.99453E+00+( 0.12593E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99453E+00+( 0.12593E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_1515_1500_jan99=fact*pe
      RETURN
      END 
      
      
      

      FUNCTION correct_pel_e1_1515_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_1515_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99480E+00+(-0.50051E-03*x)+(-0.15416E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99480E+00+(-0.50051E-03*(x2-d))+(-0.15416E-04*(x2-d)**2)
       f2= 0.99515E+00+(-0.26415E-03*(x2+d))+( 0.12705E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99515E+00+(-0.26415E-03*x)+( 0.12705E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93232E+00+( 0.54104E-02*y)+(-0.11299E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93232E+00+( 0.54104E-02*(y2-d))+(-0.11299E-03*(y2-d)**2)
       f2= 0.10041E+01+(-0.28332E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10041E+01+(-0.28332E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99123E+00+(-0.31776E-03*x)+(-0.15075E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99123E+00+(-0.31776E-03*(x2-d))+(-0.15075E-05*(x2-d)**2)
       f2= 0.99166E+00+(-0.18941E-03*(x2+d))+( 0.78313E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99166E+00+(-0.18941E-03*x)+( 0.78313E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.93421E+00+( 0.48929E-02*y)+(-0.10187E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.93421E+00+( 0.48929E-02*(y2-d))+(-0.10187E-03*(y2-d)**2)
       f2= 0.99460E+00+(-0.71288E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99460E+00+(-0.71288E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99594E+00+( 0.67596E-04*x)+( 0.27742E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99594E+00+( 0.67596E-04*(x2-d))+( 0.27742E-05*(x2-d)**2)
       f2= 0.99606E+00+( 0.61163E-04*(x2+d))+( 0.10389E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99606E+00+( 0.61163E-04*x)+( 0.10389E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97258E+00+( 0.15576E-02*y)+(-0.22705E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97258E+00+( 0.15576E-02*(y2-d))+(-0.22705E-04*(y2-d)**2)
       f2= 0.99818E+00+(-0.27000E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99818E+00+(-0.27000E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99495E+00+(-0.16044E-03*x)+(-0.21930E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99495E+00+(-0.16044E-03*(x2-d))+(-0.21930E-05*(x2-d)**2)
       f2= 0.99536E+00+( 0.60063E-04*(x2+d))+(-0.11998E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99536E+00+( 0.60063E-04*x)+(-0.11998E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95635E+00+( 0.32390E-02*y)+(-0.65342E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95635E+00+( 0.32390E-02*(y2-d))+(-0.65342E-04*(y2-d)**2)
       f2= 0.99768E+00+(-0.54881E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99768E+00+(-0.54881E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99587E+00+(-0.22965E-03*x)+( 0.32447E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99587E+00+(-0.22965E-03*(x2-d))+( 0.32447E-05*(x2-d)**2)
       f2= 0.99606E+00+(-0.80708E-04*(x2+d))+( 0.25642E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99606E+00+(-0.80708E-04*x)+( 0.25642E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92313E+00+( 0.61403E-02*y)+(-0.12674E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92313E+00+( 0.61403E-02*(y2-d))+(-0.12674E-03*(y2-d)**2)
       f2= 0.99873E+00+(-0.55082E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99873E+00+(-0.55082E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 26.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99629E+00+(-0.40582E-04*x)+( 0.44242E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99629E+00+(-0.40582E-04*(x2-d))+( 0.44242E-05*(x2-d)**2)
       f2= 0.99639E+00+( 0.93473E-04*(x2+d))+(-0.48949E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99639E+00+( 0.93473E-04*x)+(-0.48949E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.95626E+00+( 0.33530E-02*y)+(-0.69105E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.95626E+00+( 0.33530E-02*(y2-d))+(-0.69105E-04*(y2-d)**2)
       f2= 0.99551E+00+( 0.58111E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99551E+00+( 0.58111E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_1515_2250_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_pel_e1_2567_1500_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_2567_1500_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.98929E+00+(-0.10213E-02*x)+( 0.15939E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98929E+00+(-0.10213E-02*(x2-d))+( 0.15939E-04*(x2-d)**2)
       f2= 0.99022E+00+(-0.63859E-03*(x2+d))+( 0.61665E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99022E+00+(-0.63859E-03*x)+( 0.61665E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.84004E+00+( 0.16942E-01*y)+(-0.45086E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.84004E+00+( 0.16942E-01*(y2-d))+(-0.45086E-03*(y2-d)**2)
       f2= 0.10094E+01+(-0.53617E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10094E+01+(-0.53617E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98543E+00+(-0.56690E-03*x)+( 0.30922E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98543E+00+(-0.56690E-03*(x2-d))+( 0.30922E-04*(x2-d)**2)
       f2= 0.98583E+00+(-0.78940E-03*(x2+d))+( 0.57892E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98583E+00+(-0.78940E-03*x)+( 0.57892E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87348E+00+( 0.12117E-01*y)+(-0.30786E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87348E+00+( 0.12117E-01*(y2-d))+(-0.30786E-03*(y2-d)**2)
       f2= 0.10011E+01+(-0.41944E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10011E+01+(-0.41944E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99386E+00+( 0.48521E-04*x)+( 0.14338E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99386E+00+( 0.48521E-04*(x2-d))+( 0.14338E-04*(x2-d)**2)
       f2= 0.99453E+00+(-0.21200E-04*(x2+d))+( 0.48794E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99453E+00+(-0.21200E-04*x)+( 0.48794E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94100E+00+( 0.51300E-02*y)+(-0.10759E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94100E+00+( 0.51300E-02*(y2-d))+(-0.10759E-03*(y2-d)**2)
       f2= 0.10008E+01+( 0.35398E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10008E+01+( 0.35398E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98987E+00+(-0.18871E-03*x)+( 0.28448E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98987E+00+(-0.18871E-03*(x2-d))+( 0.28448E-04*(x2-d)**2)
       f2= 0.99053E+00+( 0.10135E-03*(x2+d))+( 0.41065E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99053E+00+( 0.10135E-03*x)+( 0.41065E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87771E+00+( 0.11893E-01*y)+(-0.29241E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87771E+00+( 0.11893E-01*(y2-d))+(-0.29241E-03*(y2-d)**2)
       f2= 0.10008E+01+(-0.10021E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10008E+01+(-0.10021E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98807E+00+(-0.13800E-02*x)+(-0.42142E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98807E+00+(-0.13800E-02*(x2-d))+(-0.42142E-05*(x2-d)**2)
       f2= 0.98996E+00+(-0.57027E-03*(x2+d))+( 0.84417E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98996E+00+(-0.57027E-03*x)+( 0.84417E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.86659E+00+( 0.12910E-01*y)+(-0.31083E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.86659E+00+( 0.12910E-01*(y2-d))+(-0.31083E-03*(y2-d)**2)
       f2= 0.10017E+01+(-0.57903E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10017E+01+(-0.57903E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99188E+00+( 0.19937E-03*x)+( 0.42159E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99188E+00+( 0.19937E-03*(x2-d))+( 0.42159E-04*(x2-d)**2)
       f2= 0.99320E+00+( 0.73775E-03*(x2+d))+( 0.25390E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99320E+00+( 0.73775E-03*x)+( 0.25390E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90534E+00+( 0.87790E-02*y)+(-0.20532E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90534E+00+( 0.87790E-02*(y2-d))+(-0.20532E-03*(y2-d)**2)
       f2= 0.99694E+00+( 0.10974E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99694E+00+( 0.10974E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_2567_1500_jan99=fact*pe
      RETURN
      END
      FUNCTION correct_pel_e1_2567_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_2567_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.99274E+00+(-0.81834E-03*x)+(-0.10034E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99274E+00+(-0.81834E-03*(x2-d))+(-0.10034E-05*(x2-d)**2)
       f2= 0.99332E+00+(-0.49191E-03*(x2+d))+( 0.43976E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99332E+00+(-0.49191E-03*x)+( 0.43976E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.83258E+00+( 0.16832E-01*y)+(-0.42637E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.83258E+00+( 0.16832E-01*(y2-d))+(-0.42637E-03*(y2-d)**2)
       f2= 0.10055E+01+(-0.33414E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10055E+01+(-0.33414E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99059E+00+(-0.50965E-03*x)+( 0.16964E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99059E+00+(-0.50965E-03*(x2-d))+( 0.16964E-04*(x2-d)**2)
       f2= 0.99061E+00+(-0.50617E-03*(x2+d))+( 0.20182E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99061E+00+(-0.50617E-03*x)+( 0.20182E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10617E+01+(-0.82458E-02*y)+( 0.23853E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10617E+01+(-0.82458E-02*(y2-d))+( 0.23853E-03*(y2-d)**2)
       f2= 0.99692E+00+(-0.18066E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99692E+00+(-0.18066E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99442E+00+(-0.25783E-03*x)+(-0.77472E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99442E+00+(-0.25783E-03*(x2-d))+(-0.77472E-05*(x2-d)**2)
       f2= 0.99509E+00+( 0.10588E-03*(x2+d))+( 0.28290E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99509E+00+( 0.10588E-03*x)+( 0.28290E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92951E+00+( 0.63321E-02*y)+(-0.14405E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92951E+00+( 0.63321E-02*(y2-d))+(-0.14405E-03*(y2-d)**2)
       f2= 0.99695E+00+( 0.95303E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99695E+00+( 0.95303E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99347E+00+( 0.14709E-04*x)+( 0.25219E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99347E+00+( 0.14709E-04*(x2-d))+( 0.25219E-04*(x2-d)**2)
       f2= 0.99379E+00+( 0.97314E-04*(x2+d))+( 0.19761E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99379E+00+( 0.97314E-04*x)+( 0.19761E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87811E+00+( 0.11598E-01*y)+(-0.28046E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87811E+00+( 0.11598E-01*(y2-d))+(-0.28046E-03*(y2-d)**2)
       f2= 0.99963E+00+(-0.77553E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99963E+00+(-0.77553E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99308E+00+(-0.50765E-03*x)+( 0.12877E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99308E+00+(-0.50765E-03*(x2-d))+( 0.12877E-04*(x2-d)**2)
       f2= 0.99394E+00+(-0.60055E-04*(x2+d))+( 0.20159E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99394E+00+(-0.60055E-04*x)+( 0.20159E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.86189E+00+( 0.13096E-01*y)+(-0.31349E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.86189E+00+( 0.13096E-01*(y2-d))+(-0.31349E-03*(y2-d)**2)
       f2= 0.99901E+00+(-0.15166E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99901E+00+(-0.15166E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99360E+00+(-0.19166E-03*x)+( 0.93554E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99360E+00+(-0.19166E-03*(x2-d))+( 0.93554E-05*(x2-d)**2)
       f2= 0.99491E+00+( 0.35804E-03*(x2+d))+( 0.12265E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99491E+00+( 0.35804E-03*x)+( 0.12265E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91973E+00+( 0.71910E-02*y)+(-0.16412E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91973E+00+( 0.71910E-02*(y2-d))+(-0.16412E-03*(y2-d)**2)
       f2= 0.99562E+00+( 0.12486E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99562E+00+( 0.12486E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      correct_pel_e1_2567_2250_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_pel_e1_4056_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4056_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98879E+00+(-0.10352E-02*x)+( 0.22826E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98879E+00+(-0.10352E-02*(x2-d))+( 0.22826E-04*(x2-d)**2)
       f2= 0.99126E+00+(-0.30452E-03*(x2+d))+( 0.39902E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99126E+00+(-0.30452E-03*x)+( 0.39902E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92650E+00+( 0.70558E-02*y)+(-0.17180E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92650E+00+( 0.70558E-02*(y2-d))+(-0.17180E-03*(y2-d)**2)
       f2= 0.10034E+01+(-0.25153E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10034E+01+(-0.25153E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98377E+00+(-0.11670E-02*x)+( 0.95386E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98377E+00+(-0.11670E-02*(x2-d))+( 0.95386E-05*(x2-d)**2)
       f2= 0.98485E+00+(-0.11419E-02*(x2+d))+( 0.10236E-03*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98485E+00+(-0.11419E-02*x)+( 0.10236E-03*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97716E+00+(-0.96167E-03*y)+( 0.10211E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97716E+00+(-0.96167E-03*(y2-d))+( 0.10211E-03*(y2-d)**2)
       f2= 0.99743E+00+(-0.14730E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99743E+00+(-0.14730E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99293E+00+(-0.29953E-03*x)+( 0.84805E-06*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99293E+00+(-0.29953E-03*(x2-d))+( 0.84805E-06*(x2-d)**2)
       f2= 0.99494E+00+(-0.13311E-03*(x2+d))+( 0.56717E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99494E+00+(-0.13311E-03*x)+( 0.56717E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94440E+00+( 0.50527E-02*y)+(-0.11059E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94440E+00+( 0.50527E-02*(y2-d))+(-0.11059E-03*(y2-d)**2)
       f2= 0.99998E+00+( 0.98663E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99998E+00+( 0.98663E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99004E+00+(-0.64608E-03*x)+(-0.11954E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99004E+00+(-0.64608E-03*(x2-d))+(-0.11954E-05*(x2-d)**2)
       f2= 0.99072E+00+(-0.16336E-03*(x2+d))+( 0.59020E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99072E+00+(-0.16336E-03*x)+( 0.59020E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89145E+00+( 0.10518E-01*y)+(-0.25515E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89145E+00+( 0.10518E-01*(y2-d))+(-0.25515E-03*(y2-d)**2)
       f2= 0.10005E+01+(-0.66966E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10005E+01+(-0.66966E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98838E+00+(-0.22803E-02*x)+(-0.74725E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98838E+00+(-0.22803E-02*(x2-d))+(-0.74725E-04*(x2-d)**2)
       f2= 0.98990E+00+(-0.66900E-03*(x2+d))+( 0.95871E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98990E+00+(-0.66900E-03*x)+( 0.95871E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.86696E+00+( 0.13049E-01*y)+(-0.31582E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.86696E+00+( 0.13049E-01*(y2-d))+(-0.31582E-03*(y2-d)**2)
       f2= 0.10032E+01+(-0.10833E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10032E+01+(-0.10833E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 23.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99018E+00+(-0.12737E-03*x)+( 0.30811E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99018E+00+(-0.12737E-03*(x2-d))+( 0.30811E-04*(x2-d)**2)
       f2= 0.99362E+00+( 0.10262E-02*(x2+d))+(-0.21133E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99362E+00+( 0.10262E-02*x)+(-0.21133E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89892E+00+( 0.96272E-02*y)+(-0.22791E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89892E+00+( 0.96272E-02*(y2-d))+(-0.22791E-03*(y2-d)**2)
       f2= 0.99767E+00+( 0.12376E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99767E+00+( 0.12376E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4056_2250_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_pel_e1_4247_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4247_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
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
       fact=fact*( 0.98703E+00+(-0.15768E-02*x)+(-0.11675E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98703E+00+(-0.15768E-02*(x2-d))+(-0.11675E-04*(x2-d)**2)
       f2= 0.99009E+00+(-0.16724E-04*(x2+d))+( 0.14757E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99009E+00+(-0.16724E-04*x)+( 0.14757E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.90269E+00+( 0.99216E-02*y)+(-0.25834E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.90269E+00+( 0.99216E-02*(y2-d))+(-0.25834E-03*(y2-d)**2)
       f2= 0.10098E+01+(-0.54307E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10098E+01+(-0.54307E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98396E+00+(-0.11489E-02*x)+(-0.14769E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98396E+00+(-0.11489E-02*(x2-d))+(-0.14769E-05*(x2-d)**2)
       f2= 0.98652E+00+(-0.62563E-03*(x2+d))+( 0.38008E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98652E+00+(-0.62563E-03*x)+( 0.38008E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10302E+01+(-0.55750E-02*y)+( 0.18343E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10302E+01+(-0.55750E-02*(y2-d))+( 0.18343E-03*(y2-d)**2)
       f2= 0.10015E+01+(-0.42511E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10015E+01+(-0.42511E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99408E+00+( 0.16432E-04*x)+( 0.13580E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99408E+00+( 0.16432E-04*(x2-d))+( 0.13580E-04*(x2-d)**2)
       f2= 0.99454E+00+( 0.19289E-05*(x2+d))+( 0.41582E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99454E+00+( 0.19289E-05*x)+( 0.41582E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.92781E+00+( 0.71373E-02*y)+(-0.17456E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.92781E+00+( 0.71373E-02*(y2-d))+(-0.17456E-03*(y2-d)**2)
       f2= 0.99968E+00+( 0.96557E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99968E+00+( 0.96557E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98771E+00+(-0.13452E-02*x)+(-0.44061E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98771E+00+(-0.13452E-02*(x2-d))+(-0.44061E-04*(x2-d)**2)
       f2= 0.98857E+00+( 0.16115E-03*(x2+d))+( 0.53378E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98857E+00+( 0.16115E-03*x)+( 0.53378E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87759E+00+( 0.12311E-01*y)+(-0.31133E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87759E+00+( 0.12311E-01*(y2-d))+(-0.31133E-03*(y2-d)**2)
       f2= 0.10020E+01+(-0.11086E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10020E+01+(-0.11086E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98138E+00+(-0.32661E-02*x)+(-0.11254E-03*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98138E+00+(-0.32661E-02*(x2-d))+(-0.11254E-03*(x2-d)**2)
       f2= 0.99115E+00+(-0.26282E-03*(x2+d))+( 0.53665E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99115E+00+(-0.26282E-03*x)+( 0.53665E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.85622E+00+( 0.14354E-01*y)+(-0.35590E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.85622E+00+( 0.14354E-01*(y2-d))+(-0.35590E-03*(y2-d)**2)
       f2= 0.99988E+00+( 0.46599E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99988E+00+( 0.46599E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 22.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99015E+00+( 0.27665E-03*x)+( 0.62844E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99015E+00+( 0.27665E-03*(x2-d))+( 0.62844E-04*(x2-d)**2)
       f2= 0.99403E+00+( 0.12122E-02*(x2+d))+(-0.35037E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99403E+00+( 0.12122E-02*x)+(-0.35037E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89151E+00+( 0.10655E-01*y)+(-0.26083E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89151E+00+( 0.10655E-01*(y2-d))+(-0.26083E-03*(y2-d)**2)
       f2= 0.10003E+01+( 0.58332E-05*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10003E+01+( 0.58332E-05*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4247_2250_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_pel_e1_4462_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4462_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98869E+00+(-0.12226E-02*x)+( 0.11178E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98869E+00+(-0.12226E-02*(x2-d))+( 0.11178E-04*(x2-d)**2)
       f2= 0.99094E+00+(-0.10390E-03*(x2+d))+( 0.63005E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99094E+00+(-0.10390E-03*x)+( 0.63005E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.80357E+00+( 0.22703E-01*y)+(-0.65378E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.80357E+00+( 0.22703E-01*(y2-d))+(-0.65378E-03*(y2-d)**2)
       f2= 0.10096E+01+(-0.54993E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10096E+01+(-0.54993E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98222E+00+(-0.14293E-02*x)+( 0.76271E-06*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98222E+00+(-0.14293E-02*(x2-d))+( 0.76271E-06*(x2-d)**2)
       f2= 0.98285E+00+(-0.98607E-03*(x2+d))+( 0.10482E-03*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98285E+00+(-0.98607E-03*x)+( 0.10482E-03*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.79900E+00+( 0.22173E-01*y)+(-0.63258E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.79900E+00+( 0.22173E-01*(y2-d))+(-0.63258E-03*(y2-d)**2)
       f2= 0.99854E+00+(-0.24502E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99854E+00+(-0.24502E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99448E+00+(-0.87876E-04*x)+( 0.84028E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99448E+00+(-0.87876E-04*(x2-d))+( 0.84028E-05*(x2-d)**2)
       f2= 0.99507E+00+(-0.75488E-04*(x2+d))+( 0.49586E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99507E+00+(-0.75488E-04*x)+( 0.49586E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97330E+00+( 0.10391E-02*y)+( 0.32215E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97330E+00+( 0.10391E-02*(y2-d))+( 0.32215E-04*(y2-d)**2)
       f2= 0.10016E+01+( 0.80236E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10016E+01+( 0.80236E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99113E+00+(-0.81657E-03*x)+(-0.17385E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99113E+00+(-0.81657E-03*(x2-d))+(-0.17385E-04*(x2-d)**2)
       f2= 0.99115E+00+(-0.30546E-03*(x2+d))+( 0.73645E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99115E+00+(-0.30546E-03*x)+( 0.73645E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.86087E+00+( 0.14520E-01*y)+(-0.38113E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.86087E+00+( 0.14520E-01*(y2-d))+(-0.38113E-03*(y2-d)**2)
       f2= 0.99815E+00+( 0.88851E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99815E+00+( 0.88851E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98734E+00+(-0.25988E-02*x)+(-0.88036E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98734E+00+(-0.25988E-02*(x2-d))+(-0.88036E-04*(x2-d)**2)
       f2= 0.99058E+00+(-0.39750E-03*(x2+d))+( 0.69740E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99058E+00+(-0.39750E-03*x)+( 0.69740E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.78456E+00+( 0.23348E-01*y)+(-0.63096E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.78456E+00+( 0.23348E-01*(y2-d))+(-0.63096E-03*(y2-d)**2)
       f2= 0.10023E+01+(-0.24078E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10023E+01+(-0.24078E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99288E+00+( 0.56354E-03*x)+( 0.66952E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99288E+00+( 0.56354E-03*(x2-d))+( 0.66952E-04*(x2-d)**2)
       f2= 0.99493E+00+( 0.10258E-02*(x2+d))+(-0.17564E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99493E+00+( 0.10258E-02*x)+(-0.17564E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87682E+00+( 0.12641E-01*y)+(-0.32642E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87682E+00+( 0.12641E-01*(y2-d))+(-0.32642E-03*(y2-d)**2)
       f2= 0.99447E+00+( 0.27042E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99447E+00+( 0.27042E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4462_2250_jan99=fact*pe
      RETURN
      END 


      FUNCTION correct_pel_e1_4462_3375_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4462_3375_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99215E+00+(-0.12418E-02*x)+(-0.22354E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99215E+00+(-0.12418E-02*(x2-d))+(-0.22354E-04*(x2-d)**2)
       f2= 0.99337E+00+(-0.49330E-03*(x2+d))+( 0.43769E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99337E+00+(-0.49330E-03*x)+( 0.43769E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98895E+00+(-0.12893E-02*y)+( 0.10001E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98895E+00+(-0.12893E-02*(y2-d))+( 0.10001E-03*(y2-d)**2)
       f2= 0.99923E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99923E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98786E+00+(-0.73576E-03*x)+( 0.57830E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98786E+00+(-0.73576E-03*(x2-d))+( 0.57830E-05*(x2-d)**2)
       f2= 0.98876E+00+(-0.29576E-03*(x2+d))+( 0.22245E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98876E+00+(-0.29576E-03*x)+( 0.22245E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94052E+00+( 0.37950E-02*y)+(-0.49032E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94052E+00+( 0.37950E-02*(y2-d))+(-0.49032E-04*(y2-d)**2)
       f2= 0.99481E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99481E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99516E+00+(-0.24782E-03*x)+(-0.11505E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99516E+00+(-0.24782E-03*(x2-d))+(-0.11505E-04*(x2-d)**2)
       f2= 0.99563E+00+( 0.13994E-04*(x2+d))+( 0.32670E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99563E+00+( 0.13994E-04*x)+( 0.32670E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98824E+00+(-0.22244E-03*y)+( 0.42104E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98824E+00+(-0.22244E-03*(y2-d))+( 0.42104E-04*(y2-d)**2)
       f2= 0.10002E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99241E+00+(-0.41355E-03*x)+( 0.55864E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99241E+00+(-0.41355E-03*(x2-d))+( 0.55864E-05*(x2-d)**2)
       f2= 0.99360E+00+( 0.27628E-03*(x2+d))+( 0.10840E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99360E+00+( 0.27628E-03*x)+( 0.10840E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.89554E+00+( 0.10216E-01*y)+(-0.25268E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.89554E+00+( 0.10216E-01*(y2-d))+(-0.25268E-03*(y2-d)**2)
       f2= 0.99875E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99875E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99092E+00+(-0.89474E-03*x)+( 0.32688E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99092E+00+(-0.89474E-03*(x2-d))+( 0.32688E-05*(x2-d)**2)
       f2= 0.99274E+00+( 0.16606E-03*(x2+d))+( 0.11688E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99274E+00+( 0.16606E-03*x)+( 0.11688E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.87077E+00+( 0.12779E-01*y)+(-0.31799E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.87077E+00+( 0.12779E-01*(y2-d))+(-0.31799E-03*(y2-d)**2)
       f2= 0.99948E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99948E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99384E+00+(-0.29706E-04*x)+( 0.17395E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99384E+00+(-0.29706E-04*(x2-d))+( 0.17395E-04*(x2-d)**2)
       f2= 0.99477E+00+( 0.42899E-03*(x2+d))+( 0.17642E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99477E+00+( 0.42899E-03*x)+( 0.17642E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.91953E+00+( 0.75317E-02*y)+(-0.17850E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.91953E+00+( 0.75317E-02*(y2-d))+(-0.17850E-03*(y2-d)**2)
       f2= 0.99882E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99882E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4462_3375_jan99=fact*pe
      RETURN
      END 
      FUNCTION correct_pel_e1_4817_3375_feb00(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_pel_e1_4817_3375_feb00
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 24.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10030E+01+(-0.81772E-03*x)+(-0.27991E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10030E+01+(-0.81772E-03*(x2-d))+(-0.27991E-04*(x2-d)**2)
       f2= 0.10027E+01+(-0.74268E-03*(x2+d))+( 0.99415E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10027E+01+(-0.74268E-03*x)+( 0.99415E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98987E+00+( 0.21478E-02*y)+(-0.73873E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98987E+00+( 0.21478E-02*(y2-d))+(-0.73873E-04*(y2-d)**2)
       f2= 0.10053E+01+(-0.21693E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10053E+01+(-0.21693E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99970E+00+(-0.72346E-03*x)+(-0.24546E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99970E+00+(-0.72346E-03*(x2-d))+(-0.24546E-04*(x2-d)**2)
       f2= 0.10002E+01+(-0.62001E-03*(x2+d))+( 0.95569E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+(-0.62001E-03*x)+( 0.95569E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10119E+01+(-0.15024E-02*y)+( 0.48769E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10119E+01+(-0.15024E-02*(y2-d))+( 0.48769E-04*(y2-d)**2)
       f2= 0.10064E+01+(-0.30816E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10064E+01+(-0.30816E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10063E+01+( 0.13232E-03*x)+(-0.27079E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10063E+01+( 0.13232E-03*(x2-d))+(-0.27079E-04*(x2-d)**2)
       f2= 0.10056E+01+(-0.24762E-03*(x2+d))+( 0.17585E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10056E+01+(-0.24762E-03*x)+( 0.17585E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10791E+01+(-0.82152E-02*y)+( 0.22136E-03*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10791E+01+(-0.82152E-02*(y2-d))+( 0.22136E-03*(y2-d)**2)
       f2= 0.10036E+01+(-0.21944E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10036E+01+(-0.21944E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 24.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10029E+01+(-0.71455E-06*x)+(-0.20288E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10029E+01+(-0.71455E-06*(x2-d))+(-0.20288E-04*(x2-d)**2)
       f2= 0.10030E+01+( 0.61300E-04*(x2+d))+(-0.81581E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10030E+01+( 0.61300E-04*x)+(-0.81581E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10052E+01+( 0.79282E-04*y)+(-0.11984E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10052E+01+( 0.79282E-04*(y2-d))+(-0.11984E-04*(y2-d)**2)
       f2= 0.10041E+01+(-0.18767E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10041E+01+(-0.18767E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10019E+01+(-0.41460E-03*x)+(-0.13194E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10019E+01+(-0.41460E-03*(x2-d))+(-0.13194E-04*(x2-d)**2)
       f2= 0.10021E+01+(-0.36425E-03*(x2+d))+( 0.22327E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10021E+01+(-0.36425E-03*x)+( 0.22327E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10124E+01+(-0.12992E-02*y)+( 0.38558E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10124E+01+(-0.12992E-02*(y2-d))+( 0.38558E-04*(y2-d)**2)
       f2= 0.10041E+01+(-0.13179E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10041E+01+(-0.13179E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 20.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10043E+01+( 0.15604E-03*x)+(-0.76398E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10043E+01+( 0.15604E-03*(x2-d))+(-0.76398E-05*(x2-d)**2)
       f2= 0.10040E+01+(-0.59061E-04*(x2+d))+(-0.31230E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10040E+01+(-0.59061E-04*x)+(-0.31230E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99245E+00+( 0.17472E-02*y)+(-0.65471E-04*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99245E+00+( 0.17472E-02*(y2-d))+(-0.65471E-04*(y2-d)**2)
       f2= 0.10034E+01+(-0.57339E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10034E+01+(-0.57339E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_pel_e1_4817_3375_feb00=fact*pe
      RETURN
      END 







      FUNCTION correct_ppr_e1_1645_1500_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_1645_1500_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99907E+00+( 0.34934E-03*x)+( 0.47813E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99907E+00+( 0.34934E-03*(x2-d))+( 0.47813E-04*(x2-d)**2)
       f2= 0.99766E+00+(-0.51311E-03*(x2+d))+( 0.20451E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99766E+00+(-0.51311E-03*x)+( 0.20451E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10002E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99758E+00+(-0.72384E-03*x)+(-0.41007E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99758E+00+(-0.72384E-03*(x2-d))+(-0.41007E-05*(x2-d)**2)
       f2= 0.99854E+00+(-0.30075E-03*(x2+d))+(-0.16190E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99854E+00+(-0.30075E-03*x)+(-0.16190E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10009E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10009E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10026E+01+(-0.82469E-03*x)+(-0.25370E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10026E+01+(-0.82469E-03*(x2-d))+(-0.25370E-05*(x2-d)**2)
       f2= 0.10035E+01+(-0.29490E-03*(x2+d))+(-0.26247E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10035E+01+(-0.29490E-03*x)+(-0.26247E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10050E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10050E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99646E+00+(-0.32621E-03*x)+( 0.13486E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99646E+00+(-0.32621E-03*(x2-d))+( 0.13486E-04*(x2-d)**2)
       f2= 0.99658E+00+(-0.95559E-03*(x2+d))+( 0.53832E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99658E+00+(-0.95559E-03*x)+( 0.53832E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99908E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99908E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10007E+01+( 0.18420E-03*x)+( 0.35387E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10007E+01+( 0.18420E-03*(x2-d))+( 0.35387E-04*(x2-d)**2)
       f2= 0.99960E+00+(-0.56545E-03*(x2+d))+( 0.25515E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99960E+00+(-0.56545E-03*x)+( 0.25515E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10017E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10017E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99996E+00+( 0.10953E-02*x)+( 0.41515E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99996E+00+( 0.10953E-02*(x2-d))+( 0.41515E-04*(x2-d)**2)
       f2= 0.99769E+00+( 0.21443E-03*(x2+d))+( 0.15599E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99769E+00+( 0.21443E-03*x)+( 0.15599E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99940E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99940E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_1645_1500_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_1645_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_1645_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10002E+01+( 0.94584E-05*x)+( 0.18472E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10002E+01+( 0.94584E-05*(x2-d))+( 0.18472E-04*(x2-d)**2)
       f2= 0.10002E+01+(-0.61893E-03*(x2+d))+( 0.23122E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+(-0.61893E-03*x)+( 0.23122E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99995E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99995E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10009E+01+(-0.50269E-03*x)+(-0.53086E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10009E+01+(-0.50269E-03*(x2-d))+(-0.53086E-05*(x2-d)**2)
       f2= 0.10011E+01+(-0.45582E-03*(x2+d))+( 0.75743E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10011E+01+(-0.45582E-03*x)+( 0.75743E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10003E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10003E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10031E+01+(-0.10318E-02*x)+(-0.17542E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10031E+01+(-0.10318E-02*(x2-d))+(-0.17542E-04*(x2-d)**2)
       f2= 0.10036E+01+(-0.61225E-03*(x2+d))+( 0.37275E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10036E+01+(-0.61225E-03*x)+( 0.37275E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10035E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10035E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10006E+01+(-0.12949E-03*x)+( 0.16596E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10006E+01+(-0.12949E-03*(x2-d))+( 0.16596E-04*(x2-d)**2)
       f2= 0.99945E+00+(-0.83940E-03*(x2+d))+( 0.33722E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99945E+00+(-0.83940E-03*x)+( 0.33722E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99898E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99898E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10018E+01+(-0.17026E-03*x)+( 0.10608E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10018E+01+(-0.17026E-03*(x2-d))+( 0.10608E-04*(x2-d)**2)
       f2= 0.10015E+01+(-0.32912E-03*(x2+d))+( 0.22846E-06*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10015E+01+(-0.32912E-03*x)+( 0.22846E-06*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10016E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10016E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99871E+00+( 0.20174E-03*x)+( 0.70030E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99871E+00+( 0.20174E-03*(x2-d))+( 0.70030E-05*(x2-d)**2)
       f2= 0.99872E+00+( 0.38198E-03*(x2+d))+( 0.10246E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99872E+00+( 0.38198E-03*x)+( 0.10246E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10005E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10005E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_1645_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_2445_1500_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_2445_1500_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99863E+00+(-0.58915E-03*x)+(-0.66957E-06*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99863E+00+(-0.58915E-03*(x2-d))+(-0.66957E-06*(x2-d)**2)
       f2= 0.99942E+00+(-0.54260E-03*(x2+d))+( 0.56736E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99942E+00+(-0.54260E-03*x)+( 0.56736E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99787E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99787E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10028E+01+(-0.22065E-03*x)+( 0.42317E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10028E+01+(-0.22065E-03*(x2-d))+( 0.42317E-05*(x2-d)**2)
       f2= 0.10019E+01+(-0.52251E-03*(x2+d))+( 0.18377E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10019E+01+(-0.52251E-03*x)+( 0.18377E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10002E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10077E+01+(-0.28206E-03*x)+( 0.49319E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10077E+01+(-0.28206E-03*(x2-d))+( 0.49319E-05*(x2-d)**2)
       f2= 0.10070E+01+(-0.55169E-03*(x2+d))+(-0.34057E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10070E+01+(-0.55169E-03*x)+(-0.34057E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10053E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10053E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99835E+00+(-0.26862E-03*x)+( 0.11162E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99835E+00+(-0.26862E-03*(x2-d))+( 0.11162E-04*(x2-d)**2)
       f2= 0.99757E+00+(-0.36902E-03*(x2+d))+( 0.12471E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99757E+00+(-0.36902E-03*x)+( 0.12471E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99727E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99727E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10007E+01+(-0.27331E-03*x)+( 0.82722E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10007E+01+(-0.27331E-03*(x2-d))+( 0.82722E-05*(x2-d)**2)
       f2= 0.10010E+01+(-0.40102E-03*(x2+d))+( 0.17853E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10010E+01+(-0.40102E-03*x)+( 0.17853E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10002E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99723E+00+( 0.10962E-02*x)+( 0.43432E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99723E+00+( 0.10962E-02*(x2-d))+( 0.43432E-04*(x2-d)**2)
       f2= 0.99653E+00+( 0.65454E-03*(x2+d))+(-0.10029E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99653E+00+( 0.65454E-03*x)+(-0.10029E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99770E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99770E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_2445_1500_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_2445_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_2445_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99898E+00+(-0.54724E-03*x)+(-0.35645E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99898E+00+(-0.54724E-03*(x2-d))+(-0.35645E-05*(x2-d)**2)
       f2= 0.99930E+00+(-0.22258E-03*(x2+d))+(-0.74204E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99930E+00+(-0.22258E-03*x)+(-0.74204E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99895E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99895E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10012E+01+(-0.74767E-03*x)+(-0.20976E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10012E+01+(-0.74767E-03*(x2-d))+(-0.20976E-04*(x2-d)**2)
       f2= 0.10016E+01+(-0.37410E-03*(x2+d))+(-0.19448E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10016E+01+(-0.37410E-03*x)+(-0.19448E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10005E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10005E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10041E+01+(-0.93670E-03*x)+(-0.15955E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10041E+01+(-0.93670E-03*(x2-d))+(-0.15955E-04*(x2-d)**2)
       f2= 0.10045E+01+(-0.57705E-03*(x2+d))+(-0.72145E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10045E+01+(-0.57705E-03*x)+(-0.72145E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10044E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10044E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99951E+00+(-0.35287E-03*x)+( 0.85093E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99951E+00+(-0.35287E-03*(x2-d))+( 0.85093E-05*(x2-d)**2)
       f2= 0.99946E+00+(-0.61055E-03*(x2+d))+( 0.22715E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99946E+00+(-0.61055E-03*x)+( 0.22715E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99876E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99876E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10001E+01+(-0.63180E-03*x)+(-0.61959E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10001E+01+(-0.63180E-03*(x2-d))+(-0.61959E-05*(x2-d)**2)
       f2= 0.10011E+01+(-0.33785E-03*(x2+d))+(-0.13920E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10011E+01+(-0.33785E-03*x)+(-0.13920E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10009E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10009E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99822E+00+( 0.66823E-03*x)+( 0.31288E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99822E+00+( 0.66823E-03*(x2-d))+( 0.31288E-04*(x2-d)**2)
       f2= 0.99812E+00+( 0.63055E-03*(x2+d))+(-0.16021E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99812E+00+( 0.63055E-03*x)+(-0.16021E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99920E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99920E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_2445_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_4045_2250_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4045_2250_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99985E+00+(-0.37565E-03*x)+( 0.19685E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99985E+00+(-0.37565E-03*(x2-d))+( 0.19685E-05*(x2-d)**2)
       f2= 0.99966E+00+(-0.62639E-03*(x2+d))+( 0.10193E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99966E+00+(-0.62639E-03*x)+( 0.10193E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99841E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99841E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10017E+01+(-0.73673E-03*x)+(-0.32939E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10017E+01+(-0.73673E-03*(x2-d))+(-0.32939E-04*(x2-d)**2)
       f2= 0.10025E+01+(-0.61116E-03*(x2+d))+( 0.47019E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10025E+01+(-0.61116E-03*x)+( 0.47019E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10003E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10003E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10050E+01+(-0.70193E-03*x)+(-0.74864E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10050E+01+(-0.70193E-03*(x2-d))+(-0.74864E-05*(x2-d)**2)
       f2= 0.10056E+01+(-0.88634E-04*(x2+d))+(-0.79986E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10056E+01+(-0.88634E-04*x)+(-0.79986E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10051E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10051E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99912E+00+(-0.12215E-03*x)+( 0.32007E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99912E+00+(-0.12215E-03*(x2-d))+( 0.32007E-04*(x2-d)**2)
       f2= 0.99856E+00+(-0.54568E-03*(x2+d))+( 0.25854E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99856E+00+(-0.54568E-03*x)+( 0.25854E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99771E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99771E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10005E+01+(-0.12200E-03*x)+( 0.20152E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10005E+01+(-0.12200E-03*(x2-d))+( 0.20152E-04*(x2-d)**2)
       f2= 0.99949E+00+(-0.54667E-03*(x2+d))+( 0.27299E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99949E+00+(-0.54667E-03*x)+( 0.27299E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99955E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99955E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99619E+00+( 0.35133E-03*x)+(-0.38240E-06*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99619E+00+( 0.35133E-03*(x2-d))+(-0.38240E-06*(x2-d)**2)
       f2= 0.99582E+00+( 0.54995E-03*(x2+d))+(-0.60606E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99582E+00+( 0.54995E-03*x)+(-0.60606E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99687E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99687E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4045_2250_feb98=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_4045_3375_feb98(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4045_3375_feb98
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99920E+00+(-0.51171E-03*x)+(-0.16469E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99920E+00+(-0.51171E-03*(x2-d))+(-0.16469E-04*(x2-d)**2)
       f2= 0.99867E+00+(-0.45230E-03*(x2+d))+( 0.15389E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99867E+00+(-0.45230E-03*x)+( 0.15389E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99875E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99875E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10003E+01+(-0.60089E-03*x)+(-0.12251E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10003E+01+(-0.60089E-03*(x2-d))+(-0.12251E-04*(x2-d)**2)
       f2= 0.10007E+01+(-0.45238E-03*(x2+d))+( 0.43956E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10007E+01+(-0.45238E-03*x)+( 0.43956E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99986E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99986E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10041E+01+(-0.54770E-03*x)+(-0.18224E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10041E+01+(-0.54770E-03*(x2-d))+(-0.18224E-05*(x2-d)**2)
       f2= 0.10037E+01+(-0.67450E-03*(x2+d))+(-0.78320E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10037E+01+(-0.67450E-03*x)+(-0.78320E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10032E+01+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10032E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99948E+00+(-0.50824E-03*x)+(-0.69008E-06*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99948E+00+(-0.50824E-03*(x2-d))+(-0.69008E-06*(x2-d)**2)
       f2= 0.99936E+00+(-0.65355E-03*(x2+d))+( 0.33916E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99936E+00+(-0.65355E-03*x)+( 0.33916E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99829E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99829E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99962E+00+(-0.62984E-03*x)+(-0.14367E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99962E+00+(-0.62984E-03*(x2-d))+(-0.14367E-04*(x2-d)**2)
       f2= 0.10002E+01+(-0.46012E-03*(x2+d))+( 0.14322E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10002E+01+(-0.46012E-03*x)+( 0.14322E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99985E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99985E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99664E+00+( 0.22676E-04*x)+(-0.70952E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99664E+00+( 0.22676E-04*(x2-d))+(-0.70952E-05*(x2-d)**2)
       f2= 0.99678E+00+( 0.37391E-03*(x2+d))+(-0.15429E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99678E+00+( 0.37391E-03*x)+(-0.15429E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10000E+01+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10000E+01+( 0.00000E+00*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99738E+00+( 0.00000E+00*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99738E+00+( 0.00000E+00*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4045_3375_feb98=fact*pe
      RETURN
      END

      FUNCTION correct_ppr_e1_1515_0750_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_1515_0750_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10071E+01+(-0.91977E-03*x)+(-0.34711E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10071E+01+(-0.91977E-03*(x2-d))+(-0.34711E-04*(x2-d)**2)
       f2= 0.10067E+01+(-0.91869E-03*(x2+d))+(-0.32053E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10067E+01+(-0.91869E-03*x)+(-0.32053E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10265E+01+(-0.49472E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10265E+01+(-0.49472E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97509E+00+( 0.43455E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97509E+00+( 0.43455E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10158E+01+(-0.78237E-03*x)+(-0.28723E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10158E+01+(-0.78237E-03*(x2-d))+(-0.28723E-04*(x2-d)**2)
       f2= 0.10152E+01+(-0.84656E-03*(x2+d))+(-0.30729E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10152E+01+(-0.84656E-03*x)+(-0.30729E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10377E+01+(-0.43400E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10377E+01+(-0.43400E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10782E+01+(-0.10588E-02*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10782E+01+(-0.10588E-02*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10076E+01+(-0.10251E-02*x)+(-0.90572E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10076E+01+(-0.10251E-02*(x2-d))+(-0.90572E-04*(x2-d)**2)
       f2= 0.10091E+01+( 0.19328E-03*(x2+d))+(-0.80842E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10091E+01+( 0.19328E-03*x)+(-0.80842E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.94068E+00+( 0.10037E-02*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.94068E+00+( 0.10037E-02*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.95280E+00+( 0.85014E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.95280E+00+( 0.85014E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99621E+00+(-0.11144E-02*x)+(-0.56825E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99621E+00+(-0.11144E-02*(x2-d))+(-0.56825E-04*(x2-d)**2)
       f2= 0.99588E+00+(-0.14525E-02*(x2+d))+( 0.31391E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99588E+00+(-0.14525E-02*x)+( 0.31391E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10374E+01+(-0.93595E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10374E+01+(-0.93595E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.93653E+00+( 0.90929E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.93653E+00+( 0.90929E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99777E+00+(-0.13390E-02*x)+(-0.57992E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99777E+00+(-0.13390E-02*(x2-d))+(-0.57992E-04*(x2-d)**2)
       f2= 0.99921E+00+(-0.37352E-03*(x2+d))+(-0.41054E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99921E+00+(-0.37352E-03*x)+(-0.41054E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10398E+01+(-0.95595E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10398E+01+(-0.95595E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.96452E+00+( 0.51810E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.96452E+00+( 0.51810E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99095E+00+( 0.14732E-02*x)+( 0.22084E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99095E+00+( 0.14732E-02*(x2-d))+( 0.22084E-04*(x2-d)**2)
       f2= 0.98970E+00+( 0.18868E-02*(x2+d))+(-0.47334E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98970E+00+( 0.18868E-02*x)+(-0.47334E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10146E+01+(-0.58024E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10146E+01+(-0.58024E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.93455E+00+( 0.85101E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.93455E+00+( 0.85101E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_1515_0750_jan99=fact*pe
      RETURN
      END 
 
      FUNCTION correct_ppr_e1_1515_1500_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_1515_1500_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10020E+01+(-0.85600E-03*x)+(-0.32238E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10020E+01+(-0.85600E-03*(x2-d))+(-0.32238E-04*(x2-d)**2)
       f2= 0.10025E+01+(-0.56358E-03*(x2+d))+(-0.10773E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10025E+01+(-0.56358E-03*x)+(-0.10773E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10219E+01+(-0.45913E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10219E+01+(-0.45913E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98230E+00+( 0.30528E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98230E+00+( 0.30528E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10086E+01+(-0.51630E-03*x)+(-0.25460E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10086E+01+(-0.51630E-03*(x2-d))+(-0.25460E-04*(x2-d)**2)
       f2= 0.10080E+01+(-0.74922E-03*(x2+d))+(-0.10673E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10080E+01+(-0.74922E-03*x)+(-0.10673E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10078E+01+(-0.68392E-04*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10078E+01+(-0.68392E-04*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98616E+00+( 0.31111E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98616E+00+( 0.31111E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10043E+01+(-0.18298E-03*x)+(-0.12471E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10043E+01+(-0.18298E-03*(x2-d))+(-0.12471E-04*(x2-d)**2)
       f2= 0.10038E+01+(-0.25328E-03*(x2+d))+(-0.28567E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10038E+01+(-0.25328E-03*x)+(-0.28567E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97575E+00+( 0.42672E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97575E+00+( 0.42672E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.96666E+00+( 0.59000E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.96666E+00+( 0.59000E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99745E+00+(-0.57934E-03*x)+(-0.12859E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99745E+00+(-0.57934E-03*(x2-d))+(-0.12859E-04*(x2-d)**2)
       f2= 0.99725E+00+(-0.70535E-03*(x2+d))+(-0.61957E-06*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99725E+00+(-0.70535E-03*x)+(-0.61957E-06*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10129E+01+(-0.38514E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10129E+01+(-0.38514E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97090E+00+( 0.42228E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97090E+00+( 0.42228E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10007E+01+(-0.23177E-03*x)+( 0.23148E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10007E+01+(-0.23177E-03*(x2-d))+( 0.23148E-05*(x2-d)**2)
       f2= 0.10001E+01+(-0.74320E-03*(x2+d))+( 0.83453E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10001E+01+(-0.74320E-03*x)+( 0.83453E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10202E+01+(-0.51121E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10202E+01+(-0.51121E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97088E+00+( 0.44695E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97088E+00+( 0.44695E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99327E+00+( 0.77412E-03*x)+( 0.31059E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99327E+00+( 0.77412E-03*(x2-d))+( 0.31059E-04*(x2-d)**2)
       f2= 0.99385E+00+( 0.78701E-03*(x2+d))+(-0.23072E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99385E+00+( 0.78701E-03*x)+(-0.23072E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10097E+01+(-0.34254E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10097E+01+(-0.34254E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97673E+00+( 0.26247E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97673E+00+( 0.26247E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_1515_1500_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_ppr_e1_1515_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_1515_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10004E+01+(-0.63001E-03*x)+(-0.19626E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10004E+01+(-0.63001E-03*(x2-d))+(-0.19626E-04*(x2-d)**2)
       f2= 0.10004E+01+(-0.64032E-03*(x2+d))+( 0.26398E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10004E+01+(-0.64032E-03*x)+( 0.26398E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10052E+01+(-0.14226E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10052E+01+(-0.14226E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97846E+00+( 0.35257E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97846E+00+( 0.35257E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10035E+01+(-0.67336E-03*x)+(-0.29094E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10035E+01+(-0.67336E-03*(x2-d))+(-0.29094E-04*(x2-d)**2)
       f2= 0.10039E+01+(-0.50735E-03*(x2+d))+(-0.16689E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10039E+01+(-0.50735E-03*x)+(-0.16689E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10123E+01+(-0.24041E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10123E+01+(-0.24041E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97813E+00+( 0.40351E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97813E+00+( 0.40351E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10013E+01+(-0.15712E-03*x)+(-0.53311E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10013E+01+(-0.15712E-03*(x2-d))+(-0.53311E-05*(x2-d)**2)
       f2= 0.10011E+01+(-0.37377E-03*(x2+d))+(-0.17077E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10011E+01+(-0.37377E-03*x)+(-0.17077E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98159E+00+( 0.28856E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98159E+00+( 0.28856E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.96666E+00+( 0.57189E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.96666E+00+( 0.57189E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99610E+00+(-0.57453E-03*x)+(-0.78370E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99610E+00+(-0.57453E-03*(x2-d))+(-0.78370E-05*(x2-d)**2)
       f2= 0.99596E+00+(-0.53392E-03*(x2+d))+(-0.50218E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99596E+00+(-0.53392E-03*x)+(-0.50218E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10089E+01+(-0.30122E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10089E+01+(-0.30122E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.96992E+00+( 0.44149E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.96992E+00+( 0.44149E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99762E+00+(-0.62006E-03*x)+(-0.16995E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99762E+00+(-0.62006E-03*(x2-d))+(-0.16995E-04*(x2-d)**2)
       f2= 0.99809E+00+(-0.44634E-03*(x2+d))+(-0.87844E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99809E+00+(-0.44634E-03*x)+(-0.87844E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10076E+01+(-0.25993E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10076E+01+(-0.25993E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97188E+00+( 0.42836E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97188E+00+( 0.42836E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99287E+00+( 0.18668E-03*x)+( 0.64166E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99287E+00+( 0.18668E-03*(x2-d))+( 0.64166E-05*(x2-d)**2)
       f2= 0.99298E+00+( 0.26410E-03*(x2+d))+( 0.75348E-06*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99298E+00+( 0.26410E-03*x)+( 0.75348E-06*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10027E+01+(-0.20089E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10027E+01+(-0.20089E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97723E+00+( 0.27146E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97723E+00+( 0.27146E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_1515_2250_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_ppr_e1_2567_1500_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_2567_1500_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector      
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10015E+01+(-0.29568E-03*x)+( 0.25634E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10015E+01+(-0.29568E-03*(x2-d))+( 0.25634E-04*(x2-d)**2)
       f2= 0.10009E+01+(-0.84567E-03*(x2+d))+( 0.48252E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10009E+01+(-0.84567E-03*x)+( 0.48252E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10106E+01+(-0.19703E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10106E+01+(-0.19703E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10021E+01+(-0.53898E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10021E+01+(-0.53898E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10068E+01+(-0.38922E-03*x)+( 0.10207E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10068E+01+(-0.38922E-03*(x2-d))+( 0.10207E-04*(x2-d)**2)
       f2= 0.10064E+01+(-0.81771E-03*(x2+d))+( 0.29547E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10064E+01+(-0.81771E-03*x)+( 0.29547E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10118E+01+(-0.10721E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10118E+01+(-0.10721E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10065E+01+(-0.35448E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10065E+01+(-0.35448E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10023E+01+(-0.34560E-03*x)+(-0.28846E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10023E+01+(-0.34560E-03*(x2-d))+(-0.28846E-04*(x2-d)**2)
       f2= 0.10026E+01+(-0.12978E-03*(x2+d))+(-0.60832E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10026E+01+(-0.12978E-03*x)+(-0.60832E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97567E+00+( 0.43163E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97567E+00+( 0.43163E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97306E+00+( 0.48353E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97306E+00+( 0.48353E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99436E+00+(-0.51542E-03*x)+( 0.41358E-07*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99436E+00+(-0.51542E-03*(x2-d))+( 0.41358E-07*(x2-d)**2)
       f2= 0.99393E+00+(-0.86086E-03*(x2+d))+( 0.42001E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99393E+00+(-0.86086E-03*x)+( 0.42001E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10260E+01+(-0.57509E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10260E+01+(-0.57509E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10557E+01+(-0.10378E-02*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10557E+01+(-0.10378E-02*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99628E+00+(-0.49230E-03*x)+(-0.38592E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99628E+00+(-0.49230E-03*(x2-d))+(-0.38592E-05*(x2-d)**2)
       f2= 0.99587E+00+(-0.67190E-03*(x2+d))+( 0.14124E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99587E+00+(-0.67190E-03*x)+( 0.14124E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10165E+01+(-0.44662E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10165E+01+(-0.44662E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97897E+00+( 0.27471E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97897E+00+( 0.27471E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98832E+00+( 0.70472E-03*x)+( 0.27671E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98832E+00+( 0.70472E-03*(x2-d))+( 0.27671E-04*(x2-d)**2)
       f2= 0.98861E+00+( 0.94569E-03*(x2+d))+(-0.81589E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98861E+00+( 0.94569E-03*x)+(-0.81589E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10192E+01+(-0.56170E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10192E+01+(-0.56170E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97609E+00+( 0.21343E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97609E+00+( 0.21343E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_2567_1500_jan99=fact*pe
      RETURN
      END 

      FUNCTION correct_ppr_e1_2567_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_2567_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector      
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10004E+01+(-0.38534E-03*x)+( 0.70436E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10004E+01+(-0.38534E-03*(x2-d))+( 0.70436E-05*(x2-d)**2)
       f2= 0.99989E+00+(-0.71487E-03*(x2+d))+( 0.30926E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99989E+00+(-0.71487E-03*x)+( 0.30926E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10112E+01+(-0.24751E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10112E+01+(-0.24751E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99020E+00+( 0.13827E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99020E+00+( 0.13827E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10018E+01+(-0.89395E-03*x)+(-0.44423E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10018E+01+(-0.89395E-03*(x2-d))+(-0.44423E-04*(x2-d)**2)
       f2= 0.10030E+01+(-0.25232E-03*(x2+d))+(-0.17524E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10030E+01+(-0.25232E-03*x)+(-0.17524E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10065E+01+(-0.87728E-04*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10065E+01+(-0.87728E-04*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99500E+00+( 0.12433E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99500E+00+( 0.12433E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10014E+01+(-0.62111E-04*x)+(-0.58678E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10014E+01+(-0.62111E-04*(x2-d))+(-0.58678E-05*(x2-d)**2)
       f2= 0.10012E+01+(-0.34143E-03*(x2+d))+(-0.38808E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10012E+01+(-0.34143E-03*x)+(-0.38808E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98222E+00+( 0.30581E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98222E+00+( 0.30581E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97275E+00+( 0.48144E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97275E+00+( 0.48144E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99341E+00+(-0.81662E-03*x)+(-0.21548E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99341E+00+(-0.81662E-03*(x2-d))+(-0.21548E-04*(x2-d)**2)
       f2= 0.99367E+00+(-0.59363E-03*(x2+d))+( 0.16998E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99367E+00+(-0.59363E-03*x)+( 0.16998E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10145E+01+(-0.39818E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10145E+01+(-0.39818E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98427E+00+( 0.17561E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98427E+00+( 0.17561E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99594E+00+(-0.56339E-03*x)+(-0.11061E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99594E+00+(-0.56339E-03*(x2-d))+(-0.11061E-04*(x2-d)**2)
       f2= 0.99594E+00+(-0.59264E-03*(x2+d))+( 0.10971E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99594E+00+(-0.59264E-03*x)+( 0.10971E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10128E+01+(-0.35865E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10128E+01+(-0.35865E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98230E+00+( 0.23178E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98230E+00+( 0.23178E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99084E+00+( 0.56935E-03*x)+( 0.35896E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99084E+00+( 0.56935E-03*(x2-d))+( 0.35896E-04*(x2-d)**2)
       f2= 0.99051E+00+( 0.30840E-03*(x2+d))+( 0.11197E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99051E+00+( 0.30840E-03*x)+( 0.11197E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10160E+01+(-0.47896E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10160E+01+(-0.47896E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98659E+00+( 0.69996E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98659E+00+( 0.69996E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_2567_2250_jan99=fact*pe
      RETURN
      END 


      FUNCTION correct_ppr_e1_4056_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact, correct_ppr_e1_4056_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10003E+01+(-0.53545E-03*x)+( 0.11789E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10003E+01+(-0.53545E-03*(x2-d))+( 0.11789E-04*(x2-d)**2)
       f2= 0.99999E+00+(-0.74190E-03*(x2+d))+( 0.13842E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99999E+00+(-0.74190E-03*x)+( 0.13842E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10080E+01+(-0.18263E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10080E+01+(-0.18263E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10080E+01+(-0.18263E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10080E+01+(-0.18263E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10052E+01+(-0.39692E-03*x)+(-0.24286E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10052E+01+(-0.39692E-03*(x2-d))+(-0.24286E-05*(x2-d)**2)
       f2= 0.10049E+01+(-0.69323E-03*(x2+d))+( 0.92179E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10049E+01+(-0.69323E-03*x)+( 0.92179E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10101E+01+(-0.13650E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10101E+01+(-0.13650E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10101E+01+(-0.13650E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10101E+01+(-0.13650E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99925E+00+(-0.43236E-03*x)+(-0.20951E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99925E+00+(-0.43236E-03*(x2-d))+(-0.20951E-04*(x2-d)**2)
       f2= 0.99899E+00+(-0.64233E-03*(x2+d))+(-0.19935E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99899E+00+(-0.64233E-03*x)+(-0.19935E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97697E+00+( 0.42836E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97697E+00+( 0.42836E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97697E+00+( 0.42836E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97697E+00+( 0.42836E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99263E+00+(-0.35771E-03*x)+( 0.13348E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99263E+00+(-0.35771E-03*(x2-d))+( 0.13348E-04*(x2-d)**2)
       f2= 0.99271E+00+(-0.51406E-03*(x2+d))+( 0.47934E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99271E+00+(-0.51406E-03*x)+( 0.47934E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10114E+01+(-0.32236E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10114E+01+(-0.32236E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10114E+01+(-0.32236E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10114E+01+(-0.32236E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99266E+00+(-0.93134E-03*x)+(-0.18609E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99266E+00+(-0.93134E-03*(x2-d))+(-0.18609E-04*(x2-d)**2)
       f2= 0.99355E+00+(-0.32700E-03*(x2+d))+( 0.11455E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99355E+00+(-0.32700E-03*x)+( 0.11455E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10047E+01+(-0.18707E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10047E+01+(-0.18707E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10047E+01+(-0.18707E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10047E+01+(-0.18707E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98914E+00+( 0.76081E-03*x)+( 0.45110E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98914E+00+( 0.76081E-03*(x2-d))+( 0.45110E-04*(x2-d)**2)
       f2= 0.98910E+00+( 0.59329E-03*(x2+d))+( 0.14798E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98910E+00+( 0.59329E-03*x)+( 0.14798E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10149E+01+(-0.47131E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10149E+01+(-0.47131E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10149E+01+(-0.47131E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10149E+01+(-0.47131E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4056_2250_jan99=fact*pe
      RETURN
      END 



      FUNCTION correct_ppr_e1_4247_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4247_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99745E+00+(-0.16632E-02*x)+(-0.79275E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99745E+00+(-0.16632E-02*(x2-d))+(-0.79275E-04*(x2-d)**2)
       f2= 0.99909E+00+(-0.74546E-03*(x2+d))+( 0.21382E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99909E+00+(-0.74546E-03*x)+( 0.21382E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10095E+01+(-0.22710E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10095E+01+(-0.22710E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99963E+00+(-0.38816E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99963E+00+(-0.38816E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10058E+01+(-0.47478E-03*x)+(-0.10096E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10058E+01+(-0.47478E-03*(x2-d))+(-0.10096E-04*(x2-d)**2)
       f2= 0.10050E+01+(-0.86192E-03*(x2+d))+( 0.17241E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10050E+01+(-0.86192E-03*x)+( 0.17241E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10191E+01+(-0.31618E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10191E+01+(-0.31618E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99877E+00+( 0.66226E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99877E+00+( 0.66226E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99878E+00+(-0.35131E-03*x)+(-0.92809E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99878E+00+(-0.35131E-03*(x2-d))+(-0.92809E-05*(x2-d)**2)
       f2= 0.99857E+00+(-0.44009E-03*(x2+d))+(-0.28523E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99857E+00+(-0.44009E-03*x)+(-0.28523E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97565E+00+( 0.46627E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97565E+00+( 0.46627E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98272E+00+( 0.32463E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98272E+00+( 0.32463E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99237E+00+(-0.49122E-03*x)+( 0.26292E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99237E+00+(-0.49122E-03*(x2-d))+( 0.26292E-05*(x2-d)**2)
       f2= 0.99202E+00+(-0.40830E-03*(x2+d))+( 0.47180E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99202E+00+(-0.40830E-03*x)+( 0.47180E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10160E+01+(-0.43557E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10160E+01+(-0.43557E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99976E+00+(-0.11050E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99976E+00+(-0.11050E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99543E+00+(-0.15201E-03*x)+( 0.15176E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99543E+00+(-0.15201E-03*(x2-d))+( 0.15176E-04*(x2-d)**2)
       f2= 0.99436E+00+(-0.49682E-03*(x2+d))+( 0.12917E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99436E+00+(-0.49682E-03*x)+( 0.12917E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10123E+01+(-0.35778E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10123E+01+(-0.35778E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98019E+00+( 0.25902E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98019E+00+( 0.25902E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2= 54.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98690E+00+( 0.79444E-04*x)+( 0.10798E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98690E+00+( 0.79444E-04*(x2-d))+( 0.10798E-04*(x2-d)**2)
       f2= 0.98831E+00+( 0.45155E-03*(x2+d))+( 0.25890E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98831E+00+( 0.45155E-03*x)+( 0.25890E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10055E+01+(-0.26745E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10055E+01+(-0.26745E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10288E+01+(-0.72415E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10288E+01+(-0.72415E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4247_2250_jan99=fact*pe
      RETURN
      END 


      FUNCTION correct_ppr_e1_4462_2250_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4462_2250_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10009E+01+(-0.81730E-04*x)+( 0.48084E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10009E+01+(-0.81730E-04*(x2-d))+( 0.48084E-04*(x2-d)**2)
       f2= 0.10006E+01+(-0.89966E-03*(x2+d))+( 0.16264E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10006E+01+(-0.89966E-03*x)+( 0.16264E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10100E+01+(-0.23829E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10100E+01+(-0.23829E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10100E+01+(-0.23829E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10100E+01+(-0.23829E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10031E+01+(-0.13515E-02*x)+(-0.55206E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10031E+01+(-0.13515E-02*(x2-d))+(-0.55206E-04*(x2-d)**2)
       f2= 0.10045E+01+(-0.65230E-03*(x2+d))+( 0.12098E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10045E+01+(-0.65230E-03*x)+( 0.12098E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10154E+01+(-0.24249E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10154E+01+(-0.24249E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10154E+01+(-0.24249E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10154E+01+(-0.24249E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99866E+00+(-0.54801E-03*x)+(-0.30142E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99866E+00+(-0.54801E-03*(x2-d))+(-0.30142E-04*(x2-d)**2)
       f2= 0.99809E+00+(-0.67057E-03*(x2+d))+(-0.36999E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99809E+00+(-0.67057E-03*x)+(-0.36999E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.97349E+00+( 0.50155E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.97349E+00+( 0.50155E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.97349E+00+( 0.50155E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.97349E+00+( 0.50155E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99257E+00+(-0.39439E-03*x)+( 0.34358E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99257E+00+(-0.39439E-03*(x2-d))+( 0.34358E-05*(x2-d)**2)
       f2= 0.99248E+00+(-0.35933E-03*(x2+d))+( 0.41864E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99248E+00+(-0.35933E-03*x)+( 0.41864E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10086E+01+(-0.27715E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10086E+01+(-0.27715E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10086E+01+(-0.27715E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10086E+01+(-0.27715E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99339E+00+(-0.60273E-03*x)+(-0.29190E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99339E+00+(-0.60273E-03*(x2-d))+(-0.29190E-05*(x2-d)**2)
       f2= 0.99356E+00+(-0.50834E-03*(x2+d))+( 0.23745E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99356E+00+(-0.50834E-03*x)+( 0.23745E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10075E+01+(-0.24854E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10075E+01+(-0.24854E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10075E+01+(-0.24854E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10075E+01+(-0.24854E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.98864E+00+( 0.34383E-03*x)+( 0.10138E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.98864E+00+( 0.34383E-03*(x2-d))+( 0.10138E-04*(x2-d)**2)
       f2= 0.98876E+00+( 0.25160E-03*(x2+d))+( 0.34321E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98876E+00+( 0.25160E-03*x)+( 0.34321E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10127E+01+(-0.44222E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10127E+01+(-0.44222E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10127E+01+(-0.44222E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10127E+01+(-0.44222E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4462_2250_jan99=fact*pe
      RETURN
      END 


      FUNCTION correct_ppr_e1_4462_3375_jan99(pe,theta,phi) 
      IMPLICIT none
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4462_3375_jan99
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99924E+00+(-0.56387E-04*x)+( 0.29102E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99924E+00+(-0.56387E-04*(x2-d))+( 0.29102E-04*(x2-d)**2)
       f2= 0.99842E+00+(-0.58877E-03*(x2+d))+( 0.14388E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99842E+00+(-0.58877E-03*x)+( 0.14388E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99991E+00+(-0.47199E-04*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99991E+00+(-0.47199E-04*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99991E+00+(-0.47199E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99991E+00+(-0.47199E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10010E+01+(-0.58126E-03*x)+(-0.17108E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10010E+01+(-0.58126E-03*(x2-d))+(-0.17108E-04*(x2-d)**2)
       f2= 0.10013E+01+(-0.55295E-03*(x2+d))+( 0.97972E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10013E+01+(-0.55295E-03*x)+( 0.97972E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10008E+01+( 0.84212E-06*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10008E+01+( 0.84212E-06*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10008E+01+( 0.84212E-06*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10008E+01+( 0.84212E-06*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99753E+00+(-0.17310E-03*x)+( 0.32649E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99753E+00+(-0.17310E-03*(x2-d))+( 0.32649E-05*(x2-d)**2)
       f2= 0.99713E+00+(-0.47565E-03*(x2+d))+(-0.16102E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99713E+00+(-0.47565E-03*x)+(-0.16102E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98670E+00+( 0.22758E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98670E+00+( 0.22758E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98670E+00+( 0.22758E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98670E+00+( 0.22758E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99304E+00+(-0.40527E-03*x)+( 0.15065E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99304E+00+(-0.40527E-03*(x2-d))+( 0.15065E-04*(x2-d)**2)
       f2= 0.99353E+00+(-0.24074E-03*(x2+d))+( 0.11171E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99353E+00+(-0.24074E-03*x)+( 0.11171E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10112E+01+(-0.33842E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10112E+01+(-0.33842E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10112E+01+(-0.33842E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10112E+01+(-0.33842E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99561E+00+(-0.15215E-03*x)+( 0.15499E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99561E+00+(-0.15215E-03*(x2-d))+( 0.15499E-04*(x2-d)**2)
       f2= 0.99479E+00+(-0.64493E-03*(x2+d))+( 0.31950E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99479E+00+(-0.64493E-03*x)+( 0.31950E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10065E+01+(-0.24120E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10065E+01+(-0.24120E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10065E+01+(-0.24120E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10065E+01+(-0.24120E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99071E+00+( 0.45330E-03*x)+( 0.39100E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99071E+00+( 0.45330E-03*(x2-d))+( 0.39100E-04*(x2-d)**2)
       f2= 0.99041E+00+( 0.32956E-03*(x2+d))+( 0.10891E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99041E+00+( 0.32956E-03*x)+( 0.10891E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10101E+01+(-0.38338E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10101E+01+(-0.38338E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10101E+01+(-0.38338E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10101E+01+(-0.38338E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4462_3375_jan99=fact*pe
      RETURN
      END 
      FUNCTION correct_ppr_e1_4817_3375_feb00(pe,theta,phi) 
      REAL pe,theta,phi, phisec,fact,correct_ppr_e1_4817_3375_feb00
      REAL d,x,y,x2,y2,f1,f2
      INTEGER nsector
      phisec=(mod(phi+30,60.)-30.) 
      fact=1.0 
      nsector=0 
      IF(phi.lt. 30..or. phi.gt.330.) nsector=1 
      IF(phi.gt. 30..and.phi.lt. 90.) nsector=2 
      IF(phi.gt. 90..and.phi.lt.150.) nsector=3 
      IF(phi.gt.150..and.phi.lt.210.) nsector=4 
      IF(phi.gt.210..and.phi.lt.270.) nsector=5 
      IF(phi.gt.270..and.phi.lt.330.) nsector=6 
       
      IF(nsector.eq. 1)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10035E+01+(-0.72043E-03*x)+(-0.24746E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10035E+01+(-0.72043E-03*(x2-d))+(-0.24746E-04*(x2-d)**2)
       f2= 0.10039E+01+(-0.81753E-03*(x2+d))+( 0.25940E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10039E+01+(-0.81753E-03*x)+( 0.25940E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10051E+01+(-0.39920E-04*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10051E+01+(-0.39920E-04*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10051E+01+(-0.39920E-04*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10051E+01+(-0.39920E-04*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 2)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10052E+01+(-0.32923E-03*x)+(-0.87970E-05*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10052E+01+(-0.32923E-03*(x2-d))+(-0.87970E-05*(x2-d)**2)
       f2= 0.10054E+01+(-0.37901E-03*(x2+d))+(-0.62077E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10054E+01+(-0.37901E-03*x)+(-0.62077E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.99645E+00+( 0.16794E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.99645E+00+( 0.16794E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.99645E+00+( 0.16794E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99645E+00+( 0.16794E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 3)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10027E+01+( 0.40059E-04*x)+( 0.10698E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10027E+01+( 0.40059E-04*(x2-d))+( 0.10698E-04*(x2-d)**2)
       f2= 0.10023E+01+(-0.39324E-03*(x2+d))+(-0.33945E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10023E+01+(-0.39324E-03*x)+(-0.33945E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.98161E+00+( 0.43803E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.98161E+00+( 0.43803E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.98161E+00+( 0.43803E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.98161E+00+( 0.43803E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 4)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99800E+00+(-0.11939E-03*x)+( 0.34331E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99800E+00+(-0.11939E-03*(x2-d))+( 0.34331E-04*(x2-d)**2)
       f2= 0.99928E+00+(-0.32606E-03*(x2+d))+( 0.12722E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99928E+00+(-0.32606E-03*x)+( 0.12722E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10053E+01+(-0.11807E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10053E+01+(-0.11807E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10053E+01+(-0.11807E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10053E+01+(-0.11807E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 5)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.10005E+01+(-0.56263E-03*x)+(-0.12309E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.10005E+01+(-0.56263E-03*(x2-d))+(-0.12309E-04*(x2-d)**2)
       f2= 0.10004E+01+(-0.24561E-03*(x2+d))+(-0.15636E-04*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10004E+01+(-0.24561E-03*x)+(-0.15636E-04*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10104E+01+(-0.20090E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10104E+01+(-0.20090E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10104E+01+(-0.20090E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10104E+01+(-0.20090E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
      IF(nsector.eq. 6)THEN
       x=phisec
       y=theta
       x2=  0.000
       y2=-10.000
       d=2.0
       IF    (phisec.le.x2-d)THEN
       fact=fact*( 0.99609E+00+( 0.76337E-03*x)+( 0.58227E-04*x**2))
       ELSEIF(phisec.ge.x2-d.and.phisec.le.x2+d)THEN
       f1= 0.99609E+00+( 0.76337E-03*(x2-d))+( 0.58227E-04*(x2-d)**2)
       f2= 0.99539E+00+( 0.42567E-03*(x2+d))+( 0.50114E-05*(x2+d)**2)
       fact=fact*((f2-f1)*x + ((x2+d)*f1-(x2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.99539E+00+( 0.42567E-03*x)+( 0.50114E-05*x**2))
       ENDIF
       IF    (theta.le.y2-d)THEN
       fact=fact*( 0.10055E+01+(-0.18614E-03*y)+( 0.00000E+00*y**2))
       ELSEIF(theta.ge.y2-d.and.theta.le.y2+d)THEN
       f1= 0.10055E+01+(-0.18614E-03*(y2-d))+( 0.00000E+00*(y2-d)**2)
       f2= 0.10055E+01+(-0.18614E-03*(y2+d))+( 0.00000E+00*(y2+d)**2)
       fact=fact*((f2-f1)*y + ((y2+d)*f1-(y2-d)*f2))/2./d
       ELSE
       fact=fact*( 0.10055E+01+(-0.18614E-03*y)+( 0.00000E+00*y**2))
       ENDIF
      ENDIF
       
      correct_ppr_e1_4817_3375_feb00=fact*pe
      RETURN
      END 






c
c     Burkert's momentum corrections, got 01feb1999 
c     ch1,ebeam,torcur,sector,theta,phi,q -- input
c     pmom -- input uncorrected momentum / output corrected momentum
c
      SUBROUTINE momcorr_e1_bu(ch1,ebeam,torcur,sector,pmom,theta,phi,q)
      IMPLICIT none
      CHARACTER*1 ch1
      INTEGER     sector,secte,secth
      REAL        pi,pi180,ebeam,torcur,pmom,theta,phi,q,q_h
      REAL        cx,cy,cz
      REAL        pel,phiel,phield,thetael,thetaeld
      REAL        p_h,thetah,thetahd
      INTEGER     ifirste,ifirsth
      DATA        ifirste/0/,ifirsth/0/


      pi=acos(-1.)
      pi180=acos(-1.)/180.


      IF(ch1.eq.'E') THEN
      IF(ifirste.le.1) ifirste=ifirste+1
        
        secte=sector
        pel=pmom
        cx=sin(theta*pi180)*cos(phi*pi180)
        cy=sin(theta*pi180)*sin(phi*pi180)
        cz=cos(theta*pi180)

        thetael=acos(cz)
        thetaeld=thetael*90./acos(0.)
        phiel=atan2(cy,cx)
        phield=atan2(cy,cx)*90./acos(0.)



*******************************************************************
* momentum  corrections electron start here
******************************************************************** 
** corrections for e1b data
** 1.515 GeV
*        
         if(ebeam.ge.1.5.and.ebeam.lt.1.6)then
*
** torus current 750 A
*
         if(torcur.gt.700..and.torcur.lt.800.)then  

         if(ifirste.eq.1) then
         print *,' Corrections for e1b 1.515 GeV, torus 750 A'
         endif

          if(secte.eq.1)then 
         if(cos(thetael).gt.0.95)then
           pel = pel*(1.- (cos(thetael)-0.95)*0.6)
           pel=pel*(1. - phiel*0.07 + phiel**2*0.0)
         endif 
             pel = pel*(1.-0.017)
          endif
           
           if(secte.eq.2)then
         if(cos(thetael).gt.0.96)then
           pel = pel*(1.- (cos(thetael)-0.96)*0.7)         
         endif 
          pel=pel*(1. - (phiel-pi/3.)*0.055 + (phiel-pi/3.)**2*0.05)
           pel = pel*(1.-0.017)
           endif
               
           if(secte.eq.3)then
         if(cos(thetael).gt.0.95)then
           pel = pel*(1.- (cos(thetael)-0.95)*0.60)
         endif 
             pel = pel*(1.+0.015)
          endif
         
           if(secte.eq.4)then
         if(cos(thetael).gt.0.95)then
           pel = pel*(1.- (cos(thetael)-0.95)*0.78)
         endif 
             pel = pel*(1.-0.01)
            endif
            
           if(secte.eq.5)then
         if(cos(thetael).gt.0.95)then
           pel = pel*(1.- (cos(thetael)-0.95)*0.73)
         endif 
             pel = pel*(1.-0.013) 
           endif
           
           if(secte.eq.6)then
         if(cos(thetael).gt.0.95)then
           pel = pel*(1.- (cos(thetael)-0.95)*0.68)
         endif 
             pel = pel*(1.+0.008) 
            endif  
 
          endif

*
** torus current 1500 A
*

         if(torcur.gt.1400..and.torcur.lt.1600.)then  

         if(ifirste.eq.1) then
         print *,' Corrections for e1b 1.515 GeV, torus 1500 A'
         endif

          if(secte.eq.1)then 
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.15)
           pel=pel*(1. - phiel*0.06 + phiel**2*0.06)
         endif 
             pel = pel*(1.-0.011)
          endif
           
           if(secte.eq.2)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.2)         
         endif 
          pel=pel*(1. - (phiel-pi/3.)*0.035 + (phiel-pi/3.)**2*0.04)
           pel = pel*(1.-0.011)
           endif
               
           if(secte.eq.3)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.33)
         endif 
             pel = pel*(1.+0.007)
          endif
         
           if(secte.eq.4)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.35)
         endif 
             pel = pel*(1.-0.0055)
            endif
            
           if(secte.eq.5)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.4)
         endif 
             pel = pel*(1.-0.005) 
           endif
           
           if(secte.eq.6)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.3)
         endif 
             pel = pel*(1.+0.0025) 
            endif  
 
          endif

*
** torus current 2250 A
*
         if(torcur.gt.2200..and.torcur.lt.2300.)then  

         if(ifirste.eq.1) then
         print *,' Corrections for e1b 1.515 GeV, torus 2250 A'
         endif

          if(secte.eq.1)then 
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.15)
           pel=pel*(1. - phiel*0.07 + phiel**2*0.06)
         endif 
             pel = pel*(1.-0.0105)
          endif
           
           if(secte.eq.2)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.2)         
         endif 
          pel=pel*(1. - (phiel-pi/3.)*0.025 + (phiel-pi/3.)**2*0.03)
           pel = pel*(1.-0.011)
           endif
*               
           if(secte.eq.3)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.30)
         endif 
             pel = pel*(1.+0.0014)
          endif

           if(secte.eq.4)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.35)
         endif 
             pel = pel*(1.-0.0065)
            endif
            
           if(secte.eq.5)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.35)
         endif 
             pel = pel*(1.-0.0057) 
           endif
           
           if(secte.eq.6)then
         if(cos(thetael).gt.0.93)then
           pel = pel*(1.- (cos(thetael)-0.93)*0.25)
         endif 
             pel = pel*(1.+0.0003) 
            endif  
 
          endif

         endif    



* corrections for e1a data
         if(ebeam.ge.1.6.and.ebeam.le.1.7)then
* 60% field run         
          if(torcur.gt.2200.and.torcur.lt.2300)then         

         if(ifirste.eq.1) then
         print *,' Corrections for e1a 1.6-1.7 GeV, torus 2250 A'
         endif

           if(secte.eq.1)then 
         pel=pel*(1. - phiel*0.02 + phiel**2*0.0)
             pel = pel*(1.-0.010)
          endif
           
           if(secte.eq.2)then
              pel = pel*(1. - 0.006)
           endif
              
           if(secte.eq.3)then
           pel = pel*(1.- 0.0005)
       pel = pel*(1.+(phiel-2.*pi/3.)*0.01+(phiel-2.*pi/3.)**2*0.02)         
           endif

           if(secte.eq.4)then
            pel = pel*(1. - 0.0055)
            endif

           if(secte.eq.5)then
           pel = pel*(1.-0.0050)
       pel = pel*(1.-(phiel+2.*pi/3.)*0.01-(phiel+2.*pi/3.)**2*0.03)
           endif

           if(secte.eq.6)then
            pel = pel*(1. + 0.0020)
           endif  
          endif
          
* 40% field run
         if(ifirste.eq.1) then
         print *,' Corrections for e1a 1.6-1.7 GeV, torus AGAIN '
         endif

           if(secte.eq.1)then 
            pel=pel*(1. - phiel*0.04)
            pel = pel*(1.+0.0160)
            pel = pel*(1.-0.1*(1.-cos((thetaeld-30.)*pi/180.)**2))            
            endif
           
           if(secte.eq.2)then
              pel = pel*(1.+0.016)
              pel = pel*(1. - (phiel-pi/3.)*0.01)
           endif
             
           if(secte.eq.3)then
           pel = pel*(1.+0.021)
           pel = pel*(1. + (phiel-2.*pi/3.)*0.015)
           endif

           if(secte.eq.4)then
            pel = pel*(1.+0.021)
            pel=pel*(1.-0.1*(1.-cos((thetaeld-30.)*pi/180.)**2))
            endif

           if(secte.eq.5)then
           pel = pel*(1.+0.0155)
           pel=pel*(1.-0.1*(1.-cos((thetaeld-30.)*pi/180.)**2))
           endif

           if(secte.eq.6)then
            pel = pel*(1.+0.02)
            pel=pel*(1.-0.05*(1.-cos((thetaeld-30.)*pi/180.)**2))
           endif  
 
          endif
          
          if(ebeam.ge.2.2.and.ebeam.le.2.5)then
* 60% field run             

         if(ifirste.eq.1) then
         print *,' Corrections for e1a 2.2-2.5 GeV, torus 60% '
         endif

           if(secte.eq.1)then 
            if(thetaeld.lt.23.)pel=pel*(1. - phiel*0.02 + phiel**2*0.04)
            pel = pel*(1.- 0.0085)
            if(thetaeld.le.18.)pel=pel*(1.-(18.-thetaeld)*0.0024)
           if(thetaeld.ge.20.)pel = pel*(1. + (20.-thetaeld)*0.0008)
            endif

           if(secte.eq.2)then
           pel = pel*(1.-0.0025)
           if(thetaeld.lt.23.)then
           pel = pel*(1. -(phiel-pi/3.)*0.015 + (phiel-pi/3.)**2*0.026)
           endif
           if(thetaeld.lt.18)pel=pel*(1. - (18.-thetaeld)*0.0020)
           if(thetaeld.ge.20.)pel = pel*(1. + (20.-thetaeld)*0.0006)           
           endif

           if(secte.eq.3)then
           pel = pel*(1.-0.002)
           pel = pel*(1.+(phiel-2.*pi/3.)*0.015)
          if(thetaeld.lt.23.)then
          pel = pel*(1.+(phiel-2.*pi/3.)**2*0.05)
          endif
          if(thetaeld.lt.18)pel=pel*(1. - (18. - thetaeld)*0.0010)
           endif

           if(secte.eq.4)then
            pel = pel*(1.-0.007)
            if(phiel.lt.pi.and.phiel.gt.0.)then
          if(thetaeld.lt.23.)then
          pel = pel*(1. + (phiel-pi)*0.0+(phiel - pi)**2*0.04)
          endif
            endif

            if(phiel.gt.-pi.and.phiel.lt.0.)then
          if(thetaeld.lt.23.)then
          endif
          pel = pel*(1.+(phiel+pi)*0.0 +(phiel+pi)**2*0.04)
            endif
            if(thetaeld.le.18.)pel=pel*(1.-(18.-thetaeld)*0.0015)
            endif

           if(secte.eq.5)then
        if(thetaeld.lt.23.)then
        pel = pel*(1.-(phiel+2.*pi/3.)*0.005+(phiel+2.*pi/3.)**2*0.035)
        endif
           pel = pel*(1.-0.0088)
           if(thetaeld.le.18.)pel=pel*(1.-(18.-thetaeld)*0.0025)
           endif
           
           if(secte.eq.6)then
            pel = pel*(1.+0.0025)
         if(thetaeld.lt.23.)pel = pel*(1. + (phiel+pi/3.)**2*0.03)
         if(thetaeld.le.18.)pel=pel*(1.-(18.-thetaeld)*0.0015)
         if(thetaeld.ge.20.)pel = pel*(1. + (20.-thetaeld)*0.0002)
           endif  
             
            thetaeld=thetael*90./acos(0.)   
            endif

********************************************************************
*
*eg1 run 2.565GeV
*
********************************************************************         
         if(ebeam.ge.2.4.and.ebeam.le.2.7)then
* 50% field run         

         if(ifirste.eq.1) then
         print *,' Corrections for eg1 2.4-2.7 GeV, torus 50% '
         endif

           if(secte.eq.1)then 
            pel = pel*(1. + 0.015)
             endif
           
           if(secte.eq.2)then
             pel = pel*(1. + 0.01)
           endif
             
           if(secte.eq.3)then
           pel = pel*(1. - 0.020)   
           endif

           if(secte.eq.4)then
           if(phiel.lt.pi.and.phiel.gt.0.)then
           endif
           if(phiel.gt.-pi.and.phiel.lt.0.)then
           endif           
            pel = pel*(1. + 0.006)
             endif

           if(secte.eq.5)then
          pel = pel*(1.+ 0.006)
           endif

           if(secte.eq.6)then
          pel = pel*(1. - 0.012)
           endif  
 
          endif

* 
*     4 GeV run
*          
          if(ebeam.ge.4.0.and.ebeam.le.4.5)then
* 86% field run
           if(torcur.gt.3200.and.torcur.lt.3400)then

         if(ifirste.eq.1) then
         print *,' Corrections for e1adeve 4.0-4.5 GeV, torus 3300 A '
         endif

           if(secte.eq.1)then 
           if(thetaeld.lt.45.)pel=pel*(1. - phiel*0.035 + phiel**2*0.04)
            pel = pel*(1.+ 0.001)
            if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0011)
             endif
             
           if(secte.eq.2)then
            pel = pel*(1.- 0.002)
         if(thetaeld.lt.45.)then
         pel = pel*(1. -(phiel-pi/3.)*0.018 + (phiel-pi/3.)**2*0.04)
         endif
          if(thetaeld.le.19.)pel = pel*(1.-(19.-thetaeld)*0.0018)
         endif

          if(secte.eq.3)then
      if(thetaeld.lt.45.)then
      pel = pel*(1.+(phiel-2.*pi/3.)*0.008+(phiel-2.*pi/3.)**2*0.02)
      endif
       pel = pel*(1.- 0.004)
        if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.001)
       endif     

           if(secte.eq.4)then
             pel = pel*(1.- 0.001)
             if(thetaeld.lt.45.)then
            if(phiel.lt.pi.and.phiel.gt.0.)then
               pel = pel*(1. + (phiel-pi)*0.008+(phiel - pi)**2*0.08)
            endif
            if(phiel.gt.-pi.and.phiel.lt.0.)then
               pel = pel*(1.-(phiel+pi)*0.008 +(phiel+pi)**2*0.08)
            endif
            endif
         if(thetaeld.le.18.)pel=pel*(1. - (18.-thetaeld)*0.0014)
          endif

           if(secte.eq.5)then
         if(thetaeld.lt.45.)then
         pel = pel*(1.-(phiel+2.*pi/3.)*0.005+(phiel+2.*pi/3.)**2*0.06)
         endif
          pel = pel*(1.-0.003)
        if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0022)
         endif

           if(secte.eq.6)then
            if(thetaeld.lt.45.)then 
           pel = pel*(1.+(phiel+pi/3.)*0.01+(phiel+pi/3.)**2*0.06)
           endif
            pel = pel*(1. - 0.0004)
           if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0010)
            endif
            endif  
*            
*  60% B-field run                     
*
            if(torcur.gt.2200.and.torcur.lt.2300.)then

         if(ifirste.eq.1) then
         print *,' Corrections for e1adeve 4.0-4.5 GeV, torus 2250 A '
         endif

           if(secte.eq.1)then 
           pel=pel*(1. - phiel*0.05 + phiel**2*0.05)
           pel = pel*(1.- 0.0140)
         if(thetaeld.lt.18.)pel = pel*(1. - (18.-thetaeld)*0.0015)
         if(thetaeld.lt.15.)pel = pel*(1. - (15.-thetaeld)*0.002)
             endif
             
           if(secte.eq.2)then
        pel = pel*(1.- 0.005)
         pel = pel*(1. -(phiel-pi/3.)*0.020 + (phiel-pi/3.)**2*0.10)
          if(thetaeld.le.18.)pel = pel*(1.-(18.-thetaeld)*0.0015)
          if(thetaeld.le.15.)pel = pel*(1.-(15.-thetaeld)*0.0015)
         endif

          if(secte.eq.3)then
      pel = pel*(1.+(phiel-2.*pi/3.)*0.008+(phiel-2.*pi/3.)**2*0.07)
      pel = pel*(1.+ 0.001)
        if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0015)
        if(thetaeld.le.15.)pel = pel*(1. - (15.-thetaeld)*0.000)
       endif     

           if(secte.eq.4)then
             pel = pel*(1.- 0.0105)
            if(phiel.lt.pi.and.phiel.gt.0.)then
       pel = pel*(1.+(phiel-pi)*0.0+abs((phiel-pi))**2.*0.0)
            endif
            if(phiel.gt.-pi.and.phiel.lt.0.)then
       pel = pel*(1.-(phiel+pi)*0.0 +abs((phiel+pi))**2.*0.0)
            endif
         if(thetaeld.le.18.)pel=pel*(1. - (18.-thetaeld)*0.0020)        
         if(thetaeld.le.15.)pel=pel*(1. - (15.-thetaeld)*0.0020)
            endif

           if(secte.eq.5)then
         if(thetaeld.lt.23.)then
         pel = pel*(1.-(phiel+2.*pi/3.)*0.0+(phiel+2.*pi/3.)**2*0.10)
         endif
         pel = pel*(1.-0.0130)
        if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0022)
        if(thetaeld.le.15.)pel = pel*(1. - (15.-thetaeld)*0.0025)
         endif

           if(secte.eq.6)then
            if(thetaeld.lt.23.)then 
           pel = pel*(1.+(phiel+pi/3.)*0.015+(phiel+pi/3.)**2*0.10)
           endif
            pel = pel*(1. + 0.0033)
           if(thetaeld.le.18.)pel = pel*(1. - (18.-thetaeld)*0.0012)
           if(thetaeld.le.16.5)pel = pel*(1. - (16.5-thetaeld)*0.0012)
            endif
            endif

            thetaeld=thetael*90./acos(0.)   

          endif

******************************************************************
*  corrections end here        


        pmom = pel


        ELSEIF(ch1.eq.'H') THEN
        IF(ifirsth.le.1) ifirsth=ifirsth+1

* momentum correction hadron

        cz=cos(theta*pi180)

               secth=sector
               thetah=acos(cz)              
                thetahd=thetah*90./acos(0.)

         if(ifirsth.eq.1) then
         print *,' Corrections for Hadrons '
         endif


              if(q_h.gt.0)then
                if(ebeam.gt.1.5.and.ebeam.lt.1.7)then
                if(secth.eq.1)p_h=p_h*1.012
                if(secth.eq.2)p_h=p_h*1.008
                if(secth.eq.3)p_h=p_h*0.995
                if(secth.eq.4)p_h=p_h*0.993*(1.-0.001/p_h)
                if(secth.eq.5)p_h=p_h*1.0
                if(secth.eq.6)p_h=p_h*0.991
                endif
               if(ebeam.gt.2.2.and.ebeam.lt.2.7)then
               if(secth.eq.1)p_h=p_h*1.005
               if(secth.eq.2)p_h=p_h*1.005
                if(secth.eq.3)p_h=p_h*0.995
                if(secth.eq.4)p_h=p_h*0.997*(1.-0.001/p_h)
                if(secth.eq.5)p_h=p_h*1.0
                if(secth.eq.6)p_h=p_h*0.993
               endif
               
	       if(ebeam.gt.4.0.and.ebeam.lt.4.6)then 
	  if(secth.eq.1)then
	     p_h=p_h*(1.- 0.00)
	     if(thetahd.lt.30.)p_h=p_h*(1.+(30.-thetahd)**2/100000.) 
	  endif
	  if(secth.eq.2)then
	     p_h=p_h*(1.+ 0.005)
	     if(thetahd.lt.30.)p_h=p_h*(1.+(30.-thetahd)**2/200000.) 
	  endif
	  if(secth.eq.3)then
	     p_h=p_h*(1.- 0.00)
	      if(thetahd.lt.30.)p_h=p_h*(1.-(30.-thetahd)**2/100000.) 
	  endif
	 if(secth.eq.4)then
	     p_h=p_h*(1.- 0.00)
	    if(thetahd.lt.30.)p_h=p_h*(1.-(30.-thetahd)**2/100000.) 
         endif
         if(secth.eq.5)then
            p_h=p_h*(1.- 0.005)
            if(thetahd.lt.30.)p_h=p_h*(1.-(30.-thetahd)**2/100000.)
         endif
         if(secth.eq.6)then
            p_h=p_h*(1.- 0.002)
            if(thetahd.lt.30.)p_h=p_h*(1.+(30.-thetahd)**2/200000.)
         endif
               endif     

              endif
              if(q_h.lt.0)then      
                p_h = p_h*1.000
              endif
               
        pmom = p_h

        ELSE
          print *,' ch1 should be E or H'
          stop
        ENDIF !(ch1.eq.'E')



        END

