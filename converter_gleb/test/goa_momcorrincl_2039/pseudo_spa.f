
      
      SUBROUTINE pseudo_spa2(part_t,pid,
     &                       p,theta,phi,t_current,acp,sect,acc)
      IMPLICIT none
      INCLUDE "ntpl_goa.inc"
      REAL p,theta,phi,t_current,acp 
      INTEGER pid,part_t,sect, ib     
      REAL acc
      REAL spa_e1a_2445_1500,spa_e1b_2567_1500,spa_skewed_e1a_2445_1500
      REAL spa_sknew_e1b_2567_1500,spa_sknew_e1b_4247_2250
      INTEGER n1
      DATA n1/0/

	if(beam_type) then
         ib = 1
	else
	 ib=2
	endif

      n1=n1+1
      if(part_t.eq.-1000) acc=0.
      if(part_t.eq.-1000) return

      IF( beam_type .and.
     &    abs(Eelbeam-2.445).lt.0.01.and.abs(I_torus-1500.).lt.0.01)THEN
      if(n1.eq.1) print *,'INFO: acceptance for E1A 2445-1500 to use'
c      acc = spa_e1a_2445_1500(
c     &        part_t,pid,p,theta,phi,ib,t_current,acp,sect)
      acc = spa_skewed_e1a_2445_1500(
     &        part_t,pid,p,theta,phi,ib,t_current,acp,sect)

      ELSEIF( beam_type .and.
     &    abs(Eelbeam-2.567).lt.0.01.and.abs(I_torus-1500.).lt.0.01)THEN
      if(n1.eq.1) print *,'INFO: acceptance for E1B 2567-1500 to use'
      acc = spa_sknew_e1b_2567_1500(
     &        part_t,pid,p,theta,phi,ib,t_current,acp,sect)

      ELSEIF( beam_type .and.
     &    abs(Eelbeam-4.247).lt.0.01.and.abs(I_torus-2250.).lt.0.01)THEN
      if(n1.eq.1) print *,'INFO: acceptance for E1B 4.247-2250 to use'
      acc = spa_sknew_e1b_4247_2250(
     &        part_t,pid,p,theta,phi,ib,t_current,acp,sect)

      ELSE
      if(n1.eq.1) print *,'INFO: acceptance for E1A ini(Burk) to use'
      CALL pseudo_spa_ini(part_t,p,theta,phi,ib,t_current,acp,sect,acc)


      ENDIF
      RETURN
      END




c======================================================================
c                 pseudo_spa -- initial Burkert's version
c======================================================================

      SUBROUTINE 
     &  pseudo_spa_ini(part_t,p,theta,phi,ib,t_current,acp,sect,acc)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Inputs:
c-          part_t - INTEGER 0=electron, 1=positive hadron, -1=negative hadron,
c-          p - REAL particle momentum GeV/c
c-          theta - REAL particle polar angle in degrees (0 to 180.)
c-          phi - REAL particle azimuthal angle in degrees in the sector -30. to 30.
c-          ib - ibeam type
c-          t_current - REAL torus current in A.
c-	    acp - REAL portion of phi acceptance that should be accepted.  
c-          sect - INTEGER sector number.    
c-
c-  Output:
c-          acc - REAL acceptans 0. to 1.
c-
      IMPLICIT NONE
      REAL p,theta,phi,t_current,acp    ! Input parameter
      INTEGER part_t,ib,sect                 ! 
c
      REAL acc                          ! Output parameter
c
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, delta_phi, exp, pnorm,dtheta
c-      
	if(ib.eq.1) then
	 theta0_nh = 15. ! Electron beam
         theta0_ph = 10.
	else
         theta0_nh = 10. ! Photon beam
         theta0_ph = 10.
	endif



      Acc=0.0
      pnorm=p*t_max/t_current
      IF(part_t.EQ.electron) THEN
       if(ib.eq.1) then
        if(sect.eq.2)then
          dtheta=6./pnorm
          if(theta.gt.18.+dtheta.and.theta.lt.20.+dtheta)return
        elseif(sect.eq.5)then
          dtheta=6./(p*t_max/t_current)
          if(theta.gt.19.+dtheta.and.theta.lt.21.5+dtheta)return
        endif
       endif
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.50.)then
          exp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
        elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
          exp = cel_ex_nip*(p*t_max/t_current)**pel_ex
          delta_phi = acp*phi0_el_nip
     *   *sin((theta-theta_min+theta_nip)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSEIF(part_t.EQ.pos_hadron) THEN
       if(ib.eq.1) then
        if(sect.eq.3)then
	if(pnorm.gt.0.825+1.5*(theta/35.)**6.and.
     *	(pnorm.lt.1.0125+1.35*(theta/32.)**6))return
	if(pnorm.gt.0.15+ 0.75*(theta/65.)**8.and.
     *	(pnorm.lt.0.15+0.75*(theta/50.)**6))return
        elseif(sect.eq.4)then
        if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)return
        if(pnorm.gt.0.7+1.5*(theta/41.5)**8.and.
     *  pnorm.lt.0.85+1.5*(theta/35.5)**8)return
        elseif(sect.eq.5)then
	if(pnorm.gt.1.+2.*(theta/20.)**2.and.
     *	(pnorm.lt.1.2+1.8*(theta/15.)**2))return
        if(pnorm.gt.0.5+0.7*(theta/40.)**8.and.
     *  pnorm.lt.0.7+2.*(theta/40.)**8)return
        elseif(sect.eq.6)then
         if(pnorm.gt.0.2+0.5*(theta/80.)**6.and.
     *    pnorm.lt.0.25+0.5*(theta/70.)**6)return
        endif
       endif
        theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_current/5.)**24
        if(theta.gt.theta_min)then
          exp=(p*t_max/t_current/5.)**(1./8.)
          delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSEIF(part_t.EQ.neg_hadron) THEN
        theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.130.)then
          exp = ch_ex*(p*t_max/t_current)**pim_ex
          delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
        endif
      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)WRITE(6,*)'Illegal particle type ',part_t
      		ACC=0.0
      ENDIF
      RETURN 
      END


c======================================================================
c                 pseudo_spa -- e1a_2445_1500
c======================================================================
      FUNCTION spa_e1a_2445_1500(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptans 0. to 1.
c-
      IMPLICIT NONE
      REAL acc, spa_e1a_2445_1500
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, delta_phi, exp, pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.2)then
         dtheta=6./pnorm
         if(theta.gt.18. +dtheta.and.theta.lt.20.0 +1.+dtheta)goto 9000
       elseif(sect.eq.3)then
         dtheta=6./pnorm
         if(theta.gt.39.0 +dtheta.and.theta.lt.45.1 +dtheta)goto 9000
       elseif(sect.eq.4)then
         dtheta=6./pnorm
         if(theta.gt.39.0 +dtheta.and.theta.lt.45.1 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN
      if(ibeam_type.eq.1) then
cc      IF(pid.eq.2212)THEN ! protons only
cc        if(pnorm.gt. 0.5  +(1.50    )*(theta/48.)**6)then
cc        if(pnorm.lt. 0.8  +(1.35+0.4)*(theta/43.)**6)goto 9000
cc        endif
cc      ENDIF
       if(sect.eq.3)then
c         if(pnorm.gt. 0.8250       +(1.50    )*(theta/35.)**6)then
c         if(pnorm.lt. 1.0125 +0.08 +(1.35+0.4)*(theta/32.)**6)goto 9000
c         endif
         if(pnorm.gt. 0.15 +(0.750)*(theta/65.)**8)then
         if(pnorm.lt. 0.15 +(0.750)*(theta/50.)**6)goto 9000
         endif
       elseif(sect.eq.4)then
         if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)goto 9000
c         if(pnorm.gt. 0.70      +(1.50)*(theta/41.5)**8)then
c         if(pnorm.lt. 0.85 +0.1 +(1.50)*(theta/35.5)**8)goto 9000
c         endif
       elseif(sect.eq.5)then
         if(pnorm.gt.1.00 -0.4 +(2.0 -0.3)*(theta/20.)**2)then
         if(pnorm.lt.1.20      +(1.8     )*(theta/15.)**2)goto 9000
         endif
c         if(pnorm.gt.0.50 -0.2 +(0.7)*(theta/40.)**6)then
c         if(pnorm.lt.0.70 +0.2 +(2.0)*(theta/40.)**6)goto 9000
c         endif
c       elseif(sect.eq.6)then
c         if(pnorm.gt.0.20 -0.1 +(0.50    )*(theta/80.)**6)then
c         if(pnorm.lt.0.25 -0.1 +(0.50+0.1)*(theta/70.)**6)goto 9000
c         endif
       endif
      endif
      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
        exp=(p*t_max/t_curr/5.)**(1./8.)
        delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue
      spa_e1a_2445_1500=acc
      RETURN 
      END

c======================================================================
c                 skewed pseudo_spa -- e1a_2445_1500
c======================================================================
      FUNCTION spa_skewed_e1a_2445_1500(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptance 0. to 1.
c-
      IMPLICIT NONE
      REAL acc, spa_skewed_e1a_2445_1500
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min,delta_phi,delta_phi_plus,delta_phi_minus,exp1,pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.2)then
         dtheta=6./pnorm
         if(theta.gt.18. +dtheta.and.theta.lt.20.0 +1.+dtheta)goto 9000
       elseif(sect.eq.3)then
         dtheta=6./pnorm
         if(theta.gt.39.0 +dtheta.and.theta.lt.45.1 +dtheta)goto 9000
       elseif(sect.eq.4)then
         dtheta=6./pnorm
         if(theta.gt.39.0 +dtheta.and.theta.lt.45.1 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp1 = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp1 = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN
      if(ibeam_type.eq.1) then
cc      IF(pid.eq.2212)THEN ! protons only
cc        if(pnorm.gt. 0.5  +(1.50    )*(theta/48.)**6)then
cc        if(pnorm.lt. 0.8  +(1.35+0.4)*(theta/43.)**6)goto 9000
cc        endif
cc      ENDIF
       if(sect.eq.3)then
c         if(pnorm.gt. 0.8250       +(1.50    )*(theta/35.)**6)then
c         if(pnorm.lt. 1.0125 +0.08 +(1.35+0.4)*(theta/32.)**6)goto 9000
c         endif
         if(pnorm.gt. 0.15 +(0.750)*(theta/65.)**8)then
         if(pnorm.lt. 0.15 +(0.750)*(theta/50.)**6)goto 9000
         endif
       elseif(sect.eq.4)then
         if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)goto 9000
c         if(pnorm.gt. 0.70      +(1.50)*(theta/41.5)**8)then
c         if(pnorm.lt. 0.85 +0.1 +(1.50)*(theta/35.5)**8)goto 9000
c         endif
       elseif(sect.eq.5)then
         if(pnorm.gt.1.00 -0.4 +(2.0 -0.3)*(theta/20.)**2)then
         if(pnorm.lt.1.20      +(1.8     )*(theta/15.)**2)goto 9000
         endif
c         if(pnorm.gt.0.50 -0.2 +(0.7)*(theta/40.)**6)then
c         if(pnorm.lt.0.70 +0.2 +(2.0)*(theta/40.)**6)goto 9000
c         endif
c       elseif(sect.eq.6)then
c         if(pnorm.gt.0.20 -0.1 +(0.50    )*(theta/80.)**6)then
c         if(pnorm.lt.0.25 -0.1 +(0.50+0.1)*(theta/70.)**6)goto 9000
c         endif
       endif
      endif
      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
       if(sect.eq.1)then
       delta_phi_plus=23.*(1.-exp(-0.09*(theta-9.)))
       delta_phi_minus=22.*(1.-exp(-0.12*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.2)then
       delta_phi_plus=24.*(1.-exp(-0.09*(theta-8.)))
       delta_phi_minus=21.5*(1.-exp(-0.12*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.3)then
       delta_phi_plus=23.*(1.-exp(-0.10*(theta-8.)))
       delta_phi_minus=25.*(1.-exp(-0.10*(theta-9.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.4)then
       delta_phi_plus=22.*(1.-exp(-0.12*(theta-10.)))
       delta_phi_minus=23.*(1.-exp(-0.10*(theta-9.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.5)then
       delta_phi_plus=24.*(1.-exp(-0.12*(theta-9.)))
       delta_phi_minus=24.*(1.-exp(-0.08*(theta-9.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.6)then
       delta_phi_plus=23.*(1.-exp(-0.07*(theta-8.)))
       delta_phi_minus=25.*(1.-exp(-0.09*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       endif
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp1 = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue
      spa_skewed_e1a_2445_1500=acc
      RETURN 
      END

c======================================================================
c                 pseudo_spa -- e1b_2567_1500
c======================================================================
      FUNCTION spa_e1b_2567_1500(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptans 0. to 1.
c-
        IMPLICIT NONE
      REAL acc, spa_e1b_2567_1500
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, delta_phi, exp, pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.3)then
         dtheta=6. - 2.5*pnorm
         if(theta.gt.24.0 +dtheta.and.theta.lt.28.0 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN
      if(ibeam_type.eq.1) then
cc      IF(pid.eq.2212)THEN ! protons only
cc        if(pnorm.gt. 0.5  +(1.50    )*(theta/48.)**6)then
cc        if(pnorm.lt. 0.8  +(1.35+0.4)*(theta/43.)**6)goto 9000
cc        endif
cc      ENDIF
       if(sect.eq.4)then
         if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)goto 9000
       endif
      endif
      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
       if(sect.eq.1.AND.phi.lt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.1.AND.phi.gt.0.) then
       if(theta.gt.20) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.20) exp=(p*t_max/t_curr/0.01)**(1./8.)
       endif
       if(sect.eq.2.AND.phi.lt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.2.AND.phi.gt.0.) then
       if(theta.gt.18) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.18) exp=(p*t_max/t_curr/0.05)**(1./8.)
       endif
       if(sect.eq.3.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.3.AND.phi.lt.0.) then
       if(theta.gt.16) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.16) exp=(p*t_max/t_curr/0.05)**(1./8.)
       endif
       if(sect.eq.4.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.4.AND.phi.lt.0.) then
       if(theta.gt.22) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.22) exp=(p*t_max/t_curr/0.005)**(1./8.)
       endif
       if(sect.eq.5.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.5.AND.phi.lt.0.) then
       if(theta.gt.21) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.21) exp=(p*t_max/t_curr/0.001)**(1./8.)
       endif

       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp

       if(sect.eq.6.AND.phi.lt.0.) then
       if(theta.gt.25) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.25) exp=(p*t_max/t_curr/0.01)**(1./8.)
       endif
       if(sect.eq.6.AND.phi.gt.0.) then
       if(theta.gt.30) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.30) exp=(p*t_max/t_curr/0.03)**(1./8.)
       endif
       if(sect.eq.6.AND.phi.lt.0.) then
       if(theta.gt.25) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       endif
       if(theta.lt.25) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       delta_phi=delta_phi*(0.06*theta)+2.
       endif
       endif
       if(sect.eq.6.AND.phi.gt.0.) then
       if(theta.gt.30) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       endif
       if(theta.lt.30) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       delta_phi=delta_phi*(0.04*theta)
       endif
       endif
         if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue
      spa_e1b_2567_1500=acc
      RETURN 
      END

c======================================================================
c                 pseudo_spa -- e1b_4247_2250
c======================================================================
      FUNCTION spa_e1b_4247_2250(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptans 0. to 1.
c-
        IMPLICIT NONE
      REAL acc, spa_e1b_4247_2250
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min, delta_phi, exp, pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.3)then
         dtheta=6. - 2.5*pnorm
         if(theta.gt.24.0 +dtheta.and.theta.lt.28.0 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN
      if(ibeam_type.eq.1) then
cc      IF(pid.eq.2212)THEN ! protons only
cc        if(pnorm.gt. 0.5  +(1.50    )*(theta/48.)**6)then
cc        if(pnorm.lt. 0.8  +(1.35+0.4)*(theta/43.)**6)goto 9000
cc        endif
cc      ENDIF
       if(sect.eq.4)then
         if(pnorm.lt.theta/5..and.pnorm.gt.theta/7.)goto 9000
       endif
      endif
      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
       if(sect.eq.1.AND.phi.lt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.1.AND.phi.gt.0.) then
       if(theta.gt.20) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.20) exp=(p*t_max/t_curr/0.01)**(1./8.)
       endif
       if(sect.eq.2.AND.phi.lt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.2.AND.phi.gt.0.) then
       if(theta.gt.18) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.18) exp=(p*t_max/t_curr/0.05)**(1./8.)
       endif
       if(sect.eq.3.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.3.AND.phi.lt.0.) then
       if(theta.gt.16) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.16) exp=(p*t_max/t_curr/0.05)**(1./8.)
       endif
       if(sect.eq.4.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.4.AND.phi.lt.0.) then
       if(theta.gt.22) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.22) exp=(p*t_max/t_curr/0.005)**(1./8.)
       endif
       if(sect.eq.5.AND.phi.gt.0.) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(sect.eq.5.AND.phi.lt.0.) then
       if(theta.gt.21) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.21) exp=(p*t_max/t_curr/0.001)**(1./8.)
       endif

       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp

       if(sect.eq.6.AND.phi.lt.0.) then
       if(theta.gt.25) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.25) exp=(p*t_max/t_curr/0.01)**(1./8.)
       endif
       if(sect.eq.6.AND.phi.gt.0.) then
       if(theta.gt.30) exp=(p*t_max/t_curr/5.)**(1./8.)
       if(theta.lt.30) exp=(p*t_max/t_curr/0.03)**(1./8.)
       endif
       if(sect.eq.6.AND.phi.lt.0.) then
       if(theta.gt.25) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       endif
       if(theta.lt.25) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       delta_phi=delta_phi*(0.06*theta)+2.
       endif
       endif
       if(sect.eq.6.AND.phi.gt.0.) then
       if(theta.gt.30) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       endif
       if(theta.lt.30) then
       delta_phi = acp*phi0_ph*cos((theta-theta_cut)*d2r)**exp
       delta_phi=delta_phi*(0.04*theta)
       endif
       endif
         if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue
      spa_e1b_4247_2250=acc
      RETURN 
      END

c======================================================================
c                 sknew pseudo_spa -- e1b_2567_1500
c======================================================================
      FUNCTION spa_sknew_e1b_2567_1500(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptance 0. to 1.
c-
      IMPLICIT NONE
      REAL acc, spa_sknew_e1b_2567_1500
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min,delta_phi,delta_phi_plus,delta_phi_minus,exp1,pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.3)then
         dtheta=6. - 2.5*pnorm
         if(theta.gt.24.0 +dtheta.and.theta.lt.28.0 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp1 = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp1 = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN

      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
       if(sect.eq.1)then
       delta_phi_plus=23.*(1.-exp(-0.09*(theta-10.)))
       delta_phi_minus=24.*(1.-exp(-0.10*(theta-5.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.2)then
       delta_phi_plus=24.*(1.-exp(-0.11*(theta-9.)))
       delta_phi_minus=22*(1.-exp(-0.18*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.3)then
       delta_phi_plus=23.*(1.-exp(-0.17*(theta-7.)))
       delta_phi_minus=25.*(1.-exp(-0.11*(theta-9.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.4)then
       delta_phi_plus=24.*(1.-exp(-0.11*(theta-6.)))
       delta_phi_minus=23.*(1.-exp(-0.08*(theta-11.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.5)then
       delta_phi_plus=24.*(1.-exp(-0.12*(theta-8.)))
       delta_phi_minus=23.*(1.-exp(-0.09*(theta-11.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.6)then
       delta_phi_plus=21.*(1.-exp(-0.08*(theta-10.)))
       delta_phi_minus=24.*(1.-exp(-0.10*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       endif
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp1 = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue

      if(acc.gt.0) acc=1
      spa_sknew_e1b_2567_1500=acc
      RETURN 
      END

c======================================================================
c                 sknew pseudo_spa -- e1b_4247_2250
c======================================================================
      FUNCTION spa_sknew_e1b_4247_2250(p_t,pid,p,theta,phi,ib,t_curr,acp,sect)
c-      
c-  Pseudo Single Particle acceptance function from Volker Burkert.
c-  Input:
c-        p_t = part_t - INTEGER 0=electron, 1=positive hadron, 
c-                              -1=negative hadron,
c-        pid - INTEGER -11=electron, 2212=proton, 
c-                0=undefined (corrections are applied as part_t),
c-        p - REAL particle momentum GeV/c
c-        theta - REAL particle polar angle in degrees (0 to 180.)
c-        phi - REAL particle azimuthal angle in degrees in the 
c-              sector -30. to 30.
c-        ib = ibeam_type 1-electrons; 2-photons
c-        t_curr - REAL torus current in A.
c-        acp - REAL portion of phi acceptance that should be accepted.
c-        accepted.  
c-        sect - INTEGER sector number.    
c-
c-  Output:
c-        spa - REAL acceptance 0. to 1.
c-
      IMPLICIT NONE
      REAL acc, spa_sknew_e1b_4247_2250
      REAL p,theta,phi,t_curr,acp 
      INTEGER part_t,p_t,pid,sect,ibeam_type,ib
      REAL t_max
      parameter (t_max = 3375.)
      REAL phi0_el, phi0_nh, phi0_ph
      parameter (phi0_el = 30.)
      parameter (phi0_nh = 25.)
      parameter (phi0_ph = 25.)    ! Changed for better discription.
      REAL theta0_el, theta0_nh, theta0_ph
      parameter (theta0_el = 15.5)
cc      parameter (theta0_nh = 10.) !  Different for electrons 
cc      parameter (theta0_ph = 10.) !  and photons
      REAL thetas_el, thetas_nh, thetas_ph
      parameter (thetas_el = 15.)
      parameter (thetas_nh = 15.)
      parameter (thetas_ph = 25.)
      REAL p_shift, cel_ex, pel_ex, pim_ex
      parameter (p_shift = 0.05)
      parameter (pel_ex = 0.333)
      parameter (pim_ex = 0.5)
      parameter (cel_ex = 0.35)
      REAL  ch_ex,theta_cut
      parameter (theta_cut = 75.)
      parameter (ch_ex = 0.1)
      INTEGER electron,pos_hadron, neg_hadron
      parameter (electron = 0)
      parameter (pos_hadron = 1)
      parameter (neg_hadron = -1)
c- New, for very forward electrons -
      REAL theta_nip,cel_ex_nip,phi0_el_nip
      parameter (theta_nip = 2.)
      parameter (cel_ex_nip = 0.5)
      parameter (phi0_el_nip = 20.)
c
      REAL pi,d2r,u_acc
      parameter (pi = 3.1415926)
      parameter (d2r = 0.0174533)
      parameter (u_acc = 0.20944)
      INTEGER err_count
      data err_count/0/
      REAL theta_min,delta_phi,delta_phi_plus,delta_phi_minus,exp1,pnorm,dtheta
      
      part_t = p_t
      ibeam_type = ib

      if(ibeam_type.eq.1) then
	theta0_nh = 15. ! Electron beam
        theta0_ph = 10.
      else
        theta0_nh = 10. ! Photon beam
        theta0_ph = 10.
      endif

      Acc=0.0
      pnorm=p*t_max/t_curr

c------------------------- ELECTRONS -----------------------------
      IF(part_t.EQ.electron) THEN
      IF(ibeam_type.eq.1) then
       if(sect.eq.3)then
         dtheta=6. - 2.5*pnorm
         if(theta.gt.24.0 +dtheta.and.theta.lt.28.0 +dtheta)goto 9000
       elseif(sect.eq.5)then
         dtheta=6./pnorm
         if(theta.gt.19. +dtheta.and.theta.lt.21.5 +1.+dtheta)goto 9000
       endif
      ENDIF
      theta_min = theta0_el+thetas_el/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.50.)then
        exp1 = cel_ex*(p*t_max/t_curr)**pel_ex
        delta_phi = acp*phi0_el*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
c- New, for very forward electrons -
      elseif(theta.gt.theta_min-theta_nip.and.theta.lt.50.)then
        exp1 = cel_ex_nip*(p*t_max/t_curr)**pel_ex
        delta_phi= acp*phi0_el_nip
        delta_phi= delta_phi*sin((theta-theta_min+theta_nip)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


c------------------------- POSITIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.pos_hadron) THEN
      theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_curr/5.)**24
      if(theta.gt.theta_min)then
       if(sect.eq.1)then
       delta_phi_plus=23.*(1.-exp(-0.07*(theta-11.)))
       delta_phi_minus=22.*(1.-exp(-0.18*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.2)then
       delta_phi_plus=24.*(1.-exp(-0.11*(theta-9.)))
       delta_phi_minus=21.5*(1.-exp(-0.18*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.3)then
       delta_phi_plus=23.*(1.-exp(-0.17*(theta-7.)))
       delta_phi_minus=25.*(1.-exp(-0.11*(theta-9.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.4)then
       delta_phi_plus=22.*(1.-exp(-0.18*(theta-6.)))
       delta_phi_minus=23.*(1.-exp(-0.08*(theta-11.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.5)then
       delta_phi_plus=24.*(1.-exp(-0.14*(theta-8.)))
       delta_phi_minus=24.*(1.-exp(-0.07*(theta-13.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       elseif(sect.eq.6)then
       delta_phi_plus=23.*(1.-exp(-0.06*(theta-12.)))
       delta_phi_minus=25.*(1.-exp(-0.09*(theta-8.)))
        if(phi.lt.delta_phi_plus.AND.phi.gt.(-delta_phi_minus))Acc=1
       endif
      endif


c------------------------- NEGATIVE HADRONS ---------------------------
      ELSEIF(part_t.EQ.neg_hadron) THEN
      theta_min = theta0_nh+thetas_nh/(p*t_max/t_curr+p_shift)
      if(theta.gt.theta_min.and.theta.lt.130.)then
        exp1 = ch_ex*(p*t_max/t_curr)**pim_ex
        delta_phi = acp*phi0_nh*sin((theta-theta_min)*d2r)**exp1
        if(abs(phi).lt.delta_phi)Acc=sin(theta*d2r)*delta_phi*u_acc
      endif


      ELSE
      	err_count=err_count+1
      	if(err_count.le.10)print *,'Illegal particle type ',part_t
        ACC=0.0
      ENDIF

 9000 continue

      if(acc.gt.0) acc=1
      spa_sknew_e1b_4247_2250=acc
      RETURN 
      END

