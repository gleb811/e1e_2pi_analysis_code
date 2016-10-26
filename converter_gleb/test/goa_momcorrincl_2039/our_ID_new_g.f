
      subroutine our_ID_new_g
c----------------------------------------------------------
c	"homemade" particle ID!
c----------------------------------------------------------
      implicit NONE	
      include "ntpl_goa.inc"
      include "photon.inc"
      include "pdgid.inc"	

      real*4 sctime,tr_time,dtime,p1,beta,bcor1
      real*4 t0sc,best,bsta,best1
      integer*4 isect,ipad,j
      real*4 beta_mass(3),t_exp(3,10),t_hyp(3,10)
      real*4 deltat_min,deltat_min_tot, deltat,dmt,random
      integer*4 jj,ll
      INTEGER*4  ipr, maxloop
      REAL*4  bdiff,beta_cut, b_part(9), delta_beta 
      REAL*4  mass_TG(20),mass_ST(20)
      integer ncharged


c+ Montecarlo generated events
c      if (skip_our_id) then
c         DO j=1,nEVNT
c           id_our(j)  = id_EVNT(j)
c           mass_our(j)= mass_EVNT(j)
c           beta_our(j)= beta_EVNT(j)
c         ENDDO
c         return
c      endif


	do j=2,nEVNT
	  mass_our(j) = -979      ! init mass_our
	  id_our(j) = -979        !  init id_our
c	   if(id_evnt(j).eq.2212) type *, charge_evnt(j),SCstat_EVNT(j),n_ph_in
	  if(charge_evnt(j).ne.0.and.SCstat_EVNT(j).gt.0.and.n_ph_in.gt.0) then
	   tr_time = ph_STT_HEVT(1)
           isect=sector_SCPB(SCstat_EVNT(j))
           ipad=pd_id_SCPB(SCstat_EVNT(j))
           if(ipad.eq.0.or.ipad.gt.40) goto 22
           sctime=time_SCPB(SCstat_EVNT(j))
c+to be corrected with ci(indT); now using a mean value C_tagE
           dtime=sctime - ph_TPHO_TAGR(1)
           p1=Pmom_EVNT(j)
c+to be corrected, using some target pvertex rec. inplace pf path_...
           beta= Path_SCPB(SCstat_EVNT(j))/29.97924/dtime
           bcor1=1./(1./beta+29.97924*t01/Path_SCPB(SCstat_EVNT(j)))

c+++ New particle id ctching timing only from tagger and correcting every TOF
	   t0sc= t0_i(isect,ipad)
           best1=Path_SCPB(SCstat_EVNT(j))/29.97924/(dtime+t0sc) ! Beta_best according Tagger
           best=Path_SCPB(SCstat_EVNT(j))/29.97924/(0.96+0.04*best1)/(dtime+t0sc) ! Beta_best according Tagger
           bsta=Path_SCPB(SCstat_EVNT(j))/29.97924/(sctime-tr_time-tvsst+t0sc)
	   if(best.lt.0) goto 22
           Mass_TG(J) = p1*p1* (1./best/best-1.)
           Mass_ST(j)= p1*p1* (1./bsta/bsta-1.)
	   call hfill(7208,pmom_EVNT(j),best,1.)
	   call hfill(7209,sqrt(abs(mass_TG(j))),0.,1.)
	   call hfill(7109,sqrt(abs(mass_EVNT(j))),0.,1.)

	   call hfill(7212, ((isect-1)*50+ipad)*1. ,sqrt(abs(mass_TG(j)))   ,1.)
	   call hfill(7112, ((isect-1)*50+ipad)*1. ,sqrt(abs(mass_evnt(j))) ,1.)
	   beta_our(j) =  best
	   mass_our(j) =  Mass_TG(J) !Output 
c	   beta_evnt(j) =  best !Output 

c	   type *, j, mass_evnt(j),mass_our(j), Mass_TG(J)
c+ Random coinc cleaning for 2 charged detected
c+ correcting for extra +-2ns bunch (putting the result in mass_st(j) )
	   beta_mass(1) = p1/sqrt(p1**2+.13957**2)
	   beta_mass(2) = p1/sqrt (p1**2+.494**2)
	   beta_mass(3)  = p1/sqrt(p1**2+.93828**2)
	   t_hyp(1,j)  = Path_SCPB(SCstat_EVNT(j))/(29.97924*beta_mass(1))
	   t_hyp(2,j)  = Path_SCPB(SCstat_EVNT(j))/(29.97924*beta_mass(2))
	   t_hyp(3,j)  = Path_SCPB(SCstat_EVNT(j))/(29.97924*beta_mass(3))
           t_exp(1,j)  = (dtime+t0sc)*(0.96+0.04*best1)
	   t_exp(2,j)  = (dtime+t0sc)*(0.96+0.04*best1) + 2
	   t_exp(3,j)  = (dtime+t0sc)*(0.96+0.04*best1) - 2
c	    deltat_min=10000.
c	   do jj = 1,3
c	      do ll =1,3
c	      deltat=abs(t_exp(jj)-t_hyp(ll))
c	       if(deltat.lt.deltat_min) then
c		  deltat_min = deltat
c		  j_part_good = ll
c		  j_time_good = jj
c		endif
c	       enddo
c	    enddo
c	    beta_exp_corr =Path_SCPB(SCstat_EVNT(j))/29.97924/t_exp(j_time_good)
c	    mass_st(j) = p1*p1* (1./beta_exp_corr/beta_exp_corr-1.) !writing on m_st
c- got mass_corrected in mass_st
c- ending corr for +- 2ns bunches

c+ Calculating pid number following the usual scheme (using BETA_BEST <-> Tagger)
              bdiff = 1.
              beta_cut = 1.
              maxloop = 9
              if(charge_evnt(j).lt.0) maxloop = 4
              Do  ipr = 2,maxloop
               b_part(ipr) = p1 / sqrt(p1**2+p_mass(ipr)**2)  ! Beta according to 1p 
               delta_beta = ABS(best - b_part(ipr))
c- Take the closest by beta particle
               if(delta_beta.lt.beta_cut.and.delta_beta.lt.bdiff)  then  
                bdiff=delta_beta
                id_our(j) =p_id(ipr)*charge_EVNT(j)
		if (id_our(j).eq.2112) id_our(j) = 2212                  ! correcting n->p (ch<>0)
                if (abs(id_our(j)).eq.11) id_our(j)=charge_EVNT(j)*211   ! No electron identified: ONLY PIONS
              endif
21             continue
              end do !Part id loop over known particles
	   elseif(charge_evnt(j).eq.0) then !ending our id routine
22	    continue
	    mass_our(j) = mass_EVNT(j)
	    id_our(j) = ID_EVNT(j)
	  ENDIF
	  enddo ! next particle in EVNT bank

c++ Calculating TOF/TAGGER delay assuming 1st approx pid
c-- Ending  TOF/TAGGER delay
c+ id_our(1)=  0 ->  not checked
c+ id_our(1)= -1 ->      checked TRUE
c+ id_our(1)= +1 ->      checked RANDOM
	  random=0. ! not checked: true or random
	  do jj = 1,3
	     deltat_min_tot = 0. ! Fixed T = 0 +2 -2
	     ncharged = 0
	     do j=2,nEVNT	! Looping over particles 
	      if(charge_evnt(j).ne.0.and.SCstat_EVNT(j).gt.0.and.n_ph_in.gt.0.) then
	       deltat_min=1000.
	       ncharged = ncharged +1
	       do ll =1,3	! Finding correct mass
		  if(t_exp(jj,j).le.0) then
		    ncharged=ncharged-1
                    goto 1099
		   endif
	          deltat = abs(t_exp(jj,j)-t_hyp(ll,j))
	        if(deltat.lt.deltat_min)  deltat_min = deltat
	       enddo 
	       deltat_min_tot = deltat_min_tot + deltat_min
              endif
 1099	   enddo
	     if(jj.eq.1) dmt = deltat_min_tot ! T=0 : good event
	     if(ncharged.ge.2)	random=-1 ! checked and true
	     if(deltat_min_tot.lt.dmt.and.ncharged.ge.2) random=1 ! if T=+-2 gives a better delta it's a random    
 	  enddo
	  id_our(1) = random
c--
C++  KAONS ID: forcing kaons to pions and then kaons are identified using a narrow cut around k_mass 

	  do j=2,nEVNT
	  if(id_our(j).eq. 321) then
               id_our(j)= 211
	       if (beta_OUR(j).lt.pmom_evnt(j)/sqrt(pmom_evnt(j)**2+(.4936-0.045)**2).and.
     %             beta_OUR(j).gt.pmom_evnt(j)/sqrt(pmom_evnt(j)**2+(.4936+0.085)**2))
     %         id_our(j)= 321
	  endif
	  if(id_our(j).eq.-321)  then
               id_our(j)=-211
	       if (beta_OUR(j).lt.pmom_evnt(j)/sqrt(pmom_evnt(j)**2+(.4936-0.045)**2).and.
     %             beta_OUR(j).gt.pmom_evnt(j)/sqrt(pmom_evnt(j)**2+(.4936+0.085)**2))
     %         id_our(j)= -321
	  endif
	 enddo
	     
c--
	return
	end	
