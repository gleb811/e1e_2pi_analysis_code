c
c
c List of the subroutines in this file:
c   processevent_goa      -- process every event
c   initevent_goa         -- initilyze before processing an event
c   prepfile_goa()        -- preprocessing file for electorns
c   prepfile_goa_photon() -- preprocessing file for photons
c   postfile_goa          -- postprocessing file for electrons
c
c




c======================================================================
c MCbit = 0 real data
c MCbit = 1 MC data
c====================================================================== 
      SUBROUTINE processevent_goa(Select_OK,Nstory,MCbit)

      IMPLICIT NONE
      include "ntpl_goa.inc"
      include "ntpl_clas.inc"
      include "select.inc"

      logical Select_OK, EgammaIsGood
      integer Nstory,MCbit,sect,ibuf
      integer j,point_pp,point_pm     
      real*4  m_nucleon(0:1),W2_local,pi,Eq_photon_E
      real*4  spool(50)
      integer ipool(50)
      integer max_pi,max_k,i,jj,j_l
      integer max_g,k, max_gg
      integer sector,IDP
      real*4  fcg_deltalast
      logical key_count,key_endfile,key_endblock
      integer Ifile_prev,Ifile_curr,Iblock_curr,Iblock_all
      integer ndt_ebeam,ndt_enotb,ndt_ep,ndt_p,ndt_n
      integer ndt_pip,ndt_pim,ndt_kp,ndt_km,ndt_d,ndt_g 
      integer n_e_incl,n_ep_incl,n_e_elast,n_ep_elast
      integer n_ep_delta,n_twopions
      integer n_parttot,n_evtot,n_evaccel
      integer n_partpassed1,cceffn,cceffn1
      real*4  n_partavrtot,n_partavrevnt
      real*4  dgamma
      real*4  raddeg, degrad
      real*4  fcg_prev,fcg_curr,Qdiff
      real*4  c,pp,bb,bbe,aa,mm,mme,eetot,eein,eeout
      integer key_cc
      real*4  pxyz(3),cxyz(3),nphe2,nphe2err
      character*1 ch1
      integer eloss,ip,status_eloss
      real*4  ve(3),cd(3),pm,pv,elast
      integer nbufmax_hel_ini,end_head

      integer*4  CCEFFIC
      COMMON /qqq/Iblock_all,n_e_elast,elast,n_e_incl

      data th_el,ph_EL,E_EL,P_EL,px_EL,py_EL,pz_EL,Q2/8*0./
      data pi,m_nucleon/3.1415926536,0.93955,0.93821/ !m_nucleon(0)=neutron, _(1)=proton
      data degrad/0.0174533/
      data raddeg/57.2958/
      data fcg_deltalast/0./
      data Ifile_prev,Ifile_curr,Iblock_curr/0,0,0/
      data ndt_ebeam,ndt_enotb,ndt_ep,ndt_p,ndt_n/0,0,0,0,0/
      data ndt_pip,ndt_pim,ndt_kp,ndt_km,ndt_d,ndt_g/0,0,0,0,0,0/
      data n_e_incl,n_ep_incl,n_e_elast,n_ep_elast/0,0,0,0/
      data n_ep_delta,n_twopions/0,0/
      data n_parttot,n_evtot,n_evaccel/0,0,0/
      data n_partpassed1/0/
      data n_partavrtot,n_partavrevnt/0.,0./
      data nbufmax_hel_ini/1/
      data dgamma/0./

c --- Select an event                                             ---
c --- when Select_OK is set to .false. this event will be skipped --- 
      IF(bit_mc.eq.0.and.beam_type) THEN ! No preprocess for MC or photon beam
        IF(Event_OK(Nstory)) THEN
          Select_OK = .true.
        ELSE
          Select_OK = .false.
          RETURN
        ENDIF
      ELSE
       select_OK = .true.
      ENDIF
      IF(nEVNT.gt.10) THEN
        Select_OK = .false.
        RETURN
      ENDIF       

c --- Set Qgated_total/plus Accumulating Qgated of all files ---
      DO ibuf=1,nbuf_read_fcg
      IF(indBegin_fcg(ibuf) .eq. Nstory) THEN
        IF(N_firstfile1.eq.0) THEN
        Qgated_total = Qgated_total + fcg_deltalast
        ELSE
        Qgated_total = Qgated_total + fcg_buf(ibuf)-fcg_buf(ibuf-1)
        ENDIF
        Qgated_plus  = Qgated_total + fcg_buf(ibuf+1)-fcg_buf(ibuf)
        fcg_deltalast = fcg_buf(ibuf+1)-fcg_buf(ibuf)
        N_firstfile1=1
        GOTO 1001
      ENDIF
      ENDDO
 1001 CONTINUE

c --- Set Helicitiy KEY (they were saved by buffer in preprocessing) ---
c --- beam_helkey=1 -- good;   beam_helkey=0 -- bad                  ---      
      IF(nbufmax_hel_ini.le.ihel_head) THEN
      DO i = nbufmax_hel_ini,ihel_head
         if(i.lt.ihel_head) then
            end_head=event_head(i+1)-1
         elseif(i.eq.ihel_head) then
            end_head=nevent_read
         endif
      IF(Nstory.ge.event_head(i).AND.Nstory.le.end_head) THEN
        beam_helkey = key(i)
        beam_keypair=key_pair(i)
        beam_keyboth=key_both(i)
        mult=st_mult(i)
        dose=dose_FC(i)
        if (dose.ne.0) then
           el_to_dose=mult/dose
        else
           el_to_dose=0
        endif
        ltime=tg_hevt
	nbufmax_hel_ini=i
        IF(Nstory-1.lt.event_head(i)) THEN
           fill311=.true.
           if(helicity.eq.2.OR.helicity.eq.3) then
              doplus_FC =doplus_FC +dose_FC(i)
              doplus_BPM=doplus_BPM+dose_BPM(i)
           else
              dominus_FC =dominus_FC +dose_FC(i)
              dominus_BPM=dominus_BPM+dose_BPM(i)
           endif
        ENDIF
        GOTO 1002
      ENDIF 
      ENDDO
      ENDIF
      DO i=1,ihel_head
         if(i.lt.ihel_head) then
            end_head=event_head(i+1)-1
         elseif(i.eq.ihel_head) then
            end_head=nevent_read
         endif 
      IF(Nstory.ge.event_head(i).and.Nstory.le.end_head)THEN
        beam_helkey = key(i)
        beam_keypair=key_pair(i)
        beam_keyboth=key_both(i)
        mult=st_mult(i)
        dose=dose_FC(i)
        if (dose.ne.0) then
           el_to_dose=mult/dose 
        else
           el_to_dose=0
        endif           
        ltime=tg_hevt
	nbufmax_hel_ini=i
        IF(Nstory-1.lt.event_head(i)) THEN
           fill311=.true. 
           if(helicity.eq.2.OR.helicity.eq.3) then
              doplus_FC =doplus_FC +dose_FC(i)
              doplus_BPM=doplus_BPM+dose_BPM(i)
           else
              dominus_FC =dominus_FC +dose_FC(i)
              dominus_BPM=dominus_BPM+dose_BPM(i)
           endif
        ENDIF
        GOTO 1002
      ENDIF 
      ENDDO
 1002 CONTINUE


c --- some variables initialisation ---
      if( MCbit.lt.0.or.MCbit.gt.1 ) print *,'  MCbitbad=',MCbit
      if( MCbit.lt.0.or.MCbit.gt.1 ) stop
      do i=1,8
        ChanDoneFlag(i,MCbit) = 0
      enddo
        max_pi = 3 !maximum number of pi to calculate different combinations of missing mass
        max_k  = 2 !maximum number of k to calculate different combinations of missing mass
        max_g  = 4 !maximum number of photons to calculate different combinations of invariant mass
        max_gg = 1 !maximum number of couples of two photons (max_gg=1 means we consider ONLY max_g=2
        ntot = 0
        bit = -1000
        bit1pi = -1000
        bit2pi  = -1000
        bitomega = -1000
        nNU = 0
	npi = 0
        npi_plus = 0
        npi_minus = 0
	nK  = 0
	nG  = 0
        nGG = 0
        nGGG = 0
        nGGGG = 0
	nl  = 0
        nd = 0
	neN   = 0       	
	nep  = 0         	
	nek  = 0      	
	neNp  = 0         	
	nepp  = 0       	
	neNk  = 0           	
	nepk  = 0           	
	nekk   = 0           	
	neNpp  = 0            	
	neppp    = 0          	
	neNkk   = 0           	
	neNpk   = 0           	
	neNppp    = 0
        neNG = 0
	ne0   = 0 
	neN0   = 0       	
	nep0  = 0         	
	nek0  = 0      	
	neNp0  = 0         	
	nepp0  = 0       	
	neNk0  = 0           	
	nepk0  = 0           	
	nekk0   = 0           	
	neNpp0  = 0            	
	neppp0    = 0          	
	neNkk0   = 0           	
	neNpk0   = 0           	
	neNppp0    = 0 
        ncmrest = 0
           W = -1000
           Q2 = -1000

c Some parameters
	n_ev_run = event_num
        trig_type = tg_hevt  !event_type
        trig_clas = event_clas
        fcgdose = fcg_hevt
        beam_hel = helicity
	n_ev_nt  = nevent
	n_part   = nevnt


c -------------------------- correct momenta ------------------------
       c=acos(-1.0)/180.0
       IF(mode_momcor) THEN

c --- momentum corrections for e1a ---
       IF(mode_momcortype.eq.1) THEN
         DO j = 1,nevnt
         IDP=0
         if(         j .eq. 1   ) IDP=11
         if(ID_HADR(j) .eq. 11  ) IDP=11
         if(ID_HADR(j) .eq. 2212) IDP=2212
         if(ID_HADR(j) .eq. 211)  IDP=211
         if(ID_HADR(j) .eq.-211)  IDP=-211
         if(IDP.ne.0) THEN
         CALL cormom_e1(IDP,pmom_EVNT(j),th_EVNT(j),ph_EVNT(j),Eelbeam,I_torus)
         Pmomx_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*cos(ph_EVNT(j)*c)
         Pmomy_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*sin(ph_EVNT(j)*c)
         Pmomz_EVNT(j) = pmom_EVNT(j)*cos(th_EVNT(j)*c)
         ENDIF
         ENDDO

c --- momentum corrections for eg1 ---
       ELSEIF(mode_momcortype.eq.2) THEN
         DO j = 1,nevnt
         IDP=0
         if(         j .eq. 1   ) IDP=11
         if(ID_HADR(j) .eq. 11  ) IDP=11
         if(ID_HADR(j) .eq. 2212) IDP=2212
         if(ID_HADR(j) .eq. 211)  IDP=211
         if(ID_HADR(j) .eq.-211)  IDP=-211
         IF(IDP.ne.0) THEN
         CALL cormom_eg1(IDP,pmom_EVNT(j),th_EVNT(j),ph_EVNT(j),Eelbeam,I_torus)
         Pmomx_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*cos(ph_EVNT(j)*c)
         Pmomy_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*sin(ph_EVNT(j)*c)
         Pmomz_EVNT(j) = pmom_EVNT(j)*cos(th_EVNT(j)*c)
         ENDIF
         ENDDO

c --- Burkert's corrections for eg1 ---
        ELSEIF(mode_momcortype.eq.3) THEN
         DO j = 1,nevnt
         ch1='?'
         if(j.eq.1) ch1='E'
         if(ID_HADR(j) .eq.  11 .or. ID_HADR(j) .eq.  -11 ) ch1='E'
         if(ID_HADR(j) .eq. 211 .or. ID_HADR(j) .eq. -211 ) ch1='H'
         if(ID_HADR(j) .eq. 321 .or. ID_HADR(j) .eq. -321 ) ch1='H'
         if(ID_HADR(j) .eq.2212                           ) ch1='H'
         if(ID_HADR(j) .eq.  45 .                         ) ch1='H'
         if(ID_HADR(j) .eq.  49 .or. ID_OUR(j) .eq.   47  ) ch1='H'
         IF(ch1.eq.'E'.or.ch1.eq.'H') THEN
           CALL momcorr_eg1(ch1,Eelbeam,I_torus,sector(ph_EVNT(j)),
     &          pmom_EVNT(j),th_EVNT(j),ph_EVNT(j),Charge_EVNT(j)*1.)
           Pmomx_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*cos(ph_EVNT(j)*c)
           Pmomy_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*sin(ph_EVNT(j)*c)
           Pmomz_EVNT(j) = pmom_EVNT(j)*cos(th_EVNT(j)*c)
         ENDIF
         ENDDO

c --- Burkert's corrections for e1 ---
c         DO j = 1,nevnt
c         ch1='?'
c         IF(j.eq.1) ch1='E'
c         IF(ID_HADR(j) .eq.  11 .or. ID_HADR(j) .eq.  -11 ) ch1='E'
c         IF(ID_HADR(j) .eq. 211 .or. ID_HADR(j) .eq. -211 ) ch1='H'
c         IF(ID_HADR(j) .eq. 321 .or. ID_HADR(j) .eq. -321 ) ch1='H'
c         IF(ID_HADR(j) .eq.2212                           ) ch1='H'
c         IF(ID_HADR(j) .eq.  45 .                         ) ch1='H'
c         IF(ID_HADR(j) .eq.  49 .or. ID_OUR(j) .eq.   47  ) ch1='H'
c         IF(ch1.eq.'E'.or.ch1.eq.'H') THEN
c           CALL momcorr_e1(ch1,Eelbeam,I_torus,sector(ph_EVNT(j)),
c     &          pmom_EVNT(j),th_EVNT(j),ph_EVNT(j),Charge_EVNT(j)*1.)
c           Pmomx_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*cos(ph_EVNT(j)*c)
c           Pmomy_EVNT(j) = pmom_EVNT(j)*sin(th_EVNT(j)*c)*sin(ph_EVNT(j)*c)
c           Pmomz_EVNT(j) = pmom_EVNT(j)*cos(th_EVNT(j)*c)
c         ENDIF
c         ENDDO


       ELSE
         print *,' bad mode_momcortype=',mode_momcortype
         stop
       ENDIF
       ENDIF


c Electron beam, Electron variables
        IF(beam_type)THEN
c        skip event with electron ID =0 
         IF(nEVNT.le.0.or.ID_EVNT(1).eq.0)THEN
           Select_OK = .false.
           RETURN
         ENDIF       
	 Omega = E0 - pmom_EVNT(1)
	 Q2    = 4.*E0*pmom_EVNT(1)*(sin(th_EVNT(1)*3.1416/2./180.))**2
	 W2_local = m_nucleon(1)**2+2.*m_nucleon(1)*omega-Q2
	 if(W2_local.ge.0)  W = sqrt(W2_local)
         th_EL = th_EVNT(1)
         ph_EL = ph_EVNT(1)
         E_EL = E0 - Omega
         P_EL = E_EL
         B_EL = beta_EVNT(1)
         px_EL = pmomx_EVNT(1)
         py_EL = pmomy_EVNT(1)
         pz_EL = pmomz_EVNT(1)
         x_EL = x_EVNT(1)
         y_EL = y_EVNT(1)
         z_EL = z_EVNT(1) - z_el_mid(sector(ph_EL))
         xMVRT_EL = -1000.
         yMVRT_EL = -1000.
         zMVRT_EL = -1000.
c no CC efficiency calculation now. It's not credible.
         if(input_type.eq.1.and.nDCPB.eq.1) DCstat_EVNT(1)=1
         if(DCstat_EVNT(1).gt.0) then 
           pxyz(1) = XSC_DCPB(DCstat_EVNT(1))
           pxyz(2) = YSC_DCPB(DCstat_EVNT(1))
           pxyz(3) = ZSC_DCPB(DCstat_EVNT(1))
           cxyz(1) = CXSC_DCPB(DCstat_EVNT(1))
           cxyz(2) = CYSC_DCPB(DCstat_EVNT(1))
           cxyz(3) = CZSC_DCPB(DCstat_EVNT(1))
         key_cc = cceffn( 16219, sector(ph_EL), pxyz, cxyz, 
     &                       nphe2, nphe2err, EffCC_EL )
           if(key_cc.ne.0) then
             EffCC_EL  = -1000.
*             print *,' key_cc=',key_cc
           endif
           if( .not. (EffCC_EL.ge.0. .and. EffCC_EL.le.1.0) ) then
           EffCC_EL  = -999.
           endif
         else
           EffCC_EL  = -1000.
         endif
         if(input_type.eq.1.and.nCCPB.eq.1) CCstat_EVNT(1)=1
         if(CCstat_EVNT(1).gt.0) then 
           NpheCC_EL = Nphe_CCPB(CCstat_EVNT(1))
         else
           NpheCC_EL = -1000.
         endif
         if(ECstat_EVNT(1).gt.0) then
           ECtot_EL = Etot_ECPB(ECstat_EVNT(1)) 
           ECin_EL  = Ein_ECPB (ECstat_EVNT(1)) 
           ECout_EL = Eout_ECPB(ECstat_EVNT(1))
         else
           ECtot_EL = -1000.
           ECin_EL  = -1000.
           ECout_EL = -1000.
         endif
         if(SCstat_EVNT(1).gt.0 .and. SCstat_EVNT(1).le.nSCPB) then
           PdHit_EL = ScPdHt_SCPB(SCstat_EVNT(1))
           PdHit_EL = (PdHit_EL-10000*int(PdHit_EL/10000))/100
         else
           PdHit_EL = -1000
         endif
         statEVNT_EL= Status_EVNT(1)
         trkEVNT_EL = Trk_flag_EVNT(1)

c acceptance for electron         
         call pseudo_spa2(0,11, p_EL,th_EL,(mod(ph_el+30.,60.)-30.),
     &                    I_torus,1.,sector(ph_EL),acc_EL)

c Some useful quantity related to electron
         eps   = (1.+2.*(1.+ Omega**2/Q2)*(TAN(th_EL*3.14159265/180./2.))**2)**(-1)
         eps_l = Q2*eps/Omega**2                                   
	 Eq_photon_E = (W**2-m_nucleon(1)**2)/2./m_nucleon(1) ! Equivalent photon Energy
c         Flux_vp = (1./137./2./PI**2)*(E_EL/E0)*(Eq_photon_E/Q2)/(1.-eps) 
         Flux_vp = 
     &   (1/137./4/3.1416/E0**2/m_nucleon(1)**2)*W*(W**2-m_nucleon(1)**2)/Q2/(1.-eps) ! Virtual Photon Flux

c Starting photon beam variables
       ELSEIF(n_ph_in.gt.0) THEN
          Egamma = ph_ERG_TAGR(1)
          E0 = Egamma ! to use following formula (made for electron's beam)
          Q2 = 0.
          Omega = Egamma
	  W2_local = m_nucleon(1)**2+2.*m_nucleon(1)*E0
	  if(W2_local.ge.0)  W = sqrt(W2_local)
          T_id = ph_T_id_TAGR(1)
          E_id = ph_E_id_TAGR(1)
          Tag_time =    ph_TPHO_TAGR (1)
          ST_time     = ph_STT_HEVT  (1)
          RF1_time    = ph_RF1_HEVT  (1)
          RF2_time    = ph_RF2_HEVT  (1)
          Hit_in_time = ph_hit_in_time(1)
          if(MCbit.eq.0) true = id_our(1) ! true = -1; random = 1; noinfo = 0
          if(MCbit.eq.1) then ! reading once and saving true & chan
           true =  -charge_EVNT(1) !Simulated photon data: true = -1 random =1
           chan =  Id_EVNT(1)      !Simulated photon data
          endif

          IF (event_type.eq.1
     &        .and. EgammaIsGood(Egamma) 
     &                                 ) THEN
c --- get Eid converting from energy 
          IF    (RunNumber.gt.12000.and.RunNumber.lt.13000)then !G6A
          do i=1,Etag_Nbins  ! get the number of E scint.
          if(Egamma.ge.Etag_emin(i)*Eelbeam.and.
     &       Egamma.le.Etag_emax(i)*Eelbeam)then
          if(NGammaE(i).gt.0.0) then             
            CALL HFILL(11007, i+0.5,0., 1.0/NGammaE(i))
            CALL HFILL(11009, i+0.5,0., 1.0/NGammaE(i))
            CALL HFILL(11010, i+0.5,0., 1.0)
            goto 101
          endif
          endif
          enddo
c --- get Eid from the ntuple
          ELSEIF(RunNumber.gt.19000.and.RunNumber.lt.20000)then !G6B
          i=0
          do j=1,clasg6_taghit
          if(abs(clasg6_E_gamma(j)- Egamma).lt.0.01) i=clasg6_Eid(j)
          enddo
          if(i.ne.0.and.NGammaE(i).gt.0.0)then
            CALL HFILL(11007, i+0.5,0., 1.0/NGammaE(i))
            CALL HFILL(11009, i+0.5,0., 1.0/NGammaE(i))
            CALL HFILL(11010, i+0.5,0., 1.0)
          endif
c ---
          ELSE
           print *,' select_goa: Eid: run period - ?  STOP.'
           stop
          ENDIF
          
 101      continue
          do i=1,Ttag_Nbins   ! get the number of T scint.
          if(Egamma.ge.Ttag_emin(i)*Eelbeam.and.
     &       Egamma.le.Ttag_emax(i)*Eelbeam)then
          if(NGammaT(i).gt.0.0) then
              CALL HFILL(11008, i+0.5,0., 1.0/NGammaT(i))
              goto 102
          endif
          endif
          enddo
 102      continue
          ENDIF

        ENDIF   

c+ Just electron detected (ep->eX):center of mass  variables
        Pmiss_e (0) = Omega + m_nucleon(1)       !for photon Omega = E0
        Pmiss_e (1) = - px_EL               !for photon px_EL  init @ 0.
        Pmiss_e (2) = - py_EL               !for photon py_EL  init @ 0.
        Pmiss_e (3) = E0 - pz_EL            !for photon pz_EL  init @ 0.
        Pmiss_e (5) = sqrt(Pmiss_e(1)**2 + Pmiss_e(2)**2 + Pmiss_e(3)**2) 
        Pmiss_e (4) = Pmiss_e (0)**2 - Pmiss_e (5)**2 
        call angles(0,Pmiss_e(1),Pmiss_e(2),Pmiss_e(3), Pmiss_e (6),Pmiss_e(7))        
c+ Charged particles in opposite sectors
c        if (mcbit.eq.1) call opposite_sec
         call opposite_sec(mcbit)
c-


c+ Loop over nEVNT to fill ntuple 10
	do j = 2,nevnt

c++ Nucleon
	  if(ID_hadr(j).eq.2212.or.ID_hadr(j).eq.2112) then	! proton or neutron
           ntot = ntot+1                  ! total particle number counter
           nNU=nNU+1
           if(nNU.gt.1) then
            W=-900.
            nNU = 1
            return ! Skipping events with more than 1 baryon
           endif
	   P_NU(1)  = Pmom_EVNT(j)
           B_NU(1)  = beta_hadr(j)
           M2_NU(1) = mass_HADR(j)
	   PX_NU(1) = Pmomx_EVNT(j)
	   PY_NU(1) = Pmomy_EVNT(j)
	   PZ_NU(1) = Pmomz_EVNT(j)
	   th_NU(1) = th_EVNT(j)
	   ph_NU(1) = ph_EVNT(j)
	   q_NU(1)  = Charge_EVNT(j) 
           E_NU(1)  = sqrt(Pmom_EVNT(j)**2 + m_nucleon(q_NU(1))**2)
           x_NU(1)  = x_EVNT(j)
           y_NU(1)  = y_EVNT(j)
           statEVNT_NU(1)= Status_EVNT(j)
           trkEVNT_NU(1) = Trk_flag_EVNT(j)
           thrSC_NU(1) = Thr_trigSC_EVNT(j)
             if(SCstat_EVNT(j).gt.0 .and. SCstat_EVNT(j).le.nSCPB) then 
           PdHit_NU(1) = ScPdHt_SCPB(SCstat_EVNT(j))
           PdHit_NU(1) = (PdHit_NU(1)-10000*int(PdHit_NU(1)/10000))/100
             else
           PdHit_NU(1) = -1000.
             endif
           id_SC_NU(1) = id_SC_EVNT(j)
           sec_SC_NU(1) = sec_SC_EVNT(j)
           IF(q_NU(1).gt.0.)z_NU(1)=z_EVNT(j)-z_p_mid(sector(ph_NU(1)))
           IF(q_NU(1).eq.0.)z_NU(1)=z_EVNT(j)
            if (beam_type) then 
            th1_NU(1) = E0/0.9383
            th1_NU(1) = atan(  1./(1.+th1_NU(1))/tan(th_EL*pi/180./2.)  )
            if(th1_NU(1).lt.0) th1_NU(1)=th1_NU(1)+pi
            th1_NU(1) = th_NU(1)-th1_NU(1)*180./pi
            th2_NU(1) = 1./( (0.9383/E_EL) -1. + cos(th_EL*pi/180.) )
            th2_NU(1) = atan(  1./(1.+th2_NU(1))/tan(th_EL*pi/180./2.)  )
            if(th2_NU(1).lt.0) th2_NU(1)=th2_NU(1)+pi
            th2_NU(1) = th_NU(1)-th2_NU(1)*180./pi
            endif
           if(q_NU(1).gt.0.)then
           call pseudo_spa2(+1,2212, p_NU(1),th_NU(1),
     &                     (mod(ph_NU(1)+30.,60.)-30.),I_torus,
     &                     1.,sector(ph_NU(1)),acc_NU(1))
           else
           acc_NU(1)=1.
           endif

c+++ Calling Eloss for photon data
         if(.not.beam_type.and.ID_hadr(j).eq.2212) then
          if (bit_mc.eq.0.or.(bit_mc.eq.2.and.MCbit.eq.0)) then
           ip = 1                    !  1=proton 2=pion
           ve(1) = 0.                !  x vertex
           ve(2) = 0.                !  y vertex
           if(abs(z_EVNT(j)).lt.8.5) then
            ve(3) = z_EVNT(j)        !  z vertex if z inside the 17cm long target
           else 
            ve(3) = z_EVNT(j) * 8.5 / abs(z_EVNT(j))     !  z vertex if z outside the target  
           endif    
           cd(1) = PX_NU(1)/P_NU(1)  !  px cosin
           cd(2) = PY_NU(1)/P_NU(1)  !  py cosin
           cd(3) = PZ_NU(1)/P_NU(1)  !  pz cosin
           pm    = P_NU(1)           !  p 
c          eloss = OUTPUT            !  Status flag: 1=ok 0=(dedx=0.)
c          pv = OUTPUT               ! = pcorrected = pmom + dedx
           status_eloss = eloss(ip,ve,cd,pm,pv)
           if(status_eloss.ne.0) then 
            p_nu(1) = pv
	    PX_NU(1)  = pv*cd(1)
	    Py_NU(1)  = pv*cd(2)
	    Pz_NU(1)  = pv*cd(3)
            E_NU(1)   = sqrt(P_nu(1)**2 + m_nucleon(q_NU(1))**2) 
            b_nu(1) = p_NU(1) / E_NU(1)
           endif
          endif
         endif

	         
c++ Deuterium
	  elseif(ID_hadr(j).eq.45) then	! deuteron
           IF(nd.lt.5) THEN
           ntot = ntot+1
	   nd = nd + 1
	   E_D (nd) = sqrt(Pmom_EVNT(j)**2 + 1.8756**2) 
	   P_D (nd) = Pmom_EVNT(j)
	   B_D (nd) = beta_hadr(j)
	   M2_D (nd) = mass_HADR(j)
	   PX_D(nd)  = Pmomx_EVNT(j)
	   PY_D(nd)  = Pmomy_EVNT(j)
	   PZ_D(nd)  = Pmomz_EVNT(j)
	   th_D(nd) = th_EVNT(j)
	   ph_D(nd) = ph_EVNT(j)
	   q_D (nd) = Charge_EVNT(j)	
           x_D(nd) = x_EVNT(j)
           y_D(nd) = y_EVNT(j)
           z_D(nd) = z_EVNT(j) - z_p_mid(sector(ph_D(nd)))
           statEVNT_D(nd)= Status_EVNT(j)
           trkEVNT_D(nd) = Trk_flag_EVNT(j)
            if(SCstat_EVNT(j).gt.0 .and. SCstat_EVNT(j).le.nSCPB) then 
           PdHit_D(nd) = ScPdHt_SCPB(SCstat_EVNT(j))
           PdHit_D(nd) = (PdHit_D(1)-10000*int(PdHit_D(1)/10000))/100
            else
           PdHit_D(nd) = -1000.
            endif
c WRONG ->>>      call pseudo_spa(
           ENDIF

	         
c++ Pions	  	   
	  elseif(ID_hadr(j).eq.211.OR.ID_hadr(j).eq.-211) then	! pion
           ntot = ntot+1
	   npi = npi + 1
	   E_PI (npi) = sqrt(Pmom_EVNT(j)**2 + .138**2) 
	   P_PI (npi) = Pmom_EVNT(j)
	   B_PI (npi) = beta_hadr(j)
	   M2_PI (npi) = mass_HADR(j)
	   PX_PI(npi)  = Pmomx_EVNT(j)
	   PY_PI(npi)  = Pmomy_EVNT(j)
	   PZ_PI(npi)  = Pmomz_EVNT(j)
	   th_PI(npi) = th_EVNT(j)
	   ph_PI(npi) = ph_EVNT(j)
	   q_PI (npi) = Charge_EVNT(j)	!  q=+1 pi+ q=-1 pi-
           x_PI(npi) = x_EVNT(j)
           y_PI(npi) = y_EVNT(j)
           statEVNT_PI(npi)= Status_EVNT(j)
           trkEVNT_PI(npi) = Trk_flag_EVNT(j)
           thrSC_PI(npi) = Thr_trigSC_EVNT(j)
           id_SC_PI(npi) = id_SC_EVNT(j)
           sec_SC_PI(npi) = sec_SC_EVNT(j)
           IF(q_PI(npi).le.0.) z_PI(npi)=z_EVNT(j)-z_el_mid(sector(ph_PI(npi)))
           IF(q_PI(npi).gt.0.) z_PI(npi)=z_EVNT(j)- z_p_mid(sector(ph_PI(npi)))
           ECtot_PI(npi) = -1000.
           ECin_PI(npi)  = -1000.
           ECout_PI(npi) = -1000.
           if(ECstat_EVNT(j).gt.0) then
            ECtot_PI(npi)= Etot_ECPB(ECstat_EVNT(j)) !Etot_ECPB(j)
            ECin_PI(npi) = Ein_ECPB (ECstat_EVNT(j)) !Ein_ECPB (j)
            ECout_PI(npi)= Eout_ECPB(ECstat_EVNT(j)) !Eout_ECPB(j)
           endif
           if(SCstat_EVNT(j).gt.0 .and. SCstat_EVNT(j).le.nSCPB) then 
            PdHit_PI(npi)= ScPdHt_SCPB(SCstat_EVNT(j))
            PdHit_PI(npi)= (PdHit_PI(npi)-10000*int(PdHit_PI(npi)/10000))/100
           else
            PdHit_PI(npi)= -1000.
           endif
              if(q_pi(npi).eq.1) then
              npi_plus = npi_plus+1   ! positive pions counter
              point_pp = npi          ! positive pions pointer
              endif
              if(q_pi(npi).eq.-1) then
              npi_minus = npi_minus+1 ! negative pions counter
              point_pm = npi          ! negative pions pointer
              endif
           call pseudo_spa2(int(q_PI(npi)),211*int(q_PI(npi)),p_PI(npi),
     &                     th_PI(npi),(mod(ph_PI(npi)+30.,60.)-30.),
     &                     I_torus,1.,sector(ph_PI(npi)),acc_PI(npi))

c++ Kaons
	  elseif(ID_hadr(j).eq.321.OR.ID_hadr(j).eq.-321) then	! Kaons
           ntot = ntot+1
	   nk = nk + 1
	   E_k (nk) = sqrt(Pmom_EVNT(j)**2 + .494**2) 
	   P_k (nk) = Pmom_EVNT(j)
	   B_k (nk) = beta_hadr(j)
	   M2_k (nk) = mass_HADR(j)
	   PX_K (nk)  = Pmomx_EVNT(j)
	   PY_K (nk)  = Pmomy_EVNT(j)
	   PZ_K (nk) = Pmomz_EVNT(j)
	   th_k(nk) = th_EVNT(j)
	   ph_k(nk) = ph_EVNT(j)
	   q_k (nk) = Charge_EVNT(j)	!  q=+1 k+ q=-1 k-
           x_K(nk) = x_EVNT(j)
           y_K(nk) = y_EVNT(j)
           statEVNT_K(nk)= Status_EVNT(j)
           trkEVNT_K(nk) = Trk_flag_EVNT(j)
           if(q_k(nk).lt.0) z_K(nk) = z_EVNT(j) - z_el_mid(sector(ph_k(nk)))
           if(q_k(nk).gt.0) z_K(nk) = z_EVNT(j) -  z_p_mid(sector(ph_k(nk)))
         call pseudo_spa2(int(q_K(nk)),321*int(q_K(nk)),p_K(nk),
     &                   th_K(nk),(mod(ph_K(nk)+30.,60.)-30.),
     &                   I_torus,1.,sector(ph_K(nk)),acc_K(nk))

c++ Gammas
	  elseif(ID_hadr(j).eq.22) then	! gammas
           ntot = ntot+1
	   ng = ng + 1
	   E_g (ng) =Pmom_EVNT(j)
	   PX_G (ng)  = Pmomx_EVNT(j)
	   PY_G (ng)  = Pmomy_EVNT(j)
	   PZ_G (ng)  = Pmomz_EVNT(j)
	   th_g(ng) = th_EVNT(j)
	   ph_g(ng) = ph_EVNT(j)
           x_G(ng) = x_EVNT(j)
           y_G(ng) = y_EVNT(j)
           z_G(ng) = z_EVNT(j)
	   acc_g (ng) = 1.
           if(SCstat_EVNT(j).gt.0 .and. SCstat_EVNT(j).le.nSCPB) then 
           PdHit_G(ng)= ScPdHt_SCPB(SCstat_EVNT(j))
           PdHit_G(ng)= (PdHit_G(ng)-10000*int(PdHit_G(ng)/10000))/100
           else
           PdHit_G(ng)= -1000.
           endif


c++ Leptons (others)
        else
           ntot = ntot+1
	   nl = nl + 1
             if(nl.gt.4) return   !  rejects events with more than 4 "leptons" or "unknowns"
	   E_l (nl) = -1000. 
	   P_l (nl) = Pmom_EVNT(j)
           B_l (nl) = beta_hadr(j)
	   PX_L (nl)  = Pmomx_EVNT(j)
	   PY_L (nl)  = Pmomy_EVNT(j)
	   PZ_L (nl)  = Pmomz_EVNT(j)
	   th_l(nl) = th_EVNT(j)
	   ph_l(nl) = ph_EVNT(j)
	   q_l (nl) = Charge_EVNT(j)	
	   id_l(nl) = id_hadr(j)
           x_L(nl) = x_EVNT(j)
           y_L(nl) = y_EVNT(j)
           if(q_l(nl).lt.0.)z_L(nl)=z_EVNT(j)-z_el_mid(sector(ph_L(nl)))
           if(q_l(nl).gt.0.)z_L(nl)=z_EVNT(j)- z_p_mid(sector(ph_L(nl)))
           call pseudo_spa2(int(q_l(nl)),0,
     &                     p_l(nl),th_l(nl),(mod(ph_l(nl)+30.,60.)-30.),
     &                     I_torus,1.,sector(ph_l(nl)),acc_l(nl))
        endif	   
 
	enddo ! loop over particles in 1 event

c+ in case of two pions detected ordering them
c+ 1 = pi- ; 2 = pi+
        if(npi.eq.2.and.q_pi(2).ne.1) then
           spool(1) = E_PI(1)
           spool(2) =  P_PI(1)
           spool(3) =  B_PI(1)
           spool(4) = M2_PI(1)
           spool(5) = PX_PI(1)
           spool(6) = PY_PI(1)
           spool(7) = PZ_PI(1) 
           spool(8) = th_PI(1)
           spool(9) = ph_PI(1)
           spool(10) =  q_PI(1)
           spool(11) =  acc_PI(1)
           spool(12) = x_PI(1)
           spool(13) = y_PI(1)
           spool(14) = z_PI(1)
           spool(15) = ECtot_PI(1)
           spool(16) = ECin_PI(1)
           spool(17) = ECout_PI(1)
           Ipool(18) = StatEVNT_PI(1)
           Ipool(19) = TrkEVNT_PI(1)
           spool(20) = ThrSC_PI(1)
           Ipool(21) = PdHit_PI(1)
           spool(22) = id_SC_PI(1)
           spool(23) = sec_SC_PI(1)
  
            E_PI(1)= E_PI(2)
            P_PI(1)= P_PI(2)
            B_PI(1)= B_PI(2)
            M2_PI(1)=M2_PI(2)
            PX_PI(1)=PX_PI(2)
            PY_PI(1)=PY_PI(2)
            PZ_PI(1)=PZ_PI(2)
            th_PI(1)=th_PI(2)
            ph_PI(1)=ph_PI(2)
            q_PI(1)= q_PI(2)
            acc_PI(1)= acc_PI(2)
            x_PI(1)=x_PI(2)
            y_PI(1)=y_PI(2)
            z_PI(1)=z_PI(2)
            ECtot_PI(1)=ECtot_PI(2)
            ECin_PI(1)=ECin_PI(2)
            ECout_PI(1)=ECout_PI(2)
            StatEVNT_PI(1) =StatEVNT_PI(2)
            TrkEVNT_PI(1)  =TrkEVNT_PI(2)
            ThrSC_PI(1)=ThrSC_PI(2)
            PdHit_PI(1)=PdHit_PI(2)
            id_SC_PI(1) =id_SC_PI(2)
            sec_SC_PI(1)=id_SC_PI(2)

            E_PI(2)  =spool(1)
            P_PI(2)  =spool(2)
            B_PI(2)  =spool(3)
            M2_PI(2) =spool(4)
            PX_PI(2) =spool(5)
            PY_PI(2) =spool(6)
            PZ_PI(2) =spool(7) 
            th_PI(2) =spool(8)
            ph_PI(2) =spool(9)
            q_PI(2)  =spool(10)
            acc_PI(2)=spool(11)
            x_PI(2)  =spool(12)
            y_PI(2)  =spool(13)
            z_PI(2)  =spool(14)
            ECtot_PI(2)=spool(15)
            ECin_PI(2) =spool(16)
            ECout_PI(2)=spool(17)
            StatEVNT_PI(2) =Ipool(18)
            TrkEVNT_PI(2)  =Ipool(19)
            ThrSC_PI(2)    =spool(20)
            PdHit_PI(2)    =Ipool(21)
            id_SC_PI(2)    =spool(22)
            sec_SC_PI(2)   =spool(23)

           endif



c+ Multiples gammas combinations (up to 4 gammas)
c++ G G
        if(ng.ge.2.and.ng.le.max_g) then
           jj = 0
         do j = 1, ng-1
            do i = j+1, ng
               jj = jj + 1
               nGG = jj
               call invariant_mass_2part
     %                     (E_G(j),PX_G(j),PY_G(j),PZ_G(j)
     %                     ,E_G(i),PX_G(i),PY_G(i),PZ_G(i)
     %                     ,W_GG(0,jj),W_GG(1,jj),W_GG(2,jj),W_GG(3,jj)
     %                     ,W_GG(4,jj),W_GG(5,jj),W_GG(6,jj),W_GG(7,jj))
               W2_GG(JJ) = W_GG(4,jj)
               E_GG (JJ) = W_GG(0,jj) 
               P_GG (JJ) = W_GG(5,jj)
               th_GG (JJ) = W_GG(6,jj)
               ph_GG (JJ) = W_GG(7,jj)
            enddo
         enddo
        endif

c++ G G G
        if(ng.ge.3.and.ng.le.max_g) then
         
           jj = 0
         do j = 1, ng-2
            do i = j+1, ng-1
               do k = i+1,ng 
               jj = jj + 1
               nGGG = jj    
               call invariant_mass_3part
     %                     (E_G(k),PX_G(k),PY_G(k),PZ_G(k)
     %                     ,E_G(j),PX_G(j),PY_G(j),PZ_G(j)
     %                     ,E_G(i),PX_G(i),PY_G(i),PZ_G(i)
     %                     ,W_GGG(0,jj),W_GGG(1,jj),W_GGG(2,jj),W_GGG(3,jj)
     %                     ,W_GGG(4,jj),W_GGG(5,jj),W_GGG(6,jj),W_GGG(7,jj))
               W2_GGG(JJ) = W_GGG(4,jj)
              enddo
            enddo
         enddo                 
        endif   

c++ G G G G 
        if(ng.ge.4.and.ng.le.max_g) then
         nGGGG = 1
         call invariant_mass_4part
     %        (E_G(1),PX_G(1),PY_G(1),PZ_G(1)
     %        ,E_G(2),PX_G(2),PY_G(2),PZ_G(2)
     %        ,E_G(3),PX_G(3),PY_G(3),PZ_G(3)
     %        ,E_G(4),PX_G(4),PY_G(4),PZ_G(4)
     %        ,W_GGGG(0,1),W_GGGG(1,1),W_GGGG(2,1),W_GGGG(3,1)
     %        ,W_GGGG(4,1),W_GGGG(5,1),W_GGGG(6,1),W_GGGG(7,1))
         W2_GGGG(1) = W_GGGG(4,1)
        endif   


c+ MonteCarlo corrections and forcing
c++ pi0 case (only 1pi0)
	do j = 2,nevnt   
         if(ID_hadr(j).eq.111) then ! it can be only in MC files
               nG= nG + 2 
               nGG = nGG + 1
               ntot=ntot + 2 - 1 ! pi0 is equivalent to 2 gammas - 1 'lepton'
               W_GG(0,1) = sqrt(Pmom_EVNT(j)**2 + .13497**2)
               W_GG(1,1) = Pmomx_EVNT(j)
               W_GG(2,1) = Pmomy_EVNT(j)
               W_GG(3,1) = Pmomz_EVNT(j)
               W_GG(4,1) = mass_evnt(j)
               W_GG(5,1) = Pmom_EVNT(j)
               W_GG(6,1) = th_EVNT(j)
               W_GG(7,1) = ph_EVNT(j)

               W2_GG (1) = W_GG(4,1)
               E_GG  (1) = W_GG(0,1) 
               P_GG  (1) = W_GG(5,1)
               th_GG (1) = W_GG(6,1)
               ph_GG (1) = W_GG(7,1)
          endif
         enddo




c Missing Mass calculation
c Missing mass for GG detected
c Missing and invariant mass for (GG) detected
          if (ng.eq.2) then
           ne0 = 1
           call miss_mass_2part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_e0(0),PMiss_e0(1),PMiss_e0(2),PMiss_e0(3)
     %                        ,PMiss_e0(4),PMiss_e0(5),PMiss_e0(6),PMiss_e0(7))
           M2_e0(1) = PMiss_e0(4) 
          endif   

c++ Missing mass for P detected
        if(nNU.eq.1) then
         neN = 1
         call miss_mass_2part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eN(0),PMiss_eN(1),PMiss_eN(2),PMiss_eN(3)
     %                        ,PMiss_eN(4),PMiss_eN(5),PMiss_eN(6),PMiss_eN(7))
         M2_eN (1) = PMiss_eN(4)
          E_eN (1) = PMiss_eN(0)
          P_eN (1) = PMiss_eN(5)
         th_eN (1) = PMiss_eN(6)
         ph_eN (1) = PMiss_eN(7)
c         q_eN  (1) = 

c Missing mass for P G: at maximum max_G(4) combinations (respect to 
c different gammas detected) are stored in different entries of M2_eNG vector 
         if((nG).ge.1.and.(nG).le.max_G) then
          neNG = nG
          j_l = 0
          do j=1,nG
c Missing mass
           call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_G(j),PX_G(j),PY_G(j),PZ_G(j)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNG(0,j),PMiss_eNG(1,j),PMiss_eNG(2,j),PMiss_eNG(3,j)
     %                        ,PMiss_eNG(4,j),PMiss_eNG(5,j),PMiss_eNG(6,j),PMiss_eNG(7,j))
           M2_eNG(j) = PMiss_eNG(4,j)

c+++++ Invariant mass
           call invariant_mass_2part
     %                     (E_G(j),PX_G(j),PY_G(j),PZ_G(j)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_NG(0,j),W_NG(1,j),W_NG(2,j),W_NG(3,j)
     %                     ,W_NG(4,j),W_NG(5,j),W_NG(6,j),W_NG(7,j))
           W2_NG(j) = W_NG(4,j)
          enddo
         endif

c++++ Missing and invariant mass for proton (GG) detected
          if (ng.eq.2) then
           neN0 = 1
           call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eN0(0),PMiss_eN0(1),PMiss_eN0(2),PMiss_eN0(3)
     %                        ,PMiss_eN0(4),PMiss_eN0(5),PMiss_eN0(6),PMiss_eN0(7))
           M2_eN0(1) = PMiss_eN0(4) 
           call invariant_mass_2part
     %                        (E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,W_N0(0),W_N0(1),W_N0(2),W_N0(3)
     %                        ,W_N0(4),W_N0(5),W_N0(6),W_N0(7))
           W2_N0(1) = W_N0(4)
             
          endif   
             
c Missing mass for P pi detected: at maximum max_pi(3) combinations (respect to 
c different pions detected) are stored in different entries of M2_eNp vector 
         if((npi_minus+npi_plus).ge.1.and.(npi_minus+npi_plus).le.max_pi) then
          neNp = (npi_minus+npi_plus)
          j_l = 0
          do j=1,(npi_minus+npi_plus)
c Missing mass
           call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNp(0,j),PMiss_eNp(1,j),PMiss_eNp(2,j),PMiss_eNp(3,j)
     %                        ,PMiss_eNp(4,j),PMiss_eNp(5,j),PMiss_eNp(6,j),PMiss_eNp(7,j))
     
          
c	   if((bit2pi.eq.10.or.bit2pi.eq.20).and.W.gt.0) then
	   
c	   if((W.gt.0.93).AND.(W.lt.0.95).AND.(ECout_EL/P_EL.gt.0.25-1.39*ECin_EL/P_EL)
c     %	    .AND.(NpheCC_EL.gt.25)) then
c            write(*,*) 'CMS!!!!!!', 'E_EL=', E_EL,  'PX_EL=', PX_EL, 'PY_EL=', PY_EL,
c     %	                'PZ_EL', PZ_EL, 'Q2=', Q2,
c     %                   'PMiss_eNp(0)=', PMiss_eNp(0,j), 'PMiss_eNp(1)=', PMiss_eNp(1,j),
c     %                   'PMiss_eNp(2)=', PMiss_eNp(2,j), 'PMiss_eNp(3)=', PMiss_eNp(3,j),
c     %                   'PMiss_eNp(4)=', PMiss_eNp(4,j), 'PMiss_eNp(5)=', PMiss_eNp(5,j),
c     %                   'PMiss_eNp(6)=', PMiss_eNp(6,j), 'PMiss_eNp(7)=', PMiss_eNp(7,j),
c     %                    'E_PI=', E_PI(j), 'PX_PI=', PX_PI(j), 'PY_PI=', PY_PI(j),
c     %                    'PZ_PI=', PZ_PI(j), 
c     %                   'E_NU=', E_NU(1), 'PX_NU(1)=', PX_NU(1), 'PY_NU(1)=', PY_NU(1),
c     %                    'PZ_NU(1)=', PZ_NU(1)
c     %                   
                          
c			  endif
c			  endif
			  

           M2_eNp(j) = PMiss_eNp(4,j)

c Invariant mass
           call invariant_mass_2part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Np(0,j),W_Np(1,j),W_Np(2,j),W_Np(3,j)
     %                     ,W_Np(4,j),W_Np(5,j),W_Np(6,j),W_Np(7,j))
           W2_Np(j) = W_Np(4,j)

c Missing and invariant mass for proton pi (GG) detected
            if (ng.eq.2) then
               neNp0 = neNp0+1
               j_l = j_l + 1
           call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eNp0(0,j_l),PMiss_eNp0(1,j_l),PMiss_eNp0(2,j_l),PMiss_eNp0(3,j_l)
     %                        ,PMiss_eNp0(4,j_l),PMiss_eNp0(5,j_l),PMiss_eNp0(6,j_l),PMiss_eNp0(7,j_l))
           M2_eNp0(j_l) = PMiss_eNp0(4,j_l)
           call invariant_mass_3part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_Np0(0,j_l),W_Np0(1,j_l),W_Np0(2,j_l),W_Np0(3,j_l)
     %                     ,W_Np0(4,j_l),W_Np0(5,j_l),W_Np0(6,j_l),W_Np0(7,j_l))
           W2_Np0(j_l) = W_Np0(4,j_l)
            endif   

          enddo

c+++++ Missing mass for P pi pi detected: possible combinations of max max_pi=3 pions
c+++++ detected are stored in different entries of M2_eNpp vector
          if((npi_minus+npi_plus).ge.2.and.(npi_minus+npi_plus).le.max_pi) then
           jj = 0
           j_l= 0
           do j=1,(npi_minus+npi_plus)-1
            do i=j+1,(npi_minus+npi_plus)
            jj = jj + 1
            neNpp = jj
c++++++ Missing Mass
            call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNpp(0,jj),PMiss_eNpp(1,jj),PMiss_eNpp(2,jj),PMiss_eNpp(3,jj)
     %                        ,PMiss_eNpp(4,jj),PMiss_eNpp(5,jj),PMiss_eNpp(6,jj),PMiss_eNpp(7,jj))
            M2_eNpp(jj) = PMiss_eNpp(4,jj)
            E_eNpp(jj)  = PMiss_eNpp(0,jj)

c++++++ Invariant mass
           call invariant_mass_3part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Npp(0,jj),W_Npp(1,jj),W_Npp(2,jj),W_Npp(3,jj)
     %                     ,W_Npp(4,jj),W_Npp(5,jj),W_Npp(6,jj),W_Npp(7,jj))
           W2_Npp(jj) = W_Npp(4,jj)

c++++ Missing and invariant mass for proton pi pi (GG) detected
          if (ng.eq.2) then
            neNpp0 = neNpp0 +1
            j_l = j_l + 1
            call miss_mass_5part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eNpp0(0,j_l),PMiss_eNpp0(1,j_l),PMiss_eNpp0(2,j_l),PMiss_eNpp0(3,j_l)
     %                        ,PMiss_eNpp0(4,j_l),PMiss_eNpp0(5,j_l),PMiss_eNpp0(6,j_l),PMiss_eNpp0(7,j_l))
            M2_eNpp0(j_l) = PMiss_eNpp0(4,j_l) 
           call invariant_mass_4part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_Npp0(0,j_l),W_Npp0(1,j_l),W_Npp0(2,j_l),W_Npp0(3,j_l)
     %                     ,W_Npp0(4,j_l),W_Npp0(5,j_l),W_Npp0(6,j_l),W_Npp0(7,j_l))             
     	   W2_Npp0(j_l) = W_Npp0(4,j_l)
          endif   


            enddo
           enddo

c++++++ Missing mass for P pi pi pi detected: only 1 entry for max_pi=3
           if((npi_minus+npi_plus).ge.3.and.(npi_minus+npi_plus).le.max_pi) then
            neNppp = 1
c+++++++ Missing Mass
            call miss_mass_5part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNppp(0,1),PMiss_eNppp(1,1),PMiss_eNppp(2,1),PMiss_eNppp(3,1)
     %                        ,PMiss_eNppp(4,1),PMiss_eNppp(5,1),PMiss_eNppp(6,1),PMiss_eNppp(7,1))
            M2_eNppp(1) = PMiss_eNppp(4,1)

c+++++++ Invariant mass
            call invariant_mass_4part
     %                        (E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Nppp(0,1),W_Nppp(1,1),W_Nppp(2,1),W_Nppp(3,1)
     %                     ,W_Nppp(4,1),W_Nppp(5,1),W_Nppp(6,1),W_Nppp(7,1))
            W2_Nppp(1) = W_Nppp(4,1)

c++++ Missing and invariant mass for proton pi pi pi (GG) detected
          if (ng.eq.2) then
            neNppp0 = 1
            call miss_mass_6part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNppp0(0,1),PMiss_eNppp0(1,1),PMiss_eNppp0(2,1),PMiss_eNppp0(3,1)
     %                        ,PMiss_eNppp0(4,1),PMiss_eNppp0(5,1),PMiss_eNppp0(6,1),PMiss_eNppp0(7,1))
            M2_eNppp0(1) = PMiss_eNppp0(4,1) 
            call invariant_mass_5part
     %                        (E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Nppp0(0,1),W_Nppp0(1,1),W_Nppp0(2,1),W_Nppp0(3,1)
     %                     ,W_Nppp0(4,1),W_Nppp0(5,1),W_Nppp0(6,1),W_Nppp0(7,1))
            W2_Nppp0(1) = W_Nppp0(4,1)              
          endif   
           endif
          endif
         endif


c++++  Missing mass for P k: at maximum max_k(2) combinations (respect to 
c++++  different kaons detected) are stored in different entries of M2_eNk vector 
         if((nk).ge.1.and.(nk).le.max_k) then
          neNk = nk
          j_l = 0
          do j=1,nk
c+++++ Missing mass
           call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNk(0,j),PMiss_eNk(1,j),PMiss_eNk(2,j),PMiss_eNk(3,j)
     %                        ,PMiss_eNk(4,j),PMiss_eNk(5,j),PMiss_eNk(6,j),PMiss_eNk(7,j))
           M2_eNk(j) = PMiss_eNk(4,j)

c+++++ Invariant mass
           call invariant_mass_2part
     %                     (E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Nk(0,j),W_Nk(1,j),W_Nk(2,j),W_Nk(3,j)
     %                     ,W_Nk(4,j),W_Nk(5,j),W_Nk(6,j),W_Nk(7,j))
           W2_Nk(j) = W_Nk(4,j)

c++++ Missing and invariant mass for proton k (GG) detected
          if (ng.eq.2) then
          neNk0 = neNk0 + 1
          j_l = j_l +1
           call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eNk0(0,j_l),PMiss_eNk0(1,j_l),PMiss_eNk0(2,j_l),PMiss_eNk0(3,j_l)
     %                        ,PMiss_eNk0(4,j_l),PMiss_eNk0(5,j_l),PMiss_eNk0(6,j_l),PMiss_eNk0(7,j_l))
     	   M2_eNk0(j_l) = PMiss_eNk0(4,j_l) 
           call invariant_mass_3part
     %                     (E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_Nk0(0,j_l),W_Nk0(1,j_l),W_Nk0(2,j_l),W_Nk0(3,j_l)
     %                     ,W_Nk0(4,j_l),W_Nk0(5,j_l),W_Nk0(6,j_l),W_Nk0(7,j_l))
     	  W2_Nk0(j_l) = W_Nk0(4,j_l)
          endif  

          enddo

c+++++  Missing mass for P k k: only 1 entry for max_k=2
          if(nk.ge.2.and.nk.le.max_k) then
           neNkk = 1

c+++++ Missing mass
           call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNkk(0,1),PMiss_eNkk(1,1),PMiss_eNkk(2,1),PMiss_eNkk(3,1)
     %                        ,PMiss_eNkk(4,1),PMiss_eNkk(5,1),PMiss_eNkk(6,1),PMiss_eNkk(7,1))
           M2_eNkk(1) = PMiss_eNkk(4,1)

c++++++ Invariant mass
           call invariant_mass_3part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Nkk(0,1),W_Nkk(1,1),W_Nkk(2,1),W_Nkk(3,1)
     %                     ,W_Nkk(4,1),W_Nkk(5,1),W_Nkk(6,1),W_Nkk(7,1))
           W2_Nkk(1) = W_Nkk(4,1)

c++++ Missing and invariant mass for proton K K (GG) detected
          if (ng.eq.2) then
           neNkk0 = 1
           call miss_mass_5part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eNkk0(0,1),PMiss_eNkk0(1,1),PMiss_eNkk0(2,1),PMiss_eNkk0(3,1)
     %                        ,PMiss_eNkk0(4,1),PMiss_eNkk0(5,1),PMiss_eNkk0(6,1),PMiss_eNkk0(7,1))
     	   M2_eNkk0(1) = PMiss_eNkk0(4,1) 
           call invariant_mass_4part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_Nkk0(0,1),W_Nkk0(1,1),W_Nkk0(2,1),W_Nkk0(3,1)
     %                     ,W_Nkk0(4,1),W_Nkk0(5,1),W_Nkk0(6,1),W_Nkk0(7,1))
     	   W2_Nkk0(1) = W_Nkk0(4,1)
          endif   
          endif
         endif

c++++ Missing mass for p + K + pi detected
         if(nk.eq.1.and.(npi_minus+npi_plus).eq.1) then
           neNpk = 1
c++++ Missing Mass
           call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNpk(0,1),PMiss_eNpk(1,1),PMiss_eNpk(2,1),PMiss_eNpk(3,1)
     %                        ,PMiss_eNpk(4,1),PMiss_eNpk(5,1),PMiss_eNpk(6,1),PMiss_eNpk(7,1))
           M2_eNpk(1) = PMiss_eNpk(4,1)
c++++++ Invariant mass
           call invariant_mass_3part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Npk(0,1),W_Npk(1,1),W_Npk(2,1),W_Npk(3,1)
     %                     ,W_Npk(4,1),W_Npk(5,1),W_Npk(6,1),W_Npk(7,1))
           W2_Npk(1) = W_Npk(4,1)
c++++ Missing and invariant mass for proton K PI (GG) detected
          if (ng.eq.2) then
           neNpk0 = 1
           call miss_mass_5part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                        ,PMiss_eNpk0(0,1),PMiss_eNpk0(1,1),PMiss_eNpk0(2,1),PMiss_eNpk0(3,1)
     %                        ,PMiss_eNpk0(4,1),PMiss_eNpk0(5,1),PMiss_eNpk0(6,1),PMiss_eNpk0(7,1))
     	M2_eNpk0(1) = PMiss_eNpk0(4,1) 
           call invariant_mass_4part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,E_NU(1), PX_NU(1),PY_NU(1),PZ_NU(1)
     %                     ,W_Npk0(0,1),W_Npk0(1,1),W_Npk0(2,1),W_Npk0(3,1)
     %                     ,W_Npk0(4,1),W_Npk0(5,1),W_Npk0(6,1),W_Npk0(7,1))
     	W2_Npk0(1) = W_Npk0(4,1)
          endif   
         endif           
        endif

c+++ Missing mass for p + k detected
        if(nk.eq.1.and.(npi_minus+npi_plus).eq.1) then
         nepk = 1
c++++ Missing Mass
         call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,PMiss_epk(0,1),PMiss_epk(1,1),PMiss_epk(2,1),PMiss_epk(3,1)
     %                        ,PMiss_epk(4,1),PMiss_epk(5,1),PMiss_epk(6,1),PMiss_epk(7,1))
         M2_epk(1) = PMiss_epk(4,1)
c++++++ Invariant mass
           call invariant_mass_2part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,W_pk(0,1),W_pk(1,1),W_pk(2,1),W_pk(3,1)
     %                     ,W_pk(4,1),W_pk(5,1),W_pk(6,1),W_pk(7,1))
           W2_pk(1) = W_pk(4,1)

c++++ Missing and invariant mass for pi k (GG) detected
          if (ng.eq.2) then
         nepk0 = 1             
         call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,PMiss_epk0(0,1),PMiss_epk0(1,1),PMiss_epk0(2,1),PMiss_epk0(3,1)
     %                        ,PMiss_epk0(4,1),PMiss_epk0(5,1),PMiss_epk0(6,1),PMiss_epk0(7,1))
     	   M2_epk0(1) = PMiss_epk0(4,1) 
           call invariant_mass_3part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,W_pk0(0,1),W_pk0(1,1),W_pk0(2,1),W_pk0(3,1)
     %                     ,W_pk0(4,1),W_pk0(5,1),W_pk0(6,1),W_pk0(7,1))
     	  W2_pk0(1) = W_pk0(4,1)
          endif   
       endif           

c+++ Missing mass for pi detected:different entries for different pions detected
        if((npi_minus+npi_plus).ge.1.and.(npi_minus+npi_plus).le.max_pi) then
         nep =(npi_minus+npi_plus)
         j_l = 0
         do j = 1,(npi_minus+npi_plus)
          call miss_mass_2part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,PMiss_ep(0,j),PMiss_ep(1,j),PMiss_ep(2,j),PMiss_ep(3,j)
     %                        ,PMiss_ep(4,j),PMiss_ep(5,j),PMiss_ep(6,j),PMiss_ep(7,j))
          M2_ep(j) = PMiss_ep(4,j)

c++++ Missing and invariant mass for PI  (GG) detected
          if (ng.eq.2) then
         nep0 = nep0 +1  
         j_l = j_l +1
          call invariant_mass_2part
     %                        (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,W_p0(0,j_l),W_p0(1,j_l),W_p0(2,j_l),W_p0(3,j_l)
     %                        ,W_p0(4,j_l),W_p0(5,j_l),W_p0(6,j_l),W_p0(7,j_l))
          W2_p0(j_l) = W_p0(4,j_l)     
          call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_ep0(0,j_l),PMiss_ep0(1,j_l),PMiss_ep0(2,j_l),PMiss_ep0(3,j_l)
     %                        ,PMiss_ep0(4,j_l),PMiss_ep0(5,j_l),PMiss_ep0(6,j_l),PMiss_ep0(7,j_l))
	  M2_ep0(j_l) = PMiss_ep0(4,j_l) 
          endif   

         enddo

c++++ Missing mass for pi pi detected: different entries for different combinations
c++++ of detected pions
         if((npi_minus+npi_plus).ge.2.and.(npi_minus+npi_plus).le.max_pi) then
          jj = 0
          j_l = 0
          do j=1,(npi_minus+npi_plus)-1
           do i=j+1,(npi_minus+npi_plus)
            jj = jj + 1
            nepp = jj
c+++++ Missing mass
            call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                        ,PMiss_epp(0,jj),PMiss_epp(1,jj),PMiss_epp(2,jj),PMiss_epp(3,jj)
     %                        ,PMiss_epp(4,jj),PMiss_epp(5,jj),PMiss_epp(6,jj),PMiss_epp(7,jj))
            M2_epp(jj) = PMiss_epp(4,jj)
c++++++ Invariant mass
           call invariant_mass_2part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                     ,W_pp(0,jj),W_pp(1,jj),W_pp(2,jj),W_pp(3,jj)
     %                     ,W_pp(4,jj),W_pp(5,jj),W_pp(6,jj),W_pp(7,jj))
           W2_pp(jj) = W_pp(4,jj)

c++++ Missing and invariant mass for pi pi  (GG) detected
          if (ng.eq.2) then
            nepp0 =  nepp0 +1
          j_l = j_l +1             
            call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                        ,PMiss_epp0(0,j_l),PMiss_epp0(1,j_l),PMiss_epp0(2,j_l),PMiss_epp0(3,j_l)
     %                        ,PMiss_epp0(4,j_l),PMiss_epp0(5,j_l),PMiss_epp0(6,j_l),PMiss_epp0(7,j_l))
     	   M2_epp0(j_l) = PMiss_epp0(4,j_l) 
           call invariant_mass_3part
     %                     (E_PI(j),PX_PI(j),PY_PI(j),PZ_PI(j)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,E_PI(i),PX_PI(i),PY_PI(i),PZ_PI(i)
     %                     ,W_pp0(0,j_l),W_pp0(1,j_l),W_pp0(2,j_l),W_pp0(3,j_l)
     %                     ,W_pp0(4,j_l),W_pp0(5,j_l),W_pp0(6,j_l),W_pp0(7,j_l))
     	   W2_pp0(j_l) = W_pp0(4,j_l)
          endif   
           enddo
          enddo

c+++++ Missing mass for pi pi pi detected: only 1 entry for max_pi=3
          if((npi_minus+npi_plus).ge.3.and.(npi_minus+npi_plus).le.max_pi) then
           neppp = 1
c+++++ Missing Mass
           call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,PMiss_eppp(0,1),PMiss_eppp(1,1),PMiss_eppp(2,1),PMiss_eppp(3,1)
     %                        ,PMiss_eppp(4,1),PMiss_eppp(5,1),PMiss_eppp(6,1),PMiss_eppp(7,1))
           M2_eppp(1) = PMiss_eppp(4,1)
c++++++ Invariant mass
           call invariant_mass_3part
     %                     (E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                     ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                     ,W_ppp(0,1),W_ppp(1,1),W_ppp(2,1),W_ppp(3,1)
     %                     ,W_ppp(4,1),W_ppp(5,1),W_ppp(6,1),W_ppp(7,1))
           W2_ppp(1) = W_ppp(4,1)

c++++ Missing and invariant mass for pi pi pi  (GG) detected
          if (ng.eq.2) then
           neppp0 = 1
           call miss_mass_5part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                        ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                        ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_eppp0(0,1),PMiss_eppp0(1,1),PMiss_eppp0(2,1),PMiss_eppp0(3,1)
     %                        ,PMiss_eppp0(4,1),PMiss_eppp0(5,1),PMiss_eppp0(6,1),PMiss_eppp0(7,1))
     	   M2_eppp0(1) = PMiss_eppp0(4,1) 
           call invariant_mass_4part
     %                     (E_PI(1),PX_PI(1),PY_PI(1),PZ_PI(1)
     %                     ,E_PI(2),PX_PI(2),PY_PI(2),PZ_PI(2)
     %                     ,E_PI(3),PX_PI(3),PY_PI(3),PZ_PI(3)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_ppp0(0,1),W_ppp0(1,1),W_ppp0(2,1),W_ppp0(3,1)
     %                     ,W_ppp0(4,1),W_ppp0(5,1),W_ppp0(6,1),W_ppp0(7,1))
     	   W2_ppp0(1) = W_ppp0(4,1)
          endif   
          endif
         endif           
        endif           		

c+++ Missing mass for k detected: different entries for differentk detected (max_k=2)
        if((nk).ge.1.and.(nk).le.max_k) then
         nek = nk
         j_l = 0
         do j=1,nk
          call miss_mass_2part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                        ,PMiss_ek(0,j),PMiss_ek(1,j),PMiss_ek(2,j),PMiss_ek(3,j)
     %                        ,PMiss_ek(4,j),PMiss_ek(5,j),PMiss_ek(6,j),PMiss_ek(7,j))
          M2_ek(j) = PMiss_ek(4,j)
c++++ Missing and invariant mass for K  (GG) detected
          if (ng.eq.2) then
             nek0 = nek0 + 1
             j_l = j_l + 1
          call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_ek0(0,j_l),PMiss_ek0(1,j_l),PMiss_ek0(2,j_l),PMiss_ek0(3,j_l)
     %                        ,PMiss_ek0(4,j_l),PMiss_ek0(5,j_l),PMiss_ek0(6,j_l),PMiss_ek0(7,j_l))

     	  M2_ek0(j_l) = PMiss_ek0(4,j_l) 
          call invariant_mass_2part
     %                        (E_K(j),PX_K(j),PY_K(j),PZ_K(j)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,W_k0(0,j_l),W_k0(1,j_l),W_k0(2,j_l),W_k0(3,j_l)
     %                        ,W_k0(4,j_l),W_k0(5,j_l),W_k0(6,j_l),W_k0(7,j_l))
     	  W2_k0(j_l) = W_k0(4,j_l)
          endif   
         enddo

c++++ Missing mass for k k detected: only 1 entry when max_k=2
         if(nk.ge.2.and.nk.le.max_k) then
          nekk = 1
c++++ Missing Mass
          call miss_mass_3part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                        ,PMiss_ekk(0,1),PMiss_ekk(1,1),PMiss_ekk(2,1),PMiss_ekk(3,1)
     %                        ,PMiss_ekk(4,1),PMiss_ekk(5,1),PMiss_ekk(6,1),PMiss_ekk(7,1))
          M2_ekk(1) = PMiss_ekk(4,1)

c++++++ Invariant mass
           call invariant_mass_2part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                     ,W_kk(0,1),W_kk(1,1),W_kk(2,1),W_kk(3,1)
     %                     ,W_kk(4,1),W_kk(5,1),W_kk(6,1),W_kk(7,1))
           W2_kk(1) = W_kk(4,1)

c++++ Missing and invariant mass for k k (GG) detected
          if (ng.eq.2) then
          nekk0 = 1
          call miss_mass_4part(m_nucleon(q_NU(1)),E0
     %                        ,E_EL,PX_EL,PY_EL,PZ_EL
     %                        ,E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                        ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                        ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                        ,PMiss_ekk0(0,1),PMiss_ekk0(1,1),PMiss_ekk0(2,1),PMiss_ekk0(3,1)
     %                        ,PMiss_ekk0(4,1),PMiss_ekk0(5,1),PMiss_ekk0(6,1),PMiss_ekk0(7,1))
     	   M2_ekk0(1) = PMiss_ekk0(4,1) 
           call invariant_mass_3part
     %                     (E_K(1),PX_K(1),PY_K(1),PZ_K(1)
     %                     ,E_K(2),PX_K(2),PY_K(2),PZ_K(2)
     %                     ,W_GG(0,1),W_GG(1,1),W_GG(2,1),W_GG(3,1)
     %                     ,W_kk0(0,1),W_kk0(1,1),W_kk0(2,1),W_kk0(3,1)
     %                     ,W_kk0(4,1),W_kk0(5,1),W_kk0(6,1),W_kk0(7,1))
     	   W2_kk0(1) = W_kk0(4,1)
          endif   
         endif
        endif



c---------------------------------------------------------------------
c     Quality check mode
c--------------------------------------------------------------------- 
c      IF (qc_flag) THEN

c ----- Select the events between the first and the last "blocks" -----
c ----- and Accumulate gated CHARGE                               -----
        IF(beam_type) THEN
          fcg_hevt = fcg_hevt
        ELSE
c1          dgamma = dgamma + NgammaTtotal*
c1     &             1./float(NEventTotalRun)/1.E+9
c2          dgamma = dgamma + NgammaTthisfile/float(NEventThisFile)/1.E+9
          dgamma = dgamma + 1.
          fcg_hevt = dgamma
          acc_el = 1.
        ENDIF !IF(beam_type) 
        key_endfile  = .false.
        key_endblock = .false. 
        key_count    = .false.
        Ifile_prev = Ifile_curr
        Ifile_curr = I_file         ! the current file number
        Iblock_curr = Iblock_curr   ! the current block number
	Iblock_all = Iblock_all
        IF(Ifile_prev.ne.Ifile_curr) THEN
          fcg_prev = 0.
          fcg_curr = fcg_hevt
          dgamma = 0.
          Qgated = 0.
          Iblock_curr = 1
	  Iblock_all = Iblock_all+1
          key_endfile  = .true.
          key_endblock = .false.
          key_count    = .false.
        ELSE
          fcg_prev = fcg_curr
          fcg_curr = fcg_hevt
          IF( fcg_prev .ne. fcg_curr ) THEN
            Qdiff = (fcg_curr - fcg_prev)
            Iblock_curr = Iblock_curr + 1  ! the current block number
	    Iblock_all = Iblock_all + 1
            IF(Iblock_curr.ge.3) key_endblock = .true.
            IF(Iblock_curr.ge.3) Qgated     = Qgated     + Qdiff
          ENDIF
          IF(Iblock_curr.ge.2) key_count    = .true.
        ENDIF
      ! The selected events are those that are not on the
      ! first or the last block.  In this case key_count is 
      ! set to .true. . Qgated accumulates Q for each file
      ! excluding first and last blocks.
      ! At the end of the block key_reset is set to .true.
      ! One has to:
      !   when key_endblock .eq. .true. then
      !     add to the (???) the portions of events which
      !     was taken for the passed block.
      !     And if you use some counters for a given data portions
      !     then reset them
      !   when key_endfile .eq. .true. then
      !     It means the end of file is reached
      !     you can reset some histograms/counters if you
      !     going to analyse the files separatly.
      !     If you a going to analyse all files as a hole, then do 
      !     not do anything
      !   when key_count .eq. .true. then
      !     Process events accumulating statistics somewhere
      !     for a given block. Only at the end of the block
      !     you have to add the portion of statistics to the
      !     to your hist.


        IF(beam_type) THEN
          Qgated = Qgated
        ELSE
          Qgated = NgammaTthisfile !/1.E+3
          IF( (.not.EgammaIsGood(Egamma)) .OR. (event_type.ne.1) ) THEN
            key_count = .false.
          ENDIF
        ENDIF


c--- The number of detected particles/Qgated ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(401,I_file, ndt_ebeam/Qgated,sqrt(ndt_ebeam*1.)/Qgated)
c           CALL SET_HCHANNELERR(402,I_file, ndt_enotb/Qgated,sqrt(ndt_enotb*1.)/Qgated)
c           CALL SET_HCHANNELERR(403,I_file, ndt_p/    Qgated,sqrt(ndt_p*1.)/    Qgated)
c           CALL SET_HCHANNELERR(405,I_file, ndt_ep/   Qgated,sqrt(ndt_ep*1.)/   Qgated)
c           CALL SET_HCHANNELERR(404,I_file, ndt_p/    Qgated,sqrt(ndt_p*1.)/    Qgated)
c           CALL SET_HCHANNELERR(405,I_file, ndt_n/    Qgated,sqrt(ndt_n*1.)/    Qgated)
c           CALL SET_HCHANNELERR(406,I_file, ndt_pip/  Qgated,sqrt(ndt_pip*1.)/  Qgated)
c           CALL SET_HCHANNELERR(407,I_file, ndt_pim/  Qgated,sqrt(ndt_pim*1.)/  Qgated)
c           CALL SET_HCHANNELERR(408,I_file, ndt_kp/   Qgated,sqrt(ndt_kp*1.)/   Qgated)
c           CALL SET_HCHANNELERR(409,I_file, ndt_km/   Qgated,sqrt(ndt_km*1.)/   Qgated)
c           CALL SET_HCHANNELERR(410,I_file, ndt_d/    Qgated,sqrt(ndt_d*1.)/    Qgated)
c           CALL SET_HCHANNELERR(411,I_file, ndt_g/    Qgated,sqrt(ndt_g*1.)/    Qgated)
             IF(.not.beam_type) THEN
             DO i=1,min(Erebin_Nbins,9)
             IF(Egamma.gt.Erebin_emin(i).and.
     &          Egamma.lt.Erebin_emax(i).and.
     &          Erebin_Ngamma(i).gt.0.0               )then
c           CALL SET_HCHANNELERR(800+i,I_file, ndt_ebeam/Erebin_Ngamma(i),sqrt(ndt_ebeam*1.)/Erebin_Ngamma(i))
c           CALL SET_HCHANNELERR(810+i,I_file, ndt_p/Erebin_Ngamma(i),    sqrt(ndt_p*1.)/Erebin_Ngamma(i))
c           CALL SET_HCHANNELERR(820+i,I_file, ndt_pip/Erebin_Ngamma(i),  sqrt(ndt_pip*1.)/Erebin_Ngamma(i))
c           CALL SET_HCHANNELERR(830+i,I_file, ndt_pim/Erebin_Ngamma(i),  sqrt(ndt_pim*1.)/Erebin_Ngamma(i))
             ENDIF
             ENDDO
             ENDIF
        ENDIF
        IF( key_endfile ) THEN 
          ndt_ebeam = 0
          ndt_enotb = 0
          ndt_ep    = 0
          ndt_p     = 0
          ndt_n     = 0
          ndt_pip   = 0
          ndt_pim   = 0
          ndt_kp    = 0
          ndt_km    = 0
          ndt_d     = 0
          ndt_g     = 0
        ENDIF
        IF( key_count ) THEN
          IF(nevnt.gt.0)        ndt_ebeam = ndt_ebeam+1
        DO j = 2,nevnt
          IF(ID_OUR(j).eq.  11) ndt_enotb = ndt_enotb+1
          IF(ID_OUR(j).eq. -11) ndt_ep    = ndt_ep   +1
          IF(ID_OUR(j).eq.2212) ndt_p     = ndt_p    +1
          IF(ID_OUR(j).eq.2112) ndt_n     = ndt_n    +1
          IF(ID_OUR(j).eq. 211) ndt_pip   = ndt_pip  +1
          IF(ID_OUR(j).eq.-211) ndt_pim   = ndt_pim  +1
          IF(ID_OUR(j).eq. 321) ndt_kp    = ndt_kp   +1
          IF(ID_OUR(j).eq.-321) ndt_km    = ndt_km   +1
          IF(ID_OUR(j).eq.  45) ndt_d     = ndt_d    +1
          IF(ID_OUR(j).eq.  22) ndt_g     = ndt_g    +1
        ENDDO
        ENDIF

c--- Zvertex ---
        IF( key_count) THEN
        if(beam_type .and. acc_el .ne. 0.0 ) then
          sect  = sector(ph_EL)
c          if(sect.eq.1) CALL HFILL( 161 , Z_EL,0., 1.0)
c          if(sect.eq.2) CALL HFILL( 162 , Z_EL,0., 1.0)
c          if(sect.eq.3) CALL HFILL( 163 , Z_EL,0., 1.0)
c          if(sect.eq.4) CALL HFILL( 164 , Z_EL,0., 1.0)
c          if(sect.eq.5) CALL HFILL( 165 , Z_EL,0., 1.0)
c          if(sect.eq.6) CALL HFILL( 166 , Z_EL,0., 1.0)
        elseif((.not.beam_type).and.
     &         nNU.eq.1.and.acc_nu(1).ne.0.0.and.q_nu(1).gt.0. ) then
          sect = sector(ph_NU(1))
c          if(sect.eq.1) CALL HFILL( 161 , Z_NU(1),0., 1.0)
c          if(sect.eq.2) CALL HFILL( 162 , Z_NU(1),0., 1.0)
c          if(sect.eq.3) CALL HFILL( 163 , Z_NU(1),0., 1.0)
c          if(sect.eq.4) CALL HFILL( 164 , Z_NU(1),0., 1.0)
c          if(sect.eq.5) CALL HFILL( 165 , Z_NU(1),0., 1.0)
c          if(sect.eq.6) CALL HFILL( 166 , Z_NU(1),0., 1.0)
        else
          continue
        endif
        ENDIF   

c--- Inclusive e/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(101, Iblock_all, 
c     &     n_e_incl/Qgated, 0.)
        ENDIF
        IF( key_endfile ) THEN 
          n_e_incl = 0
        ENDIF
        IF( key_count ) THEN
        IF (acc_el .gt. 0.0 ) THEN
              if ((elast.lt.46600).and.(elast.gt.42400).and.
     &  (tg_hevt.lt.0.94).and.(tg_hevt.gt.0.9)) then
          n_e_incl = n_e_incl + 1
        endif
        ENDIF
        ENDIF

c--- Inclusive ep/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(102, Iblock_all, 
c     &     n_ep_incl/Qgated, 0.)
        ENDIF
        IF(key_endfile) THEN
          n_ep_incl = 0
        ENDIF
        IF( key_count) THEN
        IF( acc_el .ne. 0.0 .and. 
     &      nNU.eq.1 .and. q_NU(1).gt.0.0 .and. acc_nu(1).ne.0.0
     &  )THEN
          n_ep_incl = n_ep_incl + 1
        ENDIF
        ENDIF   
        
c--- Elastic e(W=Mproton cut)/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(103, Iblock_all, 
c     &     n_e_elast/Qgated, 0.)
          elast = n_e_elast/Qgated
        ENDIF
        IF(key_endfile) THEN
          n_e_elast = 0
        ENDIF
        IF( key_count) THEN
        IF( acc_el .ne. 0.0 .and. w.ge.0.80 .and. w.le.1.00) THEN
          n_e_elast = n_e_elast + 1
          sect = sector(ph_EL)
c          if(sect.eq.1) CALL HFILL( 131 , W,0.,    1.0)
c          if(sect.eq.2) CALL HFILL( 132 , W,0.,    1.0)
c          if(sect.eq.3) CALL HFILL( 133 , W,0.,    1.0)
c          if(sect.eq.4) CALL HFILL( 134 , W,0.,    1.0)
c          if(sect.eq.5) CALL HFILL( 135 , W,0.,    1.0)
c          if(sect.eq.6) CALL HFILL( 136 , W,0.,    1.0)
        ENDIF
        ENDIF   

c--- Elastic ep(W=Mproton cut)/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(104, Iblock_all, 
c     &     n_ep_elast/Qgated, 0.)
c           CALL SET_HCHANNELERR(107, Iblock_all, 
c     &     tg_hevt, 0.)
            ENDIF
        IF(key_endfile) THEN
          n_ep_elast = 0
        ENDIF
        IF( key_count) THEN
        IF( acc_el .ne. 0.0 .and. 
     &      nNU.eq.1 .and. q_NU(1).gt.0.0 .and. acc_nu(1).ne.0 .and. 
     &      w.ge.0.85 .and. w.le.0.95 .and. 
     &      abs(th_NU(1)-
     &          180./3.1415*atan(1./(1.+E0/0.938)
     &          /tan(th_EL/2./180.*3.1415))).lt.10. .and.
     &      abs(abs((ph_NU(1)-ph_EL))-180.) .lt.5.
     &  )THEN
          n_ep_elast = n_ep_elast + 1
        ENDIF
        ENDIF   

c--- Delta ep(W=Delta cut)/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(105, Iblock_all, 
c     &     n_ep_delta/Qgated, 0.)
        ENDIF
        IF(key_endfile) THEN
          n_ep_delta = 0
        ENDIF
        IF( key_count) THEN
        IF( acc_el .ne. 0.0 .and. 
     &      nNU.eq.1 .and. q_NU(1).gt.0.0 .and. acc_nu(1).ne.0 .and.
     &      w.ge.1.15 .and. w.le.1.27
     &  )THEN
          n_ep_delta = n_ep_delta + 1
        ENDIF
        ENDIF   

c--- Twopions [eNp+p- or eNp+p-(mis)]/Qtot ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(106, Iblock_all, 
c     &     n_twopions/Qgated, 0.)
             IF(.not.beam_type) THEN
             DO i=1,min(Erebin_Nbins,9)
             IF(Egamma.gt.Erebin_emin(i).and.
     &          Egamma.lt.Erebin_emax(i).and.
     &          Erebin_Ngamma(i).gt.0.0               )then
           CALL SET_HCHANNELERR(840+i,I_file, 
     &     n_twopions/Erebin_Ngamma(i),sqrt(n_twopions*1.)/Erebin_Ngamma(i))
             ENDIF
             ENDDO
             ENDIF
        ENDIF
        IF(key_endfile) THEN
          n_twopions = 0
        ENDIF
        IF( key_count) THEN
        IF(
     &    ( ntot.eq.2 .and. 
     &      nNU.eq.1 .and. q_NU(1).gt.0.0 .and. acc_NU(1).ne.0.0 .and.
     &      (npi_plus.eq.1.or.npi_plus.eq.1).and.acc_pi(1).gt.0.0 .and. 
     &      neNp.ge.1 .and. 
     &      sqrt(M2_eNp(1)).ge.0.09 .and. sqrt(M2_eNp(1)).le.0.25 ) 
     &  )THEN
c     &    ( ntot.eq.3 .and.
c     &      acc_el.ne.0.0 .and. 
c     &      nNU.eq.1 .and. q_NU(1).gt.0.0 .and. acc_NU(1).ne.0.0 .and.
c     &      npi_plus.eq.1  .and. acc_pi(2).gt.0.0 .and. 
c     &      npi_minus.eq.1 .and. acc_pi(1).gt.0.0 ) 
          n_twopions = n_twopions + 1
        ENDIF
        ENDIF   

c--- The average number of particles per event ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNELERR(108, I_file, 
c     &     n_partavrtot*1., 0.)
c           CALL SET_HCHANNELERR(109, I_file, 
c     &     n_partavrevnt*1., 0.)
        ENDIF
        IF( key_endfile ) THEN 
          n_partpassed1 = 0
          n_partavrtot  = 0.
          n_partavrevnt = 0.
        ENDIF
        IF( key_count ) THEN
        IF (acc_el .ne. 0.0 ) THEN
          n_partpassed1 = n_partpassed1 + 1
          if(n_partpassed1.eq.1) then
          n_partavrtot  = ntot
          n_partavrevnt = nevnt
          else
          n_partavrtot  = ( n_partavrtot/
     &                      (n_partpassed1*1./(n_partpassed1-1.)) )
     &                    +(ntot *1./n_partpassed1*1.)
          n_partavrevnt = ( n_partavrevnt/
     &                      (n_partpassed1*1./(n_partpassed1-1.)) )
     &                    +(nevnt*1./n_partpassed1*1.)
          endif
        ENDIF
        ENDIF

c--- Just the number of events and those when acc. el. and Qgated ---
        IF(key_endblock) THEN
c           CALL SET_HCHANNEL(190, I_file, Qgated)
c           CALL SET_HCHANNEL(191, I_file, n_evtot*1.)
c           CALL SET_HCHANNEL(192, I_file, n_evaccel*1.)
        ENDIF
        IF( key_endfile ) THEN 
          n_evtot  = 0
          n_evaccel= 0
        ENDIF
        IF( key_count ) THEN 
          n_evtot = n_e_incl + 1          
        IF (acc_el .ne. 0.0 ) THEN
          n_evaccel = n_evaccel + 1          
        ENDIF
        ENDIF

c--- Mass_our and Mass_evnt distributions           --- 
c--- Pmom versus Beta plots for different particles ---
c--- in all files together                          ---
        IF( key_count ) THEN
        DO j = 1,nevnt
        pp =Pmom_EVNT(j)
        bb =Beta_our(j)
        bbe=Beta_evnt(j)
        if(Mass_our(j) .gt.0.0) mm =sqrt(Mass_our(j))
        if(Mass_our(j) .le.0.0) mm =-1000.
        if(Mass_evnt(j).gt.0.0) mme=sqrt(Mass_evnt(j))
        if(Mass_evnt(j).le.0.0) mme=-1000.
        call pseudo_spa2(int(Charge_EVNT(j)),0, Pmom_EVNT(j),
     &                  th_EVNT(1),(mod(ph_EVNT(1)+30.,60.)-30.),
     &                  I_torus,1., sector(ph_EVNT(1)), aa)
        IF( aa.ne.0.0) THEN
c          IF(j.eq.1                      ) CALL HFILL(203,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  11) CALL HFILL(204,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. -11) CALL HFILL(205,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. 211) CALL HFILL(206,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.-211) CALL HFILL(207,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. 321) CALL HFILL(208,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.-321) CALL HFILL(209,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.2212) CALL HFILL(210,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.2112) CALL HFILL(211,pp,bb,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  45) CALL HFILL(212,pp,bb,1.)
c          IF(j.ne.1.or.j.eq.1            ) CALL HFILL(213,pp,bb,1.)
c          IF(j.eq.1                       ) CALL HFILL(263,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  11) CALL HFILL(264,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. -11) CALL HFILL(265,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 211) CALL HFILL(266,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-211) CALL HFILL(267,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 321) CALL HFILL(268,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-321) CALL HFILL(269,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.2212) CALL HFILL(270,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.2112) CALL HFILL(271,pp,bbe,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  45) CALL HFILL(272,pp,bbe,1.)
c          IF(j.ne.1.or.j.eq.1             ) CALL HFILL(273,pp,bbe,1.)

c          IF(j.eq.1                      ) CALL HFILL(223,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  11) CALL HFILL(224,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. -11) CALL HFILL(225,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. 211) CALL HFILL(226,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.-211) CALL HFILL(227,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. 321) CALL HFILL(228,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.-321) CALL HFILL(229,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.2212) CALL HFILL(230,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.2112) CALL HFILL(231,mm,0.,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  45) CALL HFILL(232,mm,0.,1.)
c          IF(j.ne.1.or.j.eq.1            ) CALL HFILL(233,mm,0.,1.)
c          IF(j.eq.1                       ) CALL HFILL(243,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  11) CALL HFILL(244,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. -11) CALL HFILL(245,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 211) CALL HFILL(246,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-211) CALL HFILL(247,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 321) CALL HFILL(248,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-321) CALL HFILL(249,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.2212) CALL HFILL(250,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.2112) CALL HFILL(251,mme,0.,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  45) CALL HFILL(252,mme,0.,1.)
c          IF(j.ne.1.or.j.eq.1             ) CALL HFILL(253,mme,0.,1.)

        ENDIF
        ENDDO
        ENDIF

c--- Etot versus Pmom      --- 
c--- Ein  versus Eout      ---
c--- in all files together ---

        IF( key_count ) THEN
        DO j = 1,nevnt
        pp    = Pmom_EVNT(j)
        eetot = Etot_ECPB(j)/Pmom_EVNT(j)
        eein  = Ein_ECPB(j) /Pmom_EVNT(j)
        eeout = Eout_ECPB(j)/Pmom_EVNT(j)
        call pseudo_spa2(int(Charge_EVNT(j)),0,
     &             Pmom_EVNT(j),th_EVNT(1),(mod(ph_EVNT(1)+30.,60.)-30.),
     &             I_torus,1., sector(ph_EVNT(1)), aa)
ccc        IF(j.ne.1) THEN
ccc        write (81,*) ' QCQC==ID ,p,etot',ID_EVNT(j),Pmom_EVNT(j),Etot_ECPB(j),Etot_ECPB(j)/Pmom_EVNT(j)
ccc        ENDIF
        IF(aa.ne.0) THEN
        IF( abs(eein+eeout-eetot).lt. 0.1 ) THEN
c          IF(j.eq.1                       )CALL HFILL(301,pp,eetot,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.   11)CALL HFILL(302,pp,eetot,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  -11)CALL HFILL(303,pp,eetot,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  211)CALL HFILL(304,pp,eetot,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. -211)CALL HFILL(305,pp,eetot,1.)
c          IF(j.eq.1                       )CALL HFILL(311,eein,eeout,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.   11)CALL HFILL(312,eein,eeout,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  -11)CALL HFILL(313,eein,eeout,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq.  211)CALL HFILL(314,eein,eeout,1.)
c          IF(j.ne.1.and.ID_OUR(j).eq. -211)CALL HFILL(315,eein,eeout,1.)
c          IF(j.eq.1                       )CALL HFILL(351,pp,eetot,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  11)CALL HFILL(352,pp,eetot,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. -11)CALL HFILL(353,pp,eetot,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 211)CALL HFILL(354,pp,eetot,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-211)CALL HFILL(355,pp,eetot,1.)
c          IF(j.eq.1                       )CALL HFILL(361,eein,eeout,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.  11)CALL HFILL(362,eein,eeout,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. -11)CALL HFILL(363,eein,eeout,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq. 211)CALL HFILL(364,eein,eeout,1.)
c          IF(j.ne.1.and.ID_EVNT(j).eq.-211)CALL HFILL(365,eein,eeout,1.)
        ENDIF
        ENDIF
        ENDDO
        ENDIF


c      RETURN
c      ENDIF ! ending Run Quality checks




c+ Reaction identification
c++++++ Two pions
c++++     bit = 10 e p pi+
c+              Invariant          PMiss_NP(1)=(p pi+)   PMiss_NP(2)=(p M)     PMiss_PP(3)=(pi+ M)
c++++     bit = 11 e p pi-
c+              Invariant          PMiss_NP(1)=(p pi-)   PMiss_NP(2)=(p M)     PMiss_PP(3)=(pi- M)
c++++     bit = 12 e pi+ pi-
c+              Invariant          PMiss_NP(1)=(M pi+)   PMiss_NP(2)=(M pi-)   PMiss_PP(3)=(pi+ pi-)
c++++     bit = 13 e p pi+ pi-
c+              Invariant          PMiss_NP(1)=(p pi+)   PMiss_NP(2)=(p pi-)   PMiss_PP(3)=(pi+ pi-)
c++++++ One pion
c+        bit = 100
c++++++ Inclusive pion
c+        bit = 300
c++ TWO pions reactions RESERVED bits 0 - 99
c+++ Defining bit (exact final state)
        if (ana_chan.eq.0.or.ana_chan.eq.3) then
c+++ p pi+ pi-
c++++ the same as below but with a fixed number of particles (two or 3)
        if(ntot.eq.2.AND.nNU.eq.1.AND.q_NU(1).eq.1.AND.npi.eq.1.AND.npi_plus.eq.1) bit = 10     
        if(ntot.eq.2.AND.nNU.eq.1.AND.q_NU(1).eq.1.and.npi.eq.1.AND.npi_minus.eq.1) bit = 11   
        if(ntot.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1) bit = 12                             
        if(ntot.eq.3.AND.nNU.eq.1.and.q_NU(1).eq.1.AND.npi_plus.eq.1.AND.npi_minus.eq.1) bit = 13 
c--- Ending defining bit
c++++ Defining specific 2pion bit
c++++  the same as below but with a fixed number of particles (two or 3)
        if(nNU.eq.1.AND.q_NU(1).eq.1.AND.npi.eq.1.AND.npi_plus.eq.1) bit2pi = 20          
        if(nNU.eq.1.AND.q_NU(1).eq.1.and.npi.eq.1.AND.npi_minus.eq.1) bit2pi = 21      
        if(npi_plus.eq.1.AND.npi_minus.eq.1) bit2pi = 22                                 
        if(nNU.eq.1.and.q_NU(1).eq.1.AND.npi_plus.eq.1.AND.npi_minus.eq.1) bit2pi = 23   
c++++ Checking choosed bit
        if (bit.ge.10.and.bit.le.13)   bit2pi = bit ! put the more complete bit information
      if(bit2pi.gt.0.and.bit2pi.lt.100)  then
        call analysis_2pi(MCbit) ! Performing 2pi analysis 
        ChanDoneFlag(3,MCbit)=1
      endif

c-- Ending 2pi section
c++ ONE pion reactions RESERVED bits 100 - 299
c+++ Defining bit
      else if (ana_chan.eq.0.or.ana_chan.eq.1) then
c+++ p GG 111 -> 114
        if(ntot.eq.1.AND.nNU.eq.1.and.q_NU(1).eq.1) bit = 111
        if(ntot.eq.2.and.ngg.eq.1) bit = 112 ! to be extended
        if(ntot.eq.3.AND.nNU.eq.1.and.q_NU(1).eq.1.AND.ngg.eq.1) bit = 113
        if(ntot.eq.2.AND.nNU.eq.1.and.q_NU(1).eq.1.AND.ng.eq.1) bit = 114 ! p+1g detected
c++++ Defining specific 1pion bit 121 -> 124
        if(nNU.eq.1.and.q_NU(1).eq.1) bit1pi = 121
        if(ngg.eq.1) bit1pi = 122 ! to be extended
        if(nNU.eq.1.and.q_NU(1).eq.1.AND.ngg.eq.1) bit1pi =123 
        if(nNU.eq.1.and.q_NU(1).eq.1.AND.ng.eq.1) bit1pi = 124 ! p+1g detected
c++++ Checking choosed bit
        if (bit.ge.100.and.bit.le.120)   bit1pi = bit ! put the more complete bit information
        if(bit1pi.gt.100.and.bit2pi.lt.200)  then
          call analysis_1pi_a(MCbit) ! Performing 1pi analysis 
          ChanDoneFlag(1,MCbit)=1
        endif
      else if (ana_chan.eq.0.or.ana_chan.eq.2) then
c+++ n pi+ 211 -> 213
        if(ntot.eq.1.AND.nNU.eq.1.and.q_NU(1).eq.0) bit = 211
        if(ntot.eq.1.AND.npi.eq.1.AND.npi_plus.eq.1) bit = 212
        if(ntot.eq.2.AND.nNU.eq.1.and.q_NU(1).eq.0.AND.npi.eq.1.AND.npi_plus.eq.1) bit = 213
c++++ Defining specific 1pion bit 221 -> 223
        if(nNU.eq.1.and.q_NU(1).eq.0) bit1pi = 221
        if(npi.eq.1.AND.npi_plus.eq.1) bit1pi = 222
        if(nNU.eq.1.and.q_NU(1).eq.0.AND.npi.eq.1.AND.npi_plus.eq.1) bit1pi = 223
c--- Ending defining bit
c++++ Checking choosed bit
        if (bit.ge.200.and.bit.le.220)   bit1pi = bit ! put the more complete bit information
        if(bit1pi.gt.200.and.bit2pi.lt.300) then
          call analysis_1pi_b(MCbit) ! Performing 1pi analysis 
          ChanDoneFlag(2,MCbit)=1
        endif
c-- Ending 1pi section
c++++ Inclusive pion reserved bits 300-399
c+++ Defining bit
      else if (ana_chan.eq.0.or.ana_chan.eq.4) then
c+++at least 1pi
      if(nPI.ge.1) bit = 300
c--- Ending defining bit
      call analysis_inc_pi(MCbit) !performing the inclusive pion analysis
      ChanDoneFlag(4,MCbit)=1
c++++ Omega: reserved bits 400-499
      else if (ana_chan.eq.0.or.ana_chan.eq.5) then
       if(nNU.eq.1.and.q_NU(1).eq.1)  then
         if(ntot.eq.1)                                 bit = 400 ! P    detected
         if(ntot.eq.2.AND.ng.eq.1)                     bit = 401 ! P+1g detected
         if(ntot.eq.3.AND.ngg.eq.1)                    bit = 402 ! P+2g detected
         if(ntot.eq.2.AND.npi.eq.1.AND.npi_plus.eq.1)              bit = 405 ! P+PIp   detected    
         if(ntot.eq.3.AND.npi.eq.1.AND.npi_plus.eq.1.AND.ng.eq.1)  bit = 406 ! P+PIp+G   detected    
         if(ntot.eq.4.AND.npi.eq.1.AND.npi_plus.eq.1.AND.ngg.eq.1) bit = 407 ! P+PIp+GG   detected    
         if(ntot.eq.2.AND.npi.eq.1.AND.npi_minus.eq.1)              bit = 410 ! P+PIm   detected    
         if(ntot.eq.3.AND.npi.eq.1.AND.npi_minus.eq.1.AND.ng.eq.1)  bit = 411 ! P+PIm+G   detected    
         if(ntot.eq.4.AND.npi.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1) bit = 412 ! P+PIm+GG   detected    
         if(ntot.eq.3.AND.npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1)              bit = 415 ! P+PIp+PIm    detected    
         if(ntot.eq.4.AND.npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ng.eq.1)  bit = 416 ! P+PIp+PIm+G  detected    
         if(ntot.eq.5.AND.npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1) bit = 417 ! P+PIp+PIm+GG detected
       endif
       if(ntot.eq.4.AND.npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1)    bit = 419 ! PIp+PIm+GG detected   
       if(nNU.eq.1.and.q_NU(1).eq.1)  then 
                                                    bitomega = 430 ! P    detected
        if( ng.eq.1)                                 bitomega = 431 ! P+1g detected
        if( ngg.eq.1)                                bitomega = 432 ! P+2g detected
        if( npi.eq.1.AND.npi_plus.eq.1)              bitomega = 435 ! P+PIp   detected    
        if( npi.eq.1.AND.npi_plus.eq.1.AND.ng.eq.1)  bitomega = 436 ! P+PIp+G   detected    
        if( npi.eq.1.AND.npi_plus.eq.1.AND.ngg.eq.1) bitomega = 437 ! P+PIp+GG   detected    
        if( npi.eq.1.AND.npi_minus.eq.1)              bitomega = 440 ! P+PIm   detected    
        if( npi.eq.1.AND.npi_minus.eq.1.AND.ng.eq.1)  bitomega = 441 ! P+PIm+G   detected    
        if( npi.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1) bitomega = 442 ! P+PIm+GG   detected    
        if( npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1)              bitomega = 445 ! P+PIp+PIm    detected    
        if( npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ng.eq.1)  bitomega = 446 ! P+PIp+PIm+G  detected    
        if( npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1) bitomega = 447 ! P+PIp+PIm+GG detected  
       endif
       if(npi.eq.2.AND.npi_plus.eq.1.AND.npi_minus.eq.1.AND.ngg.eq.1)  bitomega = 449 ! PIp+PIm+GG detected    
       if (bit.ge.400.and.bit.le.419)   bitomega = bit ! put the more complete bit information
       if (bitomega.gt.0.and.bitomega.lt.499) then
         ChanDoneFlag(5,MCbit)=1
         stop !call analysis_omega(MCbit) ! Performing some analysis
       endif
      endif
c- Ending reaction id
c++++
1233   return
      end 


c======================================================================
c InitEvent_goa:
c======================================================================
      SUBROUTINE initevent_goa      
      include "ntpl_goa.inc"

c general initialisations
      indtype = 0


      RETURN 
      END




c======================================================================
c Preprocess file:
c   1) Select Good events (skipping 1st and last FCG-blocks) and fill 
c      an array Event_OK(), save FCG_buf(),... arrays
c   2) Helicities preprocessing selecting good 
c   3) Calculate medium Z coord. for Electron vertex position and save 
c      it in the z_el_mid()
c======================================================================
      SUBROUTINE prepfile_goa(Nmax)      
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "select.inc"
      LOGICAL  endread
      INTEGER*4  NMAX,ievent,ierr,i,j,ii,nn,tmp,sect,sector
      INTEGER*4  nbuf_fcg
      INTEGER*4  numtot,ngood,nbad
      CHARACTER*80 ch80
      REAL*4  w2_local 
      REAL*4  par(10),step(10),pmin(10),pmax(10),perr(10),c2(10)
      REAL*4  hhh(1000)
      REAL*4  HMAX,elast
      REAL*4   fit_boxfun2
      EXTERNAL fit_boxfun2
      REAL*4   fcg_prev,fcg_curr
      DATA     fcg_prev,fcg_curr/-1.,-1./
      REAL*4   m_nucleon(0:1)
      INTEGER*4 Ntrigger_first,Ntrigger_last
      LOGICAL   firstphysics
      integer Iblock_all,n_e_elast,n_e_incl
      DATA m_nucleon/0.93955,0.93821/ !0=neutron, 1=proton
      COMMON /qqq/Iblock_all,n_e_elast,elast,n_e_incl

      IF(Nmax.gt.nevent_max) THEN
        print *,' MMAX.gt.nevent_max==',Nmax,nevent_max
        stop
      ENDIF
      firstphysics = .false.   
      Ntrigger_file=0
      Ncebevent_file=0
      Ntrigger_first=0
      Ntrigger_last =-1
      Qgated_file=0.0

c------------------------------------------------------
c No preprocessing if MC analysis or DATA&MC together
c------------------------------------------------------

      IF(bit_mc.ne.0.or.qc_flag) THEN
        do ievent=1,Nmax
        Event_OK(ievent) = .true.
        enddo
        GOTO 1002
      ENDIF


c--------------------------------------------------
c Read the whole DATA ntupl 
c--------------------------------------------------
      endread = .false.
      ievent=0
      IF(bit_mc.eq.0) THEN
      DO WHILE (.true.)

c------------------------------------------------------
c       read one event 
c------------------------------------------------------
        ievent=ievent+1
        CALL HGNT(input_ntdata_num,ievent,ierr)
        if(input_type.eq.2) call clas2ceb
        IF(ierr.ne.0 .or. ievent.gt.Nmax) THEN
          endread = .true.
          ievent=ievent-1
          nevent_read = ievent
        ENDIF
        IF(ievent.gt.nevent_max) THEN
          print *,' ievent .gt. nevent_max ',ievent,nevent_max
          stop
        ENDIF
        Event_OK(ievent) = .true.

        IF(event_type.eq.1) THEN
          Ntrigger_last  = event_num
          NEventFileLast = Ntrigger_last 
          IF((.not.firstphysics)) THEN
          Ntrigger_first = event_num
          NEventFileFirst= Ntrigger_first 
          firstphysics = .true.
          ENDIF
        ENDIF
        Ntrigger_file = Ntrigger_last-Ntrigger_first+1


c------------------------------------------------------
c        do with Qgated, select good events 
c------------------------------------------------------

c --- save fcg values by buffers ---
        IF( .not. endread) THEN
          fcg_prev = fcg_curr
          fcg_curr = fcg_hevt  !************************ FCG
          IF(ievent.eq.1) nbuf_fcg=0
          IF(ievent.eq.1) indBegin_fcg(1)=1
          IF(fcg_curr.ne.fcg_prev .or. ievent.eq.1) THEN
            nbuf_fcg = nbuf_fcg + 1
            IF(nbuf_fcg.gt.nbufmax1) THEN
              print *,' nbuf_fcg .gt. nbufmax1==',nbuf_fcg
              stop
            ENDIF
            indBegin_fcg(nbuf_fcg)=ievent
            IF(nbuf_fcg.ge.2) indEnd_fcg(nbuf_fcg-1)=ievent-1
            fcg_buf(nbuf_fcg)   = fcg_curr
          ENDIF 
c --- select good events ---
        ELSE !( .not. endread)
          indEnd_fcg(nbuf_fcg)=ievent
          nbuf_read_fcg=nbuf_fcg
          numtot=0    
          ngood =0
          nbad  =0
          DO i=1,nbuf_read_fcg  ! 1st & last blocks & block with fcg=0
            IF(i.ne.nbuf_read_fcg) fcg_bufplus(i)=fcg_buf(i+1)
            IF(i.eq.nbuf_read_fcg) fcg_bufplus(i)=fcg_buf(i  )
            IF(i.eq.1 .or. i.eq.nbuf_read_fcg .or. 
     &         fcg_buf(i).eq.0. .or. fcg_bufplus(i).eq.0. )THEN
              DO j=indBegin_fcg(i),indEnd_fcg(i)
              Event_OK(j) = .false.            
              ENDDO
              nbad = nbad +(indEnd_fcg(i)-indBegin_fcg(i)+1)
            ELSE
              ngood= ngood+(indEnd_fcg(i)-indBegin_fcg(i)+1)
              Ncebevent_file = ngood
	      
	      
	      
	      if ((elast.lt.46600).and.(elast.gt.42400).and.
     &  (tg_hevt.lt.0.94).and.(tg_hevt.gt.0.9)) then  ! Gleb 6/16/2003
              
c          if ((z_EL.lt.-2.358).or.(z_EL.gt.1.369)) then ! Gleb 7/16/2003
              
	      
	      Qgated_file = Qgated_file + fcg_bufplus(i)-fcg_buf(i)
	      
	      endif   ! Gleb 6/16/2003
            ENDIF
          ENDDO
        write(*,74) ngood*100./(ngood+nbad),nbad*100./(ngood+nbad)
   74   format('   selected by Qgated good,bad (%)==',2f8.3)
        Ntrigger_sum  = Ntrigger_sum  + Ntrigger_file
        Ncebevent_sum = Ncebevent_sum + Ncebevent_file 
        Qgated_sum    = Qgated_sum    + Qgated_file
        write(*,76) Ntrigger_file,Ncebevent_file,Qgated_file
        write(*,*) 'n_e_incl', n_e_incl
   76   format('   Ntrig,Ncebs,Qgated: InFile: ',2i9,'  ',f10.4)
        write(*,77) Ntrigger_sum, Ncebevent_sum, Qgated_sum
   77   format('   Ntrig,Ncebs,Qgated: Total : ',2i9,'  ',f10.4)
        ENDIF !( .not. endread)
        

c-------------------------------------------------------------
c        do with Zvertex histogram, filling them then fitting 
c-------------------------------------------------------------

c --- fill Z vertex histogram taking e and p ---

        IF( mode_zcor ) THEN
        IF( .not. endread) THEN          
	  omega = E0 - pmom_EVNT(1)
	  Q2    = 4.*E0*pmom_EVNT(1)*(sin(th_EVNT(1)*3.1416/2./180.))**2
	  W2_local = m_nucleon(1)**2+2.*m_nucleon(1)*omega-Q2
	  if(W2_local.ge.0.)  W = sqrt(W2_local)
          th_EL = th_EVNT(1)
          ph_EL = ph_EVNT(1)
          E_EL  = E0 - Omega
          P_EL  = E_EL
          Z_EL  = Z_EVNT(1)
          call pseudo_spa2(0,11,p_EL,th_EL,(mod(ph_el+30.,60.)-30.),
     &                     I_torus,1.,sector(ph_EL),acc_EL)
          IF( acc_el .ne. 0.0 .and. w.ge.0.01 .and. w.le.9.99) THEN
          sect = sector(ph_EL)
          if(sect.eq.1) CALL HFILL( 161 , Z_EL,0., 1.0)
          if(sect.eq.2) CALL HFILL( 162 , Z_EL,0., 1.0)
          if(sect.eq.3) CALL HFILL( 163 , Z_EL,0., 1.0)
          if(sect.eq.4) CALL HFILL( 164 , Z_EL,0., 1.0)
          if(sect.eq.5) CALL HFILL( 165 , Z_EL,0., 1.0)
          if(sect.eq.6) CALL HFILL( 166 , Z_EL,0., 1.0)
          ENDIF
          DO j=1,nevnt
          IF(ID_EVNT(j).eq.2212) THEN
          th_NU(1) = th_EVNT(j)
          ph_NU(1) = ph_EVNT(j)
          P_NU(1)  = Pmom_EVNT(j)
          Z_NU(1)  = Z_EVNT(j)
          call pseudo_spa2(+1,2212,  p_NU(1),th_NU(1),
     &                    (mod(ph_NU(1)+30.,60.)-30.),I_torus,
     &                    1.,sector(ph_NU(1)),acc_NU(1))
          IF( acc_el .ne. 0.0 .and. w.ge.0.01 .and. w.le.9.99) THEN
          sect = sector(ph_NU(1))
          if(sect.eq.1) CALL HFILL( 171 , Z_NU(1),0., 1.0)
          if(sect.eq.2) CALL HFILL( 172 , Z_NU(1),0., 1.0)
          if(sect.eq.3) CALL HFILL( 173 , Z_NU(1),0., 1.0)
          if(sect.eq.4) CALL HFILL( 174 , Z_NU(1),0., 1.0)
          if(sect.eq.5) CALL HFILL( 175 , Z_NU(1),0., 1.0)
          if(sect.eq.6) CALL HFILL( 176 , Z_NU(1),0., 1.0)
          ENDIF
          ENDIF
          ENDDO
c --- fit Z vertex histogram and save Z vertex mid-value ---
        ELSE  ! (.not. endread)
          IF(mode_zcortype.eq.1 .or. mode_zcortype.eq.0)THEN

          DO j=1,6
          i=160+j
          par(2)=0.0
          par(3)=4.0
          par(1)=HMAX(i)/2.
          CALL HFITH (i,fit_boxfun2,'QU',3,par,step,pmin,pmax,perr,c2)
          CALL HGIVE (i,ch80,nn,tmp,tmp,tmp,tmp,tmp,tmp,tmp)
          CALL HUNPAK(i,hhh,'HIST',0)
          if(nn.gt.1000) stop 11122277
          IF(mode_zcortype.eq.1) THEN
           do ii=1,nn
           if(hhh(ii).gt.0.91*par(1)) hhh(ii)=0.9*par(1)
           if(hhh(ii).lt.0.17*par(1)) hhh(ii)=0.0
           enddo
          ENDIF
          CALL HPAK  (i,hhh)
          CALL HFITH (i,fit_boxfun2,'QU',3,par,step,pmin,pmax,perr,c2)
          z_el_mid(j) = par(2)
          ENDDO
          DO j=1,6
          i=170+j
          par(2)=0.0
          par(3)=4.0
          par(1)=HMAX(i)/2.
          CALL HFITH (i,fit_boxfun2,'QU',3,par,step,pmin,pmax,perr,c2)
          CALL HGIVE (i,ch80,nn,tmp,tmp,tmp,tmp,tmp,tmp,tmp)
          CALL HUNPAK(i,hhh,'HIST',0)
          do ii=1,nn
          if(hhh(ii).gt.0.91*par(1)) hhh(ii)=0.9*par(1)
          if(hhh(ii).lt.0.17*par(1)) hhh(ii)=0.0
          enddo
          CALL HPAK  (i,hhh)
          CALL HFITH (i,fit_boxfun2,'QU',3,par,step,pmin,pmax,perr,c2)
          z_p_mid(j) = par(2)
          ENDDO

          ELSEIF(mode_zcortype.eq.2)THEN
          
          DO j=1,6
          i=160+j
          par(1)=HMAX(i)
          par(2)=2.2
          par(3)=1.0
          par(4)=HMAX(i)/3.
          par(5)=-2.2
          par(6)=1.0
          par(7)=HMAX(i)/30.
          CALL HFITHN(i,'G+G+P0','QU',7,par,step,pmin,pmax,perr,c2)
          z_el_mid(j) = par(2)-2.0
          ENDDO
          DO j=1,6
          i=170+j
          par(1)=HMAX(i)
          par(2)=2.2
          par(3)=1.0
          par(4)=HMAX(i)/3.
          par(5)=-2.2
          par(6)=1.0
          par(7)=HMAX(i)/12.
          CALL HFITHN(i,'G+G+P0','QU',7,par,step,pmin,pmax,perr,c2)
          z_p_mid(j) = par(2)-2.0
          ENDDO
          
          ENDIF
          write(*,75) (z_el_mid(j),j=1,6),(z_p_mid(j),j=1,6)
   75     format('   zv.mid. e,p ==',12f5.1) 
          DO j=1,6
            CALL HRESET(160+j,' ')
            CALL HRESET(170+j,' ')
          ENDDO

        ENDIF !( .not. endread)
        ENDIF !( mode_zcor )

      IF( endread ) GOTO 1002
      ENDDO !DO WHILE (.true)
      ENDIF !IF(NDATA.gt.1) THEN
 1002 CONTINUE

      END



c======================================================================
c Preprocessing for photons:
c======================================================================
      SUBROUTINE prepfile_goa_photon(file_in,Nmax)      
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "photon.inc"
      CHARACTER*200 file_in
      INTEGER Nmax,i1,i2

c --- recalibrate or not recalibrate
        flag_false_g = .false.
c        tvsst = 0.       ! IS OVERWRITTEN during calibration !
c        tag_tof_corr= 0. ! IS OVERWRITTEN during calibration !
        coin_gate_ST_TAG = 3.0
        IF(input_type.eq.1 .OR. .not.do_recalibrate) THEN !no recalib
          print *,'  NO recalibration tvsst,tag_tof_corr= ',tvsst,tag_tof_corr
          do i1=1,7    ! tvsst and tag_tof_corr have been read from input
          do i2=1,48   !
          sl_i(i1,i2) = 1.
          t0_i(i1,i2) = -(tvsst+tag_tof_corr) 
          enddo
          enddo
          RETURN
        ENDIF
        IF(do_recalibrate) THEN
          print *,'  YES recalibration coin_gate_ST_TAG= ',coin_gate_ST_TAG
          tvsst = 0.        ! will be overitten
          tag_tof_corr = 0. ! will be overitten
          do i1=1,7    !Initializing delay and slope used in tof re-calibration
          do i2=1,48   !used in tof re-calibration
          sl_i(i1,i2) = 1.
          t0_i(i1,i2) = -(tvsst+tag_tof_corr) 
          enddo
          enddo
          OPEN (UNIT=1,FILE='calib_resume.txt',STATUS='UNKNOWN')
        ENDIF

c --- TAGGER RECALIBRATION  ---
        print *,' '
        print *,'Recalibration TAGGER/ST coin_gate_ST_TAG=',coin_gate_ST_TAG
        CALL TAGGER_RECALIB(file_in,5000)
        print *,'Fit TAGGER/ST histo'
        CALL TAGGER_fit_histo()

c --- TOF RECALIBRATION  ---
        print *,' '
        print *,'Recalibration TOF'
        CALL TOF_RECALIB(file_in,50000)
        print *,'Fit TOF histo'
        CALL tof_fit_histo()

        CLOSE(1)

      RETURN
      END




c======================================================================
c PostProcess file
c======================================================================
      SUBROUTINE postfile_goa
      IMPLICIT none 
      include "ntpl_goa.inc"
      integer NP,i,noen
      real    par(10),step(10),pmin(10),pmax(10),perr(10),c2(10)
      real    A(10000)
      REAL     HMAX
      REAL     fit_boxfun2
      EXTERNAL fit_boxfun2

      IF(qc_flag) THEN
      DO i=1,6
        np=3
        par(2)=0.9
        par(3)=0.05
        par(1)=HMAX(130+i)/2.
        CALL HNOENT(130+i,noen)  
          if(noen.gt.100) then
        CALL HFITHN(130+i,'G','Q',np,par,step,pmin,pmax,perr,c2)
        par(2)=abs(par(2))
        par(3)=abs(par(3))
        par(1)=abs(par(1))
        CALL HFITHN(130+i,'G','Q',np,par,step,pmin,pmax,perr,c2)
          else
        par(2)=-1.
        par(3)=-1.
        print *,' !!! EMPTY HISTO 130+i !!!'
          endif
        CALL HRESET(130+i,' ')
        CALL SET_HCHANNELERR(110+i, I_file, par(2),perr(2)) 
        CALL SET_HCHANNELERR(120+i, I_file, par(3),perr(3))
        np=3
        par(2)=0.0
        par(3)=20.0  !4.
        par(1)=HMAX(160+i)/2.
        CALL HFITH (160+i,fit_boxfun2,'QU',np,par,step,pmin,pmax,perr,c2)
        par(3)=20.0  !4. 
        par(1)=abs(par(1))
        CALL HFITH (160+i,fit_boxfun2,'QU',np,par,step,pmin,pmax,perr,c2)
        CALL SET_HCHANNELERR(140+i, I_file, par(2),perr(2)) 
        CALL SET_HCHANNELERR(150+i, I_file, par(3),perr(2))
      ENDDO

      CALL HUNPAK(11009,a,' ',0)
      do i=1,767
        if(NGammaE(i).gt.0.) a(i) = sqrt(abs(a(i))/NGammaE(i))
      enddo
      CALL HPAKE(11009,a)
      CALL HCOPY(11009, 12000+i_file, 'Nev/Ng for -- '//file_input)
      CALL HRESET(11009,' ')

      ENDIF
 
      RETURN
      END



c======================================================================
c Short Preprocess file:
c Golovach. Do not delete it, even though it's not used
c======================================================================
      SUBROUTINE shortprepfile_goa(Nmax)
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "select.inc"
      LOGICAL  endread, firstphysics
      integer Iblock_all,n_e_elast,n_e_incl
      INTEGER*4  Ntrigger_first,Ncebevent_first
      INTEGER*4  NMAX,ievent,ierr,Ntrigger,Ncebevent
      REAL*4     QQgg, Qgated_first,elast
      COMMON /qqq/Iblock_all,n_e_elast,elast,n_e_incl
      IF(Nmax.gt.nevent_max) THEN
        print *,' MMAX.gt.nevent_max==',Nmax,nevent_max
        stop
      ENDIF

      Ncebevent = 0
      firstphysics = .false.
      endread = .false.
      ievent=0

c ----- Just get Ncebevent_sum 
      call Hnoent(input_ntdata_num,Ncebevent)
      Ncebevent_sum = Ncebevent_sum + Ncebevent
      print *,' Ncebevent_sum =',Ncebevent_sum
      RETURN      


      IF(bit_mc.eq.0) THEN
      DO WHILE (.true.)

c------------------------------------------------------
c       read one event 
c------------------------------------------------------
        ievent=ievent+1
        CALL HGNT(input_ntdata_num,ievent,ierr)
        if(input_type.eq.2) call clas2ceb
        IF(ierr.ne.0 .or. ievent.gt.Nmax) THEN
          endread = .true.
          ievent=ievent-1
          nevent_read = ievent
        ENDIF
        IF(ievent.gt.nevent_max) THEN
          print *,' ievent .gt. nevent_max ',ievent,nevent_max
          stop
        ENDIF
        IF(event_type.eq.1) THEN
          Ntrigger  = event_num
          Ncebevent = Ncebevent + 1
          if(beam_type) then
          QQgg = fcg_hevt
          endif
          firstphysics = .true.
        ENDIF

c --- just Sum up Nev Qgat/Ngam ---
      IF( Ncebevent.eq.1 ) THEN
        Ntrigger_first  = Ntrigger
        Ncebevent_first = Ncebevent
        if(beam_type) then
        Qgated_first    = QQgg
        else
        Qgated_first    = 0.
        endif
      ENDIF
      IF( endread ) THEN 
        Ntrigger_sum  = Ntrigger_sum  + (Ntrigger - Ntrigger_first +1)
        Ncebevent_sum = Ncebevent_sum + (Ncebevent- Ncebevent_first+1)
	
       
       
       
*       if ((elast.lt.49566).and.(elast.gt.39240).and.
*     &  (tg_hevt.lt.0.94).and.(tg_hevt.gt.0.9)) then  ! Gleb 6/16/2003 
         
	 if ((z_EL.lt.-2.358).or.(z_EL.gt.1.369)) then ! Gleb 7/16/2003    
     
     
     
     
        if(beam_type) then
        Qgated_sum    = Qgated_sum    + (QQgg - Qgated_first)
        else
        Qgated_sum    = Qgated_sum    + (NGammaTthisfile)
        endif
	
	
	
	endif  ! Gleb 6/16/2003 

      ENDIF

      IF( endread ) GOTO 1002
      ENDDO !DO WHILE (.true)
      ENDIF !IF(NDATA.gt.1) THEN
 1002 CONTINUE

      print *,' '
      print *,' Ntrigger =',Ntrigger_sum
      print *,' Ncebevent=',Ncebevent_sum
      print *,' Qgated   =',Qgated_sum


      RETURN
      END




c======================================================================
c======================================================================
      SUBROUTINE SET_HCHANNEL(ID,ICH,VAL)
      INTEGER ID,ICH
      REAL    VAL, A(10000)
      CALL HUNPAK(id,a,' ',0)
      a(ICH) = val
      CALL HPAK(id,a)
      RETURN
      END  

      SUBROUTINE SET_HCHANNELERR(ID,ICH,VAL,ERR)
      INTEGER ID,ICH
      REAL    VAL,ERR, A(10000)
      CALL HUNPAK(id,a,' ',0)
      a(ICH) = val
      CALL HPAK(id,a)
      CALL HUNPKE(id,A,' ',0)
      a(ICH) = err
      CALL HPAKE(id,a)
      RETURN
      END

      SUBROUTINE SET_HCHANNELERRSQRT(ID,ICH,VAL)
      INTEGER ID,ICH
      REAL    VAL, A(10000)

        print *,' HCHANNELERRSQRT: ID,ICH,VAL==',ID,ICH,VAL

      CALL HUNPAK(id,a,' ',0)
      a(ICH) = val
      CALL HPAK(id,a)
      CALL HUNPKE(id,A,' ',0)
      a(ICH) = sqrt(VAL)
      CALL HPAKE(id,a)
      RETURN
      END




c======================================================================
c------------------------ FIT_BOXFUN ----------------------------------
c box like function "__|---|___"
c======================================================================
      FUNCTION fit_boxfun(x)
      implicit none
      REAL   fit_boxfun,x
      DOUBLE PRECISION dpar(24),fitfun, arg,center,w2,height
      COMMON   /hcfitd/dpar,fitfun
c - 3 parameters are:
c - dpar(1) - box height
c - dpar(2) - center
c - dpar(3) - width
      
      arg = x
      center = dpar(2)
      w2     = dpar(3)/2.d+0
      height = dpar(1)
      IF( arg .GE. center-w2 .AND. arg .LE. center+w2 ) THEN
        fitfun     = height
        fit_boxfun = height
      ELSE
        fitfun     = 0.0
        fit_boxfun = 0.0
      ENDIF

      RETURN
      END

c======================================================================
c------------------------ FIT_BOXFUN2 ---------------------------------
c fit by function like "__/---\___"
c======================================================================
      FUNCTION fit_boxfun2(x)
      implicit none
      REAL   fit_boxfun2,x
      DOUBLE PRECISION dpar(24),fitfun, arg,center,w2,height,delt
      COMMON   /hcfitd/dpar,fitfun
c - 3 parameters are:
c - dpar(1) - box height
c - dpar(2) - center
c - dpar(3) - width
      
      arg = x
      center = dpar(2)
      w2     = dpar(3)/2.d+0
      height = dpar(1)
      delt   = w2/6.
      IF    (arg .GE. center-w2+delt .AND. arg .LE. center+w2-delt)THEN
        fitfun      = height
        fit_boxfun2 = fitfun
      ELSEIF(arg .GE. center-w2-delt .AND. arg .LE. center-w2+delt)THEN
        fitfun      = (+height/2./delt) * (x - (center-w2-delt) )
        fit_boxfun2 = fitfun
      ELSEIF(arg .GE. center+w2-delt .AND. arg .LE. center+w2+delt)THEN
        fitfun      = (-height/2./delt) * (x - (center+w2+delt) )
        fit_boxfun2 = fitfun
      ELSE
        fitfun      = 0.0
        fit_boxfun2 = fitfun 
      ENDIF

      RETURN
      END
      
