
c----------------------------------------------------------------------
c                  
c----------------------------------------------------------------------

      SUBROUTINE CLAS2CEB
      IMPLICIT none
      include "ntpl_goa.inc"
      include "ntpl_clas.inc"
      include "photon.inc"
      INTEGER nevent_counter,i,jj,j_good,iii
      INTEGER NWarn1,NWarn2,NWarn3
      REAL    pigr
      DATA    nevent_counter/0/     
      DATA    NWarn1/0/
      DATA    NWarn2/0/
      DATA    NWarn3/0/

      REAL    STT_HEVT
      INTEGER sec(6),ST_pair(3),ST_all
      INTEGER hit_sc2SCPB(6,50)

      pigr=acos(-1.0)
      nevent_counter = nevent_counter + 1


c----------------------------------------------------------------------
c     HEAD   <-> HEVT
c     SEB_RE
c     SEB_Q
c----------------------------------------------------------------------
      event_num  = clas_evntid
      event_type = clas_evtype
      event_clas = clas_evntclas
      nevent     = nevent_counter
      fcg_hevt   = clas_q_l
      tg_hevt    = clas_t_l
      joppa = clas_tr_time
c----------------------------------------------------------------------
c     POL   <-> not done
c     HLS
c     FBPM
c     TGBI
c----------------------------------------------------------------------
      helicity = 0
      nHLS = 0
      do i=1,nHLS
       c_BPM(i)  =0
       c_FC(i)   =0
      enddo
      rx_FBPM=0
      ry_FBPM=0
      nTGBI=0
      do i=1,nTGBI
       helTGBI(i)=0
      enddo


c----------------------------------------------------------------------
c     for photon analysis:
c     PH_IN <-> TGPB, TAGR, HEVT
c----------------------------------------------------------------------
        IF(.not.beam_type) THEN

c --- get TGPB bank and j_good
      j_good = 0
      ph_hit_in_time(1) = 0
      nTGPB =clasg6_taggoodhit
      DO i=1,nTGPB
      tdiff_test = 100.
      flag_good_ph = .false.
      pointer_TGPB(i)= clasg6_tag_ptr(i)
      Time_TGPB(i)   = clasg6_vertex_time(i)
      Energy_TGPB(i) = clasg6_tag_energy(i)
      dt_TGPB(i)     = clasg6_dt_st_tag(i)

      STT_HEVT     = (dt_TGPB(i)+Time_TGPB(i)) !Over 
      ERG_TAGR(i)  = Energy_TGPB(i)            !Over
      TPHO_TAGR(i) = Time_TGPB(i)              !Over 

      if (i.eq.1) STT_HEVT = STT_HEVT - tvsst
      tdiff =abs(STT_HEVT-TPHO_TAGR(i)) !see ntpl_com.h to see tvsst
      if(tdiff.lt.coin_gate_ST_TAG)ph_hit_in_time(1)=ph_hit_in_time(1)+1
c     Chechking hits in ST/TAG window (see ntpl_com.h foe win value
      if (tdiff.lt.tdiff_test) then                                    
        tdiff_test = tdiff
        j_good = i     ! Recording position in TGPB  bank of good hit
        flag_good_ph = .true.
      endif
      ENDDO

      
c --- save TAGR bank only j_good
      n_ph_in=1 !clasg6_taghit ! loop over tagger hits
      if(j_good.gt.0) then
      do i=1,1
      ph_ERG_TAGR(i)  = ERG_TAGR(j_good)  ! was redifined
      ph_TPHO_TAGR(i) = TPHO_TAGR(j_good) ! was redifined
      ph_STT_HEVT(i)  = STT_HEVT          ! was redifined
      ph_RF1_HEVT(i)  = clas_rf_time1 ! from hevt bank
      ph_RF2_HEVT(i)  = clas_rf_time2 ! from hevt bank
      ph_T_id_TAGR(i) = clasg6_Tid(clasg6_tag_ptr(j_good)) ! really from TAGR bank no overwritten
      ph_E_id_TAGR(i) = clasg6_Eid(clasg6_tag_ptr(j_good)) ! really from TAGR bank no overwritten

c      ph_hit_in_time(1)=0 !defined earlier
      enddo
      endif

c----------------------------------------------------------------------
c     for photon analysis:
c     PHOTON:SEC_coin_ST <-> STPB
c----------------------------------------------------------------------

c --- initilize
      SEC_coin_ST=0
      do i=1,6
        sec(i) = 0
      enddo
      do i=1,3
        ST_pair(i) = 0
      enddo
        ST_all = 0
c --- loop over stpb, looking for coinc in ST counter
      do i=1,clasg6_st_part !number of hits in ST
      sec( clasg6_st_sector(i) ) = sec( clasg6_st_sector(i) ) + 1
      enddo  !clasg6_st_sector() is a sector number for the given hit
      if (sec(1).ne.0.or.sec(2).ne.0) ST_pair(1) =  1
      if (sec(3).ne.0.or.sec(4).ne.0) ST_pair(2) =  1
      if (sec(5).ne.0.or.sec(6).ne.0) ST_pair(3) =  1
      ST_all = ST_pair(1) + ST_pair(2) + ST_pair(3)
      sec_coin_ST = ST_all 
      

c----------------------------------------------------------------------
c     for photon analysis:
c     SC   <-> SC
c     SCPB <-> SCPB
c----------------------------------------------------------------------
      if(clasg6_nSC.gt.20) then
      nWarn3 = nWarn3 + 1
      if(nWarn1/1000*1000.eq.nWarn3 .or. nWarn3.eq.1) then
      print *,'WARNING: More than 20 hits in SC cut on 20th #',nWarn3
      endif
      endif
      nSC=min(clasg6_nSC,20)
      do i=1,nSC
       sector_SC(i)=clasg6_secSC(i)
       id_SC(i)    =clasg6_idSC(i)
       tdcl_SC(i)  =clasg6_TDCLSC(i)
       adcl_SC(i)  =clasg6_ADCLSC(i)
       tdcr_SC(i)  =clasg6_TDCRSC(i)
       adcr_SC(i)  =clasg6_ADCRSC(i)
      enddo

c----------------------------------------------------------------------
c     for photon analysis:
c     other XXPB banks are processes in the common part
c     they are equal in both el. and gamma analysis 
c----------------------------------------------------------------------

c     SEE SCPB processing for E1 

        ENDIF !IF photon analysis


c----------------------------------------------------------------------
c    EVNT bank for ELECTRONS ONLY !!!
c----------------------------------------------------------------------
      IF(beam_type) THEN      
      nEVNT = clas_gpart
      if(nEVNT.gt.10) then
      nWarn1 = nWarn1 + 1
      if(nWarn1/1000*1000.eq.nWarn1 .or. nWarn1.eq.1) then
      print *,'WARNING: More than 10 part. in EVNT cut on 10th #',nWarn1
      endif
      nEVNT = 10
      endif
      do i=1,nEVNT
      IF( clas_p(i).gt.0. .and. clas_p(i).lt.10. )THEN

      if(abs(clas_cx(i)**2+clas_cy(i)**2+clas_cz(i)**2)-1..gt.0.001)then
       th_EVNT(i)      = -1000.
       ph_EVNT(i)      = -1000.
      else
       th_EVNT(i)      = (180./pigr)*acos(clas_cz(i))
       if(clas_cx(i).ne.0) then       
       ph_EVNT(i)      = (180./pigr)*atan(clas_cy(i)/clas_cx(i))       
       else
       if(clas_cy(i).gt.0.) ph_EVNT(i) = 90.
       if(clas_cy(i).lt.0.) ph_EVNT(i) = 270.
       endif
      if(clas_cx(i).lt.0..and.clas_cy(i).gt.0.)ph_EVNT(i)=ph_EVNT(i)+180
      if(clas_cx(i).lt.0..and.clas_cy(i).lt.0.)ph_EVNT(i)=ph_EVNT(i)+180
      if(clas_cx(i).gt.0..and.clas_cy(i).lt.0.)ph_EVNT(i)=ph_EVNT(i)+360
      endif
       ID_EVNT(i)      = clas_id(i)
       ID_OUR(i)       = clas_id(i)
       Pmom_EVNT(i)    = clas_p(i)
       Pmomx_EVNT(i)   = clas_p(i)*clas_cx(i)
       Pmomy_EVNT(i)   = clas_p(i)*clas_cy(i)
       Pmomz_EVNT(i)   = clas_p(i)*clas_cz(i)
       Mass_EVNT(i)    = clas_m(i)
       Mass_OUR(i)     = clas_m(i)
       Charge_EVNT(i)  = clas_q(i)
       Beta_EVNT(i)    = clas_b(i)
       Beta_OUR(i)     = clas_b(i)
       x_EVNT(i)       = clas_vx(i)
       y_EVNT(i)       = clas_vy(i)
       z_EVNT(i)       = clas_vz(i)
       DCstat_EVNT(i)  = clas_dc(i)
       CCstat_EVNT(i)  = clas_cc(i)
       SCstat_EVNT(i)  = clas_sc(i)
       ECstat_EVNT(i)  = clas_ec(i)
       LCstat_EVNT(i)  = clas_lec(i)
       STstat_EVNT(i)  = 0
       Status_EVNT(i)  = clas_stat(i)
       if(clas_dc(i).gt.0)Trk_flag_EVNT(i) = clas_dc_stat(clas_dc(i))
       if(clas_dc(i).le.0)Trk_flag_EVNT(i) = 0
      ELSE
       nWarn2 = nWarn2 + 1
       if(nWarn1/1000*1000.eq.nWarn2 .or. nWarn2.eq.1) then
       print *,'WARNING: undefined momentum #',nWarn2
       endif
       th_EVNT(i)      = -1000.
       ph_EVNT(i)      = -1000.
       ID_EVNT(i)      = 0
       ID_OUR(i)       = 0
       Pmom_EVNT(i)    = -1000.
       Pmomx_EVNT(i)   = -1000.
       Pmomy_EVNT(i)   = -1000.
       Pmomz_EVNT(i)   = -1000.
       Mass_EVNT(i)    = -1000.
       Mass_OUR(i)     = -1000.
       Charge_EVNT(i)  = -1000.
       Beta_EVNT(i)    = -1000.
       Beta_OUR(i)     = -1000.
       x_EVNT(i)       = -1000.
       y_EVNT(i)       = -1000.
       z_EVNT(i)       = -1000.
       DCstat_EVNT(i)  = 0
       CCstat_EVNT(i)  = 0
       SCstat_EVNT(i)  = 0
       ECstat_EVNT(i)  = 0
       LCstat_EVNT(i)  = 0
       STstat_EVNT(i)  = 0
       Status_EVNT(i)  = 0
       Trk_flag_EVNT(i)= 0
      ENDIF
      enddo
      ENDIF !IF(beam_type) THEN      



c----------------------------------------------------------------------
c    EVNT bank for PHOTONS ONLY !!!
c----------------------------------------------------------------------
      IF(.not.beam_type) THEN
      nEVNT = clasg6_gpart
      if(nEVNT.gt.9) then
      nWarn1 = nWarn1 + 1
      if(nWarn1/1000*1000.eq.nWarn1 .or. nWarn1.eq.1) then
      print *,'WARNING: More than 9 part. in G6EVNT cut on 9th #',nWarn1
      endif
      nEVNT = 9
      endif
      do i=1,nEVNT
      IF( clasg6_p(i).gt.0. .and. clasg6_p(i).lt.10. )THEN

      if(abs(clasg6_cx(i)**2+clasg6_cy(i)**2+clasg6_cz(i)**2)-1..gt.0.001)then
       th_EVNT(i)      = -1000.
       ph_EVNT(i)      = -1000.
      else
       th_EVNT(i)      = (180./pigr)*acos(clasg6_cz(i))
       if(clasg6_cx(i).ne.0) then       
       ph_EVNT(i)      = (180./pigr)*atan(clasg6_cy(i)/clasg6_cx(i))       
       else
       if(clasg6_cy(i).gt.0.) ph_EVNT(i) = 90.
       if(clasg6_cy(i).lt.0.) ph_EVNT(i) = 270.
       endif
      if(clasg6_cx(i).lt.0..and.clasg6_cy(i).gt.0.)ph_EVNT(i)=ph_EVNT(i)+180
      if(clasg6_cx(i).lt.0..and.clasg6_cy(i).lt.0.)ph_EVNT(i)=ph_EVNT(i)+180
      if(clasg6_cx(i).gt.0..and.clasg6_cy(i).lt.0.)ph_EVNT(i)=ph_EVNT(i)+360
      endif
       ID_EVNT(i)      = clasg6_id(i)
       ID_OUR(i)       = clasg6_id(i)
       Pmom_EVNT(i)    = clasg6_p(i)
       Pmomx_EVNT(i)   = clasg6_p(i)*clasg6_cx(i)
       Pmomy_EVNT(i)   = clasg6_p(i)*clasg6_cy(i)
       Pmomz_EVNT(i)   = clasg6_p(i)*clasg6_cz(i)
       Mass_EVNT(i)    = clasg6_m(i)
       Mass_OUR(i)     = clasg6_m(i)
       Charge_EVNT(i)  = clasg6_q(i)
       Beta_EVNT(i)    = clasg6_b(i)
       Beta_OUR(i)     = clasg6_b(i)
       x_EVNT(i)       = clasg6_vx(i)
       y_EVNT(i)       = clasg6_vy(i)
       z_EVNT(i)       = clasg6_vz(i)
       DCstat_EVNT(i)  = clasg6_dc(i)
       CCstat_EVNT(i)  = clasg6_cc(i)
       SCstat_EVNT(i)  = clasg6_sc(i)
       ECstat_EVNT(i)  = clasg6_ec(i)
       LCstat_EVNT(i)  = clasg6_lec(i)
       STstat_EVNT(i)  = clasg6_st(i)
       Status_EVNT(i)  = clasg6_stat(i)
       if(clasg6_dc(i).gt.0)Trk_flag_EVNT(i) = clas_dc_stat(clasg6_dc(i))
       if(clasg6_dc(i).le.0)Trk_flag_EVNT(i) = 0
      ELSE
       nWarn2 = nWarn2 + 1
       if(nWarn1/1000*1000.eq.nWarn2 .or. nWarn2.eq.1) then
       print *,'WARNING: undefined momentum #',nWarn2
       endif
       th_EVNT(i)      = -1000.
       ph_EVNT(i)      = -1000.
       ID_EVNT(i)      = 0
       ID_OUR(i)       = 0
       Pmom_EVNT(i)    = -1000.
       Pmomx_EVNT(i)   = -1000.
       Pmomy_EVNT(i)   = -1000.
       Pmomz_EVNT(i)   = -1000.
       Mass_EVNT(i)    = -1000.
       Mass_OUR(i)     = -1000.
       Charge_EVNT(i)  = -1000.
       Beta_EVNT(i)    = -1000.
       Beta_OUR(i)     = -1000.
       x_EVNT(i)       = -1000.
       y_EVNT(i)       = -1000.
       z_EVNT(i)       = -1000.
       DCstat_EVNT(i)  = 0
       CCstat_EVNT(i)  = 0
       SCstat_EVNT(i)  = 0
       ECstat_EVNT(i)  = 0
       LCstat_EVNT(i)  = 0
       STstat_EVNT(i)  = 0
       Status_EVNT(i)  = 0
       Trk_flag_EVNT(i)= 0
      ENDIF
      enddo
      ENDIF !IF(beam_type) THEN      

c----------------------------------------------------------------------
c    shift EVNT bank for PHOTONS ONLY !!!
c----------------------------------------------------------------------
      IF(.not.beam_type) THEN      
      if(.not.flag_false_g.and.nEVNT.gt.0) then
         do i=1,nEVNT
          ID_EVNT    (nEVNT+2-i) =     ID_EVNT(nEVNT+1-i)
          Pmom_EVNT  (nEVNT+2-i) =   Pmom_EVNT(nEVNT+1-i)
          Pmomx_EVNT (nEVNT+2-i) =  Pmomx_EVNT(nEVNT+1-i)
          Pmomy_EVNT (nEVNT+2-i) =  Pmomy_EVNT(nEVNT+1-i)
          Pmomz_EVNT (nEVNT+2-i) =  Pmomz_EVNT(nEVNT+1-i)
          th_EVNT    (nEVNT+2-i) =     th_EVNT(nEVNT+1-i)
          ph_EVNT    (nEVNT+2-i) =     ph_EVNT(nEVNT+1-i)
          ID_OUR     (nEVNT+2-i) =      ID_OUR(nEVNT+1-i)
          Mass_OUR   (nEVNT+2-i) =    Mass_OUR(nEVNT+1-i)
          Beta_OUR   (nEVNT+2-i) =    Beta_OUR(nEVNT+1-i)
          Mass_EVNT  (nEVNT+2-i) =   Mass_EVNT(nEVNT+1-i)
          Charge_EVNT(nEVNT+2-i) = Charge_EVNT(nEVNT+1-i)
          Beta_EVNT  (nEVNT+2-i) =   Beta_EVNT(nEVNT+1-i)
          X_EVNT     (nEVNT+2-i) =      X_EVNT(nEVNT+1-i)
          Y_EVNT     (nEVNT+2-i) =      Y_EVNT(nEVNT+1-i)
          Z_EVNT     (nEVNT+2-i) =      Z_EVNT(nEVNT+1-i)
          DCstat_EVNT(nEVNT+2-i) = DCstat_EVNT(nEVNT+1-i)
          CCstat_EVNT(nEVNT+2-i) = CCstat_EVNT(nEVNT+1-i)
          SCstat_EVNT(nEVNT+2-i) = SCstat_EVNT(nEVNT+1-i)
          ECstat_EVNT(nEVNT+2-i) = ECstat_EVNT(nEVNT+1-i)
          LCstat_EVNT(nEVNT+2-i) = LCstat_EVNT(nEVNT+1-i)
          STstat_EVNT(nEVNT+2-i) = STstat_EVNT(nEVNT+1-i)
          Status_EVNT(nEVNT+2-i) = Status_EVNT(nEVNT+1-i)
         enddo
          nEVNT = nEVNT + 1   ! if last is not 0 , adding a new record 

      endif
      ENDIF !IF(beam_type) THEN      






c----------------------------------------------------------------------
c    MVRT Vertex bank is not used. Use EVNT vertex instead 
c----------------------------------------------------------------------
      ntrk_MVRT = -1000 ! use vertex coordinates from EVNT
      x_MVRT    = -1000
      y_MVRT    = -1000
      z_MVRT    = -1000


c----------------------------------------------------------------------
c    SCPB
c----------------------------------------------------------------------
      do i = 1,6   !sector
      do jj = 1,50 !id
         hit_sc2SCPB(i,jj) = 0
      enddo
      enddo
      nSCPB = clas_sc_part
      if(nSCPB.gt.20) then
      print *,' INFO: More than 20 particles in SCPB cut on 20th'
      nSCPB = 20
      endif
      do i=1,nSCPB
       SCstat_SCPB(i) = clas_sc_stat(i)
       ScPdHt_SCPB(i) = 10000*clas_sc_sect(i) + 100*clas_sc_pd(i)
     &                 +clas_sc_hit(i)
       Time_SCPB(i)   = clas_sc_t(i)
       Path_SCPB(i)   = clas_sc_r(i)
         Status_SCPB(i) = clas_sc_stat(i)  ! no in the final ntupl
         Edep_SCPB(i)   = clas_edep(i)     !
         Chi2SC_SCPB(i) = clas_sc_c2(i)    !
         sector_SCPB(i) = clas_sc_sect(i)  !
         pd_id_SCPB(i)  = clas_sc_pd(i)    !
         hit_id_SCPB(i) = clas_sc_hit(i)   !
c --- Matching SC bank hits
      SCstat_SCPB(i) = 0
      do jj=1,nSC
        if(adcl_SC(jj).ne.0.and.adcr_SC(jj).ne.0
     &  .and.tdcl_SC(jj).ne.0.and.tdcr_SC(jj).ne.0) then ! 1 request for matching: SC adc tdc RL <> 0
        if(id_sc(jj).eq.pd_id_SCPB(i)
     &  .and.sector_SC(jj).eq.sector_SCPB(i)) then       ! 2 reqest: id sec matching
        if(hit_sc2SCPB(sector_SC(jj),id_sc(jj)).eq.0) then
           hit_sc2SCPB(sector_SC(jj),id_sc(jj)) = 1 ! this sec/id was hit 
           SCstat_SCPB(i) = jj 
        else
           SCstat_SCPB(i) = -1 ! if multiple hits in SC no match
        endif
        endif
        endif
      enddo

      enddo


c----------------------------------------------------------------------
c    DCPB
c----------------------------------------------------------------------
      nDCPB = clas_dc_part
      if(nDCPB.gt.20) then
      print *,' INFO: More than 20 particles in DCPB cut on 20th'
      nDCPB = 10
      endif
      do i=1,nDCPB
       ScTr_DCPB(i)  = 100*clas_dc_sect(i) + clas_dc_trk(i)
       XSC_DCPB(i)   = clas_dc_xsc(i)
       YSC_DCPB(i)   = clas_dc_ysc(i)
       ZSC_DCPB(i)   = clas_dc_zsc(i)
       CXSC_DCPB(i)  = clas_dc_cxsc(i)
       CYSC_DCPB(i)  = clas_dc_cysc(i)
       CZSC_DCPB(i)  = clas_dc_czsc(i)
       STATUS_DCPB(i)= clas_dc_stat(i)
      enddo


c----------------------------------------------------------------------
c    CCPB
c----------------------------------------------------------------------
      nCCPB = clas_cc_part
      if(nCCPB.gt.20) then
      print *,' INFO: More than 20 particles in CCPB cut on 20th'
      nCCPB = 10
      endif
      do i=1,nCCPB
       ScSgHt_CCPB(i) = 1000*clas_cc_sect(i) + 100*clas_cc_segm(i) 
     &                  + clas_cc_hit(i)
       NPHE_CCPB(i)   = clas_nphe(i)
       STATUS_CCPB(i) = -1000
      enddo


c----------------------------------------------------------------------
c    ECPB
c----------------------------------------------------------------------
      nECPB=clas_ec_part
      if(nECPB.gt.20) then
      print *,' INFO: More than 20 particles in ECPB cut on 20th'
      nECPB = 10
      endif
      do i=1,nECPB
       ETOT_ECPB(i)   = clas_etot(i)
       EIN_ECPB(i)    = clas_ec_ei(i)
       EOUT_ECPB(i)   = clas_ec_eo(i)
       STATUS_ECPB(i) = clas_ec_stat(i)
       M2_hit_ECPB(i) = clas_ec_m2(i)
       M3_hit_ECPB(i) = clas_ec_m3(i)
       M4_hit_ECPB(i) = clas_ec_m4(i)
      enddo


c      IF( LC_analysis ) THEN
c      nLCPB = clas_lac_part
c      do i=1,nLCPB
c       ScHt_LCPB(i) = 100*clas_lec_sect(i) + clas_lec_hit(i)
c       ETOT_LCPB(i) = clas_lec_etot(i)
c       X_LCPB(i)    = clas_lec_x(i)
c       Y_LCPB(i)    = clas_lec_y(i)
c       Z_LCPB(i)    = clas_lec_z(i)
c      enddo
c      nEC1R=0
c      do i=1,nEC1R
c       sec_EC1R(i)     =0
c       E_tot_EC1R(i)   =0
c       x_m_EC1R(i)     =0
c       y_m_EC1R(i)     =0
c       z_m_EC1R(i)     =0
c       x_in_EC1R(i)    =0
c       y_in_EC1R(i)    =0
c       ibl_EC1R(i)     =0
c       ncluster_EC1R(i)=0
c      enddo
c      ENDIF


c      IF( LC_analysis ) THEN 
c       nTDPL=0
c       do i=1,nTDPL
c       trkins_TDPL(i)
c       XSC_TDPL(i)
c       YSC_TDPL(i)
c       ZSC_TDPL(i)
c       CXSC_TDPL(i)  
c       CYSC_TDPL(i) 
c       CZSC_TDPL(i)  
c       XLC_TDPL(i)
c       YLC_TDPL(i)
c       ZLC_TDPL(i)
c       CXLC_TDPL(i)  
c       CYLC_TDPL(i)  
c       CZLC_TDPL(i)  
c      ENDIF


      RETURN
      END



