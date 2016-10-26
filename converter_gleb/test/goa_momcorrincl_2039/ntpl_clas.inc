c---------------------------------------------------------------------
c                           CLAS Ntuple
c---------------------------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  E1 and COMMON ntuples
c                  HEVT
c                  EVNT
c                  ??PB
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL
     +  !HEVT
     +  clas_q_l,clas_t_l,clas_tr_time,clas_rf_time1,clas_rf_time2
     +  !EVNT
     + ,clas_p(20),clas_m(20),clas_b(20)
     + ,clas_cx(20),clas_cy(20),clas_cz(20)
     + ,clas_vx(20),clas_vy(20),clas_vz(20)
     +  !DCPB
     + ,clas_dc_xsc(20),clas_dc_ysc(20),clas_dc_zsc(20)
     + ,clas_dc_cxsc(20),clas_dc_cysc(20),clas_dc_czsc(20)
     + ,clas_dc_xec(20),clas_dc_yec(20),clas_dc_zec(20)
     + ,clas_dc_thcc(20),clas_dc_c2(20)
     +  !ECPB
     + ,clas_etot(20),clas_ec_ei(20),clas_ec_eo(20),clas_ec_t(20)
     + ,clas_ec_r(20),clas_ech_x(20),clas_ech_y(20),clas_ech_z(20)
     + ,clas_ec_m2(20),clas_ec_m3(20),clas_ec_m4(20)
     + ,clas_ec_c2(20),clas_edep(20)
     +  !SCPB
     + ,clas_sc_t(20),clas_sc_r(20),clas_sc_c2(20)
     +  !CCPB
     + ,clas_cc_t(20),clas_cc_r(20),clas_cc_c2(20)
     +  !LCPB
     + ,clas_lec_etot(20),clas_lec_t(20),clas_lec_r(20)
     + ,clas_lec_x(20),clas_lec_y(20),clas_lec_z(20),clas_lec_c2(20)

      INTEGER
     +  !HEVT 
     +  clas_npart,clas_evstat,clas_evntid,clas_evtype,clas_evntclas
     +  !EVNT
     + ,clas_gpart,clas_id(20),clas_stat(20)
     + ,clas_dc(20),clas_cc(20),clas_sc(20),clas_ec(20),clas_lec(20)
     + ,clas_q(20)
     +  !DCPB
     + ,clas_dc_part,clas_dc_sect(20)
     + ,clas_dc_trk(20),clas_dc_stat(20)
     +  !ECPB
     + ,clas_ec_part,clas_ec_stat(20),clas_ec_sect(20)
     + ,clas_ec_whol(20),clas_ec_inst(20),clas_ec_oust(20)
     +  !SCPB
     + ,clas_sc_part,clas_sc_sect(20)
     + ,clas_sc_hit(20),clas_sc_pd(20),clas_sc_stat(20)
     +  !CCPB
     + ,clas_cc_part,clas_cc_sect(20),clas_cc_hit(20)
     + ,clas_cc_segm(20),clas_nphe(20)
     +  !LCPB
     + ,clas_lac_part,clas_lec_sect(20),clas_lec_hit(20)
     + ,clas_lec_stat(20)


c --- E1 ntuples
      COMMON /CLAS_HEVT/ 
     +  clas_NPART,clas_evstat,clas_evntid,clas_evtype,clas_evntclas,clas_q_l,clas_t_l
     + ,clas_tr_time,clas_rf_time1,clas_rf_time2
      COMMON /CLAS_EVNT/
     +  clas_GPART,clas_id,clas_stat,clas_dc,clas_cc,clas_sc,clas_ec,clas_lec,clas_p,clas_m,clas_q,clas_b
     + ,clas_cx,clas_cy,clas_cz,clas_vx,clas_vy,clas_vz
      COMMON /CLAS_DCPB/
     +  clas_DC_PART,clas_dc_sect,clas_dc_trk,clas_dc_stat,clas_dc_xsc,clas_dc_ysc
     + ,clas_dc_zsc,clas_dc_cxsc,clas_dc_cysc,clas_dc_czsc,clas_dc_xec,clas_dc_yec,clas_dc_zec,clas_dc_thcc
     + ,clas_dc_c2
      COMMON /CLAS_ECPB/
     +  clas_EC_PART,clas_ec_stat,clas_ec_sect,clas_ec_whol,clas_ec_inst,clas_ec_oust,clas_etot,clas_ec_ei
     + ,clas_ec_eo,clas_ec_t,clas_ec_r,clas_ech_x,clas_ech_y,clas_ech_z,clas_ec_m2,clas_ec_m3,clas_ec_m4,clas_ec_c2
      COMMON /CLAS_SCPB/
     +  clas_SC_PART,clas_sc_sect,clas_sc_hit,clas_sc_pd,clas_sc_stat,clas_edep,clas_sc_t,clas_sc_r,clas_sc_c2
      COMMON /CLAS_CCPB/
     +  clas_CC_PART,clas_cc_sect,clas_cc_hit,clas_cc_segm,clas_nphe,clas_cc_t,clas_cc_r,clas_cc_c2
      COMMON /CLAS_LCPB/
     +  clas_LAC_PART,clas_lec_sect,clas_lec_hit,clas_lec_stat,clas_lec_etot,clas_lec_t,clas_lec_r,clas_lec_x,clas_lec_y,clas_lec_z
     + ,clas_lec_c2

 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  G6 ntuples
c                  PART - not used
c                  EVNT - different from E1 (the variable 'st')
c                  STPB, TGPB
c                  TAGR
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL
     +  !PART
     +  clas_xPART(10),clas_yPART(10),clas_zPART(10)
     + ,clas_ePART(10),clas_pxPART(10),clas_pyPART(10),clas_pzPART(10)
     +  !EVNT G6
     + ,clasg6_p(20),clasg6_m(20),clasg6_b(20)
     + ,clasg6_cx(20),clasg6_cy(20),clasg6_cz(20)
     + ,clasg6_vx(20),clasg6_vy(20),clasg6_vz(20)
     +  !STPB G6
     + ,clasg6_st_time(20),clasg6_st_rtrk(20)
     +  !TAGR G6
     + ,clasg6_E_gamma(20),clasg6_T_gamma_noRF(20),clasg6_T_gamma(20)
     +  !TGPB G6
     + ,clasg6_vertex_time(20),clasg6_tag_energy(20),clasg6_dt_st_tag(20)

      INTEGER 
     +  !PART
     +  clas_nPRT,clas_pidPART(10),clas_qPART(10)
     + ,clas_flagsPART(10),clas_iPART(10)
     +  !SC G6
     + ,clasg6_nSC,clasg6_secSC(288),clasg6_idSC(288)
     + ,clasg6_TDCLSC(288),clasg6_ADCLSC(288)
     + ,clasg6_TDCRSC(288),clasg6_ADCRSC(288)
     +  !EVNT G6
     + ,clasg6_gpart,clasg6_id(20),clasg6_stat(20)
     + ,clasg6_dc(20),clasg6_cc(20),clasg6_sc(20),clasg6_ec(20),clasg6_lec(20),clasg6_st(20)
     + ,clasg6_q(20)
     +  !STPB G6
     + ,clasg6_st_part,clasg6_st_sector(20),clasg6_st_ihit(20)
     + ,clasg6_st_trkno(20),clasg6_st_status(20)
     +  !TAGR G6
     + ,clasg6_taghit,clasg6_tagstat(20),clasg6_Tid(20),clasg6_Eid(20)
     +  !TGPB G6
     + ,clasg6_taggoodhit,clasg6_tag_ptr(20)


      COMMON /CLAS_PART/ 
     +  clas_nPRT,clas_pidPART,clas_xPART,clas_yPART,clas_zPART
     + ,clas_ePART,clas_pxPART,clas_pyPART,clas_pzPART
     + ,clas_qPART,clas_flagsPART,clas_iPART
      COMMON /CLASG6_SC/ 
     +  clasg6_nSC,clasg6_secSC,clasg6_idSC
     + ,clasg6_TDCLSC,clasg6_ADCLSC,clasg6_TDCRSC,clasg6_ADCRSC
      COMMON /CLASG6_EVNT/ 
     +  clasg6_gpart,clasg6_id,clasg6_stat
     + ,clasg6_dc,clasg6_cc,clasg6_sc,clasg6_ec,clasg6_lec,clasg6_st
     + ,clasg6_p,clasg6_m,clasg6_q,clasg6_b
     + ,clasg6_cx,clasg6_cy,clasg6_cz,clasg6_vx,clasg6_vy,clasg6_vz
      COMMON /CLASG6_STPB/ 
     +  clasg6_st_part,clasg6_st_sector,clasg6_st_ihit,clasg6_st_trkno
     + ,clasg6_st_time,clasg6_st_rtrk,clasg6_st_status
      COMMON /CLASG6_TAGR/ 
     +  clasg6_taghit,clasg6_E_gamma,clasg6_T_gamma_noRF,clasg6_T_gamma
     + ,clasg6_tagstat
     + ,clasg6_Tid,clasg6_Eid
      COMMON /CLASG6_TGPB/ 
     +  clasg6_taggoodhit,clasg6_tag_ptr
     + ,clasg6_vertex_time,clasg6_tag_energy,clasg6_dt_st_tag

