

      SUBROUTINE book_ntpl_goa      
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INTEGER*4 NN,i,n_sec,n_pad
      INTEGER Iblock_all,n_e_elast,n_e_incl
      CHARACTER*1 ch1
      real*4  elast
      COMMON   /qqq/Iblock_all,n_e_elast,elast,n_e_incl
c------------------------------------------------------
c Genral RUNDATA and Ngamma ntuple 51 and histograms
c------------------------------------------------------
      call hbnt  (51,'RUNDATA',' ')

      call hbname(51,'RUNNUMER', RunNumber,     'RunNumber:I')

      call hbname(51,'RUNGEN', NEventFRUN11,    'NEventFRUN11:I')
      call hbname(51,'RUNGEN', NEventFileFirst, 'NEventFileFirst:I')
      call hbname(51,'RUNGEN', NEventFileLast,  'NEventFileLast:I')
      call hbname(51,'RUNGEN', NEventTotalRun,  'NEventTotalRun:I')
      call hbname(51,'RUNGEN', NEventThisFile,  'NEventThisFile:I')

      call hbname(51,'RUNGAMMA', NGammaE,    'NGammaE(767):R')
      call hbname(51,'RUNGAMMA', NGammaE_2,  'NGammaE_2(767):R')
      call hbname(51,'RUNGAMMA', NGammaT,    'NGammaT(61):R')
      call hbname(51,'RUNGAMMA', NGammaT_2,  'NGammaT_2(61):R')

        if(.not.beam_type) then
      CALL HBOOK1(11002,'NGammaE',             767, 0., 767., 0.0)
      CALL HBOOK1(11003,'NGammaE_2',           767, 0., 767., 0.0)
      CALL HBOOK1(11004,'NGammaT',             61,  0.,  61., 0.0)
      CALL HBOOK1(11005,'NGammaT_2',           61,  0.,  61., 0.0)
      CALL HBOOK1(11006,'NGamma in Energy intervals',   800,  0.,  61., 0.0)

      CALL HBOOK1(11007,'NeventsE/NgammaE All files',    767, 0., 767., 0.0)
      CALL HBOOK1(11008,'NeventsT/NgammaT All files',    61,  0.,  61., 0.0)
      CALL HBOOK1(11009,'NeventsE/NgammaE 1 file',       767, 0., 767., 0.0)
      CALL HBOOK1(11010,'NeventsE',                      767, 0., 767., 0.0)
      CALL HBARX (11007)
      CALL HBARX (11008)
      CALL HBARX (11009)
      CALL HBARX (11010)
        endif


c---------------------------------------------------------------------
c  General DATA and MC ntuples 60 and 61
c---------------------------------------------------------------------

      DO NN = 60,61
      IF(ana_gen.eq.1) THEN

      call hbnt(NN,'PART',' ')

c ----- general -----
      call hbname(NN,'GENERAL',n_ev_run,   'n_ev_run:I')
      call hbname(NN,'GENERAL',trig_type,  'trig_type:I')
      call hbname(NN,'GENERAL',trig_clas,  'trig_clas:I')
      call hbname(NN,'GENERAL',fcgdose,    'fcgdose:R')
      call hbname(NN,'GENERAL',event_clas, 'event_clas:I')
              IF(hel_ana.eq.1) THEN
      call hbname(NN,'GENERAL',beam_hel,   'beam_hel:I')
      call hbname(NN,'GENERAL',beam_helkey,'beam_helkey:I')
      call hbname(NN,'GENERAL',beam_keypair,'beam_keypair:I')
      call hbname(NN,'GENERAL',beam_keyboth,'beam_keyboth:I')
      call hbname(NN,'GENERAL',doplus_BPM, 'doplus_BPM:I')
      call hbname(NN,'GENERAL',dominus_BPM, 'dminus_BPM:I')
      call hbname(NN,'GENERAL',doplus_FC,  'doplus_FC:I')
      call hbname(NN,'GENERAL',dominus_FC, 'dominus_FC:I')
        ENDIF
      call hbname(NN,'GENERAL',n_ev_nt,    'n_ev_nt:I')
      call hbname(NN,'GENERAL',n_part,     'n_part:I')
      call hbname(NN,'GENERAL',ntot,       'ntot:I')
      call hbname(NN,'GENERAL',bit,        'bit:I')
      call hbname(NN,'GENERAL',sec_coin_SC, 'sec_SC:I')
      call hbname(NN,'GENERAL',sec_coin_ST, 'sec_ST:I')
      call hbname(NN,'GENERAL',Iblock_all,  'Iblock_all:I')
      call hbname(NN,'GENERAL',n_e_incl,  'n_e_incl:I')
      call hbname(NN,'GENERAL',elast,  'elast:R')
      call hbname(NN,'GENERAL',tg_hevt,     'tg_hevt:R')
        IF(beam_type) THEN
*      call hbname(NN,'GENERAL',fcg_hevt, 'fcg_hevt:R')
*      call hbname(NN,'GENERAL',Qgated, 'Qgated:R')
*      call hbname(NN,'GENERAL',Qgated_plus,  'Qgated_plus:R')
        ENDIF

        IF(beam_type) THEN
      call hbname(NN,'ELECTRON', Omega,      'omega:R')      
      call hbname(NN,'ELECTRON', Q2,         'Q2:R')
      call hbname(NN,'ELECTRON', W,          'W:R')
      call hbname(NN,'ELECTRON', eps,        'eps:R')
      call hbname(NN,'ELECTRON', eps_l,      'eps_l:R')
      call hbname(NN,'ELECTRON', Flux_vp,    'Flux_vp:R')
      call hbname(NN,'ELECTRON', E_EL,       'E_EL:R')
      call hbname(NN,'ELECTRON', P_EL,       'P_EL:R')
      call hbname(NN,'ELECTRON', B_EL,       ' B_EL:R')
      call hbname(NN,'ELECTRON', th_EL,      'th_EL:R')      
      call hbname(NN,'ELECTRON', ph_EL,      'ph_EL:R')
      call hbname(NN,'ELECTRON', acc_EL,     'acc_EL:R')
      call hbname(NN,'ELECTRON', X_EL,       'X_EL:R')
      call hbname(NN,'ELECTRON', Y_EL,       'Y_EL:R')
      call hbname(NN,'ELECTRON', Z_EL,       'Z_EL:R')
      call hbname(NN,'ELECTRON', XMVRT_EL,   'XMVRT_EL:R')
      call hbname(NN,'ELECTRON', YMVRT_EL,   'YMVRT_EL:R')
      call hbname(NN,'ELECTRON', ZMVRT_EL,   'ZMVRT_EL:R')
      call hbname(NN,'ELECTRON', ECtot_EL,   'ECtot_EL:R')
      call hbname(NN,'ELECTRON', ECin_EL,    'ECin_EL:R')
      call hbname(NN,'ELECTRON', ECout_EL,   'ECout_EL:R')
      call hbname(NN,'ELECTRON', NpheCC_EL,  'NpheCC_EL:R')
      call hbname(NN,'ELECTRON', EffCC_EL,   'EffCC_EL:R')
      call hbname(NN,'ELECTRON', statEVNT_EL,'statEVNT_EL:I')
      call hbname(NN,'ELECTRON', trkEVNT_EL, 'trkEVNT_EL:I')
      call hbname(NN,'ELECTRON', PdHit_EL,   'PdHit_EL:I')
        ELSE
      call hbname(NN,'PH_BEAM',Egamma,     'Egamma:R') 
      call hbname(NN,'PH_BEAM',W,          'W:R')
      call hbname(NN,'PH_BEAM',Tag_time,   'Tag_time:R')
      call hbname(NN,'PH_BEAM',ST_time,    'ST_time:R')
      call hbname(NN,'PH_BEAM',RF1_time,   'RF1_time:R')
      call hbname(NN,'PH_BEAM',RF2_time,   'RF2_time:R')
      call hbname(NN,'PH_BEAM',hit_in_time,'hit_in_time:I')
      call hbname(NN,'PH_BEAM',chan,       'chanMC:I')
      call hbname(NN,'PH_BEAM',true,       'rand:I')
        ENDIF

      call hbname(NN,'NUCLEON',nNU,    'nNU[0,1]:I')
      call hbname(NN,'NUCLEON',E_NU,   'E_NU(nNU):R')
      call hbname(NN,'NUCLEON',P_NU,   'P_NU(nNU):R')
      call hbname(NN,'NUCLEON',B_NU,   'B_NU(nNU):R')
      call hbname(NN,'NUCLEON',M2_NU,  'M2_NU(nNU):R')
      call hbname(NN,'NUCLEON',th_NU,  'th_NU(nNU):R')
      call hbname(NN,'NUCLEON',ph_NU,  'ph_NU(nNU):R')
      call hbname(NN,'NUCLEON',q_NU,   'q_NU(nNU):R')
      call hbname(NN,'NUCLEON',acc_NU, 'acc_NU(nNU):R')
      call hbname(NN,'NUCLEON',x_NU,   'x_NU(nNU):R')
      call hbname(NN,'NUCLEON',y_NU,   'y_NU(nNU):R')
      call hbname(NN,'NUCLEON',z_NU,   'z_NU(nNU):R')
      call hbname(NN,'NUCLEON',statEVNT_NU,'statEVNT_NU(nNU):I')
      call hbname(NN,'NUCLEON',trkEVNT_NU, 'trkEVNT_NU(nNU):I')
      if (.not.beam_type) then 
      call hbname(NN,'NUCLEON',thrSC_NU,   'thrSC_NU(nNU):R')
      endif
      if (beam_type) then  
      call hbname(NN,'NUCLEON',th1_NU,     'th1_NU(nNU):R')
      call hbname(NN,'NUCLEON',th2_NU,     'th2_NU(nNU):R')
      endif
      call hbname(NN,'NUCLEON',PdHit_NU,   'PdHit_NU(nNU):I')

      call hbname(NN,'DEUTERON',nD,        'nD[0,5]:I')
      call hbname(NN,'DEUTERON',E_D,       'E_D(nD):R')
      call hbname(NN,'DEUTERON',P_D,       'P_D(nD):R')
      call hbname(NN,'DEUTERON',B_D,       'B_D(nD):R')
      call hbname(NN,'DEUTERON',M2_D,      'M2_D(nD):R')
      call hbname(NN,'DEUTERON',th_D,      'th_D(nD):R')
      call hbname(NN,'DEUTERON',ph_D,      'ph_D(nD):R')
      call hbname(NN,'DEUTERON',q_D,       'q_D(nD):R')
      call hbname(NN,'DEUTERON',acc_D,     'acc_D(nD):R')
      call hbname(NN,'DEUTERON',x_D,       'x_D(nD):R')
      call hbname(NN,'DEUTERON',y_D,       'y_D(nD):R')
      call hbname(NN,'DEUTERON',z_D,       'z_D(nD):R')
      call hbname(NN,'DEUTERON',statEVNT_D,'statEVNT_D(nD):I')
      call hbname(NN,'DEUTERON',trkEVNT_D, 'trkEVNT_D(nD):I')
      call hbname(NN,'DEUTERON',PdHit_D,   'PdHit_D(nD):I')

      call hbname(NN,'PION',npi_plus,   'npi_plus:I')
      call hbname(NN,'PION',npi_minus,  'npi_minus:I')
      call hbname(NN,'PION',nPI,        'nPI[0,10]:I')
      call hbname(NN,'PION',E_PI,       'E_PI(nPI):R')
      call hbname(NN,'PION',P_PI,       'P_PI(nPI):R')
      call hbname(NN,'PION',B_PI,       'B_PI(nPI):R')
      call hbname(NN,'PION',M2_PI,      'M2_PI(nPI):R')
      call hbname(NN,'PION',th_PI,      'th_PI(nPI):R')
      call hbname(NN,'PION',ph_PI,      'ph_PI(nPI):R')
      call hbname(NN,'PION',q_PI,       'q_PI(nPI):R')
      call hbname(NN,'PION',acc_PI,     'acc_PI(nPI):R')
      call hbname(NN,'PION',x_PI,       'x_PI(nPI):R')
      call hbname(NN,'PION',y_PI,       'y_PI(nPI):R')
      call hbname(NN,'PION',z_PI,       'z_PI(nPI):R')
      call hbname(NN,'PION',ECtot_PI,   'ECtot_PI(nPI):R')
      call hbname(NN,'PION',ECin_PI,    'ECin_PI(nPI):R')
      call hbname(NN,'PION',ECout_PI,   'ECout_PI(nPI):R')
      call hbname(NN,'PION',statEVNT_PI,'statEVNT_PI(nPI):I')
      call hbname(NN,'PION',trkEVNT_PI, 'trkEVNT_PI(nPI):I')
      
      
            if (.not.beam_type) then
      call hbname(NN,'PION',thrSC_PI,   'thrSC_PI(nPI):R')
      endif
      call hbname(NN,'PION',PdHit_PI,   'PdHit_PI(nPI):I')

      call hbname(NN,'KAON',nK,        'nK[0,10]:I')
      call hbname(NN,'KAON',E_K,       'E_K(nK):R')
      call hbname(NN,'KAON',P_K,       'P_K(nK):R')
      call hbname(NN,'KAON',B_K,       'B_K(nK):R')
      call hbname(NN,'KAON',M2_K,      'M2_K(nK):R')
      call hbname(NN,'KAON',th_K,      'th_K(nK):R')
      call hbname(NN,'KAON',ph_K,      'ph_K(nK):R')
      call hbname(NN,'KAON',q_K,       'q_K(nK):R')
      call hbname(NN,'KAON',acc_K,     'acc_K(nK):R')
      call hbname(NN,'KAON',x_K,       'x_K(nK):R')
      call hbname(NN,'KAON',y_K,       'y_K(nK):R')
      call hbname(NN,'KAON',z_K,       'z_K(nK):R')
      call hbname(NN,'KAON',statEVNT_K,'statEVNT_K(nK):I')
      call hbname(NN,'KAON',trkEVNT_K, 'trkEVNT_K(nK):I')

      call hbname(NN,'GAMMA',nG,      'nG[0,10]:I')
      call hbname(NN,'GAMMA',E_G,     'E_G(nG):R')
      call hbname(NN,'GAMMA',th_G,    'th_G(nG):R')
      call hbname(NN,'GAMMA',ph_G,    'ph_G(nG):R')
      call hbname(NN,'GAMMA',acc_G,   'acc_G(nG):R')
      call hbname(NN,'GAMMA',x_G,     'x_G(nG):R')
      call hbname(NN,'GAMMA',y_G,     'y_G(nG):R')
      call hbname(NN,'GAMMA',z_G,     'z_G(nG):R')
      call hbname(NN,'GAMMA',PdHit_G, 'PdHit_G(nG):I')

c      call hbname(NN,'LEPTON',nL,   'nL[0,4]:I')
c      call hbname(NN,'LEPTON',id_L, 'id_L(nL):R')
c      call hbname(NN,'LEPTON',E_L,  'E_L(nL):R')
c      call hbname(NN,'LEPTON',P_L,  'P_L(nL):R')
c      call hbname(NN,'LEPTON',B_L,  'B_L(nL):R')
c      call hbname(NN,'LEPTON',th_L, 'th_L(nL):R')
c      call hbname(NN,'LEPTON',ph_L, 'ph_L(nL):R')
c      call hbname(NN,'LEPTON',q_L,  'q_L(nL):R')
c      call hbname(NN,'LEPTON',acc_L,'acc_L(nL):R')
c      call hbname(NN,'LEPTON',x_L,  'x_L(nL):R')
c      call hbname(NN,'LEPTON',y_L,  'y_L(nL):R')
c      call hbname(NN,'LEPTON',z_L,  'z_L(nL):R')

       call hbname(NN,'GAGA',nGG,   'nGG[0,6]:I')
       call hbname(NN,'GAGA',W2_GG, 'W2_GG(nGG):R')
       call hbname(NN,'GAGA',E_GG,  'E_GG(nGG):R')
       call hbname(NN,'GAGA',P_GG,  'P_GG(nGG):R')

c      call hbname(NN,'GAGAGA',nGGG,       'nGGG[0,4]:I')
c      call hbname(NN,'GAGAGA',W2_GGG,     'W2_GGG(nGGG):R')

c      call hbname(NN,'GAGAGAGA',nGGGG,   'nGGGG[0,1]:I')
c      call hbname(NN,'GAGAGAGA',W2_GGGG, 'W2_GGGG(nGGGG):R')


c ----- DETected groups of particles -----

      call hbname(NN,'DT_ELNU',neN,   'neN[0,1]:I')
      call hbname(NN,'DT_ELNU',M2_eN, 'M2_eN(neN):R')
c      call hbname(NN,'DT_ELNU',E_eN,  'E_eN(neN):R')
c      call hbname(NN,'DT_ELNU',p_eN,  'p_eN(neN):R')
c      call hbname(NN,'DT_ELNU',th_eN, 'th_eN(neN):R')
c      call hbname(NN,'DT_ELNU',ph_eN, 'ph_eN(neN):R')
c      call hbname(NN,'DT_ELNU',q_eN,  'q_eN(neN):R')

      call hbname(NN,'DT_ep',nep,     'nep[0,3]:I')
      call hbname(NN,'DT_ep',M2_ep,   'M2_ep(nep):R')

      call hbname(NN,'DT_ek',nek,     'nek[0,3]:I')
      call hbname(NN,'DT_ek',M2_ek,   'M2_ek(nek):R')

      call hbname(NN,'DT_eNp',neNp,   'neNp[0,3]:I')
      call hbname(NN,'DT_eNp',M2_eNp, 'M2_eNp(neNp):R')
      call hbname(NN,'DT_eNp',W2_Np,  'W2_Np(neNp):R')

      call hbname(NN,'DT_epp',nepp,   'nepp[0,3]:I')
      call hbname(NN,'DT_epp',M2_epp, 'M2_epp(nepp):R')
      call hbname(NN,'DT_epp',W2_pp,  'W2_pp(nepp):R')

      call hbname(NN,'DT_eNk',neNk,   'neNk[0,3]:I')
      call hbname(NN,'DT_eNk',M2_eNk, 'M2_eNk(neNk):R')
      call hbname(NN,'DT_eNk',W2_Nk,  'W2_Nk(neNk):R')

      call hbname(NN,'DT_epk',nepk,   'nepk[0,3]:I')
      call hbname(NN,'DT_epk',M2_epk, 'M2_epk(nepk):R')
      call hbname(NN,'DT_epk',W2_pk,  'W2_pk(nepk):R')

      call hbname(NN,'DT_ekk',nekk,   'nekk[0,3]:I')
      call hbname(NN,'DT_ekk',M2_ekk, 'M2_ekk(nekk):R')
      call hbname(NN,'DT_ekk',W2_kk,  'W2_kk(nekk):R')

      call hbname(NN,'DT_eNpp',neNpp,   'neNpp[0,3]:I')
      call hbname(NN,'DT_eNpp',M2_eNpp, 'M2_eNpp(neNpp):R')
      call hbname(NN,'DT_eNpp',E_eNpp,  'E_eNpp(neNpp):R')
      call hbname(NN,'DT_eNpp',W2_Npp,  'W2_Npp(neNpp):R')

      call hbname(NN,'DT_eppp',neppp,   'neppp[0,3]:I')
      call hbname(NN,'DT_eppp',M2_eppp, 'M2_eppp(neppp):R')
      call hbname(NN,'DT_eppp',W2_ppp,  'W2_ppp(neppp):R')

      call hbname(NN,'DT_eNkk',neNkk,   'neNkk[0,3]:I')
      call hbname(NN,'DT_eNkk',M2_eNkk, 'M2_eNkk(neNkk):R')
      call hbname(NN,'DT_eNkk',W2_Nkk,  'W2_Nkk(neNkk):R')

      call hbname(NN,'DT_eNpk',neNpk,   'neNpk[0,3]:I')
      call hbname(NN,'DT_eNpk',M2_eNpk, 'M2_eNpk(neNpk):R')
      call hbname(NN,'DT_eNpk',W2_Npk,  'W2_Npk(neNpk):R')

      call hbname(NN,'DT_eNppp',neNppp,   'neNppp[0,3]:I')
      call hbname(NN,'DT_eNppp',M2_eNppp, 'M2_eNppp(neNppp):R')
      call hbname(NN,'DT_eNppp',W2_Nppp,  'W2_Nppp(neNppp):R')

c      call hbname(NN,'DT_pipi',ncmrest,   'ncmrest[0,1]:I')
c      call hbname(NN,'DT_pipi',E_lab,     'E_lab(ncmrest):R')
c      call hbname(NN,'DT_pipi',p_lab,     'p_lab(ncmrest):R')
c      call hbname(NN,'DT_pipi',th_lab,    'th_lab(ncmrest):R')
c      call hbname(NN,'DT_pipi',ph_lab,    'ph_lab(ncmrest):R')
c      call hbname(NN,'DT_pipi',E_cm,      'E_cm(ncmrest):R')
c      call hbname(NN,'DT_pipi',p_cm,      'p_cm(ncmrest):R')
c      call hbname(NN,'DT_pipi',th_cm,     'th_cm(ncmrest):R')
c      call hbname(NN,'DT_pipi',ph_cm,     'ph_cm(ncmrest):R')
c      call hbname(NN,'DT_pipi',E_rest,      'E_rest(ncmrest):R')
c      call hbname(NN,'DT_pipi',p_rest,      'p_rest(ncmrest):R')
c      call hbname(NN,'DT_pipi',th_rest,     'th_rest(ncmrest):R')
c      call hbname(NN,'DT_pipi',ph_rest,     'ph_rest(ncmrest):R')

c      call hbname(NN,'DT_fit',nfit,    'nfit[0,3]:I')
c      call hbname(NN,'DT_fit',fit_Eg,  'fit_Eg(nfit):R')
c      call hbname(NN,'DT_fit',fit_dEg, 'fit_dEg(nfit):R')
c      call hbname(NN,'DT_fit',fit_p,   'fit_p(nfit):R')
c      call hbname(NN,'DT_fit',fit_dp,  'fit_dp(nfit):R')
c      call hbname(NN,'DT_fit',fit_th,  'fit_th(nfit):R')
c      call hbname(NN,'DT_fit',fit_dt,  'fit_dt(nfit):R')
c      call hbname(NN,'DT_fit',fit_fi,  'fit_fi(nfit):R')
c      call hbname(NN,'DT_fit',fit_df,  'fit_df(nfit):R')
c      call hbname(NN,'DT_fit',fit_W,   'fit_W(nfit):R')
c      call hbname(NN,'DT_fit',fit_ch2, 'fit_ch2(nfit):R')
c      call hbname(NN,'DT_fit',fit_err, 'fit_err(nfit):I')


c ----- DETected groups of particles (neutral/charged combinations)

c 1 gamma
      call hbname(NN,'DT_eNG',neNG,   'neNG[0,4]:I')
      call hbname(NN,'DT_eNG',M2_eNG, 'M2_eNG(neNG):R')
      call hbname(NN,'DT_eNG',W2_NG,  'W2_NG(neNG):R')

c 2 Gammas
c      call hbname(NN,'DT_e0',ne0,     'ne0[0,1]:I')
c      call hbname(NN,'DT_e0',M2_e0,   'M2_e0(ne0):R')

c      call hbname(NN,'DT_eN0',neN0,   'neN0[0,1]:I')
c      call hbname(NN,'DT_eN0',M2_eN0, 'M2_eN0(neN0):R')
c      call hbname(NN,'DT_eN0',W2_N0,  'W2_N0(neN0):R')
c      call hbname(NN,'DT_eN0',E_eN0,  'E_eN0(neN0):R')
c      call hbname(NN,'DT_eN0',p_eN0,  'p_eN0(neN0):R')
c      call hbname(NN,'DT_eN0',th_eN0, 'th_eN0(neN0):R')
c      call hbname(NN,'DT_eN0',ph_eN0, 'ph_eN0(neN0):R')
c      call hbname(NN,'DT_eN0',q_eN0,  'q_eN0(neN0):R')

      call hbname(NN,'DT_ep0',nep0,   'nep0[0,3]:I')
      call hbname(NN,'DT_ep0',M2_ep0, 'M2_ep0(nep0):R')
      call hbname(NN,'DT_ep0',W2_p0,  'W2_p0(nep0):R')

c      call hbname(NN,'DT_ek0',nek0,   'nek0[0,3]:I')
c      call hbname(NN,'DT_ek0',M2_ek0, 'M2_ek0(nek0):R')
c      call hbname(NN,'DT_ek0',W2_k0,  'W2_k0(nek0):R')

c      call hbname(NN,'DT_eNp0',neNp0,  'neNp0[0,3]:I')
c      call hbname(NN,'DT_eNp0',M2_eNp0,'M2_eNp0(neNp0):R')
c      call hbname(NN,'DT_eNp0',W2_Np0 ,'W2_Np0(neNp0):R')

c      call hbname(NN,'DT_eNk0',neNk0,  'neNk0[0,3]:I')
c      call hbname(NN,'DT_eNk0',M2_eNk0,'M2_eNk0(neNk0):R')
c      call hbname(NN,'DT_eNk0',W2_Nk0, 'W2_Nk0(neNk0):R')

c      call hbname(NN,'DT_epk0',nepk0,  'nepk0[0,3]:I')
c      call hbname(NN,'DT_epk0',M2_epk0,'M2_epk0(nepk0):R')
c      call hbname(NN,'DT_epk0',W2_pk0, 'W2_pk0(nepk0):R')

c      call hbname(NN,'DT_ekk0',nekk0,  'nekk0[0,3]:I')
c      call hbname(NN,'DT_ekk0',M2_ekk0,'M2_ekk0(nekk0):R')
c      call hbname(NN,'DT_ekk0',W2_kk0, 'W2_kk0(nekk0):R')

      call hbname(NN,'DT_epp0',nepp0,  'nepp0[0,3]:I')
      call hbname(NN,'DT_epp0',M2_epp0,'M2_epp0(nepp0):R')
      call hbname(NN,'DT_epp0',W2_pp0, 'W2_pp0(nepp0):R')

c      call hbname(NN,'DT_eppp0',neppp0,  'neppp0[0,3]:I')
c      call hbname(NN,'DT_eppp0',M2_eppp0,'M2_eppp0(neppp0):R')
c      call hbname(NN,'DT_eppp0',W2_ppp0, 'W2_ppp0(neppp0):R')

c      call hbname(NN,'DT_eNkk0',neNkk0,  'neNkk0[0,3]:I')
c      call hbname(NN,'DT_eNkk0',M2_eNkk0,'M2_eNkk0(nekk0):R')
c      call hbname(NN,'DT_eNkk0',W2_Nkk0, 'W2_Nkk0(nekk0):R')

c      call hbname(NN,'DT_eNpk0',neNpk0,  'neNpk0[0,3]:I')
c      call hbname(NN,'DT_eNpkk0',M2_eNpk0,'M2_eNpk0(nekk0):R')
c      call hbname(NN,'DT_eNpk0',W2_Npk0, 'W2_Npk0(nekk0):R')

c      call hbname(NN,'DT_eNppp0',neNppp0,  'neNppp0[0,3]:I')
c      call hbname(NN,'DT_eNppp0',M2_eNppp0,'M2_eNppp0(neppp0):R')
c      call hbname(NN,'DT_eNppp0',W2_Nppp0, 'W2_Nppp0(neppp0):R')


      ENDIF !IF(ana_gen.eq.1)
      ENDDO !DO NN=10,11




c---------------------------------------------------------------------
c (P pi0)  Specific channel ntuples #25
c---------------------------------------------------------------------

      IF (ana_chan.eq.1) THEN
      
      call hbnt  (25,'1PION',' ')

      call hbname(25,'ana_1pia',indtype,        'indtype[0,5]:I')
      call hbname(25,'ana_1pia',n_ev_run_a, 'n_ev_run(indtype):I')
      call hbname(25,'ana_1pia',trig_type_a,'trig_type(indtype):I')
      call hbname(25,'ana_1pia',trig_clas_a,'trig_clas(indtype):I')
      call hbname(25,'ana_1pia',n_ev_nt_a,  'n_ev_nt(indtype):I')
      call hbname(25,'ana_1pia',n_part_a,   'n_part(indtype):I')
      call hbname(25,'ana_1pia',ntot_a,     'ntot(indtype):I')
      call hbname(25,'ana_1pia',bit_a,      'bit(indtype):I')
      call hbname(25,'ana_1pia',bit1pi_a,   'bit1pi(indtype):I')
      call hbname(25,'ana_1pia',bit2pi_a,   'bit2pi(indtype):I')
      call hbname(25,'ana_1pia',sec_coin_SC_a,      'sec_SC:I')
      call hbname(25,'ana_1pia',sec_coin_ST_a,      'sec_ST:I')

      IF(beam_type) THEN
      call hbname(25,'ana_1pia',Qgated_total,  'Qgated_total:R')
      call hbname(25,'ana_1pia',Qgated_plus,   'Qgated_plus:R')
      call hbname(25,'ana_1pia',omega_a,       'omega(indtype):R')
      call hbname(25,'ana_1pia',Q2_a,          'Q2(indtype):R')
      call hbname(25,'ana_1pia',W_el_a,        'W(indtype):R')
      call hbname(25,'ana_1pia',eps_a,         'eps(indtype):R')
      call hbname(25,'ana_1pia',eps_l_a,       'eps_l(indtype):R')
      call hbname(25,'ana_1pia',Flux_vp_a,     'Flux_vp(indtype):R')
      call hbname(25,'ana_1pia',E_EL_a,        'E_EL(indtype):R')
      call hbname(25,'ana_1pia',P_EL_a,        'P_EL(indtype):R')
      call hbname(25,'ana_1pia',B_EL_a,        'B_EL(indtype):R')
      call hbname(25,'ana_1pia',th_EL_a,       'th_EL(indtype):R')
      call hbname(25,'ana_1pia',ph_EL_a,       'ph_EL(indtype):R')
      call hbname(25,'ana_1pia',acc_EL_a,      'acc_EL(indtype):R')
      call hbname(25,'ana_1pia',x_EL_a,        'x_EL(indtype):R')
      call hbname(25,'ana_1pia',y_EL_a,        'y_EL(indtype):R')
      call hbname(25,'ana_1pia',z_EL_a,        'z_EL(indtype):R')
      call hbname(25,'ana_1pia',ECtot_EL_a,    'ECtot_EL(indtype):R')
      call hbname(25,'ana_1pia',ECin_EL_a,     'ECin_EL(indtype):R')
      call hbname(25,'ana_1pia',ECout_EL_a,    'ECout_EL(indtype):R')
      
      ELSE
      call hbname(25,'ana_1pia',Egamma_a,      'Egamma(indtype):R')
      call hbname(25,'ana_1pia',W_ph_a,           'W(indtype):R')
      call hbname(25,'ana_1pia',Tag_time_a,    'Tag_time(indtype):R')
      call hbname(25,'ana_1pia',ST_time_a,     'ST_time(indtype):R')
      call hbname(25,'ana_1pia',RF1_time_a,    'RF1_time(indtype):R')
      call hbname(25,'ana_1pia',RF2_time_a,    'RF2_time(indtype):R')
      call hbname(25,'ana_1pia',hit_in_time_a, 'hit_in_time(indtype):I')
      ENDIF

c Proton  in lab
      call hbname(25,'ana_1pia',det_P,      'det_P(indtype):I')
      call hbname(25,'ana_1pia',E_P_lab,    'E_P_lab(indtype):R')
      call hbname(25,'ana_1pia',P_P_lab,    'P_P_lab(indtype):R')
      call hbname(25,'ana_1pia',B_P_lab,    'B_P_lab(indtype):R')
      call hbname(25,'ana_1pia',M2_P_lab,   'M2_P_lab(indtype):R')
      call hbname(25,'ana_1pia',th_P_lab,   'th_P_lab(indtype):R')
      call hbname(25,'ana_1pia',ph_P_lab,   'ph_P_lab(indtype):R')
      call hbname(25,'ana_1pia',acc_P_lab,  'acc_P_lab(indtype):R')
      call hbname(25,'ana_1pia',x_P_lab,    'x_P_lab(indtype):R')
      call hbname(25,'ana_1pia',y_P_lab,    'y_P_lab(indtype):R')
      call hbname(25,'ana_1pia',z_P_lab,    'z_P_lab(indtype):R')

      if (beam_type) then
      call hbname(25,'ana_1pia',th1_P_lab,  'th1_P_lab(indtype):R')
      call hbname(25,'ana_1pia',th2_P_lab,  'th2_P_lab(indtype):R')
      endif
      call hbname(25,'ana_1pia',M2_eN_a,  'M2_eN(indtype):R')
      call hbname(25,'ana_1pia',E_eN_a,   'E_eN(indtype):R')
      call hbname(25,'ana_1pia',P_eN_a,   'P_eN(indtype):R')
      call hbname(25,'ana_1pia',th_eN_a,  'th_eN(indtype):R')
      call hbname(25,'ana_1pia',ph_eN_a,  'ph_eN(indtype):R')

c pi0 in lab
      call hbname(25,'ana_1pia',det_PI0,      'det_PI0(indtype):I')
      call hbname(25,'ana_1pia',E_PI0_lab,    'E_PI0_lab(indtype):R')
      call hbname(25,'ana_1pia',P_PI0_lab,    'P_PI0_lab(indtype):R')
      call hbname(25,'ana_1pia',B_PI0_lab,    'B_PI0_lab(indtype):R')
      call hbname(25,'ana_1pia',M2_PI0_lab,   'M2_PI0_lab(indtype):R')
      call hbname(25,'ana_1pia',th_PI0_lab,   'th_PI0_lab(indtype):R')
      call hbname(25,'ana_1pia',ph_PI0_lab,   'ph_PI0_lab(indtype):R')

c gamma 1 in lab
      call hbname(25,'ana_1pia',det_G1,       'det_G1(indtype):I')
      call hbname(25,'ana_1pia',E_G1_lab,     'E_G1_lab(indtype):R')
      call hbname(25,'ana_1pia',th_G1_lab,    'th_G1_lab(indtype):R')
      call hbname(25,'ana_1pia',ph_G1_lab,    'ph_G1_lab(indtype):R')
c      call hbname(25,'ana_1pia',acc_G1_lab,   'acc_G1_lab(indtype):R')
c      call hbname(25,'ana_1pia',x_G1_lab,     'x_G1_lab(indtype):R')
c      call hbname(25,'ana_1pia',y_G1_lab,     'y_G1_lab(indtype):R')
c      call hbname(25,'ana_1pia',z_G1_lab,     'z_G1_lab(indtype):R')

c gamma 2 in lab
      call hbname(25,'ana_1pia',det_G2,       'det_G2(indtype):I')
      call hbname(25,'ana_1pia',E_G2_lab,     'E_G2_lab(indtype):R')
      call hbname(25,'ana_1pia',th_G2_lab,    'th_G2_lab(indtype):R')
      call hbname(25,'ana_1pia',ph_G2_lab,    'ph_G2_lab(indtype):R')
c      call hbname(25,'ana_1pia',acc_G2_lab,   'acc_G2_lab(indtype):R')
c      call hbname(25,'ana_1pia',x_G2_lab,     'x_G2_lab(indtype):R')
c      call hbname(25,'ana_1pia',y_G2_lab,     'y_G2_lab(indtype):R')
c      call hbname(25,'ana_1pia',z_G2_lab,     'z_G2_lab(indtype):R')

c Missing and invariant eNp combinations
      call hbname(25,'ana_1pia',M2_eNp_a,   'M2_eNp(indtype):R')
      call hbname(25,'ana_1pia',W2_Np_a,    'W2_Np(indtype):R')

c Delta+ (Ppi0) in lab
      call hbname(25,'ana_1pia',M2_Dp,      'M2_Dp(indtype):R')
      call hbname(25,'ana_1pia',E_Dp_lab,   'E_Dp_lab(indtype):R')
      call hbname(25,'ana_1pia',P_Dp_lab,   'P_Dp_lab(indtype):R')
      call hbname(25,'ana_1pia',th_Dp_lab,  'th_Dp_lab(indtype):R')
      call hbname(25,'ana_1pia',ph_Dp_lab,  'ph_Dp_lab(indtype):R')

c pi0 in hcm  frame  
      call hbname(25,'ana_1pia',E_PI0_hcm,   'E_PI0_hcm(indtype):R')
      call hbname(25,'ana_1pia',P_PI0_hcm,   'P_PI0_hcm(indtype):R')
      call hbname(25,'ana_1pia',th_PI0_hcm,  'th_PI0_hcm(indtype):R')
      call hbname(25,'ana_1pia',ph_PI0_hcm,  'ph_PI0_hcm(indtype):R')



c---------------------------------------------------------------------
c (n pi+)  Specific channel ntuples #26
c---------------------------------------------------------------------

      ELSEIF(ana_chan.eq.2) THEN

      call hbnt  (26,'1PION',' ')
         
      call hbname(26,'ana_1pib',indtype,        'indtype[0,5]:I')
      call hbname(26,'ana_1pib',n_ev_run_a, 'n_ev_run(indtype):I')
      call hbname(26,'ana_1pib',trig_type_a,'trig_type(indtype):I')
      call hbname(26,'ana_1pib',trig_clas_a,'trig_clas(indtype):I')
      call hbname(26,'ana_1pib',n_ev_nt_a,  'n_ev_nt(indtype):I')
      call hbname(26,'ana_1pib',n_part_a,   'n_part(indtype):I')
      call hbname(26,'ana_1pib',ntot_a,     'ntot(indtype):I')
      call hbname(26,'ana_1pib',bit_a,      'bit(indtype):I')
      call hbname(26,'ana_1pib',bit1pi_a,   'bit1pi(indtype):I')
      call hbname(26,'ana_1pib',bit2pi_a,   'bit2pi(indtype):I')
      call hbname(26,'ana_1pib',sec_coin_SC_a,      'sec_SC:I')
      call hbname(26,'ana_1pib',sec_coin_ST_a,      'sec_ST:I')

      IF(beam_type) THEN
      call hbname(26,'ana_1pib',omega_a,       'omega(indtype):R')
      call hbname(26,'ana_1pib',Q2_a,          'Q2(indtype):R')
      call hbname(26,'ana_1pib',W_el_a,        'W(indtype):R')
      call hbname(26,'ana_1pib',eps_a,         'eps(indtype):R')
      call hbname(26,'ana_1pib',eps_l_a,       'eps_l(indtype):R')
      call hbname(26,'ana_1pib',Flux_vp_a,     'Flux_vp(indtype):R')
      call hbname(26,'ana_1pib',E_EL_a,        'E_EL(indtype):R')
      call hbname(26,'ana_1pib',P_EL_a,        'P_EL(indtype):R')
      call hbname(26,'ana_1pib',B_EL_a,        'B_EL(indtype):R')
      call hbname(26,'ana_1pib',th_EL_a,       'th_EL(indtype):R')
      call hbname(26,'ana_1pib',ph_EL_a,       'ph_EL(indtype):R')
      call hbname(26,'ana_1pib',acc_EL_a,      'acc_EL(indtype):R')
      call hbname(26,'ana_1pib',x_EL_a,        'x_EL(indtype):R')
      call hbname(26,'ana_1pib',y_EL_a,        'y_EL(indtype):R')
      call hbname(26,'ana_1pib',z_EL_a,        'z_EL(indtype):R')
      call hbname(26,'ana_1pib',ECtot_EL_a,    'ECtot_EL(indtype):R')
      call hbname(26,'ana_1pib',ECin_EL_a,     'ECin_EL(indtype):R')
      call hbname(26,'ana_1pib',ECout_EL_a,    'ECout_EL(indtype):R')
       ELSE
      call hbname(26,'ana_1pib',Egamma_a,      'Egamma(indtype):R')
      call hbname(26,'ana_1pib',W_ph_a,        'W(indtype):R')
      call hbname(26,'ana_1pib',Tag_time_a,    'Tag_time(indtype):R')
      call hbname(26,'ana_1pib',ST_time_a,     'ST_time(indtype):R')
      call hbname(26,'ana_1pib',RF1_time_a,    'RF1_time(indtype):R')
      call hbname(26,'ana_1pib',RF2_time_a,    'RF2_time(indtype):R')
      call hbname(26,'ana_1pib',hit_in_time_a, 'hit_in_time(indtype):I')
      ENDIF

c Neutron in lab
      call hbname(26,'ana_1pib',det_n,      'det_n(indtype):I')
      call hbname(26,'ana_1pib',E_n_lab,    'E_n_lab(indtype):R')
      call hbname(26,'ana_1pib',P_n_lab,    'P_n_lab(indtype):R')
      call hbname(26,'ana_1pib',B_n_lab,    'B_n_lab(indtype):R')
      call hbname(26,'ana_1pib',M2_n_lab,   'M2_n_lab(indtype):R')
      call hbname(26,'ana_1pib',th_n_lab,   'th_n_lab(indtype):R')
      call hbname(26,'ana_1pib',ph_n_lab,   'ph_n_lab(indtype):R')
      call hbname(26,'ana_1pib',acc_n_lab,  'acc_n_lab(indtype):R')
      call hbname(26,'ana_1pib',x_n_lab,    'x_n_lab(indtype):R')
      call hbname(26,'ana_1pib',y_n_lab,    'y_n_lab(indtype):R')
      call hbname(26,'ana_1pib',z_n_lab,    'z_n_lab(indtype):R')

c pi+ in lab
      call hbname(26,'ana_1pib',det_PIp,      'det_PIp(indtype):I')
      call hbname(26,'ana_1pib',E_PIp_lab,    'E_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',P_PIp_lab,    'P_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',B_PIp_lab,    'B_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',M2_PIp_lab,   'M2_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',th_PIp_lab,   'th_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',ph_PIp_lab,   'ph_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',acc_PIp_lab,  'acc_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',x_PIp_lab,    'x_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',y_PIp_lab,    'y_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',z_PIp_lab,    'z_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',ECtot_PIp_lab,'ECtot_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',ECin_PIp_lab, 'ECin_PIp_lab(indtype):R')
      call hbname(26,'ana_1pib',ECout_PIp_lab,'ECout_PIp_lab(indtype):R')

c Missing and invariant eNp combinations
      call hbname(26,'ana_1pib',M2_eNp_a,  'M2_eNp(indtype):R')
      call hbname(26,'ana_1pib',W2_Np_a,   'W2_Np(indtype):R')

c Delta+ (Npi+) in lab
      call hbname(26,'ana_1pib',M2_Dp,     'M2_Dp(indtype):R')
      call hbname(26,'ana_1pib',E_Dp_lab,  'E_Dp_lab(indtype):R')
      call hbname(26,'ana_1pib',P_Dp_lab,  'P_Dp_lab(indtype):R')
      call hbname(26,'ana_1pib',th_Dp_lab, 'th_Dp_lab(indtype):R')
      call hbname(26,'ana_1pib',ph_Dp_lab, 'ph_Dp_lab(indtype):R')

c pi+ in hcm rest frame
      call hbname(26,'ana_1pib',E_PIp_hcm, 'E_PIp_hcm(indtype):R')
      call hbname(26,'ana_1pib',P_PIp_hcm, 'P_PIp_hcm(indtype):R')
      call hbname(26,'ana_1pib',th_PIp_hcm,'th_PIp_hcm(indtype):R')
      call hbname(26,'ana_1pib',ph_PIp_hcm,'ph_PIp_hcm(indtype):R')




c---------------------------------------------------------------------
c (pi+ pi-)  Specific channel ntuples #20
c---------------------------------------------------------------------

      ELSEIF (ana_chan.eq.3) THEN
      call hbnt  (20,'2PION',' ')
      call hbname(20,'ana_2pi',elast,  'elast:R')
*      call hbname(20,'ana_2pi',fcg_hevt, 'fcg_hevt:R')
*      call hbname(20,'ana_2pi',Qgated, 'Qgated:R')
*      call hbname(20,'ana_2pi',Qgated_plus, 'Qgated_plus:R') 
      call hbname(20,'ana_2pi',Iblock_all,  'Iblock_all:I')
*      call hbname(20,'ana_2pi',n_e_incl,  'n_e_incl:I')
      call hbname(20,'ana_2pi',tg_hevt,     'tg_hevt:R')
      call hbname(20,'ana_2pi',event_clas, 'event_clas:I')
      call hbname(20,'ana_2pi',trkEVNT_PI, 'trkEVNT_PI:I')
      call hbname(20,'ana_2pi',trkEVNT_NU, 'trkEVNT_NU:I')
      call hbname(20,'ana_2pi', trkEVNT_EL, 'trkEVNT_EL:I')
      call hbname(20,'ana_2pi',indtype,        'indtype[0,5]:I')
      call hbname(20,'ana_2pi',n_ev_run_a, 'n_ev_run(indtype):I')
      call hbname(20,'ana_2pi',trig_type_a,'trig_type(indtype):I')
      call hbname(20,'ana_2pi',trig_clas_a,'trig_clas(indtype):I')
      call hbname(20,'ana_2pi',n_ev_nt_a,  'n_ev_nt(indtype):I')
      call hbname(20,'ana_2pi',n_part_a,   'n_part(indtype):I')
      call hbname(20,'ana_2pi',ntot_a,     'ntot(indtype):I')
      call hbname(20,'ana_2pi',bit_a,      'bit(indtype):I')
      call hbname(20,'ana_2pi',bit1pi_a,   'bit1pi(indtype):I')
      call hbname(20,'ana_2pi',bit2pi_a,   'bit2pi(indtype):I')
      call hbname(20,'ana_2pi',sec_coin_SC_a,    'sec_SC(indtype):I')
      call hbname(20,'ana_2pi',sec_coin_ST_a,    'sec_ST(indtype):I')
      call hbname(20,'ana_2pi',sec_coin_SC_thr_a,'sec_SC_thr(indtype):I')

      IF(beam_type) THEN
      call hbname(20,'ana_2pi',omega_a,       'omega(indtype):R')
      call hbname(20,'ana_2pi',Q2_a,          'Q2(indtype):R')
      call hbname(20,'ana_2pi',W_EL_a,        'W(indtype):R')
      call hbname(20,'ana_2pi',eps_a,         'eps(indtype):R')
      call hbname(20,'ana_2pi',eps_l_a,       'eps_l(indtype):R')
      call hbname(20,'ana_2pi',Flux_vp_a,     'Flux_vp(indtype):R')
      call hbname(20,'ana_2pi',E_EL_a,        'E_EL(indtype):R')
      call hbname(20,'ana_2pi',P_EL_a,        'P_EL(indtype):R')
      call hbname(20,'ana_2pi',B_EL_a,        'B_EL(indtype):R')
      call hbname(20,'ana_2pi',th_EL_a,       'th_EL(indtype):R')
      call hbname(20,'ana_2pi',ph_EL_a,       'ph_EL(indtype):R')
      call hbname(20,'ana_2pi',acc_EL_a,      'acc_EL(indtype):R')
      call hbname(20,'ana_2pi',x_EL_a,        'x_EL(indtype):R')
      call hbname(20,'ana_2pi',y_EL_a,        'y_EL(indtype):R')
      call hbname(20,'ana_2pi',z_EL_a,        'z_EL(indtype):R')
      call hbname(20,'ana_2pi',ECtot_EL_a,    'ECtot_EL(indtype):R')
      call hbname(20,'ana_2pi',ECin_EL_a,     'ECin_EL(indtype):R')
      call hbname(20,'ana_2pi',ECout_EL_a,    'ECout_EL(indtype):R')
      call hbname(20,'ana_2pi',NpheCC_EL_a,   'NpheCC_EL(indtype):R')
      call hbname(20,'ana_2pi',EffCC_EL_a,    'EffCC_EL(indtype):R')
      call hbname(20,'ana_2pi',PdHit_EL_a,    'PdHit_EL(indtype):I')
      call hbname(20,'ana_2pi',jopa,          'joap(indtype):R')
      
      ELSE
      call hbname(20,'ana_2pi',Egamma_a,      'Egamma(indtype):R')
      call hbname(20,'ana_2pi',W_ph_a,        'W(indtype):R')
      call hbname(20,'ana_2pi',Tag_time_a,    'Tag_time(indtype):R')
      call hbname(20,'ana_2pi',ST_time_a,     'ST_time(indtype):R')
      call hbname(20,'ana_2pi',RF1_time_a,    'RF1_time(indtype):R')
      call hbname(20,'ana_2pi',RF2_time_a,    'RF2_time(indtype):R')
      call hbname(20,'ana_2pi',hit_in_time_a, 'hit_in_time(indtype):I')
      call hbname(20,'ana_2pi',chan_a,        'chanMC(indtype):I')
      call hbname(20,'ana_2pi',true_a,        'rand(indtype):I')
      call hbname(20,'ana_2pi',trig_SC,       'trig_SC(indtype):R')
      call hbname(20,'ana_2pi',T_id_a,       'T_Id(indtype):I')
      call hbname(20,'ana_2pi',E_id_a,       'E_Id(indtype):I')
      call hbname(20,'ana_2pi',Eg_fit,     'Eg_fit(indtype):R')
      call hbname(20,'ana_2pi',DEg_fit,    'DEg_fit(indtype):R')
      call hbname(20,'ana_2pi',DP_P_fit,   'DP_P_fit(indtype):R')
      call hbname(20,'ana_2pi',DP_PIp_fit, 'DP_PIp_fit(indtype):R')
      call hbname(20,'ana_2pi',DP_PIm_fit, 'DP_PIm_fit(indtype):R')
      call hbname(20,'ana_2pi',DTH_P_fit,  'DTH_P_fit(indtype):R')
      call hbname(20,'ana_2pi',DTH_PIp_fit,'DTH_PIp_fit(indtype):R')
      call hbname(20,'ana_2pi',DTH_PIm_fit,'DTH_PIm_fit(indtype):R')
      call hbname(20,'ana_2pi',DPH_P_fit,  'DPH_P_fit(indtype):R')
      call hbname(20,'ana_2pi',DPH_PIp_fit,'DPH_PIp_fit(indtype):R')
      call hbname(20,'ana_2pi',DPH_PIm_fit,'DPH_PIm_fit(indtype):R')
      call hbname(20,'ana_2pi',CHI_fit,    'CHI_fit(indtype):R')
      call hbname(20,'ana_2pi',Err_fit,    'Err_fit(indtype):I')
      ENDIF

c Missing masses e+part
      if (beam_type) call hbname(20,'ana_2pi',M2_eN_a2,      'M2_eN(indtype):R')
      if (beam_type) call hbname(20,'ana_2pi',M2_epip_a2,    'M2_epip(indtype):R')
      if (beam_type) call hbname(20,'ana_2pi',M2_epim_a2,    'M2_epim(indtype):R')
      if (beam_type) call hbname(20,'ana_2pi',M2_eNpip_a2,   'M2_eNpip(indtype):R')
      if (beam_type) call hbname(20,'ana_2pi',M2_eNpim_a2,   'M2_eNpim(indtype):R')

c Proton in lab
      call hbname(20,'ana_2pi',det_P,      'det_P(indtype):I')
      call hbname(20,'ana_2pi',E_P_lab,    'E_P_lab(indtype):R')
      call hbname(20,'ana_2pi',P_P_lab,    'P_P_lab(indtype):R')
      call hbname(20,'ana_2pi',B_P_lab,    'B_P_lab(indtype):R')
      call hbname(20,'ana_2pi',M2_P_lab,   'M2_P_lab(indtype):R')
      call hbname(20,'ana_2pi',th_P_lab,   'th_P_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_P_lab,   'ph_P_lab(indtype):R')
      call hbname(20,'ana_2pi',acc_P_lab,  'acc_P_lab(indtype):R')
      call hbname(20,'ana_2pi',x_P_lab,    'x_P_lab(indtype):R')
      call hbname(20,'ana_2pi',y_P_lab,    'y_P_lab(indtype):R')
      call hbname(20,'ana_2pi',z_P_lab,    'z_P_lab(indtype):R')
      if (.not.beam_type) call hbname(20,'ana_2pi',id_SC_p,'id_SC_p(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',sec_SC_p,'sec_SC_p(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',thrSC_p,'thrSC_p(indtype):R')
      if (beam_type) then
      call hbname(20,'ana_2pi',th1_P_lab,  'th1_P_lab(indtype):R')
      call hbname(20,'ana_2pi',th2_P_lab,  'th2_P_lab(indtype):R')
      endif
      if (beam_type)   call hbname(20,'ana_2pi',PdHit_P_lab,'PdHit_P_lab(indtype):I')
      
c pi+ in lab
      call hbname(20,'ana_2pi',det_PIp,      'det_PIp(indtype):I')
      call hbname(20,'ana_2pi',E_PIp_lab,    'E_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',P_PIp_lab,    'P_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',B_PIp_lab,    'B_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',M2_PIp_lab,   'M2_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',th_PIp_lab,   'th_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_PIp_lab,   'ph_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',acc_PIp_lab,  'acc_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',x_PIp_lab,    'x_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',y_PIp_lab,    'y_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',z_PIp_lab,    'z_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',ECtot_PIp_lab,'ECtot_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',ECin_PIp_lab, 'ECin_PIp_lab(indtype):R')
      call hbname(20,'ana_2pi',ECout_PIp_lab,'ECout_PIp_lab(indtype):R')
      if (beam_type) call hbname(20,'ana_2pi',PdHit_PIp_lab,'PdHit_PIp_lab(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',id_SC_pip,'id_SC_pip(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',sec_SC_pip,'sec_SC_pip(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',thrSC_pip,'thrSC_pip(indtype):R')

c pi- in lab
      call hbname(20,'ana_2pi',det_PIm,      'det_PIm(indtype):I')
      call hbname(20,'ana_2pi',E_PIm_lab,    'E_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',P_PIm_lab,    'P_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',B_PIm_lab,    'B_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',M2_PIm_lab,   'M2_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',th_PIm_lab,   'th_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_PIm_lab,   'ph_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',acc_PIm_lab,  'acc_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',x_PIm_lab,    'x_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',y_PIm_lab,    'y_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',z_PIm_lab,    'z_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',ECtot_PIm_lab,'ECtot_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',ECin_PIm_lab, 'ECin_PIm_lab(indtype):R')
      call hbname(20,'ana_2pi',ECout_PIm_lab,'ECout_PIm_lab(indtype):R')
      if (beam_type)   call hbname(20,'ana_2pi',PdHit_PIm_lab,'PdHit_PIm_lab(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',id_SC_pim,'id_SC_pim(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',sec_SC_pim,'sec_SC_pim(indtype):I')
      if (.not.beam_type) call hbname(20,'ana_2pi',thrSC_pim,'thrSC_pim(indtype):R')

c Missing and invariant eNpp combinations
      call hbname(20,'ana_2pi',M2_eNpp_a,  'M2_eNpp(indtype):R')
      call hbname(20,'ana_2pi',E_eNpp_a,   'E_eNpp(indtype):R')
      call hbname(20,'ana_2pi',W2_Npp_a,   'W2_Npp(indtype):R')

c Rho (pi+pi-) in lab and in hcm
      call hbname(20,'ana_2pi',M2_rho,     'M2_rho(indtype):R')
      call hbname(20,'ana_2pi',E_rho_lab,  'E_rho_lab(indtype):R')
      call hbname(20,'ana_2pi',P_rho_lab,  'P_rho_lab(indtype):R')
      call hbname(20,'ana_2pi',th_rho_lab, 'th_rho_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_rho_lab, 'ph_rho_lab(indtype):R')
      call hbname(20,'ana_2pi',E_rho_hcm,  'E_rho_hcm(indtype):R')
      call hbname(20,'ana_2pi',P_rho_hcm,  'P_rho_hcm(indtype):R')
      call hbname(20,'ana_2pi',th_rho_hcm, 'th_rho_hcm(indtype):R')
      call hbname(20,'ana_2pi',ph_rho_hcm, 'ph_rho_hcm(indtype):R')
c Delta++ (Ppi+) in lab and in hcm
      call hbname(20,'ana_2pi',M2_Dpp,    'M2_Dpp(indtype):R')
      call hbname(20,'ana_2pi',E_Dpp_lab, 'E_Dpp_lab(indtype):R')
      call hbname(20,'ana_2pi',P_Dpp_lab, 'P_Dpp_lab(indtype):R')
      call hbname(20,'ana_2pi',th_Dpp_lab,'th_Dpp_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_Dpp_lab,'ph_Dpp_lab(indtype):R')
      call hbname(20,'ana_2pi',E_Dpp_hcm, 'E_Dpp_hcm(indtype):R')
      call hbname(20,'ana_2pi',P_Dpp_hcm, 'P_Dpp_hcm(indtype):R')
      call hbname(20,'ana_2pi',th_Dpp_hcm,'th_Dpp_hcm(indtype):R')
      call hbname(20,'ana_2pi',ph_Dpp_hcm,'ph_Dpp_hcm(indtype):R')
c Delta0 (Ppi-) in lab and in hcm
      call hbname(20,'ana_2pi',M2_D0,     'M2_D0(indtype):R')
      call hbname(20,'ana_2pi',E_D0_lab,  'E_D0_lab(indtype):R')
      call hbname(20,'ana_2pi',P_D0_lab,  'P_D0_lab(indtype):R')
      call hbname(20,'ana_2pi',th_D0_lab, 'th_D0_lab(indtype):R')
      call hbname(20,'ana_2pi',ph_D0_lab, 'ph_D0_lab(indtype):R')
      call hbname(20,'ana_2pi',E_D0_hcm,  'E_D0_hcm(indtype):R')
      call hbname(20,'ana_2pi',P_D0_hcm,  'P_D0_hcm(indtype):R')
      call hbname(20,'ana_2pi',th_D0_hcm, 'th_D0_hcm(indtype):R')
      call hbname(20,'ana_2pi',ph_D0_hcm, 'ph_D0_hcm(indtype):R')
c pi+ in rho rest frame
      call hbname(20,'ana_2pi',E_PIp_rho, 'E_PIp_rho(indtype):R')
      call hbname(20,'ana_2pi',P_PIp_rho, 'P_PIp_rho(indtype):R')
      call hbname(20,'ana_2pi',th_PIp_rho,'th_PIp_rho(indtype):R')
      call hbname(20,'ana_2pi',ph_PIp_rho,'ph_PIp_rho(indtype):R')
      call hbname(20,'ana_2pi',psi_rho,   'psi_rho(indtype):R')
c pi+ in Delta++ rest frame
      call hbname(20,'ana_2pi',E_PIp_Dpp, 'E_PIp_Dpp(indtype):R')
      call hbname(20,'ana_2pi',P_PIp_Dpp, 'P_PIp_Dpp(indtype):R')
      call hbname(20,'ana_2pi',th_PIp_Dpp,'th_PIp_Dpp(indtype):R')
      call hbname(20,'ana_2pi',ph_PIp_Dpp,'ph_PIp_Dpp(indtype):R')
      call hbname(20,'ana_2pi',psi_Dpp,   'psi_Dpp(indtype):R')
c pi+ in Delta0 rest frame
      call hbname(20,'ana_2pi',E_PIm_D0,  'E_PIm_D0(indtype):R')
      call hbname(20,'ana_2pi',P_PIm_D0,  'P_PIm_D0(indtype):R')
      call hbname(20,'ana_2pi',th_PIm_D0, 'th_PIm_D0(indtype):R')
      call hbname(20,'ana_2pi',ph_PIm_D0, 'ph_PIm_D0(indtype):R')
      call hbname(20,'ana_2pi',psi_D0,    'psi_D0(indtype):R')

c---------------------------------------------------------------------
c Inclusive pion  Specific channel ntuple #31
c---------------------------------------------------------------------

      ELSEIF (ana_chan.eq.4) THEN

      call hbnt  (31,'Inc_pi',' ')
      
c ----- general -----
      call hbname(31,'GENERAL',n_ev_run,   'n_ev_run:I')
      call hbname(31,'GENERAL',trig_type,  'trig_type[-10,10]:I')
      call hbname(31,'GENERAL',trig_clas,  'trig_clas[0,15]:I')
      call hbname(31,'GENERAL',beam_hel,   'beam_hel:I')
      call hbname(31,'GENERAL',beam_helkey,'beam_helkey:I')
      call hbname(31,'GENERAL',beam_keypair,'beam_keypair:I')
      call hbname(31,'GENERAL',beam_keyboth,'beam_keyboth:I')
      call hbname(31,'GENERAL',doplus_BPM, 'doplus_BPM:I')
      call hbname(31,'GENERAL',dominus_BPM, 'dminus_BPM:I')
      call hbname(31,'GENERAL',doplus_FC,  'doplus_FC:I')
      call hbname(31,'GENERAL',dominus_FC, 'dominus_FC:I')
      call hbname(31,'GENERAL',n_ev_nt,    'n_ev_nt:I')
      call hbname(31,'GENERAL',n_part,     'n_part:I')
      call hbname(31,'GENERAL',ntot,       'ntot:I')
      call hbname(31,'GENERAL',bit,        'bit:I')
        IF(beam_type) THEN
      call hbname(31,'GENERAL',Qgated_total, 'Qgated_total:R')
      call hbname(31,'GENERAL',Qgated_plus,  'Qgated_plus:R')
        ENDIF

        IF(beam_type) THEN
      call hbname(31,'ELECTRON', Omega,   'omega:R')      
      call hbname(31,'ELECTRON', Q2,      'Q2:R')
      call hbname(31,'ELECTRON', W,       'W:R')
      call hbname(31,'ELECTRON', eps,     'eps:R')
      call hbname(31,'ELECTRON', eps_l,   'eps_l:R')
      call hbname(31,'ELECTRON', Flux_vp, 'Flux_vp:R')
      call hbname(31,'ana_inc_pi',cosTH_v,   'cosTH_v:R')
      call hbname(31,'ELECTRON', E_EL,    'E_EL:R')
      call hbname(31,'ELECTRON', P_EL,    'P_EL:R')
      call hbname(31,'ELECTRON', B_EL,    'B_EL:R')
      call hbname(31,'ELECTRON', th_EL,   'th_EL:R')      
      call hbname(31,'ELECTRON', ph_EL,   'ph_EL:R')
      call hbname(31,'ELECTRON', acc_EL,  'acc_EL:R')
      call hbname(31,'ELECTRON', X_EL,    'X_EL:R')
      call hbname(31,'ELECTRON', Y_EL,    'Y_EL:R')
      call hbname(31,'ELECTRON', Z_EL,    'Z_EL:R')
      call hbname(31,'ELECTRON', ECtot_EL,'ECtot_EL:R')
      call hbname(31,'ELECTRON', ECin_EL, 'ECin_EL:R')
      call hbname(31,'ELECTRON', ECout_EL,'ECout_EL:R')
        ELSE
      call hbname(31,'PH_BEAM',Egamma,     'Egamma:R') 
      call hbname(31,'PH_BEAM',W,          'W:R')
      call hbname(31,'PH_BEAM',Tag_time,   'Tag_time:R')
      call hbname(31,'PH_BEAM',ST_time,    'ST_time:R')
      call hbname(31,'PH_BEAM',RF1_time,   'RF1_time:R')
      call hbname(31,'PH_BEAM',RF2_time,   'RF2_time:R')
      call hbname(31,'PH_BEAM',hit_in_time,'hit_in_time:I')
      call hbname(31,'PH_BEAM',chan,       'chanMC:I')
      call hbname(31,'PH_BEAM',true,       'rand:I')
        ENDIF

      call hbname(31,'NUCLEON',nNU,    'nNU[0,1]:I')
      call hbname(31,'NUCLEON',E_NU,   'E_NU(nNU):R')
      call hbname(31,'NUCLEON',P_NU,   'P_NU(nNU):R')
      call hbname(31,'NUCLEON',B_NU,   'B_NU(nNU):R')
      call hbname(31,'NUCLEON',M2_NU,  'M2_NU(nNU):R')
      call hbname(31,'NUCLEON',th_NU,  'th_NU(nNU):R')
      call hbname(31,'NUCLEON',ph_NU,  'ph_NU(nNU):R')
      call hbname(31,'NUCLEON',q_NU,   'q_NU(nNU):R')
      call hbname(31,'NUCLEON',acc_NU, 'acc_NU(nNU):R')
      call hbname(31,'NUCLEON',x_NU,   'x_NU(nNU):R')
      call hbname(31,'NUCLEON',y_NU,   'y_NU(nNU):R')
      call hbname(31,'NUCLEON',z_NU,   'z_NU(nNU):R')
      if (beam_type) then
      call hbname(31,'NUCLEON',th1_NU, 'th1_NU(nNU):R')
      call hbname(31,'NUCLEON',th2_NU, 'th2_NU(nNU):R')
      endif
      call hbname(31,'DEUTERON',nD,    'nD[0,5]:I')
      call hbname(31,'DEUTERON',E_D,   'E_D(nD):R')
      call hbname(31,'DEUTERON',P_D,   'P_D(nD):R')
      call hbname(31,'DEUTERON',B_D,   'B_D(nD):R')
      call hbname(31,'DEUTERON',M2_D,  'M2_D(nD):R')
      call hbname(31,'DEUTERON',th_D,  'th_D(nD):R')
      call hbname(31,'DEUTERON',ph_D,  'ph_D(nD):R')
      call hbname(31,'DEUTERON',q_D,   'q_D(nD):R')
      call hbname(31,'DEUTERON',acc_D, 'acc_D(nD):R')
      call hbname(31,'DEUTERON',x_D,   'x_D(nD):R')
      call hbname(31,'DEUTERON',y_D,   'y_D(nD):R')
      call hbname(31,'DEUTERON',z_D,   'z_D(nD):R')

      call hbname(31,'GAMMA',nG,      'nG[0,10]:I')

      call hbname(31,'PION',npi_plus, 'npi_plus:I')
      call hbname(31,'PION',npi_minus,'npi_minus:I')
      call hbname(31,'PION',nPI,      'nPI[0,10]:I')
      call hbname(31,'PION',E_PI,     'E_PI(nPI):R')
      call hbname(31,'PION',P_PI,     'P_PI(nPI):R')
      call hbname(31,'PION',B_PI,     'B_PI(nPI):R')
      call hbname(31,'PION',M2_PI,    'M2_PI(nPI):R')
      call hbname(31,'PION',th_PI,    'th_PI(nPI):R')
      call hbname(31,'PION',ph_PI,    'ph_PI(nPI):R')
      call hbname(31,'PION',q_PI,     'q_PI(nPI):R')
      call hbname(31,'PION',acc_PI,   'acc_PI(nPI):R')
      call hbname(31,'PION',x_PI,     'x_PI(nPI):R')
      call hbname(31,'PION',y_PI,     'y_PI(nPI):R')
      call hbname(31,'PION',z_PI,     'z_PI(nPI):R')
      call hbname(31,'PION',ECtot_PI, 'ECtot_PI(nPI):R')
      call hbname(31,'PION',ECin_PI,  'ECin_PI(nPI):R')
      call hbname(31,'PION',ECout_PI, 'ECout_PI(nPI):R')
      call hbname(31,'ana_inc_pi',nPI_CM,   'nPI_CM[0,10]:I')
      call hbname(31,'ana_inc_pi',E_PI_CM,  'E_PI_CM(nPI_CM):R')
      call hbname(31,'ana_inc_pi',P_PI_CM,  'P_PI_CM(nPI_CM):R')
      call hbname(31,'ana_inc_pi',th_PI_CM, 'th_PI_CM(nPI_CM):R')
      call hbname(31,'ana_inc_pi',ph_PI_CM, 'ph_PI_CM(nPI_CM):R')
      call hbname(31,'ana_inc_pi',P_hort,   'P_hort(nPI_CM):R')
      call hbname(31,'DT_ELNU',neN,   'neN[0,1]:I')
      call hbname(31,'DT_ELNU',M2_eN, 'M2_eN(neN):R')
      call hbname(31,'DT_ep',nep,     'nep[0,3]:I')
      call hbname(31,'DT_ep',M2_ep,   'M2_ep (nep):R')

ccccccccccccccccccc
      if(hel_ana.eq.1) then
         call hbnt  (311,'HEL',' ')  
         call hbname(311,'GENERAL',n_ev_run,   'n_ev_run:I')
         call hbname(311,'GENERAL',beam_hel,   'beam_hel:I')
         call hbname(311,'hel_match',mult,   'mult:I')
         call hbname(311,'hel_match',dose,   'dose:I')
         call hbname(311,'hel_match',el_to_dose,'el_to_dose:R')
         call hbname(311,'hel_match',ltime,   'ltime:R')
      endif

c---------------------------------------------------------------------
c Omega Specific channel ntuples #30
c---------------------------------------------------------------------

      ELSEIF (ana_chan.eq.5) THEN

      call hbnt  (30,'OMEGA',' ')
      
      call hbname(30,'ana_omg',indtype,        'indtype[0,5]:I')
      call hbname(30,'ana_omg',n_ev_run_a, 'n_ev_run(indtype):I')
      call hbname(30,'ana_omg',trig_type_a,'trig_type(indtype):I')
      call hbname(30,'ana_omg',trig_clas_a,'trig_clas(indtype):I')
      call hbname(30,'ana_omg',n_ev_nt_a,  'n_ev_nt(indtype):I')
      call hbname(30,'ana_omg',n_part_a,   'n_part(indtype):I')
      call hbname(30,'ana_omg',ntot_a,     'ntot(indtype):I')
      call hbname(30,'ana_omg',bit_a,      'bit(indtype):I')
      call hbname(30,'ana_omg',bitomega_a,   'bitomega(indtype):I')
      call hbname(30,'ana_omg',sec_coin_SC_a,  'sec_SC(indtype):I')
      call hbname(30,'ana_omg',sec_coin_ST_a,  'sec_ST(indtype):I')

      IF(beam_type) THEN
      call hbname(30,'ana_omg',omega_a,       'omega(indtype):R')
      call hbname(30,'ana_omg',Q2_a,          'Q2(indtype):R')
      call hbname(30,'ana_omg',W_EL_a,        'W(indtype):R')
      call hbname(30,'ana_omg',eps_a,         'eps(indtype):R')
      call hbname(30,'ana_omg',eps_l_a,       'eps_l(indtype):R')
      call hbname(30,'ana_omg',Flux_vp_a,     'Flux_vp(indtype):R')
      call hbname(30,'ana_omg',E_EL_a,        'E_EL(indtype):R')
      call hbname(30,'ana_omg',P_EL_a,        'P_EL(indtype):R')
      call hbname(30,'ana_omg',B_EL_a,        'B_EL(indtype):R')
      call hbname(30,'ana_omg',th_EL_a,       'th_EL(indtype):R')
      call hbname(30,'ana_omg',ph_EL_a,       'ph_EL(indtype):R')
      call hbname(30,'ana_omg',acc_EL_a,      'acc_EL(indtype):R')
      call hbname(30,'ana_omg',x_EL_a,        'x_EL(indtype):R')
      call hbname(30,'ana_omg',y_EL_a,        'y_EL(indtype):R')
      call hbname(30,'ana_omg',z_EL_a,        'z_EL(indtype):R')
      call hbname(30,'ana_omg',ECtot_EL_a,    'ECtot_EL(indtype):R')
      call hbname(30,'ana_omg',ECin_EL_a,     'ECin_EL(indtype):R')
      call hbname(30,'ana_omg',ECout_EL_a,    'ECout_EL(indtype):R')
      ELSE
      call hbname(30,'ana_omg',Egamma_a,      'Egamma(indtype):R')
      call hbname(30,'ana_omg',W_ph_a,        'W(indtype):R')
      call hbname(30,'ana_omg',Tag_time_a,    'Tag_time(indtype):R')
      call hbname(30,'ana_omg',ST_time_a,     'ST_time(indtype):R')
      call hbname(30,'ana_omg',RF1_time_a,    'RF1_time(indtype):R')
      call hbname(30,'ana_omg',RF2_time_a,    'RF2_time(indtype):R')
      call hbname(30,'ana_omg',hit_in_time_a, 'hit_in_time(indtype):I')
      call hbname(30,'ana_omg',chan_a,'chanMC(indtype):I')
      call hbname(30,'ana_omg',true_a,'rand(indtype):I')
      ENDIF

c Proton in lab
      call hbname(30,'ana_omg',det_P,      'det_P(indtype):I')
      call hbname(30,'ana_omg',E_P_lab,    'E_P_lab(indtype):R')
      call hbname(30,'ana_omg',P_P_lab,    'P_P_lab(indtype):R')
      call hbname(30,'ana_omg',B_P_lab,    'B_P_lab(indtype):R')
      call hbname(30,'ana_omg',M2_P_lab,   'M2_P_lab(indtype):R')
      call hbname(30,'ana_omg',th_P_lab,   'th_P_lab(indtype):R')
      call hbname(30,'ana_omg',ph_P_lab,   'ph_P_lab(indtype):R')
      call hbname(30,'ana_omg',acc_P_lab,  'acc_P_lab(indtype):R')
c      call hbname(30,'ana_omg',x_P_lab,    'x_P_lab(indtype):R')
c      call hbname(30,'ana_omg',y_P_lab,    'y_P_lab(indtype):R')
      call hbname(30,'ana_omg',z_P_lab,    'z_P_lab(indtype):R')
      if (beam_type) then
      call hbname(30,'ana_omg',th1_P_lab,  'th1_P_lab(indtype):R')
      call hbname(30,'ana_omg',th2_P_lab,  'th2_P_lab(indtype):R')
      endif
c gamma 1 in lab
      call hbname(30,'ana_omg',det_G1,       'det_G1(indtype):I')
      call hbname(30,'ana_omg',E_G1_lab,     'E_G1_lab(indtype):R')
      call hbname(30,'ana_omg',th_G1_lab,    'th_G1_lab(indtype):R')
      call hbname(30,'ana_omg',ph_G1_lab,    'ph_G1_lab(indtype):R')
c      call hbname(30,'ana_omg',acc_G1_lab,   'acc_G1_lab(indtype):R')
c      call hbname(30,'ana_omg',x_G1_lab,     'x_G1_lab(indtype):R')
c      call hbname(30,'ana_omg',y_G1_lab,     'y_G1_lab(indtype):R')
c      call hbname(30,'ana_omg',z_G1_lab,     'z_G1_lab(indtype):R')

c gamma 2 in lab
      call hbname(30,'ana_omg',det_G2,       'det_G2(indtype):I')
      call hbname(30,'ana_omg',E_G2_lab,     'E_G2_lab(indtype):R')
      call hbname(30,'ana_omg',th_G2_lab,    'th_G2_lab(indtype):R')
      call hbname(30,'ana_omg',ph_G2_lab,    'ph_G2_lab(indtype):R')
c      call hbname(30,'ana_omg',acc_G2_lab,   'acc_G2_lab(indtype):R')
c      call hbname(30,'ana_omg',x_G2_lab,     'x_G2_lab(indtype):R')
c      call hbname(30,'ana_omg',y_G2_lab,     'y_G2_lab(indtype):R')
c      call hbname(30,'ana_omg',z_G2_lab,     'z_G2_lab(indtype):R')
c pi0 in lab
      call hbname(30,'ana_omg',det_PI0,      'det_PI0(indtype):I')
      call hbname(30,'ana_omg',E_PI0_lab,    'E_PI0_lab(indtype):R')
      call hbname(30,'ana_omg',P_PI0_lab,    'P_PI0_lab(indtype):R')
      call hbname(30,'ana_omg',B_PI0_lab,    'B_PI0_lab(indtype):R')
      call hbname(30,'ana_omg',M2_PI0_lab,   'M2_PI0_lab(indtype):R')
      call hbname(30,'ana_omg',th_PI0_lab,   'th_PI0_lab(indtype):R')
      call hbname(30,'ana_omg',ph_PI0_lab,   'ph_PI0_lab(indtype):R')

c pi+ in lab
      call hbname(30,'ana_omg',det_PIp,      'det_PIp(indtype):I')
      call hbname(30,'ana_omg',E_PIp_lab,    'E_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',P_PIp_lab,    'P_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',B_PIp_lab,    'B_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',M2_PIp_lab,   'M2_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',th_PIp_lab,   'th_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',ph_PIp_lab,   'ph_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',acc_PIp_lab,  'acc_PIp_lab(indtype):R')
c      call hbname(30,'ana_omg',x_PIp_lab,    'x_PIp_lab(indtype):R')
c      call hbname(30,'ana_omg',y_PIp_lab,    'y_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',z_PIp_lab,    'z_PIp_lab(indtype):R')
      call hbname(30,'ana_omg',ECtot_PIp_lab,'ECtot_PIp_lab(indtype):R')
c      call hbname(30,'ana_omg',ECin_PIp_lab, 'ECin_PIp_lab(indtype):R')
c      call hbname(30,'ana_omg',ECout_PIp_lab,'ECout_PIp_lab(indtype):R')

c pi- in lab
      call hbname(30,'ana_omg',det_PIm,      'det_PIm(indtype):I')
      call hbname(30,'ana_omg',E_PIm_lab,    'E_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',P_PIm_lab,    'P_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',B_PIm_lab,    'B_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',M2_PIm_lab,   'M2_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',th_PIm_lab,   'th_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',ph_PIm_lab,   'ph_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',acc_PIm_lab,  'acc_PIm_lab(indtype):R')
c      call hbname(30,'ana_omg',x_PIm_lab,    'x_PIm_lab(indtype):R')
c      call hbname(30,'ana_omg',y_PIm_lab,    'y_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',z_PIm_lab,    'z_PIm_lab(indtype):R')
      call hbname(30,'ana_omg',ECtot_PIm_lab,'ECtot_PIm_lab(indtype):R')
c      call hbname(30,'ana_omg',ECin_PIm_lab, 'ECin_PIm_lab(indtype):R')
c      call hbname(30,'ana_omg',ECout_PIm_lab,'ECout_PIm_lab(indtype):R')
c Missing and invariant eNp combinations
      call hbname(30,'ana_omg',M2_eNp_a,   'M2_eNp(indtype):R')
      call hbname(30,'ana_omg',W2_Np_a,    'W2_Np(indtype):R')
c Missing and invariant eNpp combinations
      call hbname(30,'ana_omg',M2_eNpp_a,  'M2_eNpp(indtype):R')
      call hbname(30,'ana_omg',E_eNpp_a,   'E_eNpp(indtype):R')
      call hbname(30,'ana_omg',W2_Npp_a,   'W2_Npp(indtype):R')

c Omega  in lab and in hcm
      call hbname(30,'ana_omg',M2_omega,     'M2_omega(indtype):R')
      call hbname(30,'ana_omg',E_omega_lab,  'E_omega_lab(indtype):R')
      call hbname(30,'ana_omg',P_omega_lab,  'P_omega_lab(indtype):R')
      call hbname(30,'ana_omg',th_omega_lab, 'th_omega_lab(indtype):R')
      call hbname(30,'ana_omg',ph_omega_lab, 'ph_omega_lab(indtype):R')
      call hbname(30,'ana_omg',E_omega_hcm,  'E_omega_hcm(indtype):R')
      call hbname(30,'ana_omg',P_omega_hcm,  'P_omega_hcm(indtype):R')
      call hbname(30,'ana_omg',th_omega_hcm, 'th_omega_hcm(indtype):R')
      call hbname(30,'ana_omg',ph_omega_hcm, 'ph_omega_hcm(indtype):R')
c pi+ in omega rest frame
      call hbname(30,'ana_omg',E_PIp_omega, 'E_PIp_omega(indtype):R')
      call hbname(30,'ana_omg',P_PIp_omega, 'P_PIp_omega(indtype):R')
      call hbname(30,'ana_omg',th_PIp_omega,'th_PIp_omega(indtype):R')
      call hbname(30,'ana_omg',ph_PIp_omega,'ph_PIp_omega(indtype):R')
      call hbname(30,'ana_omg',psi_omega,   'psi_omega(indtype):R')

ccccccccccccccccccc
      ENDIF ! ending specific analysis


c----------------------------------------------------------------------
c Photon histograms histograms 
c----------------------------------------------------------------------

c+ TOF calibration Histos
      IF(.not.beam_type)THEN
      CALL HBOOK1(7001,'TAGGER recalib: tvsst ', 200, -20.0, 20.0, 0.0)
      CALL HBOOK1(7002,'Overall TOF recalib: tag_tof_corr ', 200, -20.0, 20.0, 0.0)
      
      CALL HBOOK2(7110,'TOF rec: SEB pions  Beta % p ',     100,0.,3.0,  100,  0.,1.4,  0.0)
      CALL HBOOK1(7111,'TOF rec: SEB pions  Mass2 ',        100, -0.1, 0.2, 0.0)
      CALL HBOOK2(7112,'TOF rec: SEB  Mass % pad ',     300,0.,300.,  100,  0.,1.2,  0.0)
      CALL HBOOK1(7109,'TOF rec: SEB  Mass ',        100, 0., 1.2, 0.0)
      CALL HBOOK2(7208,'TOF rec: OUR  Beta % p ',     100,0.,3.0,  100,  0.,1.4,  0.0)
      CALL HBOOK1(7209,'TOF rec: OUR  Mass ',        100, 0., 1.2, 0.0)
      CALL HBOOK2(7212,'TOF rec: OUR  Mass % pad ',      300,0.,300.,  100,  0.,1.2,  0.0)
      CALL HBOOK2(7101,'TOF rec: SEB t0 pi % pad SECTOR 1',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7102,'TOF rec: SEB t0 pi % pad SECTOR 2',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7103,'TOF rec: SEB t0 pi % pad SECTOR 3',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7104,'TOF rec: SEB t0 pi % pad SECTOR 4',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7105,'TOF rec: SEB t0 pi % pad SECTOR 5',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7106,'TOF rec: SEB t0 pi % pad SECTOR 6',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7121,'TOF rec: SEB <t0> pi % pad SECTOR 1',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7122,'TOF rec: SEB <t0> pi % pad SECTOR 2',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7123,'TOF rec: SEB <t0> pi % pad SECTOR 3',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7124,'TOF rec: SEB <t0> pi % pad SECTOR 4',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7125,'TOF rec: SEB <t0> pi % pad SECTOR 5',  48,0.,48.,  100,-10.,10,  0.0)
      CALL HBOOK2(7126,'TOF rec: SEB <t0> pi % pad SECTOR 6',  48,0.,48.,  100,-10.,10,  0.0)
      DO n_sec = 1,6
      DO n_pad = 1,48
        CALL HBOOK1(n_sec*1000+n_pad,'Timing Tof by Tof', 100,-10.,10,  0.0) 
      ENDDO
      ENDDO   
      ENDIF
      


c----------------------------------------------------------------------
c Quality check histograms 
c----------------------------------------------------------------------

      IF(mode_zcortype.eq.1 .or. mode_zcortype.eq.0) THEN
        CALL HBOOK1(161,'Zvertex- (Sect1)',60, -30.,30., 0.0)
        CALL HBOOK1(162,'Zvertex- (Sect2)',60, -30.,30., 0.0)
        CALL HBOOK1(163,'Zvertex- (Sect3)',60, -30.,30., 0.0)
        CALL HBOOK1(164,'Zvertex- (Sect4)',60, -30.,30., 0.0)
        CALL HBOOK1(165,'Zvertex- (Sect5)',60, -30.,30., 0.0)
        CALL HBOOK1(166,'Zvertex- (Sect6)',60, -30.,30., 0.0)
        CALL HBOOK1(171,'Zvertex+ (Sect1)',60, -30.,30., 0.0)
        CALL HBOOK1(172,'Zvertex+ (Sect2)',60, -30.,30., 0.0)
        CALL HBOOK1(173,'Zvertex+ (Sect3)',60, -30.,30., 0.0)
        CALL HBOOK1(174,'Zvertex+ (Sect4)',60, -30.,30., 0.0)
        CALL HBOOK1(175,'Zvertex+ (Sect5)',60, -30.,30., 0.0)
        CALL HBOOK1(176,'Zvertex+ (Sect6)',60, -30.,30., 0.0)
      ELSEIF(mode_zcortype.eq.2) THEN
        CALL HBOOK1(161,'Zvertex- (Sect1)',60, -5.,5., 0.0)
        CALL HBOOK1(162,'Zvertex- (Sect2)',60, -5.,5., 0.0)
        CALL HBOOK1(163,'Zvertex- (Sect3)',60, -5.,5., 0.0)
        CALL HBOOK1(164,'Zvertex- (Sect4)',60, -5.,5., 0.0)
        CALL HBOOK1(165,'Zvertex- (Sect5)',60, -5.,5., 0.0)
        CALL HBOOK1(166,'Zvertex- (Sect6)',60, -5.,5., 0.0)
        CALL HBOOK1(171,'Zvertex+ (Sect1)',60, -5.,5., 0.0)
        CALL HBOOK1(172,'Zvertex+ (Sect2)',60, -5.,5., 0.0)
        CALL HBOOK1(173,'Zvertex+ (Sect3)',60, -5.,5., 0.0)
        CALL HBOOK1(174,'Zvertex+ (Sect4)',60, -5.,5., 0.0)
        CALL HBOOK1(175,'Zvertex+ (Sect5)',60, -5.,5., 0.0)
        CALL HBOOK1(176,'Zvertex+ (Sect6)',60, -5.,5., 0.0)
      ELSE
        print *,' bad mode_zcortype=',mode_zcortype
        stop
      ENDIF


      IF (qc_flag) THEN

        CALL HBOOK1(401,'det(beam e)/Qtot OUR',   Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(402,'det(other e)/Qtot OUR',  Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(403,'det(other e+)/Qtot OUR', Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(404,'det(p)/Qtot OUR',        Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(405,'det(n)/Qtot OUR',        Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(406,'det(Pi+)/Qtot OUR',      Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(407,'det(Pi-)/Qtot OUR',      Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(408,'det(K+)/Qtot OUR',       Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(409,'det(K-)/Qtot OUR',       Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(410,'det(D-)/Qtot OUR',       Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(411,'det(G)/Qtot OUR',       Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (401)
        CALL HBARX (402)
        CALL HBARX (403)
        CALL HBARX (404)
        CALL HBARX (405)
        CALL HBARX (406)
        CALL HBARX (407)
        CALL HBARX (408)
        CALL HBARX (409)
        CALL HBARX (410)
        CALL HBARX (411)


        CALL HBOOK1(101,'Inclusive e/Qtot',20000, 0.,20000.,0.0)
        CALL HBOOK1(102,'Inclusive ep/Qtot',20000, 0.,20000.,0.0)
        CALL HBOOK1(103,'Elastic e/Qtot (W cut)',20000, 0.,20000.,0.0)
        CALL HBOOK1(104,'Elastic ep/Qtot (ep coinc, W ph th cuts)',20000, 0.,20000.,0.0)
        CALL HBOOK1(105,'Delta ep/Qtot (ep coinc, W cut)',20000, 0.,20000.,0.0)
        CALL HBOOK1(106,'2pions/Qtot(eppi+ pi-(miss))',20000, 0.,20000.,0.0)
	CALL HBOOK1(107,'Live Time',20000, 0.,20000.,0.0)
        CALL HBOOK1(108,'Ntot per event',Num_files, 0.,Num_files*1.,0.0)
	CALL HBOOK1(109,'Nevnt per event',Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (101)
        CALL HBARX (102)
	CALL HBARX (103)
        CALL HBARX (104)
	CALL HBARX (105)
        CALL HBARX (106)
        CALL HBARX (108)
        CALL HBARX (109)

        CALL HBOOK1(190,'Qgated',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(191,'Nevents',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(192,'Nevents acc. el',Num_files, 0.,Num_files*1.,0.0)

        CALL HBOOK2(203,'ElectronBeam (P-B-OUR)',  100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(204,'ElectronOther (P-B-OUR)', 100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(205,'PositronOther (P-B-OUR)', 100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(206,'PiPlus (P-B-OUR)',        100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(207,'PiMinus (P-B-OUR)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(208,'KPlus (P-B-OUR)',         100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(209,'KMinus (P-B-OUR)',        100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(210,'Proton (P-B-OUR)',        100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(211,'Neutron (P-B-OUR)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(212,'Deutron (P-B-OUR)',       100,0.,5.,  100,0.,1.2,  0.0)      
        CALL HBOOK2(213,'All part(P-B-OUR)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(263,'ElectronBeam (P-B-EVNT)', 100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(264,'ElectronOther (P-B-EVNT)',100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(265,'PositronOther (P-B-EVNT)',100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(266,'PiPlus (P-B-EVNT)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(267,'PiMinus (P-B-EVNT)',      100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(268,'KPlus (P-B-EVNT)',        100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(269,'KMinus (P-B-EVNT)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(270,'Proton (P-B-EVNT)',       100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(271,'Neutron (P-B-EVNT)',      100,0.,5.,  100,0.,1.2,  0.0)
        CALL HBOOK2(272,'Deutron (P-B-EVNT)',      100,0.,5.,  100,0.,1.2,  0.0)      
        CALL HBOOK2(273,'All part(P-B-EVNT)',      100,0.,5.,  100,0.,1.2,  0.0)

        CALL HBOOK1(223,'ElectronBeam Mass-OUR',  100,0.,0.001,  0.0)
        CALL HBOOK1(224,'ElectronOther Mass-OUR', 100,0.,0.001,  0.0)
        CALL HBOOK1(225,'PositronOther Mass-OUR', 100,0.,0.001,  0.0)
        CALL HBOOK1(226,'PiPlus Mass-OUR',        100,0.,0.25,   0.0)
        CALL HBOOK1(227,'PiMinus Mass-OUR',       100,0.,0.25,   0.0)
        CALL HBOOK1(228,'KPlus Mass-OUR',         100,0.,1.,     0.0)
        CALL HBOOK1(229,'KMinus Mass-OUR',        100,0.,1.,     0.0)
        CALL HBOOK1(230,'Proton Mass-OUR',        100,0.6,1.4,   0.0)
        CALL HBOOK1(231,'Neutron Mass-OUR',       100,0.6,1.4,   0.0)
        CALL HBOOK1(232,'Deutron Mass-OUR',       100,0.,4.,     0.0)      
        CALL HBOOK1(233,'All part Mass-OUR',      100,0.05,2.2,  0.0)
        CALL HBOOK1(243,'ElectronBeam Mass-EVNT', 100,0.,0.001,  0.0)
        CALL HBOOK1(244,'ElectronOther Mass-EVNT',100,0.,0.001,  0.0)
        CALL HBOOK1(245,'PositronOther Mass-EVNT',100,0.,0.001,  0.0)
        CALL HBOOK1(246,'PiPlus Mass-EVNT',       100,0.,0.25,   0.0)
        CALL HBOOK1(247,'PiMinus Mass-EVNT',      100,0.,0.25,   0.0)
        CALL HBOOK1(248,'KPlus Mass-EVNT',        100,0.,1.,     0.0)
        CALL HBOOK1(249,'KMinus Mass-EVNT',       100,0.,1.,     0.0)
        CALL HBOOK1(250,'Proton Mass-EVNT',       100,0.6,1.4,   0.0)
        CALL HBOOK1(251,'Neutron Mass-EVNT',      100,0.6,1.4,   0.0)
        CALL HBOOK1(252,'Deutron Mass-EVNT',      100,0.,4.,     0.0)      
        CALL HBOOK1(253,'All part Mass-EVNT',     100,0.05,2.2,  0.0)


        CALL HBOOK1(111,'Elastic Peak Position (Sect1)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(112,'Elastic Peak Position (Sect2)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(113,'Elastic Peak Position (Sect3)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(114,'Elastic Peak Position (Sect4)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(115,'Elastic Peak Position (Sect5)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(116,'Elastic Peak Position (Sect6)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (111)
        CALL HBARX (112)
        CALL HBARX (113)
        CALL HBARX (114)
        CALL HBARX (115)
        CALL HBARX (116)
       
        CALL HBOOK1(121,'Elastic Peak Width (Sect1)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(122,'Elastic Peak Width (Sect2)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(123,'Elastic Peak Width (Sect3)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(124,'Elastic Peak Width (Sect4)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(125,'Elastic Peak Width (Sect5)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(126,'Elastic Peak Width (Sect6)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (121)
        CALL HBARX (122)
        CALL HBARX (123)
        CALL HBARX (124)
        CALL HBARX (125)
        CALL HBARX (126)
        
        CALL HBOOK1(131,'Elastic Peak (Sect1)',30, 0.8,1.0, 0.0)
        CALL HBOOK1(132,'Elastic Peak (Sect2)',30, 0.8,1.0, 0.0)
        CALL HBOOK1(133,'Elastic Peak (Sect3)',30, 0.8,1.0, 0.0)
        CALL HBOOK1(134,'Elastic Peak (Sect4)',30, 0.8,1.0, 0.0)
        CALL HBOOK1(135,'Elastic Peak (Sect5)',30, 0.8,1.0, 0.0)
        CALL HBOOK1(136,'Elastic Peak (Sect6)',30, 0.8,1.0, 0.0)

        CALL HBOOK1(141,'Zvert Position (Sect1)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(142,'Zvert Position (Sect2)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(143,'Zvert Position (Sect3)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(144,'Zvert Position (Sect4)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(145,'Zvert Position (Sect5)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(146,'Zvert Position (Sect6)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (141)
        CALL HBARX (142)
        CALL HBARX (143)
        CALL HBARX (144)
        CALL HBARX (145)
        CALL HBARX (146)

        CALL HBOOK1(151,'Zvert Width (Sect1)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(152,'Zvert Width (Sect2)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(153,'Zvert Width (Sect3)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(154,'Zvert Width (Sect4)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(155,'Zvert Width (Sect5)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(156,'Zvert Width (Sect6)',Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (151)
        CALL HBARX (152)
        CALL HBARX (153)
        CALL HBARX (154)
        CALL HBARX (155)
        CALL HBARX (156)

c        CALL HBOOK1(161,'Zvertex (Sect1)',30, -7.,7., 0.0)
c        CALL HBOOK1(162,'Zvertex (Sect2)',30, -7.,7., 0.0)
c        CALL HBOOK1(163,'Zvertex (Sect3)',30, -7.,7., 0.0)
c        CALL HBOOK1(164,'Zvertex (Sect4)',30, -7.,7., 0.0)
c        CALL HBOOK1(165,'Zvertex (Sect5)',30, -7.,7., 0.0)
c        CALL HBOOK1(166,'Zvertex (Sect6)',30, -7.,7., 0.0)


        CALL HBOOK2(301,'ELECbeam-OUR (P-Etot)',   100,0.,5.,   100,0.,1.2, 0.0)
        CALL HBOOK2(302,'ELECother-OUR (P-Etot)',  100,0.,5.,   100,0.,1.2, 0.0)
        CALL HBOOK2(303,'POSIT-OUR (P-Etot)',      100,0.,5.,   100,0.,1.2, 0.0)
        CALL HBOOK2(304,'PiP-OUR  (P-Etot)',       100,0.,5.,   100,0.,1.2, 0.0)
        CALL HBOOK2(305,'PiM-OUR  (P-Etot)',       100,0.,5.,   100,0.,1.2, 0.0)
        CALL HBOOK2(311,'ELECbeam-OUR (Ein-Eout)', 100,0.,0.8,  100,0.,0.8, 0.0)
        CALL HBOOK2(312,'ELECother-OUR (Ein-Eout)',100,0.,0.8,  100,0.,0.8, 0.0)
        CALL HBOOK2(313,'POSIT-OUR (Ein-Eout)',    100,0.,0.8,  100,0.,0.8, 0.0)
        CALL HBOOK2(314,'PiP-OUR  (Ein-Eout)',     100,0.,0.8,  100,0.,0.8, 0.0)
        CALL HBOOK2(315,'PiM-OUR  (Ein-Eout)',     100,0.,0.8,  100,0.,0.8, 0.0)

        CALL HBOOK2(351,'ELECbeam-EVNT (P-Etot)',   100,0.,5.,  100,0.,1.2, 0.0)
        CALL HBOOK2(352,'ELECother-EVNT (P-Etot)',  100,0.,5.,  100,0.,1.2, 0.0)
        CALL HBOOK2(353,'POSIT-EVNT (P-Etot)',      100,0.,5.,  100,0.,1.2, 0.0)
        CALL HBOOK2(354,'PiP-EVNT  (P-Etot)',       100,0.,5.,  100,0.,1.2, 0.0)
        CALL HBOOK2(355,'PiM-EVNT  (P-Etot)',       100,0.,5.,  100,0.,1.2, 0.0)
        CALL HBOOK2(361,'ELECbeam-EVNT (Ein-Eout)', 100,0.,0.8, 100,0.,0.8, 0.0)
        CALL HBOOK2(362,'ELECother-EVNT (Ein-Eout)',100,0.,0.8, 100,0.,0.8, 0.0)
        CALL HBOOK2(363,'POSIT-EVNT (Ein-Eout)',    100,0.,0.8, 100,0.,0.8, 0.0)
        CALL HBOOK2(364,'PiP-EVNT  (Ein-Eout)',     100,0.,0.8, 100,0.,0.8, 0.0)
        CALL HBOOK2(365,'PiM-EVNT  (Ein-Eout)',     100,0.,0.8, 100,0.,0.8, 0.0)

        DO i=1,9
        write(ch1,701) i
  701   format(i1)
        CALL HBOOK1(800+i,'Inclusive/Ngamma, interval E'//ch1,Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(810+i,'Nprotons/Ngamma, interval E'//ch1, Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(820+i,'Npiplus/Ngamma, interval E'//ch1, Num_files,  0.,Num_files*1.,0.0)
        CALL HBOOK1(830+i,'Npiminus/Ngamma, interval E'//ch1, Num_files, 0.,Num_files*1.,0.0)
        CALL HBOOK1(840+i,'Ntwopions/Ngamma, interval E'//ch1,Num_files, 0.,Num_files*1.,0.0)
        CALL HBARX (800+i)
        CALL HBARX (810+i)
        CALL HBARX (820+i)
        CALL HBARX (830+i)
        CALL HBARX (840+i)
        ENDDO

      ENDIF ! ending Run Quality checks


      END


