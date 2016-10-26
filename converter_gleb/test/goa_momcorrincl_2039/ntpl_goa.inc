

      character*10 MaxEventNum
      integer*4 maxpages
      parameter ( maxpages= 100000)
      parameter ( MaxEventNum='100000000')
     
      integer*4 IQUEST
      integer*4 ipawc
      common /QUEST/IQUEST(100)
      common /pawc/ipawc(maxpages*128)


c---------------------------------------------------------------------
c                    Input Ntuple from CEB
c     !!! ATENTION: The varriavbles ordering in these commons    !!!
c     !!! must identically coincide with the variables ordering  !!!
c     !!! in the calls of HBNAME() in CEB analyser               !!! 
c     !!! (note: it is correspondance with HBNAME() calls in CEB !!!
c     !!!        NOT with the commons in CEB)                    !!!
c---------------------------------------------------------------------

c RUNNUMBER
      INTEGER*4            RUNNumber
      common/com_runnumber/RUNNumber

c RUNGEN
      INTEGER*4     NEventFRUN11
      INTEGER*4     NEventFileFirst, NEventFileLast
      INTEGER*4     NEventTotalRun, NEventThisFile
      REAL*4        xxTMPxx(2000)
      common/g6mapdata/NEventFRUN11
     &                ,NEventFileFirst, NEventFileLast
     &                ,NEventTotalRun,  NEventThisFile
     &                ,xxTMPxx

c SHORTPREP
      INTEGER*4     Ntrigger_file,Ntrigger_sum
      INTEGER*4     Ncebevent_file,Ncebevent_sum
      REAL*4        Qgated_file,Qgated_sum
      common/SHORTPREP/ Ntrigger_file,Ntrigger_sum
     &                 ,Ncebevent_file,Ncebevent_sum
     &                 ,Qgated_file,Qgated_sum

c RUNGAMMA
      REAL*4   NGammaE(767),   NGammaE_2(767)
      REAL*4   NGammaT(61),    NGammaT_2(61)
      REAL*4   ped_SC_R(288),ped_SC_L(288)
      common/com_rungamma/
     &         NGammaE, NGammaE_2
     &        ,NGammaT, NGammaT_2
     &        ,ped_SC_R,ped_SC_L

c NGamma data for g6
      INTEGER*4 Ttag_Nbins, Etag_Nbins
      REAL*4    Ttag_emin(61), Ttag_emax(61)
      REAL*4    Etag_emin(767),Etag_emax(767)
      COMMON/tagger_boundaries/Ttag_Nbins,Ttag_emin,Ttag_emax
     &                        ,Etag_Nbins,Etag_emin,Etag_emax

      INTEGER*4 Trebin_Nbins
      REAL*4    Trebin_emin(61), Trebin_emax(61), Trebin_Ngamma(61)
      INTEGER*4 Erebin_Nbins
      REAL*4    Erebin_emin(767),Erebin_emax(767),Erebin_Ngamma(767)
      COMMON/tagger_rebin/Trebin_Nbins,Trebin_emin,Trebin_emax
     &                   ,Trebin_Ngamma
     &                   ,Erebin_Nbins,Erebin_emin,Erebin_emax
     &                   ,Erebin_Ngamma

      REAL*4                     NgammaTtotal,NgammaTthisfile
      COMMON/ttagger_NgammaTtota/NgammaTtotal,NgammaTthisfile


c head
      INTEGER*4    event_num,event_type,event_clas
      COMMON /HEAD/event_num,event_type,event_clas

      INTEGER*4      nevent
      COMMON /SEB_RE/nevent

      REAL*4        fcg_hevt,tg_hevt
      COMMON /SEB_Q/fcg_hevt,tg_hevt

      INTEGER    helicity
      COMMON/POL/helicity

      INTEGER*4      nHLS,c_BPM,c_FC
      COMMON /HLSnt/ nHLS,c_BPM(40),c_FC(40)

      INTEGER      rx_FBPM,ry_FBPM
      COMMON /FBPM/rx_FBPM,ry_FBPM

      INTEGER*4      nTGBI,heltgbi
      COMMON /TGBInt/nTGBI,heltgbi(4)

      INTEGER*4 n_ph_in,ph_hit_in_time,ph_T_Id_TAGR,ph_E_Id_TAGR
      REAL*4    ph_ERG_TAGR,ph_TPHO_TAGR
      REAL*4    ph_STT_HEVT,ph_RF1_HEVT,ph_RF2_HEVT
      COMMON/PH_IN/
     &  n_ph_in
     & ,ph_ERG_TAGR(1),ph_TPHO_TAGR(1)
     & ,ph_STT_HEVT(1),ph_RF1_HEVT(1),ph_RF2_HEVT(1)
     & ,ph_hit_in_time(1),ph_T_Id_TAGR(1),ph_E_Id_TAGR(1)

      INTEGER       sec_coin_SC
      real*4        thr_trigsc_evnt,id_SC_EVNT,sec_SC_EVNT
      INTEGER       sec_coin_ST
      COMMON/ph_trigger/sec_coin_SC,thr_trigsc_evnt(10),id_SC_EVNT(10),sec_SC_EVNT(10)

      COMMON/PHOTON/ sec_coin_ST

c TAGR
      integer*4 nTAGR,STATUS_TAGR,T_Id_TAGR,E_Id_TAGR 
      real*4    ERG_TAGR,TTAG_TAGR,TPHO_TAGR
      common/TAGRnt/nTAGR 
     &      ,ERG_TAGR(50),TTAG_TAGR(50),TPHO_TAGR(50)
     &      ,STATUS_TAGR(50),T_Id_TAGR(50),E_Id_TAGR(50)

c TGPB      
      REAL*4    Time_TGPB, Energy_TGPB ,dt_TGPB
      integer*4 nTGPB,pointer_TGPB
      common   /TGPBnt/nTGPB,pointer_TGPB(20),Time_TGPB(20)
     &                ,Energy_TGPB(20), dt_TGPB(20) 

      INTEGER*4      ntrk_MVRT
      REAL*4         x_MVRT,y_MVRT,z_MVRT
      COMMON/MVRTnt/ ntrk_MVRT,x_MVRT,y_MVRT,z_MVRT
      

      INTEGER*4 nEVNT
      INTEGER*4 ID_EVNT,ID_OUR
      INTEGER*4 Charge_EVNT
      INTEGER*4 DCstat_EVNT,CCstat_EVNT,SCstat_EVNT
      INTEGER*4 ECstat_EVNT,LCstat_EVNT,STstat_EVNT
      INTEGER*4 Status_EVNT,Trk_flag_EVNT
      REAL*4    th_EVNT,ph_EVNT
      REAL*4    Pmom_EVNT, Pmomx_EVNT, Pmomy_EVNT, Pmomz_EVNT 
      REAL*4    Mass_EVNT, Mass_our, Beta_EVNT, Beta_our
      REAL*4    X_EVNT, Y_EVNT, Z_EVNT

      COMMON /EVNT/nEVNT
     &            ,th_EVNT(10),ph_EVNT(10)
     &            ,ID_EVNT(10),ID_OUR(10)
     &            ,Pmom_EVNT(10)
     &            ,Pmomx_EVNT(10),Pmomy_EVNT(10),Pmomz_EVNT(10)
     &	          ,Mass_EVNT(10),Mass_OUR(10),Charge_EVNT(10)
     &	          ,Beta_EVNT(10),Beta_OUR(10)
     &		  ,X_EVNT(10),Y_EVNT(10),Z_EVNT(10)
     &            ,DCstat_EVNT(10),CCstat_EVNT(10),SCstat_EVNT(10)
     &		  ,ECstat_EVNT(10),LCstat_EVNT(10),STstat_EVNT(10)
     &            ,Status_EVNT(10),Trk_flag_EVNT(10)
  	

      INTEGER  nDCPB
      INTEGER  ScTr_DCPB
      REAL*4   XSC_DCPB, YSC_DCPB, ZSC_DCPB
      REAL*4   CXSC_DCPB,CYSC_DCPB,CZSC_DCPB
      INTEGER  Status_DCPB
      common   /DCPBnt/
     &           nDCPB
     &          ,ScTr_DCPB(20)
     &          ,XSC_DCPB(20),  YSC_DCPB(20),  ZSC_DCPB(20)
     &          ,CXSC_DCPB(20), CYSC_DCPB(20) ,CZSC_DCPB(20)
     &          ,Status_DCPB(20)

      INTEGER  nCCPB
      INTEGER  ScSgHt_CCPB
      REAL*4   Nphe_CCPB
      INTEGER  Status_CCPB
      common   /CCPBnt/
     &          nCCPB
     &         ,ScSgHt_CCPB(20)
     &         ,Nphe_CCPB(20)
     &         ,Status_CCPB(20)

c ECPB
      INTEGER*4 nECPB
      REAL*4    ETOT_ECPB, EIN_ECPB, EOUT_ECPB
      INTEGER*4 Status_ECPB
      REAL*4    M2_hit_ECPB, M3_hit_ECPB, M4_hit_ECPB
      common   /ECPBnt/
     &         nECPB
     &        ,ETOT_ECPB(20) ,EIN_ECPB(20), EOUT_ECPB(20)
     &        ,Status_ECPB(20)
     &        ,M2_hit_ECPB(20), M3_hit_ECPB(20), M4_hit_ECPB(20)

c SCPB
      INTEGER*4 nSCPB      
      REAL*4    Edep_SCPB , Time_SCPB,  Path_SCPB ,Chi2SC_SCPB         
      INTEGER*4 ScPdHt_SCPB,Status_SCPB
      INTEGER*4 sector_SCPB, pd_id_SCPB, hit_id_SCPB,SCstat_SCPB 
      common/SCPBnt/nSCPB
     &             ,SCstat_SCPB(20)
     &             ,ScPdHt_SCPB(20)
     &             ,Time_SCPB(20),Path_SCPB(20)
     &             ,Status_SCPB(20)
     &             ,Edep_SCPB(20)
     &             ,Chi2SC_SCPB(20)
     &             ,sector_SCPB(20), pd_id_SCPB(20), hit_id_SCPB(20)



c SC
      INTEGER*4 nSC,sector_SC
      INTEGER*4 ID_SC,ADCL_SC,ADCR_SC,TDCL_SC,TDCR_SC
      common/SCRnt/nSC,sector_SC(20)
     &            ,ID_SC(20)
     &            ,TDCL_SC(20),ADCL_SC(20),TDCR_SC(20),ADCR_SC(20)


c---------------------------------------------------------------------
c                 Common Parameters
c---------------------------------------------------------------------
c Parameters
      INTEGER*4 trig_event,n_events,n_trigs,ana_gen,ana_chan
      INTEGER*4 num_files,i_file,N_firstfile1
      INTEGER*4 input_type,input_ntdata_num,input_ntmc_num
      REAL*4    E0,Eelbeam,I_torus
      LOGICAL   choice,beam_type,filter_flag,do_recalibrate
      LOGICAL   qc_flag,mode_momcor,mode_zcor
      INTEGER*4 mode_momcortype, mode_zcortype
      REAL*4    Qgated
      INTEGER*4 hel_ana
      INTEGER*4 ChanDoneFlag(1:8,0:1)
      COMMON /el_parms/E0,Eelbeam,n_events,trig_event,n_trigs,choice
     &           ,input_type,input_ntdata_num,input_ntmc_num
     &           ,beam_type,ana_gen,ana_chan,filter_flag,do_recalibrate
     &           ,I_torus,qc_flag,mode_momcor,mode_zcor
     &           ,mode_momcortype,mode_zcortype
     &           ,num_files,i_file,N_firstfile1
     &           ,Qgated
     &	         ,hel_ana
     &           ,ChanDoneFlag
      CHARACTER*100 file_input,file_inputnum
      COMMON /com_file_input/file_input,file_inputnum

c Particle Id
      INTEGER*4 ID_hadr
      REAL*4    mass_hadr,beta_hadr
      COMMON    ID_hadr(10),mass_hadr(10),beta_hadr(10)


c bit_mc
      integer    bit_mc
      common/gpv/bit_mc

c eee
      INTEGER*4  nevent_max
      PARAMETER (nevent_max=1500000)
      LOGICAL   Event_OK(nevent_max)
      INTEGER*4  nbufmax1
      PARAMETER (nbufmax1=50000)
      INTEGER*4 nevent_read,nbuf_read_fcg
      INTEGER*4 indBegin_fcg(nbufmax1), indEnd_fcg(nbufmax1)
      REAL*4    fcg_buf(nbufmax1), fcg_bufplus(nbufmax1)      
      REAL*4    z_el_mid(6),z_p_mid(6)
      INTEGER*4	head(5000),key(5000),key_pair(5000)
      INTEGER*4 key_both(5000),event_head(5000),st_mult(5000)
      INTEGER*4	dose_BPM(5000),dose_FC(5000)
      INTEGER*4 ihel_head
      

      COMMON/events_arrays/ 
     &    nevent_read,nbuf_read_fcg
     &   ,Event_OK
     &   ,indBegin_fcg,indEnd_fcg,fcg_buf,fcg_bufplus
     &   ,z_el_mid,z_p_mid
     &   ,head,key,key_pair,key_both,event_head,st_mult
     &   ,dose_BPM,dose_FC
     &   ,ihel_head



      INTEGER*4 mult,dose
      REAL*4    el_to_dose,ltime
      LOGICAL fill311

      COMMON/hel_match/
     &     mult,dose,el_to_dose,ltime
     &    ,fill311

c---------------------------------------------------------------------
c                 Ntuple 9000 used as common parameters
c---------------------------------------------------------------------
c Common kinematics
      REAL*4    W_EVNT,omega_EVNT,Q2_EVNT,Tp_el_EVNT,thp_el_EVNT
     &         ,Mm2a_EVNT,Mm2b_EVNT,Mm2c_EVNT
     &	       ,s1_EVNT,s2_EVNT,s3_EVNT 
     &	       ,p_miss_EVNT,th_miss_EVNT,ph_miss_EVNT
      COMMON /kinema/  
     &	        W_EVNT,omega_EVNT,Q2_EVNT,Tp_el_EVNT,thp_el_EVNT
     &         ,Mm2a_EVNT,Mm2b_EVNT,Mm2c_EVNT
     &	       ,s1_EVNT,s2_EVNT,s3_EVNT 
     &	       ,p_miss_EVNT,th_miss_EVNT,ph_miss_EVNT
     

	
c---------------------------------------------------------------------
c                 General Ntuple 60
c---------------------------------------------------------------------
c general
      REAL*4        FCGDOSE
      INTEGER*4     N_EV_RUN,TRIG_TYPE,TRIG_CLAS,beam_hel,beam_helkey
     &             ,beam_keypair,beam_keyboth
     &		   ,doplus_BPM,dominus_BPM,doplus_FC,dominus_FC
     &             ,N_EV_NT,N_PART, ntot, bit,sec_coin
      COMMON/N_EV/  N_EV_RUN,TRIG_TYPE,TRIG_CLAS,FCGDOSE
     &             ,beam_hel,beam_helkey,beam_keypair,beam_keyboth
     &		   ,doplus_BPM,dominus_BPM,doplus_FC,dominus_FC
     &             ,N_EV_NT, N_PART, ntot, bit,sec_coin

      

c Qgated = FCGdose for selected events 
      REAL*4        Qgated_total,Qgated_plus
      COMMON/Qgated/Qgated_total,Qgated_plus

c Electron
	REAL*4     Omega, Q2, W, E_EL, P_EL, B_EL, PX_EL, PY_EL, PZ_EL
     &            ,th_EL, ph_EL ,acc_EL
     &            ,X_EL,Y_EL,Z_EL, XMVRT_EL,YMVRT_EL,ZMVRT_EL
     &            ,ECtot_EL,ECin_EL,ECout_EL
     &            ,NpheCC_EL, EffCC_EL
        INTEGER*4  statEVNT_EL,trkEVNT_EL,PdHit_EL
	COMMON/ELECTRON/
     &	           Omega, Q2, W, E_EL, P_EL, B_EL, th_EL, ph_EL
     &            ,acc_EL, PX_EL, PY_EL, PZ_EL
     &            ,X_EL,Y_EL,Z_EL, XMVRT_EL,YMVRT_EL,ZMVRT_EL
     &            ,ECtot_EL,ECin_EL,ECout_EL
     &            ,NpheCC_EL, EffCC_EL
     &            ,statEVNT_EL,trkEVNT_EL
     &            ,PdHit_EL
        REAL*4                eps, eps_l, Flux_vp
        COMMON/Electron_more/ eps, eps_l, Flux_vp

c Nucleon
	INTEGER*4  nNU,id_SC_NU,sec_SC_NU  
	REAL*4     E_NU, P_NU, B_NU, M2_NU ,PX_NU, PY_NU, PZ_NU, th_NU
     &            ,ph_NU, q_NU ,acc_NU,X_NU,Y_NU,Z_NU
     &            ,th1_NU, th2_NU,thrSC_NU
        INTEGER*4  statEVNT_NU,trkEVNT_NU,PdHit_NU
	COMMON/NUCLEON/
     &	           nNU
     &            ,E_NU(1), P_NU(1),B_NU(1), M2_NU(1)
     &            ,th_NU(1), ph_NU(1), q_NU(1) ,acc_NU(1)
     &            ,PX_NU(1), PY_NU(1), PZ_NU(1)
     &            ,X_NU(1),Y_NU(1),Z_NU(1)
     &            ,statEVNT_NU(1),trkEVNT_NU(1)
     &            ,th1_NU(1), th2_NU(1),thrSC_NU(1)
     &            ,PdHit_NU(1),id_SC_NU(1),sec_SC_NU(1)

c Deuteron
	INTEGER*4  nD     
	REAL*4     E_D, P_D, B_D, M2_D ,PX_D, PY_D, PZ_D, th_D
     &            ,ph_D, q_D ,acc_D,X_D,Y_D,Z_D
        INTEGER*4  statEVNT_D,trkEVNT_D,PdHit_D
	COMMON/DEUTERON/
     &	           nD
     &            ,E_D(5), P_D(5), B_D(5), M2_D(5)
     &            ,th_D(5), ph_D(5), q_D(5) ,acc_D(5)
     &            ,PX_D(5), PY_D(5), PZ_D(5)
     &            ,X_D(5),Y_D(5),Z_D(5)
     &            ,statEVNT_D(5),trkEVNT_D(5)
     &            ,PdHit_D(5)

c Pion
        INTEGER*4   npi_plus,npi_minus,nPI,id_SC_PI,sec_SC_PI
        REAL*4      E_PI, P_PI, B_PI,M2_PI, PX_PI, PY_PI, PZ_PI
     &             ,th_PI, ph_PI, q_PI ,acc_PI, E_Cal_PI
     &             ,x_PI,y_PI,z_PI
     &             ,ECtot_PI,ECin_PI,ECout_PI,thrSC_PI
        INTEGER*4   statEVNT_PI,trkEVNT_PI,PdHit_PI
	COMMON/PI/  npi_plus,npi_minus, nPI
     &             ,E_PI(10), P_PI(10),B_PI(10),M2_PI(10)
     &             ,th_PI(10),ph_PI(10),q_PI(10),acc_PI(10)
     &             ,E_Cal_PI(10)
     &             ,PX_PI(10), PY_PI(10), PZ_PI(10)
     &             ,x_PI(10),y_PI(10),z_PI(10)
     &             ,ECtot_PI(10),ECin_PI(10),ECout_PI(10)
     &             ,statEVNT_PI(10),trkEVNT_PI(10)
     &             ,thrSC_PI(10)
     &             ,PdHit_PI(10)
     &             ,id_SC_PI(10),sec_SC_PI(10)

c Kaon
        INTEGER*4 nK
	REAL*4    E_K, P_K,  B_K ,M2_K , PX_K, PY_K, PZ_K, th_K,ph_K
     &           ,q_K, acc_K 
     &           ,x_K,y_K,z_K
        INTEGER*4 statEVNT_K,trkEVNT_K
	COMMON/KAONS/ nK, E_K(10), P_K(10),B_K(10),M2_K(10)
     &               ,th_K(10), ph_K(10), q_K(10) ,acc_K(10)
     &               ,PX_K(10), PY_K(10), PZ_K(10)
     &               ,x_K(10),y_K(10),z_K(10)
     &               ,statEVNT_K(10),trkEVNT_K(10)

c Gamma
        INTEGER*4 nG
	REAL*4    E_G, PX_G, PY_G, PZ_G, th_G, ph_G, acc_G
     &          , x_G,y_G,z_G
        INTEGER*4 PdHit_G
	COMMON/GAMMA/  nG, E_G(10), th_G(10), ph_G(10), acc_G(10)
     &                ,PX_G(10), PY_G(10), PZ_G(10)
     &                ,x_G(10),y_G(10),z_G(10)
     &                ,PdHit_G(10)

c Lepton
        INTEGER*4 nL
	REAL*4    E_L, P_L , B_L,  PX_L, PY_L, PZ_L, th_L , ph_L
     &           ,q_L, id_L, acc_L
     &           ,x_L,y_L,z_L
	COMMON/LEPTONS/ nL, id_L(4), E_L(4), P_L(4), B_L(4)
     &                 ,th_L(4), ph_L(4), q_L(4), acc_L(4)
     &                 ,PX_L(4), PY_L(4), PZ_L(4)
     &                 ,x_L(4),y_L(4),z_L(4)

c GA_GA
c GA_GA
c GA_GA_GA
      INTEGER*4 nGG,nGGG,nGGGG
      REAL*4    W2_GG  ,E_GG , P_GG,  th_GG,  ph_GG 
     &         ,W2_GGG
     &         ,W2_GGGG
      COMMON/Ga_Ga/       nGG,  W2_GG(6),E_GG(6),P_GG(6)
     &                   ,th_GG(6), ph_GG(6)         	
      COMMON/Ga_Ga_Ga/    nGGG, W2_GGG(4)  
      COMMON/Ga_Ga_Ga_Ga/ nGGGG,W2_GGGG(1)  
     
c DETected particles combinations
c missing respect to electron + ... (e.g. e + Nucleon + pi  --> _eNp)
c invariant respect to ...          (e.g.     Nucleon + pi  --> _Np)	     
	INTEGER  neN,nep,nek,neNp,nepp,neNk,nepk,nekk,neNpp,neppp
     &    	,neNkk, neNpk, neNppp
	REAL*4   M2_eN, E_eN , P_eN,  th_eN,  ph_eN,  q_eN 
     &          ,M2_ep, M2_ek, M2_eNp, M2_epp, M2_eNk, M2_epk, M2_ekk 
     &  	,M2_eNpp,E_eNpp
     &  	,M2_eppp, M2_eNkk, M2_eNpk, M2_eNppp
	REAL*4   W2_Np, W2_pp, W2_Nk, W2_pk, W2_kk, W2_Npp, W2_ppp
     &          ,W2_Nkk, W2_Npk, W2_Nppp

	COMMON/DET_eN/    neN,   M2_eN(1),
     &                    E_eN(1), P_eN(1),th_eN(1),ph_eN(1),q_eN(1)
	COMMON/DET_ep/    nep,   M2_ep(3)             	
	COMMON/DET_ek/    nek,   M2_ek(3)             	
	COMMON/DET_eNp/   neNp,  M2_eNp(3),  W2_Np(3)             	
	COMMON/DET_epp/   nepp,  M2_epp(3),  W2_pp(3)           	
	COMMON/DET_eNk/   neNk,  M2_eNk(3),  W2_Nk(3)           	
	COMMON/DET_epk/   nepk,  M2_epk(3),  W2_pk(3)           	
	COMMON/DET_ekk/   nekk,  M2_ekk(3),  W2_kk(3)           	
	COMMON/DET_eNpp/  neNpp, M2_eNpp(3), E_eNpp(3), W2_Npp(3)
	COMMON/DET_eppp/  neppp, M2_eppp(3), W2_ppp(3)            	
	COMMON/DET_eNkk/  neNkk, M2_eNkk(3), W2_Nkk(3)           	
	COMMON/DET_eNpk/  neNpk, M2_eNpk(3), W2_Npk(3)           	
	COMMON/DET_eNppp/ neNppp,M2_eNppp(3),W2_Nppp(3)           	

c Combinations with neutrals particles
	REAL*4   
     &           M2_eNG 
     &          ,M2_e0, M2_eN0
     &          ,M2_ep0,M2_ek0,M2_eNp0,M2_epp0,M2_eNk0,M2_epk0,M2_ekk0 
     &  	,M2_eNpp0, M2_eppp0, M2_eNkk0, M2_eNpk0, M2_eNppp0
	REAL*4   
     &           W2_NG
     &          ,W2_N0, W2_p0, W2_k0
     &          ,W2_Np0, W2_pp0, W2_Nk0, W2_pk0, W2_kk0, W2_Npp0 
     &          ,W2_ppp0, W2_Nkk0, W2_Npk0, W2_Nppp0
	INTEGER*4  
     &    	  neNG
     &    	 ,ne0, neN0, nep0,  nek0, neNp0, nepp0, neNk0, nepk0
     &           ,nekk0, neNpp0, neppp0, neNkk0, neNpk0, neNppp0

c 1 Gamma detected
	COMMON/DET_eNG/     neNG,   M2_eNG(4),   W2_NG(4)
c 2 Gamma detected
	COMMON/DET_e0/      ne0,    M2_e0(1)
	COMMON/DET_eN0/     neN0,   M2_eN0(1),   W2_N0(1)	
	COMMON/DET_ep0/     nep0,   M2_ep0(3),   W2_p0(3)            	
	COMMON/DET_ek0/     nek0,   M2_ek0(3),   W2_k0(3)           	
	COMMON/DET_eNp0/    neNp0,  M2_eNp0(3),  W2_Np0(3)             	
	COMMON/DET_epp0/    nepp0,  M2_epp0(3),  W2_pp0(3)           	
	COMMON/DET_eNk0/    neNk0,  M2_eNk0(3),  W2_Nk0(3)           	
	COMMON/DET_epk0/    nepk0,  M2_epk0(3),  W2_pk0(3)           	
	COMMON/DET_ekk0/    nekk0,  M2_ekk0(3),  W2_kk0(3)           	
	COMMON/DET_eNpp0/   neNpp0, M2_eNpp0(3), W2_Npp0(3)           	
	COMMON/DET_eppp0/   neppp0, M2_eppp0(3), W2_ppp0(3)            	
	COMMON/DET_eNkk0/   neNkk0, M2_eNkk0(3), W2_Nkk0(3)           	
	COMMON/DET_eNpk0/   neNpk0, M2_eNpk0(3), W2_Npk0(3)           	
	COMMON/DET_eNppp0/  neNppp0,M2_eNppp0(3),W2_Nppp0(3)           	

c 2Pions in Lab, Center of mass, Rest frame variables
       integer     ncmrest
       real*4 
     &             E_lab,  P_lab, th_lab ,ph_lab
     &            ,E_cm,P_cm,th_cm,ph_cm
     &            ,E_rest, P_rest, th_rest, ph_rest
       COMMON/twopions/ncmrest
     &            ,E_lab(1),  P_lab(1), th_lab(1) ,ph_lab(1)
     &            ,E_cm(1),   P_cm(1),  th_cm(1)  ,ph_cm(1)
     &            ,E_rest(1), P_rest(1),th_rest(1),ph_rest(1)

c Efficiency
	REAL*4     Mm2_eff,p_miss_eff,th_miss_eff,ph_miss_eff
	COMMON/Efficiency/
     &	           Mm2_eff,p_miss_eff,th_miss_eff,ph_miss_eff 
          
c PHOTON ANALYSIS
        REAL*4         Egamma,Tag_time,ST_time
     &                ,RF1_time, RF2_time 
        INTEGER hit_in_time,chan,true,T_id, E_id
        COMMON/PH_BEAM/Egamma,Tag_time, ST_time
     &                ,RF1_time, RF2_time,hit_in_time,chan,true,T_id, E_id

c Kinemaic fit result
       integer*4 nfit,fit_err
       real*4 fit_p, fit_dp, fit_th, fit_dt, fit_fi, fit_df
     &       ,fit_Eg,fit_DEg,fit_ch2,fit_W
       common/fit/nfit
     &       ,fit_Eg(3),fit_DEg(3)
     &       ,fit_p(3),fit_dp(3), fit_th(3), fit_dt(3), fit_fi(3)
     &       ,fit_df(3)
     &       ,fit_W(3)
     &       ,fit_ch2(3),fit_err(3)




c----------------------------------------------------------------------
c Specific  analysis
c----------------------------------------------------------------------

c General run info
      integer*4 indtype, n_ev_run_a,trig_type_a,trig_clas_a,n_ev_nt_a,n_part_a,ntot_a
     $                 ,sec_coin_SC_a,sec_coin_ST_a,sec_coin_SC_thr_a
     &         ,bit_a,bit1pi_a,bit2pi_a,bitomega_a
      common/ana_genrun/indtype
     &                 ,n_ev_run_a(5), trig_type_a(5),trig_clas_a(5) ,n_ev_nt_a(5)
     &                 ,n_part_a(5) ,ntot_a(5)
     &                 ,bit_a(5), bit1pi_a(5),bit2pi_a(5),bitomega_a(5),sec_coin_SC_a(5)
     &                 ,sec_coin_ST_a(5),sec_coin_SC_thr_a(5)

c Electron beam
       real*4       Omega_a ,Q2_a , W_el_a , eps_a, eps_l_a ,Flux_vp_a
     &             ,E_EL_a, P_EL_a, B_EL_a, th_EL_a, ph_EL_a
     &             ,acc_EL_a
     &             ,x_EL_a,y_EL_a,z_EL_a
     &             ,ECtot_EL_a,ECin_EL_a,ECout_EL_a
     &             ,NpheCC_EL_a,EffCC_EL_a,jopa,joppa
       integer*4    PdHit_EL_a
       common /ana_el/
     &      Omega_a(5) ,Q2_a(5) , W_el_a(5) , eps_a(5), eps_l_a(5)
     &     ,Flux_vp_a(5)
     &     ,E_EL_a(5), P_EL_a(5), B_EL_a(5), th_EL_a(5), ph_EL_a(5)
     &     ,acc_EL_a(5)
     &     ,x_EL_a(5),y_EL_a(5),z_EL_a(5)
     &     ,ECtot_EL_a(5),ECin_EL_a(5),ECout_EL_a(5)
     &     ,NpheCC_EL_a(5),EffCC_EL_a(5)
     &     ,PdHit_EL_a(5),jopa,joppa

c Photon beam
       integer*4   hit_in_time_a,chan_a,true_a,T_id_a, E_id_a
       real*4     Egamma_a,W_ph_a,Tag_time_a,ST_time_a
     &           ,RF1_time_a,RF2_time_a
     #       , trig_SC
       common /ana_ph/
     &     Egamma_a(5),W_ph_a(5),Tag_time_a(5)
     &    ,ST_time_a(5)
     &    ,RF1_time_a(5),RF2_time_a(5), hit_in_time_a(5)
     &    ,chan_a(5),true_a(5)
     #       , trig_SC(5)
     #     ,  T_id_a(5), E_id_a(5)

c Single particles in the lab
       integer det_P,det_n, det_PIp, det_PIm, det_PI0,det_g1,det_g2
     #       ,id_SC_p,sec_SC_P,id_SC_pip,sec_SC_pip,id_SC_pim,sec_SC_pim
       real*4  E_P_lab,   P_P_lab,   B_P_lab,   M2_P_lab  
     #       , th_P_lab,   ph_P_lab
     #       , E_n_lab,   P_n_lab,   B_n_lab,   M2_n_lab
     #       , th_n_lab,   ph_n_lab
     #       , E_PIp_lab, P_PIp_lab, B_PIp_lab, M2_PIp_lab
     #       , th_PIp_lab, ph_PIp_lab
     #       , E_PIm_lab, P_PIm_lab, B_PIm_lab, M2_PIm_lab
     #       , th_PIm_lab, ph_PIm_lab
     #       , E_PI0_lab, P_PI0_lab, B_PI0_lab, M2_PI0_lab
     #       , th_PI0_lab, ph_PI0_lab
     #       , E_G1_lab, th_G1_lab, ph_G1_lab 
     #       , E_G2_lab, th_G2_lab, ph_G2_lab 
     #       , M2_eNpp_a, E_eNpp_a,W2_Npp_a
     #       , M2_eNp_a, W2_Np_a
     #       , acc_P_lab
     #       , acc_n_lab
     #       , acc_PIp_lab
     #       , acc_PIm_lab
     #       , acc_PI0_lab
     #       , acc_G1_lab
     #       , acc_G2_lab
     #       , x_P_lab, y_P_lab, z_P_lab
     #       , x_n_lab, y_n_lab, z_n_lab
     #       , x_PIp_lab, y_PIp_lab, z_PIp_lab
     #       , x_PIm_lab, y_PIm_lab, z_PIm_lab
     #       , x_PI0_lab, y_PI0_lab, z_PI0_lab
     #       , x_G1_lab, y_G1_lab, z_G1_lab
     #       , x_G2_lab, y_G2_lab, z_G2_lab
     #       , ECtot_PIp_lab, ECin_PIp_lab, ECout_PIp_lab
     #       , ECtot_PIm_lab, ECin_PIm_lab, ECout_PIm_lab
     #       , th1_P_lab,th2_P_lab
     #       , M2_eN_a, E_eN_a, P_eN_a, th_eN_a, ph_eN_a
     #       ,thrSC_p,thrSC_pip,thrSC_pim
      integer*4  PdHit_P_lab,PdHit_PIp_lab,PdHit_PIm_lab

       common /ana_p_lab/          
     #   det_P(5)    
     #,  E_P_lab(5), P_P_lab(5), B_P_lab(5), M2_P_lab(5)
     #,  th_P_lab(5),  ph_P_lab(5)
     #,  acc_P_lab(5)
     #,  x_P_lab(5), y_P_lab(5), z_P_lab(5)
     #,  th1_P_lab(5), th2_P_lab(5)
     #,  PdHit_P_lab(5)
     #,  M2_eN_a(5), E_eN_a(5), P_eN_a(5), th_eN_a(5), ph_eN_a(5)
     #,  id_SC_p(5),sec_SC_P(5),thrSC_p(5)
       common /ana_n_lab/           
     #   det_n(5)    
     #,  E_n_lab(5), P_n_lab(5), B_n_lab(5), M2_n_lab(5)
     #,  th_n_lab(5), ph_n_lab(5)
     #,  acc_n_lab(5)
     #,  x_n_lab(5), y_n_lab(5), z_n_lab(5)

       common /ana_pip_lab/ 
     #   det_PIp(5)
     #,  E_PIp_lab(5), P_PIp_lab(5), B_PIp_lab(5), M2_PIp_lab(5)
     #,  th_PIp_lab(5), ph_PIp_lab(5)
     #,  acc_PIp_lab(5)
     #,  x_PIp_lab(5), y_PIp_lab(5), z_PIp_lab(5)
     #,  ECtot_PIp_lab(5), ECin_PIp_lab(5), ECout_PIp_lab(5)
     #,  PdHit_PIp_lab(5)
     #,  id_SC_pip(5),sec_SC_Pip(5),thrSC_pip(5)
       common /ana_pim_lab/ 
     #   det_PIm(5)
     #,  E_PIm_lab(5), P_PIm_lab(5), B_PIm_lab(5), M2_PIm_lab(5)
     #,  th_PIm_lab(5), ph_PIm_lab(5)
     #,  acc_PIm_lab(5)
     #,  x_PIm_lab(5), y_PIm_lab(5), z_PIm_lab(5)
     #,  ECtot_PIm_lab(5), ECin_PIm_lab(5), ECout_PIm_lab(5)
     #,  PdHit_PIm_lab(5)
     #,  id_SC_pim(5),sec_SC_Pim(5),thrSC_pim(5)

       common /ana_pi0_lab/  
     #   det_PI0(5)
     #,  E_PI0_lab(5), P_PI0_lab(5), B_PI0_lab(5), M2_PI0_lab(5)
     #,  th_PI0_lab(5), ph_PI0_lab(5)
     #,  acc_PI0_lab(5)
     #,  x_PI0_lab(5), y_PI0_lab(5), z_PI0_lab(5)


       common /ana_g1_lab/  
     #   det_g1(5)
     #,  E_G1_lab(5), th_G1_lab(5), ph_G1_lab(5)
     #,  acc_G1_lab(5)
     #,  x_G1_lab(5), y_G1_lab(5), z_G1_lab(5)
       common /ana_g2_lab/ 
     #   det_g2(5)
     #,  E_G2_lab(5), th_G2_lab(5), ph_G2_lab(5)
     #,  acc_G2_lab(5)
     #,  x_G2_lab(5), y_G2_lab(5), z_G2_lab(5)

       common /ana_miss_Npp/  
     #  M2_eNpp_a(5), E_eNpp_a(5), W2_Npp_a(5)
       common /ana_miss_Np/ 
     #  M2_eNp_a(5), W2_Np_a(5)
       real*4
     #  M2_eN_a2
     # ,M2_epip_a2
     # ,M2_epim_a2
     # ,M2_eNpip_a2
     # ,M2_eNpim_a2
       common /ana_miss_2a/ 
     #  M2_eN_a2(5)
     # ,M2_epip_a2(5)
     # ,M2_epim_a2(5)
     # ,M2_eNpip_a2(5)
     # ,M2_eNpim_a2(5)

c Single particles in hcm
       real*4  
     # E_PI0_hcm,  P_PI0_hcm,   th_PI0_hcm,  ph_Pi0_hcm
     #,E_PIp_hcm,  P_PIp_hcm,   th_PIp_hcm,  ph_Pip_hcm
       common/ana_pi0_hcm/ 
     #  E_PI0_hcm(5),  P_PI0_hcm(5),   th_PI0_hcm(5),  ph_Pi0_hcm(5)
       common/ana_pip_hcm/ 
     #  E_PIp_hcm(5),  P_PIp_hcm(5),   th_PIp_hcm(5),  ph_Pip_hcm(5)


c 1 pion specific analysis
       integer*4 bit1pi
       common/general_1pi/bit1pi

       real*4  
     # M2_Dp,    E_Dp_lab,    P_Dp_lab,    th_Dp_lab,    ph_Dp_lab
       common/ana_1pi_res/
     # M2_Dp(5), E_Dp_lab(5), P_Dp_lab(5), th_Dp_lab(5), ph_Dp_lab(5)

c 2 pions specific analysis
       integer*4 bit2pi
              common/general_2pi/bit2pi

c omega specific analysis
       integer*4 bitomega
       common/general_omega/bitomega

       real*4  
     # M2_rho, E_rho_lab, P_rho_lab, th_rho_lab, ph_rho_lab    
     #       , E_rho_hcm, P_rho_hcm, th_rho_hcm, ph_rho_hcm
     #,M2_Dpp, E_Dpp_lab, P_Dpp_lab, th_Dpp_lab, ph_Dpp_lab
     #        ,E_Dpp_hcm, P_Dpp_hcm, th_Dpp_hcm, ph_Dpp_hcm
     #,M2_D0,  E_D0_lab,   P_D0_lab,  th_D0_lab,  ph_D0_lab
     #       , E_D0_hcm,   P_D0_hcm,  th_D0_hcm,  ph_D0_hcm 
     #       , E_PIp_rho, P_PIp_rho,  th_PIp_rho, ph_Pip_rho,psi_rho
     #       , E_PIp_Dpp, P_PIp_Dpp,  th_PIp_Dpp, ph_Pip_Dpp,psi_Dpp
     #       , E_PIm_D0,  P_PIm_D0,   th_PIm_D0,  ph_Pim_D0,psi_D0
     #       , E_PIp_omega, P_PIp_omega,  th_PIp_omega, ph_Pip_omega,psi_omega
     #, M2_omega, E_omega_lab, P_omega_lab, th_omega_lab, ph_omega_lab    
     #       , E_omega_hcm, P_omega_hcm, th_omega_hcm, ph_omega_hcm

       common /ana_2pi_res/   
     # M2_rho(5), E_rho_lab(5),P_rho_lab(5),th_rho_lab(5),ph_rho_lab(5)
     #          , E_rho_hcm(5),P_rho_hcm(5),th_rho_hcm(5),ph_rho_hcm(5)
     #,M2_Dpp(5), E_Dpp_lab(5),P_Dpp_lab(5),th_Dpp_lab(5),ph_Dpp_lab(5)
     #           ,E_Dpp_hcm(5),P_Dpp_hcm(5),th_Dpp_hcm(5),ph_Dpp_hcm(5)
     #,M2_D0(5),  E_D0_lab(5),  P_D0_lab(5), th_D0_lab(5),ph_D0_lab(5)
     #          , E_D0_hcm(5),  P_D0_hcm(5), th_D0_hcm(5),ph_D0_hcm(5)
     #          , E_PIp_rho(5),P_PIp_rho(5),th_PIp_rho(5),ph_Pip_rho(5),psi_rho(5)
     #          , E_PIp_Dpp(5),P_PIp_Dpp(5),th_PIp_Dpp(5),ph_Pip_Dpp(5),psi_Dpp(5)
     #          , E_PIm_D0(5), P_PIm_D0(5),  th_PIm_D0(5),ph_Pim_D0(5),psi_D0(5)
     #       , E_PIp_omega(5), P_PIp_omega(5),  th_PIp_omega(5), ph_Pip_omega(5),psi_omega(5)
     #,M2_omega(5), E_omega_lab(5), P_omega_lab(5), th_omega_lab(5), ph_omega_lab(5)    
     #       , E_omega_hcm(5), P_omega_hcm(5), th_omega_hcm(5), ph_omega_hcm(5)

        INTEGER*4 Err_fit
        REAL*4    Eg_fit, DEg_fit
     #          , DP_P_fit, DP_PIp_fit, DP_PIm_fit
     #          , DTH_P_fit, DTH_PIp_fit, DTH_PIm_fit
     #          , DPH_P_fit, DPH_PIp_fit, DPH_PIm_fit
     #          , CHI_fit
        common /ana_2pi_fit/ 
     #            Eg_fit(5) 
     #          , dEg_fit(5) 
     #          , dp_p_fit(5)
     #          , DP_PIp_fit(5), DP_PIm_fit(5)
     #          , DTH_P_fit(5), DTH_PIp_fit(5), DTH_PIm_fit(5)
     #          , DPH_P_fit(5), DPH_PIp_fit(5), DPH_PIm_fit(5)
     #          , CHI_fit(5), Err_fit(5)

c Inlcusive pion specific analysis
	integer*4   nPI_CM
        REAL*4      costh_v,E_PI_CM, P_PI_CM
     &             ,th_PI_CM, ph_PI_CM,P_hort

	COMMON/ana_inc_pi/nPI_CM,cosTH_v
     &		   ,E_PI_CM(10), P_PI_CM(10)
     &             ,th_PI_CM(10), ph_PI_CM(10),P_hort(10)


