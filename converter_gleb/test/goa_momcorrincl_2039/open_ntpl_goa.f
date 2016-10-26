
      SUBROUTINE open_ntpl_goa(InFile)
      
      IMPLICIT NONE  
      include "ntpl_goa.inc"
      include "ntpl_clas.inc"
      integer*4 istat,lrecl,i
      CHARACTER*(*) InFile
      data lrecl/8192/


c ----- OPEN input nt of NEXT file and set -----
      if(input_type.eq.1) lrecl = 8192
      if(input_type.eq.2) lrecl = 0
      call hropen(12,'ESCA',InFile,'P',lrecl,istat)
        IF(istat.ne.0) print *,' ***Can not open InFile:',InFile
        IF(istat.ne.0) print *,' ***istat,lrecl=',istat,lrecl
        IF(istat.ne.0) print *,' ***trying one more time'
        IF(istat.ne.0) lrecl=0
        IF(istat.ne.0) call hrend('ESCA')
        IF(istat.ne.0) call hropen(12,'ESCA',InFile,'P',lrecl,istat)


c-------------- CEB nt ------------------------------
      IF(input_type.eq.1) THEN

      call hrin  (81,9999,0)
      call hbname(81,' ',0,'$CLEAR')
      call hbname(81,'G6',NEventFRUN11,'$SET')     

      call hrin  (99,9999,0)
      call hbname(99,' ',0,'$CLEAR')
      call hbname(99,'HEAD',event_num,'$SET')
      call hbname(99,'POL',helicity,'$SET')
      call hbname(99,'SEB_RE',nevent,'$SET')
      call hbname(99,'SEB_Q',fcg_hevt,'$SET')
      call hbname(99,'EVNT',nEVNT,'$SET')
      call hbname(99,'SEB_DC',nDCPB,'$SET')
      call hbname(99,'SEB_CC',nCCPB,'$SET')
      call hbname(99,'SEB_EC',nECPB,'$SET')
        if(beam_type) then 
      call hbname(99,'SEB_SC',nSCPB,'$SET') ! att. 'SCPB'
        endif
        if(beam_type) then
      call hbname(99,'FBPM',rx_FBPM,'$SET')
      call hbname(99,'TGBI',nTGBI,'$SET')
      call hbname(99,'HLS',nHLS,'$SET')
        endif
        if(.not.beam_type) then 
      call hbname(99,'PH_IN',n_ph_in,'$SET') ! photon case
      call hbname(99,'PHOTON',SEC_COIN_ST,'$SET')
      call hbname(99,'SC',nSC,'$SET')
      call hbname(99,'SCPB',  nSCPB,'$SET') ! att. 'SEB_SC'
cc      call hbname(99,'TAGR',  nTAGR,'$SET')
cc      call hbname(99,'SEB_TG',nTGPB,'$SET')
        endif

      call hrin  (98,9999,0)
      call hbname(98,' ',0,'$CLEAR')
      call hbname(98,'MC_HEAD',event_num,'$SET')
      call hbname(98,'POL',helicity,'$SET')
      call hbname(98,'MC_SEB_R',nevent,'$SET')
      call hbname(98,'MC_EVNT',nEVNT,'$SET')
      call hbname(98,'MC_SEB_E',nECPB,'$SET')
        if(.not.beam_type) then
      call hbname(98,'MC_PH_IN',n_ph_in,'$SET') ! photon case
        endif

c---------------- CLAS nt ------------------------          
      ELSEIF(input_type.eq.2) THEN

      call hrin  (input_ntdata_num,9999,0)
      call hbname(input_ntdata_num,' ',0,'$CLEAR')
      call hbname(input_ntdata_num,'HEVT', clas_npart,   '$SET')
      if(beam_type)      call hbname(input_ntdata_num,'EVNT', clas_gpart,  '$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'EVNT', clasg6_gpart,'$SET')
      call hbname(input_ntdata_num,'DCPB', clas_dc_part, '$SET')
      call hbname(input_ntdata_num,'ECPB', clas_ec_part, '$SET')
      call hbname(input_ntdata_num,'SCPB', clas_sc_part, '$SET')
      call hbname(input_ntdata_num,'CCPB', clas_cc_part, '$SET')
      call hbname(input_ntdata_num,'LECPB',clas_lac_part,'$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'PART',clas_nPRT,        '$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'SC',  clasg6_nSC,       '$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'STPB',clasg6_st_part,   '$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'TAGR',clasg6_taghit,    '$SET')
      if(.not.beam_type) call hbname(input_ntdata_num,'TGPB',clasg6_taggoodhit,'$SET')

      ELSE
      print *,' input_type is bad=',input_type
      stop
      ENDIF


c ----- set some variables -----
      DO i=1,6
      z_el_mid(i)=0.0
      z_p_mid(i)=0.0
      ENDDO

      RETURN
      END
