ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	This program reads BOS format preprocessed CLAS event file,
c	then it does some physics analysis
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      PROGRAM analyser_goa
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"

      CHARACTER*200 Fname(999), file_in,OutFile,OutFile2
      INTEGER*4 ierr98,ierr99,ierr81
      INTEGER*4 ifile,icycle,i,i1,i2
      INTEGER*4 nfiles,runnum
      INTEGER*4 nstory,nstory_MC,nstory_81, NDATA,NMC,ND81
      INTEGER*4 NDATAmax,NevMC
      CHARACTER*1   ch1 
      LOGICAL   firsttime,eof,filter_ok_data,filter_ok_MC,Select_OK

      LOGICAL   HEXIST

c----------------------------------------------------------------------
c Initialize
c----------------------------------------------------------------------      

      CALL init_par_goa(Nfiles,Fname,OutFile,OutFile2)
      CALL book_ntpl_goa

      filter_ok_data = .true.
      filter_ok_mc = .true.
      Qgated_total = 0.0
      Qgated_plus  = 0.0
      doplus_BPM   = 0.0
      dominus_BPM  = 0.0
      doplus_FC    = 0.0
      dominus_BPM  = 0.0

c----------------------------------------------------------------------
c Loop over input files 
c read data from BOS file until eof reached
c close file and end run upon eof
c----------------------------------------------------------------------

      write(*,721)
  721 format(///'======================= START =======================')

      DO Ifile=1, Nfiles

c --- Some settings ---
      I_file =Ifile
      file_in   =Fname(Ifile)
      file_input=Fname(Ifile)
      N_firstfile1 = 0
      eof       = .false.
      firstTime = .true.
      Nstory    = 0
      Nstory_MC = 0
      ierr98 = 0
      ierr99 = 0
      write(*,720) Ifile,Nfiles,file_in
  720 format(/' Starting to Analyze data FILE (',I2,'of',I2,'): '/' ',A78)
      CALL open_ntpl_goa(file_in)

   
      
c --- Get the number of events in the 2 ntuples ---
      ND81=0
      NMC=0
      NDATA=0
      if(HEXIST(81))               call Hnoent(81,ND81)
      if(HEXIST(input_ntmc_num))   call Hnoent(input_ntmc_num,  NMC)
      if(HEXIST(input_ntdata_num)) call Hnoent(input_ntdata_num,NDATA)
      bit_mc=-1
      if(NDATA.gt.1)              bit_mc= 0 !DATA ntuple exists
      if(NMC  .gt.1)              bit_mc= 1 !MC ntuple exists
      if(NDATA.gt.1.and.NMC.gt.1) bit_mc= 2 !DATA & MC ntuples exist

c-------------------------------------------------------------------
c  ntuple 81
c-------------------------------------------------------------------
      IF(ND81.gt.0.and.bit_mc.eq.0)THEN
cc         write(*,*) 'TEMPORARY CHANGED NTUPLE 81 PROCESSING !!!!'
cc      IF(ND81.gt.0.and.bit_mc.eq.5)THEN
      print *,'Processing ntuple 81'
      Nstory_81 = 0
      DO WHILE (.true.)
c ----- read and write ntuple 81 ----- 
        Nstory_81 =  Nstory_81 + 1
        call HGNT(81,Nstory_81,ierr81)
        IF(ierr81.ne.0) THEN
          Nstory_81 = Nstory_81 - 1
          GOTO 2001
        ENDIF
      ENDDO
      print *,'No ntuple 81'
      ENDIF
 2001 CONTINUE

c-------------------------------------------------------------------
c  If MC and DATA ntuples are presents I will syncronize them 
c  MonteCarlo loop (ntuple 98/9) Data looop (ntuple 99/10 )                                 
c-------------------------------------------------------------------

      IF(bit_mc.eq.2) THEN
      print *,'NO PreProcessing'
      print *,'Processing MC and DATA banks in syncro'
      Nstory_MC = 0
      DO WHILE (.true.)

	Nstory_MC =  Nstory_MC + 1
        call HGNT(input_ntmc_num,Nstory_MC,ierr98)
        IF(ierr98.ne.0 .or.(Nstory_MC.gt.n_events.and.n_events.ne.0)) THEN
          Nstory_MC = Nstory_MC - 1
          GOTO 1002
        ENDIF
        call initevent_goa
        call id_choice
        call processevent_goa(Select_OK,Nstory_MC,1)
        if(.not.select_OK) goto 1021
        if(filter_flag) call filter(filter_ok_mc)
	if(filter_ok_mc.and.ana_gen.eq.1) call hfnt(61)

        nevMC = nevent !checking the MC-out event number
c Trying to get the event from GSIM-out (data)
	Nstory = Nstory + 1
        call HGNT(input_ntdata_num,Nstory,ierr99)
        if(input_type.eq.2) call clas2ceb
        if(input_type.eq.2.and.beam_type)      call our_ID_new_e1
        if(input_type.eq.2.and..not.beam_type) call our_ID_new_g
        IF(ierr99.ne.0 .or. (Nstory.gt.n_events.and.n_events.ne.0)) THEN
          Nstory_MC = Nstory_MC - 1
          GOTO 1002
        ENDIF
        if (nevent.ne.nevMC) then
          Nstory =  Nstory - 1
          goto 1011
        else
          call id_choice
          call processevent_goa(Select_OK,Nstory,0)
          if(.not.select_OK) goto 1011
          if(filter_flag) call filter(filter_ok_data)
	  if(filter_ok_data.and.ana_gen.eq.1)  CALL hfnt(60)
        endif
 1011 CONTINUE
c+ Filling specific channel ntuples
c  ChanDoneFlag( ,1) - flag GENERATI: 0- analysis_xx non e' stat eseguita; 1- ... si
c  ChanDoneFlag( ,0) - flag RECOSTRu: 0- analysis_xx non e' stat eseguita; 1- ... si
        if(indtype.ne.0.and.ana_chan.gt.0.and.ana_chan.le.8) then
        if( (filter_ok_MC  .and.ChanDoneFlag(ana_chan,1).eq.1) .or.
     &      (filter_ok_data.and.ChanDoneFlag(ana_chan,0).eq.1) ) then 
          if(ChanDoneFlag(ana_chan,1).eq.0 .and.
     &       ChanDoneFlag(ana_chan,0).eq.1 ) then
          call sas_init(1)
          endif
          if(ana_chan.eq.1) call hfnt(25)
          if(ana_chan.eq.2) call hfnt(26)
          if(ana_chan.eq.3) call hfnt(20)
          if(ana_chan.eq.5) call hfnt(30)
        endif
        endif
c+ Writing events passed
 1021   continue
        if (Nstory_MC/20000*20000.eq.Nstory_MC) then
         write(*,711) Ifile,Nfiles,Nstory_MC,NMC,'MC/DATA passed      '
  711    format(' ','  File# ',i4,' (',i4,')   ',i8,' (',i8,') ',a20)
        end if

      ENDDO  ! next event 
      ENDIF  ! If MC and DATA ntuples
 1002 CONTINUE 


c------------------------------------------------------------------- 
c  MonteCarlo loop (ntuple 98)                  
c-------------------------------------------------------------------

      IF(bit_mc.eq.1) THEN
      print *,'NO PreProcessing'
      print*, 'Processing MC bank'
      DO WHILE (.true.)
	Nstory_MC =  Nstory_MC + 1
        call HGNT(input_ntmc_num,Nstory_MC,ierr98)
        IF(ierr98.ne.0 .or.(Nstory_MC.gt.n_events.and.n_events.ne.0)) THEN
          Nstory_MC = Nstory_MC - 1
          GOTO 8000
        ENDIF
        call initevent_goa
        call id_choice               ! choose between our or seb id
        call processevent_goa(Select_OK,Nstory_MC,1)
        IF( .not. Select_OK ) GOTO 1012
        if(filter_flag) call filter(filter_ok_mc)
	if(filter_ok_mc.and.ana_gen.eq.1) CALL hfnt(61)
c+ Filling specific channel ntuples
        if(filter_ok_mc) then
        if(indtype.ne.0.and.ana_chan.gt.0.and.ana_chan.le.8) then
        if(ChanDoneFlag(ana_chan,1).eq.1) then  !1-MC
          if(ana_chan.eq.1) call hfnt(25)
          if(ana_chan.eq.2) call hfnt(26)
          if(ana_chan.eq.3) call hfnt(20)
          if(ana_chan.eq.5) call hfnt(30)
        endif
        endif
        endif
c+ Writing events passed     
1012  continue
      if (Nstory_MC/20000*20000.eq.Nstory_MC) then
        write(*,711) Ifile,Nfiles,Nstory_MC,NMC,'MC/DATA passed      '
      end if

      ENDDO  ! next event
      ENDIF
8000  continue


c------------------------------------------------------------------- 
c  Data loop (ntuple 99/10)                  
c-------------------------------------------------------------------

      IF(bit_mc.eq.0) THEN
      print *,'PreProcessing DATA bank'
      NDATAmax=NDATA
ccc        CALL shortprepfile_goa(NDATAmax)!don't delete, leave commented
ccc        GOTO 9001                       !don't delete, leave commented
      if(beam_type.AND.hel_ana.eq.1) then
         call new_pre_hel(Ifile,file_in,outfile,NDATAmax)
      endif
      if(n_events.gt.0) NDATAmax=min(NDATA,n_events)
c --- preprocessing ---
      if(beam_type)      call prepfile_goa(NDATAmax)
      if(.not.beam_type) call prepfile_goa_photon(file_in,NDATAmax)
      if(beam_type.AND.hel_ana.eq.1) then
         call match_sel(file_in)
      endif
c --- nt51 for Gamma analysis ---
      if(.not.beam_type) then 
        CALL get_51goa
        CALL Ttag_rebin
        call hfnt(51)
      endif
c --- main data loop ---
      print *,'Processing    DATA bank'
      DO WHILE (.true.)
	 Nstory =  Nstory + 1
         call HGNT(input_ntdata_num,Nstory,ierr99)
         if(input_type.eq.2) call clas2ceb
         if(input_type.eq.2.and.beam_type)      call our_ID_new_e1
         if(input_type.eq.2.and..not.beam_type) call our_ID_new_g
         IF(ierr99.ne.0.or.(Nstory.gt.n_events.and.n_events.ne.0)) THEN
           Nstory = Nstory - 1
           GOTO 9000
         ENDIF	
         call initevent_goa
         call id_choice
         call processevent_goa(Select_OK,Nstory,0)
         IF( .not. Select_OK ) GOTO 1010     
         if(qc_flag) goto 1010
         if(filter_flag) call filter(filter_ok_data)
	 if(filter_ok_data.and.ana_gen.eq.1)	 CALL hfnt(60)
c+ Filling specific channel ntuples
        if(filter_ok_data) then
        if(indtype.ne.0.and.ana_chan.gt.0.and.ana_chan.le.8) then
        if(ChanDoneFlag(ana_chan,0).eq.1) then  !1-MC
          if(ana_chan.eq.1) call hfnt(25)
          if(ana_chan.eq.2) call hfnt(26)
          if(ana_chan.eq.3) call hfnt(20)
          if(ana_chan.eq.4) call hfnt(31)
          if(ana_chan.eq.4.and.hel_ana.eq.1.and.fill311) then
            call hfnt(311)
            fill311=.false.
          endif
          if(ana_chan.eq.5) call hfnt(30)
        endif
        endif
        endif
c-       
1010 	 continue  
c+ Writing events passed     
      if (Nstory/20000*20000.eq.Nstory.or.Nstory.eq.100) then
         write(*,711) Ifile,Nfiles,Nstory,NDATA,'DATA passed         '
      end if

      ENDDO  ! next event 
      ENDIF
      
C----------------------------------------------------------
 9000   continue
	print *,'FILE completed. ',
     &          ' Events written to nt DATA/MC = ',Nstory,Nstory_MC

        print *,'PostProcessing data for the file '
        CALL postfile_goa
        print *,'Do hrend'
	CALL hrend('ESCA')
        IF(.not.qc_flag) Qgated=0.0
        WRITE(71,701) Ifile,Fname(Ifile),
     &                Ntrigger_file,Ncebevent_file,Qgated_file
 701    FORMAT(' ',i5,' ',a100,' ',2i9,e12.5)
 9001   continue

      ENDDO   ! next data file
      WRITE(71,702) Ntrigger_sum,Ncebevent_sum,Qgated_sum
 702  FORMAT(/' Total ',2i9,e12.5)
      
C----------------------------------------------------------
	
c        CALL eff_pi_minus(0,0)     ! option 0 means final efficiency calculation

 9999 CONTINUE

        DO i=1,6
        CALL hdelet(160+i)
        CALL hdelet(170+i)
        ENDDO
        IF(qc_flag) THEN
        CALL hdelet(131)
        CALL hdelet(132)
        CALL hdelet(133)
        CALL hdelet(134)
        CALL hdelet(135)
        CALL hdelet(136)
        ENDIF



        CLOSE(71)
        CALL hcdir('//DST',' ')
 	CALL hrout(0,icycle,' ')
	CALL hrend('DST')
 

	STOP
	END

