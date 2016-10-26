

      SUBROUTINE init_par_goa(Nfiles,Fname,OutFile,OutFile2)      
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "photon.inc"
      INTEGER*4     istat,lrecl,ierr
      INTEGER*4     Nfiles,i,ichoice
      CHARACTER*200  Fname(999), OutFile,OutFile2
      DATA lrecl/8192/                !record lenght in machine words


      ana_chan    = 0
      filter_flag = .false.
      Ntrigger_sum  = 0
      Ncebevent_sum = 0
      Qgated_sum    = 0.


c----------------------------------------------------------------------
c         Enter input information
c----------------------------------------------------------------------


c --- enter input file type ---
 11   WRITE (*,'(" Enter Input files type (1=ceb nt, 2=clas nt, 3=clas nt_filt): ",$)')
      READ  (*,*) input_type
      if    ( input_type.eq.1 ) then
        input_ntdata_num = 99
        input_ntmc_num   = 98
      elseif( input_type.eq.2 ) then
        input_ntdata_num = 10
        input_ntmc_num   =  9
      elseif( input_type.eq.3 ) then
        input_type=2 ! restoring previus input_type, used in other sub 
        input_ntdata_num = 97
        input_ntmc_num   =  9
      else
        print *,' ERROR: input files type should be 1 or 2'
        goto 1
      endif


c --- read input files names ---
 1    WRITE (*,'(" Enter NUMBER of files you wish to analyze: ",$)')
      READ *, Nfiles
      num_files=Nfiles
      WRITE (*,'(" Enter INPUT file names one by one line: ")')
      do i=1,Nfiles
        READ  (*,'(a200)') Fname(i)
        Fname(i) = Fname(i)(1:index(Fname(i),' ')-1)
      enddo
      do i=1,Nfiles
        open (unit=1,file=Fname(i),status='old',err=2)
        goto 3
  2     continue
        print *,'Input file ',Fname(i),' does NOT exist'
        goto 1
  3     continue
        close (1)  
      enddo
      write (*,'(" Input files are: ")')
      do i=1,Nfiles
        write (*,'(bna)') Fname(i)
      end do

      
c --- enter beam type ---
      WRITE (*,'(" Enter BEAM type (0 = photon, 1 = electron) ",$)')
      READ  (*,*) ichoice
      if (ichoice.eq.1) then
      beam_type= .true.      ! (electron)
      else
      beam_type = .false.    ! (photon)
      endif


c --- enter beam energy ---
      WRITE (*,'(" Enter BEAM energy (GeV) (it has sense also for gamma!) ",$)')
      READ  (*,*) E0
      Eelbeam=E0
      WRITE (*,'(" Enter torus current ",$)')
      READ  (*,*) I_torus


c --- Calibration ---
      IF(input_type.eq.2 .and. .not.beam_type)THEN
      WRITE (*,'(" Want Tagger and TOF  recalibration (0-no, 1-yes) ",$)')
      READ  (*,*) ichoice
      if (ichoice.eq.1) then
      do_recalibrate = .true.      
      else
      do_recalibrate = .false.    
      endif
      IF(.not.do_recalibrate)THEN
      WRITE (*,'(" Enter  tvsst       > ",$)')
      READ  (*,*) tvsst
      WRITE (*,'(" Enter tag_tof_corr > ",$)')
      READ  (*,*) tag_tof_corr
      ENDIF
      ELSE
      do_recalibrate = .false.
      ENDIF


c --- set Particle ID mode ---
      WRITE (*,'(" Enter particle ID mode ",$)')
      WRITE (*,'(" (1-homemade, 2-SEB) ",$)')
      READ  (*,'(i2)') ichoice
      if (ichoice.eq.1) then
        choice= .true.
      else
        choice= .false.
      endif


c --- set momentum correction mode ---
      WRITE (*,'(" Enter momentum correction mode ",$)')
      WRITE (*,'(" (0-do not correct, 1-correct) ",$)')
      READ  (*,'(i2)') ichoice
      if (ichoice.eq.1) then
        mode_momcor = .true.
      else
        mode_momcor = .false.
      end if
      IF(mode_momcor) THEN
        WRITE (*,'(" Enter momentum correction type")')
        WRITE (*,'("  1 --> e1A-B  ")')
        WRITE (*,'("  2 --> eg1  ")')
        WRITE (*,'("  3 --> eg1-Burketr ")')
        WRITE (*,'("  ?> ",$)')
        READ  (*,'(i2)') mode_momcortype
      ENDIF

c --- set Zvertix position correction mode ---
      WRITE (*,'(" Enter Zvertix correction mode  ",$)')
      WRITE (*,'(" (0-do not correct, 1-correct) ",$)')
      READ  (*,'(i2)') ichoice
      if     (ichoice.eq.1) then
        mode_zcor    = .true.
      else
        mode_zcor    = .false.
        mode_zcortype= 0
      endif
      IF(mode_zcor) THEN
       WRITE (*,'("  1 --> E1 full  target (tested: e1a)")')
       WRITE (*,'("  2 --> E1 empty target (tested: e1a)")')
       WRITE (*,'("  ?> ",$)')
       READ  (*,'(i2)') mode_zcortype
      ENDIF


c --- set quality check mode ---
      WRITE (*,'(" Enter quality check mode ",$)')
      WRITE (*,'(" (quality check mode?  0=not 1=yes) ",$)')
      READ  (*,'(i2)') ichoice
      if (ichoice.eq.1) then
        qc_flag = .true.
      else
        qc_flag = .false.
      endif


c --- set filter mode ---
      if( .not. qc_flag ) then
      WRITE (*,'(" Enter filter mode ",$)')
      WRITE (*,'(" (0-nofilter, 1-filtered) ",$)')
      READ  (*,'(i2)') ichoice
      if (ichoice.eq.1)  filter_flag= .true.


c --- set general ntuple mode ---
      WRITE (*,'(" (general ntuple ? 0=not 1=yes) ",$)')
      READ  (*,'(i2)') ana_gen


c --- set Specific channel mode ---
      WRITE (*,'(" Specific channel? 0=not 1=yes ",$)')
      READ  (*,'(i2)') ichoice
      if (ichoice.eq.1) then
      WRITE (*,'(" ( 1 p pi0 ) ",$)')
      WRITE (*,'(" ( 2 n pi+ ) ",$)')
      WRITE (*,'(" ( 3 p pi+ pi-) ",$)')
      WRITE (*,'(" ( 4 inclusive pi) ",$)')
      WRITE (*,'(" ( 5 p omega) ",$)')
      READ  (*,'(i2)') ana_chan
      endif
      endif

c---- set beam polarization analysis mode
      hel_ana=0
      if(beam_type) then
         WRITE (*,'(" Beam helicity analysis? 0=not, 1=yes ",$)')
         READ  (*,'(i2)') hel_ana
      endif

c --- read output files names ---      
      WRITE (*,'(" Enter OUTPUT file name: ",$)')
      READ  (*,'(a80)') OutFile
      OutFile2= OutFile(1:index(OutFile,' ')-1)//'_files'
      OutFile = OutFile(1:index(OutFile,' ')-1)//'.hbook'
      WRITE (*,'(" Output file name is   : ",$)')
      WRITE (*,'(bna)') OutFile
      WRITE (*,'(bna)') OutFile2


c --- set the number of events to process ---
      WRITE (*,'(" Enter number of events you want to procees (0 - all file): ",$)')
      READ  *, n_events
      if(n_events.le.0) then
        n_events=0
        PRINT *,'Allora, facciamo tutto'
      else
        PRINT *,'Allora, facciamo solo ',n_events,' events'
      endif


c----------------------------------------------------------------------
c              Initialize CERN "system" 
c----------------------------------------------------------------------

      call hlimit(128*maxpages)       !
      call hbset ('BSIZE',lrecl,ierr) !define buffersize for CWN
      !to increase max number of records in output hbook file to 64K
      IQUEST(10) = 65000

c Just to initialize minuit...
      call tmpmininit           

c Open output hbook file
      call hropen(2,'DST',OutFile,'NQ',lrecl,istat)

c Open output file containing the input files names
      OPEN(71,File=OutFile2,status='UNKNOWN')

c some vars. initialization


      RETURN
      END 




      SUBROUTINE tmpmininit
      REAL     *4 a1,a2,a3,a4,a5,a6
      call hbook1(1,' ',3,  0.0,0.0,0.0)
      call hfill (1, 1.0,  0.0,1.0)
      call hfill (1, 2.0,  0.0,1.0)
      call hfill (1, 3.0,  0.0,1.0)
      a1=0.0
      call hfithn (1, 'P0','Q', 1, a1,a2,a3,a4,a5,a6)
      call hdelet (1)
      RETURN
      END

