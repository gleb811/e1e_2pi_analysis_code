


c======================================================================
c  maps from the file
c======================================================================
      SUBROUTINE get_51goa
      IMPLICIT none
      INCLUDE "ntpl_goa.inc"
      cHARACTER*99    clas_map

      CHARACTER*99  namefilemap,namefilenevents
      CHARACTER*1   ch1 
      INTEGER*4     runnum,vtime,nc, MODA_NGREAD, i, i1, NN, NNN
      INTEGER*4     ierr
      REAL*4   width_ratio(61),livetimeNorm,livetimeProd
      REAL*4   G2SL(61),G4SL(61),PSnorm(61),PSprod(61)
      REAL*4   NGammaTNorm(61),  NGammaTtmp(61)
      REAL*4   NGammaENorm(767), NGammaEtmp(767)
      REAL*4   PP, add1,add2


c ----- determine run number from the file name -----
      do while(.true.)
      if(index(file_input,'/').ne.0) then
        file_input = file_input(index(file_input,'/')+1:99)
      else
        goto 3001
      endif
      enddo
 3001 continue
      i1=999
      do i=0,9
      write(ch1,'(I1)') i
      if(index(file_input,ch1).ne.0) i1=min( index(file_input,ch1), i1)
      enddo
      if(i1.eq.999) goto 1001
      file_inputnum=file_input(i1:i1+5-1)
      read(file_inputnum,*,ERR=1001) RunNumber
      goto 1002
 1001 continue
      print *,' Can not get RunNumber from the Input file Name=',file_input
      stop
 1002 continue
      print *,'  INFO: RunNumber=',RunNumber
      runnum=RunNumber 


      IF(beam_type) RETURN


c ----- get The Number Of Events for the whole run -----
      call revinm('GENOVA_PARMS','goa_nevens_in_run.dat',namefilenevents)
      open(9,FILE=namefilenevents,STATUS='OLD',ERR=2001)
      DO WHILE (.true.)      
        read(9,*,ERR=2003,END=2004) i1,NEventTotalRun
        if(i1.eq.RunNumber) goto 2002
      ENDDO
 2004 print *,' The number of Events not found for the run# ',RunNumber
      stop
 2001 print *,' Cannot open file goa_nevens_in_run.dat (from GENOVA_PARMS)'
      stop      
 2003 print *,' Error reading goa_nevens_in_run.dat file'
      stop      
 2002 continue


c ----- get First and Last Trig Number and Nevents for this file ---
      CALL Hnoent(input_ntdata_num,NNN)
      CALL HGNT(input_ntdata_num,1,  ierr)
      if(input_type.eq.2) call clas2ceb
      NEventFileFirst = event_num
      CALL HGNT(input_ntdata_num,NNN,ierr)
      if(input_type.eq.2) call clas2ceb
      NEventFileLast  = event_num
      NEventThisFile = NEventFileLast - NEventFileFirst + 1 !i1-i2
      print *,'  NEvent: TotalRun,ThisFile=',NEventTotalRun,NEventThisFile


        if(.not.beam_type) then
c+++ Getting SC ADC threshold from the MAP
     	call revinm('CLAS_PARMS','Maps/SC_CALIBRATIONS.map',clas_map)
     	call map_get_float(clas_map,'pedestals','right',
     &                 288,ped_SC_R, runnum, vtime)
     	call map_get_float(clas_map,'pedestals','left',
     &                 288,ped_SC_L, runnum, vtime)



c ----- g6 normalisation map file -----
c ----- MODA_NGREAD = 1 -> read NGamma from g6norm map. -----
c ----- MODA_NGREAD = 2 -> calculate NGamma          -----
      if    (RunNumber.gt.12000.and.RunNumber.lt.13000)then !G6A
        call revinm('CLAS_PARMS','Maps/g6normalisation.map',namefilemap)
        print *,'  INFO: g6normalisation.map is used'
      elseif(RunNumber.gt.19000.and.RunNumber.lt.20000)then !G6B
        call revinm('CLAS_PARMS','Maps/g6Bnormalisation.map',namefilemap)
        print *,'  INFO: g6Bnormalisation.map is used'
      else
        print *,' I do not know which norm. map to use. Stop.'
        stop
      endif

      MODA_NGREAD = 1
      
c ----- Read the number of photons per E/T counters ---
      IF (MODA_NGREAD.eq.1) THEN
c      call map_get_int(namefilemap,'PRODUCTION','NEvent',
c     &                 1,NEventMap, runnum, vtime )
      if(nc.ne.0) print *,' Error while g6 map reading for NEventMap'

      call map_get_float(namefilemap,'PRODUCTION','NGammaE',
     &                 767,NGammaE, runnum, vtime )
      if(nc.ne.0)  print *,' Error while map reading for NGammaE'
      
      call map_get_float(namefilemap,'PRODUCTION','NGammaE_2',
     &                 767,NGammaE_2, runnum, vtime )
      if(nc.ne.0) print *,' Error while map reading for NGammaE_2'

      call map_get_float(namefilemap,'PRODUCTION','NGammaT',
     &                 61,NGammaT, runnum, vtime )
      if(nc.ne.0) print *,' Error while map reading for NGammaT'

      call map_get_float(namefilemap,'PRODUCTION','NGammaT_2',
     &                 61,NGammaE_2, runnum, vtime )
      if(nc.ne.0) print *,' Error while map reading for NGammaT_2'
      ENDIF

c ----- CALCULATE the number of photons per E/T counters ---
      IF (MODA_NGREAD.eq.2) THEN

      call map_get_float(namefilemap,'GENERAL','width_ratio',
     &                 61,width_ratio, runnum, vtime )

      call map_get_float(namefilemap,'NORMALISATION','livetime',
     &                 1,livetimeNorm, runnum, vtime )
      call map_get_float(namefilemap,'NORMALISATION','PS',
     &                 61,PSnorm, runnum, vtime )
      call map_get_float(namefilemap,'NORMALISATION','NGammaT',
     &                 61,NGammaTNorm, runnum, vtime )
      call map_get_float(namefilemap,'NORMALISATION','NGamma',
     &                767,NGammaENorm, runnum, vtime )

      call map_get_float(namefilemap,'PRODUCTION','livetime',
     &                 1,livetimeProd, runnum, vtime )
      call map_get_float(namefilemap,'PRODUCTION','G2SL',
     &                 61,G2SL, runnum, vtime )
      call map_get_float(namefilemap,'PRODUCTION','G4SL',
     &                 61,G4SL, runnum, vtime )
      DO i=1,61
        PSprod(i) = G2SL(i) - width_ratio(i)*G4SL(i)
        PSnorm(i) = PSnorm(i)
      ENDDO
      NN=0
      PP=0.
      DO i=1,61
      IF(i.le.18) THEN
      IF(PSprod(i).ne.0. .and. PSnorm(i).ne.0.) THEN
        NN=NN+1
        PP=PP + PSprod(i)/PSnorm(i)
      ENDIF
      ENDIF
      ENDDO
      PP=PP/NN
      DO i=1,61
        NGammaT(i) = NGammaTNorm(i)*PP*(livetimeProd/livetimeNorm)
      ENDDO
      DO i=1,767
        NGammaE(i) = NGammaENorm(i)*PP*(livetimeProd/livetimeNorm)
      ENDDO
     
      ENDIF


c ----- Scale the number of photons by the num. of trig. in this file 
c ----- calculate and print the total numbers of photons 
c ----- (they should be equal for E and T counters)
c ----- the resulting NgammaTthisfile is set
      add1 = 0.
      DO i=1,61
        NGammaT(i) = NGammaT(i)*NEventThisFile/float(NEventTotalRun)
        add1=add1+NGammaT(i)
      ENDDO
      add2 = 0.
      DO i=1,767
        NGammaE(i) = NGammaE(i)*NEventThisFile/float(NEventTotalRun)
        add2=add2+NGammaE(i)
      ENDDO
      NgammaTthisfile = add1 ! add1 or add2
      print *,'  NGammaTthisfile(sumE,sumT) = ',add1,add2

c ----- save the number of photons into the histograms 11004,11002 -----
c ----- summing them up to the Ngamma in previous files            -----
      CALL HUNPAK(11004,NGammaTtmp,'HIST',0)
      DO i=1,61
      NGammaTtmp(i)=NGammaTtmp(i) + NGammaT(i)
      ENDDO
      CALL HPAK(11004,NGammaTtmp)

      CALL HUNPAK(11002,NGammaEtmp,'HIST',0)
      DO i=1,767
      NGammaEtmp(i)=NGammaEtmp(i) + NGammaE(i)
      ENDDO
      CALL HPAK(11002,NGammaEtmp)


        endif !if(.not.beam_type) then

      RETURN
      END




c======================================================================
c     wewewe
c======================================================================
      SUBROUTINE Ttag_rebin
      IMPLICIT NONE
      include "ntpl_goa.inc"
      CHARACTER*99 map_Tboundaries,map_Eboundaries
      INTEGER*4 i,i2,idum
      REAL*4    rdum,E2max,E2min,emid,add0,add1,add2
      REAL*4    emax_T,emin_T,de_T,  emax_E,emin_E,de_E
      REAL*4    emax_gev,emin_gev
      REAL*4    Emin(800),Emax(800),Namma(800)
      LOGICAL EgammaIsGood


c ----- read tagT-boundaries.dat -----      
      call revinm('CLAS_PARMS','tagT-boundaries.dat',map_Tboundaries)
      open(9,FILE=map_Tboundaries,ERR=9001)
      DO i=1,61
        read(9,*,END=9002,ERR=9002)idum,Ttag_emin(i),rdum,Ttag_emax(i)
      ENDDO
      close(9)

c ----- read tagT-boundaries.dat -----      
      call revinm('CLAS_PARMS','tagE-boundaries.dat',map_Eboundaries)
      open(9,FILE=map_Eboundaries,ERR=9003)
      DO i=1,767
      IF( mod(i+2,2).eq.1 ) THEN
        read(9,*,END=9004,ERR=9004) idum,idum,idum,idum,idum,
     &      rdum,rdum,     Etag_emin(i),Etag_emax(i)
      ELSE
        read(9,*,END=9004,ERR=9004) idum,idum,idum,idum,idum,
     &      rdum,rdum,rdum,Etag_emin(i),Etag_emax(i)
      ENDIF
      ENDDO
      close(9)
      
c --- set Bins ---
      Ttag_Nbins = 60  !the number of T counter to use
      Etag_Nbins = 700 !the number of E counter to use   
      if(Ttag_Nbins.gt.800) print *,' Ttag_rebin: .gt.800=',Ttag_Nbins
      if(Ttag_Nbins.gt.800) stop
      if(Etag_Nbins.gt.800) print *,' Etag_rebin: .gt.800=',Etag_Nbins
      if(Etag_Nbins.gt.800) stop
      
c --- Energy binning ---
      if    (RunNumber.gt.12000.and.RunNumber.lt.13000)then !G6A
        Erebin_Nbins = 4 
        Trebin_Nbins = 4 
        Emin_gev = 3.3
        Emax_gev = 3.9
        print *,'  INFO: PHOTONS ENERGY limits: ',Emin_gev,Emax_gev
      elseif(RunNumber.gt.19000.and.RunNumber.lt.20000)then !G6B
        Erebin_Nbins = 8 
        Trebin_Nbins = 8 
        Emin_gev = 2.9
        Emax_gev = 5.2
        print *,'  INFO: PHOTONS ENERGY limits: ',Emin_gev,Emax_gev
      else
        print *,' I do not know which Energy Limits to use. Stop.'
        stop
      endif

      Emax_T = Emax_gev/Eelbeam             !in % of Ebeam 
      Emin_T = Emin_gev/Eelbeam             !in % of Ebeam 
      DE_T   = (Emax_T-Emin_T)/Trebin_Nbins !in % of Ebeam 
      Emax_E = Emax_gev/Eelbeam             !in % of Ebeam 
      Emin_E = Emin_gev/Eelbeam             !in % of Ebeam 
      DE_E   = (Emax_E-Emin_E)/Erebin_Nbins !in % of Ebeam 


c --- reorder Etaggers ---
      add0=0.
      DO i=1,Etag_Nbins
        Emin(i)  =Etag_Emin(Etag_Nbins-i+1)
        Emax(i)  =Etag_Emax(Etag_Nbins-i+1)
        Namma(i) =NGammaE  (Etag_Nbins-i+1)
        add0=add0+Namma(i)
      ENDDO

c --- rebin for Etaggers ---
      DO i2=1,Erebin_Nbins
        E2min = Emin_E + DE_E*(i2-1)
        E2max = Emin_E + DE_E*(i2)
        Erebin_emin(i2)   = E2min
        Erebin_emax(i2)   = E2max
        Erebin_Ngamma(i2) = 0.0
        DO i=1,Etag_Nbins
        IF(  EgammaIsGood( Eelbeam*(Emin(i)+Emax(i))/2. )  ) THEN
        if     (Emin(i).ge.E2min.and.Emax(i).le.E2max) then
          Erebin_Ngamma(i2)=Erebin_Ngamma(i2) + 
     &                      Namma(i)
        else if(Emin(i).lt.E2min.and.Emax(i).gt.E2min) then
          Erebin_Ngamma(i2)=Erebin_Ngamma(i2) + 
     &                      Namma(i)*((Emax(i)-E2min)/(Emax(i)-Emin(i)))
        else if(Emin(i).lt.E2max.and.Emax(i).gt.E2max) then
          Erebin_Ngamma(i2)=Erebin_Ngamma(i2) + 
     &                      Namma(i)*((E2max-Emin(i))/(Emax(i)-Emin(i)))
        else if(Emin(i).le.E2min.and.Emax(i).ge.E2max) then
          print *,' Ttag_rebin: too small Energy bin'
          print *,' STOP'
          stop 
        endif
        ENDIF
        ENDDO
      ENDDO

c --- multiplay ---
      DO i=1,Erebin_Nbins
        Erebin_emax(i)=Erebin_emax(i)*Eelbeam
        Erebin_emin(i)=Erebin_emin(i)*Eelbeam
      ENDDO

c --- check total number of Gammas ---
      add1=0.
      do i=1,Etag_Nbins
      add1 = add1 + NGammaE(i)
      enddo
      add2=0.
      do i=1,Erebin_Nbins
      add2 = add2 + Erebin_Ngamma(i)
      enddo
      NgammaTtotal=add1
      if(abs((NgammaTthisfile-add0)/NgammaTthisfile).gt.0.3 .or.
     &   abs((NgammaTthisfile-add1)/NgammaTthisfile).gt.0.3  ) then
        print *,' !!! NgammaTthisfile=',NgammaTthisfile
        print *,' !!! add0,1,2==',add0,add1,add2
        print *,' !!! STOP'
        stop
      endif
      print *,'  NgammaE(1-200),NGRebinned  = ',NgammaTtotal,add2


c ----- Fill Ngamma in Energy bins histogram -----
      CALL HDELET(11006)
      CALL HBOOK1(11006,'NGamma of Energy', Erebin_Nbins,  
     &            Erebin_emin(1),Erebin_emax(Erebin_Nbins), 0.0)
      DO i=1,Erebin_Nbins
        emid = (Erebin_emax(i)+Erebin_emin(i))/2.
        CALL HFILL(11006, emid,0.0, Erebin_Ngamma(i))
      ENDDO


      RETURN
 9001 print *,' Error in opening tagT-bound=',map_Tboundaries
      stop
 9002 print *,' Error or EOF while reading tag-Tboundries.dat='
      stop
 9003 print *,' Error in opening tagE-bound=',map_Eboundaries
      stop
 9004 print *,' Error or EOF while reading tag-Eboundries.dat='
      stop
      END



c======================================================================
c     wewewe
c======================================================================
      FUNCTION EgammaIsGood(Eg)
      IMPLICIT NONE
      include "ntpl_goa.inc"
      LOGICAL  EgammaIsGood
      INTEGER*4 i,moda,ifirstzero
      REAL*4    Eg
      DATA ifirstzero/0/


      MODA = 0
      if(ifirstzero.eq.0)then
      ifirstzero=1
      print *,'  INFO: EgammaIsGood MODA =',MODA
      endif
      
      IF (MODA.eq.0) THEN ! All good
      EgammaIsGood = .true.
      RETURN
      ENDIF


      EgammaIsGood = .false.
      do i=1,Etag_Nbins  ! get the number of E scint.
      if(Eg.ge.Etag_emin(i)*Eelbeam.and.
     &   Eg.le.Etag_emax(i)*Eelbeam)then
      if(NGammaE(i).gt.0.0) then                       

         IF (MODA.eq.1) THEN ! 13off-13onoff
         if(i.ge.187) return
         if(i.ge. 17 .and. i.le. 22) return
         if(i.ge. 34 .and. i.le. 37) return
         if(i.ge. 50 .and. i.le. 52) return
         if(i.ge. 65 .and. i.le. 71) return
         if(i.ge. 96 .and. i.le. 99) return
         if(i.ge.130 .and. i.le.135) return
         if(i.ge.145 .and. i.le.150) return
         if(i.ge.175 .and. i.le.180) return
         ENDIF
         IF (MODA.eq.2) THEN ! 23off-23on-23onoff
         if(i.ge.187) return
         if(i.ge.  1 .and. i.le.  9) return
         if(i.ge. 17 .and. i.le. 22) return
         if(i.ge. 31 .and. i.le. 34) return
         if(i.ge. 49 .and. i.le. 51) return
         if(i.ge. 65 .and. i.le. 70) return
         if(i.ge. 96 .and. i.le. 99) return
         if(i.ge.129 .and. i.le.134) return
         if(i.ge.145 .and. i.le.150) return
         ENDIF

         EgammaIsGood = .true.
         RETURN
      endif
      endif
      enddo

      RETURN
      END


c======================================================================
c     wewewe
c======================================================================
      FUNCTION CorFact_multhits_23(Eg)
      IMPLICIT NONE
      REAL*4  CorFact_multhits_23,Eg,Ebin,Tbin
      Ebin = 985.38 - 1029.4*Eg/4.114
      Tbin = 0.1017*Ebin + 0.205
      CorFact_multhits_23= 1. + 0.0507 + 0.005623*Tbin
      RETURN
      END
      FUNCTION CorFact_multhits_13(Eg)
      IMPLICIT NONE
      REAL*4  CorFact_multhits_13,Eg,Ebin,Tbin
      Ebin = 985.38 - 1029.4*Eg/4.114
      Tbin = 0.1017*Ebin + 0.205
      CorFact_multhits_13= 1. + 0.1871 - 0.004488*Tbin +0.6204E-3*Tbin*Tbin
      RETURN
      END

