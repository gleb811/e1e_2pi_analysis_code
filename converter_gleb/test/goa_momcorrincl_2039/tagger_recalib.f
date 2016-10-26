      subroutine tagger_recalib(file_name,nu_ev)
      IMPLICIT NONE
      INCLUDE  "ntpl_goa.inc"
      INCLUDE  "photon.inc"

      INTEGER n_tot,nu_ev
      integer*4 ierr,Nstory
      CHARACTER*200 file_name
      LOGICAL   eof,flag_print
      real  dtime
      integer j,j_scpb,s_scpb,p_SCPB
      real c,m_pi,t0_pi
      data c,m_pi/29.9792438,0.139/

      n_tot = 0

c+ Opening file
cceebb        CALL FPARM(
cceebb     >'OPEN UNIT=12 FILE="'//file_name//'" '//
cceebb     >'ACTION=READ RECL=32760 STATUS=OLD FORM=BINARY')
c     file is already open
      ierr = 0
      Nstory = 0
      EOF = .false.
cc      Print *,' ' 
cc      Print *,'----  Starting to re-calibrate TAGGER. Run: ',file_name(3:20)


c --- READ EVENT ---
      DO WHILE (.true.)
cceebb           CALL read_event(EOF)
        Nstory =  Nstory + 1
        call HGNT(input_ntdata_num,Nstory,ierr)
        if(input_type.eq.2) call clas2ceb
        IF(ierr.ne.0) THEN
          Nstory = Nstory - 1
          EOF = .true.
          GOTO 1015
        ENDIF
        flag_print=.false.

c+++ Making TAGGER claculations
        if(.not.flag_false_g.and.n_ph_in.eq.1) then 
          if (n_tot.ge.nu_ev) goto 1015  
          n_tot = n_tot+1
          dtime = ph_TPHO_TAGR(1)-ph_STT_HEVT(1)
          call hfill(7001,dtime,0.,1.)    
          flag_print=.true.   

c++++ Making TOF pre-analysis with PI+ finding
          if(nEVNT.gt.0) then 
          do j=2,nEVNT
              if (ID_EVNT(j).eq.211.and.mass_evnt(j).lt.0.2) then
                j_SCPB = SCstat_EVNT(j)
                if (j_SCPB.eq.0) goto 1014
                s_SCPB = sector_SCPB(j_SCPB)
                p_SCPB = pd_id_SCPB(j_SCPB)
                t0_pi = path_SCPB(j_SCPB)*Sqrt(pmom_EVNT(j)**2+m_pi**2)/
c     #         (pmom_EVNT(j)*c*(0.96+0.04*beta_EVNT(j)))
     #          (pmom_EVNT(j)*c)
     #          -(Time_SCPB(j_SCPB) - ph_TPHO_TAGR(1))
                CALL HFILL(7002,t0_pi,0.,1.)
              endif
          enddo
          endif
c---  
          endif
          IF(EOF) goto 1015  
1014      continue
          if (n_tot/1000*1000.eq.n_tot.and.flag_print) then
           flag_print = .false.
           write(*,*) n_tot,Nstory,' event passed for TAGGER recalibration'
          end if  
	ENDDO  ! next event 


c-- Ending loop on events
1015    continue  
	print *,'  --- TAGGER recalib. performed using ',n_tot,' events from ',Nstory
        write(1,*) 'TAGGER recalibration performed using ',n_tot,' events'
        write(1,*) 'over ', Nstory,' events scanned'
        write(1,*) ' ' 


      RETURN
      END


