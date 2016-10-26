      subroutine tof_recalib(file_name,nu_ev)
      IMPLICIT NONE
      include  "ntpl_goa.inc"
      include  "photon.inc"

      INTEGER i,j,nu_ev,j_SCPB,s_SCPB, p_SCPB,n_tot
      integer*4 Nstory,ierr
      CHARACTER *80 file_name
      LOGICAL   eof,flag_print
      real c,m_pi,t0_pi
      data c,m_pi/29.9792438,0.139/

c+ Time constant and slope initialization (file by file)
      do i=1,7
      do j=1,48
        t0_i(i,j) = -(tvsst+tag_tof_corr)
        sl_i(i,j) = 1. 
        n_pi(i,j) = 0
      enddo
      enddo
      n_tot = 0


c+ Opening file
c        CALL FPARM(
c     >'OPEN UNIT=12 FILE="'//file_name//'" '//
c     >'ACTION=READ RECL=32760 STATUS=OLD FORM=BINARY')
c     file is already open
      ierr = 0
      Nstory = 0
      EOF = .false.
cc      Print *,' ' 
cc      Print *,'----  Starting to re-calibrate TOF. Run: ',file_name(3:20)


c --- READ EVENT ---
        DO WHILE (.true.)
c         CALL read_event(EOF)
        Nstory =  Nstory + 1
        call HGNT(input_ntdata_num,Nstory,ierr)
        if(input_type.eq.2) call clas2ceb
        IF(ierr.ne.0) THEN
          Nstory = Nstory - 1
          EOF = .true.
          GOTO 1015
        ENDIF
        call our_id_new_g           ! To avoid random coinc contamination
        flag_print=.false.

cc           print *,' 000 Nstory,flag_false_g=',Nstory,flag_false_g
cc           print *,' 000 nEVNT,n_ph_in=',nEVNT,n_ph_in
cc           print *,' 000 id_our(1)=',id_our(1)

c+++ Making t0_i calculation
          if(.not.flag_false_g.and.nEVNT.gt.0.and.n_ph_in.eq.1
     %       .and.id_our(1).eq.-1 ! Only checked good events
     %     ) then 
c++++ PI+ finding
         do j=2,nEVNT

cc           print *,' 111 ID_EVNT(j),mass_evnt(j)=',ID_EVNT(j),mass_evnt(j)
cc           print *,' 111 ph_STT_HEVT(1),ph_TPHO_TAGR(1)=',ph_STT_HEVT(1),ph_TPHO_TAGR(1),abs(ph_STT_HEVT(1)-ph_TPHO_TAGR(1))
cc           print *,' 111 pmom_EVNT(j)=',pmom_EVNT(j)


          if (ID_EVNT(j).eq.211 !
     #       .and.mass_evnt(j).lt.0.05.and.mass_evnt(j).gt.-0.2     ! pions SEB_MASS_2: 0.01 - 0.05  
     #       .and.abs(ph_STT_HEVT(1)-ph_TPHO_TAGR(1)).lt.1.0
     #       .and.pmom_EVNT(j).lt.1.5
     #  ) then ! TOF recal using pions in main RF peak +-1ns
           if (n_tot.ge.nu_ev) goto 1015
           if(beta_EVNT(j).gt.1.0) goto 1014  

c           print *,' 222 tut'

           j_SCPB = SCstat_EVNT(j)
           if (j_SCPB.eq.0) goto 1014
           n_tot = n_tot+1         
           s_SCPB = sector_SCPB(j_SCPB)
           p_SCPB = pd_id_SCPB(j_SCPB)
           if (s_SCPB.eq.0.or.p_SCPB.eq.0) goto 1014
           t0_pi = path_SCPB(j_SCPB)*Sqrt(pmom_EVNT(j)**2+m_pi**2)/
c     #             (pmom_EVNT(j)*c*(0.96+0.04*beta_EVNT(j)))
     #             (pmom_EVNT(j)*c)
     #            -(Time_SCPB(j_SCPB) - ph_TPHO_TAGR(1))
           if(abs(t0_pi).lt.5) then 
            n_pi(s_SCPB,p_SCPB) = n_pi(s_SCPB,p_SCPB) + 1
	    t0_i(s_SCPB,p_SCPB) = (t0_i(s_SCPB,p_SCPB) *
     # 			     (n_pi(s_SCPB,p_SCPB)-1)+t0_pi)/n_pi(s_SCPB,p_SCPB)
            CALL HFILL(7110,pmom_EVNT(j),beta_evnt(j),1.)
            CALL HFILL(7111,mass_EVNT(j),0.,1.)
            CALL HFILL(7100+s_SCPB,p_SCPB*1.,t0_pi,1.)
            CALL HFILL(1000*s_SCPB+p_SCPB,t0_pi,0.,1.)
            flag_print=.true.
           endif
          endif
         enddo
        endif
c---  

        IF(EOF) goto 1015  
1014    continue
         if (n_tot/5000*5000.eq.n_tot.and.flag_print) then
          flag_print = .false.
          write(*,*) n_tot,' pions passed for TOF recalibration. From ',Nstory,' events'
         end if  
	END do  ! next event 
c-- Ending loop on events
1015    continue  
	print *,'  ---  TOF recalib. performed using ',n_tot,' pions from ',Nstory,' events'
        write(1,*) 'TOF recalibration performed using ',n_tot,' pions'
        write(1,*) 'over ', Nstory,' events scanned'
        write(1,*) ' ' 



      RETURN
      END
