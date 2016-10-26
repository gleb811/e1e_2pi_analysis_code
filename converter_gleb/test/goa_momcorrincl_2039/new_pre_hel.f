c======================================================================
c Preprocess helicity:
c======================================================================
      SUBROUTINE new_pre_hel(Ifile,file_in,outfile,NMax)
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "select.inc"



      CHARACTER*(*) file_in,outfile
      INTEGER*4  Ifile
      CHARACTER*(80) out_head,out_tgbi,out_hls,out,out_key
      INTEGER*4  NMAX,ievent,ierr,i,j,ind,icycle
      INTEGER*4  ihel_tgbi
      INTEGER*4  j_hls,j_head
      INTEGER*4  helpre_head,helcur_head
      INTEGER*4  helpre_tgbi,helcur_tgbi
      INTEGER*4  hel_hls,ihel_hls
      INTEGER*4  uno,due,last,pattern
      INTEGER*4  FC,BPM,bit29,bit30,bit31,oldbit31,newbit31
      INTEGER*4  vBPM(5000),vFC(5000)
      INTEGER*4  buffer_hlsb_old(30),buffer_hlsb_new(30)
      INTEGER*4  buffer_hlsf_old(30),buffer_hlsf_new(30)
      INTEGER*4  buffer_bit30_old(30),buffer_bit30_new(30)
      INTEGER*4  nrow_hls(5000),event_hls(5000)
      INTEGER*4  tgbi(5000),hls(5000),count(5000)
      INTEGER*4  head_new(5000),nev_head(5000)
      INTEGER*4  tgbi_new(5000),count_new(5000)
      INTEGER*4  hls_new(5000),nev_hls(5000)
      INTEGER*4  lcount(2),dub_bef,dub_aft
      INTEGER*4  sync1,sync2
      INTEGER*4  offset,new_offset,off_diff,off_tot
      INTEGER*4  ncut_hls,off_hls(100),pos_hls(100)      
      INTEGER*4  ncut,off(100),pos(100)
      INTEGER*4  ngood,done,done20,prec_hls_bank,n_inv,old_pair
      INTEGER*4  npair,h1,h2,h3,h4,next_ev_num
      INTEGER*4  start_point,end_point
      INTEGER*4  numrun
      REAL*4     run_number,old_run
      REAL*4     mean_dose
      INTEGER*4  i_head
      INTEGER*4  counter
      LOGICAL endread,match,flag20,beam,same_hls_bank
      LOGICAL forward
      save run_number
      save lcount

c------------------------------------------------------
c No preprocessing if MC analysis or DATA&MC together
c------------------------------------------------------
      IF(bit_mc.ne.0.or.qc_flag) THEN
        GOTO 1002
      ENDIF



c------------------------------------------------------
c Output file names
c------------------------------------------------------
      file_in=file_in(1:index(file_in,' ')-1)
      out_head=file_in(1:index(file_in,'h')-2)//'.head'
      out_tgbi=file_in(1:index(file_in,'h')-2)//'.tgbi'
      out_hls =file_in(1:index(file_in,'h')-2)//'.hls'
      out     =file_in(1:index(file_in,'h')-2)//'.hel'
      out_key =file_in(1:index(file_in,'h')-2)//'.key'

c------------------------------------------------------
c Look for run number: this is needed to check if the 
c present file belongs to the same run of the previous.
c In that case the tgbi helicity scaler information is 
c connected to the previous file
c------------------------------------------------------
      old_run=run_number
c      ind=index(file_in,'_')
c      run_number=file_in(ind+2:ind+7)

      run_number=RunNumber

      print *,'run number is:',run_number


      
c------------------------------------------------------
c Open files that will contain helicity info
c------------------------------------------------------
      OPEN(UNIT=15,FILE=out_head,FORM='formatted')
      OPEN(UNIT=25,FILE=out_tgbi,FORM='formatted')
      OPEN(UNIT=35,FILE=out_hls ,FORM='formatted')
      OPEN(UNIT=45,FILE=out     ,FORM='formatted')
      OPEN(UNIT=55,FILE=out_key ,FORM='formatted')
      


      print *,'reading helicity info'
     

c---- select matching algorithm as a function of run number

      if(run_number.ge.14992.AND.run_number.lt.15000) then
         goto 10001
      elseif(run_number.ge.15070) then
         goto 20001
      endif

c=============================================================
c=============================================================
c=============================================================
c First matching algorithm
c=============================================================
c=============================================================
c=============================================================
         
c--------------------------------------------------
c Read the whole DATA ntupl 
c--------------------------------------------------
      endread = .false.
      ievent=0
      ihel_head =0
      ihel_tgbi=0
      ihel_hls=0
      FC=0
      BPM=0
      uno=0
      due=0
      pattern=0
      bit31=0
      bit30=0
      bit29=0
      helcur_tgbi=-1
      helcur_head=-1      
      DO j=1,30
         buffer_hlsb_new(j)=0
         buffer_hlsf_new(j)=0
         buffer_bit30_new(j)=0
      ENDDO
      IF(bit_mc.eq.0) THEN
      DO WHILE (.true.)

c------------------------------------------------------
c       read one event 
c------------------------------------------------------
        ievent=ievent+1
        CALL HGNT(99,ievent,ierr)
        IF(ierr.ne.0 .or. ievent.gt.Nmax) THEN
          endread = .true.
          ievent=ievent-1
          nevent_read = ievent
        ENDIF




c------------------------------------------------------
c Read Helicities and save information in arrays:
c
c head        helicity and sync read from head bank
c tgbi        helicity and sync read from tgbi bank
c count       helicity scaler read from tgbi bank
c event_head  event number in the ceb ntuple
c hls         helicity and sync read from hls bank
c------------------------------------------------------

      IF( .not. endread) THEN
        

        helpre_head = helcur_head
        helcur_head = helicity
        IF(helcur_head.ne.helpre_head) THEN
           ihel_head=ihel_head+1
           write(15,*)  ievent, helcur_head
           head(ihel_head)=helcur_head
           event_head(ihel_head)=ievent
        ENDIF

        helpre_tgbi=helcur_tgbi
        helcur_tgbi=0
        call mvbits(helTGBI(1),14,2,helcur_tgbi,0)
        IF(helcur_tgbi.ne.helpre_tgbi) THEN
           ihel_tgbi=ihel_tgbi+1
           write(25,*)ievent, helcur_tgbi,heltgbi(2)
           tgbi(ihel_tgbi)=helcur_tgbi
           count(ihel_tgbi)=heltgbi(2)
        ENDIF

        IF(nHLS.ne.0) THEN
           DO j=1,30
              buffer_hlsb_old(j)=buffer_hlsb_new(j)
              buffer_hlsb_new(j)=0
              buffer_hlsf_old(j)=buffer_hlsf_new(j)
              buffer_hlsf_new(j)=0
              buffer_bit30_old(j)=buffer_bit30_new(j)
              buffer_bit30_new(j)=0
           ENDDO
           DO j=1,NHLS
              call mvbits(c_bpm(j),0,24,buffer_hlsb_new(j),0)
              call mvbits(c_fc (j),0,24,buffer_hlsf_new(j),0)
              call mvbits(c_fc(j),30,1,buffer_bit30_new(j),0)
           ENDDO
           same_hls_bank=.true.
           DO j=1,NHLS
              same_hls_bank=same_hls_bank.AND.
     &                 (buffer_hlsb_new(j).eq.buffer_hlsb_old(j)).AND.
     &                 (buffer_hlsf_new(j).eq.buffer_hlsf_old(j)).AND. 
     &                 (buffer_bit30_new(j).eq.buffer_bit30_old(j))
           ENDDO
           if(.not.same_hls_bank) then
              do i=1,4
                 call mvbits(c_fc(i),30,1,uno,0)
                 call mvbits(c_fc(i+1),30,1,due,0)
                 if(uno.ne.due) call mvbits(c_fc(i+1),29,1,pattern,0)
              enddo
              DO j=1,NHLS
                 ihel_hls=ihel_hls+1
                 call mvbits(c_fc(j),0,24,FC,0)
                 call mvbits(c_bpm(j),0,24,BPM,0)
                 call mvbits(c_bpm(j),29,1,bit29,0)
                 call mvbits(c_bpm(j),30,1,bit30,0)
                 call mvbits(c_bpm(j),31,1,bit31,0)
                 
                 write(35,*)  ievent,bit31,bit30,bit 29,bpm,fc
                 if(pattern.eq.1) then
                    hel_hls=2*bit30+1-bit29
c     print *,'inverted sync bit'
                 else
                    hel_hls=2*bit30+bit29
                 endif
                 hls(ihel_hls)=hel_hls
                 nrow_hls(ihel_hls)=NHLS
                 event_hls(ihel_hls)=ievent
                 vBPM(ihel_hls)=BPM
                 vFC (ihel_hls)=FC
              ENDDO
           endif
        ENDIF
      ENDIF

      IF( endread ) GOTO 1000
      ENDDO !DO WHILE (.true)
      ENDIF 
 1000 CONTINUE

c--------------------------------------------------------------
c  Close files that contain helicity info
c--------------------------------------------------------------
      CLOSE(15)
      CLOSE(25)
      CLOSE(35)


c--------------------------------------------------------------
c Start helicity matching
c--------------------------------------------------------------
c ihel_HEAD                     dimension of head and tgbi vectors
c ihel_hls                      dimension of hls vector
      

c-----------------------------------------------------------------
c Define the helicity scaler offset:use recorded info if this file 
c and the previous one belong to same run 
c----------------------------------------------------------------- 
     
      IF(run_number.ne.old_run) THEN
         print *,'first file of a new run...'
         lcount(1)=0
         lcount(2)=0
      ENDIF
      offset=(count(1)-lcount(1))*2
      IF(count(1).ne.count(2)) THEN
         offset=offset-lcount(2)
      ELSE
         offset=offset-1-lcount(2)
      ENDIF 
      print *,'offset=',offset
      print *,count(1),count(2)


c------------------------------------------------------------------
c Initialize final vectors
c------------------------------------------------------------------
      DO j=1,5000
         head_new(j)=-1
         tgbi_new(j)=-1
         count_new(j)=-1
         nev_head(j)=-1
         dose_BPM(j)=0
         dose_FC (j)=0
         key(j)=0
         key_pair(j)=0
         key_both(j)=0
      ENDDO

c--- the final vectors will contain the helicity information after 
c--- the matching. That will be the helicity if match was found or 
c--- -1 otherwise. The no match problem can occur when that helicity 
c--- state was missing in the tgbi and head bank or when no agreement 
c--- with the hls was obtained.



c------------------------------------------------------------------
c  Initialize variables to start matching
c------------------------------------------------------------------
c pos   vector containing the position of the missing scaler counts
c off   vector containing the number of missing count at the 
c       corresponding position


      i=1
      j=1
      pos(1)=0
      off(1)=offset
      
      print *,'start matching helicity info'
c----------------------------------------------------------------
c Look for missing count first
c----------------------------------------------------------------
      do j=1,ihel_head-1
         if(count(j+1)-count(j).gt.1) then
            i=i+1
            pos(i)=j !----- last index before missing count
            if(mod(head(j),2).eq.0) then
               dub_bef=1
            else
               dub_bef=0
            endif               
            if(mod(head(j+1),2).eq.0) then
               dub_aft=1
            else
               dub_aft=0
            endif               
            offset=(count(j+1)-count(j))*2
            offset=offset-2+(1-dub_bef)+dub_aft
            off(i)=offset
         elseif(count(j+1)-count(j).eq.1) then
            sync1=0
            sync2=0
            call mvbits(head(j),0,1,sync1,0)
            call mvbits(head(j+1),0,1,sync2,0)
            if(sync2-sync1.ne.1) then
               i=i+1
               pos(i)=j
               off(i)=sync1-sync2+1 
            endif
         endif
      enddo
      ncut=i
      pos(ncut+1)=ihel_head

c------------------------------------------------------------------
c Start matching.....
c------------------------------------------------------------------
      off_diff=0
      new_offset=0
      prec_hls_bank=20
      beam=.true.
      do i=1,ncut
         done=0
         done20=0
         if(beam) then
            new_offset=new_offset+mod(off(i),20)
         elseif(off(i).gt.40-prec_hls_bank) then
            new_offset=new_offset+mod(off(i)-(40-prec_hls_bank),20)
         else
            new_offset=new_offset+off(i) 
         endif
         if(off(i).ge.20) then
            flag20=.true.
         else
            flag20=.false.
         endif
  101    j=pos(i)+1
         if(pos(i)+1+new_offset+off_diff.gt.ihel_hls) goto 202 
         off_tot=new_offset+off_diff
         match=.true.
c         if(vFC(j+off_tot).eq.0.AND.(.not.beam)) then
c            match=.false.
c            beam=.false.
c            goto 102
c         endif
c         beam =.true.
         do while(j.le.pos(i+1).AND.j+off_tot.le.ihel_hls)
            if(j+off_tot.gt.1) then
               if(event_hls(j+off_tot).ne.event_hls(j+off_tot-1)) then
                  prec_hls_bank=nrow_hls(j+off_tot-1)
               endif
               if(vFC(j+off_tot).eq.0.AND.
     &            vFC(min(pos(i+1)+off_tot,ihel_hls)).eq.0) then
                  beam=.false.
                  print *,'no beam'
                  goto 102
               endif
            endif
            match=match.AND.head(j).eq.tgbi(j)
     &                 .AND.tgbi(j).eq.hls(j+off_tot)
            if(match.AND.j.gt.pos(i)+1.AND.j.lt.pos(i+1)) then
               key(j)=1
               dose_BPM(j)=vBPM(j+off_tot)
               dose_FC (j)=vFC (j+off_tot)
               head_new(j+off_tot) =head(j)
               tgbi_new(j+off_tot) =tgbi(j)
               count_new(j+off_tot)=count(j)
               nev_head(j+off_tot)=event_head(j)
            endif
c            print *,j,j+off_tot,head_new(j+off_tot),tgbi_new(j+off_tot),
c     &              hls(j+off_tot)     
            j=j+1 
         enddo


c----- If mismatch is found try to change offset to recover
  102    if(match) then
            print *,'match found'
            beam=.true.
         elseif(.not.match.AND.done.lt.5) then
            j=pos(i)+1
            do while(j.le.pos(i+1).AND.j+off(i).le.ihel_hls)
               key(j)=0
               dose_BPM(j)=0
               dose_FC (j)=0
               head_new(j+off_tot) =-1
               tgbi_new(j+off_tot) =-1
               count_new(j+off_tot)=-1
               nev_head(j+off_tot) =-1
               j=j+1
            enddo
            done=done+1
            off_diff=off_diff-1
         print *,'no match found: retrying with offset changed of',done
         goto 101
         elseif(.not.match.AND.done.eq.5.AND.flag20) then
            done=0
            off_diff=off_diff+5+20
            done20=done20+1
            if(off_diff.gt.off(i)-mod(off(i),20)) then
               done=6
               off_diff=off_diff-done20*20
               flag20=.false.
            endif               
         print *, 'missing count # greater than 20; trying offset+20'
         goto 101
         elseif(.not.match.AND.done.eq.5.) then
            done=done+1
            off_diff=off_diff+5	
         print *,'not match found after 4 attempt'
         print *,'gone back to initial offset and proceeding forward'
         goto 101       
         endif
      enddo

c---- After 4 attempt whereas succeded or not go to following helicity group

c---- Concludes writing matched helicity in text file
  202 continue
      DO i=1,ihel_hls
         write(45,*) nev_head(i),head_new(i),tgbi_new(i),hls(i),
     &               count_new (i),event_hls(i),nrow_hls(i),
     &               vbpm(i),vfc(i)
      ENDDO
      CLOSE(45)

c-------------------------------------------------------------------
c find pairs of helicity states
c-------------------------------------------------------------------

      n_inv=0
      ind=0
      npair=0

c --- find a cut ----
      j=0


 2001 CONTINUE
      j=j+1
      start_point=pos(j)
      end_point  =pos(j+1)
      DO i=start_point+1,end_point-3 
         h1=head(i+0)
         h2=head(i+1)
         h3=head(i+2)
         h4=head(i+3)
         IF( (h1.eq.0 .and. h2.eq.1 .and. h3.eq.0 .and. h4.eq.1) .OR.
     &        (h1.eq.2 .and. h2.eq.3 .and. h3.eq.2 .and. h4.eq.3) ) THEN
            n_inv=i+2
            GOTO 3001 
         ENDIF
      ENDDO
      GOTO 5001                 ! no cut has been found
c --- if the cut was found then go back from ncut to nbufstart ---
c --- taking good events until bad one occures                 ---
 3001 CONTINUE
      i=n_inv-1+4
      DO WHILE(i-4-3 .gt. start_point)
         i=i-4
         h1=head(i-3)
         h2=head(i-2)
         h3=head(i-1)
         h4=head(i-0)
         IF(h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3 .OR.
     &        h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)THEN
            npair=npair+1
            key_pair(i-3)=npair
            key_pair(i-2)=npair
            key_pair(i-1)=npair
            key_pair(i-0)=npair
         ELSE
            GOTO 4001 
         ENDIF
      ENDDO
c --- after the cut was found then go forth from ncut ---
c --- taking good events until bad one occures        ---
 4001 CONTINUE
      i=n_inv-4
      DO WHILE( i+4+3.le.end_point)
        i=i+4 
        h1=head(i+0)
        h2=head(i+1)
        h3=head(i+2)
        h4=head(i+3)
        IF((h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3) .OR.
     &     (h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)  )THEN
           npair=npair+1
           key_pair(i+0)=npair
           key_pair(i+1)=npair
           key_pair(i+2)=npair
           key_pair(i+3)=npair           
        ELSE
          GOTO 2001 
        ENDIF
      ENDDO
 5001 CONTINUE
      IF(j.lt.ncut) GOTO 2001
c-----------------------------------------------------------------------


 6001 ngood=0
      DO i=1,ihel_head
         if(i.eq.1) then
            old_pair=0
         else
            old_pair=key_pair(i-1)
         endif
         if(i.eq.ihel_head) then
            next_ev_num=nevent_read
         else
            next_ev_num=event_head(i+1)
         endif
         st_mult(i)=next_ev_num-event_head(i)
         if(key_pair(i).ne.old_pair.AND.key_pair(i).ne.0) then
            key_both(i+0)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+1)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+2)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+3)=key(i)*key(i+1)*key(i+2)*key(i+3)
         endif
         if(key_both(i).eq.1) then
            ngood=ngood+1
         else
            dose_fc(i)=0
            dose_BPM(i)=0
         endif
         write(55,*) event_head(i), st_mult(i)
     &              , head(i), key(i), key_pair(i+0)
     &              , key_both(i), dose_fc(i),dose_BPM(i)

      ENDDO
      print *,'helicity matching succeded at',ngood*100/ihel_head,'%'
      CLOSE(55)


c---  information about the helicity scaler for the last matched event
      last=ihel_hls
      lcount(1)=count_new(last)
      IF(mod(head_new(last),2).eq.0) THEN
         lcount(2)=1
      ELSE
         lcount(2)=0
      ENDIF
      IF(count_new(ihel_hls).eq.-1) THEN
         DO WHILE(count_new(last).eq.-1)
            last=last-1
         ENDDO
         IF(mod(head_new(last),2).eq.0) THEN
            lcount(1)=count_new(last)+int((ihel_hls-last+1)/2)
            lcount(2)=1-mod(ihel_hls-last,2)
         ELSE
            lcount(1)=count_new(last)+int((ihel_hls-last)/2)
            lcount(2)=mod(ihel_hls-last,2)
         ENDIF
      ENDIF
      print *,'last count is',lcount(1),' with mode',lcount(2)

      print *, 'helicity analysis concluded'

      goto 1002 


c===================================================================
c===================================================================
c===================================================================
c Second matching algorithm
c===================================================================
c===================================================================
c===================================================================
10001 continue



         
c--------------------------------------------------
c Read the whole DATA ntupl 
c--------------------------------------------------
      endread = .false.
      ievent=0
      ihel_head =0
      ihel_tgbi=0
      ihel_hls=0
      FC=0
      BPM=0
      uno=0
      due=0
      pattern=0
      bit31=0
      bit30=0
      bit29=0
      oldbit31=0
      newbit31=0
      helcur_tgbi=-1
      helcur_head=-1  
      DO j=1,30
         buffer_hlsb_new(j)=0
         buffer_hlsf_new(j)=0
         buffer_bit30_new(j)=0
      ENDDO
      IF(bit_mc.eq.0) THEN
      DO WHILE (.true.)

c------------------------------------------------------
c       read one event 
c------------------------------------------------------
        ievent=ievent+1
        CALL HGNT(99,ievent,ierr)
        IF(ierr.ne.0 .or. ievent.gt.Nmax) THEN
          endread = .true.
          ievent=ievent-1
          nevent_read = ievent
        ENDIF




c------------------------------------------------------
c Read Helicities and save information in arrays:
c
c head        helicity and sync read from head bank
c tgbi        helicity and sync read from tgbi bank
c count       helicity scaler read from tgbi bank
c event_head  event number in the ceb ntuple
c hls         helicity and sync read from hls bank
c------------------------------------------------------

      IF( .not. endread) THEN
        

        helpre_head = helcur_head
        helcur_head = helicity
        IF(helcur_head.ne.helpre_head) THEN
           ihel_head=ihel_head+1
           write(15,*)  ievent, helcur_head
           head(ihel_head)=helcur_head
           event_head(ihel_head)=ievent
        ENDIF

        helpre_tgbi=helcur_tgbi
        helcur_tgbi=0
        call mvbits(helTGBI(1),14,2,helcur_tgbi,0)
        IF(helcur_tgbi.ne.helpre_tgbi) THEN
           ihel_tgbi=ihel_tgbi+1
           write(25,*)ievent, helcur_tgbi,heltgbi(2)
           tgbi(ihel_tgbi)=helcur_tgbi
           count(ihel_tgbi)=heltgbi(2)
        ENDIF

        IF(nHLS.ne.0) THEN
           DO j=1,30
              buffer_hlsb_old(j)=buffer_hlsb_new(j)
              buffer_hlsb_new(j)=0
              buffer_hlsf_old(j)=buffer_hlsf_new(j)
              buffer_hlsf_new(j)=0
              buffer_bit30_old(j)=buffer_bit30_new(j)
              buffer_bit30_new(j)=0
           ENDDO
           DO j=1,NHLS
              call mvbits(c_bpm(j),0,24,buffer_hlsb_new(j),0)
              call mvbits(c_fc (j),0,24,buffer_hlsf_new(j),0)
              call mvbits(c_fc(j),30,1,buffer_bit30_new(j),0)
           ENDDO
           same_hls_bank=.true.
           DO j=1,NHLS
              same_hls_bank=same_hls_bank.AND.
     &             (buffer_hlsb_new(j).eq.buffer_hlsb_old(j)).AND.
     &             (buffer_hlsf_new(j).eq.buffer_hlsf_old(j)).AND.
     &             (buffer_bit30_new(j).eq.buffer_bit30_old(j))
           ENDDO
           if(.not.same_hls_bank) then
              DO j=1,NHLS
                 ihel_hls=ihel_hls+1
                 oldbit31=bit31
                 call mvbits(c_fc(j),0,24,FC,0)
                 call mvbits(c_bpm(j),0,24,BPM,0)
                 call mvbits(c_bpm(j),29,1,bit29,0)
                 call mvbits(c_bpm(j),30,1,bit30,0)
                 call mvbits(c_bpm(j),31,1,bit31,0) 
             		
                 write(35,*)  ievent,bit31,bit30,bit29,bpm,fc
c----- check to find wrong bit 31              
                 if(j.lt.NHLS.AND.j.gt.1) then
                    call mvbits(c_bpm(j+1),31,1,newbit31,0)
                    if(oldbit31.eq.bit31.AND.bit31.eq.newbit31) then
                       bit31=1-bit31
                    endif
                 endif
c----------------------------------
                 hel_hls=2*bit30+bit31               
                 hls(ihel_hls)=hel_hls
                 nrow_hls(ihel_hls)=NHLS
                 event_hls(ihel_hls)=ievent
                 vBPM(ihel_hls)=BPM
                 vFC (ihel_hls)=FC
              ENDDO
           endif
        ENDIF
      ENDIF

      IF( endread ) GOTO 11000
      ENDDO !DO WHILE (.true)
      ENDIF 
11000 CONTINUE

c--------------------------------------------------------------
c  Close files that contain helicity info
c--------------------------------------------------------------
      CLOSE(15)
      CLOSE(25)
      CLOSE(35)


      
c-----------------------------------------------------------------
c Define the helicity scaler offset:use recorded info if this file 
c and the previous one belong to same run 
c----------------------------------------------------------------- 
     
      IF(run_number.ne.old_run) THEN
         print *,'first file of a new run...'
         lcount(1)=0
         lcount(2)=0
      ENDIF
      offset=(count(1)-lcount(1))*2
      IF(count(1).ne.count(2)) THEN
         offset=offset-lcount(2)
      ELSE
         offset=offset-1-lcount(2)
      ENDIF 

      print *,'offset=',offset
      print *,count(1),count(2)


c------------------------------------------------------------------
c Initialize final vectors
c------------------------------------------------------------------
      DO j=1,5000
         head_new(j) =-1
         tgbi_new(j) =-1
         count_new(j)=-1
         nev_head(j) =-1
         hls_new(j)  =-1
         nev_hls(j)  =-1
         dose_BPM(j) = 0
         dose_FC (j) = 0
         key(j)      = 0
         key_pair(j) = 0
         key_both(j) = 0
      ENDDO

c--- the final vectors will contain the helicity information after 
c--- the matching. That will be the helicity if match was found or 
c--- -1 otherwise. 



c----- calculate mean dose (dose uncertanty is assumed to be one count)
      mean_dose=0
      DO i=1,ihel_hls
         mean_dose=mean_dose+vFC(i)
      ENDDO
      mean_dose=mean_dose/ihel_hls 
      mean_dose=mean_dose*0.95

c---- find "integrated" hls events
      
      pos_hls(1)=1
      off_hls(1)=0
      i=1
      DO j=1,ihel_hls
         if(vFC(j).gt.mean_dose*1.5) then
            print *,'integrated hls event found'
            i=i+1
            pos_hls(i)=j !----- index of integrated event
            off_hls(i)=nint(vFC(j)/mean_dose)-1
         endif
      ENDDO
      ncut_hls=i
      pos_hls(ncut_hls+1)=ihel_hls+1

c---- look for missing counts
      pos(1)=1
      off(1)=offset
      i=1
      do j=1,ihel_head-1
         if(count(j+1)-count(j).gt.1) then
            i=i+1
            pos(i)=j+1 !----- first index after missing count
            if(mod(head(j),2).eq.0) then
               dub_bef=1
            else
               dub_bef=0
            endif               
            if(mod(head(j+1),2).eq.0) then
               dub_aft=1
            else
               dub_aft=0
            endif               
            offset=(count(j+1)-count(j))*2
            offset=offset-2+(1-dub_bef)+dub_aft
            off(i)=offset
         endif
      enddo
      ncut=i
      pos(ncut+1)=ihel_head+1



c---- if beginning of run....
      if(count(1).lt.20) then
         off_hls(1)=1
         off(1)    =0
      endif


      print*,'start matching helicity info'

c------------------------------------------------------------------
c Start matching.....
c------------------------------------------------------------------
      off_diff=0
      off_tot=0


      i_head=2
      i=0
      DO WHILE(i.lt.ncut_hls)
         i=i+1
         done=0

c---- starting....
10101    j=pos_hls(i)

         if(j.eq.1.AND.off(1).eq.0) then
            j=1
            off_tot=0
         elseif(j.eq.1.AND.off(1).ne.0) then
            DO WHILE(pos_hls(i+1).lt.off(1)+1)
               i=i+1
            ENDDO
            pos_hls(i)=off(1)+1            
            j=off(1)+1
            off_tot=-off(1)
         elseif(j+off_tot+1.eq.pos(i_head)) then
            off_tot=off_tot-off(i_head)
            i_head=i_head+1 
         endif

         off_tot=off_tot+off_hls(i)
         match=.true.
         do while(j.ge.pos_hls(i).AND.j.le.pos_hls(i+1)-1)
            if(j.ne.pos_hls(i)) then
               match=match.AND.(head(j+off_tot).eq.tgbi(j+off_tot))
     &              .AND.(tgbi(j+off_tot).eq.hls(j))
            endif
            if(match.AND.j.ne.pos_hls(i)) then
               key(j+off_tot)=1
               dose_BPM(j+off_tot)=vBPM(j)
               dose_FC (j+off_tot)=vFC (j)
               head_new(j+off_tot) =head(j+off_tot)
               tgbi_new(j+off_tot) =tgbi(j+off_tot)
               count_new(j+off_tot)=count(j+off_tot)
               nev_head(j+off_tot)=event_head(j+off_tot)
               hls_new(j+off_tot)=hls(j)
               nev_hls(j+off_tot)=event_hls(j)
            endif
c            print *,j,j+off_tot,head_new(j+off_tot),tgbi_new(j+off_tot),
c     &              hls(j),vBPM(j),vFC (j),count_new(j+off_tot),
c     &              nev_head(j+off_tot)

            j=j+1 
         enddo    
         
            
c----- If mismatch is found try to change offset to recover
10102    if(match) then
            print *,'match found'
         
         elseif(.not.match.AND.done.lt.3) then
            j=pos_hls(i)+1            
            if(j+off_tot+1.eq.pos(i_head-1)) then
               off_tot=off_tot+off(i_head-1)
               i_head=i_head-1
            endif 
            do while(j.le.pos_hls(i+1)-1)              
               key(j+off_tot)=0
               dose_BPM(j+off_tot)=0
               dose_FC (j+off_tot)=0
               head_new(j+off_tot) =-1
               tgbi_new(j+off_tot) =-1
               count_new(j+off_tot)=-1
               nev_head(j+off_tot)=-1
               nev_hls(j+off_tot)=-1
               j=j+1
            enddo
            done=done+1
            off_tot=off_tot-off_hls(i)+1
            off_diff=-1
         print *,'no match found: retrying with offset changed of',done
         goto 10101

         elseif(.not.match.AND.done.eq.3) then
            done=done+1
            off_tot=off_tot-4
         print *,'not match found after 4 attempt'
         print *,'gone back to initial offset and proceeding forward'
         goto 10101      
         endif
      enddo

c---- After 4 attempt whereas succeded or not go to following helicity group

c---- Concludes writing matched helicity in text file
10202 continue
      DO i=1,ihel_head
         write(45,*) nev_head(i),head_new(i),tgbi_new(i),hls_new(i),
     &               count_new (i),nev_hls(i),
     &               dose_BPM(i),dose_FC(i)
      ENDDO
      CLOSE(45)


c-------------------------------------------------------------------
c find pairs of helicity states
c-------------------------------------------------------------------

      n_inv=0
      ind=0
      npair=0

c --- find a cut ----
      j=1


12001 CONTINUE
      DO WHILE(key(j).ne.1.AND.j.le.ihel_head)
         j=j+1
      ENDDO
      start_point=j
      if(start_point.eq.ihel_head) goto 16001
      DO WHILE(key(j).eq.1.AND.j.lt.ihel_head)
         j=j+1
      ENDDO
      if(key(j).eq.1) then
         end_point=j-1
      else
         end_point=j
      endif
      DO i=start_point,end_point-3 
         h1=head(i+0)
         h2=head(i+1)
         h3=head(i+2)
         h4=head(i+3)
         IF( (h1.eq.0 .and. h2.eq.1 .and. h3.eq.0 .and. h4.eq.1) .OR.
     &       (h1.eq.2 .and. h2.eq.3 .and. h3.eq.2 .and. h4.eq.3) ) THEN
            n_inv=i+2
            GOTO 13001 
         ENDIF
      ENDDO
      GOTO 15001                 ! no cut has been found
c --- if the cut was found then go back from ncut to nbufstart ---
c --- taking good events until bad one occures                 ---
13001 CONTINUE
      i=n_inv-1+4
      DO WHILE(i-4-3 .ge. start_point)
         i=i-4
         h1=head(i-3)
         h2=head(i-2)
         h3=head(i-1)
         h4=head(i-0)
         IF(h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3 .OR.
     &        h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)THEN
            npair=npair+1
            key_pair(i-3)=npair
            key_pair(i-2)=npair
            key_pair(i-1)=npair
            key_pair(i-0)=npair
         ELSE
            GOTO 14001 
         ENDIF
      ENDDO
c --- after the cut was found then go forth from ncut ---
c --- taking good events until bad one occures        ---
14001 CONTINUE
      i=n_inv-4
      DO WHILE( i+4+3.le.end_point)
        i=i+4 
        h1=head(i+0)
        h2=head(i+1)
        h3=head(i+2)
        h4=head(i+3)
        IF((h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3) .OR.
     &     (h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)  )THEN
           npair=npair+1
           key_pair(i+0)=npair
           key_pair(i+1)=npair
           key_pair(i+2)=npair
           key_pair(i+3)=npair           
        ELSE
          GOTO 12001 
        ENDIF
      ENDDO
15001 CONTINUE
      IF(end_point.lt.ihel_head) GOTO 12001
c-----------------------------------------------------------------------


16001 ngood=0
      DO i=1,ihel_head
         if(i.eq.1) then
            old_pair=0
         else
            old_pair=key_pair(i-1)
         endif
         if(i.eq.ihel_head) then
            next_ev_num=nevent_read
         else
            next_ev_num=event_head(i+1)
         endif
         st_mult(i)=next_ev_num-event_head(i)
         if(key_pair(i).ne.old_pair.AND.key_pair(i).ne.0) then
            key_both(i+0)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+1)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+2)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+3)=key(i)*key(i+1)*key(i+2)*key(i+3)
         endif
         if(key_both(i).eq.1) then
            ngood=ngood+1
         else
            dose_fc(i)=0
            dose_BPM(i)=0
         endif
         write(55,*) event_head(i), st_mult(i)
     &              , head(i), key(i), key_pair(i+0)
     &              , key_both(i), dose_fc(i),dose_BPM(i)

      ENDDO
      print *,'helicity matching succeded at',ngood*100/ihel_head,'%'
      CLOSE(55)

 


c---  information about the helicity scaler for the last matched event
      last=ihel_hls+off_tot
      lcount(1)=count(last)
      IF(mod(head(last),2).eq.0) THEN
         lcount(2)=1
      ELSE
         lcount(2)=0
      ENDIF

      print *,'last count is',lcount(1),' with mode',lcount(2)

      print *, 'helicity analysis concluded'

      goto 1002 

c=============================================================
c=============================================================
c=============================================================
c Third matching algorithm
c=============================================================
c=============================================================
c=============================================================

20001  CONTINUE

c--------------------------------------------------
c Read the whole DATA ntupl 
c--------------------------------------------------
      endread = .false.
      ievent=0
      ihel_head =0
      ihel_tgbi=0
      ihel_hls=0
      FC=0
      BPM=0
      uno=0
      due=0
      pattern=0
      bit31=0
      bit30=0
      bit29=0
      helcur_tgbi=-1
      helcur_head=-1      
      DO j=1,30
         buffer_hlsb_new(j)=0
         buffer_hlsf_new(j)=0
         buffer_bit30_new(j)=0
      ENDDO
      IF(bit_mc.eq.0) THEN
      DO WHILE (.true.)

c------------------------------------------------------
c       read one event 
c------------------------------------------------------
        ievent=ievent+1
        CALL HGNT(99,ievent,ierr)
        IF(ierr.ne.0 .or. ievent.gt.Nmax) THEN
          endread = .true.
          ievent=ievent-1
          nevent_read = ievent
        ENDIF




c------------------------------------------------------
c Read Helicities and save information in arrays:
c
c head        helicity and sync read from head bank
c tgbi        helicity and sync read from tgbi bank
c count       helicity scaler read from tgbi bank
c event_head  event number in the ceb ntuple
c hls         helicity and sync read from hls bank
c------------------------------------------------------

      IF( .not. endread) THEN
        

        helpre_head = helcur_head
        helcur_head = helicity
        IF(helcur_head.ne.helpre_head) THEN
           ihel_head=ihel_head+1
           head(ihel_head)=helcur_head
           event_head(ihel_head)=ievent
           write(15,*)  ievent, helcur_head 
       ENDIF

        helpre_tgbi=helcur_tgbi
        helcur_tgbi=0
        call mvbits(helTGBI(1),14,2,helcur_tgbi,0)
        IF(helcur_tgbi.ne.helpre_tgbi) THEN
           ihel_tgbi=ihel_tgbi+1
           write(25,*)ievent, helcur_tgbi,heltgbi(2)
           tgbi(ihel_tgbi)=helcur_tgbi
           count(ihel_tgbi)=heltgbi(2)
        ENDIF

        IF(nHLS.ne.0) THEN
           DO j=1,30
              buffer_hlsb_old(j)=buffer_hlsb_new(j)
              buffer_hlsb_new(j)=0
              buffer_hlsf_old(j)=buffer_hlsf_new(j)
              buffer_hlsf_new(j)=0
              buffer_bit30_old(j)=buffer_bit30_new(j)
              buffer_bit30_new(j)=0
           ENDDO
           DO j=1,NHLS
              call mvbits(c_bpm(j),0,24,buffer_hlsb_new(j),0)
              call mvbits(c_fc (j),0,24,buffer_hlsf_new(j),0)              
              call mvbits(c_fc(j),30,1,buffer_bit30_new(j),0)
           ENDDO
           same_hls_bank=.true.
           DO j=1,NHLS
              same_hls_bank=same_hls_bank.AND.
     &             (buffer_hlsb_new(j).eq.buffer_hlsb_old(j)).AND.
     &             (buffer_hlsf_new(j).eq.buffer_hlsf_old(j)).AND.       
     &             (buffer_bit30_new(j).eq.buffer_bit30_old(j))
           ENDDO
           if(.not.same_hls_bank) then
              DO j=1,NHLS
                 ihel_hls=ihel_hls+1
                 oldbit31=bit31
                 call mvbits(c_fc(j),0,24,FC,0)
                 call mvbits(c_bpm(j),0,24,BPM,0)
                 call mvbits(c_bpm(j),29,1,bit29,0)
                 call mvbits(c_bpm(j),30,1,bit30,0)
                 call mvbits(c_bpm(j),31,1,bit31,0) 
             		
                 write(35,*)  ievent,bit31,bit30,bit29,bpm,fc
c----- check to find wrong bit 31              
                 if(j.lt.NHLS.AND.j.gt.1) then
                    call mvbits(c_bpm(j+1),31,1,newbit31,0)
                    if(oldbit31.eq.bit31.AND.bit31.eq.newbit31) then
                       bit31=1-bit31
                    endif
                 endif
c----------------------------------
c----- in the third algorithm bit 31 is ALWAYS used as sync bit

                 hel_hls=2*bit30+bit31               
                 hls(ihel_hls)=hel_hls
                 nrow_hls(ihel_hls)=NHLS
                 event_hls(ihel_hls)=ievent
                 vBPM(ihel_hls)=BPM
                 vFC (ihel_hls)=FC
              ENDDO
           endif
        ENDIF
      ENDIF

      IF( endread ) GOTO 21000
      ENDDO !DO WHILE (.true)
      ENDIF 
21000 CONTINUE

c--------------------------------------------------------------
c  Close files that contain helicity info
c--------------------------------------------------------------
      CLOSE(15)
      CLOSE(25)
      CLOSE(35)


c--------------------------------------------------------------
c Start helicity matching
c--------------------------------------------------------------
c ihel_HEAD                     dimension of head and tgbi vectors
c ihel_hls                      dimension of hls vector
      

c-----------------------------------------------------------------
c Define the helicity scaler offset:use recorded info if this file 
c and the previous one belong to same run 
c----------------------------------------------------------------- 
     
      IF(run_number.ne.old_run) THEN
         print *,'first file of a new run...'
         lcount(1)=0
         lcount(2)=0
      ENDIF
      offset=(count(1)-lcount(1))*2
      IF(count(1).ne.count(2)) THEN
         offset=offset-lcount(2)
      ELSE
         offset=offset-1-lcount(2)
      ENDIF 
      print *,'offset=',offset
      print *,count(1),count(2)


c------------------------------------------------------------------
c Initialize final vectors
c------------------------------------------------------------------
      DO j=1,5000
         head_new(j)=-1
         tgbi_new(j)=-1
         count_new(j)=-1
         nev_head(j)=-1
         dose_BPM(j)=0
         dose_FC (j)=0
         key(j)=0
         key_pair(j)=0
         key_both(j)=0
      ENDDO

c--- the final vectors will contain the helicity information after 
c--- the matching. That will be the helicity if match was found or 
c--- -1 otherwise. The no match problem can occur when that helicity 
c--- state was missing in the tgbi and head bank or when no agreement 
c--- with the hls was obtained.



c------------------------------------------------------------------
c  Initialize variables to start matching
c------------------------------------------------------------------
c pos   vector containing the position of the missing scaler counts
c off   vector containing the number of missing count at the 
c       corresponding position


      i=1
      j=1
      pos(1)=0
      off(1)=offset
      
      print *,'start matching helicity info'
c----------------------------------------------------------------
c Look for missing count first
c----------------------------------------------------------------
      do j=1,ihel_head-1
         if(count(j+1)-count(j).gt.1) then
            i=i+1
            pos(i)=j !----- last index before missing count
            if(mod(head(j),2).eq.0) then
               dub_bef=1
            else
               dub_bef=0
            endif               
            if(mod(head(j+1),2).eq.0) then
               dub_aft=1
            else
               dub_aft=0
            endif               
            offset=(count(j+1)-count(j))*2
            offset=offset-2+(1-dub_bef)+dub_aft
            off(i)=offset
         elseif(count(j+1)-count(j).eq.1) then
            sync1=0
            sync2=0
            call mvbits(head(j),0,1,sync1,0)
            call mvbits(head(j+1),0,1,sync2,0)
            if(sync2-sync1.ne.1) then
               i=i+1
               pos(i)=j
               off(i)=sync1-sync2+1 
            endif
         endif
      enddo
      ncut=i
      pos(ncut+1)=ihel_head
      if(count(1).lt.20) then
         pos(1)=1
         off(1)=-1
      endif

c------------------------------------------------------------------
c Start matching.....
c------------------------------------------------------------------
      off_diff=0
      new_offset=0
      prec_hls_bank=20
      beam=.true.
      i=0
      do while(i.lt.ncut)
         i=i+1
         done=0
         done20=0
         if(beam.AND.i.gt.1) then
            new_offset=new_offset+mod(off(i),20)
         elseif(beam.AND.i.eq.1) then
            new_offset=new_offset+off(i)
         elseif(off(i).gt.40-prec_hls_bank) then
            new_offset=new_offset+mod(off(i)-(40-prec_hls_bank),20)
         else
            new_offset=new_offset+off(i) 
         endif
         if(off(i).ge.20) then
            flag20=.true.
         else
            flag20=.false.
         endif
20101    j=pos(i)+1
         if(i.eq.1.AND.off(i).gt.30) then
            print*,'helicity states missing before this file'
            done=6
            forward=.true.
            match=.false.
            off_tot=0
            new_offset=0
            goto 20103
         endif
20102    if(pos(i)+1+new_offset+off_diff.gt.ihel_hls) goto 20202 
         off_tot=new_offset+off_diff
         match=.true.
         counter=0
c         if(vFC(j+off_tot).eq.0.AND.(.not.beam)) then
c            match=.false.
c            beam=.false.
c            goto 20102
c         endif
c         beam =.true.
         do while(j.le.pos(i+1).AND.j+off_tot.le.ihel_hls)
            if(j+off_tot-1.gt.0) then
               if(event_hls(j+off_tot).ne.event_hls(j+off_tot-1)) then
                  prec_hls_bank=nrow_hls(j+off_tot-1)
               endif
            endif
            if(vFC(j+off_tot).eq.0.AND.
     &         vFC(min(pos(i+1)+off_tot-1,ihel_hls)).eq.0) then
               beam=.false.
               print *,'no beam'
               goto 20103
            endif
            match=match.AND.head(j).eq.tgbi(j)
     &                 .AND.tgbi(j).eq.hls(j+off_tot)
            if(match.AND.j.gt.pos(i)+1.AND.j.lt.pos(i+1)) then
               counter=counter+1
               key(j)=1
               dose_BPM(j)=vBPM(j+off_tot)
               dose_FC (j)=vFC (j+off_tot)
               head_new(j+off_tot) =head(j)
               tgbi_new(j+off_tot) =tgbi(j)
               count_new(j+off_tot)=count(j)
               nev_head(j+off_tot)=event_head(j)
            endif
            if(.not.match.AND.counter.gt.30) then
               done=6
               forward=.true.
               goto 20103
            endif
            if(.not.match.AND.beam.AND.forward) goto 20103
c            print *,j,j+off_tot,head_new(j+off_tot),tgbi_new(j+off_tot),
c     &              hls(j+off_tot)     
            j=j+1 
         enddo


c----- If mismatch is found try to change offset to recover
20103    if(match) then
            print *,'match found'
            beam=.true.
            forward=.false.
         elseif(.not.match.AND.done.lt.5) then
            j=pos(i)+1
            do while(j.le.pos(i+1).AND.j+off(i).le.ihel_hls)
               key(j)=0
               dose_BPM(j)=0
               dose_FC (j)=0
               head_new(j+off_tot) =-1
               tgbi_new(j+off_tot) =-1
               count_new(j+off_tot)=-1
               nev_head(j+off_tot) =-1
               j=j+1
            enddo
            done=done+1
            off_diff=off_diff-1
         print *,'no match found: retrying with offset changed of',done
         goto 20101
         elseif(.not.match.AND.done.eq.5.AND.flag20) then
            done=0
            off_diff=off_diff+5+20
            done20=done20+1
            if(off_diff.gt.off(i)-mod(off(i),20).OR.done20.gt.2) then
               done=6
               off_diff=off_diff-done20*20
               flag20=.false.
            endif               
         print *, 'missing count # greater than 20; trying offset+20'
         goto 20101
         elseif(.not.match.AND.done.eq.5.AND.(.not.forward)) then
            done=done+1
            off_diff=off_diff+5	
         print *,'not match found after 4 attempt'
         print *,'gone back to initial offset and proceeding forward'
         forward=.true.
         goto 20101
         elseif(.not.match.AND.done.eq.6.AND.forward) then
            print *,'trying to find new alignement....'
c---- look for new HLS bank 
            j_hls=j+off_tot
            j_head=j
            do while(event_hls(j_hls).eq.event_hls(j_hls+1)
     &                .OR.vFC(j_hls+1).eq.0) 
               j_hls=j_hls+1
            enddo
            do while(event_head(j_head).le.event_hls(j_hls))
               j_head=j_head+1
            enddo
            j=j_head-1
            new_offset=j_hls+1-(j_head-1)
            off_diff=0
            do while(pos(i+1).le.j)
               i=i+1
            enddo
            forward=.false.
            done=0
            goto 20102
         endif
      enddo

c---- After 4 attempt whereas succeded or not go to following helicity group

c---- Concludes writing matched helicity in text file
20202 continue
      DO i=1,ihel_hls
         write(45,*) nev_head(i),head_new(i),tgbi_new(i),hls(i),
     &               count_new (i),event_hls(i),nrow_hls(i),
     &               vbpm(i),vfc(i)
      ENDDO
      CLOSE(45)


c-------------------------------------------------------------------
c find pairs of helicity states
c-------------------------------------------------------------------

      n_inv=0
      ind=0
      npair=0

c --- find a cut ----
      j=1


22001 CONTINUE
      DO WHILE(key(j).ne.1.AND.j.le.ihel_head)
         j=j+1
      ENDDO
      start_point=j
      if(start_point.eq.ihel_head) goto 26001
      DO WHILE(key(j).eq.1.AND.j.lt.ihel_head)
         j=j+1
      ENDDO
      if(key(j).eq.1) then
         end_point=j-1
      else
         end_point=j
      endif
      DO i=start_point,end_point-3 
         h1=head(i+0)
         h2=head(i+1)
         h3=head(i+2)
         h4=head(i+3)
         IF( (h1.eq.0 .and. h2.eq.1 .and. h3.eq.0 .and. h4.eq.1) .OR.
     &       (h1.eq.2 .and. h2.eq.3 .and. h3.eq.2 .and. h4.eq.3) ) THEN
            n_inv=i+2
            GOTO 23001 
         ENDIF
      ENDDO
      GOTO 25001                 ! no cut has been found
c --- if the cut was found then go back from ncut to nbufstart ---
c --- taking good events until bad one occures                 ---
23001 CONTINUE
      i=n_inv-1+4
      DO WHILE(i-4-3 .ge. start_point)
         i=i-4
         h1=head(i-3)
         h2=head(i-2)
         h3=head(i-1)
         h4=head(i-0)
         IF(h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3 .OR.
     &        h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)THEN
            npair=npair+1
            key_pair(i-3)=npair
            key_pair(i-2)=npair
            key_pair(i-1)=npair
            key_pair(i-0)=npair
         ELSE
            GOTO 24001 
         ENDIF
      ENDDO
c --- after the cut was found then go forth from ncut ---
c --- taking good events until bad one occures        ---
24001 CONTINUE
      i=n_inv-4
      DO WHILE( i+4+3.le.end_point)
        i=i+4 
        h1=head(i+0)
        h2=head(i+1)
        h3=head(i+2)
        h4=head(i+3)
        IF((h1.eq.0 .and. h2.eq.1 .and. h3.eq.2 .and. h4.eq.3) .OR.
     &     (h1.eq.2 .and. h2.eq.3 .and. h3.eq.0 .and. h4.eq.1)  )THEN
           npair=npair+1
           key_pair(i+0)=npair
           key_pair(i+1)=npair
           key_pair(i+2)=npair
           key_pair(i+3)=npair           
        ELSE
          GOTO 22001 
        ENDIF
      ENDDO
25001 CONTINUE
      IF(end_point.lt.ihel_head) GOTO 22001
c-----------------------------------------------------------------------


26001 ngood=0
      DO i=1,ihel_head
         if(i.eq.1) then
            old_pair=0
         else
            old_pair=key_pair(i-1)
         endif
         if(i.eq.ihel_head) then
            next_ev_num=nevent_read
         else
            next_ev_num=event_head(i+1)
         endif
         st_mult(i)=next_ev_num-event_head(i)
         if(key_pair(i).ne.old_pair.AND.key_pair(i).ne.0) then
            key_both(i+0)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+1)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+2)=key(i)*key(i+1)*key(i+2)*key(i+3)
            key_both(i+3)=key(i)*key(i+1)*key(i+2)*key(i+3)
         endif
         if(key_both(i).eq.1) then
            ngood=ngood+1
         else
            dose_fc(i)=0
            dose_BPM(i)=0
         endif
         write(55,*) event_head(i), st_mult(i)
     &              , head(i), key(i), key_pair(i+0)
     &              , key_both(i), dose_fc(i),dose_BPM(i)

      ENDDO
      print *,'helicity matching succeded at',ngood*100/ihel_head,'%'
      CLOSE(55)




c---  information about the helicity scaler for the last matched event
      last=ihel_hls
      lcount(1)=count_new(last)
      IF(mod(head_new(last),2).eq.0) THEN
         lcount(2)=1
      ELSE
         lcount(2)=0
      ENDIF
      IF(count_new(ihel_hls).eq.-1) THEN
         DO WHILE(count_new(last).eq.-1)
            last=last-1
         ENDDO
         IF(mod(head_new(last),2).eq.0) THEN
            lcount(1)=count_new(last)+int((ihel_hls-last+1)/2)
            lcount(2)=1-mod(ihel_hls-last,2)
         ELSE
            lcount(1)=count_new(last)+int((ihel_hls-last)/2)
            lcount(2)=mod(ihel_hls-last,2)
         ENDIF
      ENDIF
      print *,'last count is',lcount(1),' with mode',lcount(2)

      print *, 'helicity analysis concluded'

      goto 1002 
 







 1002 END













     
c======================================================================
c Match Faraday cup and helicity selection
c======================================================================      


      SUBROUTINE match_sel(file_in)
      IMPLICIT NONE
      INCLUDE "ntpl_goa.inc"
      INCLUDE "select.inc"




      INTEGER*4 i
      INTEGER*4 begin_FCG,begin_hel,begin_pair,pair_index
      INTEGER*4 end_FCG,end_hel,end_pair
      CHARACTER*(*) file_in
      CHARACTER*(80) out



c---- define name of output file
      out     =file_in(1:index(file_in,'h')-2)//'.match'
      OPEN(UNIT=45,FILE=out     ,FORM='formatted')


c---- identify start and stop point for faraday cup gated      
      DO i=1,nevent_read         
         if(Event_OK(i+1).AND.(.not.Event_OK(i))) then
           begin_FCG=i+1
           goto 1001
        endif
      ENDDO
 1001 DO i=begin_FCG,nevent_read
         if(Event_OK(i).AND.(.not.Event_OK(i+1))) then
           end_FCG=i
           goto 1002
        endif
      ENDDO

c----- indentify start point for dose key
      
 1002 i=1
      DO WHILE(event_head(i).lt.begin_FCG)
         i=i+1
      ENDDO
      begin_hel=i


c----- indentify end point for dose key
      DO WHILE(event_head(i)-1.lt.end_FCG)
         i=i+1
      ENDDO
      end_hel=i-2

c----- define key to 0 for events before start and after end
      DO i=1,begin_hel-1
         key(i)=0
      ENDDO
      DO i=end_hel+1,ihel_head
         key(i)=0
      ENDDO

c---- identify start point for pairs      
      i=1
      begin_pair=1
      DO WHILE(begin_pair.lt.begin_hel.AND.i.lt.ihel_head)
         If(key_pair(i).ne.key_pair(i+1)) begin_pair=i+1
         i=i+1
      ENDDO
      if(i.eq.ihel_head) begin_pair=ihel_head

c----- identify end point for pairs
      i=ihel_head
      end_pair=ihel_head
      DO WHILE(end_pair.gt.end_hel.AND.i.gt.1)
         IF(key_pair(i).ne.key_pair(i-1)) end_pair=i-1
         i=i-1      
      ENDDO
      if(i.eq.1) end_pair=begin_pair

c----- define pair key and both key to 0 before start and after end
      DO i=1,begin_pair-1
         key_pair(i)=0
         key_both(i)=0
         dose_BPM(i)=0
         dose_FC(i) =0
      ENDDO
      DO i=end_pair+1,ihel_head
         key_pair(i)=0
         key_both(i)=0
         dose_BPM(i)=0
         dose_FC(i) =0
      ENDDO

c---- save result of final match in output file
      write(45,*) begin_fcg
      do i=begin_hel,end_hel
         write(45,*) event_head(i)-begin_fcg,st_mult(i)
     &              ,head(i),key(i),key_pair(i),key_both(i)
     &              ,dose_BPM(i),dose_FC(i)
      enddo
      write(45,*) end_fcg
      CLOSE(45)
      END
