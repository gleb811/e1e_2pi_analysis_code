      subroutine book_trgs_nt(id)
      IMPLICIT NONE
      
#include "trgs_nt.inc"

      integer id,nwrite
      
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRAUTH='Avto Tkabladze')
      PARAMETER (CRNAME='BOOK_TRGS_NT')
      
      IF (NWRITE .LT. 1) THEN
        NWRITE = NWRITE + 1
        CRMESS='Booking ntuple for TRGS'
        CALL RECMES(CRNAME,'I',CRMESS)
      ENDIF
      
      call hbname(id,'TRGS',clock_ug,'clock_ug:I'//
     1',clock_rg:I'//
     1',clock_lg:I'//
     1',fcup_ug:I'// 
     1',fcup_rg:I'// 
     1',fcup_lg:I'// 
     1',trig_ug:I'// 
     1',trig_ac:I') 

      return
      end

