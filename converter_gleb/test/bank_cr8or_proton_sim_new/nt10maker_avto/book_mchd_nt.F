      subroutine book_mchd_nt(id)
      IMPLICIT NONE
      
#include "mchd_nt.inc"

      integer id,nwrite
      
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRAUTH='Avto Tkabladze')
      PARAMETER (CRNAME='BOOK_MCHD_NT')
      
      IF (NWRITE .LT. 1) THEN
        NWRITE = NWRITE + 1
        CRMESS='Booking ntuple for MCHD'
        CALL RECMES(CRNAME,'I',CRMESS)
      ENDIF
      
      call hbname(id,'MCHD',weight,'weight:R')

      return
      end


