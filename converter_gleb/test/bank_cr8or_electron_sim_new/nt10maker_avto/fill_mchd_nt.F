      subroutine fill_mchd_nt

      IMPLICIT NONE
      SAVE
      
#include "bcs.inc"
#include "mchd_nt.inc"

      integer mamind, nami
      integer ind
      data nami/0/ 
         
      If (nami.eq.0) nami = mamind(iw,'MCHD')
      
      ind	= nami+1           

      ind 	= IW(ind-1)
      
      if(ind.eq.0) go to 99
      
      weight = RW(ind+6) 

99    return

      end




