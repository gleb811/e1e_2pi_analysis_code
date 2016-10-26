      subroutine fill_trgs_nt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fill_trgs_nt -- fill the TRGS block of the N-tuple, if   c
c it exists.                                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c A Tkabladze, OHIO University                             c
c 12 July 2004                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      SAVE
#include "bcs.inc"
#include "trgs_nt.inc"

      integer mamind, nami
      integer ind
      data nami/0/ 
         
      If (nami.eq.0) nami = mamind(iw,'TRGS')
      
      ind	= nami+1           
      ind 	= IW(ind-1)
      if(ind.eq.0) go to 99

      clock_ug = IW(ind+15)
      clock_rg = IW(ind+17)
      clock_lg = IW(ind+33) 
      fcup_ug  = IW(ind+16)
      fcup_rg  = IW(ind+18)
      fcup_lg  = IW(ind+34) 
      trig_ug  = IW(ind+61)
      trig_ac  = IW(ind+62)
99    return

      end
