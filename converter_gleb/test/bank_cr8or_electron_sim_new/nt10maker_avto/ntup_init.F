      subroutine ntup_init(data_type, beam_type)
      
      implicit none
      save

#include "ntup.inc"

      integer data_type
      character beam_type
      
      lmctk_nt_do	= .false.
      lpart_nt_do	= .false.
      lhead_nt_do	= .false.
      lbeam_nt_do	= .false.
      lec_nt_do		= .false.
      lsc_nt_do		= .false.
      lcc_nt_do		= .false.
      lec1_nt_do	= .false.
      lechb_nt_do	= .false.
      lccrc_nt_do	= .false.
      ltrk_nt_do	= .false.
      ltgpb_nt_do       = .false.
      ltagr_nt_do       = .false.      
      lmchd_nt_do       = .false.   
      lseb_nt_do	= .true.   ! we need only SEB 
*      ltrgs_nt_do       = .true.   ! lumi information
      ltrgs_nt_do       = .false.   ! lumi information

      if(data_type.eq.2) then 
         if(beam_type.eq.'g') then !  photon  MC and data
            lmctk_nt_do = .true.
            ltgpb_nt_do = .true.
            ltagr_nt_do = .true.      
         endif
         ltrgs_nt_do = .false.
         lpart_nt_do	= .true.   ! eprod MC
      elseif(data_type.eq.3) then 
         lpart_nt_do	= .true.   ! eprod MC
         lmchd_nt_do    = .true.   ! weight from MCHD bank
      endif

      return
      end






