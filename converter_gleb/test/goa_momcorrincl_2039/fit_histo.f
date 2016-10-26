c###########################################
       subroutine tof_fit_histo()
c###########################################
       IMPLICIT NONE
       include  "ntpl_goa.inc"
       include  "photon.inc"

       INTEGER*4 n_sec, n_pad,num_hist 
       INTEGER*4 au_ran, au_wid,au_par
       REAL*4    min_f,max_f,wid_f,a_p,b_p,c_p
       REAL*4    mean,mean_e,sigma,sigma_e,fit_ok

c             num_hist = see after
              au_ran = 0  ! NO auto range option  
              min_f = 0.9 ! Assuming a time resolution of 300 ps 
              max_f = 0.9 !I fit only Mean +- 3*sigma(=0.9 ns)
              au_wid =0   ! NO auto fit option  
              wid_f = 0.3 ! Assuming a time resolution of 300 ps 
              au_par= 0   ! NO auto fit par fix  option
              a_p = 3.
              b_p = 4.
              c_p = 1.5
        write(1,*) '---------------------------'
        write(1,*) ' TOF  re-calibration'
        write(1,*) ' '
        write(1,*) 'Sec pad      t0_i'

      DO n_sec = 1,6
       DO n_pad = 1,48
	num_hist = n_sec*1000+n_pad
        CALL fit_histo(num_hist,au_ran,min_f,max_f,au_wid,wid_f,au_par,a_p,b_p,c_p,
     %                 mean,mean_e,sigma,sigma_e,fit_ok)
c If mean_val = 0 I assume the initialization value -(tvsst-tag_tof_corr)
        if (mean.ne.0)  t0_i(n_sec,n_pad) = mean  
        CALL HFILL(7120+n_sec,n_pad*1.,t0_i(n_sec,n_pad),1.) 
        write(1,*) n_sec,n_pad,t0_i(n_sec,n_pad)
       END DO
      END DO
        write(1,*) ' Ended tagger re-calibration'
        write(1,*) '----------------------------'
      return
      end

c###########################################
       subroutine tagger_fit_histo()
c###########################################
       IMPLICIT NONE
       include "photon.inc"
       include "ntpl_goa.inc"

       INTEGER*4 num_hist 
       INTEGER*4 au_ran, au_wid,au_par
       REAL*4    min_f,max_f,wid_f,a_p,b_p,c_p
       REAL*4    mean,mean_e,sigma,sigma_e,fit_ok
c+ Fitting ST/Tagger difference : OVERWRITINGn tvsst
              num_hist = 7001
              au_ran = 0  ! NO auto range option  
              min_f = 0.9 ! Assuming a time resolution of 300 ps 
              max_f = 0.9 !I fit only Mean +- 3*sigma(=0.9 ns)
              au_wid =0   ! NO auto fit option  
              wid_f = 0.15 ! Assuming a time resolution of 150 ps 
              au_par= 0   ! NO auto fit par fix  option
              a_p = 3.
              b_p = 4.
              c_p = 1.5
        CALL fit_histo(num_hist,au_ran,min_f,max_f,au_wid,wid_f,au_par,a_p,b_p,c_p,
     %                 mean,mean_e,sigma,sigma_e,fit_ok)
       tvsst = -mean ! new tagger/st calibration
c-
c+  Fitting TOF overall : OVERWRITING tag_tof_corr
              num_hist = 7002
              au_ran = 0  ! NO auto range option  
              min_f = 0.9 ! Assuming a time resolution of 300 ps 
              max_f = 0.9 !I fit only Mean +- 3*sigma(=0.9 ns)
              au_wid =0   ! NO auto fit option  
              wid_f = 0.3 ! Assuming a time resolution of 300 ps 
              au_par= 0   ! NO auto fit par fix  option
              a_p = 3.
              b_p = 4.
              c_p = 1.5

        CALL fit_histo(num_hist,au_ran,min_f,max_f,au_wid,wid_f,au_par,a_p,b_p,c_p,
     %                 mean,mean_e,sigma,sigma_e,fit_ok)
      tag_tof_corr  = -(tvsst+mean) ! new tag_tof_corr  calibration
c-
        write(1,*) '---------------------------------'
        write(1,*) ' Tagger/ST  re-calibration'
        write(1,*) ' '
        write(1,*) ' TAG/ST                    tvsst = ',tvsst
        write(1,*) ' Overall TOF/TG     tag_tof_corr = ',tag_tof_corr
        write(1,*) ' Ended Tagger/ST  re-calibration'
        write(1,*) '---------------------------------'


            return
      end



c#########################################
       subroutine fit_histo(nhist,auto_range,min_to_fit,max_to_fit,
     %            auto_width,width_fit,auto_par,a,b,c,
     %            mean_val,mean_err,sigma_val,sigma_err,ok_val)    
c#########################################
c
c Subroutine to fit single histogram
c 
c INPUT: 
c       nhist = # of histo to fit
c  auto_range = flag to specifie if to use automatic fit range 
c                1 = +- 1.45 hist_devia
c                0 = use the values below
c  min_to_fit = min value of fitted interval 
c  max_to_fit = max value of fitted interval
c  auto_width = flag to specifie if to use automatic gaussian widht 
c                1 = hist_devia/3.
c                0 = use the values below
c  width_fit = min value of fitted interval 
c  auto_par  = flag to specifie if to use automatic gaussian parameter fix:
c                min_mean, max_mean, min_widht, max_width
c                1 = center hist_mean +- 3.0*hist_devia
c                    width  hist_devia    /4 *2.5  
c                0 = center hist_mean +- a * width_fit 
c                    width  width_fit /b  *c
c OUTPUT:
c                  mean_val  
c                  mean_err 
c                  sigma_val 
c                  sigma_err 
c                  ok_val
c
c####################################################

       IMPLICIT NONE
       include  "ntpl_goa.inc"
       include  "photon.inc"
	
      CHARACTER*80 hist_title
      INTEGER*4 hist_ne,hist_nch,hist_nleft,hist_nright
      REAL*4    hist_xlower,hist_xupper,hist_dx,hist_mean,hist_devia
      REAL*4    hist_min,hist_max,hist_int
      REAL*4    mean_val,mean_err,sigma_val,sigma_err,ok_val
      INTEGER*4  tmp
      REAL*4    parfit(10),parerr(10),parstp(10),parmin(10),parmax(10)
      REAL*4    chichi
      REAL*4    hstati,hmin,hmax,hsum
      INTEGER*4 nhist
      INTEGER*4 auto_range, auto_width, auto_par
      REAL*4    min_to_fit,max_to_fit,width_fit,a,b,c
c ***** This loop inside:                                           ***
c ***** The one-cell histograms are fitted by the gaussian function ***
c ***** the result of the fit: peak value and the chisquare error   ***
c ***** of the peak value are stored in the arrays amid(), achi()   ***
c ***** the array amok() contains 1.0 IF the hist. has been fitted  ***
c ***** and 0.0 IF the midvalue used instead of peak value          ***

       
              CALL hnoent(nhist,hist_ne)

c **** IF the hist. is empty to skip it

	      IF (hist_ne .LE. 0) THEN
	        mean_val  =  0.0
	        sigma_val =  0.0
	        mean_err  =  0.0
	        sigma_err =  0.0
	        ok_val    = -1.0
	        GOTO 599
	      ENDIF
	      
c **** get some info about the histogram

	      CALL hgive(nhist,hist_title,hist_nch,
     &                   hist_xlower,hist_xupper,tmp,tmp,tmp,tmp,tmp)
	      hist_dx   = (hist_xupper-hist_xlower)/hist_nch
	      hist_mean = hstati (nhist,1,' ',1)
	      hist_devia= hstati (nhist,2,' ',1)

c **** IF the statistics is low ( .le. 11 entries ) but it is not empty

	      IF (hist_ne .LE. 11) THEN
	        mean_val  = hist_mean  
	        sigma_val = abs(hist_mean)+hist_dx
	        mean_err  =  hist_mean/sqrt(1.*hist_ne)
	        sigma_err =  sigma_val
	        ok_val    = 1.0
	        GOTO 599
	      ENDIF

c **** IF Statistics is good but all the entries are in one bin

	      IF (hist_devia .LE. 0.0) THEN
	        mean_val  = hist_mean
	        sigma_val = abs(hist_mean)+hist_dx
	        mean_err  =  hist_mean/sqrt(1.*hist_ne)
	        sigma_err =  sigma_val	        
	        ok_val    = 1.0
	        GOTO 599
	      ENDIF

c * IF statist. is great enough and the entries are not in one bin *
c * THEN determine the region  (mean +/- xxx) (nn1,nn1+1,...nn2)   *
c * in which the hist will be fitted                               *
c * and get min,max,int values                                     *

              if(auto_range.eq.1) then
 	       CALL hxi(nhist, hist_mean - 1.8*hist_devia, hist_nleft )
	       CALL hxi(nhist, hist_mean + 1.4*hist_devia, hist_nright)
	      else 
               CALL hxi(nhist, hist_mean - min_to_fit, hist_nleft ) 
	       CALL hxi(nhist, hist_mean + max_to_fit, hist_nright)
              endif 
             IF(hist_nleft .LT.1        ) hist_nleft =1
              IF(hist_nright.GT.hist_nch ) hist_nright=hist_nch
	      hist_max = hmax(nhist)
	      hist_min = hmin(nhist)
	      hist_int = hsum(nhist)

c * IF all the entries of the region are now in one bin of the hist. * 
	        
	      IF (hist_nleft .GE. hist_nright) THEN
	        mean_val  = hist_mean
	        sigma_val = abs(hist_mean)+hist_dx
	        mean_err  =  hist_mean/sqrt(1.*hist_ne)
	        sigma_err =  sigma_val	        
	        ok_val    =  1.0
	        GOTO 599
	      ENDIF

c * The time have come to fit the histogram by gaussian *
c * We will fit them by Gaussian + Constant function    *
c * in the restricted area.                             *

ccc parfit (2)  = hist_mean + 0.6*hist_devia
ccc parfit (3)  = 0.4*hist_devia 

	      iquest (11) = hist_nleft
	      iquest (12) = hist_nright
	      parfit (1)  = hist_max/2.0                ! factor
	      parfit (2)  = hist_mean                   ! center
              if (auto_width.eq.1) then 
	       parfit (3)  =  hist_devia/3.0             ! width 
              else
 	       parfit (3)  = width_fit                   ! fixed width 
              endif 
	      parfit (4)  = 0.                          ! c0 * 1
	      parfit (5)  = 0.                          ! c1 * x
	      parstp (1)  = parfit(1)/10.0
	      parstp (2)  = parfit(2)/10.0
	      parstp (3)  = parfit(3)/10.0
	      parstp (4)  = hist_max/20.
	      parstp (5)  = hist_max/20.
	      parmin (1)  = hist_max/2.
              if (auto_par.eq.1) then
	       parmin (2)  = hist_mean-3.0*hist_devia   ! center
	       parmax (2)  = hist_mean+3.0*hist_devia
	       parmin (3)  = hist_devia/4.0             ! width 
	       parmax (3)  = hist_devia*2.5
              else
	       parmin (2)  = hist_mean- a * width_fit ! center
	       parmax (2)  = hist_mean+ a * width_fit 
	       parmin (3)  = width_fit / b               ! width 
	       parmax (3)  = c * width_fit 
              endif

	      parmin (3)  = hist_devia/4.0
	      parmin (4)  = -hist_max/12.
	      parmin (5)  = -hist_max
	      parmax (1)  = hist_max*3.0
	      parmax (4)  = hist_max/12.
	      parmax (5)  = hist_max
	  
	      CALL hfithn (nhist,'G','BFQR',3,
     &                     parfit,parstp,parmin,parmax,parerr,chichi)

c * IF the difference between the fitted peak value (parfit(2)) *
c * and the mean value is greater than hist's deviation THEN    *
c * we use the hist-meanvalue as a peak value                   *
c * IN BOTH CASE THE VALUE WILL BE TAKEN INTO ACCOUN IN         *
c * THE STACK HISTOGRAMS FIT (because the ok_val is set to 1)   * 
	
	      IF (abs(hist_mean-parfit(2)) .GT. 1.0*hist_devia) THEN
	        mean_val = hist_mean
	        sigma_val= hist_devia
	        mean_err = sqrt(hist_devia**2+hist_dx**2)
                sigma_err= sqrt(hist_devia**2+hist_dx**2)
	        ok_val = 1.0
	      ELSE
	        mean_val  = parfit(2)
	        sigma_val = parfit(3)
	        mean_err  = parerr(2)
	        sigma_err = parerr(3)
	        ok_val = 1.0
	      ENDIF
	      

 599  CONTINUE

      return
      end
