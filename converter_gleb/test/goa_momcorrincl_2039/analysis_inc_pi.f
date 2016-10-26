
      SUBROUTINE analysis_inc_pi(MC_ana)
c-------------------------------------------------------------------------      
      IMPLICIT NONE
      include "ntpl_goa.inc"
      include "select.inc"
      integer MC_ana,ntype
        real raddeg, degrad,psi(3)
        integer i
        data degrad/0.0174533/
        data raddeg/57.2958/
c+ Variables for CM
c+ i index 
c+         0 = Energy of cm-rest particle
c+         1 = px of cm-rest   particle
c+         2 = py of  cm-rest  particle
c+         3 = pz of  cm-rest  particle
c+         4 = sqared mass of cm-rest   particle
c+         5 = momentum of cm-rest  particle
c+         6 = theta of cm-rest  particle
c+         7 = phi of cm-rest  particle
c+ j index runs over possibles combination (e.g. j=1 (p pi1) j=2 (pi1 pi2) j=3 (p pi2) j=4 pi1 for singol meson) 
      REAL*4 p_vec(0:7,3),p_mes(0:7,4)

c+ Some general setting
c++ Choosing wich entry in the ntuple to fill
c++ Photon beam:
c++ MC only          : 1 <=> MC;  
c++ DATA only        : 1 <=> DA;  
c++ MC.and.DATA      : 1 <=> MC;  2 <=> DA 
c++ Electron beam:
c++ MC only          : 1 <=> MC;  
c++ DATA only        : 1 <=> DA;  
c++ MC.and.DATA      : 1 <=> MC; 2 <=> DA
            ntype = 1 ! By default it fills 1 variable (DATA.or.MC)
             if (bit_mc.eq.2) then ! If both are present: 1st = MC 2nd = data
               if(MC_ana.eq.0) ntype = 2 
             endif 
             indtype = ntype
c-- Ending Choosing wich entry in the ntuple to fill

c+ Filling some general variable
             cosTH_v= (6.-P_EL*cos(th_EL/57.3))/sqrt(q2+omega**2)
             
c+ Starting lab to CM transformation for the pions
         DO i=1,10
            E_PI_CM(i)  = -1000.
            P_PI_CM(i)  = -1000.
            TH_PI_CM(i) = -1000.
            PH_PI_CM(i) = -1000.
            P_hort(i)   = -1000.
         ENDDO

	 nPI_CM=nPI
         DO i=1,nPI
            call lab2cms( 9
     %                   ,E0,E_el,th_el,ph_el,W,Q2
     %                   ,E_pi(i),Px_pi(i),Py_pi(i),Pz_pi(i)
     %                   ,0.,0.,0.,0.
     %                   ,0.,0.,0.,0.
     %                   ,p_vec,p_mes,psi)
            E_PI_CM(i)  = p_mes(0,4) 
            P_PI_CM(i)  = p_mes(5,4)
            TH_PI_CM(i) = p_mes(6,4)
            PH_PI_CM(i) = p_mes(7,4)
            P_hort(i)   = P_EL*sin(TH_EL/57.3)/sqrt(q2+omega**2)*
     %                    sin((ph_pi(i)-ph_EL-180)/57.3)*sin(th_pi(i)/57.3)
         ENDDO
   

      return
      end 


