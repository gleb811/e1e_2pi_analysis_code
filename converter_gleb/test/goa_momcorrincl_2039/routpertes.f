c*******************************************************************

        integer function eloss(ipart,vertex,cdir,pmom,pvert)

        IMPLICIT NONE
        integer ipart
        real vertex(3),cdir(3),pmom,pvert
        real energy,evert,pcut,dp,mass(2)
        data mass/0.938,0.139/

        if(ipart.ne.1.and.ipart.ne.2) then
           eloss = 0
           return
        endif
        eloss=1
        energy=sqrt(pmom**2+mass(ipart)**2)-mass(ipart)
        if(energy.ge.0.58)then
          call perte(ipart,vertex,cdir,0.58,evert)
          if(evert.eq.0.) then
            eloss=0
            return
          endif
          if(evert.gt.0.6) then
             eloss = 0
             pvert = 0
             return
          endif
          pcut=sqrt((0.58+mass(ipart))**2-mass(ipart)**2)
          pvert=sqrt((evert+mass(ipart))**2-mass(ipart)**2)
          dp=0.58*(pvert-pcut)/energy
          pvert=pmom+dp
          evert=sqrt(pvert**2+mass(ipart)**2)-mass(ipart)
        else
          call perte(ipart,vertex,cdir,energy,evert)
          if(evert.eq.0.) then
            eloss=0
            return
          endif
        endif
        pvert=sqrt((evert+mass(ipart))**2-mass(ipart)**2)

        return
        end




	subroutine perte(TypePart,vertex,cdir,EnerFin,EnerIni)

        implicit NONE

	real extremtube

	real vertex(3),aupoint(3),cdir(3),point(3),start(3),np(3)
        real EnerFin,EnerIni,d_leg
	integer TypePart

	real centre(3),dir(3),rayon,long,alpha,phi,theta
	real pos1(3),pos2(3),dist,intersec(3),norme,csinus
	integer i,w,extrm
	integer milieu(3)
	real epaisseur(3),deltaE

c
c-------------Intersection avec le Start Counter
c

        d_leg=10.2

        alpha=atan2(cdir(2),cdir(1))
        alpha=3.*alpha/3.1416
        w=nint(alpha)
        phi=w*60.

        point(1)=vertex(1)+10.*cdir(1)		
        point(2)=vertex(2)+10.*cdir(2)		
        point(3)=vertex(3)+10.*cdir(3)

        np(1)=cos(phi*3.14159265/180.)
        np(2)=sin(phi*3.14159265/180.)
        np(3)=0.

        start(1)=d_leg*np(1)
        start(2)=d_leg*np(2)
        start(3)=d_leg*np(3)

        call IntersectionPlan(vertex,point,start,np,aupoint,w)


c ---------- Traversee de cible
c	print *,'++++++++++++++++++++Cible'
c		endif

	long=17.
	centre(1)=0.
	centre(2)=0.
	centre(3)=0.
	rayon=6./2.
	dir(1)=0.
	dir(2)=0.
	dir(3)=1.

	call IntersectionCylindre(vertex,aupoint,centre,dir,rayon,pos1,pos2,w)

	if (w.eq.0.or.w.eq.2.or.((vertex(1)-pos1(1))*(vertex(1)-pos2(1))
     +     +(vertex(2)-pos1(2))*(vertex(2)-pos2(2))).gt.0.) then
c		print *,'Intercepte pas la cible !'
c                milieu(1)=1
                enerini=0.
		return
c                epaisseur(1)=0.
c		goto 567
	endif
        if((vertex(3)-pos1(3))*(vertex(3)-pos2(3)).gt.0.) then
                enerini=0.
		return
	endif

c--- ca sort par le cylindre ou les spheres ?
	extrm=0
	if (pos1(3).gt.(long/2.-rayon)) extrm=1
	if (-(long/2.-rayon).gt.pos1(3)) extrm=-1

c--- si ca sort par une sphere on recalcule l'intersection :
	if (extrm.ne.0) then 
c		print *,'Ca sort par une sphere'
		centre(3)=1.*extrm*(long/2.-rayon)
		call IntersectionSphere(vertex,aupoint,centre,rayon,pos1,pos2,w)
		if (w.eq.0) then 
c			print *,'Vertex hors de la cible'
                        enerini=0.
			return
			endif
		if (w.eq.2) then
			pos1(1)=pos2(1)
			pos1(2)=pos2(2)
			pos1(3)=pos2(3)
			endif
		endif

c	print *,'----------------'
c	print *,'Sortie de la cible calculee en ',pos1(1),pos1(2),pos1(3)
	dist=0
	do i=1,3
	dist=dist+(pos1(i)-vertex(i))*(pos1(i)-vertex(i))
	enddo
	dist=sqrt(dist)
c	print *,'Distance parcourue dans la cible : ',dist
	milieu(1)=1
	epaisseur(1)=dist*10.
c				le sous-programme pertes marche avec des mm


c ------------ Traversee du tube en carbone
c	print *,'++++++++++++++++tube carbone'
 567    continue
	rayon=9.4
        extremtube=17.2
	call IntersectionCylindre(vertex,aupoint,centre,dir,rayon,pos1,pos2,w)
	if (w.eq.0.or.w.eq.2.or.((vertex(1)-pos1(1))*(vertex(1)-pos2(1))
     +     +(vertex(2)-pos1(2))*(vertex(2)-pos2(2))).gt.0.) then
			print *,'Pas dans le tube de carbone !!!!!'
                        enerini=0.
			return
			endif
c--- On teste si ca sort par le cylindre, ou par le cone


	if (pos1(3).gt.extremtube.or.w.eq.-1) then
c			print *,'Sortie par le cone avant'
			centre(1)=0.
			centre(2)=0.
			centre(3)=extremtube+9.4/tan(30.*3.14159265/180.)
			alpha=30.
			call IntersectionCone(vertex,aupoint,centre,dir,alpha,
     +							pos1,pos2,w)
			intersec(1)=pos1(1)
			intersec(2)=pos1(2)
			intersec(3)=pos1(3)
			centre(3)=centre(3)+0.15/sin(30.*3.14159265/180.)
			call IntersectionCone(vertex,aupoint,centre,dir,alpha,
     +							pos1,pos2,w)
			else
			intersec(1)=pos1(1)
			intersec(2)=pos1(2)
			intersec(3)=pos1(3)
			rayon=rayon+.15
			call IntersectionCylindre(vertex,aupoint,centre,dir,
     +						rayon,pos1,pos2,w)
			endif
c	print *,'----------------'
c	print *,'Traversee du carbone de '
c	print *,intersec(1),intersec(2),intersec(3)
c	print *,'a'
c	print *,pos1(1),pos1(2),pos1(3)
	dist=0
	do i=1,3
	dist=dist+(pos1(i)-intersec(i))*(pos1(i)-intersec(i))
	enddo
	dist=sqrt(dist)	
c	Print *,'Distance parcourue dans le carbone :',dist
	milieu(2)=2
	epaisseur(2)=dist*10.
c				le sous-programme pertes marche avec des mm
	

c ------------- Traversee des scintillateurs

c----- On commence par calculer dans quel cadran la particule est sortie:
c	print *,'++++++++++++++++++++++++scintillateur'
	alpha=atan2(aupoint(2),aupoint(1))
c	print *,'en x :',aupoint(1),' et en y : ',aupoint(2)
c	print *,'alpha : ',alpha
	alpha=3.*alpha/3.14159265359
c	print *,'Secteur reel : ',alpha
	w=nint(alpha)	
c	print *,'La particule est sortie dans le cadran',w
	phi=w*180./3.
c	print *,'phi :',phi
	theta=90.
	
	if (aupoint(3).gt.22) then
c			print *,'Sortie par un scintillateur avant'
			theta=60.
			endif
c	print *,'theta = ',theta
	dir(1)=sin(theta*3.14159265/180.)*cos(phi*3.14159265/180.)
	dir(2)=sin(theta*3.14159265/180.)*sin(phi*3.14159265/180.)
	dir(3)=cos(theta*3.14159265/180.)

c	do i=1,3
c	print *,'dir(',i,')=',dir(i)
c	enddo
c  ....................dir contient  un vecteur unitaire perpendiculaire
c					au scintillateur

	norme=0.
	do i=1,3
	norme=norme+(aupoint(i)-vertex(i))*(aupoint(i)-vertex(i))
	enddo
	norme=sqrt(norme)
c	print *,'norme =',norme
	
	csinus=0.
	do i=1,3
	csinus=csinus+(aupoint(i)-vertex(i))*dir(i)
	enddo
c	print *,'produit scalaire =',csinus

	csinus=csinus/norme
c	print *,'cosinus = ',csinus
	
	rayon=0.3
c			epaisseur des scintillateurs
	dist=rayon/csinus

c	print *,'Distance parcourue dans les scintillateurs',dist
	milieu(3)=3
	epaisseur(3)=dist*10.
c				le sous-programme pertes marche avec des mm

	

	call InitRange()
c				On initialise la table des ranges

	EnerIni=EnerFin*1000.
c				le s-prog pertes marche avec des Mev
	do i=3,1,-1
	call pertes(TypePart,milieu(i),epaisseur(i),EnerIni,deltaE)
	EnerIni=EnerIni-deltaE
	enddo
	EnerIni=EnerIni/1000.
	return
	end

c************************************************************************
	subroutine pertes(particule,milieu,epaisseur,EnergieFinale,perts)

        implicit NONE

c
c	particule :	1	proton
c			2	pion
c
c	milieu :	1	H1
c			2	carbone
c			3	scintillateur
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	EnergieFinale est donne en Mev   !
c	epaisseur donnee en milimetres   !
c	pertes rendues en Mev		 !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	common/ran_data/sc_pr(150),h1_pr(150),sc_pi(150),h1_pi(150)
        real sc_pr,h1_pr,sc_pi,h1_pi

	integer particule,milieu
	real	epaisseur,EnergieFinale
	real perts

	double precision rangei,rangef,rangem,coefCsurSC
	integer er

	double precision Emin,Emax,Emilieu

	coefCsurSC=1.035/1.7
	
	Emin=EnergieFinale

	if (milieu.eq.1) then
		if (particule.eq.1) call interpole(h1_pr,150,4.,4.,Emin,rangef,er)
		if (particule.eq.2) call interpole(h1_pi,150,4.,4.,Emin,rangef,er)
		endif
	if (milieu.eq.2.or.milieu.eq.3) then
		if (particule.eq.1) call interpole(sc_pr,150,4.,4.,Emin,rangef,er)
		if (particule.eq.2) call interpole(sc_pi,150,4.,4.,Emin,rangef,er)
		endif
	if (milieu.eq.2) rangef=rangef*coefCsurSC
 
	if (er.eq.-1) then
c			print *,'Hors limites'
			perts=0
			return
			endif

c	print *,'range final : ',rangef

	rangei=rangef+epaisseur

c	print *,'range initial ',rangei

c....... petite dichotomie

	Emin=EnergieFinale
	Emax=600
	
2	Emilieu=(Emin+Emax)/2.
c	print *,'Energies :',Emin,'  <>   ',Emax
	if (milieu.eq.1) then
		if (particule.eq.1) call 
     +			interpole(h1_pr,150,4.,4.,Emilieu,rangem,er)
		if (particule.eq.2) call
     +			interpole(h1_pi,150,4.,4.,Emilieu,rangem,er)
		endif
	if (milieu.eq.2.or.milieu.eq.3) then
		if (particule.eq.1) call
     +			interpole(sc_pr,150,4.,4.,Emilieu,rangem,er)
		if (particule.eq.2) call
     +			interpole(sc_pi,150,4.,4.,Emilieu,rangem,er)
		endif
	if (milieu.eq.2) rangem=rangem*coefCsurSC
	if (rangem.gt.rangei) then
				Emax=Emilieu
				else
				Emin=Emilieu
				endif
	if (abs(rangem-rangei).gt.(epaisseur*1e-4).and.
     +				(Emax-Emin).gt.1e-2) goto 2

	perts=EnergieFinale-Emilieu
	if (milieu.eq.2) perts=perts*1.096
	return
	end


c************************************************************************
	subroutine interpole(tab,t,xmin,pas,x,resul,err)


        implicit NONE

c	Cette fonction retourne la valeur interpolee d'une fonction
c	dont on ne connait qu'un echantillonage discret
c	Les valeurs echantillonnees sont contenues dans le tableau tab
c	de taille t.
c	xmin est la valeur minimum de l'echantillonage
c	pas est le pas d'echantillonage 
c		x=xmin  ->   tab(1)
c		x=xmin+pas -> tab(2)
c		etc...
c		
c		x est la valeur de la variable 
c		err prend la valeur -1 si x hors limites

	real xmin,pas,xi(3),yi(3)
	double precision x,resul,dx
	real tab,resultat,parabole
	integer err,t
	dimension tab(t)
	
	integer i


	err=0
	if (xmin.gt.x.or.x.gt.(xmin+(t-1)*pas)) then
				err=-1
				return
				endif
c	print *,'x : ',x
	dx=(x-xmin)/pas
c	print *,'indice non arrondi',dx
	i=(x-xmin)/pas
c	print *,'Indice arrondi',i

	dx=(x-(xmin+pas*float(i)))/pas
c					x=xmin+i*pas+dx
c	print *,'dx ',dx
	
	xi(1)=-pas
c	print *,'i=',i

	if (i.eq.0) then
			xi(1)=0
			i=1
			endif
	xi(2)=xi(1)+pas
	xi(3)=xi(2)+pas
	yi(1)=tab(i)
	yi(2)=tab(i+1)
	yi(3)=tab(i+2)

c	do i=1,3
c	print *,'tab(',xi(i),')=',yi(i)
c	enddo

	resultat=parabole(yi,xi,sngl(dx*pas))
c 				on interpole par une parabole
	
	

c	resul=tab(i+1)+(tab(i+2)-tab(i+1))*dx
c					on interpole par une bete droite
c	print *,'resultat ',resultat
	resul=dble(resultat)

	return
	end


c****************************************************************************

	real function parabole(yi,xi,x)

c****************************************************************************


        implicit NONE

	real yi(3),xi(3),yc(5),xc(5),x
	integer i
	double precision a,b,c,det

	do i=1,3
	yc(i)=yi(i)
	xc(i)=xi(i)
	enddo
	do i=4,5
	yc(i)=yi(i-3)
	xc(i)=xi(i-3)
	enddo


	a=0
	b=0
	c=0
	do i=1,3
	c=c+xc(i+1)*xc(i+2)*(xc(i+2)-xc(i+1))*yc(i)
	b=b-(xc(i+2)+xc(i+1))*(xc(i+2)-xc(i+1))*yc(i)
	a=a+(xc(i+2)-xc(i+1))*yc(i)
	enddo

	det=(xc(3)-xc(1))*(xc(3)-xc(2))*(xc(2)-xc(1))
	a=a/det
	b=b/det
	c=c/det
	
	parabole=a*x**2+b*x+c
	return
	end

c**********************************************************************

	subroutine IntersectionCylindre(m1,m2,m0,d,r,m,mm,rep)

c**********************************************************************


        implicit NONE

	real m1(3),m2(3),m0(3),d(3),m(3),mm(3),r
	real mm1(4),mm2(4),mm0(4),dd(4)
	integer rep,i
	real a,b,c,k1,k2,kt,delta,norme

	norme=0
	do i=1,3
	norme=norme+d(i)*d(i)
	enddo
	
	if (1e-10.gt.norme) then 
				rep=-3
				return
				endif

	do i=1,3
	d(i)=d(i)/norme
	enddo

	do i=1,3
	mm1(i)=m1(i)
	mm2(i)=m2(i)
	mm0(i)=m0(i)
	dd(i)=d(i)
	enddo

	mm1(4)=m1(1)
	mm2(4)=m2(1)
	mm0(4)=m0(1)
	dd(4)=d(1)

c	do i=1,4
c	enddo


	a=0
	do i=1,3

	a=a+(m2(i)-m1(i))*(m2(i)-m1(i))*(1-d(i)*d(i))
     +			-2*(m2(i)-m1(i))*(mm2(i+1)-mm1(i+1))*d(i)*dd(i+1)
	enddo
	
	if (1e-10.gt.abs(a)) then
			rep=0
			return
			endif

	b=0
	do i=1,3
	b=b+(m2(i)-m1(i))*(m1(i)-m0(i))*(1-d(i)*d(i))
	b=b-(m2(i)-m1(i))*(mm1(i+1)-mm0(i+1))*d(i)*dd(i+1)
	b=b-(m1(i)-m0(i))*(mm2(i+1)-mm1(i+1))*d(i)*dd(i+1)
	enddo

	c=0
	do i=1,3
	c=c+(m1(i)-m0(i))*(m1(i)-m0(i))*(1-d(i)*d(i))
	c=c-2*(m1(i)-m0(i))*(mm1(i+1)-mm0(i+1))*d(i)*dd(i+1)
	enddo
	c=c-r*r
	
	delta=b*b-a*c
	if (0.gt.delta) then
			rep=0
			return
			endif

	k1=(-b+sqrt(delta))/a
	k2=(-b-sqrt(delta))/a


	rep=-1
	if (k1.gt.0.and.1.gt.k1) rep=1
	if (k2.gt.0.and.1.gt.k2) then
			if (rep.eq.1) then
					rep=2
					else
					rep=1
					kt=k1
					k1=k2
					k2=kt
					endif
			endif


	if (rep.eq.2.or.rep.eq.-1) then
			if (abs(k1).gt.abs(k2)) then
					kt=k1
					k1=k2
					k2=kt
					endif
			endif
							
	do i=1,3
	m(i)=m1(i)+k1*(m2(i)-m1(i))
	mm(i)=m1(i)+k2*(m2(i)-m1(i))
	enddo


	return
	end	
	

c ******************************************************************

	subroutine IntersectionPlan(m1,m2,m0,n,m,rep)

c ******************************************************************


        implicit NONE

	real m1(3),m2(3),m0(3),n(3),m(3)
	integer rep,i
	real kd,k

	if (n(1).eq.0.and.n(2).eq.0.and.n(3).eq.0) then
					rep=-3
					return
					endif	
	rep=0
	k=0
	do i=1,3 
	k=k+(m1(i)-m0(i))*n(i)
	enddo
	
	if (1e-10.gt.abs(k)) rep=-2

	kd=0
	do i=1,3
	kd=kd+(m2(i)-m1(i))*n(i)
	enddo

	if (1e-10.gt.kd) then
		if (rep.ne.-2) rep=0
		return
		endif
	k=0
	do i=1,3
	k=k+(m0(i)-m1(i))*n(i)
	enddo
	k=k/kd

	if (k.gt.0.and.1.gt.k)	then 
		rep=1 
		else 
		rep=-1
		endif

	do i=1,3
	m(i)=m1(i)+k*(m2(i)-m1(i))
	enddo
	return
	end	
	
c ******************************************************************

	subroutine IntersectionSphere(m1,m2,m0,r,m,mm,rep)

c ******************************************************************


        implicit NONE

	real m1(3),m2(3),m0(3),m(3),mm(3),r
	integer rep,i
	real a,b,c,k1,k2,kt,delta

	a=0
	do i=1,3
	a=a+(m2(i)-m1(i))*(m2(i)-m1(i))
	enddo

	if (1e-10.gt.abs(a)) then
			rep=0
			return
			endif

	b=0
	do i=1,3
	b=b+(m2(i)-m1(i))*(m1(i)-m0(i))
	enddo

	c=0
	do i=1,3
	c=c+(m1(i)-m0(i))*(m1(i)-m0(i))
	enddo
	c=c-r*r

	delta=b*b-a*c
	if (0.gt.delta) then
			rep=0
			return
			endif

	k1=(-b+sqrt(delta))/a
	k2=(-b-sqrt(delta))/a

	rep=-1
	if (k1.gt.0.and.1.gt.k1) rep=1
	if (k2.gt.0.and.1.gt.k2) then
			if (rep.eq.1) then
					rep=2
					else
					rep=1
					kt=k1
					k1=k2
					k2=kt
					endif
			endif


	if (rep.eq.2.or.rep.eq.-1) then
			if (abs(k1).gt.abs(k2)) then
					kt=k1
					k1=k2
					k2=kt
					endif
			endif
							
	do i=1,3
	m(i)=m1(i)+k1*(m2(i)-m1(i))
	mm(i)=m1(i)+k2*(m2(i)-m1(i))
	enddo

	return
	end	

c ******************************************************************

	subroutine IntersectionCone(m1,m2,m0,d,alpha,m,mm,rep)

c ******************************************************************



        implicit NONE

	real m1(3),m2(3),m0(3),d(3),m(3),mm(3),alpha
	real mm1(4),mm2(4),mm0(4),dd(4)
	integer rep,i
	real a,b,c,k1,k2,kt,delta,norme,c2

	c2=cos(alpha/180.*3.14159265359)
	c2=c2*c2

	norme=0
	do i=1,3
	norme=norme+d(i)*d(i)
	enddo
	
	if (1e-10.gt.norme) then 
				rep=-3
				return
				endif

	do i=1,3
	d(i)=d(i)/norme
	enddo

	do i=1,3
	mm1(i)=m1(i)
	mm2(i)=m2(i)
	mm0(i)=m0(i)
	dd(i)=d(i)
	enddo

	mm1(4)=m1(1)
	mm2(4)=m2(1)
	mm0(4)=m0(1)
	dd(4)=d(1)

	a=0
	do i=1,3

	a=a+(m2(i)-m1(i))*(m2(i)-m1(i))*(d(i)*d(i)-c2)
     +			+2*(m2(i)-m1(i))*(mm2(i+1)-mm1(i+1))*d(i)*dd(i+1)
	enddo
	
	if (1e-10.gt.abs(a)) then
			rep=0
			return
			endif

	b=0
	do i=1,3
	b=b+(m2(i)-m1(i))*(m1(i)-m0(i))*(d(i)*d(i)-c2)
	b=b+(m2(i)-m1(i))*(mm1(i+1)-mm0(i+1))*d(i)*dd(i+1)
	b=b+(m1(i)-m0(i))*(mm2(i+1)-mm1(i+1))*d(i)*dd(i+1)
	enddo

	c=0
	do i=1,3
	c=c+(m1(i)-m0(i))*(m1(i)-m0(i))*(d(i)*d(i)-c2)
	c=c+2*(m1(i)-m0(i))*(mm1(i+1)-mm0(i+1))*d(i)*dd(i+1)
	enddo
	
	delta=b*b-a*c
	if (0.gt.delta) then
			rep=0
			return
			endif

	k1=(-b+sqrt(delta))/a
	k2=(-b-sqrt(delta))/a

	rep=-1
	if (k1.gt.0.and.1.gt.k1) rep=1
	if (k2.gt.0.and.1.gt.k2) then
			if (rep.eq.1) then
					rep=2
					else
					rep=1
					kt=k1
					k1=k2
					k2=kt
					endif
			endif


	if (rep.eq.2.or.rep.eq.-1) then
			if (abs(k1).gt.abs(k2)) then
					kt=k1
					k1=k2
					k2=kt
					endif
			endif
							
	do i=1,3
	m(i)=m1(i)+k1*(m2(i)-m1(i))
	mm(i)=m1(i)+k2*(m2(i)-m1(i))
	enddo


	return
	end	
	


c **************************************************************************
	subroutine InitRange()

        implicit NONE

c *
c *  - for initialization of data table
c *  - max kinetic energy = 600 MeV both for protons and pions
c * 
c *	energie des particules 4 Mev -> 600 Mev par pas de 4Mev
c *	      
c *       XX_pr(i) = range for protons in material XX ( mm )
c *       XX_pi(i) = range for pions   in material XX ( mm )
c *                  where XX can be:
c *                    sc = scintillator
c*		       h1 = hydrogen
c *
c **************************************************************************
c *
	common/ran_data/sc_pr(150),h1_pr(150),sc_pi(150),h1_pi(150)
        real sc_pr,h1_pr,sc_pi,h1_pi

c
	data sc_pr/
     +   0.2301,   0.7882,   1.6362,   2.7515,   4.1142,
     +   5.7374,   7.5904,   9.6697,  11.9568,  14.4648,
     +  17.3025,  20.2227,  23.3460,  26.6687,  30.1872,
     +  33.8977,  37.7969,  41.8815,  46.1481,  50.5936,
     +  55.2148,  60.0087,  64.9724,  70.1030,  75.3977,
     +  80.8538,  86.4687,  92.2398,  98.1646, 104.2407,
     + 110.4656, 116.8371, 123.3529, 130.0109, 136.8089,
     + 143.7448, 150.8166, 158.0223, 165.3600, 172.8279,
     + 180.4241, 188.1467, 195.9942, 203.9647, 212.0567,
     + 220.2685, 228.5985, 237.0451, 245.6069, 254.2825,
     + 263.0703, 271.9689, 280.9769, 290.0930, 299.3160,
     + 308.6444, 318.0770, 327.6126, 337.2498, 346.9876,
     + 356.8247, 366.7599, 376.7922, 386.9203, 397.1433,
     + 407.4600, 417.8691, 428.3699, 438.9612, 449.6418,
     + 460.4111, 471.2676, 482.2105, 493.2391, 504.3520,
     + 515.5483, 526.8273, 538.1880, 549.6293, 561.1504,
     + 572.7501, 584.4280, 596.1827, 608.0136, 619.9197,
     + 631.9003, 643.9542, 656.0808, 668.2791, 680.5483,
     + 692.8875, 705.2961, 717.7729, 730.3174, 742.9286,
     + 755.6056, 768.3479, 781.1544, 794.0247, 806.9575,
     + 819.9524, 833.0084, 846.1252, 859.3015, 872.5369,
     + 885.8304, 899.1818, 912.5897, 926.0540, 939.5739,
     + 953.1487, 966.7776, 980.4599, 994.1953,1007.9833,
     +1021.8230,1035.7139,1049.6553,1063.6473,1077.6885,
     +1091.7793,1105.9183,1120.1060,1134.3416,1148.6245,
     +1162.9543,1177.3317,1191.7552,1206.2246,1220.7405,
     +1235.3019,1249.9089,1264.5619,1279.2600,1294.0029,
     +1308.7922,1323.6265,1338.5059,1353.4309,1368.4016,
     +1383.4185,1398.4814,1413.5908,1428.7473,1443.9508,
     +1459.2012,1474.5009,1489.8494,1505.2472,1520.6951/
c
	data H1_PR/
     +    1.373,    4.879,   10.302,   17.510,   26.376,
     +   36.997,   49.171,   62.877,   77.994,   94.611,
     +  113.085,  132.514,  153.332,  175.514,  199.036,
     +  223.874,  250.006,  277.409,  306.062,  335.943,
     +  367.032,  399.308,  432.753,  467.347,  503.072,
     +  539.910,  577.843,  616.854,  656.926,  698.043,
     +  740.189,  783.350,  827.509,  872.652,  918.765,
     +  965.834, 1013.846, 1062.786, 1112.643, 1163.403,
     + 1215.054, 1267.585, 1320.983, 1375.238, 1430.337,
     + 1486.270, 1543.026, 1600.595, 1658.967, 1718.130,
     + 1778.077, 1838.796, 1900.279, 1962.517, 2025.500,
     + 2089.219, 2153.666, 2218.832, 2284.709, 2351.289,
     + 2418.562, 2486.523, 2555.162, 2624.471, 2694.444,
     + 2765.073, 2836.349, 2908.266, 2980.818, 3053.996,
     + 3127.793, 3202.203, 3277.218, 3352.833, 3429.039,
     + 3505.831, 3583.203, 3661.146, 3739.656, 3818.724,
     + 3898.347, 3978.515, 4059.226, 4140.470, 4222.242,
     + 4304.537, 4387.348, 4470.669, 4554.495, 4638.819,
     + 4723.637, 4808.939, 4894.724, 4980.984, 5067.714,
     + 5154.908, 5242.561, 5330.667, 5419.221, 5508.217,
     + 5597.649, 5687.515, 5777.807, 5868.521, 5959.650,
     + 6051.193, 6143.141, 6235.491, 6328.240, 6421.379,
     + 6514.908, 6608.821, 6703.113, 6797.778, 6892.817,
     + 6988.223, 7083.990, 7180.118, 7276.602, 7373.438,
     + 7470.623, 7568.155, 7666.028, 7764.244, 7862.796,
     + 7961.684, 8060.903, 8160.457, 8260.334, 8360.541,
     + 8461.073, 8561.931, 8663.113, 8764.615, 8866.439,
     + 8968.592, 9071.062, 9173.855, 9276.972, 9380.414,
     + 9484.182, 9588.277, 9692.698, 9797.456, 9902.548,
     +10007.974,10113.740,10219.855,10326.317,10433.134/
c
c
	data sc_pi/
     +   1.0489,   3.6855,   7.6072,  12.7846,  18.9594,
     +  25.9762,  34.0466,  42.8675,  52.2316,  62.5075,
     +  73.2553,  84.3481,  96.5312, 108.7142, 121.1779,
     + 134.7961, 148.4144, 162.0327, 176.1488, 191.1414,
     + 206.1341, 221.1267, 236.1194, 251.7828, 268.0259,
     + 284.2689, 300.5120, 316.7550, 332.9980, 349.7382,
     + 367.0592, 384.3802, 401.7012, 419.0222, 436.3432,
     + 453.6642, 470.9852, 489.0165, 507.2115, 525.4066,
     + 543.6015, 561.7965, 579.9915, 598.1864, 616.3814,
     + 634.5764, 652.8911, 671.7372, 690.5833, 709.4293,
     + 728.2753, 747.1215, 765.9673, 784.8134, 803.6594,
     + 822.5055, 841.3514, 860.1974, 879.0702, 898.3461,
     + 917.6221, 936.8979, 956.1740, 975.4499, 994.7256,
     +1014.0015,1033.2775,1052.5535,1071.8293,1091.1053,
     +1110.3812,1129.6570,1148.9330,1168.2090,1187.6671,
     +1207.1704,1226.6737,1246.1771,1265.6805,1285.1838,
     +1304.6873,1324.1907,1343.6938,1363.1973,1382.7007,
     +1402.2040,1421.7073,1441.2107,1460.7140,1480.2173,
     +1499.7206,1519.2239,1538.7274,1558.2635,1577.8225,
     +1597.3817,1616.9408,1636.4998,1656.0590,1675.6179,
     +1695.1770,1714.7362,1734.2952,1753.8542,1773.4132,
     +1792.9723,1812.5312,1832.0905,1851.6494,1871.2085,
     +1890.7677,1910.3267,1929.8855,1949.4446,1969.0039,
     +1988.5627,2008.1219,2027.6488,2047.1282,2066.6069,
     +2086.0862,2105.5654,2125.0447,2144.5237,2164.0029,
     +2183.4819,2202.9614,2222.4404,2241.9194,2261.3987,
     +2280.8779,2300.3572,2319.8362,2339.3154,2358.7947,
     +2378.2742,2397.7529,2417.2319,2436.7114,2456.1907,
     +2475.6697,2495.1489,2514.6284,2534.1069,2553.5862,
     +2573.0654,2592.5444,2612.0239,2631.3318,2650.6313/
c
c
	data H1_PI/
     +     6.78,    24.22,    50.85,    85.48,   127.25,
     +   175.86,   230.30,   290.22,   354.93,   423.86,
     +   498.62,   576.05,   657.61,   742.56,   830.57,
     +   921.20,  1015.61,  1110.02,  1210.44,  1311.10,
     +  1412.54,  1519.18,  1625.82,  1732.46,  1843.89,
     +  1956.13,  2068.38,  2180.95,  2298.41,  2415.88,
     +  2533.34,  2650.80,  2771.58,  2893.94,  3016.30,
     +  3138.65,  3261.01,  3386.07,  3512.89,  3639.71,
     +  3766.53,  3893.35,  4020.17,  4149.82,  4280.50,
     +  4411.17,  4541.84,  4672.52,  4803.19,  4933.87,
     +  5067.58,  5201.48,  5335.37,  5469.27,  5603.17,
     +  5737.06,  5870.96,  6005.08,  6141.57,  6278.06,
     +  6414.56,  6551.05,  6687.54,  6824.03,  6960.53,
     +  7097.02,  7233.50,  7371.98,  7510.46,  7648.95,
     +  7787.43,  7925.91,  8064.39,  8202.88,  8341.36,
     +  8479.84,  8618.32,  8757.48,  8897.39,  9037.30,
     +  9177.21,  9317.11,  9457.02,  9596.92,  9736.83,
     +  9876.74, 10016.64, 10156.55, 10296.46, 10436.65,
     + 10577.46, 10718.28, 10859.10, 10999.91, 11140.73,
     + 11281.55, 11422.37, 11563.18, 11704.00, 11844.82,
     + 11985.63, 12126.45, 12267.27, 12408.19, 12549.46,
     + 12690.73, 12832.01, 12973.28, 13114.56, 13255.83,
     + 13397.10, 13538.38, 13679.65, 13820.92, 13962.20,
     + 14103.47, 14244.74, 14386.02, 14527.29, 14668.56,
     + 14809.89, 14951.23, 15092.57, 15233.91, 15375.25,
     + 15516.59, 15657.92, 15799.26, 15940.60, 16081.94,
     + 16223.28, 16364.61, 16505.95, 16647.29, 16788.63,
     + 16929.96, 17071.30, 17212.64, 17353.96, 17495.03,
     + 17636.10, 17777.17, 17918.25, 18059.32, 18200.39,
     + 18341.47, 18482.54, 18623.61, 18764.69, 18905.76/
c
	
	return	
	end
