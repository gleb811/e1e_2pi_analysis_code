C ****************************************************************************
	SUBROUTINE KINFIT(INUC)
C *
C * kinematical fit : Y + Nuc --> 1 + 2 + 3
C * CONLES package
C *
C * additional routines from CERNLIB
C * - CONLES package
C * - VZERO
C * - REQN
C * - ATG
C * - PROB
C *
C * 10 kinematical variables : 3 for each particle + photon energy
C *
C * ---------- COMMON BLOCK
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
C *
C * ---------- INPUT
C * INUC       nucleus type :
C *              1 = Hydrogen
C *              2 = Deuterium
C *              3 = He-3
C *              4 = neutron (quasi-free)
C * ID_PART(3) particle type :
C *              = 1 proton
C *              = 2 pi +/-
C *              = 3 kaon
C *              = 4 pi 0
C * PM_KIN(3) momentum value for each particle (MeV/c) ...
C * EP_KIN(3) ... error
C * TH_KIN(3) theta value for each particle (rad) ...
C * ET_KIN(3) ... error
C * FI_KIN(3) phi value for each particle (rad) ...
C * EF_KIN(3) ... error
C * PM_GAM    photon energy (MeV) ...
C * EP_GAM    ... error
C *
C * ---------- NOTE
C *            - units: Mev , MeV/c , radiants
C *            - set errors =0 for unmeasured variables
C *            - after fit the previous variables contain the new values if
C *              convergence was reached or old value if not (flag IER_KIN)
C *
C * ---------- OUTPUT
C * PP_INV(3)  invariant masses 1-2 , 2-3 , 1-3 (MeV)
C * NST_KIN    nr of steps
C * C4_KIN     constraints value = SUM{C_i}
C * CHI_KIN    chi-square
C * CL_KIN     Confidence Level (< 0.001)
C * IER_KIN    = 0 convergence reached. if not:
C *            = 1 probability (C.L.) less than CST(1) (except first iteration)
C *            = 2 nr of iterations greater than NST(1)
C *            = 3 total nr of cutsteps greater than NST(2)
C *            = 4 matrix singular
C *            = 5 unphysical values at first iteration
C *            = 6 nr of variable (NY) or nr of constraints (NF) too large
C *            = 7 not enough informations
C *            = 8 error in particle or nucleus number
C *            = 9 error in recalculated missing informations
C *            =10 1 miss p
C *            =11 2 miss p
C *            =12 3 miss p : matrix inversion error
C *            =13 3 miss p : p<0
C *            =14 1 miss p + 1 miss part : E<0
C *            =15 1 miss p + 1 miss part : cos(theta)>1
C *            =16 2 miss p + 1 miss part : index error
C *            =17 2 miss p + 1 miss part : solve error
C *
C ****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
C
      REAL*8 V(1275),DELY(50),F(20),B(1000),PL(50),CHSQ,PR
      COMMON/CONF1/V,DELY,F,B,PL,
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      REAL*8 DET,FTEST,FTESTL,CHSQV,DELYP(50),CST(4)
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP,
     1   CST,NSTEP,NCUT,NCST,NST(3)
C
      DIMENSION TAB_PAR(4),TAB_NUC(4)
      DATA TAB_PAR /938.27231,139.568,493.677,134.974/
      DATA TAB_NUC /938.27231,1876.092,2809.436,939.56563/
      DATA PI /3.141592653/
C
        do i=1,10
	  SP_VAR(I)=0.
	  ER_VAR(I)=0.
	  FT_VAR(I)=0.
	ENDDO
	DO 1 I=1,3
1	PP_INV(I)=0.
C
	NST_KIN = 0
	C4_KIN  = 0.
	CHI_KIN = 0.
	CL_KIN  = 0.
	IER_KIN = 0
C
	IF (0.GE.INUC.OR.INUC.GT.4) THEN
	  IER_KIN=8
	  GOTO 114
	ENDIF
	AMNUC=TAB_NUC(INUC)
C
	DO I=1,3
	  IPT=ID_PART(I)
	  IF (0.GE.IPT.OR.IPT.GT.4) THEN
	    IER_KIN=8
	    GOTO 114
	  ENDIF
	  AMASS(I)=TAB_PAR(IPT)
	ENDDO
C
	IU1=0
	IU2=0
	IU3=0
	DO I=1,3
	  IF (EP_KIN(I).EQ.0.) IU1=IU1+1
	  IF (ET_KIN(I).EQ.0.) IU2=IU2+1
	  IF (EF_KIN(I).EQ.0.) IU3=IU3+1
	ENDDO
	ITOT=IU1+IU2+IU3
C ============================================================================
C                                                   1) not enough informations
	IF (ITOT.GT.4) THEN
	  IER_KIN=7
	  GOTO 114
	ENDIF
C ============================================================================
C                                                   2) no missing informations
	IF (ITOT.EQ.0) GOTO 100
C ============================================================================
C                                                3) max 3 missing informations
	IMIS=0
	IF (ITOT.LT.4) THEN
C ----------------------------------------------------------------------------
C                                                    a) 1 to 3 missing momenta
	  IF (IU1.GT.0.AND.IU2.EQ.0.AND.IU3.EQ.0) THEN
	    IMIS=1
	    IF (IU1.EQ.1) CALL RECALC1_1
	    IF (IU1.EQ.2) CALL RECALC1_2
	    IF (IU1.EQ.3) CALL RECALC1_3
	  ENDIF
C ----------------------------------------------------------------------------
C                                          b) 1 particle missing (p,theta,phi)
	  IF (IU1.EQ.1.AND.IU2.EQ.1.AND.IU3.EQ.1) THEN
	    IOUT=0
	    DO I=1,3
	      IF ((EP_KIN(I).EQ.0.).AND.(ET_KIN(I).EQ.0.).AND.
     +            (EF_KIN(I).EQ.0.)) IOUT=I
	    ENDDO
	    IF (IOUT.GT.0) THEN
	      IMIS=2
	      CALL RECALC2(IOUT)
	    ENDIF
	  ENDIF
C ============================================================================
C                            4) 4 missing info: 1 particle + 1 other momentum.
C                                in this case no fit
	ELSE
C
	  IF (IU1.EQ.2.AND.IU2.EQ.1.AND.IU3.EQ.1) THEN
	    IMIS=3
	    CALL RECALC3
C -ale 4/11/93 : no fit if 4 missing info
	    IF (IER_KIN.EQ.0) GOTO 113
C -ale 4/11/93
	  ENDIF
C
	ENDIF
C
	IF (IMIS.EQ.0) IER_KIN=9
	IF (IER_KIN.NE.0) GOTO 114
C	
C ----------------------------------------------------------------------------
C                                      fill vectors for fit
C                                       - variables : 1/p , pi/2 - theta , phi
100	DO I=1,3
	  SP_VAR(I)=1./PM_KIN(I)
	  SP_VAR(I+3)=PI/2.-TH_KIN(I)
	  SP_VAR(I+6)=FI_KIN(I)
C
	  ER_VAR(I)=SP_VAR(I)**2*EP_KIN(I)
	  ER_VAR(I+3)=ET_KIN(I)
	  ER_VAR(I+6)=EF_KIN(I)
	ENDDO
	SP_VAR(10)=PM_GAM
	ER_VAR(10)=EP_GAM
C ----------------------------------------------------------------------------
C                                                                  fit routine
	CALL CONKIN
C
	NST_KIN = NSTEP
	C4_KIN  = FTEST
	CHI_KIN = CHSQ
	CL_KIN  = PR
	IER_KIN = IEF
C
	IF (IER_KIN.NE.0) THEN
c...	  DO I=1,3
c...	    IF(EP_KIN(I).EQ.0.) PM_KIN(I)=0.
c...	    IF(ET_KIN(I).EQ.0.) TH_KIN(I)=0.
c...	    IF(EF_KIN(I).EQ.0.) FI_KIN(I)=0.
c...	  ENDDO
	  RESULT=IER_KIN
	  GOTO 114
	ENDIF
C ----------------------------------------------------------------------------
C                                                    conversion for humans:
C                                                      - 1/p        --> p
C                                                      - pi/2-theta --> theta
	DO I=1,3
	  PM_KIN(I)=1./FT_VAR(I)
	  EP_KIN(I)=PM_KIN(I)**2*ER_VAR(I)
	  TH_KIN(I)=PI/2. - FT_VAR(I+3)
	  ET_KIN(I)=ER_VAR(I+3)
	  FI_KIN(I)=FT_VAR(I+6)
	  EF_KIN(I)=ER_VAR(I+6)
	ENDDO
	PM_GAM=FT_VAR(10)
	EP_GAM=ER_VAR(10)
C ----------------------------------------------------------------------------
C                           Invariant masses (MeV) : PP_INV = M_12, M_23, M_31
C
113	DO J=1,3
	  K1=J				      
	  K2=MAX(1,MOD(J+1,4))	      
	  P1=PM_KIN(K1)
	  P2=PM_KIN(K2)
          Z1=COS(TH_KIN(K1))*COS(TH_KIN(K2))+
     *    SIN(TH_KIN(K1))*SIN(TH_KIN(K2))*COS(FI_KIN(K1)-
     *    FI_KIN(K2))
          PP_INV(J)=SQRT((SQRT(P1**2+AMASS(K1)**2)+
     *    SQRT(P2**2+AMASS(K2)**2))**2-
     *    (P1**2+P2**2+2*P1*P2*Z1))
	END DO
C
c thierry
c	call consta
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC1_1
C *
C * estimates initial values for 1 umeasured momenta
C * uses energy-momenta conservation eqs
C *
C ****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      DIMENSION P(4),E_BIL(3)
C
	DO 1 I=1,3
1	IF (EP_KIN(I).EQ.0) NOMEA=I
C
	DO 2 I=1,2
2	P(I) = 0.
	P(3) = PM_GAM
C ----------------------------------------------------------------------------
C                                           using momenta conservation gives
C                                           three differrent 
C                                           estimation of the missing momentum
	DO I=1,3
	  IF (I.NE.NOMEA) THEN
	    P(1) = P(1)-PM_KIN(I)*SIN(TH_KIN(I))*COS(FI_KIN(I))
	    P(2) = P(2)-PM_KIN(I)*SIN(TH_KIN(I))*SIN(FI_KIN(I))
	    P(3) = P(3)-PM_KIN(I)*COS(TH_KIN(I))
	  ENDIF
	ENDDO
	P(1) = P(1) / (SIN(TH_KIN(NOMEA))*COS(FI_KIN(NOMEA)))
	P(2) = P(2) / (SIN(TH_KIN(NOMEA))*SIN(FI_KIN(NOMEA)))
	P(3) = P(3) / COS(TH_KIN(NOMEA))
C ----------------------------------------------------------------------------
C                                      only from energy conservation calculate
C                                      missing momentum
C
	P2 = (E0 - E1)**2 - AMASS(NOMEA)**2
	IF (P2.GT.0) THEN
	  P(4) = SQRT(P2)
	ELSE
	  P(4) = -100000.
	ENDIF
C ----------------------------------------------------------------------------
C                                         test all solutions and take the best
	E0 = PM_GAM + AMNUC
	E1 = 0.
	DO I=1,3
	  IF (I.NE.NOMEA) E1=E1+SQRT(PM_KIN(I)**2+AMASS(I)**2)
	ENDDO
C
	DO I=1,3
	  IF(P(I).GT.0) THEN
	    E_BIL(I) = ABS(E0-E1-SQRT(P(I)**2+AMASS(NOMEA)**2))
	  ELSE
	    E_BIL(I) = 100000.
	  ENDIF
	ENDDO
	MIN_SOL = MINFZE(E_BIL(1),3)
C	
	IOK=0
	IF (P(4).LT.0) THEN
	  IF (E_BIL(MIN_SOL).LT.80.) THEN
	    PM_KIN(NOMEA) = P(MIN_SOL)
	  ELSE
	    IER_KIN = 10
	    GOTO 114
	  ENDIF
	ELSE
	  DIFF = ABS(P(MIN_SOL)-P(4))
	  IF (DIFF.LT.50) THEN
	    PM_KIN(NOMEA) = (P(MIN_SOL)+P(4))/2.
	  ELSE
	    PM_KIN(NOMEA) = P(4)
	  ENDIF
	ENDIF
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC1_2
C *
C * estimates initial values for 2 umeasured momenta
C * 
****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      DIMENSION A(2,2),B(2,1),WORK(10)
      DIMENSION A1(3,3),C(3),AM(2),PM(2,3),EBIL(3)
C
	DO I=1,3
	  EBIL(I) = 100000.
	  IF (EP_KIN(I).NE.0) MEA=I
	ENDDO
C
	C(1) = -PM_KIN(MEA)*SIN(TH_KIN(MEA))*COS(FI_KIN(MEA))
	C(2) = -PM_KIN(MEA)*SIN(TH_KIN(MEA))*SIN(FI_KIN(MEA))
	C(3) = -PM_KIN(MEA)*COS(TH_KIN(MEA)) + PM_GAM
C
	N0=0
	DO I=1,3
	  IF (I.NE.MEA) THEN
	    N0=N0+1
	    A1(1,N0)=SIN(TH_KIN(I))*COS(FI_KIN(I))
	    A1(2,N0)=SIN(TH_KIN(I))*SIN(FI_KIN(I))
	    A1(3,N0)=COS(TH_KIN(I))
	    AM(N0)=AMASS(I)
	  ENDIF
	ENDDO
C
	KSOL=0
	DO 4 I=1,2
	  DO 5 J=I+1,3
	    DO K=1,2
	      A(1,K)=A1(I,K)
	      A(2,K)=A1(J,K)
	    ENDDO
	    B(1,1)=C(I)
	    B(2,1)=C(J)
C
	    IFAIL=0
	    CALL REQN(2,A,2,WORK,IFAIL,1,B)
	    IF (IFAIL.NE.0) GOTO 5
	    DO K=1,2
	      IF (B(K,1).LT.0) GOTO 5
	    ENDDO
C
	    KSOL=KSOL+1
	    DO K=1,2
	      PM(K,KSOL)=B(K,1)
	    ENDDO
	    EBIL(KSOL) = SQRT(PM_KIN(MEA)**2 + AMASS(MEA)**2)-
     *                   (PM_GAM+AMNUC)
	    DO K=1,2
	      EBIL(KSOL)=EBIL(KSOL)+SQRT(PM(K,KSOL)**2+AM(K)**2)
	    ENDDO
C
5	  ENDDO
4	ENDDO
C
	IF (KSOL.EQ.0) THEN
	  IER_KIN=11
	  GOTO 114
	ENDIF
C
	MIN_SOL=1
	IF (KSOL.GT.1) THEN
	  DO K=2,KSOL
	    IF (ABS(EBIL(K)).LT.ABS(EBIL(MIN_SOL))) MIN_SOL=K
	  ENDDO
	ENDIF
C
	N0=0
	DO K=1,3
	  IF (EP_KIN(K).EQ.0.) THEN
	    N0=N0+1
	    PM_KIN(K)=PM(N0,MIN_SOL)
	  ENDIF
	ENDDO
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC1_3
C *
C * estimates initial values for unmeasured momenta
C * uses momenta conservation equations:
C *  i   - sum(P_xi)=0.       P_xi/Pi = sin(theta) cos(phi)
C *  ii  - sum(P_yi)=0.       P_yi/Pi = sin(theta) sin(phi)
C *  iii - sum(P_zi)=P_Y      P_zi/Pi = cos(theta)
C * Matrix notation
C *
C *  | P_x1/P1 P_x2/P2 P_x3/P3 |   |P1|   | 0 |
C *  | P_y1/P1 P_y2/P2 P_y3/P3 | x |P2| = | 0 |
C *  | P_z1/P1 P_z2/P2 P_z3/P3 |   |P3|   |P_Y|
C * 
****************************************************************************
C *
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      DIMENSION A(3,3),B(3,1),WORK(10)
C
	IFAIL=0
C
	DO 1 I=1,2
1	B(I,1)=0.
	B(3,1)=PM_GAM
	DO I=1,3
	  A(1,I)=SIN(TH_KIN(I))*COS(FI_KIN(I))
	  A(2,I)=SIN(TH_KIN(I))*SIN(FI_KIN(I))
	  A(3,I)=COS(TH_KIN(I))
	ENDDO
C
	CALL REQN(3,A,3,WORK,IFAIL,1,B)
	IF (IFAIL.NE.0) THEN
	  IER_KIN=12
	  GOTO 114
	ENDIF
C
	DO I=1,3
	  IF (B(I,1).GT.0) THEN
	    PM_KIN(I)=B(I,1)
	  ELSE
	    IER_KIN=13
	    GOTO 114
	  ENDIF
	ENDDO
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC2(IOUT)
C *
C * 31/05/94
C *
C * estimates initial values for 1 unmeasured particle
C * IOUT : unmeasured particle
C * uses momentum equations for Y + N --> (12) + 3
C * 
C ****************************************************************************
C *
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      DATA PI /3.141592653/
C
	PX=0.
	PY=0.
	PZ=0.
	DO I=1,3
	  IF (I.NE.IOUT) THEN
	    PX=PX+PM_KIN(I)*SIN(TH_KIN(I))*COS(FI_KIN(I))
	    PY=PY+PM_KIN(I)*SIN(TH_KIN(I))*SIN(FI_KIN(I))
	    PZ=PZ+PM_KIN(I)*COS(TH_KIN(I))
	  ENDIF
	ENDDO
C
	PM_PP=SQRT(PX**2+PY**2+PZ**2)
	TH_PP=ACOS(PZ/PM_PP)
	IF ((PX.EQ.0.).AND.(PY.EQ.0.)) THEN
	  FI_PP=0.
	ELSE
	  FI_PP=ATAN2(PY,PX)
	END IF
	IF (FI_PP.LT.0.) FI_PP=2.*PI+FI_PP
C
	X1=SIN(TH_PP)
	X2=PM_GAM/PM_PP-COS(TH_PP)
	TH_N=ATAN2(X1,X2)
	IF (TH_N.LT.0.) TH_N=2.*PI+TH_N
C
	PM_N=PM_PP*SIN(TH_PP)/SIN(TH_N)
	FI_N=FI_PP+PI
	IF (FI_N.GT.2.*PI) FI_N=FI_N-2.*PI
C
	PM_KIN(IOUT)=PM_N
	TH_KIN(IOUT)=TH_N
	FI_KIN(IOUT)=FI_N
C
	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC2_MERDE(IOUT)
C *
C * estimates initial values for 1 unmeasured particle
C * IOUT : unmeasured particle
C * uses 4-contrainsts equations
C *  i   - sum(P_xi)=0.       P_xi/Pi = sin(theta) cos(phi)
C *  ii  - sum(P_yi)=0.       P_yi/Pi = sin(theta) sin(phi)
C *  iii - sum(P_zi)=P_Y      P_zi/Pi = cos(theta)
C *  iv  - sum sqrt(Pi^2+Mi^2)=P_Y + M0
C *
C * ii/i      ==> tan(phi_out)
C * iv        ==> P_out
C * iii/P_out ==> cos(theta_out)
C * 
C ****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
C
	F1=0.
	F2=0.
	T1=PM_GAM
	E1=PM_GAM+AMNUC
	DO I=1,3
	  IF (I.NE.IOUT) THEN
cthierry   def completely different here of the phi angle. All the rest
C of the program has another one...
	    F1=F1-PM_KIN(I)*SIN(TH_KIN(I))*SIN(FI_KIN(I))
	    F2=F2-PM_KIN(I)*SIN(TH_KIN(I))*COS(FI_KIN(I))
	    T1=T1-PM_KIN(I)*COS(TH_KIN(I))
	    E1=E1-SQRT(PM_KIN(I)**2+AMASS(I)**2)
	  ENDIF
	ENDDO
C
	FI_KIN(IOUT)=ATG(F1,F2)
C
	E1=E1**2-AMASS(IOUT)**2
C
	IF (E1.LE.0) THEN
	  IER_KIN=14
	  GOTO 114
	ENDIF
C
	PM_KIN(IOUT)=SQRT(E1)
C
	T1=T1/PM_KIN(IOUT)
C
	IF (ABS(T1).GT.1) THEN
	  IER_KIN=15
	  GOTO 114
	ENDIF
C
	TH_KIN(IOUT)=ACOS(T1)
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE RECALC3
C *
C * from P.A.W. mod. by Ale (05.1992)
C * see Luc B. for detailes
C *
C * calculate values for 1 unmeasured particle + 1 other momentum.
C * NO fit ! No constraints = No miss.
C * part a (MOM)   : particle completily measured
C * part b (NOMOM) : particle with measured angles and unmeasured momentum
C * part c (NOTRK) : particle not measured at all
C *
C * uses energy and total momentum  conservation :
C *        _       _     _    _
C *  i   - Pc^2 = [P_Y - Pa - Pb]^2
C *  ii  - sqrt(Pc^2+Mc^2) = P_Y + M0 - sqrt(Pa^2+Ma^2) sqrt(Pa^2+Ma^2) 
C *
C ****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      REAL P3(4)
C
C	sort out which momentum has been measured, and which tracks.
	MOM   = 0
	NOMOM = 0
	NOTRK = 0
	DO I = 1,3
		IF (EP_KIN(I).GT.0.) THEN 
		MOM = I
		GOTO 10
		END IF
		END DO
10	DO I = 1,3
		IF ((ET_KIN(I).GT.0.).AND.(I.NE.MOM)) THEN
		NOMOM = I
		GOTO 20
		END IF
		END DO
20	NOTRK = 6 - MOM - NOMOM
C
	IF (NOTRK.EQ.0.OR.MOM.EQ.0.OR.NOMOM.EQ.0) THEN
	  IER_KIN=16
	  GOTO 114
	ENDIF
C
C	method using 4 angles, 1 momentum and the photon energy
	CALL SOLVE(MOM,NOMOM,NOTRK,R1,R2)
C
	IF(IER_KIN.NE.0) GOTO 114
C
	RR = MAX(R1,R2)
	PM_KIN(NOMOM) = RR
C
	P3(1) = - (PM_KIN(MOM)*SIN(TH_KIN(MOM))*COS(FI_KIN(MOM))+
     *          RR*SIN(TH_KIN(NOMOM))*COS(FI_KIN(NOMOM)))
	P3(2) = - (PM_KIN(MOM)*SIN(TH_KIN(MOM))*SIN(FI_KIN(MOM))+
     *          RR*SIN(TH_KIN(NOMOM))*SIN(FI_KIN(NOMOM)))
	P3(3) = PM_GAM - PM_KIN(MOM)*COS(TH_KIN(MOM))-
     *          RR*COS(TH_KIN(NOMOM))
	P3(4) = SQRT(P3(1)**2+P3(2)**2+P3(3)**2+AMASS(NOTRK)**2)
C
	CALL POL_CONV(P3,TH_KIN(NOTRK),FI_KIN(NOTRK),PM_KIN(NOTRK))
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE SOLVE(MOM,NOMOM,NOTRK,R1,R2)
C *
C * from P.A.W. mod by Ale (05.1992)
C *
C ****************************************************************************
C *
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      COMMON/COM_KINEM/ID_PART(3),PM_KIN(3),TH_KIN(3),FI_KIN(3),
     +      EP_KIN(3),ET_KIN(3),EF_KIN(3),PM_GAM,EP_GAM,
     +      PP_INV(3),NST_KIN,C4_KIN,CHI_KIN,CL_KIN,IER_KIN
      REAL*8 P1(4),A,B,C,F1,F2,F3,TH2,PH2,DELTA2,DELTA
C
C write eqn. as F1*E + F2*P + F3
C whence (P**2)*(F2**2-F1**2) + 2*F2*F3*P + (F3**2 - (F1**2)*M**2)
C
	TH2=DBLE(TH_KIN(NOMOM))
	PH2=DBLE(FI_KIN(NOMOM))
C
	P1(1)=DBLE(PM_KIN(MOM)*SIN(TH_KIN(MOM))*COS(FI_KIN(MOM)))
	P1(2)=DBLE(PM_KIN(MOM)*SIN(TH_KIN(MOM))*SIN(FI_KIN(MOM)))
	P1(3)=DBLE(PM_KIN(MOM)*COS(TH_KIN(MOM)))
	P1(4)=DBLE(SQRT(P1(1)**2+P1(2)**2+P1(3)**2+AMASS(MOM)**2))
C F1
	F1=DBLE(2.*(P1(4)-AMNUC-PM_GAM))
C F2
	F2=DBLE(2.*(PM_GAM*COS(TH2)-P1(1)*SIN(TH2)*COS(PH2)-
     *     P1(2)*SIN(TH2)*SIN(PH2)-P1(3)*COS(TH2)))
C F3
	F3=DBLE(2.*P1(4)*(-AMNUC-PM_GAM) + 2.*AMNUC*PM_GAM +
     *     2.*P1(3)*PM_GAM + AMASS(MOM)**2 + AMASS(NOMOM)**2 +
     *     AMNUC**2 - AMASS(NOTRK)**2)
C Quadratic root formula aX2 + bX + c
	A = (F2**2 - F1**2)
	B = 2.*F2*F3
	C = F3**2 - (F1**2)*(DBLE(AMASS(NOMOM))**2)
C
	DELTA2 = B**2 - 4.*A*C
	IF (DELTA2.GE.0.) THEN
		DELTA = SQRT(DELTA2)
	ELSE
	 	IER_KIN=17
		GOTO 114
	END IF
	R1 = (-B + DELTA)/(2.*A)
	R2 = (-B - DELTA)/(2.*A)
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE POL_CONV(P,THETA,PHI,PTOT)
C *
C * from P.A.W. mod by Ale (05.1992)
C * convert cartesian 3-momentum P(3) to polar coords PTOT,THETA,PHI
C *
C ****************************************************************************
C *
      REAL P(4)
      DATA PI /3.141592653/
C
	PTOT  = SQRT(P(1)**2+P(2)**2+P(3)**2)
	THETA = ACOS(P(3)/PTOT)
	IF ((P(1).EQ.0.).AND.(P(2).EQ.0.)) THEN
		PHI = 0.
	ELSE
		PHI = ATAN2(P(2),P(1))
	END IF
	IF (PHI.LT.0.) PHI = 2.*PI + PHI
C
114	RETURN
	END
C ****************************************************************************
	SUBROUTINE CONKIN
C *
C * Interface between kinfit and conles package
C *
C * --------- common block
C *      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
C * --------- input
C * SP_VAR  experimental values
C *         1/p_1,1/p_2,1/p_3,lambda_1,lambda_2,lambda_3,phi_1,phi_2,phi_3,p_Y
C *         (momentum in MeV/c , angles in rad , energies in MeV)
C * ER_VAR  experimental errors
C * --------- output
C * ER_VAR  errors
C * FT_VAR  fitted values
C *
C * --------- notes
C *           parameters used in fit:
C *             - 1/p
C *             - lambda = pi/2 - theta
C *             - phi
C *           units:
C *             - MeV
C *             - MeV/c
C *             - rad
C *
C ****************************************************************************
C *
C
      COMMON/INKIN/SP_VAR(10),ER_VAR(10),AMASS(3),AMNUC,FT_VAR(10)
      REAL*8 V(1275),DELY(50),F(20),B(1000),PL(50),CHSQ,PR
      COMMON/CONF1/V,DELY,F,B,PL,
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      REAL*8 DET,FTEST,FTESTL,CHSQV,DELYP(50),CST(4)
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP,
     1   CST,NSTEP,NCUT,NCST,NST(3)
C
C
      REAL*8 Y(50),YF(50),ETEMP(4),CC(4)
      DATA PI /3.141592653/
C
C ----------------------------------------------------------------------------
C                                                    solo per inizializzazione
	ICASE=1
	IDIAG=2
	CC(1)=0.1E-5
	CC(2)=0.01
	CC(3)=0.1
	CC(4)=1.0
	CALL CONDEF(ICASE,10,10,5,CC(1),CC(2),CC(3),CC(4))
	CALL CONDIA(ICASE,2,2)
cthierry modif tempo	CALL CONDIA(ICASE,10,2)
	CALL CONPAR(ICASE)
C ----------------------------------------------------------------------------
C                                         num variabili e num constraint
C                                         The covariance matrix is set to zero
	NY=10
	NF=4
C
c...	CALL VZERO(F,20)
c...	CALL VZERO(B,1000)
c...	CALL VZERO(DELY,50)
c...	CALL VZERO(PL,50)
	CALL CONLES(1)
C ----------------------------------------------------------------------------
C                                           matrice errori in forma vettoriale
	K=0
	DO I=1,NY
	  K=K+I
	  V(K)=ER_VAR(I)**2
	ENDDO
C
	CALL CONLES(2)
	DO I=1,NY
	  Y(I)=SP_VAR(I)
	ENDDO
C ----------------------------------------------------------------------------
C                                                            MAIN FITTING LOOP
	IBR = 0
	DO WHILE (IBR.LE.1)
C ----------------------------------------------------------------------------
C                           The correction DELY(I) for each variable is added
11	DO 12 I=1,NY
12	YF(I) = Y(I) + DELY(I)
C ----------------------------------------------------------------------------
C                                         calcolo valori F(i) dei constrainsts
	ET = 0.
	PX = 0.
	PY = 0.
	PZ = 0.
	DO I=1,3
	  PM = 1./(YF(I))
	  ETEMP(I)=SQRT(PM**2+AMASS(I)**2)
	  PX = PM*COS(YF(I+3))*COS(YF(I+6)) + PX
	  PY = PM*COS(YF(I+3))*SIN(YF(I+6)) + PY
	  PZ = PM*SIN(YF(I+3)) + PZ
	  ET = ET + ETEMP(I)
	END DO
	F(1) = PX 
	F(2) = PY 
	F(3) = PZ - YF(10)
	F(4) = ET - YF(10) - AMNUC
C ----------------------------------------------------------------------------
C                                              Check of convergence conditions
	CALL CONCHK(IBR)
	IF(IBR.EQ.0) THEN
C ----------------------------------------------------------------------------
C                                Derivatives of the constraint functions 
C                                with respect to variables Y(I) are calculated
c...	CALL VZERO(B,NY*NF)
	DO I=1,3
	  I1=0
	  I2=NY
	  I3=2*NY
	  I4=3*NY
	  PM = 1./(YF(I))
	  AL = YF(I+3)
	  AF = YF(I+6)
	  B(I1+I)   = -PM**2*COS(AL)*COS(AF)
	  B(I1+I+3) = -PM*SIN(AL)*COS(AF)
	  B(I1+I+6) = -PM*COS(AL)*SIN(AF)
	  B(I1+10)  = 0.
C
	  B(I2+I)   = -PM**2*COS(AL)*SIN(AF)
	  B(I2+I+3) = -PM*SIN(AL)*SIN(AF)
	  B(I2+I+6) =  PM*COS(AL)*COS(AF)
	  B(I2+10)  = 0.
C
	  B(I3+I)   = -PM**2*SIN(AL)
	  B(I3+I+3) = PM*COS(AL)
	  B(I3+I+6) = 0.
	  B(I3+10)  = -1.
C
	  B(I4+I)   = -PM**3/ETEMP(I)
	  B(I4+I+3) = 0.
	  B(I4+I+6) = 0.
	  B(I4+10)  = -1.
	ENDDO
C ----------------------------------------------------------------------------
C                            One step of least square calculation is performed
	CALL CONLES(3)
	END IF
	END DO
C
	IF(IBR.NE.2) GOTO 114
C ----------------------------------------------------------------------------
C                          Convergence reached.The output variables are filled
	CALL CONLES(4)
	DO I=1,NY
	  FT_VAR(I)=YF(I)
	ENDDO
C ----------------------------------------------------------------------------
C                              errori dagli elementi diagonali della matrice V
	K=0
	DO I=1,NY
	  K=K+I
	  ER_VAR(I)=SQRT(ABS(V(K)))
	ENDDO
C
114	RETURN
	END

      SUBROUTINE CONCOM
C.        THIS PACKAGE OF SUBROUTINES CONTAINS ALL THE STEPS
C.        NECESSARY FOR LEAST SQUARES FITS WITH CONSTRAINT EQUATIONS.
C.        THE USER HAS COMPLETE CONTROL OVER THE STEPS TO EXECUTE,
C.        ANALYTICAL/MUMERICAL DIFFERENTIATION, CONVERGENCE CRITERIA ETC
C.        MEASURED VALUES AND THEIR ERROR MATRIX, STARTING VALUES
C.        FOR UNKNOWN VARIABLES, AND THE CONSTRAINT EQUATIONS
C.        HAVE TO BE SUPPLIED.
C.            DOCUMENTATION = IN LINE
C.            ORIGIN = V. BLOBEL '80, STATUS = FIELD PROVEN (DESY)
C
C
C     PROGRAM PACKAGE FOR CONSTRAINED LEAST SQUARE FIT
C     ------------------------------------------------
C
C     SUBROUTINES
C     CONPAR   (ADD. ENTRIES  CONDEF  CONDIA )
C     CONLES  CONCHK  CONPRI  CONPMT  CONSTA  CONDER
C
C     ADDITIONAL SUBROUTINES
C     SMINV   PROB
C     THE PACKAGE DESCRIBED BELOW IS A RATHER GENERAL SET OF
C     SUBROUTINES TO PERFORM CONSTRAINED LEAST SQUARE FITS.
C     FOR EACH APPLICATION THE USER HAS TO WRITE A CALLING
C     PROGRAMM.
C
C     THE SET OF SUBROUTINES IS WRITTEN BY
C
C             VOLKER BLOBEL
C             II. INSTITUT FUER EXPERIMENTALPHYSIK
C             UNIVERSITY OF HAMBURG
C
C     THE CONLES PACKAGE IS USED IN THE ANALYSIS PROGRAMS FOR
C     THE PLUTO DETECTOR AT DESY FOR THE FOLLOWING PROBLEMS
C
C     1.-  KINEMATICAL FIT
C
C     INPUT ARE THE TRACK PARAMETERS 1/P, LAMBDA AND PHI AND
C     SOME HYPOTHESIS ON THE MASS ASSIGNMENT. CONSTRAINTS ARE
C     THE 4 ENERGY AND MOMENTUM EQUATIONS.
C
C
C     2.-  VERTEX FIT
C
C     HAVING A WELL DEFINED INTERACTION POINT, THE PARAMETERS
C     OF TRACKS ORIGINATING FROM THIS POINT ARE MODIFIED BY
C     THE CONSTRAINT, THAT THEY HAVE TO EXTRAPOLATE TO THE
C     INTERACTION POINT.
C
C
C     3.-  K ZERO FINDING
C
C     FOR SELECTED PAIRS OF (+-) TRACKS, CANDIDATES FOR K ZEROS
C     ARE FOUND BY FINDING A DECAY POINT IN SPACE, CONSTRAINING
C     THE TWO TRACKS TO PASS THROUGH THIS INTERACTION POINT,
C     WITH TRANSVERSE MOMENTUM BALANCE IN THE TWO DIRECTIONS,
C     ORTHOGONAL TO THE LINE FROM THE INTERACTION POINT TO THE
C     DECAY POINT.
C     THERE ARE 13 MEASURED PARAMETERS (2*5 TRACK PARAMETERS
C     AND 3 COORDINATES OF THE INTERACTION POINT),
C     4 UNMEASURED PARAMETERS (K ZERO MOMENTUM AND 3 COORDI-
C     NATES OF THE DECAY POINT) AND 7 CONSTRAINT EQUATIONS
C     (4 GEOMETRICAL CONSTRAINTS AND 3 KINEMATICAL CONSTRAINTS.)
C
C
C
C
C     ********************   USAGE   ***********************
C
C
C     USER COMMON /CONF1/
C
C     COMMON/CONF1/NY,NF,V(1275),DELY(50),F(20),B(1000),PL(50),
C    1   IEF,ND,CHSQ,PR,IDIAG,IUNPH
C
C     NY      = NR OF VARIABLES NY (MAX. 50)
C     NF      = NR OF CONSTRAINT EQUATIONS (MAX. 20)
C     V( )    = COVARIANCE MATRIX OF Y( )
C     DELY( ) = CORRECTIONS TO ORIGINAL VALUES Y( )
C     F( )    = VALUES OF CONSTRAINT EQUATIONS
C     B( )    = NY-BY-NF MATRIX OF DERIVATIVES
C     PL( )   = PULLS FOR VARIABLES Y( )
C     IEF     = ERROR CODE
C     ND      = NR OF DEGREES OF FREEDOM
C     CHSQ    = CHI SQUARE VALUE
C     PR      = PROBABILITY
C     IDIAG   = DIAGNOSTIC MARKER
C     IUNPH   = MARKER FOR UNPHYSICAL VALUES
C
C
C     STEP
C        1   CALL CONDEF(ICASE,NST1,NST2,NST3,CST1,CST2,CST3,CST4)
C        2   CALL CONDIA(ICASE,NPRT,NPRA)
C
C        3   CALL CONPAR(ICASE)
C        4   NY= NR OF VARIABLES Y(I)
C        5   NF= NR OF CONSTRAINTS F(J)
C        6   CALL CONLES(1)
C        7   V( )= . . .   DEFINE COVARIANCE MATRIX V( )
C        8   CALL CONLES(2)
C        9   Y(I)= . . .   ESTIMATE INITIAL VALUES FOR UNMEASURED Y(I)
C       10   CALL CONDES(I,STP)   DEFINE STEP SIZES FOR UNMEASURED Y(I)
C       11   YF(I)=Y(I)+DELY(I)   ADD CORRECTIONS TO ORIGINAL Y(I)
C       12   F(J)= . . .   COMPUTE CONSTRAINTS WITH CORRECTED YF(I)
C       13   CALL CONDER(IREP)
C       14   IF(IREP.NE.0) GOTO 11
C       15   CALL CONCHK(IBR)
C       16   IF(IBR.NE.0) GOTO (11,20,21),IBR
C       17   B( )= . . .    COMPUTE MATRIX OF DERIVATIVES
C       18   CALL CONLES(3)
C       19   GOTO 11
C       20   CALL CONLES(4)   END OF SUCCESSFUL FIT
C       21   CONTINUE         END OF UNSUCCESSSFUL FIT
C
C       22   CALL CONSTA      PRINT FINAL STATISTICS
C
C
C     EXPLANATION
C     -----------
C
C        1   CALL CONDEF(ICASE,NST1,NST2,NST3,CST1,CST2,CST3,CST4)
C        2   CALL CONDIA(ICASE,NPRT,NPRA)
C
C     THESE STEPS ARE OPTIONAL. IF USED, THEY ARE EXECUTED ONCE
C     AT THE BEGINNING OF THE PROGRAM. EACH DIFFERENT APPLICATION
C     OF THE PROGRAM CAN BE GIVEN A CASE-NUMBER ICASE (STEP 3), WHICH
C     DETERMINES THE SET OF CONSTANTS USED IN CONCHK FOR THE
C     RECOGNITION OF CONVERGENCE. THE DEFAULT VALUES OF THE CONSTANTS
C     ARE THE FOLLOWING (COMPARE STEP 15).
C            NST1 = 10                     CST1 = 0.1E-3
C            NST2 =  4                     CST2 = 0.01
C            NST3 =  4                     CST3 = 0.1
C                                          CST4 = 1.0
C     BY CALLING CONDIA SHORT DIAGNOSTIC (IDIAG=1) IS PRINTED FOR
C     THE FIRST NPRT TIMES AND LONG DIAGNOSTIC (IDIAG=2) FOR THE
C     FIRST NPRA TIMES FOR CASE-NUMBER ICASE.
C
C
C        3   CALL CONPAR(ICASE)
C
C     THE CASE-NUMBER ICASE IS ASSIGNED TO THE CURRENT APPLICATION.
C     ALLOWED VALUES ARE 1 TO 10. THE SET OF CONSTANTS USED IN
C     CONCHK IS DEFINED.
C
C
C        4   NY= NR OF VARIABLES Y(I)
C        5   NF= NR OF CONSTRAINTS F(J)
C        6   CALL CONLES(1)
C
C     AFTER DEFINITION OF THE NUMBERS OF VARIABLES AND OF CONSTRAINTS,
C     THE COVARIANCE MATRIX IS SET TO ZERO BY CALLING CONLES(1).
C
C
C        7   V( )= . . .   DEFINE COVARIANCE MATRIX V( )
C        8   CALL CONLES(2)
C
C     THE USER HAS TO STORE THE COVARIANCE MATRIX OF THE Y(I) IN
C     MATRIX V( ) IN SYMMETRIC STORAGE MODE. ELEMENTS OF V( )
C     REFERRING TO UNMEASURED VARIABLES Y(I) HAVE TO BE ZERO.
C     SUBROUTINE SSVW CAN BE USED TO STORE SUBMATRICES INTO THE
C     MATRIX V( ). ASSUME THAT THE COVARIANCE MATRIX FOR THE
C     VARIABLES Y(I),Y(I+1) . . . Y(I+N-1) IS STORED IN A
C     SYMMETRIC N-BY-N MATRIX W( ). THIS SUBMATRIX CAN BE STORED
C     IN THE SYMMETRIC MATRIX V( ) AT ROWS/COLS I, I+1 . . . I+N-1
C     BY THE FOLLOWING CALL.
C                  - - - -
C        CALL SSVW(V,I,W,N)
C                  -
C     THE NUMBER OF DEGREES OF FREEDOM FOR THE FIT IS NF MINUS
C     RANK DEFECT OF THE MATRIX V( ). THE RANK DEFECT OF THE MATRIX V( )
C     IS EQUAL TO THE NUMBER OF ZERO DIAGONAL ELEMENTS OF V( ).
C     THE CORRECTIONS DELY(I) AND OTHER PARAMETERS ARE SET ZERO.
C
C
C        9   Y(I)= . . .   ESTIMATE INITIAL VALUES FOR UNMEASURED Y(I)
C
C     FOR UNMEASURED VARIABLES Y(I) REASONABLE GOOD INITIAL VALUES
C     HAVE TO BE ESTIMATED. THIS IS USUALLY DONE BY USING
C     SELECTED CONSTRAINTS.
C
C
C       10   CALL CONDES(I,STP)   DEFINE STEP SIZES FOR UNMEASURED Y(I)
C
C     STEPS 10 AND 13,14 ARE OPTIONAL. IF USED, STEP 17 IS NOT
C     NECESSARY. THEY ARE USED TO OBTAIN THE MATIRX B( ) OF DERIVATIVES
C     BY NUMERICAL METHODS. BECAUSE STEP SIZES ARE NECESSARY, THEY
C     HAVE TO BE DEFINED BY THE USER FOR UNMEASURED VARIABLES. THE
C     ORDER OF MAGNITUDE OF THE STEP SIZES SHOULD BE ABOUT ONE
C     (EXPECTED) STANDARD DEVIATION OF THE CORRESPONDING VARIABLE.
C     FOR MEASURED VARIABLES THE STEP SIZES ARE TAKEN FROM THE
C     COVARIANCE MATRIX.
C     FOR IDIAG=2, THE CALCULATED MATRIX IS PRINTED. IF THE MATRIX
C     B IS CALCULATED ALSO AFTER CONCHK ANALYTICALLY, THE PRINTOUT
C     MAY BE COMPARED, TO FIND ERRORS IN THE ANALYTIC DERIVATIVES.
C
C
C       11   YF(I)=Y(I)+DELY(I)   ADD CORRECTIONS TO ORIGINAL Y(I)
C
C     THE CORRECTION DELY(I) FOR EACH VARIABLE IS COMPUTED BY
C     CONLES(3). THIS HAS TO BE ADDED TO THE MEASURED OR INITIAL
C     VALUES Y(I) TO ABTAIN THE CORRECTED VALUES YF(I). Y(I) AND
C     YF(I) HAVE TO BE DIFFERENT ARRAYS.
C     IF SOME PARAMETER YF( ) HAVE AN UNPHYSICAL VALUE, MARKER IUNPH
C     SHOULD BE SET TO 1, GO TO STEP 15.
C
C
C       12   F(J)= . . .   COMPUTE CONSTRAINTS WITH CORRECTED YF(I)
C
C     THE VALUES OF THE CONSTRAINT EQUATIONS F(J) HAVE TO BE COMPUTED
C     USING THE CORRECTED VALUES YF(I).
C
C
C       13   CALL CONDER(IREP)
C       14   IF(IREP.NE.0) GOTO 11
C
C     SEE REMARK AT STEP 10. IN CONDER CORRECTIONS DELY(I) ARE MODIFIED
C     FOR ONE PARAMETER AT A TIME, THE STEPS 11 AND 12 HAVE TO BE
C     EXECUTED AGAIN, UNTIL CONDER RETURNS WITH IREP = 0. THEN THE
C     MATRIX B( ) IS DETERMINED, THE CORRECTIONS HAVE THE ORIGINAL
C     VALUES AGAIN.
C
C
C       15   CALL CONCHK(IBR)
C       16   IF(IBR.NE.0) GOTO (11,20,21),IBR
C
C     THIS ROUTINE DETERMINES CONVERGENCE OR NON-CONVERGENCE. FOR
C     IBR = 0, FURTHER NORMAL ITERATIONS ARE NECESSARY. SPECIAL
C     CONDITIONS AND CONVERGENCE IS SIGNALLED BY THE ARGUMENT
C     IBR AND BY THE ERROR CODE IEF.
C     IBR = 0   FURTHER ITERATIONS ARE NECCESSARY, GO TO STEP 17,18.
C     IBR = 1   A SO CALLED CUTSTEP IS MADE, THE CORRECTIONS ARE
C               HALF THE PREVIOUS ONES, GO BACK TO STEP 11.
C               TAKEN, IF THE CORRECTED VARIABLES HAVE UNPHYSICAL
C               VALUES (IUNPH=1) OR IF THE SUM OF THE ABSOLUTE
C               VALUES OF THE CONSTRAINTS IS LARGER THAN IN THE
C               PREVIOUS STEP.
C     IBR = 2   CONVERGENCE REACHED, GO TO STEP 20. TAKEN, IF
C               SUM OF ABSOLUTE VALUES OF THE CONSTRAINTS LESS
C               THAN CST(2) AND LAST CHISQUARE DIFFERENCE LESS
C               THAN CST(3).
C     IBR = 3   FIT UNSUCCESSFULL, GO TO STEP 21. THE REASON FOR
C               NON-CONVERGENCE IS GIVEN BY THE ERROR CODE IEF.
C     IEF = 1   PROBABILITY LESS THAN CST(1) (EXCEPT FIRST ITERATION)
C     IEF = 2   NR OF ITERATIONS GREATER THAN NST(1).
C     IEF = 3   TOTAL NR OF CUTSTEPS GREATER THAN NST(2).
C     IEF = 4   MATRIX SINGULAR.
C     IEF = 5   UNPHYSICAL VALUES AT FIRST ITERATION.
C     IEF = 6   NY OR NF TOO LARGE.
C
C
C       17   B( )= . . .    COMPUTE MATRIX OF DERIVATIVES
C     THE NY-BY-NF MATRIX B( ) OF DERIVATIVES OF THE CONSTRAINT
C     FUNCTIONS F(J) W.R.T. THE VARIABLES Y(I) HAS TO BE CALCULATED.
C     THIS STEP IS NOT NECCESSARY, IF STEPS 10,13,14 (NUMERICAL
C     DIFFERENTIATION) ARE USED.
C     THE ORDER OF THE MATRIX ELEMENTS IS AS FOLLOWS.
C     B(1)    = DERIVATIVE OF F(1) WRT Y(1)
C     B(2)    = DERIVATIVE OF F(1) WRT Y(2)
C     ...
C     B(NF)   = DERIVATIVE OF F(1) WRT Y(NF)
C     B(NF+1) = DERIVATIVE OF F(2) WRT Y(1)
C     ...
C
C
C       18   CALL CONLES(3)
C       19   GOTO 11
C
C     ONE STEP OF LEAST SQUARE CALCULATION IS PERFORMED. THE
C     CORRECTIONS DELY(I) ARE CALCULATED FROM THE VALUES OF THE
C     CONSTRAINT FUNCTIONS F(J) USING THE COVARIANCE MATRIX V( )
C     (THE INVERSE OF V IS G) AND THE MATRIX OF DERIVATIVES B( ).
C
C                               -1
C        (  DELY  )   ( G  B T )   (        0               )
C        (        ) = (        ) * (                        )
C        ( LAMBDA )   ( B   0  )   ( -F + B * PREVIOUS DELY )
C
C     THE METHOD USED FOR THE SOLUTION OF THE EQUATION ABOVE
C     REALLY DOES NOT WORK WITH THE INVERSE MATRIX G, BUT WORKS
C     IN A SPECIAL WAY WITH PARTITIONED MATRICES.
C     THE CHISQUARE, DEFINED BY
C                         T
C              CHSQ = DELY * G * DELY
C     IS CALCULATED USING THE SIMPLER FORMULA
C                            T
C              CHSQ = -LAMBDA * (-F+B*PREVIOUS DELY)
C
C
C       20   CALL CONLES(4)   END OF SUCCESSSFUL FIT
C
C     THE COVARIANCE MATRIX OF FITTED VARIABLES YF(I) IS STORED
C     IN V( ). PULLS, DEFINED BY PULL(I) = DELY(I) DIVIDED BY THE
C     SQRT OF THE DIFFERENCE BETWEEN THE CORRESPONDING V-MATRIX-
C     ELEMENTS BEFORE AND AFTER THE FIT, ARE STORED IN PL(I).
C     FOR NORMAL DISTRIBUTED MEASUREMENTS THE PULLS ARE EACH
C     DISTRIBUTED ACCORDING TO THE NORMAL DISTRIBUTION WITH
C     MEAN ZERO AND VARIANCE ONE.
C
C       21   CONTINUE        END OF UNSUCCESSSFUL FIT
C
C       22   CALL CONSTA     PRINT FINAL STATISTICS
C
C     STATISTIC ON USED CONSTANTS, TOTAL NUMBER OF APPLICATIONS,
C     USED NUMBERS OF ITERATIONS AND THE DISTRIBUTION OF THE
C     CHISQUARE PROBABILITY IS PRINTED FOR EACH CASE.
C
C
       END
      SUBROUTINE CONPAR(ICASE)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF1/V(1275),DELY(50),F(20),B(1000),PL(50),
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP(50),
     1   CST(4),NSTEP,NCUT,NCST,NST(3)
      COMMON/CONF4/CSTD(10,4),NPROB(10,10),LCASE,NSTD(10,3),ISTD(10,2),
     1   NCASE(10,8),NITER(10,10)
      INTEGER NSTF(3)/10,4,4/
C
C LUC CHANGE
C
c      REAL*8 CSTF(4)/.1E-3,.01,.1,1.0/
      REAL*8 CSTF(4)/.1,.01,.1,1.0/
      INTEGER IG/0/
      LOGICAL OUT
      OUT(I)=I.LT.1.OR.I.GT.10
      JCASE=ICASE
      IF(OUT(JCASE)) JCASE=10
      IF(IG.EQ.0) GOTO 81
   10 DO 12 J=1,3
      NST(J)=NSTD(JCASE,J)
   12 CST(J)=CSTD(JCASE,J)
      CST(4)=CSTD(JCASE,4)
      LCASE=JCASE
      IDIAG=0
      NCASE(JCASE,1)=NCASE(JCASE,1)+1
      IF(ISTD(JCASE,1).GE.NCASE(JCASE,1)) IDIAG=1
      IF(ISTD(JCASE,2).GE.NCASE(JCASE,1)) IDIAG=2
      GOTO 100
C
      ENTRY CONDEF(ICASE,NST1,NST2,NST3,CST1,CST2,CST3,CST4)
      JCASE=ICASE
      IF(OUT(ICASE)) JCASE=10
      IF(IG.EQ.0) GOTO 82
   20 NSTD(JCASE,1)=NST1
      NSTD(JCASE,2)=NST2
      NSTD(JCASE,3)=NST3
      CSTD(JCASE,1)=CST1
      CSTD(JCASE,2)=CST2
      CSTD(JCASE,3)=CST3
      CSTD(JCASE,4)=CST4
      GOTO 100
C
      ENTRY CONDIA(ICASE,NPRT,NPRA)
      JCASE=ICASE
      IF(OUT(ICASE)) JCASE=10
      IF(IG.EQ.0) GOTO 83
   30 ISTD(JCASE,1)=NPRT
      ISTD(JCASE,2)=NPRA
      GOTO 100
C
   83 IG=IG+1
   82 IG=IG+1
   81 IG=IG+1
      DO 98 I=1,10
      DO 90 J=1,3
      NSTD(I,J)=NSTF(J)
   90 CSTD(I,J)=CSTF(J)
      CSTD(I,4)=CSTF(4)
      DO 92 J=1,2
   92 ISTD(I,J)=0
      DO 94 J=1,8
   94 NCASE(I,J)=0
      DO 96 J=1,10
      NITER(I,J)=0
   96 NPROB(I,J)=0
   98 CONTINUE
      GOTO (10,20,30),IG
  100 RETURN
      END
      SUBROUTINE CONSTA
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF4/CSTD(10,4),NPROB(10,10),LCASE,NSTD(10,3),ISTD(10,2),
     1   NCASE(10,8),NITER(10,10)
      INTEGER NRS(10)/1,2,3,4,5,6,7,8,9,10/
      CHARACTER*8 TNST(3),TCST(4),TCASE(8),TITER(10),TPROB(10)
      CHARACTER*8 SCONS(3),SCASE(3),SITER(3),SPROB(3)
      DATA TNST/'  NST(1)','  NST(2)','  NST(3)'/,TCST/'  CST(1)',
     1   '  CST(2)','  CST(3)','  CST(4)'/,TCASE/'   TOTAL','NO ERROR',
     2   ' IEF = 1',' IEF = 2',' IEF = 3',' IEF = 4',' IEF = 5',
     3   ' IEF = 6'/
      DATA TITER/'  IT = 1','       2','       3','       4',
     1   '       5','       6','       7','       8','       9',
     2   '  GE  10'/
      DATA TPROB/' 0 -  10','10 -  20','20 -  30','30 -  40',
     1   '40 -  50','50 -  60','60 -  70','70 -  80','80 -  90',
     2   '90 - 100'/
      DATA SCONS/'CONS','TANT','S   '/,SCASE/3*'    '/,
     1   SITER/'ITER','ATIO','NS  '/,SPROB/'PROB',' IN ','P.C.'/
      IF(LCASE.LT.1.OR.LCASE.GT.10) GOTO 100
      WRITE(6,101)
      WRITE(6,102) NRS,SCONS
      DO 10 J=1,3
   10 WRITE(6,103) TNST(J),(NSTD(I,J),I=1,10)
      DO 20 J=1,4
   20 WRITE(6,104) TCST(J),(CSTD(I,J),I=1,10)
      WRITE(6,106) (TCASE(J),J=3,8)
      WRITE(6,102) NRS,SCASE
      DO 30 J=1,8
   30 WRITE(6,103) TCASE(J),(NCASE(I,J),I=1,10)
      WRITE(6,105)
      WRITE(6,102) NRS,SITER
      DO 40 J=1,10
   40 WRITE(6,103) TITER(J),(NITER(I,J),I=1,10)
      WRITE(6,102) NRS,SPROB
      DO 50 J=1,10
   50 WRITE(6,103) TPROB(J),(NPROB(I,J),I=1,10)
  100 RETURN
  101 FORMAT('0CONLES - STATISTICS'/1X,4('----'),'---')
  102 FORMAT('0',12X,'ICASE =',I4,9I11/1X,3A4)
  103 FORMAT(1X,A8,4X,10I11)
  104 FORMAT(1X,A8,4X,10E11.3)
  105 FORMAT('0',12X,'THE FOLLOWING REFERS TO SUCCESSFUL FITS ONLY')
  106 FORMAT('0',A8,10X,'PROB LT CST(1)'/
     1        1X,A8,10X,'NR OF ITERATIONS GT NST(1)'/
     2        1X,A8,10X,'NR OF CUTSTEPS GT NST(2)'/
     3        1X,A8,10X,'MATRIX SINGULAR'/
     4        1X,A8,10X,'UNPHYSICAL PARAMETER AT FIRST ITERATION'/
     5        1X,A8,10X,'NY OR NF TOO LARGE'/)
      END

      SUBROUTINE CONLES(IGOTO)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF1/V(1275),DELY(50),F(20),B(1000),PL(50),
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP(50),
     1   CST(4),NSTEP,NCUT,NCST,NST(3)
      COMMON/CONF3/W(2485),R(70),DG(50)
C
      GOTO (10,20,30,70),IGOTO
C
C     THE COVARIANCE MATRIX V IS SET ZERO.
   10 IF(NY.GT.50.OR.NF.GT.20) GOTO 100
      NS=NF+NY
      KY=(NY*NY+NY)/2
      DO 12 I=1,KY
   12 V(I)=0.0
      GOTO 100
C     THE NUMBER OF DEGREES OF FREEDOM (ND) IS TAKEN FROM THE
C     RANK OF THE MATRIX V
   20 IF(NY.GT.50.OR.NF.GT.20) GOTO 100
      DET=0.0
      K=0
C     ND=NF-RANK DEFECT OF V
C     STORE DIAGONAL OF V IN DG
      ND=NF
      DO 26 I=1,NY
      DELY(I)=0.0
      PL(I)=0.0
      K=K+I
      DG(I)=0.0
      IF(V(K).LE.0.0) GOTO 22
      DG(I)=V(K)
      GOTO 26
   22 ND=ND-1
      IJ=(I*I-I)/2
      DO 24 J=1,NY
      IF(J.LE.I) IJ=IJ+1
      V(IJ)=0.0
      IF(J.GE.I) IJ=IJ+J
   24 CONTINUE
   26 CONTINUE
C     NSTEP, NCUT, NCST, CHSQ , IUNPH ARE SET SEZO
      NSTEP=0
      NCUT=0
      NCST=0
      CHSQ=0.0
      IUNPH=0
      DO 28 I=1,NF
   28 F(I)=0.0
      GOTO 100
C
C     ONE STEP OF LEAST SQUARE CALCULATION IS PERFORMED
   30 IF(NSTEP.NE.0.OR.IDIAG.EQ.0) GOTO 40
C     DIAGNOSTIC OUTPUT
      WRITE(6,101) NY,NF,ND
      IF(IDIAG.NE.2) GOTO 34
      CALL CONPRI
      WRITE(6,102)
      KA=1
      DO 32 J=1,NF
      KB=KA+NY-1
      WRITE(6,103) J,(B(K),K=KA,KB)
   32 KA=KB+1
   34 WRITE(6,104) FTEST,(F(J),J=1,NF)
      WRITE(6,105)
   40 NSTEP=NSTEP+1
      CHSQV=CHSQ
C     BUILD MATRIX
      do 41 i=1,2485
41    w(i)=0.
      DO 42 I=1,KY
   42 W(I)=-V(I)
      L=0
      K=KY
      DO 48 J=1,NF
      DO 44 I=1,NY
      K=K+1
      L=L+1
   44 W(K)=B(L)
      DO 46 I=1,J
      K=K+1
   46 W(K)=0.0
   48 CONTINUE
C     BUILD VECTOR
      DO 50 I=1,NY
   50 R(I)=0.0
      K=NY
      DO 52 I=1,NF
      K=K+1
   52 R(K)=-F(I)
      K=0
      DO 54 I=1,NF
      DO 54 J=1,NY
      K=K+1
   54 R(NY+I)=R(NY+I)+B(K)*DELY(J)
      DO 56 I=1,NF
   56 F(I)=R(NY+I)
C     SOLVE SYSTEM
      CALL CONPMT(W,R,NY,NF,DET)
      IF(DET.EQ.0.0) GOTO 62
      DO 58 I=1,NY
      DELYP(I)=DELY(I)
   58 DELY(I)=R(I)
C     CALCULATE CHISQUARE
      CHSQ=0.0
      DO 60 I=1,NF
   60 CHSQ=CHSQ-F(I)*R(NY+I)
      GOTO 100
C     FOR SINGULAR MATRIX SET DELY = 0
   62 DO 64 I=1,NY
   64 DELY(I)=0.0
      GOTO 100
   70 DO 72 I=1,KY
   72 V(I)=W(I)
      K=0
      DO 74 I=1,NY
      PL(I)=0.0
      K=K+I
      IF(DG(I).EQ.0.0) GOTO 74
      SQ=DG(I)-V(K)
      IF(SQ.LE.0.0) GOTO 74
      PL(I)=DELY(I)/SQRT(SQ)
   74 CONTINUE
      IF(IDIAG.EQ.2) CALL CONPRI
      IEF=0
C
  100 RETURN
  101 FORMAT('0',7('---')/' DIAGNOSTIC OUTPUT  -  NY,NF,ND = '
     1 ,I2,2(',',I2))
  102 FORMAT('0DERIVATIVE MATRIX B'/1X,'F(K) / VAR.')
  103 FORMAT(I4,6X,10G12.4/(10X,10G12.4))
  104 FORMAT('0ITERATION',18X,G9.2,3X,10G9.2/(40X,10G9.2))
  105 FORMAT(' NSTEP  NCUT  DET',7X,'CHSQ  FTEST',7X,'F(K) = . . .')
      END

      SUBROUTINE CONPRI
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF1/V(1275),DELY(50),F(20),B(1000),PL(50),
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      REAL*8 ER(50),RS(50)
C     PURPOSE
C        PRINT DIAGNOSTIC OUTPUT (CALLED BY CONLES)
C
      WRITE(6,101)
      K=0
      DO 10 I=1,NY
      K=K+I
      ER(I)=0.0
      IF(V(K).LE.0.0) GOTO 10
      ER(I)=SQRT(V(K))
   10 CONTINUE
      K=0
      DO 40 I=1,NY
      K=K+I
      IF(ER(I).EQ.0.0) GOTO 30
      DO 20 J=1,I
      RS(J)=0.0
      IF(ER(J).NE.0.0) RS(J)=V(K-I+J)/(ER(I)*ER(J))
   20 CONTINUE
      WRITE(6,102) I,PL(I),ER(I),(RS(J),J=1,I)
      GOTO 40
   30 WRITE(6,103) I
   40 CONTINUE
  100 RETURN
  101 FORMAT('0VAR.  PULL   ERROR',6X,'CORRELATION MATRIX')
  102 FORMAT(I3,F8.2,G12.4,1X,15F6.2/(24X,15F6.2))
  103 FORMAT(I3, 6X,'UNMEASURED')
      END
      SUBROUTINE CONCHK(IBR)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF1/V(1275),DELY(50),F(20),B(1000),PL(50),
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP(50),
     1   CST(4),NSTEP,NCUT,NCST,NST(3)
      COMMON/CONF4/CSTD(10,4),NPROB(10,10),LCASE,NSTD(10,3),ISTD(10,2),
     1   NCASE(10,8),NITER(10,10)
c - ale
      EXTERNAL PROB
      REAL*4 VPR,CPR,PROB
c - ale
C     PURPOSE
C        TEST CONVERGENCE OF A LEAST SQUARE FIT BY CONCHK
C
      IBR=0
      IEF=0
      IF(NY.GT.50.OR.NF.GT.20) GOTO 306
      FTEST=0.0
      DO 10 I=1,NF
   10 FTEST=FTEST+ABS(F(I))
      IF(NSTEP.EQ.0) FTESTL=FTEST
      FTESTL=MIN(FTEST,FTESTL)
      IF(NSTEP.EQ.0.AND.IUNPH.NE.0) GOTO 305
      IF(NSTEP.EQ.0) GOTO 40
      IF(IDIAG.NE.0)
     1 WRITE(6,101) NSTEP,NCUT,NCST,DET,CHSQ,FTEST,(F(K),K=1,NF)
      IF(DET.EQ.0.0) GOTO 304
      IF(NCST.EQ.0) GOTO 20
C
C     LAST STEP WAS A CUTSTEP
C     TEST UNPHYSICAL DATA OR DIVERGENCE
      IF(IUNPH.NE.0) GOTO 50
      IF(FTEST.GT.1.1*FTESTL+CST(2)) GOTO 50
      NCST=0
      GOTO 40
C
C     LAST STEP WAS A NORMAL STEP
   20 DTEST=ABS(CHSQ-CHSQV)
      PR=1.0
c - ale
      VPR=CST(4)*CHSQ
      IF (CST(4)*CHSQ.LE.1.E-32) VPR=1.E-10
      IF (ND.GT.0) THEN
        CPR=PROB(VPR,ND)
        PR =DBLE(CPR)
      ENDIF
c - ale
C
C     TEST UNPHYSICAL DATA OR DIVERGENCE
      IF(IUNPH.NE.0) GOTO 50
      IF(FTEST.GT.1.1*FTESTL+CST(2)) GOTO 50
C
C     TEST PROBABILITY/CHISQUARE
      IF(NSTEP.GE.2.AND.PR.LT.CST(1)) GOTO 301
      IF(NSTEP.LT.2.AND.PR.LT.CST(1)) GOTO  40
C
C     TEST CONVERGENCE
   30 IF(FTEST.GT.CST(2)) GOTO 40
      IF(DTEST.GT.CST(3)) GOTO 40
      GOTO 200
C
C     TEST NO OF ITERATIONS
   40 IF(NSTEP.GE.NST(1)) GOTO 302
C
C     PERFORM NEXT ITERATION (NORMAL RETURN)
      GOTO 100
C
C     DIVERGENCE, USE CUTSTEP MODE
   50 IF(NSTEP.GE.NST(1)) GOTO 303
      IF(NCST.NE.0) GOTO 52
      IF(NCUT.GE.NST(2)) GOTO 303
      NCUT=NCUT+1
      FTESTB=FTEST
      NCST=0
   52 IF(NCST.GE.NST(3)) GOTO 56
      NCST=NCST+1
      IUNPH=0
C     CUTSTEP MODE
C     USE HALF THE PREVIOUS CORRECTION
      DO 54 I=1,NY
   54 DELY(I)=0.5*(DELY(I)+DELYP(I))
C
C     CUTSTEP (RETURN 1)
      IBR=1
      GOTO 100
C
C     NEXT STEP, ALTHOUGH LARGER FTEST
   56 NCST=0
      FTESTL=FTESTB
      GOTO 100
C
C     CONVERGENCE, STOP ITERATION (RETURN 2)
  200 IF(IDIAG.NE.0) WRITE(6,102)
      NCASE(LCASE,2)=NCASE(LCASE,2)+1
      NSTP=MIN0(10,NSTEP)
      NITER(LCASE,NSTP)=NITER(LCASE,NSTP)+1
      NPR=1.0+10.0*PR
      NPR=MAX0( 1,NPR)
      NPR=MIN0(10,NPR)
      NPROB(LCASE,NPR)=NPROB(LCASE,NPR)+1
      IBR=2
      GOTO 100
C
C     NO CONVERGENCE, STOP ITERATION (RETURN 3 WITH IER = FAIL CODE)
  306 IEF=IEF+1
  305 IEF=IEF+1
  304 IEF=IEF+1
  303 IEF=IEF+1
  302 IEF=IEF+1
  301 IEF=IEF+1
      NCASE(LCASE,IEF+2)=NCASE(LCASE,IEF+2)+1
      IF(IDIAG.NE.0) WRITE(6,103) IEF
      IBR=3
  100 RETURN
C
  101 FORMAT(I6,I3,'.',I2,G9.2,F7.1,G9.2,3X,10G9.2/(40X,10G9.2))
  102 FORMAT('0CONVERGENCE')
  103 FORMAT('0FIT UNSUCCESSFULL  -  ERROR CODE IER =',I2)
      END
      SUBROUTINE CONDER(IREP)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      COMMON/CONF1/V(1275),DELY(50),F(20),B(1000),PL(50),
     1   CHSQ,PR,NY,NF,IEF,ND,IDIAG,IUNPH
      COMMON/CONF2/DET,FTEST,FTESTL,CHSQV,DELYP(50),
     1   CST(4),NSTEP,NCUT,NCST,NST(3)
      COMMON/CONF3/W(2485),R(70),DG(50)
      INTEGER NFT/0/
      REAL*8 STEP(50),FF(20)
C     PURPOSE
C        OBTAIN MATRIX B( ) = DERIVATIVE MATRIX
C        BY NUMERICAL DIFFERENTIATION
      IREP=0
      IF(NY.GT.50.OR.NF.GT.20) GOTO 100
      IF(NFT.NE.0) GOTO 40
      IF(NSTEP.NE.0) GOTO 15
C     FIRST TIME, STORE STEPS FROM DIAGONAL OF V
      K=0
      DO 10 I=1,NY
      K=K+I
      DGL=V(K)
      IF(DGL.LE.0.0) GOTO 10
      STEP(I)=0.5*SQRT(DGL)
   10 CONTINUE
      GOTO 25
C     NOT THE FIRST TIME, TAKE STEP FROM MATRIX W
   15 K=0
      DO 20 I=1,NY
      K=K+I
      DGL=W(K)*DG(I)*DG(I)
      IF(DGL.LE.0.0) GOTO 20
      STEP(I)=0.5*SQRT(DGL)
   20 CONTINUE
   25 NYF=NY*NF
      DO 30 I=1,NYF
   30 B(I)=0.0
      GOTO 60
C
   40 IF(NFT.GT.NY) GOTO 90
      IF(ISW.EQ.2) GOTO 50
      DO 45 J=1,NF
   45 FF(J)=F(J)
C     MINUS STEP
      DELY(NFT)=SAVE-STEP(NFT)
      ISW=2
      GOTO 80
   50 K=NFT
C     NUMERICAL DIFFERENTIATION
      DO 55 J=1,NF
      B(K)=0.5*(FF(J)-F(J))/STEP(NFT)
   55 K=K+NY
      DELY(NFT)=SAVE
   60 NFT=NFT+1
      IF(NFT.GT.NY) GOTO 80
      IF(STEP(NFT).EQ.0.0) GOTO 60
      SAVE=DELY(NFT)
C     PLUS STEP
      DELY(NFT)=SAVE+STEP(NFT)
      ISW=1
   80 IREP=1
      GOTO 100
   90 NFT=0
      IF(NSTEP.NE.0.OR.IDIAG.EQ.0) GOTO 100
      WRITE(6,101)
      IF(IDIAG.NE.2) GOTO 100
      KA=1
      DO 95 J=1,NF
      KB=KA+NY-1
      WRITE(6,102) J,(B(K),K=KA,KB)
   95 KA=KB+1
  100 RETURN
      ENTRY CONDES(L,STP)
      STEP(L)=ABS(STP)
      RETURN
  101 FORMAT('0MATRIX B BY NUMERICAL DIFFERENTIATION')
  102 FORMAT(I4,6X,10G12.4/(10X,10G12.4))
      END
      SUBROUTINE CONPMT(A,B,NY,NF,DET)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      REAL*8 A(*),B(*),DET
      COMMON/MATCOM/H(50),LR(200),NRD(50)
C     PURPOSE
C     PREPARE MATRIX FOR INVERSION (CALLED BY CONLES)
C
      EXTERNAL SMINV
      LOGICAL*1 LR
      REAL*8 HJ,AJK
      NYF=NY+NF
      NM =0
      K  =0
      DO 10 I=1,NYF
      IF(I.GT.NY) GOTO 2
      K=K+I
      IF(A(K).NE.0.0) GOTO 4
    2 LR(I)=.FALSE.
      GOTO 10
    4 LR(I)=.TRUE.
      NM=NM+1
      NRD(NM)=I
   10 CONTINUE
C
      NY1=NY+1
      DO 50 I=NY1,NYF
      II=(I*I-I)/2
C
      DO 20 M=1,NM
      J=NRD(M)
      HJ=0.0D0
      JK=(J*J-J)/2
      DO 18 K=1,NY
      IF(K.LE.J) JK=JK+1
      IF(.NOT.LR(K)) GOTO 16
      HJ=HJ+A(JK)*A(II+K)
   16 IF(K.GE.J) JK=JK+K
   18 CONTINUE
   20 H(J)=HJ
C
      DO 30 K=I,NYF
      IK=(K*K-K)/2
      AJK=0.0D0
      DO 25 M=1,NM
      J=NRD(M)
   25 AJK=AJK+A(IK+J)*H(J)
      A(IK+I)=A(IK+I)+AJK
   30 CONTINUE
C
      DO 40 M=1,NM
      J=NRD(M)
   40 A(II+J)=-H(J)
C
   50 CONTINUE
C
      CALL XMINV(A,B,NYF,1,DET)
      RETURN
      END
      SUBROUTINE SMINV(A,B,N,MV,DET)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      REAL*8 A(*),B(*),DET,AM,D,E
      REAL*8 EPS/1.E-6/
      COMMON/MATCOM/H(50),LR(200),NRD(50)
      DIMENSION DIAG(200)
      LOGICAL*1 LR
C
C     PURPOSE
C        OBTAIN SOLUTION OF A SYSTEM OF LINEAR EQUTIONS A * X = B WITH
C        SYMMETRIC MATRIX A AND INVERSE (FOR MV = 1) OR MATRIX
C        INVERSION ONLY (FOR MV = 0)
C
C     USAGE
C                   - - - --
C        CALL SMINV(A,B,N,MV,DET)
C                   - -      ---
C
C           A = SYMMETRIC N-BY-N MATRIX IN SYMMETRIC STORAGE MODE
C               A(1) = A11, A(2) = A12, A(3) = A22, A(4) = A13, . . .
C               REPLACED BY INVERSE MATRIX
C           B = N-VECTOR   (FOR MV = 0 USE A DUMMY ARGUMENT)
C               REPLACED BY SOLUTION VECTOR
C          MV = SEE ABOVE
C         DET = DETERMINANT OF A
C
C
C     RESTRICTION   N LE 200
C     METHOD OF SOLUTION IS BY ELIMINATION USING THE LARGEST PIVOTAL
C     DIVISOR AT EACH STAGE. A DETERMINANT OF ZERO INDICATES
C     A SINGULAR MATRIX. IN THIS CASE ALL REMAINING ROWS AND COLS OF
C     MATRIX A ARE SET TO ZERO.
C
C          CONSTRUCT TABLE
      DO 10 I=1,N
   10 LR(I)=.FALSE.
      NI=N
      GOTO 14
      ENTRY XMINV(A,B,N,MV,DET)
      NI=0
      DO 12 I=1,N
      IF(LR(I)) GOTO 12
      NI=NI+1
   12 CONTINUE
   14 DET=1.0
      M=0
      DO 16 I=1,N
      M=M+I
   16 DIAG(I)=ABS(A(M))
C          LOOP BEGIN
      DO 60 I=1,NI
C          SEARCH FOR PIVOT
      K=0
      M=0
      AM=0.0
      DO 20 J=1,N
      M=M+J
      IF(LR(J)) GOTO 20
      IF(ABS(A(M)).LE.AM) GOTO 20
C          TEST FOR LINEARITY
      IF(ABS(A(M)).LT.EPS*DIAG(J)) GOTO 20
      AM=ABS(A(M))
      K=J
      KK=M
   20 CONTINUE
C          TEST FOR ZERO MATRIX
      IF(K.EQ.0) GOTO 90
C          PREPARATION FOR ELIMINATION
      LR(K)=.TRUE.
      DET=DET*A(KK)
      D=1.0/A(KK)
      A(KK)=-D
      IF(MV.EQ.1) B(K)=B(K)*D
      JK=KK-K
      JL=0
C          ELIMINATION
      DO 50 J=1,N
      IF(J-K) 24,22,26
   22 JK=KK
      JL=JL+J
      GOTO 50
   24 JK=JK+1
      GOTO 28
   26 JK=JK+J-1
   28 E=A(JK)
      A(JK)=D*E
      IF(MV.EQ.1) B(J)=B(J)-B(K)*E
      LK=KK-K
      DO 40 L=1,J
      JL=JL+1
      IF(L-K) 34,32,36
   32 LK=KK
      GOTO 40
   34 LK=LK+1
      GOTO 38
   36 LK=LK+L-1
   38 A(JL)=A(JL)-A(LK)*E
   40 CONTINUE
   50 CONTINUE
C          LOOP END
   60 CONTINUE
C          CHANGE SIGN
      M=0
      DO 80 I=1,N
      DO 70 J=1,I
      M=M+1
   70 A(M)=-A(M)
   80 CONTINUE
      GOTO 100
   90 DET=0.0
C          CLEAR REST OF MATRIX
      M=0
      DO 95 I=1,N
      DO 95 J=1,I
      M=M+1
      IF(.NOT.LR(I)) A(M)=0.0
      IF(.NOT.LR(J)) A(M)=0.0
   95 A(M)=-A(M)
  100 RETURN
      END
      FUNCTION INDS(I,J)
C
C     PURPOSE
C        OBTAIN LINEAR INDEX FOR ELEMENT IJ OF A MATRIX IN SYMMETRIC
C        STORAGE MODE 11 12 22 13 23 33 14 . . .
C
C     USAGE
C                     - -
C        INDEX = INDS(I,J)
C        -----
C        WHERE I,J = INDICES OF MATRIX ELEMENT
C
C     EXAMPLE
C        INDS(1,4) = INDS(4,1) = 7
C
      IF(I-J) 10,10,20
   10 INDS=I+(J*J-J)/2
      GOTO 100
   20 INDS=J+(I*I-I)/2
  100 RETURN
      END
      SUBROUTINE SSVW(V,I,W,N)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      REAL*8 V(*),W(*)
C
C     PURPOSE
C        STORE SYMM. N-BY-N SUBMATRIX W IN SYMM. MATRIX V
C        AT ROWS/COLS   I, I+1, . . . I+N-1
C
C     USAGE
C                  - - - -
C        CALL SSVW(V,I,W,N)
C                  -
C
      IV=(I*I+I)/2-1
      IW=0
      DO 20 K=1,N
      DO 10 L=1,K
      IV=IV+1
      IW=IW+1
   10 V(IV)=W(IW)
   20 IV=IV+I-1
      RETURN
      END
      SUBROUTINE SSWV(V,I,W,N)
c - ale
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c - ale
      REAL*8 V(*),W(*)
C
C     PURPOSE
C        GET SYMM. N-BY-N SUBMATRIX OUT OF
C        SYMM. MATRIX V AT ROWS/COLS I,I+1, . . . I+N-1
C
C     USAGE
C                  - -   -
C        CALL SSWV(V,I,W,N)
C                      -
C
      IV=(I*I+I)/2-1
      IW=0
      DO 20 K=1,N
      DO 10 L=1,K
      IV=IV+1
      IW=IW+1
   10 W(IW)=V(IV)
   20 IV=IV+I-1
      RETURN
      END
