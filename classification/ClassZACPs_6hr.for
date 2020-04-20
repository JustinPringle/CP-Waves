C *******************************************************************
C Program for the optimization of fuzzy rules for CP Classification
C    Goals -  Flood relevant meteorological situations
C             using surface pressure data
C
C
C    Version 0.1    Lausanne September 12 2001
C                   Andras Bardossy
C                   
C *******************************************************************
      COMMON /TPERIOD/ NBY,NEY,NBM,NEM,NBD,NED
      COMMON /FIRSTYEAR/ NRBY
      
      CHARACTER INF*80
      CHARACTER FILNA*80
C
C Read control file
C
      OPEN(81,FILE='ZAclasscps.con')
      WRITE(*,*)  ' You are using a program for Fuzzy Classification'
      WRITE(*,*)  '           Andras Bardossy               '
      WRITE(*,*)
      
C
C     Get name of input file
C
      READ (81,'(A)') FILNA
      WRITE(*,'(A)') FILNA
      READ(81,*) NRBY
      
       
      OPEN(30,File=FILNA,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=1500)
      READ(81,*) NBY,NEY      
      DO 1 I=1,1
	  NBM=1
      NBD=1
      NBHr =0
	  NEM=12
	  NED=31
      NEHr=18

      CALL INITIAL
      CALL CLASS
      CLOSE(81)
1     CONTINUE
C	GOTO 1
99    STOP
      END

	SUBROUTINE JDATE(NY,NM,ND,JTAG)
C *******************************************************************
C Julian date
C *******************************************************************
	DIMENSION MLEN(12),MSLEN(13)
	DO 1 I=1,12
	  MLEN(I)=4*31
1     CONTINUE
      MLEN(4)=4*30
	MLEN(6)=4*30
	MLEN(9)=4*30
	MLEN(11)=4*30
	MLEN(2)=4*28
	IY=NY/4
	IY=IY*4-NY
	IF(IY.EQ.0) THEN
	  MLEN(2)=4*29
      ENDIF
	MSLEN(1)=0
      DO 2 I=1,12
	MSLEN(I+1)=MSLEN(I)+MLEN(I)
2     CONTINUE
C     4 readings in one day (one evry 6 hrs
      IF(NHr.eq.0) THEN
        Time = 0
      ELSEIF(NHr.eq.6) THEN
        Time = 1
      ELSEIF(NHr.eq.12) THEN
        Time = 2
      ELSEIF(NHr.eq.18) THEN
        Time = 3
      ENDIF
      JTAG=MSLEN(NM)+ 4*(ND-1)+1 + Time
	RETURN
	END

      SUBROUTINE INITIAL
C ******************************************************************
C Read Pressure data, rules, etc.
C ******************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /ATLAG/ TNATL,XMAX1,XMIN1
      COMMON /NOIDAT/ NXI,NYI,NXA,NYA,VALISL,VALAZO
      COMMON /DAILY/ CPR(7400,25,26)
      COMMON /CLASSES/ CDOF(4700,30),NTAG
      COMMON /MONAT/ MON(4700)
	COMMON /IDO/ NAP(4700),IEV(4700)
      COMMON /OFIL/ IF1,IF2,FILNA
	CHARACTER FILNA*80
      CHARACTER*1  CFC, answer
      CHARACTER*50 INF, AVGNAM,SIGNAM,OUTNAM,GCFILE
      character*1 ichr(3000)
      real E(25,26)
      COMMON /AVERA/ AVG(4*366,25,26)
      COMMON /SIGRA/ SIG(4*366,25,26)
      COMMON /TPERIOD/ NBY,NEY,NBM,NEM,NBD,NED
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
      
      DIMENSION BA(400),BS(400)
C
C     Get name of output file
C
c      READ (81,'(A)') INF
c      OPEN(21,FILE=INF)
c      WRITE(*,'(A)') INF

c45    READ (81,'(A)') OUTNAM
c      IF(OUTNAM.EQ.' ')  OUTNAM='TEST.CLS'
C      OPEN(9,FILE=OUTNAM)

10    CONTINUE
C
C     Determine which grid user is interested in
C
c      READ(81,*) nby,nbm,nbd,nbh
c      READ(81,*) ney,nem,ned,neh

C
C Read filename for average pressure data
C
      READ (81,'(A)') AVGNAM
      write(*,*) AVGNAM
      READ (81,'(A)') SIGNAM
      write(*,*) SIGNAM
      OPEN(28,FILE=AVGNAM)
      OPEN(29,FILE=SIGNAM)
C Grid descrition
	LOND=21
	LATD=17
	
	
C
C Calculate seasonal averages in required resolution
C
      write(*,*) ' Reading averages'


C Daily smoothed means are used
       DO 703 J=1,4*365
	   READ(28,*) KDUM
	   READ(29,*) KDUM
         READ(28,*) (BA(K),K=1,LOND*LATD)
         READ(29,*) (BS(K),K=1,LOND*LATD)
         
	   LL=1
	  DO 704 I=1,LATD
	  DO 704 K=1,LOND
	    AVG(J,K,I)=BA(LL)
	    SIG(J,K,I)=BS(LL)
	    LL=LL+1
704     CONTINUE
703    CONTINUE
        J=4*366
	  DO 712 I=1,LATD
	   DO 712 K=1,LOND
	     AVG(J,K,I)=AVG(4*365,K,I)
	     SIG(J,K,I)=SIG(4*365,K,I)
712    CONTINUE
      CLOSE(28)
      CLOSE(29)

C Grid descrition

	XMINLT=-10.0
	XMINLD=0.0
	
      DLAT=2.5
	DLON=2.5

      READ (81,*) XS,XW
      READ (81,*) XN,XE
C
      CALL READRULE(XS,XW,XN,XE)
C
C Convert rule locations to indices
C
      CALL CONVRULE(LOND,LATD,XMINLT,XMINLD,DLAT,DLON)
        READ(81,*) FILNA
	  READ(81,*) IF1,IF2
      CLOSE(81)

C      OPEN(77,FILE='d:\nmc\fors\testdate')
C Read selected area  (SW-NE) corners
C
C Find four corners as indices
      XMIN1=99999.9
      DO 901 I=1,LOND
         XH=ABS(XW+FLOAT(I-1)*DLON-XMINLD)
         IF(XH.LT.XMIN1) THEN
           INW=I
           XMIN1=XH
         ENDIF
901   CONTINUE
      XMIN1=99999.9
      DO 902 I=1,LOND
         XH=ABS(XE+FLOAT(I-1)*DLON-XMINLD)
         IF(XH.LT.XMIN1) THEN
           INE=I
           XMIN1=XH
         ENDIF
902   CONTINUE
      XMIN1=99999.9
      DO 903 I=1,LATD
         XH=ABS(XS+FLOAT(I-1)*DLON-XMINLT)
         IF(XH.LT.XMIN1) THEN
           INS=I
           XMIN1=XH
         ENDIF
903   CONTINUE
      XMIN1=99999.9
      DO 904 I=1,LATD
         XH=ABS(XN+FLOAT(I-1)*DLON-XMINLT)
         IF(XH.LT.XMIN1) THEN
           INN=I
           XMIN1=XH
         ENDIF
904   CONTINUE

C
C Begin day by day reading
C
C     Read the specified grid into array 
C
      IDELTHR=6
      NTAG=0
      NR=1
	IYR=NBY
	IMO=NBM
	IDY=NBD
	IHR=NBHr
322   CALL cdslpZA (NR,IYR,IMO,IDY,IHR,E,NERR)
      IF(IYR.LT.100) IYR=IYR+1900
      IF(NERR.EQ.-1) THEN
        WRITE(*,*) 'PRESSURE DATA READ'
        WRITE(*,*) ' TOTAL ',NTAG,' DAYS'
        WRITE(*,*) 'Last day =', IYR,IMO,IDY
        RETURN
      ENDIF

      if(iyr.lt.nby) goto 765
      if(iyr.eq.nby. and .imo.lt.nbm ) goto 765
      if(iyr.eq.nby. and .imo.eq.nbm . and .idy.lt.nbd) goto 765
c      if(ihr.lt.3) goto 765
607   CONTINUE
      IF(NTAG.GT.0) THEN
      IF((KYR.NE.IYR).OR.(KMO.NE.IMO).OR.(KDY.NE.IDY)) THEN
       NTAG=NTAG+1
       WRITE(77,*) IYR,IMO,IDY,NTAG
       WRITE(*,*) IYR,IMO,IDY,NTAG
       CPR(NTAG,1,1)=-1.0
       MYR=KYR
       MN=KMO
       MDY=KDY
       
       CALL CHDATE (MYR,MN,MDY,MHR,KYR,KMO,KDY,KHR,IDELTHR)
       GOTO 607
      ENDIF
      ENDIF
C
C Subrtact mean values and divide with std
      NTAG=NTAG+1
	CALL JDATE(IYR,IMO,IDY,JTAG)
	IF(E(1,1).GE.0.0) THEN
      DO 601 I=1,LOND
       DO 602 J=1,LATD
         EHE=(E(I,J)-AVG(JTAG,I,J))/SIG(JTAG,I,J)
         IF(ABS(EHE).GT.6.0) THEN
           LOLO=1
         ENDIF
         
C Change norming  first only
         EHE=0.33*EHE+0.5
C        EHE=0.5*EHE+0.5
        IF(EHE.GE.1.0) EHE=1.0
        IF(EHE.LE.0.) EHE=0.0
        CPR(NTAG,I,J)=EHE
602    CONTINUE
601   CONTINUE
      ENDIF
      WRITE(77,*) IYR,IMO,IDY,NTAG
      MON(NTAG)=IMO
      NAP(NTAG)=JTAG
	IEV(NTAG)=IYR
2019  format(4I4,I6)
2001  format(11F7.1)
      IDELTHR=6
      MYR=IYR
      MN=IMO
      MDY=IDY
      MHR=IHR
      CALL CHDATE (MYR,MN,MDY,MHR,KYR,KMO,KDY,KHR,IDELTHR)
	IYR=KYR
	IMO=KMO
	IDY=KDY
765   CONTINUE
      if(iyr.gt.ney) THEN
        WRITE(*,*) 'PRESSURE DATA READ'
        WRITE(*,*) ' TOTAL ',NTAG,' DAYS'
        RETURN
      ENDIF
      GOTO 322
      END

      SUBROUTINE CHDATE (MYR,MN,MDY,MHR,MYRCH,MNCH,MDYCH,MHRCH,IDELTHR)
C
C         MYR = Original year
C          MN = Original month
C         MDY = Original day of month
C         MHR = Original hour of day
C       MYRCH = Changed year
C        MNCH = Changed month
C       MDYCH = Changed day of month
C       MHRCH = Changed hour of day
C     IDELTHR = Number of hours by which to change original date/time
C

      INTEGER*2 NDAYS(12)
      DATA NDAYS /31,28,31,30,31,30,31,31,30,31,30,31/
      NDAYS(2)=28
      IF (MOD(MYR,4).EQ.0) NDAYS(2) = 29

      MHRCH = MHR + IDELTHR
      MDYCH = MDY
      MNCH = MN
      MYRCH = MYR

10    IF (MHRCH.GE.0 .AND. MHRCH.LT.24) RETURN
      IF (MHRCH.LT.0) THEN
         MHRCH = MHRCH + 24
         MDYCH = MDYCH - 1
         IF (MDYCH.EQ.0) THEN
            MNCH = MNCH - 1
            IF (MNCH.EQ.0) THEN
               MNCH =  12
               MYRCH = MYRCH - 1
               IF (MOD(MYRCH,4) .EQ. 0) NDAYS(2)=29
            ENDIF
            MDYCH = NDAYS(MNCH)
         ENDIF
      ENDIF
      IF (MHRCH.GT.23) THEN
         MHRCH = MHRCH - 24
         MDYCH = MDYCH + 1
         IF (MDYCH.GT.NDAYS(MNCH)) THEN
            MNCH = MNCH + 1
            IF (MNCH.EQ.13) THEN
               MNCH  = 1
               MYRCH = MYRCH + 1
               IF (MOD(MYRCH,4) .EQ. 0) NDAYS(2)=29
            ENDIF
            MDYCH = 1
         ENDIF
      ENDIF
      GO TO 10

      END
      
      SUBROUTINE READRULE(XS,XW,XN,XE)
C *****************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
	COMMON /REALRULE/ XLOC1(30,40),XLOC2(30,40)
      DIMENSION NHELP(30)
      character OUTNAM*80
      READ (81,'(A)') OUTNAM
      OPEN(50,FILE=OUTNAM)
      IOLD=0
      ISUM=0
1     READ(50,*,END=2) I1,I2,I3,VI4,VI5
      WRITE(*,*) I1,I2,I3,VI4,VI5
      I5=INT(VI5+0.1)
      IF(VI4.GT.0.0) then
      I4=INT(VI4+0.1)
      else
      I4=INT(VI4-0.1)
      endif
      HEW=(VI5-XW)*(VI4-XE)
C	IF(HEW.GT.0.0) GOTO 1
      HNS=(VI4-XS)*(VI5-XN)
	IF(HNS.GT.0.0) GOTO 1

c      write(*,*) I1,I2,I3,VI4,VI5

      IF(I1.EQ.99) GOTO 2
      IF(I1.NE.IOLD) THEN
        IF(ISUM.GT.0) NSUM(ISUM)=NADD
        ISUM=ISUM+1
        NADD=1
        IOLD=I1
      ELSE
        NADD=NADD+1
      ENDIF
      NTYPE(ISUM,NADD)=I2
      XLOC1(ISUM,NADD)=VI4
      XLOC2(ISUM,NADD)=VI5
      GO TO 1
2     IF(ISUM.GT.0)  NSUM(ISUM)=NADD
      NTOTR=ISUM
C Find number of types represented 1 ... 4
      DO 3 I=1,NTOTR
        DO 4 J=1,4
         NHELP(J)=0
4       CONTINUE
         DO 5 J=1,NSUM(I)
           NHELP(NTYPE(I,J))=1
5        CONTINUE
        NHH=0
        DO 6 J=1,4
           NHH=NHH+NHELP(J)
6       CONTINUE
      NUTY(I)=NHH
3     CONTINUE
      CLOSE(50)
      RETURN
      END

      SUBROUTINE CONVRULE(LOND,LATD,XMINLT,XMINLD,DLAT,DLON)
C ************************************************************************
C Convert rules to actual grid
C ************************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /RULES2/ NRMAT(30,25,26)
      COMMON /RULEOLD/ NOLDMAT(30,25,26),MOLD
	COMMON /REALRULE/ XLOC1(30,40),XLOC2(30,40)
      REAL EE(100),FF(100)
      NRMAT=0
      DO 1 I=1,LOND
         EE(I)=-FLOAT(I-1)*DLON+XMINLD
c	   IF(EE(I).LT.0.0) EE(I)=360.0+EE(I)
1     CONTINUE
      DO 2 J=1,LATD
         FF(J)=-FLOAT(J-1)*DLAT+XMINLT
2     CONTINUE
      DO 3 I=1,NTOTR
      DO 4 J=1,NSUM(I)
        XX=XLOC2(I,J)

         IND1=0
         XMI=99999.9
         DO 5 I1=1,LOND
          IF(ABS(XX-EE(I1)).LT.XMI) THEN
            XMI=ABS(XX-EE(I1))
            IND1=I1
          ENDIF
5        CONTINUE
        NLOC1(I,J)=IND1
        XX=XLOC1(I,J)
         IND2=0
         XMI=99999.9
         DO 6 I1=1,LATD
          IF(ABS(XX-FF(I1)).LT.XMI) THEN
            XMI=ABS(XX-FF(I1))
            IND2=I1
          ENDIF
6        CONTINUE
        NLOC2(I,J)=IND2
        
        NRMAT(I,NLOC1(I,J),NLOC2(I,J))=NTYPE(I,J)
4       CONTINUE
3       CONTINUE
       OPEN(41,FILE='Cporig')
c       OPEN(41,FILE='Cporig')
       DO 21 I=1,NTOTR
       WRITE(41,*) I
        DO 22 J=1,LATD
         WRITE(41,'(80I2)')(NRMAT(I,K,J),K=1,LOND)
         DO 22 K=1,LOND
          NOLDMAT(I,K,J)=NRMAT(I,K,J)
22      CONTINUE
21     CONTINUE
       CLOSE(41)
       MOLD=0
      RETURN
      END
      


      SUBROUTINE CLASS
C ****************************************************************
C Subroutine Class
C ****************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /ATLAG/ TNATL,XMAX2,XMIN2
      COMMON /NOIDAT/ NXI,NYI,NXA,NYA,VALISL,VALAZO
      COMMON /DAILY/ CPR(7400,25,26)
      COMMON /CLASSES/ CDOF(4700,30),NTAG
      COMMON /OBJECTI/ NOBD,NNBD,NCOB
      COMMON /PRCS/ PRCVAL(4700,30),NPSTAT,NPCTAG
      COMMON /MONAT/ MON(4700)
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /PRCATL/ LDRYAVG(30,12),LWETAVG(30,12),PSAVG(30,12),DRYLIM
      COMMON /PRCAVG2/ LDRYAVG2(30,12),LWETAVG2(30,12),EXTLIM
	COMMON /IDO/ NAP(4700),IEV(4700)
      COMMON /OFIL/ IF1,IF2,FILNA
      REAL DOF(30),ATM(10),AIDB(10),PEX(5)
      REAL AAND(10),AOR(10)
	CHARACTER FILNA*80
C Open output file No. 1
        IF(IEV(1).GT.1000) THEN
	    IVV=IEV(1)
	  ELSE
	    IVV=IEV(1)+1900
	  ENDIF
        IF(IVV.GT.9) THEN
	    WRITE(FILNA(IF1:IF2),'(I4)') IVV
        ELSE
	    WRITE(FILNA(IF1:IF1),'(A)') '0'
	    WRITE(FILNA(IF2:IF2),'(I1)') IVV
        ENDIF
	  OPEN(67,FILE=FILNA)
        JOTAG=0
C Calculate DOF-s for E  Lp distance
      PEX(1)=1.25
      PEX(2)=5.0
      PEX(3)=5.0
      PEX(4)=1.25
C     PEX(1)=2.5
C     PEX(2)=2.5
C     PEX(3)=2.5
C     PEX(4)=2.5
      NOBD=0
      NCOB=0
      DO 101 ITG=1,NTAG
      IF(CPR(ITG,1,1).LT.0.) THEN
        DO 102 KL=1,NTOTR
         CDOF(ITG,KL)=0.0
102     CONTINUE
        CDOF(ITG,NTOTR+1)=1.0
        GOTO 101
      ENDIF
      DO 5 I=1,NTOTR
         DOF(I)=0.0
         DO 6 J=1,4
           ATM(J)=0.0
           AIDB(J)=0.0
           AAND(J)=0.0
           AOR(J)=0.0
6        CONTINUE
         DO 7 J=1,NSUM(I)
           XNORM=CPR(ITG,NLOC1(I,J),NLOC2(I,J))
           XNORMP1 = CPR(ITG+1,NLOC1(I,J),NLOC2(I,J)) 
            PEXP=PEX(NTYPE(I,J))
            IF(NTYPE(I,J).EQ.1) THEN
              AMB=XMEM1(1.0-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0) ATM(1)=ATM(1)+AMB**PEXP
              AIDB(1)=AIDB(1)+1.0
c              IF(AIDB(1).LT.1.5) THEN
c                AAND(1)=AMB
c                AOR(1)=AMB
c              ELSE
c                AAND(1)=AAND(1)*AMB
c                AOR(1)=AOR(1)+AMB-AOR(1)*AMB
c              ENDIF
            ENDIF
            IF(NTYPE(I,J).EQ.2) THEN
              AMB=XMEM2(1.-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(2)=ATM(2)+AMB**PEXP
                AIDB(2)=AIDB(2)+1.0
              IF(AIDB(2).LT.1.5) THEN
                AAND(2)=AMB
                AOR(2)=AMB
              ELSE
                AAND(2)=AAND(2)*AMB
                AOR(2)=AOR(2)+AMB-AOR(2)*AMB
              ENDIF
            ENDIF
            IF(NTYPE(I,J).EQ.3) THEN
              AMB=XMEM2(XNORM,0)
              AMB=1.0-AMB
                IF(AMB.GT.0.0) ATM(3)=ATM(3)+AMB**PEXP
                AIDB(3)=AIDB(3)+1.0
              IF(AIDB(3).LT.1.5) THEN
                AAND(3)=AMB
                AOR(3)=AMB
              ELSE
                AAND(3)=AAND(3)*AMB
                AOR(3)=AOR(3)+AMB-AOR(3)*AMB
              ENDIF
            ENDIF
            IF(NTYPE(I,J).EQ.4) THEN
              AMB=XMEM1(XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(4)=ATM(4)+AMB**PEXP
                AIDB(4)=AIDB(4)+1.0
              IF(AIDB(4).LT.1.5) THEN
                AAND(4)=AMB
                AOR(4)=AMB
              ELSE
                AAND(4)=AAND(4)*AMB
                AOR(4)=AOR(4)+AMB-AOR(4)*AMB
              ENDIF
            ENDIF
C        write(9,*) NTYPE(I,J),XNORM,AMB,AIDB(NTYPE(I,J))
7        CONTINUE
         AEXP=1.0
         AMULT=1.0
         DO 8 J=1,4
            PEXP=PEX(J)
           IF(AIDB(J).GT.0.0) THEN
             IF(ATM(J).GT.0.0) THEN
               ATM(J)=1.0-(ATM(J)/AIDB(J))**(1./PEXP)
             ELSE
                ATM(J)=1.0
             ENDIF
C new   AND-OR
c            IF((J.EQ.1).OR.(J.EQ.4)) ATM(J)=AOR(J)*0.6+AAND(J)*0.4
c            IF((J.EQ.2).OR.(J.EQ.3)) ATM(J)=AOR(J)*0.1+AAND(J)*0.9
C new end
           ELSE
             ATM(J)=1.0
             AEXP=4./3.
           ENDIF
           AMULT=AMULT*ATM(J)
8        CONTINUE
         IF(AMULT.GT.0.0) AMULT=AMULT**AEXP
         DOF(I)=AMULT
C       write(*,1112) I,DOF(I),ATM(1),ATM(2),ATM(3),ATM(4)
1112  FORMAT(I6,5F10.5)
5     CONTINUE
      DMAX=0.0
      NCL=NTOTR+1
      NCL=99
      DO 9 I=1,NTOTR
C       Change the degree of fit
        CDOF(ITG,I)=alp*CDOF(ITG-1,I)+alp*CDOF(ITG+1,I)+(1-2*alp)*DOF(I)
        IF(DOF(I).GT.DMAX) THEN
          DMAX=DOF(I)
          NCL=I
        ENDIF
9     CONTINUE
      WRITE(61,'(3I6,30F10.6)') IEV(ITG),NAP(ITG),NCL
     x   ,(DOF(I1),I1=1,NTOTR)
109   IF(NAP(ITG).EQ.JOTAG+1) THEN
        IF(NCL.LT.10) THEN
	    WRITE(67,1900) NCL
	  ELSE
	    WRITE(67,1901) NCL
	  ENDIF
        JOTAG=JOTAG+1
	ELSE
	  IF(NAP(ITG).LT.JOTAG) THEN
          CLOSE(67)
	    IF(IEV(ITG).GT.1000) THEN
	       IVV=IEV(ITG)
	    ELSE
	       IVV=IEV(ITG)+1900
	    ENDIF
	    IF(IVV.GT.9) THEN
	      WRITE(FILNA(IF1:IF2),'(I4)') IVV
	    ELSE
	      WRITE(FILNA(IF1:IF1),'(I1)') 0
	      WRITE(FILNA(IF2:IF2),'(I1)') IVV
          ENDIF	      
		OPEN(67,FILE=FILNA)
          IF(NCL.LT.10) THEN
	      WRITE(67,1900) NCL
	    ELSE
	      WRITE(67,1901) NCL
	    ENDIF
          JOTAG=1
	  ELSE
	    WRITE(67,1902)
		JOTAG=JOTAG+1
		GOTO 109
	  ENDIF
	ENDIF
101   CONTINUE
1900  FORMAT('CP0',I1)
1901  FORMAT('CP',I2)
1902  FORMAT('CP99')
      RETURN
      END

      REAL FUNCTION XMEM1(X,I1)
C *********************************************************************
C Calculates membership values classes 1-4
C *********************************************************************
      XMEM1=0.0
      IF(X.LE.0.0) XMEM1=1.0
      IF(X.GE.0.5) XMEM1=0.0
      IF((X.GT.0.0).AND.(X.LT.0.5)) XMEM1=1.-X/0.5
      IF(I1.EQ.1) XMEM1=1.0-XMEM1
      RETURN
      END

      REAL FUNCTION XMEM2(X,I1)
C *********************************************************************
C Calculates membership values classes 2-3
C *********************************************************************
      XMEM2=0.0
      IF(X.LE.0.2) XMEM2=(X+0.3)/0.5
      IF(X.GE.0.6) XMEM2=0.0
      IF((X.GT.0.2).AND.(X.LT.0.6)) XMEM2=1.-(X-0.2)/0.4
      IF(I1.EQ.1) XMEM2=1.0-XMEM2
      RETURN
      END

      SUBROUTINE CDSLPNEWOLD (NR,IYA,IMO,ITA,IHO,A,NERR)
C ****************************************************************
C Read ASCII CD
C ****************************************************************
      DIMENSION A(72,15)
	CHARACTER DUMMY*20
	DO 2 J=1,3
	  READ(21,'(A)',END=999) DUMMY
2     CONTINUE
      READ(21,*) IDU,IYA,IMO,ITA,IHO
	DO 4 J=1,2
	  READ(21,'(A)',END=999) DUMMY
4     CONTINUE
      DO 3 L1=1,72
	   READ(21,'(7X,15F8.1)') (A(L1,L2),L2=1,15)
3     CONTINUE
      NERR=0
	RETURN
999   NERR=-1
      CLOSE(21)
      RETURN
      END
      SUBROUTINE CDSLPZA (NR,IYA,IMO,ITA,IHO,A,NERR)
C ****************************************************************
C Read Binary data
C ****************************************************************
      COMMON /AVERA/ AVG(4*366,25,26)
      COMMON /SIGRA/ SIG(4*366,25,26)
      COMMON /FIRSTYEAR/ NRBY
      
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
      
      DIMENSION A(25,26),B(1080)
	CHARACTER DUMMY*20
C
      NERR=0
	CALL JDATE(IYA,IMO,ITA,JTAG)
c	IREC=(IYA-1979)*366+JTAG
	IREC=4*(IYA-NRBY)*366+JTAG
      READ(30,REC=IREC) KYA,KMO,KTA,(B(L1),L1=1,LOND*LATD)
      L3=1
      DO 3 L2=1,LATD
        DO 3 L1=1,LOND
	   A(L1,L2)=B(L3)
c	   A(L1,L2)=(B(L3)-AVG(JTAG,L1,L2))/SIG(JTAG,L1,L2)
	   L3=L3+1
3     CONTINUE
      NERR=0
      IF(IYA.NE.KYA) THEN
      DO 31 L1=1,LOND
        DO 31 L2=1,LATD
	   A(L2,L1)=-9.0
31    CONTINUE
	WRITE(*,*) IYA,KYA,IMO,KMO,IREC
	  NERR=-9
	ENDIF
      IF(IMO.NE.KMO) NERR=-9
      IF(ITA.NE.KTA) NERR=-9
	RETURN
999   NERR=-1
     
      RETURN
      END
