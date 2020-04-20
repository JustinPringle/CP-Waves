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
c
c     
      CHARACTER FILNA*100
      COMMON /FIRSTYEAR/ NRBY
      
      OPEN(81,FILE='ZAoptcps.con')
      WRITE(*,*)  ' You are using a program for Fuzzy Classification'
      WRITE(*,*)  '           Andras Bardossy               '
      WRITE(*,*)
C
C     Get name of output file
C
      READ (81,'(A)') FILNA
      WRITE(*,'(A)') FILNA
       
      OPEN(30,File=FILNA,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=1500)
      READ(81,*) NRBY
      write(82,'(I4)')NRBY
      CALL INITIAL
      CLOSE(30)
      WRITE(60,*)'THE STARTING NUMBER OF RULES ARE:',NTOTR     
      CALL INITPREC      
      CALL ICLASS
      CALL OPTIMIZE
      CALL WRITERULE
      STOP
      END

	SUBROUTINE JDATE(NY,NM,ND,NHr,JTAG)
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

      SUBROUTINE INITPREC
C ******************************************************************
C Read Precipitation data
C ******************************************************************
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /TPERIOD/ NBY,NEY
      DIMENSION KPRC(300)
      DIMENSION NSEL(30)
      CHARACTER GCFILE*60
C
C  Read precipitation data (max 30 stations according PRCVAL(.,.))
C
	READ(81,*) NPALL
	READ(81,*) NPSTAT
C
C  NSEL = SELECTED STATIONS
C
      DO 2 I=1,NPSTAT
      READ(81,*) NSEL(I)
2     CONTINUE
C  FILE-NAMES
      READ(81,'(A)') GCFILE
	READ(81,*) IA1,IE1
      ITGS=1
      DO 620 I=NBY,NEY
	IYT=I
c      IF(I.GT.1000) THEN
c	  IYT=I-1900
c	ELSE
c	  IYT=I
c	ENDIF
      WRITE(GCFILE(IA1:IE1),'(I4)') IYT
      OPEN(79,FILE=GCFILE,STATUS='OLD')
621    READ(79,'(7X,288I6)',END=629) (KPRC(L),L=1,NPALL)
       DO 1 L=1,NPSTAT
         PRCVAL(ITGS,L)=FLOAT(KPRC(NSEL(L)))/10.
1      CONTINUE
	WRITE(33,*)PRCVAL(ITGS,1)
       ITGS=ITGS+1
       GOTO 621
629   CLOSE(79)
620   CONTINUE
      NPCTAG=ITGS-1
      WRITE(*,*) 'DAYS with precip data' ,NPCTAG
C
C Initialize Precipitation indices
C
      CLOSE(81)
      CALL INITPIDX
      RETURN
      END

      SUBROUTINE INITPIDX
C ******************************************************************
C Initialize precipitation indices for original classification
C ******************************************************************
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /MONAT/ MON(7400)
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /PRCATL/ LDRYAVG(30,12),LWETAVG(30,12),PSAVG(30,12),DRYLIM
      COMMON /PRCAVG2/ LDRYAVG2(30,12),LWETAVG2(30,12),EXTLIM
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA
C
      DRYLIM=25
      EXTLIM=35
      DO 1 M=1,12
      DO 1 I=1,30
        LDRYAVG(I,M)=0
        LWETAVG(I,M)=0
        LDRYAVG2(I,M)=0
        LWETAVG2(I,M)=0
        PSAVG(I,M)=0.0
        DO 1 J=1,30
          LWET(I,J,M)=0
          LWET2(I,J,M)=0
          LDRY(I,J,M)=0
          LDRY2(I,J,M)=0
          PS(I,J,M)=0.0
1     CONTINUE
      DO 2 I=1,NPCTAG
C        IF(MON(I).EQ.12) MON(I)=0
C        MS=MON(I)/3+1
         MS=MSEA(MON(I))
c        MS=1
c        IF((MON(I).LT.5).OR.(MON(I).GT.10)) MS=2
        DO 3 J=1,NPSTAT
          IF(PRCVAL(I,J).GT.-0.1) THEN
            IF(PRCVAL(I,J).GT.DRYLIM) THEN
              LWETAVG(J,MS)=LWETAVG(J,MS)+1
            ELSE
              LDRYAVG(J,MS)=LDRYAVG(J,MS)+1
            ENDIF
            IF(PRCVAL(I,J).GE.EXTLIM) THEN
              LWETAVG2(J,MS)=LWETAVG2(J,MS)+1
            ELSE
              LDRYAVG2(J,MS)=LDRYAVG2(J,MS)+1
            ENDIF
            PSAVG(J,MS)=PSAVG(J,MS)+PRCVAL(I,J)
          ENDIF
3       CONTINUE
2     CONTINUE
      RETURN
      END

      SUBROUTINE ZEROPREC
C ******************************************************************
C Initialize precipitation indices for original classification
C ******************************************************************
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /Variance/ X2(30,25,26),URules(25,26),RULEX2(25,26),
     X UANOM(30,25,26)
      COMMON /SIZE/ LOND,LATD
      COMMON /RULES/ NTOTR
      DO 1 M=1,4
      DO 2 I=1,30
      DO 3 J=1,NPSTAT
          LWET(I,J,M)=0
          LDRY(I,J,M)=0
          LWET2(I,J,M)=0
          LDRY2(I,J,M)=0
          PS(I,J,M)=0.0
3     CONTINUE
2     CONTINUE
1     CONTINUE
      DO 4 I = 1,NTOTR
      DO 5 J = 1,LOND
      DO 6 K = 1,LATD
      X2(I,J,K) = 0.0
      UANOM(I,J,K) = 0.0
6     CONTINUE
5     CONTINUE
4     CONTINUE
      
      DO 7 I = 1,LATD
      DO 8 J = 1,LOND
      URules(I,J) = 0.0
      RULEX2(I,J) = 0.0
8     CONTINUE
7     CONTINUE
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
      COMMON /DAILY/ CPR(4*7400,25,26),ANOM(4*7400,25,26)
      COMMON /CLASSES/ CDOF(4*7400,30),NTAG
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA      
      COMMON /MONAT/ MON(4*7400)
      CHARACTER*1  CFC, answer
      CHARACTER*100 INF, AVGNAM,SIGNAM,OUTNAM,GCFILE
      character*1 ichr(3000)
      real E(25,26)
      COMMON /AVERA/ AVG(4*366,25,26)
      COMMON /SIGRA/ SIG(4*366,25,26)
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
       
      COMMON /TPERIOD/ NBY,NEY
      COMMON /OBJ_WEI/ WW1,WW2,WW3,WW4,WW5,WW6
      DIMENSION BA(400),BS(400),MA(12)
C
C Read control file
C
C   Reads in the weights of the different parts of the objective func...
      READ (81,*) WW1,WW2,WW3,WW4,WW5,WW6
      OUTNAM='TEST.CLS'
      OPEN(9,FILE=OUTNAM)

10    CONTINUE
C
C     Determine which grid user is interested in
C   Reads in the beginning and end date
      READ(81,*) nby,nbm,nbd,nbh
      READ(81,*) ney,nem,ned,neh
C      write(82,'(8I6)') nby,nbm,nbd,nbh,ney,nem,ned,neh 

C
C Read filename for average pressure data
C
      READ (81,'(A)') AVGNAM
      write(*,*) AVGNAM
      READ (81,'(A)') SIGNAM
      write(*,*) SIGNAM
      OPEN(28,FILE=AVGNAM)
      OPEN(29,FILE=SIGNAM)
	  READ(81,*) NMSEA
	  DO 954 IK=1,NMSEA
	  READ(81,*) WSEA(IK)
	  READ(81,*) (MA(LL),LL=1,12)
	  DO 953 LL=1,12
	    IF(MA(LL).GT.0) THEN
	      MSEA(MA(LL))=IK
	    ENDIF
953     CONTINUE	  
954   CONTINUE	  
C
C Calculate seasonal averages in required resolution
C
C Grid descrition
	  LOND=21
	  LATD=17

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
        J=366
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
c      XW=360.-XW
c      XE=360.-XE
c       IF(XW.GE.360.) XW=XW-360.
c       IF(XE.GE.360.) XE=XE-360.0
C
c       IF(XW.LT.0.0) XW=360.0+XW
c       IF(XE.LT.0.0) XE=360.0+XE
C
C Read initial rules
C
      CALL READRULE(XS,XW,XN,XE)

C
C Convert rule locations to indices
C
      CALL CONVRULE(LOND,LATD,XMINLT,XMINLD,DLAT,DLON)
      

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
         IF(XH.LE.XMIN1) THEN
           INE=I
           XMIN1=XH
         ENDIF
902   CONTINUE
      XMIN1=99999.9
      DO 903 I=1,LATD
         XH=ABS(XS+FLOAT(I-1)*DLAT-XMINLT)
         IF(XH.LT.XMIN1) THEN
           INS=I
           XMIN1=XH
         ENDIF
903   CONTINUE
      XMIN1=99999.9
      DO 904 I=1,LATD
         XH=ABS(XN+FLOAT(I-1)*DLAT-XMINLT)
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
       KYR=1992
	KMO=1
	KDY=1
	IHR=0
       NR=1
322   CALL cdslpZA(NR,KYR,KMO,KDY,IHR,E,NERR)
C      write(*,*)'HEllO'
      IF(NERR.EQ.-1) THEN
        WRITE(*,*) 'PRESSURE DATA READ'
        WRITE(*,*) ' TOTAL ',NTAG,' DAYS'
        WRITE(*,*) 'Last day =', KYR,KMO,KDY
        RETURN
      ENDIF

      if(kyr.lt.nby) goto 765
      if(kyr.eq.nby. and .kmo.lt.nbm ) goto 765
      if(kyr.eq.nby. and .kmo.eq.nbm . and .kdy.lt.nbd) goto 765
c      if(ihr.lt.3) goto 765
607   CONTINUE
      IF(NTAG.GT.0) THEN
c      IF((KYR.NE.IYR).OR.(KMO.NE.IMO).OR.(KDY.NE.IDY)) THEN
       IF(NERR.EQ.-9) THEN
       NTAG=NTAG+1
       WRITE(*,*) KYR,KMO,KDY,NTAG
       write(10,*) KYR,KMO,KDY,NTAG
C       write(*,*)'Hello'
       CPR(NTAG,1,1)=-1.0
       MYR=KYR
       MN=KMO
       MDY=KDY
       MHR = KHR
       CALL CHDATE (MYR,MN,MDY,MHR,KYR,KMO,KDY,KHR,IDELTHR)
       GOTO 322
      ENDIF
      ENDIF
C      write(*,*)KHR
C
C Subrtact mean values and divide with std
      NTAG=NTAG+1
C      write(*,*)NTAG
	CALL JDATE(kYR,kMO,kDY,KHR,JTAG)
C      CALL TESTSMALL(E,INS,INN,INW,INE,KNOT)
	IF(E(1,1).GE.0.0) THEN
      DO 601 I=1,LOND
       DO 602 J=1,LATD
        EHE=(E(I,J)-AVG(JTAG,I,J))/SIG(JTAG,I,J)
        ANOM(NTAG,I,J)=EHE
        IF(ABS(EHE).GT.5.0) THEN
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
c      WRITE(*,*) KYR,KMO,KDY,NTAG
      MON(NTAG)=KMO
2019  format(4I4,I6)
2001  format(11F7.1)
      IDELTHR=6
      MYR=KYR
      MN=KMO
      MDY=KDY
      MHR=KHR
      CALL CHDATE (MYR,MN,MDY,MHR,KYR,KMO,KDY,KHR,IDELTHR)
765   CONTINUE
      if(kyr.gt.ney) THEN
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
      character OUTNAM*50
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
C	Checks whether the grid points are within specified grid
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
      WRITE(60,*)'BUT HERE IS:',NTOTR
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

      SUBROUTINE ICLASS
C ****************************************************************
C Subroutine Class
C ****************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /ATLAG/ TNATL,XMAX2,XMIN2
      COMMON /NOIDAT/ NXI,NYI,NXA,NYA,VALISL,VALAZO
      COMMON /DAILY/ CPR(4*7400,25,26),ANOM(4*7400,25,26)
      COMMON /CLASSES/ CDOF(4*7400,30),NTAG
      COMMON /OBJECTI/ NOBD,NNBD,NCOB
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /MONAT/ MON(4*7400)
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /PRCATL/ LDRYAVG(30,12),LWETAVG(30,12),PSAVG(30,12),DRYLIM
      COMMON /PRCAVG2/ LDRYAVG2(30,12),LWETAVG2(30,12),EXTLIM
      COMMON /CLASSFREQ/ NUMCPS(30),NOMCPS(30)
      COMMON /Variance/ X2(30,25,26),URules(25,26),RULEX2(25,26),
     X UANOM(30,25,26)
      COMMON /SIZE/ LOND,LATD
      REAL DOF(30),ATM(10),AIDB(10),PEX(5)
      REAL AAND(10),AOR(10)
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
C Absolute frequency of old CPs is NOMCPS()
	DO 301 I=1,30
	  NOMCPS(I)=0
301   CONTINUE
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
      DO 9 I=1,NTOTR
        CDOF(ITG,I)=DOF(I)
        IF(DOF(I).GT.DMAX) THEN
          DMAX=DOF(I)
          NCL=I
        ENDIF
9     CONTINUE
C Absolute frequency count
      IF(NCL.LT.30) THEN
        NOMCPS(NCL)=NOMCPS(NCL)+1
	ENDIF
C Classification done - looking for precipitation
C      MS=1
C      IF((MON(ITG).LT.5).OR.(MON(ITG).GT.10)) MS=2
C       MS=MON(ITG)/3+1
       MS=MSEA(MON(ITG))
      DO 201 JST=1,NPSTAT
        IF(PRCVAL(ITG,JST).GT.-0.5) THEN
          IF(PRCVAL(ITG,JST).GT.DRYLIM) THEN
            LWET(NCL,JST,MS)=LWET(NCL,JST,MS)+1
          ELSE
            LDRY(NCL,JST,MS)=LDRY(NCL,JST,MS)+1
          ENDIF
          IF(PRCVAL(ITG,JST).GE.EXTLIM) THEN
            LWET2(NCL,JST,MS)=LWET2(NCL,JST,MS)+1
          ELSE
            LDRY2(NCL,JST,MS)=LDRY2(NCL,JST,MS)+1
          ENDIF
          PS(NCL,JST,MS)=PS(NCL,JST,MS)+PRCVAL(ITG,JST)
        ENDIF
201   CONTINUE
C     This is Where the Variance is stored
      DO 10 I=1,LOND
       DO 11 J = 1,LATD
        UANOM(NCL,I,J) = UANOM(NCL,I,J) + ANOM(ITG,I,J)
        X2(NCL,I,J) = X2(NCL,I,J) + ANOM(ITG,I,J)**2
11      CONTINUE
10      CONTINUE
      
101   CONTINUE
      
C     Average the values (variance within)
      
      DO 12 L = 1,NTOTR          
       DO 13 J = 1,LOND
        DO 14 K = 1,LATD
        IF (NOMCPS(L).GT.0) THEN        	 
        UANOM(L,J,K) = UANOM(L,J,K)/NOMCPS(L)
        URules(J,K) = URules(J,K) + UANOM(L,J,K)
        RULEX2(J,K) = RULEX2(J,K) + UANOM(L,J,K)**2
        X2(L,J,K) = X2(L,J,K)/NOMCPS(L)
        ENDIF
14       CONTINUE
13       CONTINUE
       
12       CONTINUE
C      Variance Between 
        DO 15 J = 1,LOND
         DO 16 K = 1,LATD
         URules(J,K) =URules(J,K)/FLOAT(NTOTR)
         RULEX2(J,K) = RULEX2(J,K)/FLOAT(NTOTR)
16       CONTINUE
15       CONTINUE
       
       
C      WRITE(9,'(30F5.2)')(DOF(I),I=1,NTOTR)
C
C Calculate precipitation indices
C
      RETURN
      END

      SUBROUTINE PRCIDX(PRI1,PRI2,PRI3,PRI4,PRI5,PRI6,PRI7)
C *********************************************************************
C Calculate precipitation indices
C I have changed the original objective funcs
C PRI1 = P(theta|CP)-Pbar - 2.5m
C PRI2 = av(H(CP)) - Hbar
C PRI3 = H(CP)/H - 1
C PRI4 = P(theta|CP)-Pbar - 3.5m
C PRI5 = var within class
C PRI6 = var between classes
C PRI7 = Linked to P(CP|theta) - perhaps use entropy
C *********************************************************************
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /PRCATL/ LDRYAVG(30,12),LWETAVG(30,12),PSAVG(30,12),DRYLIM
      COMMON /PRCAVG2/ LDRYAVG2(30,12),LWETAVG2(30,12),EXTLIM
      COMMON /MONAT/ MON(4*7400)
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA
      COMMON /Variance/ X2(30,25,26),URules(25,26),RULEX2(25,26),
     X UANOM(30,25,26)
      COMMON /SIZE/LOND,LATD
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      DIMENSION PRI1(35,12),PRI2(35,12),PRI3(35,12),PRI4(35,12),
     X PRI7(35,12)
c      REAL PRI5,PRI6,VAR,RuleVar
      PRI5 = 0.0
      PRI6 = 0.0
	  DO 1 M=1,4
      DO 9 I=1,NPSTAT
      PRI1(I,M)=0.0
      PRI2(I,M)=0.0
      PRI3(I,M)=0.0
      PRI4(I,M)=0.0
      PRI7(I,M)=0.0

C Loop over CPs
        MSU=LWETAVG(I,M)+LDRYAVG(I,M)
        IF(MSU.GT.0) THEN
          DO 2 J=1,NTOTR
            MCP=LWET(J,I,M)+LDRY(J,I,M)
            IF(MCP.GT.0) THEN
              HELP=(FLOAT(LWET(J,I,M))/FLOAT(MCP)-FLOAT(LWETAVG(I,M))/
     x          FLOAT(MSU))**2*(FLOAT(MCP)/FLOAT(MSU))
              PRI1(I,M)=PRI1(I,M)+HELP
              MCB=LWET(J,I,M)
              MSB=LWETAVG(I,M)
              IF(MCB.GT.0) THEN
                HELP=(PS(J,I,M)/FLOAT(MCB)-PSAVG(I,M)/FLOAT(MSB))**2
                HELP=HELP*(FLOAT(MCB)/FLOAT(MSB))
                PRI2(I,M)=PRI2(I,M)+HELP
              ENDIF
C Consider completely dry situations
	        IF(PSAVG(I,M).GT.0.01) THEN
               T1=ABS(PS(J,I,M)*FLOAT(MSU)/(FLOAT(MCP)*PSAVG(I,M))-1.0)
c	         T1=ALOG(T1)
c               T1=ABS(PS(J,I,M)*FLOAT(MSB)/(FLOAT(MCB)*PSAVG(I,M))-1.0)
C Concentrate on extremes take squared difference as goal
C
       T1=ABS(PS(J,I,M)*FLOAT(MSU)/(FLOAT(MCP)*PSAVG(I,M))-1.0)**1.25
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	        ELSE
	         T1=1.0
	        ENDIF
              PRI3(I,M)=PRI3(I,M)+T1*FLOAT(MCP)/FLOAT(MSU)
            ENDIF
2         CONTINUE
        ENDIF
        IF(PRI1(I,M).GT.0.0) THEN
          PRI1(I,M)=SQRT(PRI1(I,M))
        ENDIF
        MSU=LWETAVG2(I,M)+LDRYAVG2(I,M)
        IF(MSU.GT.0) THEN
          DO 3 J=1,NTOTR
            MCP=LWET2(J,I,M)+LDRY2(J,I,M)
            IF(MCP.GT.0) THEN
      HELP=(FLOAT(LWET2(J,I,M))/FLOAT(MCP)-FLOAT(LWETAVG2(I,M)
     x          )/FLOAT(MSU))**2*(FLOAT(MCP)/FLOAT(MSU))
              PRI4(I,M)=PRI4(I,M)+HELP
            ENDIF
           IF (LWETAVG2(I,M).GT.0) THEN
               HELP = (FLOAT(LWET2(J,I,M))/FLOAT(LWETAVG2(I,M)))
     x*(FLOAT(MCP)/FLOAT(MSU))
               PRI7(I,M) =  PRI7(I,M) + HELP
         ENDIF
3        CONTINUE
        ENDIF
        PRI4(I,M)=SQRT(PRI4(I,M))
9     CONTINUE
1     CONTINUE
	
C Variance within
      DO 4 L=1,NTOTR
       VAR = 0.0
       DO 5 K = 1,LOND
        DO 6 J = 1,LATD
        HELP = X2(L,K,J) - UANOM(L,K,J)**2
        VAR = VAR + HELP
	 
6       CONTINUE
5       CONTINUE
        
        PRI5 = PRI5 + VAR/FLOAT(LOND*LATD)

4     CONTINUE
       
	RuleVar = 0.0       
C Variance between
      DO 7 K = 1,LOND
       DO 8 J = 1,LATD
       HELP = RULEX2(K,J)-URules(K,J)**2
       RuleVar = RuleVar + HELP
8      CONTINUE
7      CONTINUE
       PRI6 = PRI6 + RuleVar/FLOAT(LOND*LATD)
      
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
	CALL JDATE(IYA,IMO,ITA,IHO,JTAG)
c	IREC=(IYA-1979)*366+JTAG
	IREC=4*(IYA-NRBY)*366+JTAG
      	
C	WRITE(41,*)IYA,IMO,ITA
C	WRITE(41,*)KYA,KMO,KTA
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
	WRITE(82,*) IYA,KYA,IMO,KMO,IREC
	  NERR=-9
	ENDIF
      IF(IMO.NE.KMO) NERR=-9
      IF(ITA.NE.KTA) NERR=-9
	RETURN
999   NERR=-1
     
      RETURN
      END

      SUBROUTINE OPTIMIZE
C ****************************************************************
C Optimize fuzzy definitions
C ****************************************************************
      COMMON /DAILY/ CPR(4*7400,25,26),ANOM(4*7400,25,26)
      COMMON /CLASSES/ CDOF(4*7400,30),NTAG
      COMMON /OBJECTI/ NOBD,NNBD,NCOB
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /RULES2/ NRMAT(30,25,26)
      COMMON /RULEOLD/ NOLDMAT(30,25,26),MOLD
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
      
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
	  COMMON /OBJ_WEI/ WW1,WW2,WW3,WW4,WW5,WW6
      COMMON /CLASSFREQ/ NUMCPS(30),NOMCPS(30)
      DIMENSION DOFNEW(4*7400)
      DIMENSION PRI1(35,12),PRI2(35,12),PRI3(35,12),PRI4(35,12),
     1 PRI7(35,12)
      DOUBLE PRECISION XRAND,RRAND
      WRITE(*,*) NOBD,MOLD,NCOB

            
c      OPEN(57,FILE='protko')
      OPEN(57,FILE='protko')
C Initialize precipitation objective
      RETVO=-10.0

      CALL PRCIDX(PRI1,PRI2,PRI3,PRI4,PRI5,PRI6,PRI7)
      PROBJ=0.0
      PROB2=0.0
      PROB3=0.0
      PROB4=0.0
      PROB5=0.0
      PROB6=0.0
      PROB7=0.0
      DO 101 I=1,NPSTAT
       DO 107 IHL=1,NMSEA
        PROBJ=PROBJ+(0.5-PRI1(I,IHL))*WSEA(IHL)
        PROB2=PROB2+(5.0-PRI2(I,IHL))*WSEA(IHL)
        PROB3=PROB3+(1.0-PRI3(I,IHL))*WSEA(IHL)
        PROB4=PROB4+(0.5-PRI4(I,IHL))*WSEA(IHL)        
        PROB7=PROB7+(0.5-PRI7(I,IHL))*WSEA(IHL)
107    CONTINUE
        WRITE(*,'(12F8.3)') (PRI1(I,LL),LL=1,NMSEA) 
        WRITE(*,'(12F8.3)') (PRI2(I,LL),LL=1,NMSEA) 
        WRITE(*,'(12F8.3)') (PRI3(I,LL),LL=1,NMSEA) 
        WRITE(*,'(12F8.3)') (PRI4(I,LL),LL=1,NMSEA) 
 101   CONTINUE
       PROB5=PROB5+PRI5
       PROB6=PROB6+(5.0-PRI6)
      WRITE(*,*) PROBJ,PROB2,PROB3,PROB4,PROB5,PROB6,PROB7
C Start annealing
      XRAND=1.873
      T1=7.85*FLOAT(NPSTAT)
      T0=T1
        LBRE=INE-INW+1
	  LCRE=INN-INS+1
	KITER=5000

      DO 1 KT=1,500

       IPOS=0
       INEG=0
       DO 2 K=1,KITER
c	  WRITE(*,*) K
11       IRU=INT(RRAND(XRAND)*FLOAT(NTOTR))+1
         IF(IRU.GT.NTOTR) IRU=NTOTR
C3        IX=INT(DRAND(XRAND)*FLOAT(LOND))+1
3        IX=INT(RRAND(XRAND)*FLOAT(LBRE))+1
         IX=IX+INW-1
         IF((IX.LT.INW).OR.(IX.GT.INE)) GOTO 3
4        IY=INT(RRAND(XRAND)*FLOAT(LCRE))+1
         IY=IY+INS-1
         IF((IY.LT.INS).OR.(IY.GT.INN)) GOTO 4

         IBX=INT(RRAND(XRAND)*18.99999)-14
         IF(IBX.LT.0) IBX=0
         IF(IBX.GT.4) IBX=4
       
	   IF(IBX.EQ.NRMAT(IRU,IX,IY)) GOTO 3
         MRULES=NSUM(IRU)
         IF((IBX.GT.0).AND.(NRMAT(IRU,IX,IY).EQ.0)) THEN
           MRULES=NSUM(IRU)+1
         ENDIF
         IF((IBX.EQ.0).AND.(NRMAT(IRU,IX,IY).GT.0)) THEN
           MRULES=NSUM(IRU)-1
         ENDIF
         IF ((MRULES.GT.25).AND.(IBX.NE.0)) GOTO 11
         IF((MRULES.LT.8).AND.(IBX.EQ.0)) GOTO 11
C Check neighbours
C Left
         IF(IX.GT.1) THEN
           IF(IABS(NRMAT(IRU,IX-1,IY)).GT.0) THEN
             MO=IABS(NRMAT(IRU,IX-1,IY)-IBX)
             IF(MO.GT.1) GOTO 11
           ENDIF
         ELSE
c           IF(IABS(NRMAT(IRU,72,IY)).GT.0) THEN
c             MO=IABS(NRMAT(IRU,72,IY)-IBX)
c             IF(MO.GT.1) GOTO 11
c           ENDIF
         ENDIF
C Right
         IF(IX.LT.25) THEN
           IF(IABS(NRMAT(IRU,IX+1,IY)).GT.0) THEN
             MO=IABS(NRMAT(IRU,IX+1,IY)-IBX)
             IF(MO.GT.1) GOTO 11
           ENDIF
         ELSE
           IF(IABS(NRMAT(IRU,1,IY)).GT.0) THEN
             MO=IABS(NRMAT(IRU,1,IY)-IBX)
             IF(MO.GT.1) GOTO 11
           ENDIF
         ENDIF
C Up
         IF(IY.LT.25) THEN
           IF(IABS(NRMAT(IRU,IX,IY+1)).GT.0) THEN
             MO=IABS(NRMAT(IRU,IX,IY+1)-IBX)
             IF(MO.GT.1) GOTO 11
           ENDIF
         ENDIF
C Down
         IF(IY.GT.1) THEN
           IF(IABS(NRMAT(IRU,IX,IY-1)).GT.0) THEN
             MO=IABS(NRMAT(IRU,IX,IY-1)-IBX)
             IF(MO.GT.1) GOTO 11
           ENDIF
         ENDIF
         CALL CHA(DOFNEW,IRU,IX,IY,IBX,IUJO,MNEW,ICOB,RETVN)
         IF(RETVO.LT.-1.) THEN
           RETVO=RETVN
           WRITE(*,*) RETVO
         ENDIF
C Precipitation indices
         CALL PRCIDX(PRI1,PRI2,PRI3,PRI4,PRI5,PRI6,PRI7)
      PRUBJ=0.0
      PRUB2=0.0
      PRUB3=0.0
      PRUB4=0.0
      PRUB5=0.0
      PRUB6=0.0
      PRUB7=0.0
C Sum for stations - and seasons
      DO 102 I=1,NPSTAT
       DO 102 IHL=1,NMSEA
        PRUBJ=PRUBJ+(0.5-PRI1(I,IHL))*WSEA(IHL)
        PRUB2=PRUB2+(5.0-PRI2(I,IHL))*WSEA(IHL)
        PRUB3=PRUB3+(1.0-PRI3(I,IHL))*WSEA(IHL)
        PRUB4=PRUB4+(0.5-PRI4(I,IHL))*WSEA(IHL)
        PRUB7=PRUB7+(0.5-PRI7(I,IHL))*WSEA(IHL)
102   CONTINUE
       PRUB5=PRUB5+PRI5
       
       PRUB6=PRUB6+(5.0-PRI6)
C Calculate objective function value
C  UJO = New value (after change)
C  OBD = Old value
C  Float = Integer wird in single-precision real Zahl umgewandelt
C  Anfangswerte
C        UJO=FLOAT(IUJO)+0.25*FLOAT(MNEW)+FLOAT(ICOB)+1000.0*PRUBJ
C        OBD=FLOAT(NOBD)+0.25*FLOAT(MOLD)+FLOAT(NCOB)+1000.0*PROBJ
c      UJO=0.*FLOAT(IUJO)+0.*FLOAT(MNEW)+0.*FLOAT(ICOB)+100.*PRUBJ+
c     x   0.*PRUB2+100.0*PRUB3+100.0*PRUB4+0.010*RETVN
c      OBD=0.*FLOAT(NOBD)+0.*FLOAT(MOLD)+0.*FLOAT(NCOB)+100.*PROBJ+
c     x   0.*PROB2+100.0*PROB3+100.0*PROB4+0.010*RETVO
      UJO=0.*FLOAT(IUJO)+0.*FLOAT(MNEW)+0.*FLOAT(ICOB)+WW1*PRUBJ+
     x   0.*PRUB2+WW2*PRUB3+WW3*PRUB4+WW4*PRUB5+WW5*PRUB6+WW6*PRUB7
      OBD=0.*FLOAT(NOBD)+0.*FLOAT(MOLD)+0.*FLOAT(NCOB)+WW1*PROBJ+
     x   0.*PROB2+WW2*PROB3+WW3*PROB4+WW4*PROB5+WW5*PROB6+WW6*PROB7


C Multiplier according frequencies
C      CALL EVALNUNO(UMUO,UMUU)
C	 IF(UJO.GT.0.0) THEN
C	   UJO=UJO*UMUU
C	 ELSE
C	   UJO=UJO/UMUU
C	 ENDIF
C	 IF(OBD.GT.0.0) THEN
C	   OBD=OBD*UMUO
C	 ELSE
C	   OBD=OBD/UMUO
C	 ENDIF
c      WRITE(*,*)'old Obj'
c      WRITE(*,*)OBD
c      WRITE(*,*)'new Obj'
c      WRITE(*,*)UJO
        IF(UJO.LT.OBD) THEN
           ISW=1
           IPOS=IPOS+1
         ELSE
           DX=(OBD-UJO)/T0
           PX=RRAND(XRAND)
C ***************
           IF(PX.LT.(EXP(DX))) THEN
             INEG=INEG+1
             ISW=1
           ELSE
             ISW=0
           ENDIF
         ENDIF
         IF(ISW.EQ.1) THEN
          NSUM(IRU)=MRULES
          NOBD=IUJO
          MOLD=MNEW
          NCOB=ICOB
          OBD=UJO
          PROBJ=PRUBJ
          PROB2=PRUB2
          PROB3=PRUB3
          PROB4=PRUB4
          PROB5=PRUB5
          PROB6=PRUB6
          PROB7=PRUB7
          RETVO=RETVN
          NRMAT(IRU,IX,IY)=IBX
          DO 10 ITG=1,NTAG
           CDOF(ITG,IRU)=DOFNEW(ITG)
10        CONTINUE
          DO 202 IGG=1,30
	      NOMCPS(IGG)=NUMCPS(IGG)
202       CONTINUE
         ENDIF
C       WRITE(*,*) IPOS,INEG,NOBD
2      CONTINUE
       WRITE(*,'(3I6,6F7.3,2F8.1)') KT,IPOS,INEG,PROBJ,
     x   PROB3,PROB4,PROB5,PROB6,PROB7,OBD,RETVO
       WRITE(57,'(3I6,6F7.3,2F8.1)') KT,IPOS,INEG,PROBJ,
     x   PROB3,PROB4,PROB5,PROB6,PROB7,OBD,RETVO
       IF((IPOS+INEG).LT.3) GOTO 29
C       T0=T1/ALOG(FLOAT(KT+1))
C
C Correct initial temperature if needed
C

        IF(KT.LE.10) THEN
          ARANY=FLOAT(IPOS+INEG)/FLOAT(KITER)
	    IF(ARANY.GT.0.80) THEN
	      T0=T0*0.8
	    ENDIF
	    IF(ARANY.LT.0.40) THEN
	      T0=T0*1.25
	    ENDIF
        ENDIF
        T0=T0*0.97
c       OPEN(41,FILE='CPACT')
c       OPEN(55,FILE='CPACT.def')
	 CALL WRITERULE

1     CONTINUE
29    CONTINUE
       OPEN(41,FILE='10_CPNEW')
c       OPEN(41,FILE='10_CPNEW')
       DO 21 I=1,NTOTR
       WRITE(41,*) I
        DO 22 J=LATD,1,-1
          WRITE(41,'(80I1)')(NRMAT(I,K,J),K=1,LOND)
22      CONTINUE          
21     CONTINUE
       CLOSE(41)
      RETURN
      END

      SUBROUTINE CHA(DOFNEW,IRU,IX,IY,IBX,IUJO,MNEW,ICOB,RETVN)
C *********************************************************************
C Calculate objective function
C *********************************************************************
      COMMON /DAILY/ CPR(4*7400,25,26),ANOM(4*7400,25,26)
      COMMON /CLASSES/ CDOF(4*7400,30),NTAG
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /RULES2/ NRMAT(30,25,26)
      COMMON /RULEOLD/ NOLDMAT(30,25,26),MOLD
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
      
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /MONAT/ MON(4*7400)
      COMMON /MONSEAS/ WSEA(12),MSEA(12),NMSEA
      
      COMMON /PRCS/ PRCVAL(4*7400,30),NPSTAT,NPCTAG
      COMMON /PRECIDX/ LWET(30,30,12),LDRY(30,30,12),PS(30,30,12)
      COMMON /PREIDX2/ LWET2(30,30,12),LDRY2(30,30,12)
      COMMON /PRCATL/ LDRYAVG(30,12),LWETAVG(30,12),PSAVG(30,12),DRYLIM
      COMMON /PRCAVG2/ LDRYAVG2(30,12),LWETAVG2(30,12),EXTLIM
      COMMON /CLASSFREQ/ NUMCPS(30),NOMCPS(30)
      COMMON /Variance/ X2(30,25,26),URules(25,26),RULEX2(25,26),
     1 UANOM(30,25,26)      
      DIMENSION DOFNEW(4*7400)
      DIMENSION ATM(10),AIDB(10),PEX(4)
C Zero precipitation counters
      CALL ZEROPREC
c      WRITE(*,*)UANOM(1,1,1)
C Absolute frequency of CPs is NUMCPS()
	DO 301 I=1,30
	  NUMCPS(I)=0
301   CONTINUE
C Define exponents
      PEX(1)=1.25
      PEX(2)=5.0
      PEX(3)=5.0
      PEX(4)=1.25
      IUJO=0
      ICOB=0
      DO 1 ITG=1,NTAG
       IF(CPR(ITG,1,1).LT.0.) GOTO 1
       DO 10 KL=1,4
        ATM(KL)=0.0
        AIDB(KL)=0.0
10     CONTINUE
        DO 2 I=1,LOND
C Check if in area
         IF(INW.LT.INE) THEN
           IF((IX.LT.INW).OR.(IX.GT.INE)) GOTO 2
         ELSE
           IF((IX.LT.INW).AND.(IX.GT.INE)) GOTO 2
         ENDIF
C
          DO 3 J=1,LATD
           XNORM=CPR(ITG,I,J)
           IF((I.EQ.IX).AND.(J.EQ.IY)) THEN
            IF(IBX.EQ.0) GOTO 3
            IF(IBX.EQ.1) THEN
              AMB=XMEM1(1.0-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0) ATM(1)=ATM(1)+AMB**PEX(1)
              AIDB(1)=AIDB(1)+1.0
            ENDIF
            IF(IBX.EQ.2) THEN
              AMB=XMEM2(1.-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(2)=ATM(2)+AMB**PEX(2)
              AIDB(2)=AIDB(2)+1.0
            ENDIF
            IF(IBX.EQ.3) THEN
              AMB=XMEM2(XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0) ATM(3)=ATM(3)+AMB**PEX(3)
              AIDB(3)=AIDB(3)+1.0
            ENDIF
            IF(IBX.EQ.4) THEN
              AMB=XMEM1(XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(4)=ATM(4)+AMB**PEX(4)
              AIDB(4)=AIDB(4)+1.0
            ENDIF
           ELSE
            IF(NRMAT(IRU,I,J).EQ.0) GOTO 3
            IF(NRMAT(IRU,I,J).EQ.1) THEN
              AMB=XMEM1(1.0-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0) ATM(1)=ATM(1)+AMB**PEX(1)
              AIDB(1)=AIDB(1)+1.0
            ENDIF
            IF(NRMAT(IRU,I,J).EQ.2) THEN
              AMB=XMEM2(1.-XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(2)=ATM(2)+AMB**PEX(2)
              AIDB(2)=AIDB(2)+1.0
            ENDIF
            IF(NRMAT(IRU,I,J).EQ.3) THEN
              AMB=XMEM2(XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0) ATM(3)=ATM(3)+AMB**PEX(3)
              AIDB(3)=AIDB(3)+1.0
            ENDIF
            IF(NRMAT(IRU,I,J).EQ.4) THEN
              AMB=XMEM1(XNORM,0)
              AMB=1.0-AMB
              IF(AMB.GT.0.0)  ATM(4)=ATM(4)+AMB**PEX(4)
              AIDB(4)=AIDB(4)+1.0
            ENDIF
           ENDIF
3         CONTINUE
2       CONTINUE
C Evaluate categories
        AMULT=1.0
        AEXP=1.0
        DO 5 KL=1,4
           IF(AIDB(KL).GT.0.0) THEN
             IF(ATM(KL).GT.0.0) THEN
               ATM(KL)=1.0-(ATM(KL)/AIDB(KL))**(1./PEX(KL))
             ELSE
                ATM(KL)=1.0
             ENDIF
           ELSE
             ATM(KL)=1.0
             AEXP=4./3.
           ENDIF
           AMULT=AMULT*ATM(KL)
5       CONTINUE
        IF(AMULT.GT.0.0) AMULT=AMULT**AEXP
        DOFNEW(ITG)=AMULT
        NCL=NTOTR+1
        DMAX=0.0
        DO 6 K=1,NTOTR
         IF(K.NE.IRU) THEN
           IF(CDOF(ITG,K).GT.DMAX) THEN
             NCL=K
             DMAX=CDOF(ITG,K)
           ENDIF
         ELSE
           IF(AMULT.GT.DMAX) THEN
             NCL=K
             DMAX=AMULT
           ENDIF
         ENDIF
6       CONTINUE
C Absolute frequencies
      IF(NCL.LT.30) THEN
        NUMCPS(NCL)=NUMCPS(NCL)+1
	ENDIF
C
C Classification done - looking for precipitation
C
c       MS=1
c       IF((MON(ITG).LT.5).OR.(MON(ITG).GT.10)) MS=2
C        MS=MON(ITG)/3+1
         MS=MSEA(MON(ITG))       
      DO 201 JST=1,NPSTAT
        IF(PRCVAL(ITG,JST).GT.-0.5) THEN
          IF(PRCVAL(ITG,JST).GT.DRYLIM) THEN
            LWET(NCL,JST,MS)=LWET(NCL,JST,MS)+1
          ELSE
            LDRY(NCL,JST,MS)=LDRY(NCL,JST,MS)+1
          ENDIF
          IF(PRCVAL(ITG,JST).GT.EXTLIM) THEN
            LWET2(NCL,JST,MS)=LWET2(NCL,JST,MS)+1
          ELSE
            LDRY2(NCL,JST,MS)=LDRY2(NCL,JST,MS)+1
          ENDIF
          PS(NCL,JST,MS)=PS(NCL,JST,MS)+PRCVAL(ITG,JST)
        ENDIF
201   CONTINUE
C     This is Where the Variance is stored
      DO 7 I=1,LOND
       DO 8 J = 1,LATD
        UANOM(NCL,I,J) = UANOM(NCL,I,J) + ANOM(ITG,I,J)
        X2(NCL,I,J) = X2(NCL,I,J) + ANOM(ITG,I,J)**2
8      CONTINUE
7      CONTINUE
1     CONTINUE
C     Average the values (variance within)
      DO 11 I = 1,NTOTR
       DO 12 J = 1,LOND
        DO 13 K = 1,LATD
        IF(NUMCPS(I).GT.0)THEN
        UANOM(I,J,K) = UANOM(I,J,K)/FLOAT(NUMCPS(I))
        URules(J,K) = URules(J,K) + UANOM(I,J,K)
        RULEX2(J,K) = RULEX2(J,K) + UANOM(I,J,K)**2
        X2(I,J,K) = X2(I,J,K)/FLOAT(NUMCPS(I))
        ENDIF
13       CONTINUE
12       CONTINUE
11       CONTINUE
C       WRITE(*,*)(NUMCPS(I),I=1,NTOTR)
C      Variance Between 
        DO 14 J = 1,LOND
         DO 15 K = 1,LATD
         URules(J,K) =URules(J,K)/FLOAT(NTOTR)
         RULEX2(J,K) = RULEX2(J,K)/FLOAT(NTOTR)
15       CONTINUE
14       CONTINUE
C        WRITE(*,*)UANOM(1,1,1)
        MNEW=MOLD
C Alternative 1 IF sequence
         GOTO 87
        IF(IBX.EQ.NOLDMAT(IRU,IX,IY)) THEN
          IF(IBX.EQ.0) THEN
            MNEW=MOLD-1
          ENDIF
          IF(IBX.NE.0) THEN
            IF(NRMAT(IRU,IX,IY).EQ.0) THEN
              MNEW=MOLD-1
            ELSE
              IF(NRMAT(IRU,IX,IY).GT.IBX) THEN
                MNEW=MNEW-NRMAT(IRU,IX,IY)+IBX
              ELSE
                MNEW=MNEW+NRMAT(IRU,IX,IY)-IBX
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IF(0.EQ.NOLDMAT(IRU,IX,IY)) THEN
            IF(NRMAT(IRU,IX,IY).EQ.0) THEN
              MNEW=MOLD+1
            ENDIF
          ELSE
            IF(IBX.EQ.0) THEN
              IF(NRMAT(IRU,IX,IY).GE.NOLDMAT(IRU,IX,IY)) THEN
                MNEW=MOLD-NRMAT(IRU,IX,IY)+NOLDMAT(IRU,IX,IY)+1
              ELSE
                MNEW=MOLD+NRMAT(IRU,IX,IY)-NOLDMAT(IRU,IX,IY)+1
              ENDIF
            ELSE
              IF(NRMAT(IRU,IX,IY).EQ.0) THEN
                MNEW=MOLD-1+IABS(NOLDMAT(IRU,IX,IY)-IBX)
              ELSE
                MNEW=MOLD+IABS(NOLDMAT(IRU,IX,IY)-IBX)-
     x               IABS(NOLDMAT(IRU,IX,IY)-NRMAT(IRU,IX,IY))
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C Alternative 2 multiplication sequence
87    CONTINUE
      INO=0
      IF(NOLDMAT(IRU,IX,IY).GT.0) INO=1
      INU=0
      IF(NRMAT(IRU,IX,IY).GT.0) INU=1
      INA=0
      IF(IBX.GT.0) INA=1
      MNEW=MOLD-IABS(INO-INU)+IABS(INO-INA)
      MNEW=MNEW-INO*INU*IABS(NOLDMAT(IRU,IX,IY)-NRMAT(IRU,IX,IY))+
     x     INO*INA*IABS(NOLDMAT(IRU,IX,IY)-IBX)
      CALL RULEVAL(RETVN,IRU,IX,IY,IBX)
      RETURN
      END

      SUBROUTINE RULEVAL(RETVN,IRU,IX,IY,IBX)
C *********************************************************************
C Subroutine to evaluate rule - connectivity etc
C *********************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /RULES2/ NRMAT(30,25,26)
      COMMON /RULEOLD/ NOLDMAT(30,25,26),MOLD
      COMMON /AREA/ INS,INW,INN,INE
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
      DIMENSION MSUM(30),MTYPE(30,40),ML1(30,40),ML2(30,40),MDB(30,4)
      DIMENSION AVGPO1(30,4),AVGPO2(30,4),DISPO1(30,4),DISPO2(30,4)
      RETVN=0.0
      IKEEP=NRMAT(IRU,IX,IY)
      NRMAT(IRU,IX,IY)=IBX
C  4 - 1 distance
      DO 1 IR=1,NTOTR
       RETVL=0.0
       MSUM(IR)=0
       DO 11 I=1,LOND
C Check for area
        DO 21 J=1,LATD
         IF(NRMAT(IR,I,J).GT.0) THEN
           MSUM(IR)=MSUM(IR)+1
           MTYPE(IR,MSUM(IR))=NRMAT(IR,I,J)
           ML1(IR,MSUM(IR))=I
           ML2(IR,MSUM(IR))=J
         ENDIF
21      CONTINUE
11     CONTINUE
1     CONTINUE
      DO 2 IR=1,NTOTR
        DO 3 L=1,4
          AVGPO1(IR,L)=0
          AVGPO2(IR,L)=0
          DISPO1(IR,L)=0
          DISPO2(IR,L)=0
          MDB(IR,L)=0
3       CONTINUE
        DO 4 K=1,MSUM(IR)
          AVGPO1(IR,MTYPE(IR,K))=AVGPO1(IR,MTYPE(IR,K))+FLOAT(ML1(IR,K))
          AVGPO2(IR,MTYPE(IR,K))=AVGPO2(IR,MTYPE(IR,K))+FLOAT(ML2(IR,K))
          DISPO1(IR,MTYPE(IR,K))=DISPO1(IR,MTYPE(IR,K))
     x       +FLOAT(ML1(IR,K)**2)
          DISPO2(IR,MTYPE(IR,K))=DISPO2(IR,MTYPE(IR,K))
     x       +FLOAT(ML2(IR,K)**2)
          MDB(IR,MTYPE(IR,K))=MDB(IR,MTYPE(IR,K))+1
4       CONTINUE
        NHI=0
        NORMAX=0
        DO 5 L=1,4
          MD=MDB(IR,L)
          IF(MD.GT.NORMAX) NORMAX=MD
          IF(MD.GT.0) THEN
            RHELP1=AVGPO1(IR,L)/FLOAT(MD)
            RHELP1=DISPO1(IR,L)/FLOAT(MD)-RHELP1*RHELP1
            RHELP2=AVGPO2(IR,L)/FLOAT(MD)
            RHELP2=DISPO2(IR,L)/FLOAT(MD)-RHELP2*RHELP2
            RTOT=RHELP1+RHELP2
            IF((L.LE.2).AND.(RTOT.GT.10.)) THEN
              RETVL=RETVL+RHELP1+RHELP2
            ENDIF
            IF((L.GT.2).AND.(RTOT.GT.5.)) THEN
              RETVL=RETVL+RHELP1+RHELP2
            ENDIF
          ELSE
            NHI=NHI+1
          ENDIF
5       CONTINUE
      IF(RETVL.GT.10.) THEN
        RETVN=RETVN+RETVL
      ENDIF
      IF(NORMAX.GT.10) THEN
        RETVN=RETVN+FLOAT(NORMAX-10)*500.0
      ENDIF
      IF(NHI.GT.0) THEN
        RETVN=RETVN+30000.0
      ENDIF
      IF(NHI.GT.1) THEN
        RETVN=RETVN+90000.0
      ENDIF
c      WRITE(*,*) IR,RETVN
2     CONTINUE
c      READ(*,*) TETU
      NRMAT(IRU,IX,IY)=IKEEP
      RETURN
      END


      FUNCTION DRAND(IX)
C *********************************************************************
C
C PORTABLE RANDOM NUMBER GENERATOR 
C  REFERENCE: SCRAGE,L.,"A MORE PORTABLE FORTRAN RANDOM
C                        NUMBER GENERATOR."
C        IN  ACM TRANS. ON MATH. SOFTWARE, V.5.,N.2.,1979
C
C  IX = IX * A MOD P  RECURSION USED
C
C *********************************************************************
      DOUBLE PRECISION A,P,IX,B15,B16,XHI,XALO,LEFTLO,FHI,K
      A    = 16807.D0
      B15  = 32768.D0
      B16  = 65536.D0
      P    = 2147483647.D0
C
C GET 15 HIGH ORDER BITS OF IX
C
      XHI= IX/B16
      XHI= XHI-DMOD(XHI,1.D0)
C
C GET 16 LOWER BITS OF IX AND LOW PRODUCT
C
      XALO = (IX-XHI*B16)*A
C
C GET 15 HIGH ORDER BITS OF LOW PRODUCT
C
      LEFTLO = XALO/B16
      LEFTLO = LEFTLO - DMOD(LEFTLO,1.D0)
C
C FORM THE 31 HIGHEST BITS OF FULL PRODUCT
C
      FHI = XHI*A+LEFTLO
C
C GET OVERFLOW PAST 31ST BIT OF FULL PRODUCT
C
      K = FHI/B15
      K = K-DMOD(K,1.D0)
C
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P. THE PARANTHESES ARE
C ESSENTIAL
C
      IX = (((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
C
C ADD P BACK IF NECESSARY
C
      IF(IX.LT.0.D0) IX=IX+P
C
C MULTIPLY BY (1/(2**31-1))
C
      DRAND = IX*4.656612875D-10
      RETURN
      END

	SUBROUTINE WRITERULE
C ****************************************************************************
C Write rules in DEF format
C ****************************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /RULES2/ NRMAT(30,25,26)
      COMMON /RULEOLD/ NOLDMAT(30,25,26),MOLD
      COMMON /SIZE/ LOND,LATD,XMINLT,XMINLD
	COMMON /REALRULE/ XLOC1(30,40),XLOC2(30,40)
      REAL EE(100),FF(100)
      
      OPEN(41,FILE='CPACT')
      OPEN(55,FILE='CPACT.def')
      
      DLON=2.50
      DLAT=2.50

      DO 1 I=1,LOND
         EE(I)=-FLOAT(I-1)*DLON+XMINLD
c	   IF(EE(I).LT.0.0) EE(I)=360.0+EE(I)
1     CONTINUE
      DO 2 J=1,LATD
         FF(J)=-FLOAT(J-1)*DLAT+XMINLT
2     CONTINUE

      DO 3 I=1,NTOTR
	KL=0
        DO 4 J=1,4
	    DO 5 J1=1,LATD
	      DO 6 I1=1,LOND
	        IF(NRMAT(I,I1,J1).EQ.J) THEN
	          KL=KL+1
	          NTYPE(I,KL)=J
                XLOC2(I,KL)=EE(I1)
	          XLOC1(I,KL)=FF(J1)
	        ENDIF
6           CONTINUE
5         CONTINUE
4       CONTINUE
C WRITE Rule
     	KUOLD=0
      DO 7 IL=1,KL
	IF(KUOLD.NE.NTYPE(I,IL)) THEN
	  KU=1
	  KUOLD=NTYPE(I,IL)
	ELSE
	  KU=KU+1
	ENDIF
	WRITE(55,'(3I3,2F10.2)') I,NTYPE(I,IL),KU,
     x    XLOC1(I,IL),XLOC2(I,IL)
7     CONTINUE
3     CONTINUE
	 CLOSE(55)
       DO 91 I=1,NTOTR
       WRITE(41,*) I
        DO 92 J=1,LATD
         WRITE(41,'(80I2)')(NRMAT(I,K,J),K=1,LOND)
92      CONTINUE
91     CONTINUE
       CLOSE(41)

      RETURN
	END


      SUBROUTINE EVALNUNO(UMUO,UMUU)
C **************************************************************************
C * Evaluate CP frequencies to have sufficien number of cases
C **************************************************************************
      COMMON /RULES/ NSUM(30),NTYPE(30,40),NLOC1(30,40),NLOC2(30,40),
     1  NUTY(30),NTOTR
      COMMON /CLASSFREQ/ NUMCPS(30),NOMCPS(30)
      COMMON /CLASSES/ CDOF(4*7400,30),NTAG
      UMUU=1.0
	UMUO=1.0
c	XLIMIT=0.025
	XLIMIT=0.05
	DO 1 I=1,NTOTR
	  HUF=FLOAT(NUMCPS(I))/FLOAT(NTAG)
        IF(HUF.LT.XLIMIT) THEN
	    DADD=(XLIMIT-HUF)/XLIMIT
	  ELSE
	    DADD=0.0
	  ENDIF
	  UMUU=UMUU*(1.+4.*DADD)
	  HOF=FLOAT(NOMCPS(I))/FLOAT(NTAG)
        IF(HOF.LT.XLIMIT) THEN
	    DADD=(XLIMIT-HOF)/XLIMIT
	  ELSE
	    DADD=0.0
	  ENDIF
	  UMUO=UMUO*(1.+4.*DADD)
1     CONTINUE
      RETURN
	END

	SUBROUTINE TESTSMALL(E1,INS,INN,INW,INE,KNOT)
C ****************************************************************
C Reduce matrix
C ****************************************************************
      DIMENSION E1(72,15)
	K1=0
	KNOT=0
	DO 1 I=INW,72
        K1=K1+1
	  K2=0
	  DO 2 J=INS,INN
          K2=K2+1
	    IF(E1(I,J).LT.100.0) KNOT=1
2       CONTINUE
1     CONTINUE        
	DO 3 I=1,INE
        K1=K1+1
	  K2=0
	  DO 4 J=INS,INN
          K2=K2+1
	    IF(E1(I,J).LT.100.0) KNOT=1
4       CONTINUE
3     CONTINUE 
      IF(E1(1,1).LE.0.0) E1(1,1)=0.0
      IF(KNOT.EQ.1) E1(1,1)=-9.0       
      RETURN
	END

      DOUBLE PRECISION FUNCTION RRAND(XXX)
C ********************************************************
C Random number generator using Random.org files
C ********************************************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /NRANDS/ RR(10000),NIND,ISACT
      DIMENSION IRA(10000)
      CHARACTER FILN*120
66    CONTINUE      
      IF(ISACT.NE.99) THEN        
        OPEN(99,FILE='/home/justin/CP_Analysis/data/RANDOMORG/Rnames2')
        NIND=99999
        ISACT=99
      ENDIF
      NIND=NIND+1
      IF(NIND.GT.10000) THEN
        READ(99,'(A)',END=99) FILN
        FILN='/home/justin/CP_Analysis/data/RANDOMORG/'//FILN
        OPEN(98,FILE=FILN)
        DO 1 I=1,2000
          K=(I-1)*5+1
          READ(98,*) (IRA(J),J=K,K+4)
1       CONTINUE
        CLOSE(98)
        DO 2 I=1,10000
          RR(I)=FLOAT(IRA(I))/1000000.0
2       CONTINUE
        NIND=1    
        IF(NIND.GT.20000) NIND=5        
      ENDIF
      RRAND=RR(NIND)
      RETURN
99    CLOSE(99)
      ISACT=0
      GOTO 66      
      END        
