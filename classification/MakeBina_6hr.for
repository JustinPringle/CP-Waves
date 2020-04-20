      DIMENSION AA(17,21),BB(357)
      DIMENSION AV(4*366,357),STD(4*366,357),NTOT(4*366)
      CHARACTER DUMMY*80

      ISEL=1
      IF(ISEL.EQ.1) THEN
        IYA=1979
        IMO=1
        ITA=1
        IHr = 0
        NTOT=0

       OPEN(20,
     xFile='../data/ECMWF/ERAI/ERAI-700hPa-geopotential-height.txt')
      OPEN(30,File='../data/ECMWF/ERAI/ERAI_7bina.dat'
     x,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=1500)
        OPEN(31,File='../data/ECMWF/ERAI/ZaERAIAvg.dat')
        OPEN(32,File='../data/ECMWF/ERAI/ZaERAIStd.dat')
      ENDIF

      IF(ISEL.EQ.2) THEN
        IYA=1958
        IMO=1
        ITA=1
        NTOT=0

       OPEN(20,
     xFile='../data/ECMWF/ERA40/ERA40-700hPa-geopotential-height.txt')
      OPEN(30,File='../data/ECMWF/ERA40/ERA40_7bina.dat'
     x,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=1500)

        OPEN(31,File='../data/ECMWF/ERA40/ZaERA40Avg.dat')
        OPEN(32,File='../data/ECMWF/ERA40/ZaERA40Std.dat')
      ENDIF
      IBA=IYA
1     READ(20,1001,END=99) NY,NM,NT        
1001  FORMAT(I4,1X,I2,1X,I2)
      READ(20,*) (BB(K),K=1,357)
      CALL JDATE(IYA,IMO,ITA,IHr,JTAG)
c there are 4 readings per day at 6 hrly intervals - 0,6,12,18
      IREC=4*(IYA-IBA)*366+JTAG
      IF(NY.LT.IBA) GOTO 101
      WRITE(87,*) NY,NM,NT,IYA,IMO,ITA,IHr,JTAG,IREC
      WRITE(30,REC=IREC) IYA,IMO,ITA,(BB(L1),L1=1,357)
      IF(BB(1).LT.0.0) THEN
        WRITE(78,*) IYA,IMO,ITA
        MHR=0
        IDELTHR=6
        CALL CHDATE (IYA,IMO,ITA,MHR,MYRCH,MNCH,MDYCH,MHRCH,IDELTHR)
        IYA=MYRCH
        IMO=MNCH
        ITA=MDYCH
        READ(20,'(A)') DUMMY
        READ(20,'(A)') DUMMY
        GOTO 1
      ENDIF
c     changed JTAG to ITA
      NTOT(JTAG)=NTOT(JTAG)+1
      DO 2 J=1,357
        AV(JTAG,J)=AV(JTAG,J)+BB(J)
2     CONTINUE
      MHR=IHr
      IDELTHR=6
      CALL CHDATE (IYA,IMO,ITA,MHR,MYRCH,MNCH,MDYCH,MHRCH,IDELTHR)
      IYA=MYRCH
      IMO=MNCH
      ITA=MDYCH
      IHr=MHRCH
101      READ(20,'(A)') DUMMY
      READ(20,'(A)') DUMMY
      GOTO 1
99      CLOSE(20)

      DO 3 JTAG=1,4*365
      DO 3 J=1,357
        AV(JTAG,J)=AV(JTAG,J)/FLOAT(NTOT(JTAG))
3     CONTINUE

      IF(ISEL.EQ.1) THEN
        IYA=1979
        IMO=1
        ITA=1
        IHr=0
        NTOT=0

       OPEN(20,
     xFile='../data/ECMWF/ERAI/ERAI-700hPa-geopotential-height.txt')
      ENDIF
      IF(ISEL.EQ.2) THEN
        IYA=1958
        IMO=1
        ITA=1
        NTOT=0

      OPEN(20,
     xFile='../data/ECMWF/ERA40/ERA40-700hPa-geopotential-height.txt')
      ENDIF
41    READ(20,1001,END=999) NY,NM,NT
      READ(20,*) (BB(K),K=1,357)
      CALL JDATE(IYA,IMO,ITA,JTAG)
      IF(BB(1).LT.0.0) THEN
        WRITE(78,*) IYA,IMO,ITA
        MHR=0
        IDELTHR=6
        CALL CHDATE (IYA,IMO,ITA,MHR,MYRCH,MNCH,MDYCH,MHRCH,IDELTHR)
        IYA=MYRCH
        IMO=MNCH
        ITA=MDYCH
        READ(20,'(A)') DUMMY
        READ(20,'(A)') DUMMY
        GOTO 41
      ENDIF
      NTOT(JTAG)=NTOT(JTAG)+1
      DO 21 J=1,357
        STD(JTAG,J)=STD(JTAG,J)+(BB(J)-AV(JTAG,J))**2
21     CONTINUE
      MHR=IHr
      IDELTHR=6
      CALL CHDATE (IYA,IMO,ITA,MHR,MYRCH,MNCH,MDYCH,MHRCH,IDELTHR)
      IYA=MYRCH
      IMO=MNCH
      ITA=MDYCH
      IHr=MHRCH
      READ(20,'(A)') DUMMY
      READ(20,'(A)') DUMMY
      GOTO 41
999      CLOSE(20)


      WRITE(*,*) IREC
      DO 31 JTAG=1,4*365
      DO 31 J=1,357
        STD(JTAG,J)=STD(JTAG,J)/FLOAT(NTOT(JTAG))
        STD(JTAG,J)=SQRT(STD(JTAG,J))
31    CONTINUE
      LAG=4*15
      NDB=357
      CALL EFILTER(AV,NDB,LAG)
      CALL EFILTER(STD,NDB,LAG)
      DO 5 I=1,4*365
        WRITE(31,'(I6)') I
        WRITE(31,'(21F10.3)') (AV(I,J),J=1,NDB)
        WRITE(32,'(I6)') I
        WRITE(32,'(21F10.3)') (STD(I,J),J=1,NDB)
5     CONTINUE
      STOP
      END

      SUBROUTINE EFILTER(AV,NDB,LAG)
C ********************************************************************
C Calculate smoothes time series
C ********************************************************************
      DIMENSION AV(4*366,357),AVA(4*366,357),WSUM(4*365)
      AVA=0.0
      DO 1 I=1,NDB
        WSUM=0.0
        DO 2 J=1,4*365
          DO 3 K=J-LAG,J+LAG
            WW=1.-ABS(FLOAT(K-J))/FLOAT(LAG)
            K1=K
            IF(K.LT.1) K1=4*365+K
            IF(K.GT.365) K1=K-4*365
            AVA(J,I)=AVA(J,I)+AV(K1,I)*WW
            WSUM(J)=WSUM(J)+WW
3         CONTINUE
2       CONTINUE
        DO 4 J=1,4*365
          AVA(J,I)=AVA(J,I)/WSUM(J)
4       CONTINUE
1     CONTINUE
      AV=AVA
      RETURN
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
      IF(NHr.eq.0)THEN
        Time = 0
      ELSEIF(NHr.eq.6)THEN
        Time = 1
      ELSEIF(NHr.eq.12)THEN
        Time = 2
      ELSEIF(NHr.eq.18)THEN
        Time = 3
      ENDIF
      JTAG=MSLEN(NM)+ 4*(ND-1)+ND + Time
      RETURN
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

