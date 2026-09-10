      PROGRAM DIHEDKIN
C     Analyze peptide kinetics based on dihedral angle transitions
C     Dihedral types
C     IT=1 - (PHI,PSI) pair - use alpha and beta basins
C     IT=2 - 3-modal, g-,t,g+
C     IT=3 - 2-modal, a,b
C     INPUT:
C       1. N1,N2,N3 - number of cordinates of types 1,2,3
C                     (e.g. N1=2 for one (phi,psi) pair)
C       2. DT,TCUT - time step per ps and STRAN cutoff
C       3. file-name - data file name
C       4. IMASK(K),K=1,NS  - 0/1 to select POP and RATE output
C     OUTPUT: full data to popk.dat
C             POP(K),KPATT(K) , K=1,NSTATE
C             RATE(K,L)       , K,L=1,NSTATE
C           :  masked data to popkm.dat
C             POPO(K),KPATTO(K) , K=1,NSO
C             RATEO(K,L)       , K,L=1,NSO
C
C     K. Kuczera  Lawrence, KS.  18-Apr-2017
C     MAXT - number of columns in data file = no of time points
C     MAXC - number of columns in data file = no of dihedrals
C     MAXS - number of microstates
C
      PARAMETER(MAXT=1000000,MAXC=20,MAXS=1000)
      REAL*4 TIME(MAXT),DIHE(MAXT,MAXC),RESTIM(MAXS),REST(MAXS)
      REAL*4 RATE(MAXS,MAXS),POP(MAXS)
      INTEGER IS(MAXC),MS(MAXT),KPOW(MAXC),KPATT(MAXS,MAXC)
      INTEGER IFOUND(MAXS),IPATT(MAXS,MAXC),NTRAN(MAXS,MAXS)
      INTEGER IMASK(MAXS),KPATTO(MAXS,MAXC)
      INTEGER ICON(MAXS,MAXS),JCON(MAXS,MAXS),JJC(MAXS)
      REAL*4  STRAN(MAXS,MAXS),RATEO(MAXS,MAXS),POPO(MAXS)
      REAL*4  RATC(MAXS,MAXS)
      INTEGER I,J,N1,N2,N3,NT,NH1
      INTEGER I1,I2,I3,IH1,IH2,IH3
      CHARACTER*72 FILE1
C
C...Define parameters for coordinate classes
C...This is for NATA/CHARMM36/TIP4P
C...A=alpha, B=C7eq, C=beta
      PHIA =  -90.0
      PSIA =  -50.0
      DPPA =   20.0
      PHIB = -100.0
      PSIB =   70.0
      DPPB =   20.0
      PHIC = -100.0
      PSIC =  150.0
      DPPC =   30.0
      DI3A =  -60.0
      DI3B =  180.0
      DI3C =   60.0
      DIS3 =   30.0
      DI2A =  -95.0
      DI2B =   95.0
      DIS2 =   30.0
C
      WRITE(*,*) '...testing Phi-Psi distances RAB,RAC,RBC'
      DX = ABS(PHIB-PHIA)
      IF(DX.GT.180.0) DX = DX-360.0
      DY = ABS(PSIB-PSIA)
      IF(DY.GT.180.0) DY = DY-360.0
      RAB = SQRT(DX**2 + DY**2)
      WRITE(*,'(A,3F10.2)') '...dphi,dpsi,RAB=',DX,DY,RAB
      DX = ABS(PHIC-PHIA)
      IF(DX.GT.180.0) DX = DX-360.0
      DY = ABS(PSIC-PSIA)
      IF(DY.GT.180.0) DY = DY-360.0
      RAC = SQRT(DX**2 + DY**2)
      WRITE(*,'(A,3F10.2)') '...dphi,dpsi,RAC=',DX,DY,RAC
      DX = ABS(PHIB-PHIC)
      IF(DX.GT.180.0) DX = DX-360.0
      DY = ABS(PSIB-PSIC)
      IF(DY.GT.180.0) DY = DY-360.0
      RBC = SQRT(DX**2 + DY**2)
      WRITE(*,'(A,3F10.2)') '...dphi,dpsi,RBC=',DX,DY,RBC
C... read basic counts
      READ(5,*) N1,N2,N3
      WRITE(*,*) ' >>>>>>>>>>>>Dihedral kinetics analysis<<<<<<<<<<<<<'
      READ(5,*) DT,TCUT
      WRITE(*,'(A,3I5,2F8.4)') ' ...N1,N2,N3,DT,TCUT=',N1,N2,N3,DT,TCUT
C
C... calculate indices NT - tot no of dihed data cols
C    NC - tot no of coordinates
C    I1 I2 I3 = high coor ranges for types 1,2,3
C    IH1 IH2 IH3 = high data ranges for types 1,2,3 - J in DIHE(I,J)
      NT=N1+N2+N3
      N1H=N1/2
      NC=N1H+N2+N3
C
      I1=N1H
      I2=I1+N2
      I3=NC
C
      IH1=N1
      IH2=IH1+N2
      IH3=IH2+N3
      WRITE(*,*) '...I1,I2,I3,IH1,IH2,IH3=',I1,I2,I3,IH1,IH2,IH3
      WRITE(*,*) '...Number of cordinates NC =',NC
C
C... Calculate power coefficients to generate decimal state numbering
C... Here PP = 3 , 3-fold = 3 and 2-fold = 2 state
C      
      NSTATE=1
      K=1
      KPOW(1)=1
      IF(N1H.GT.0) THEN
        DO I=1,N1H
            K=K+1
            KPOW(K)=KPOW(K-1)*3
          NSTATE=NSTATE*3
        END DO
      END IF
      IF(N2.GT.0) THEN
        DO I=1,N2
            K=K+1
            KPOW(K)=KPOW(K-1)*3
          NSTATE=NSTATE*3
        END DO
      END IF
      IF(N3.GT.0) THEN
        DO I=1,N3
            K=K+1
            KPOW(K)=KPOW(K-1)*2
          NSTATE=NSTATE*2
        END DO
      END IF
C...test it
       WRITE(*,*) '...Total microstates =',NSTATE
       WRITE(*,*) '...The KPOWs...='
       WRITE(*,'(20I8)') (KPOW(K),K=1,NC)
       IF(NSTATE.GT.MAXS) STOP ' Increase MAXS'
C
       DO K=1,NC
         RESTIM(K) = 0.0
         REST(K) = 0.0
       END DO
C...initialize IPATT paterns at 9
C   the IPATT values will be found from data
       DO K=1,NSTATE
         DO J=1,NC
           IPATT(K,J)=9
         END DO
         IFOUND(K)=0
       END DO
C
C...For computing fun, generate patterns from numbers
C...the KPATT values are calculated from algorithm
C 
      DO I=1,NSTATE
        II = I-1
        NN = NC
        DO J=1,NC-1
          IS(NN) = INT(II/KPOW(NN))
          II = II - IS(NN)*KPOW(NN)
          NN = NN-1
        END DO
        IS(1) = II
        DO J=1,NC
          KPATT(I,J) = IS(J)
        END DO
      END DO
C...test
       WRITE(*,*) '...Patterns IPATT and KPATT'
       DO K=1,NSTATE
         WRITE(*,'(I4,4X,3I4,4X,3I4)') K,(IPATT(K,J),J=1,NC),
     $         (KPATT(K,J),J=1,NC)
       END DO
C
	READ(5,901) FILE1
	OPEN (UNIT=1,FILE=FILE1,FORM='FORMATTED',STATUS='OLD')
        WRITE(6,*) ' ...Data will be read from ',FILE1
901	FORMAT(A40)
C
C...read in output mask for POP and RATE
      READ(5,*) (IMASK(I),I=1,NSTATE)
      WRITE(*,*) '...output IMASK:'
      WRITE(*,*) (IMASK(I),I=1,NSTATE)
C
C
C ...Read in data = dihedral angle time series
      ND=1
  10  READ(1,900,END=20,ERR=20) TIME(ND),(DIHE(ND,I),I=1,NT)
        ND=ND+1
        GO TO 10
  20  CONTINUE
      ND=ND-1
      WRITE(*,903) ND
 900  FORMAT(F10.1,4F9.3)
 903  FORMAT(2X,'No. OF DATA READ N=',I12)
C
C ...Test read
      WRITE(*,*) '... Testing data read ...'
      DO I=1,10
        WRITE(*,904) TIME(I),(DIHE(I,J),J=1,NT)

      END DO
 904  FORMAT(1X,F10.1,50F8.2)
C
C... Analyze data:  
C    Algorithm: for each coord determine micro state = min or transition
C               then find global state
C               and assign residence times and analyze transitions
C               
C ...Loop over time steps
      ITOTS=0
      ITRANS=0
      MSIP=0
      MSIPIN=0
      TTRAN=0.0
      DO J=1,NC
       IS(J)=-1
      END DO
C
      DO I=1,ND
C ... Analyze individual coordinates : separate for types 1,2,3
C ... phi/psi : 3 states - A=alpha B=C7eq C=beta
C ... J loop over coordinates
        DO J=1,I1
          JPHI=2*J-1
          JPSI=2*J
          IS(J)=-1
          DPHI=ABS(DIHE(I,JPHI) - PHIA)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
          DPSI=ABS(DIHE(I,JPSI) - PSIA)
          IF(DPSI.GT.180.0) DPSI=DPSI-360.0
            DIST=SQRT(DPHI**2 + DPSI**2)
            IF(DIST.LE.DPPA) IS(J)=0
          DPHI=ABS(DIHE(I,JPHI) - PHIB)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
          DPSI=ABS(DIHE(I,JPSI) - PSIB)
          IF(DPSI.GT.180.0) DPSI=DPSI-360.0
            DIST=SQRT(DPHI**2 + DPSI**2)
            IF(DIST.LE.DPPB) IS(J)=1
          DPHI=ABS(DIHE(I,JPHI) - PHIC)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
          DPSI=ABS(DIHE(I,JPSI) - PSIC)
          IF(DPSI.GT.180.0) DPSI=DPSI-360.0
            DIST=SQRT(DPHI**2 + DPSI**2)
            IF(DIST.LE.DPPC) IS(J)=2
        END DO
C
C...single 3-modal, states A,B,C=(g-,t,g+)
        JDIH=IH1
        DO J=I1+1,I2
          JDIH=JDIH+1
          IS(J)=-1
          DPHI=-999.00
          DPHI=ABS(DIHE(I,JDIH) - DI3A)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
            IF(ABS(DPHI).LE.DIS3) IS(J)=0
          DPHI=ABS(DIHE(I,JDIH) - DI3B)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
            IF(ABS(DPHI).LE.DIS3) IS(J)=1
          DPHI=ABS(DIHE(I,JDIH) - DI3C)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
            IF(ABS(DPHI).LE.DIS3) IS(J)=2
C...test
C          WRITE(*,*), ' test3: I,J,JDIH,IS(J),DIHE(I,JDIH)',
C     $       I,J,JDIH,IS(J),DIHE(I,JDIH)
        END DO
C
C...single 2-modal, states A,B
        JDIH=IH2
        DO J=I2+1,I3
          JDIH=JDIH+1
          IS(J)=-1
          DPHI=ABS(DIHE(I,JDIH) - DI2A)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
            IF(ABS(DPHI).LE.DIS2) IS(J)=0
          DPHI=ABS(DIHE(I,JDIH) - DI2B)
          IF(DPHI.GT.180.0) DPHI=DPHI-360.0
            IF(ABS(DPHI).LE.DIS2) IS(J)=1
        END DO
C ...        
C ... analyze global state MS = microstate number
C ... .LT.0 = transition state, GT.0 - one of microstates
      ISTT = 0
      DO J=1,NC
       IF(IS(J).LT.0) ISTT=ISTT-1
      END DO
      IF(ISTT.LT.0) THEN
        MS(I) = -1
      ELSE
        MS(I)=1
        DO J=1,NC
         MS(I) = MS(I) + IS(J)*KPOW(J)
        END DO
      END IF
      MSI=MS(I)
C
C...Count direct residence time and save explored patterns
      IF(MSI.GT.0) THEN
        REST(MSI) = REST(MSI) + 1.0
        ITOTS = ITOTS + 1
        IF(IFOUND(MSI).LE.0) THEN
          DO J=1,NC
            IPATT(MSI,J) = IS(J)
          END DO
          IFOUND(MSI) = 1
        END IF
      END IF
C
C...Analyze transitions MSI = current microstate  
C...                   MSIP = previous microstate
C...                 MSIPIN = previous microstate with > 0 no
C... MSI,MSIP > 0 means we are inside a microstate
C... MSI,MSIP < 0 means we are in transition region
C... Start: with first step - if MSI>0, tag it
      IF(I.EQ.1) THEN
        IF(MSI.GT.0) THEN
          MSIPIN = MSI
          RESTIM(MSI)=RESTIM(MSI)+1.0
        ELSE
          TTRAN=1.0
        END IF
        MSIP = MSI
      ELSE
C...all other steps for I>1
C... branch #1. Both negative - no change, count time
        IF(MSI.LT.0 .AND. MSIP.LT.0) THEN
          TTRAN=TTRAN+1.0
          MSIP=MSI
        END IF
C...branch #2. Remain in same pos state, no change
        IF(MSIP.GT.0 .AND. MSI.GT.0 .AND. MSIP.EQ.MSI) THEN
           RESTIM(MSI)=RESTIM(MSI)+1.0
           MSIPIN=MSI
           MSIP=MSI
        END IF
C...state changes now
C...branch #3. pos to pos direct jump
        IF(MSIP.GT.0 .AND. MSI.GT.0 .AND. MSIP.NE.MSI) THEN
C...test
          WRITE(*,*) '    ----Transition ',MSIP,' -> ', MSI
          NTRAN(MSI,MSIP) = NTRAN(MSI,MSIP) + 1
          RESTIM(MSI)=RESTIM(MSI)+1.0
          MSIPIN=MSI
          MSIP=MSI
        END IF
C...branch #4. pos to neg - start a transition, set MSIPIN
        IF(MSIP.GT.0 .AND. MSI.LT.0) THEN
          TTRAN=1.0
          MSIPIN=MSIP
          MSIP=MSI
        END IF
C...branch #5. neg to pos  - three options
C... 
        IF(MSIP.LT.0 .AND. MSI.GT.0) THEN
          IF(MSIPIN.GT.0) THEN
C.........Option 5A: revisit same pos state w/o tran
            IF(MSI.EQ.MSIPIN) THEN
              RESTIM(MSI)=RESTIM(MSI)+TTRAN+1.0
              MSIP=MSI
            ELSE
C.........Option 5B: transition between different microstates
C...test
              WRITE(*,*) '    ----Transition ',MSIPIN,' -> ', MSI
              NTRAN(MSI,MSIPIN) = NTRAN(MSI,MSIPIN) + 1
              RESTIM(MSI)=RESTIM(MSI)+ (TTRAN/2.0) + 1.0
              RESTIM(MSIPIN)=RESTIM(MSIPIN)+(TTRAN)/2.0
              MSIPIN=MSI
              MSIP=MSI
              TTRAN=0.0
            END IF !  MSI.EQ.MSIPIN
          END IF !  MSIPIN.GT.0
C........Option C: special for first microstate found - no tran
C...     TPREV = the time before first state is found
          IF(MSIPIN.EQ.0 .AND. MSI.GT.0) THEN
            RESTIM(MSI)=RESTIM(MSI)+1.0
            MSIPIN=MSI
            MSIP=MSI
            TPREV=TTRAN
            TTRAN=0.0
          END IF
        END IF  ! MSIP.LT.0 .AND. MSI.GT.0
      END IF  ! I.EQ.1
C
C... testing
      WRITE(*,'(A,I10,10I4)') '+++I,IS(J),ISTT,MS(I)=',
     $      I,(IS(J),J=1,NC),ISTT,MS(I)

C ... end of data loop
      END DO  ! DO I=1,ND
C
C...output  
C...convention: k(i,j) = rate of j->i transition
C...FRAC = part of time spent inside microstate boundaries
C
       DO K=1,NSTATE
         REST(K) = REST(K)*DT
         RESTIM(K) = RESTIM(K)*DT
       END DO
       FRAC=REAL(ITOTS)/REAL(ND)
       WRITE(*,*) ' ...ND,ITOTS,ITOTS/ND=',ND,ITOTS,FRAC
       WRITE(*,*) ' ...Raw microstate residence times [ps], DT=',DT
       DO K=1,NSTATE
       WRITE(*,'(1X,I8,F10.2,10I2)') K,REST(K),(IPATT(K,J),J=1,NC)
       END DO
C
       WRITE(*,*) ' ...Corrected microstate residence times [ps]'
       WRITE(*,*) '    = raw + ttran/2 '
       DO K=1,NSTATE
       WRITE(*,'(1X,I8,F10.2,3I2,4X,3I2)') K,RESTIM(K),
     $       (IPATT(K,J),J=1,NC),(KPATT(K,J),J=1,NC)
       END DO
C
C...This checks if corrected residence times add up to total time
C...TPREV and TTRAN are times before first state and after final state
       TOTTIM=0.0
       DO K=1,NSTATE
         TOTTIM=TOTTIM+RESTIM(K)
       END DO
       DO K=1,NSTATE
         POP(K)=RESTIM(K)/TOTTIM
       END DO
       TPREV = TPREV*DT
       TTRAN = TTRAN*DT
       WRITE(*,'(A,F12.1)') ' ...Testing: SUM of TIMEs =',TOTTIM
       WRITE(*,'(A,2F12.1)') ' ...edges: TPREV,TTRAN=',TPREV,TTRAN
       WRITE(*,*)
C
       WRITE(*,*) ' ...Microstate populations POP=RESTIM/TOTTIM'
       DO K=1,NSTATE
       WRITE(*,'(1X,I8,F10.6,2X,10I2)') K,POP(K),
     $       (KPATT(K,J),J=1,NC)
       END DO
C
C
       WRITE(*,*) ' ...Non-zero transition counts'
       DO K=1,NSTATE
       DO J=1,NSTATE
         IF(NTRAN(J,K).GT.0) THEN
           WRITE(*,'(1X,I4," -> ",I4,I8)') K,J,NTRAN(J,K)
         ENDIF
       END DO
       END DO
C
       WRITE(*,*) ' ...Transition count matrix'
       DO K=1,NSTATE
       WRITE(*,'(1X,I4,4X,18I4)') K,(NTRAN(K,J),J=1,NSTATE)
       END DO
C
       DO K=1,NSTATE
       DO L=K,NSTATE
         S = (NTRAN(K,L) + NTRAN(L,K))/2.0
         STRAN(K,L) = S
         STRAN(L,K) = S
       END DO
       END DO
       WRITE(*,*) ' ...Symmetrized transition count matrix'
       DO K=1,NSTATE
       WRITE(*,'(1X,I4,2X,18F6.1)') K,(STRAN(K,J),J=1,NSTATE)
       END DO
C
       WRITE(*,*) ' ... Rates per ns '
       WRITE(*,*) '     rate(i,j) = 1000*stran(i,j)/restim(j)'
       DO K=1,NSTATE
       DO J=1,NSTATE
         IF(RESTIM(J).GT.0) THEN
           RATE(K,J) = 1000*STRAN(K,J)/RESTIM(J)
         ELSE 
           RATE(K,J) = 0.0
         END IF
       END DO
       END DO
C
       DO K=1,NSTATE
       WRITE(*,'(1X,I4,2X,18F6.1)') K,(RATE(K,J),J=1,NSTATE)
       END DO
C
      OPEN(UNIT=8,FILE="popk.dat",FORM="FORMATTED",STATUS="NEW")
C
C ... Write output file
C
       DO K=1,NSTATE
         WRITE(8,'(E16.8,2X,10I2)') POP(K),(KPATT(K,J),J=1,NC)
       END DO
       DO K=1,NSTATE
         WRITE(8,'(I4,2X,18E16.8)') K,(RATE(K,J),J=1,NSTATE)
       END DO
       WRITE(*,*) '... populations and rates written to popk.dat'
C...start connectivity calculatino
C...Check connectivity using counts
C...ICON(K,L) = 0/1 - connectivity table
C...NICON - tot no of connections for state I
C...KK - iteration counter
C... JCON(KK,JJ) - JJ-th connection at iteration KK
C...  JJC(KK) - total no of connections at iteration KK
C... for each coord perform NSTATE iterations
C... ICUT - cutoff for ignored transitions (=0 - use all)
C... IFOUND(I) - tag found connections to speed up search
C
       ICUT=0
       DO K=1,NSTATE
         DO L=1,NSTATE
           ICON(K,L) = 0
         END DO
       END DO
C      
       DO I=1,NSTATE
           DO L=1,NSTATE
             IFOUND(L)=0
           END DO
C...test
C         WRITE(*,*) '+++Connect for state I=',I
         NICON=0
         JJ=0
         KK=1
C...direct connections in column I of NTRAN
         DO J=1,NSTATE
           IF(NTRAN(J,I).GT.ICUT) THEN 
             ICON(J,I)=1
             NICON=NICON+1
             JJ=JJ+1
             JCON(KK,JJ)=J
             IFOUND(J)=1
           END IF
         END DO
C...test
C             WRITE(*,*) '...direct: NICON,JJ=',NICON,JJ
C             WRITE(*,*) '...JCON(1,*)=',(JCON(1,K),K=1,JJ)
C...iterate through connected columns till end
100      IF(KK.LE.NSTATE) THEN
C...test
C             WRITE(*,*) '......iteration: KK=',KK
             JJC(KK)=JJ
C...test
             JTEST=JJC(KK)
C             WRITE(*,*) '......JJC(KK)=',JTEST
C             WRITE(*,*) '......JCON(KK,I)=',(JCON(KK,K),K=1,JTEST)
           JJ=0
           DO K=1,JTEST
           JJCK=JCON(KK,K)
           DO L=1,NSTATE
              IF(NTRAN(L,JJCK).GT.ICUT .AND. IFOUND(L).LT.1) THEN
                ICON(L,I)=1
                JJ=JJ+1
                JCON(KK+1,JJ)=L
                NICON=NICON+1
                IFOUND(L)=1
C...test
C              WRITE(*,*) '......K,L,JJCK,JJ,NICON=',K,L,JJCK,JJ,NICON
              END IF
           END DO
           END DO  ! K=1,JTEST
           KK=KK+1
           GO TO 100
         END IF  ! KK.LT.NSTATE
       END DO    ! I=1,NSTATE
C
C...Symmetrize
C       DO K=1,NSTATE
C       DO L=K+1,NSTATE
C         IF(ICON(L,K).EQ.1 .OR. ICON(K,L).EQ.1)  THEN
C            ICON(K,L)=1
C            ICON(L,K)=1
C         END IF
C       END DO
C       END DO
C
       WRITE(*,*)  '...Connectivity matrix for ICUT=',ICUT
       DO K=1,NSTATE
         WRITE(*,'(I4,4X,32I2)') K,(ICON(L,K),L=1,NSTATE)
       END DO
C
C ...Re-calculate connectivity with ICUT>0
       ICUT=2
       DO K=1,NSTATE
         DO L=1,NSTATE
           ICON(K,L) = 0
         END DO
       END DO
C      
       DO I=1,NSTATE
           DO L=1,NSTATE
             IFOUND(L)=0
           END DO
         NICON=0
         JJ=0
         KK=1
C...direct connections in column I of NTRAN
         DO J=1,NSTATE
           IF(NTRAN(J,I).GT.ICUT) THEN 
             ICON(J,I)=1
             NICON=NICON+1
             JJ=JJ+1
             JCON(KK,JJ)=J
             IFOUND(J)=1
           END IF
         END DO
C...iterate through connected columns till end
200      IF(KK.LE.NSTATE) THEN
             JJC(KK)=JJ
             JTEST=JJC(KK)
           JJ=0
           DO K=1,JTEST
           JJCK=JCON(KK,K)
           DO L=1,NSTATE
              IF(NTRAN(L,JJCK).GT.ICUT .AND. IFOUND(L).LT.1) THEN
                ICON(L,I)=1
                JJ=JJ+1
                JCON(KK+1,JJ)=L
                NICON=NICON+1
                IFOUND(L)=1
              END IF
           END DO
           END DO  ! K=1,JTEST
           KK=KK+1
           GO TO 200
         END IF  ! KK.LT.NSTATE
       END DO    ! I=1,NSTATE
C
C ...print again
       WRITE(*,*)  '...Connectivity matrix for ICUT=',ICUT
       DO K=1,NSTATE
         WRITE(*,'(I4,4X,32I2)') K,(ICON(L,K),L=1,NSTATE)
       END DO
C...end connectivity calculation
C
C... generate masked population and rates
C
      NSO = 0
      DO K=1,NSTATE
        NSO = NSO + IMASK(K)
      END DO
      WRITE(*,*) '... number of output states after mask=',NSO
      POPT=0.0
      KK=0
      DO K=1,NSTATE
        IF(IMASK(K).GT.0) THEN
          KK=KK+1
          POPO(KK) = POP(K)
          POPT = POPT + POP(K)
          DO J=1,NC
            KPATTO(KK,J)=KPATT(K,J)
          END DO
        END IF
      END DO
C...Normalize POPO
      DO K=1,NSO
        POPO(K) = POPO(K)/POPT
      END DO
       WRITE(*,*) ' ...Microstate populations - masked'
       DO K=1,NSO
       WRITE(*,'(1X,I8,F10.6,2X,10I2)') K,POPO(K),
     $       (KPATTO(K,J),J=1,NC)
       END DO
       WRITE(*,*) '... total masked state pop POPT=',POPT
C
C...eliminate infrequent transitions - below TCUT
       DO K=1,NSTATE
       DO J=1,NSTATE
         IF(RESTIM(J).GT.0 .AND. STRAN(K,J).GT.TCUT) THEN
           RATC(K,J) = 1000*STRAN(K,J)/RESTIM(J)
         ELSE 
           RATC(K,J) = 0.0
         END IF
       END DO
       END DO
C
C...   Fill out masked rate matrix
       KK=0
       DO K=1,NSTATE
         IF(IMASK(K).GT.0) THEN
           KK=KK+1
           LL=0
           DO L=1,NSTATE
             IF(IMASK(L).GT.0) THEN
               LL=LL+1
               RATEO(KK,LL)=RATC(K,L)
             END IF
           END DO
         END IF    ! IMASK(K).GT.0
       END DO      ! DO K=1,NSTATE
       WRITE(*,*) ' ...masked rates per ns '
C
       DO K=1,NSO
       WRITE(*,'(1X,I4,2X,18F6.1)') K,(RATEO(K,J),J=1,NSO)
       END DO
C
C
      OPEN(UNIT=9,FILE="popkm.dat",FORM="FORMATTED",STATUS="NEW")
C
C
C ... Write masked output file
C
       DO K=1,NSO
         WRITE(9,'(E16.8,2X,10I2)') POPO(K),(KPATTO(K,J),J=1,NC)
       END DO
       DO K=1,NSO
         WRITE(9,'(I4,2X,18E16.8)') K,(RATEO(K,J),J=1,NSO)
       END DO
       WRITE(*,*) '... masked pops and rates written to popkm.dat'
      STOP 'What was that?'
      END
