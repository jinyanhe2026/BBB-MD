      PROGRAM REDDIM2
C     Generate reduced dimesionality rate matrix RHS from full size RATE
C     matrix, using Hummer & Szabo JPCB 2014  Eq. (12)
C       added: Calculate times=1/rates and fluxes=rates*probs
C     PEQ(I),I=1,NS - population of microstates
C     PCG(K),K=1,NCG - pop of aggregated states
C     A - aggegration matrix
C     
C     INPUT:
C       1. NS,NCG,RCUT - NS=number of microstates and NCG=aggregates
C              RCUT - [ns-1] for rates < RCUT: TIME=1/RATE -> -1
C       2. IA(I,J) - NCG lines of NS 0/1 - aggregation mask
C       3. file-name - data file name; file contains
C          PEQ(I),I=1,NS - populations from kinetics
C          I,RATE(I,J),I,J=1,NS - rate matrix from trajectory
C     OUTPUT: to redi.dat
C             PCG(K),K=1,NCG - aggregate pops
C             RHS(K,L),K,L=1,NCG - aggregate rate matrix
C             EV,VECT - eigenvalues and eigenvect of RHS
C
C     compilation: gfortran reddim2.f hshldr.f matinv.f -o reddim2
C
C     K.Kuczera Lawrence,KS. 9-May-2017; updated 04/19/22,01/25/24
C     MAXS - maximum number of microstates
C
      PARAMETER(MAXS=5000)
      REAL*8 RATE(MAXS,MAXS),RHS(MAXS,MAXS),VECT(MAXS,MAXS)
      REAL*8 PEQ(MAXS),PCG(MAXS),EV(MAXS),W(MAXS),A(MAXS,MAXS)
      REAL*8 WORK(MAXS,MAXS),RHSS(MAXS,MAXS),TOT
      REAL*8 VECL(MAXS,MAXS),VECR(MAXS,MAXS)
      INTEGER IA(MAXS,MAXS)
      INTEGER INDX(MAXS)
      INTEGER I,J,K,NC,NCG
      REAL*8  RCUT
      CHARACTER*72 FILE1
C
C... read input
      READ(5,*) NS,NCG,RCUT
      WRITE(*,*) ' >>>>>>>>>Kinetic dimensionality reduction<<<<<<<<<<'
      WRITE(*,*) ' ...# microstates NS =',NS,' aggregates NCG =',NCG
      WRITE(*,'(A,E12.4)') ' ...  rate cutoff RCUT=',RCUT
      IF(NS.GT.MAXS) STOP ' Increase MAXS'
      DO K=1,NCG
        READ(*,*) (IA(I,K),I=1,NS)
      END DO
      DO I=1,NS
      DO K=1,NCG
        A(I,K) = 0.0
        IF(IA(I,K).GT.0) A(I,K)=1.0
      END DO
      END DO
      WRITE(*,*) '... Testing A matrix'
      DO I=1,NS
        WRITE(*,'(I2,2X,18F4.1)') I,(A(I,K),K=1,NCG)
      END DO
C
	READ(5,900) FILE1
	OPEN (UNIT=1,FILE=FILE1,FORM='FORMATTED',STATUS='OLD')
        WRITE(6,*) ' ...Data will be read from ',FILE1
900	FORMAT(A72)
C
C ...Read in data
      DO I=1,NS
        READ(1,*) PEQ(I)
      END DO
      DO I=1,NS
        READ(1,*) (RATE(I,J),J=1,NS)
      END DO
C
C ...Test read
      WRITE(*,*) '... Populations '
      TOT=0.0
      DO I=1,NS
        WRITE(*,'(I4,F12.6,20I2)') I,PEQ(I)
        TOT=TOT+PEQ(I)
      END DO
      WRITE(*,'(A,F12.6)') '...... Testing: TOT =',TOT
      WRITE(*,*) '... Input rates '
      DO I=1,NS
        WRITE(*,901) I,(RATE(I,J),J=1,NS)
      END DO
901   FORMAT(I4,18F8.3,15(/4X,18F8.3))
C
C... 1. Calculate PCG = At*PEQ
      DO K=1,NCG
        TOT=0.0D0
        DO J=1,NS
          TOT=TOT+A(J,K)*PEQ(J)
        END DO
        PCG(K) = TOT
      END DO
      WRITE(*,*) '... Aggregate populations PCG'
      DO K=1,NCG
        WRITE(*,'(I4,F12.6)') K,PCG(K)
      END DO
C
C... 2. WORK = PEQ*1t - RATE  NSxNS
C
      DO I=1,NS
        DO J=1,NS
          WORK(I,J) = PEQ(I) - RATE(I,J)
          VECL(I,J) = WORK(I,J)
        END DO
      END DO
      WRITE(*,*) '... WORK at step 2 '
      DO I=1,NS
        WRITE(*,901) I,(WORK(I,J),J=1,NS)
      END DO
C
C... 3. Invert : MATINV destroys input, output = VECT
C
      CALL MATINV(WORK,VECT,W,INDX,NS,MAXS)
C
      WRITE(*,*) '... Inverse at step 3'
      DO I=1,NS
        WRITE(*,901) I,(VECT(I,J),J=1,NS)
      END DO
C
C... testing inverse 
C
      EPS = 1.0D-6
      ICNT=0
      DO I=1,NS
      DO J=1,NS
        TOT=0
        DO K=1,NS
          TOT=TOT+VECL(I,K)*VECT(K,J)
        END DO
        IF(I.EQ.J) THEN
          IF(ABS(TOT-1.0D0).GT.EPS) ICNT=ICNT+1
        ELSE
          IF(ABS(TOT).GT.EPS) ICNT=ICNT+1
        END IF
      END DO
      END DO
      WRITE(*,*) '... 3: testing inverse: ICNT=0 is good'
      WRITE(*,*) '    ICNT=',ICNT
C
C... 4. WORK = VECT*Dn*A    NSxNCG
C 
      DO I=1,NS
      DO J=1,NCG
        TOT=0.0
        DO K=1,NS
          TOT=TOT+VECT(I,K)*PEQ(K)*A(K,J)
        END DO
        WORK(I,J)=TOT
      END DO
      END DO
C
C... 5. VECT = At*WORK   NCGxNCG
C 
      DO K=1,NCG
      DO L=1,NCG
        TOT=0.0
        DO I=1,NS
          TOT=TOT+A(I,K)*WORK(I,L)
        END DO
        VECT(K,L)=TOT
      END DO
      END DO
C
      WRITE(*,*) '... Matrix at step 5'
      DO I=1,NCG
        WRITE(*,901) I,(VECT(I,J),J=1,NCG)
      END DO
      DO I=1,NCG
      DO J=1,NCG
        VECL(I,J) = VECT(I,J)
      END DO
      END DO
C
C
C... 6. Invert VECT - output in WORK
C
      CALL MATINV(VECT,WORK,W,INDX,NCG,MAXS)
C
C... testing inverse 
C
      EPS = 1.0D-6
      ICNT=0
      DO I=1,NCG
      DO J=1,NCG
        TOT=0
        DO K=1,NCG
          TOT=TOT+VECL(I,K)*WORK(K,J)
        END DO
        IF(I.EQ.J) THEN
          IF(ABS(TOT-1.0D0).GT.EPS) ICNT=ICNT+1
        ELSE
          IF(ABS(TOT).GT.EPS) ICNT=ICNT+1
        END IF
      END DO
      END DO
      WRITE(*,*) '... 6: testing inverse: ICNT=0 is good'
      WRITE(*,*) '    ICNT=',ICNT
C
C
C... 7. RHS = PCG*1t - DN*WORK
C
      DO K=1,NCG
      DO L=1,NCG
        RHS(K,L) = PCG(K) - PCG(K)*WORK(K,L)
      END DO
      END DO
C
      WRITE(*,*) '... Aggregate rates RHS [ns-1]'
      DO K=1,NCG
        WRITE(*,902) K,(RHS(K,L),L=1,NCG)
      END DO
902   FORMAT(I4,10F10.6)
C
      WRITE(*,*) '... Aggregate rates RHS [microseconds-1]'
      FACT=1000.0
      DO K=1,NCG
        WRITE(*,912) K,(FACT*RHS(K,L),L=1,NCG)
      END DO
912   FORMAT(I4,10F10.3)
C
C...  Convert rates to times
      WRITE(*,'(A,F10.6)') '... Times=1/rates [ns] with RCUT=',RCUT
      DO K=1,NCG
        DO L=1,NCG
          EV(L)=0.0
          TOT=ABS(RHS(K,L))
          IF(TOT.GT.RCUT) THEN
            EV(L) = 1.0/TOT
          ELSE 
            EV(L) = 0.0
          ENDIF
        END DO
        WRITE(*,910) K,(EV(L),L=1,NCG)
      END DO
910   FORMAT(I4,18F12.2,15(/4X,18F12.2))
C
C...  Get fluxes - set diagonal to zero
      WRITE(*,*) '... Fluxes F(I,J)=1000*R(I,J)*PCG(J)[us-1]'
      DO K=1,NCG
        DO L=1,NCG
          EV(L)=1000.0*RHS(K,L)*PCG(L)
        END DO
        EV(K) = 0.0
        WRITE(*,911) K,(EV(L),L=1,NCG)
      END DO
911   FORMAT(I4,18F12.6,15(/4X,18F12.6))
C
C
C...symetrize
      DO K=1,NCG
        K1=K+1
        DO L=K1,NCG
          TOT=SQRT(RHS(K,L)*RHS(L,K))
          RHSS(K,L)=TOT
          RHSS(L,K)=TOT
        END DO
        RHSS(K,K) = RHS(K,K)     
      END DO
C
      WRITE(*,*) '... Symmetrized rates RHSS'
      DO K=1,NCG
        WRITE(*,902) K,(RHSS(K,L),L=1,NCG)
      END DO
C
C...  Make a copy - HSHLDR overwrites input
      DO I=1,NCG
      DO J=1,NCG
        VECT(I,J)=RHSS(I,J)
      END DO
      END DO
C... 8. Diagonalize VECT=RHSS
C       output VECT=eigenvectors, EV= eigenvalues
C
      CALL HSHLDR(VECT,EV,W,NCG,MAXS)
C
C...
      WRITE(*,*) '... Eigenvalues as rates -EV(I), I=1,NS [ns-1]'
      WRITE(*,'(18E14.6)') (-EV(I),I=1,NCG)
      WRITE(*,*) '... Lifetimes -1/EV(I) [ns], I=1,NS-1'
      DO I=1,NCG-1
        WRITE(*,'(F14.6)') -1.0D0/EV(I)
      END DO
      WRITE(*,*) '... Eigenvectors of RHSS in rows '
      DO I=1,NCG
        WRITE(*,903) I,(VECT(J,I),J=1,NCG)
      END DO
903   FORMAT(I4,18E14.6,15(/I4,18E14.6))
C
C... Find Left and Right eigenvectors of RHS
C    VECT(I,NS) = SQRT(Peq(I))
C
      DO I=1,NCG
        W(I) = VECT(I,NCG)
      END DO
      DO I=1,NCG
      DO J=1,NCG
        VECR(I,J) = VECT(I,J)*W(I)
        VECL(I,J) = VECT(I,J)/W(I)
      END DO
      END DO
      WRITE(*,*) '... L, R & sym EV with slowest rate'
      DO J=1,NCG
        WRITE(*,'(3F12.6)') VECL(J,NCG-1),VECR(J,NCG-1),
     $           VECT(J,NCG-1)
      END DO
C
C...Peq analysis
      WRITE(*,*) '... Peq in CG from zero eigenvector: normalized'
      TOT=0.0
      DO I=1,NCG
        TOT=TOT+VECR(I,NCG)
      END DO
      WRITE(*,*) '... TOT =',TOT
      WRITE(*,*) '...  VECR(I,NCG)/TOT  PCG '
      DO I=1,NCG
        W(I)=VECR(I,NCG)/TOT
      END DO
      DO I=1,NCG
        WRITE(*,'(I4,2F12.6)') I,W(I),PCG(I)
      END DO
C
C...testing diagonalization
      WRITE(*,*) '... Testing diagonalization, ICNT=0 is good'
      EPS = 1.0D-6
      ICNT=0
      DO I=1,NCG
        DO J=1,NCG
          TOT=0.0
          DO K=1,NCG
          TOT=TOT+RHSS(J,K)*VECT(K,I)
          END DO
          IF(ABS(TOT-EV(I)*VECT(J,I)).GT.EPS) ICNT=ICNT+1
        END DO
      END DO
      WRITE(*,*) '... ICNT=',ICNT
C
C...   Output to redi.dat: PCG, RHS, EV/EVECT
       OPEN(UNIT=8,FILE="redi.dat",FORM="FORMATTED",STATUS="NEW")
       DO I=1,NCG
        WRITE(8,'(E16.8)') PCG(I)
       END DO
       DO I=1,NCG
        WRITE(8,'(I4,2X,18E16.8)') I,(RHS(I,J),J=1,NCG)
       END DO
       DO K=1,NCG
        WRITE(8,'(E14.6,18E16.8)') EV(K),(VECT(J,K),J=1,NCG)
       END DO
       WRITE(*,*) '... PCG, RHS and eigenv written to redi,dat'
C
      STOP 'What was that'
      END
