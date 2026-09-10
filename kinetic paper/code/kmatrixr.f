      PROGRAM KMATRIXR
C     Process kinetic rate matrix K : only POP & RATE, no patterns
C        I/O: no integers
C        Added sorting of second EV
C     1. Add diagonal elements and symmetrize
C     2. Diagonalize
C     
C     INPUT:
C       1. NS - number of states
C       2. file-name - data file name; file contains
C          POP(I),I=1,NS 
C             - state populations
C           RATE(I,J),I,J=1,NS - rate matrix from trajectory
C     OUTPUT: to peqk.dat
C             PEQ(I) from kinetics
C             K matrix with filled out diagonal
C           : to ev.dat
C             - eigenvalues and eigenvectors
C
C     External routines: hshldr, simsor
C
C     compilation: gfortran kmatrixr.f hshldr.f simsor.f -o kmatrixr
C
C     K. Kuczera  Lawrence, KS.  5-Feb-2018, 7-May-2018
C     MAXS - maximum number of microstates
C
      PARAMETER(MAXS=7000)
      REAL*8 RATE(MAXS,MAXS),RATS(MAXS,MAXS),VECT(MAXS,MAXS)
      REAL*8 POP(MAXS),PEQ(MAXS),EV(MAXS),W(MAXS)
      REAL*8 TOT,SUM,VECR(MAXS,MAXS),VECL(MAXS,MAXS)
      INTEGER IPTA(MAXS)
      INTEGER I,J,K,NS
      CHARACTER*72 FILE1
C
C... read input
      READ(5,*) NS
      WRITE(*,*) ' >>>>>>>>>>>>Kinetic matrix analysis<<<<<<<<<<<<<'
      WRITE(*,*) ' ...number of states NS =',NS
      IF(NS.GT.MAXS) STOP ' Increase MAXS'
C
	READ(5,901) FILE1
	OPEN (UNIT=1,FILE=FILE1,FORM='FORMATTED',STATUS='OLD')
        WRITE(6,*) ' ...Data will be read from ',FILE1
901	FORMAT(A72)
C
C ...Read in data
      DO I=1,NS
        READ(1,'(E16.8)') POP(I)
      END DO
      DO I=1,NS
        READ(1,'(5E16.8)') (RATE(I,J),J=1,NS)
      END DO
C
C ...Test read
      WRITE(*,*) '... Populations '
      TOT=0.0
      DO I=1,NS
        WRITE(*,'(I4,F12.6)') I,POP(I)
        TOT=TOT+POP(I)
      END DO
      WRITE(*,'(A,F12.6)') '...... Testing: TOT =',TOT
      WRITE(*,*) '... Input rates '
      DO I=1,NS
        WRITE(*,'(I4,18F8.2,500(/4X,18F8.2))') I,(RATE(I,J),J=1,NS)
      END DO
C
C... 1. Add diagonal elements
      
      DO I=1,NS
        TOT=0.0D0
        DO J=1,NS
          TOT=TOT+RATE(J,I)
        END DO
        RATE(I,I) = -TOT
      END DO
      WRITE(*,*) '... diagonal elements added to RATE(I,J) '
      DO I=1,NS
        WRITE(*,'(I4,18F8.2,500(/4X,18F8.2))') I,(RATE(I,J),J=1,NS)
      END DO
C
C... 2. Symmetrize
C
      DO I=1,NS
        I1=I+1
        DO J=I1,NS
          TOT=SQRT(RATE(I,J)*RATE(J,I))
          RATS(I,J) = TOT
          RATS(J,I) = TOT
        END DO
        RATS(I,I) = RATE(I,I)
      END DO
      WRITE(*,*) '... Symmetrized rate matrix '
      DO I=1,NS
        WRITE(*,'(I4,18F8.2,500(/4X,18F8.2))') I,(RATS(I,J),J=1,NS)
      END DO
C 
C...  Make a copy - HSHLDR overwrites input
      DO I=1,NS
      DO J=1,NS
        VECT(I,J)=RATS(I,J)
      END DO
      END DO
C
C... 3. Diagonalize symmetrized matrix RATS
C       output VECT=eigenvectors, EV= eigenvalues
C       order = by increasing EV, zero is last
C
      CALL HSHLDR(VECT,EV,W,NS,MAXS)
C
C...
      WRITE(*,*) '... Eigenvalues as rates -EV(I), I=1,NS [ns-1]'
      WRITE(*,'(18E14.6)') (-EV(I),I=1,NS)
      WRITE(*,*) '... Lifetimes -1/EV(I) [ns], I=1,NS-1'
      DO I=1,NS-1
        WRITE(*,'(F14.3)') -1.0D0/EV(I)
      END DO
      WRITE(*,*) '... Eigenvectors of symmetrized RATS in rows '
      DO I=1,NS
        WRITE(*,'(I4,18E14.6,500(/4X,18E14.6))') I,(VECT(J,I),J=1,NS)
      END DO
C
C... Find Left and Right eigenvectors of RATE
C    VECT(I,NS) = SQRT(Peq(I))
C
      DO I=1,NS
        W(I) = VECT(I,NS)
      END DO
      DO I=1,NS
      DO J=1,NS
        VECR(I,J) = VECT(I,J)*W(I)
        VECL(I,J) = VECT(I,J)/W(I)
      END DO
      END DO
      WRITE(*,*) '... L,R & sym EV with slowest rate'
      DO J=1,NS
        WRITE(*,'(I4,2X,3F12.6)') J,VECL(J,NS-1),VECR(J,NS-1),
     $           VECT(J,NS-1)
      END DO
C
C ...Print sorted VEL(J,NS-1) with pointers
C
      DO I=1,NS
        W(I) = VECL(I,NS-1)
      END DO
      CALL SIMSOR(W,IPTA,NS)
      WRITE(*,*) '... Left EV with slowest rate, sorted with pointers'
      DO J=1,NS
        WRITE(*,'(F12.6,I6)') W(J),IPTA(J)
      END DO
C
C...Peq analysis
      WRITE(*,*) '... Peq from zero eigenvector: normalized'
      TOT=0.0
      DO I=1,NS
        TOT=TOT+VECR(I,NS)
      END DO
      WRITE(*,*) '... TOT =',TOT
      WRITE(*,*) '...  VECR(I,NS)/TOT  POP '
      DO I=1,NS
        W(I)=VECR(I,NS)/TOT
      END DO
      DO I=1,NS
        WRITE(*,'(I4,2F12.6)') I,W(I),POP(I)
      END DO
C
C...testing diagonalization
      WRITE(*,*) '... Testing diagonalization, ICNT=0 is good'
      EPS = 1.0D-6
      ICNT=0
      DO I=1,NS
        DO J=1,NS
          TOT=0.0
          DO K=1,NS
          TOT=TOT+RATS(J,K)*VECT(K,I)
          END DO
          IF(ABS(TOT-EV(I)*VECT(J,I)).GT.EPS) ICNT=ICNT+1
        END DO
      END DO
      WRITE(*,*) '... ICNT=',ICNT
C
C...   Output to peqk.dat
       OPEN(UNIT=8,FILE="peqk.dat",FORM="FORMATTED",STATUS="NEW")
       DO I=1,NS
        WRITE(8,'(E16.8)') W(I)
       END DO
       DO I=1,NS
        WRITE(8,'(5E16.8)') (RATE(I,J),J=1,NS)
       END DO
       WRITE(*,*) '...kinetic Peq and K with diag written to peqk.dat'
       
C...   Output to ev.dat
       OPEN(UNIT=9,FILE="ev.dat",FORM="FORMATTED",STATUS="NEW")
       DO K=1,NS
        WRITE(9,910) EV(K),(VECT(J,K),J=1,NS)
       END DO
 910   FORMAT(E16.8,18E16.8,500(/18E16.8))
       WRITE(*,*) '...rows=eigenvalues+eigenvectors written to ev.dat'
       
C
      STOP 'What was that'
      END
