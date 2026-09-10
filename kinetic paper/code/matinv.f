      SUBROUTINE MATINV(A,Y,VV,INDX,N,NP)
C
C    Matrix inversion using LU decomposition
C    from "Numerical recipes", W.H.Press, B.P.Flannery,
C    S.A.Teukolsky and W.T.Vetterling, Cambridge University
C    Press, 1986/1992.
C
C Modified to REAL*8
C
C    A - input matrix, destroyed on return
C    Y - output matrix = inverse of A
C    N - order of A,B (actual dimension)
C    NP - first dimension of A in calling routine
C
      INTEGER I,J,N,NP
      INTEGER INDX(NP)
      REAL*8 A(NP,NP),Y(NP,NP),VV(NP)
      REAL*8 D
C
C ... Setup Y as identity matrix
C
      DO 20 I=1,N
      DO 10 J=1,N
      Y(I,J)=0.0
10    CONTINUE
      Y(I,I)=1.0
20    CONTINUE
C
C ... Perform LU decomposition
C
      CALL LUDCM(A,VV,N,NP,INDX,D)
C
C ... Find inverse by columns
C
      DO 30 J=1,N
      CALL LUBKS(A,N,NP,INDX,Y(1,J))
30    CONTINUE
      RETURN
      END

      SUBROUTINE LUDCM(A,VV,N,NP,INDX,D)
C
C Modified to REAL*8
C
      REAL*8 TINY
      PARAMETER (TINY=1.0D-20)
      INTEGER N,NP,INDX(NP)
      REAL*8 A(NP,NP),VV(NP)
      INTEGER I,IMAX,J,K
      REAL*8 AAMAX,D,DUM,SUM
C
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF(AAMAX.EQ.0.) THEN
          WRITE(*,*) ' **** INVERSION: SINGULAR MATRIX ****'
          RETURN
        ENDIF
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END

      SUBROUTINE LUBKS(A,N,NP,INDX,B)
C
C Modified to REAL*8
C
      INTEGER N,NP
      REAL*8 A(NP,NP),B(N)
      INTEGER INDX(N)
      INTEGER I,II,J,LL
      REAL*8 SUM
C
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
C
      RETURN
      END
