      SUBROUTINE HSHLDR(A,DD,VV,N,NM)
C   
C   Householder diagonalization of a symmetric matrix, using
C   code from "Numerical Recipes" by W.H.Press, B.P.Flannery,
C   S.A.Teukolsky and W.T.Vetterling, Cambridge University Press,
C   Cambridge, 1986.
C    
C    A - on input: matrix to be diagonalized
C        on output: eigenvectors in columns, ordered by ascending
C                   eigenvalues
C    DD - eigenvalues
C    VV - work vector
C    N  - current dimension
C    NM - dimension of A in calling routine
C
      REAL*8 A,DD,VV
      DIMENSION A(NM,NM),DD(NM),VV(NM)
C
      CALL TRED2(A,N,NM,DD,VV)
      CALL TQLI(DD,VV,N,NM,A)
      CALL EIGSRT(DD,A,N,NM)
C
      RETURN
      END

      SUBROUTINE EIGSRT(D,V,N,NP)
C
C  Modified to give ascending order of eigenvalues
C
      IMPLICIT REAL*8 (A-H)
      IMPLICIT REAL*8 (O-Z)
      DIMENSION D(NP),V(NP,NP)
C
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END

      SUBROUTINE TQLI(D,E,N,NP,Z)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT REAL*8 (O-Z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
C
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.20000) STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
C     Omit lines from here >>>
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
C     >>> to here when finding only eigenvalues.
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE TRED2(A,N,NP,D,E)
      IMPLICIT REAL*8 (A-H)
      IMPLICIT REAL*8 (O-Z)
      DIMENSION A(NP,NP),D(NP),E(NP)
C
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1  
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
C     Omit following line if finding only eigenvalues
                A(J,I)=A(I,J)/H
C
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
C     Omit following line if finding only eigenvalues.
      D(1)=0.
C
      E(1)=0.
      DO 23 I=1,N
C     Delete lines from here >>>
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
C     >>> to here when finding only eigenvalues.
        D(I)=A(I,I)
C     Also delete lines from here >>>
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
C     >>> to here when finding only eigenvalues.
23    CONTINUE
      RETURN
      END
