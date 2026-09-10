      SUBROUTINE SIMSOR(A,IPTA,N)
C
C ... Simple sorting routine with pointers
C ... On output Elements of A are sorted in ascending order
C ... with IPTA pointing to original elements
C
      REAL*8 A(N),XM,T
      INTEGER IPTA(N)
      INTEGER N,JM,JT
C
      DO I=1,N
        IPTA(I)=I
      END DO
C
      DO I=1,N-1
        XM = A(I)
        JM = I
C
        DO J=I+1,N 
          IF(A(J).LE.XM) THEN
            JM = J
            XM = A(J)
          END IF
        END DO
C
        IF(JM.NE.I) THEN
          T = A(I)
          A(I) = XM
          A(JM) = T
          JT = IPTA(I)
          IPTA(I) = IPTA(JM)
          IPTA(JM) = JT
        END IF
C
      END DO  ! I=1,N-1
C
      RETURN
      END
