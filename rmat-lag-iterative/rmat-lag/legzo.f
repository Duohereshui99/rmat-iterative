ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [0,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
C Author: J. M. Jin
C Downloaded from http://jin.ece.illinois.edu/routines/routines.html
          IMPLICIT REAL*8 (A-H,O-Z)
          DIMENSION X(N),W(N)
          data pi/3.1415926535898d0/,one/1/
          N0=(N+1)/2
          DO 45 NR=1,N0
c            Z=COS(pi*(NR-0.25D0)/N)
             Z=COS(pi*(NR-0.25D0)/(N+0.5d0))
  10         Z0=Z
             P=1
             DO 15 I=1,NR-1
  15            P=P*(Z-X(I))
             F0=1
             IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0
             F1=Z
             DO 20 K=2,N
                PF=(2-one/K)*Z*F1-(1-one/K)*F0
                PD=K*(F1-Z*PF)/(1-Z*Z)
                F0=F1
  20            F1=PF
             IF (Z.EQ.0) GO TO 40
             FD=PF/P
             Q=0
             DO 35 I=1,NR-1
                WP=1
                DO 30 J=1,NR-1
                   IF (J.NE.I) WP=WP*(Z-X(J))
  30            CONTINUE
  35            Q=Q+WP
             GD=(PD-Q*FD)/P
             Z=Z-FD/GD
             IF (ABS(Z-Z0).GT.ABS(Z)*1.0D-15) GO TO 10
  40         X(NR)=Z
             X(N+1-NR)=-Z
             W(NR)=2/((1-Z*Z)*PD*PD)
  45         W(N+1-NR)=W(NR)
          x(1:n)=(1+x(n:1:-1))/2
          w(1:n)=w(n:1:-1)/2
          RETURN
          END
