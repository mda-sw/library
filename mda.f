C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Carry out a MULTIPLE DISCRIMINANT ANALYSIS                           C
C              (DISCRIMINANT FACTOR ANALYSIS,                           C
C               CANONICAL DISCRININANT ANALYSIS).                       C
C                                                                       C
C                                                                       C
C  To call:   CALL MDA(N,M,NG,DATA,GP,IPRINT,NOG,MEAN,MGP,TOTAL,        C
C                      BETWEEN,BETW2,CPROJ,W1,W2,IERR)     where        C
C                                                                       C
C                                                                       C
C  N, M  : integer dimensions of ...                                    C
C  DATA  : input data (real).                                           C
C          On output the first NG-1 columns of DATA contain projections C
C          of the N items on the discriminant factors.                  C
C  NG    : (integer) number of groups.                                  C
C  GP    : Integer vector of length N giving group assignments.         C
C          Must be specified correctly - no 0s or values > NG.          C
C  IPRINT: integer; print options (= 3: full; otherwise none).          C
C  NOG   : integer vector of length NG (to contain group cardinalities).C
C  MEAN  : real vector of length M (number of attributes or variables). C
C  MGP   : real array of dimensions 2 by M.                             C
C  TOTAL : real array of dimensions M by M; on output contains inverse  C
C          of total variance/covariance matrix.                         C
C  BETWEEN: real array of dimensions NG by NG.                          C
C  BETW2 : real array of dimensions NG by NG.                           C
C  CPROJ : real array of dimensions M by NG; on output contains the     C
C          coefficients of the discriminant factors in terms of the     C
C          original variables.                                          C
C  W1, W2: real vectors of length M.                                    C
C  IERR  : initially 0; = 1 if there is no convergence in the TQL2      C
C          eigenroutine; = 2 if the total variance-covariance to be     C
C          inverted is singular, - in this case, check that there are   C
C          no columns with identical values, that N > M > NG, etc.      C
C                                                                       C
C                                                                       C
C  Inputs here are N, M, NG, DATA, GP, IPRINT (and IERR).               C
C  The principle output information is contained in DATA and CPROJ; and C
C  W1(NG-1), W1(NG-2) ... contain the eigenvalues in decreasing order.  C
C                                                                       C
C  Notes: we require that N > M > NG (otherwise, an error is likely in  C
C         the matrix inversion routine due to singularity).             C
C         NG-1 eigenvalues eigenvectors are output.                     C
C                                                                       C
C                                                                       C
C  F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.     C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE MDA(N,M,NG,DATA,GP,IPRINT,NOG,MEAN,MGP,TOTAL,
     X                                  BETWEEN,BETW2,CPROJ,W1,W2,IERR)
        REAL    DATA(N,M),TOTAL(M,M),MEAN(M),MGP(NG,M),BETWEEN(NG,NG) 
        REAL    W1(M), W2(M), BETW2(NG,NG), CPROJ(M,NG)
        INTEGER GP(N), NOG(NG)
C
C          Form global mean.
C
        DO 200 J = 1, M
           MEAN(J) = 0.0
           DO 100 I = 1, N
	      DATA(I,J)=DATA(I,J)*15.0     ! *********** BEMERKEN SIE!
              MEAN(J) = MEAN(J) + DATA(I,J)
  100      CONTINUE
           MEAN(J) = MEAN(J)/FLOAT(N)
  200   CONTINUE
C
C          Form (total) variance-covariance matrix.
C
        DO 500 J1 = 1, M
           DO 400 J2 = 1, M
              TOTAL(J1,J2) = 0.0
              DO 300 I = 1, N
                 TOTAL(J1,J2) = TOTAL(J1,J2) + 
     X           (DATA(I,J1)-MEAN(J1))*(DATA(I,J2)-MEAN(J2))
  300         CONTINUE
              TOTAL(J1,J2) = TOTAL(J1,J2)/FLOAT(N)
  400      CONTINUE
  500   CONTINUE
C
        IMAT = 1
C       CALL OUTMAT(IMAT,M,TOTAL)
C
C          Form group means.
C
        DO 700 J = 1, M
           DO 600 K = 1, NG
              MGP(K,J) = 0.0
  600      CONTINUE
  700   CONTINUE
C
        DO 900 I = 1, N
           G = GP(I)
           IF (G.EQ.0) GOTO 9000
           NOG(G) = NOG(G) + 1
           DO 800 J = 1, M
              MGP(G,J) = MGP(G,J) + DATA(I,J)
  800      CONTINUE
  900   CONTINUE
C
        DO 1100 K = 1, NG
           DO 1000 J = 1, M
              MGP(K,J) = MGP(K,J)/NOG(K)
 1000      CONTINUE
 1100   CONTINUE
C
C          Invert variance-covariance matrix.
C
        CALL MATINV(M,TOTAL,D,W1,W2)
        IF (D.GT.0.000000001) GOTO 1150
           IERR = 2
           GOTO 9000
 1150   CONTINUE
        IMAT = 2
C       CALL OUTMAT(IMAT,M,TOTAL)
C
C          Form the symmetric variant of the between-groups 
C          variance-covariance matrix for diagonalization.
C
        DO 1200 K1 = 1, NG
         DO 1200 K2 = 1, NG
         BETWEEN(K1,K2) = 0.0
          DO 1200 J1 = 1, M
           DO 1200 J2 = 1, M
            D1 = MGP(K1,J1) - MEAN(J1)
            D2 = FLOAT(NOG(K1))/FLOAT(N)
            D2 = SQRT(D2)
            D3 = MGP(K2,J2) - MEAN(J2)
            D4 = FLOAT(NOG(K2))/FLOAT(N)
            D4 = SQRT(D4)
            BETWEEN(K1,K2) = BETWEEN(K1,K2) +
     X                       (D1*D2)*TOTAL(J1,J2)*(D3*D4)
 1200   CONTINUE
C
        IMAT = 4
C       CALL OUTMAT(IMAT,M,TOTAL)
C
C          Carry out eigenreduction.
C
        NG2 = NG
        CALL TRED2(NG,NG2,BETWEEN,W1,W2,BETW2)
        CALL TQL2(NG,NG2,W1,W2,BETW2,IERR)
        IF (IERR.NE.0) GOTO 9000
C
C          Output eigenvalues and eigenvectors.
C
        IF (IPRINT.GT.1) CALL OUTEVL(N,M,NG,W1)
        IF (IPRINT.GT.1) CALL OUTEVC(N,M,NG,BETW2,NG-1)
C
C          Convert eigenvectors in NG-space to those in M-space.
C
        DO 1300 J = 1, M
           DO 1300 K = 1, NG
              CPROJ(J,K) = 0.0
              DO 1300 J2 = 1, M
                 DO 1300 K2 = 1, NG
                    D1 = MGP(K2,J2) - MEAN(J2)
                    D2 = FLOAT(NOG(K2))/FLOAT(N)
                    D1 = D1*SQRT(D2)
                    CPROJ(J,K)=CPROJ(J,K)+
     X              TOTAL(J,J2)*D1*BETW2(K2,NG-K+1)
 1300   CONTINUE
        IF (IPRINT.GT.1) CALL OUTEVC(N,NG,M,CPROJ,NG-1)
C
C          Determine projections and output them.
C
        CALL PROJX(N,M,NG,DATA,MEAN,CPROJ,W2,TOTAL)
        IF (IPRINT.EQ.3) CALL OUTPRX(N,M,NG,DATA)
C
C
C       
 9000   CONTINUE
        RETURN
        END

        SUBROUTINE OUTMAT(IMAT,M,ARRAY)
        DIMENSION ARRAY(M,M)
C
        WRITE (6,900) IMAT
        DO 100 K1 = 1, M
           WRITE (6,1000) (ARRAY(K1,K2),K2=1,K1)
  100   CONTINUE
C
  900   FORMAT(' IMAT =',I6)
 1000   FORMAT(10(2X,F8.4))
        RETURN
        END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                      C
C  Invert a symmetric matrix and calculate its determinant.            C
C                                                                      C
C                                                                      C
C  To call:      CALL MATINV(M,ARRAY,DET,W1,W2)   where                C
C                                                                      C
C                                                                      C
C  M       : dimension of ...                                          C
C  ARRAY   : input matrix which is replaced by its inverse.            C
C  NORDER  : degree of matrix (order of determinant)                   C
C  DET     : determinant of input matrix.                              C
C  W1, W2  : work vectors of dimension M.                              C
C                                                                      C
C                                                                      C
C  Reference: Philip B Bevington, "Data Reduction and Error Analysis   C
C             for the Physical Sciences", McGraw-Hill, New York, 1969, C
C             pp. 300-303.                                             C
C                                                                      C
C                                                                      C
C  F. Murtagh, ST-ECF, Garching-bei-Muenchen, January 1986.            C
C                                                                      C
C----------------------------------------------------------------------C
        SUBROUTINE MATINV(M,ARRAY,DET,IK,JK)
        REAL    ARRAY(M,M), IK(M), JK(M)
C
   10   DET = 1.0
   11   DO 100 K = 1, M
C       Find largest element ARRAY(I,J) in rest of matrix.
        AMAX = 0.0
   21      DO 30 I = K, M
              DO 30 J = K, M
   23            IF (ABS(AMAX)-ABS(ARRAY(I,J))) 24,24,30
   24            AMAX = ARRAY(I,J)
                 IK(K) = I
                 JK(K) = J
   30      CONTINUE
C          Interchange rows and columns to put AMAX in ARRAY(K,K).
   31      IF (AMAX) 41,32,41
   32      DET = 0.0
           GOTO 140
   41      I = IK(K)
           IF (I-K) 21,51,43
   43      DO 50 J = 1, M
              SAVE = ARRAY(K,J)
              ARRAY(K,J) = ARRAY(I,J)
   50      ARRAY(I,J) = -SAVE
   51      J = JK(K)
           IF (J-K) 21,61,53
   53      DO 60 I = 1, M
              SAVE = ARRAY(I,K)
              ARRAY(I,K) = ARRAY(I,J)
   60      ARRAY(I,J) = -SAVE
C          Accumulate elements of inverse matrix.
   61      DO 70 I = 1, M
              IF (I-K) 63,70,63
   63         ARRAY(I,K) = -ARRAY(I,K)/AMAX
   70      CONTINUE
   71      DO 80 I = 1, M
              DO 80 J = 1, M
                 IF (I-K) 74,80,74
   74            IF (J-K) 75,80,75
   75            ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
   80      CONTINUE
   81      DO 90 J = 1, M
              IF (J-K) 83,90,83
   83         ARRAY(K,J) = ARRAY(K,J)/AMAX
   90      CONTINUE
           ARRAY(K,K) = 1.0/AMAX        
  100   DET = DET * AMAX
C       Restore ordering of matrix.
  101   DO 130 L = 1, M
           K = M - L + 1
           J = IK(K)
           IF (J-K) 111,111,105
  105      DO 110 I = 1, M
              SAVE = ARRAY(I,K)
              ARRAY(I,K) = -ARRAY(I,J)
  110      ARRAY(I,J) = SAVE
  111      I = JK(K)
           IF (I-K) 130,130,113
  113      DO 120 J = 1, M
              SAVE = ARRAY(K,J)
              ARRAY(K,J) = -ARRAY(I,J)
  120      ARRAY(I,J) = SAVE
  130   CONTINUE
  140   RETURN
        END           
        
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C 
C Reduce a real, symmetric matrix to a symmetric, tridiagonal                C
C matrix.                                                                    C 
C                                                                            C
C To call:    CALL TRED2(NM,N,A,D,E,Z)    where                              C
C                                                                            C
C NM = row dimension of A and Z;                                             C
C N = order of matrix A (will always be <= NM);                              C
C A = symmetric matrix of order N to be reduced to tridiagonal form;         C
C D = vector of dim. N containing, on output, diagonal elts. of trid. matrix;C
C E = working vector of dim. at least N-1 to contain subdiagonal elts.;      C
C Z = matrix of dims. NM by N containing, on output, orthogonal              C
C                    transformation matrix producting the reduction.         C
C                                                                            C
C Normally a call to TQL2 will follow the call to TRED2 in order to          C
C produce all eigenvectors and eigenvalues of matrix A.                      C
C                                                                            C
C Algorithm used: Martin et al., Num. Math. 11, 181-195, 1968.               C
C                                                                            C
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK             C
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 1976,         C
C pp. 489-494.                                                               C
C                                                                            C
C F. Murtagh, ST-ECF Garching-bei-Muenchen, January 1986.                    C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
        REAL A(NM,N),D(N),E(N),Z(NM,N)
C
        DO 100 I = 1, N
           DO 100 J = 1, I
              Z(I,J) = A(I,J)
  100   CONTINUE
        IF (N.EQ.1) GOTO 320
        DO 300 II = 2, N
           I = N + 2 - II
           L = I - 1
           H = 0.0
           SCALE = 0.0
           IF (L.LT.2) GOTO 130
           DO 120 K = 1, L
              SCALE = SCALE + ABS(Z(I,K))
  120      CONTINUE
           IF (SCALE.NE.0.0) GOTO 140
  130      E(I) = Z(I,L)
           GOTO 290
  140      DO 150 K = 1, L
              Z(I,K) = Z(I,K)/SCALE
              H = H + Z(I,K)*Z(I,K)
  150      CONTINUE
C
           F = Z(I,L)
           G = -SIGN(SQRT(H),F)
           E(I) = SCALE * G
           H = H - F * G
           Z(I,L) = F - G
           F = 0.0
C
           DO 240 J = 1, L
              Z(J,I) = Z(I,J)/H
              G = 0.0
C             Form element of A*U.
              DO 180 K = 1, J
                 G = G + Z(J,K)*Z(I,K)
  180         CONTINUE
              JP1 = J + 1
              IF (L.LT.JP1) GOTO 220
              DO 200 K = JP1, L
                 G = G + Z(K,J)*Z(I,K)
  200         CONTINUE
C             Form element of P where P = I - U U' / H .
  220         E(J) = G/H
              F = F + E(J) * Z(I,J)
  240      CONTINUE
           HH = F/(H + H)
C          Form reduced A.
           DO 260 J = 1, L
              F = Z(I,J)
              G = E(J) - HH * F
              E(J) = G
              DO 250 K = 1, J
                 Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  250         CONTINUE
  260      CONTINUE
  290      D(I) = H
  300   CONTINUE
  320   D(1) = 0.0
        E(1) = 0.0
C       Accumulation of transformation matrices.
        DO 500 I = 1, N
           L = I - 1
           IF (D(I).EQ.0.0) GOTO 380
           DO 360 J = 1, L
              G = 0.0
              DO 340 K = 1, L
                 G = G + Z(I,K) * Z(K,J)
  340         CONTINUE
              DO 350 K = 1, L
                 Z(K,J) = Z(K,J) - G * Z(K,I)
  350         CONTINUE
  360      CONTINUE
  380      D(I) = Z(I,I)
           Z(I,I) = 1.0
           IF (L.LT.1) GOTO 500
           DO 400 J = 1, L
              Z(I,J) = 0.0
              Z(J,I) = 0.0
  400      CONTINUE
  500   CONTINUE
C
        RETURN
        END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C Determine eigenvalues and eigenvectors of a symmetric,                     C
C tridiagonal matrix.                                                        C
C                                                                            C
C To call:    CALL TQL2(NM,N,D,E,Z,IERR)    where                            C
C                                                                            C
C NM = row dimension of Z;                                                   C
C N = order of matrix Z;                                                     C
C D = vector of dim. N containing, on output, eigenvalues;                   C
C E = working vector of dim. at least N-1;                                   C
C Z = matrix of dims. NM by N containing, on output, eigenvectors;           C
C IERR = error, normally 0, but 1 if no convergence.                         C
C                                                                            C
C Normally the call to TQL2 will be preceded by a call to TRED2 in           C
C order to set up the tridiagonal matrix.                                    C
C                                                                            C
C Algorithm used: QL method of Bowdler et al., Num. Math. 11,                C
C 293-306, 1968.                                                             C
C                                                                            C
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK             C
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 1976,         C
C pp. 468-474.                                                               C
C                                                                            C
C F. Murtagh, ST-ECF Garching-bei-Muenchen, January 1986.                    C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
        REAL    D(N), E(N), Z(NM,N)
        DATA    EPS/1.E-12/
C
        IERR = 0
        IF (N.EQ.1) GOTO 1001
        DO 100 I = 2, N
           E(I-1) = E(I)
  100   CONTINUE
        F = 0.0
        B = 0.0
        E(N) = 0.0
C
        DO 240 L = 1, N
           J = 0
           H = EPS * (ABS(D(L)) + ABS(E(L)))
           IF (B.LT.H) B = H
C          Look for small sub-diagonal element.
           DO 110 M = L, N
              IF (ABS(E(M)).LE.B) GOTO 120
C             E(N) is always 0, so there is no exit through the bottom 
C             of the loop.
  110      CONTINUE
  120      IF (M.EQ.L) GOTO 220
  130      IF (J.EQ.30) GOTO 1000
           J = J + 1
C          Form shift.
           L1 = L + 1
           G = D(L)
           P = (D(L1)-G)/(2.0*E(L))
           R = SQRT(P*P+1.0)
           D(L) = E(L)/(P+SIGN(R,P))
           H = G-D(L)
C
           DO 140 I = L1, N
              D(I) = D(I) - H
  140      CONTINUE
C
           F = F + H
C          QL transformation.
           P = D(M)
           C = 1.0
           S = 0.0
           MML = M - L
C
           DO 200 II = 1, MML
              I = M - II
              G = C * E(I)
              H = C * P
              IF (ABS(P).LT.ABS(E(I))) GOTO 150
              C = E(I)/P
              R = SQRT(C*C+1.0)
              E(I+1) = S * P * R
              S = C/R
              C = 1.0/R
              GOTO 160
  150         C = P/E(I)
              R = SQRT(C*C+1.0)
              E(I+1) = S * E(I) * R
              S = 1.0/R
              C = C * S
  160         P = C * D(I) - S * G
              D(I+1) = H + S * (C * G + S * D(I))
C             Form vector.
              DO 180 K = 1, N
                 H = Z(K,I+1)
                 Z(K,I+1) = S * Z(K,I) + C * H
                 Z(K,I) = C * Z(K,I) - S * H
  180         CONTINUE
  200      CONTINUE
           E(L) = S * P
           D(L) = C * P
           IF (ABS(E(L)).GT.B) GOTO 130
  220      D(L) = D(L) + F
  240   CONTINUE
C
C       Order eigenvectors and eigenvalues.
        DO 300 II = 2, N
           I = II - 1
           K = I
           P = D(I)
           DO 260 J = II, N
              IF (D(J).GE.P) GOTO 260
              K = J
              P = D(J)
  260      CONTINUE
           IF (K.EQ.I) GOTO 300
           D(K) = D(I)
           D(I) = P
           DO 280 J = 1, N
              P = Z(J,I)
              Z(J,I) = Z(J,K)
              Z(J,K) = P
  280      CONTINUE
  300   CONTINUE
C
        GOTO 1001
C       Set error - no convergence to an eigenvalue after 30 iterations.
 1000   IERR = 1
 1001   RETURN
        END  

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output eigenvalues in order of decreasing value.                     C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE OUTEVL(N,M,NG,VALS)
        DIMENSION       VALS(NG)
C
        TOT = 0.0
        DO 100 K = 2, NG
           TOT = TOT + VALS(K)
  100   CONTINUE
C
        WRITE (6,1000)
        CUM = 0.0
        K = NG + 1 
        WRITE (6,1010)
        WRITE (6,1020)
  200   CONTINUE
        K = K - 1
        CUM = CUM + VALS(K)
        VPC = VALS(K) * 100.0 / TOT
        VCPC = CUM * 100.0 / TOT
        WRITE (6,1030) VALS(K),VPC,VCPC
        IF (K.GT.2) GOTO 200
C
        RETURN
 1000 FORMAT(1H0,'EIGENVALUES FOLLOW.',/)
 1010 FORMAT(' Eigenvalues       As Percentages    Cumul. Percentages')
 1020 FORMAT(' -----------       --------------    ------------------')
 1030 FORMAT(F10.4,9X,F10.4,10X,F10.4)
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output FIRST SEVEN eigenvectors associated with eigenvalues in       C
C  decreasing order.                                                    C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE OUTEVC(N1,N2,N3,VECS,N4)
        DIMENSION       VECS(N3,N3)
C
        NUM = MIN0(N4,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K1 = 1, N3
        WRITE (6,1030) K1,(VECS(K1,K2),K2=1,NUM)
  100   CONTINUE
C
        RETURN
 1000   FORMAT(1H0,'EIGENVECTORS FOLLOW.',/)
 1010   FORMAT('  VBLE.   EV-1    EV-2    EV-3    EV-4    EV-5    EV-6 
     X   EV-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------  
     X------')
 1030   FORMAT(I5,2X,7F8.4)
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output projections on discriminant factors.                          C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE OUTPRX(N,M,NG,PRJN)
        REAL    PRJN(N,M)
C
        NUM = MIN0(N,M,NG,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, N
           WRITE (6,1030) K,(PRJN(K,J),J=1,NUM-1)
  100   CONTINUE
C
 1000   FORMAT(1H0,'PROJECTIONS OF ROW-POINTS FOLLOW.',/)
 1010   FORMAT(' OBJECT  PROJ-1  PROJ-2  PROJ-3  PROJ-4  PROJ-5  PROJ-6
     X  PROJ-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------
     X  ------')
 1030   FORMAT(I5,2X,7F8.4)
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output projections of column points on up to first 7 discriminant    C
C  axes.                                                                C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE OUTPRY(N,M,NG,PRJNS)
        REAL    PRJNS(NG,NG)
C
        NUM = MIN0(N,M,MG,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, M
           WRITE (6,1030) K,(PRJNS(K,J),J=1,NUM)
  100   CONTINUE
C
 1000   FORMAT(1H0,'PROJECTIONS OF COLUMN-POINTS FOLLOW.',/)
 1010   FORMAT('  VBLE.  PROJ-1  PROJ-2  PROJ-3  PROJ-4  PROJ-5  PROJ-6
     X  PROJ-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------
     X  ------')
 1030   FORMAT(I5,2X,7F8.4)
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Form projections of row-points on (up to) first 7 factors.           C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE PROJX(N,M,NG,DATA,MEAN,EVEC,VEC,TOTINV)
        REAL    DATA(N,M), EVEC(M,M), VEC(M), TOTINV(M,M), MEAN(M)
C
        NUM = MIN0(N,M,NG,7)
        DO 300 K = 1, N
           DO 50 L = 1, M
              VEC(L) = DATA(K,L)
   50      CONTINUE
           DO 200 I = 1, NUM
              DATA(K,I) = 0.0
              DO 100 J1 = 1, M
C                DO 75 J2 = 1, M
                    DATA(K,I) = DATA(K,I) + (VEC(J1) - MEAN(J1))*
     X                          EVEC(J1,I)
   75            CONTINUE
  100         CONTINUE
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Determine projections of column points on (up to) 7 factors.         C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE PROJY(N,M,NG,EVALS,A,Z,VEC)
        REAL    EVALS(M), A(M,M), Z(M,M), VEC(M)
C
        NUM = MIN0(N,M,NG,7)
        DO 300 J1 = 1, M
           DO 50 L = 1, M
              VEC(L) = A(J1,L)
   50      CONTINUE
           DO 200 J2 = 1, NUM
              A(J1,J2) = 0.0
              DO 100 J3 = 1, M
                 A(J1,J2) = A(J1,J2) + VEC(J3)*Z(J3,M-J2+1)
  100         CONTINUE
              IF (EVALS(M-J2+1).GT.0.0) A(J1,J2) = 
     X                                  A(J1,J2)/SQRT(EVALS(M-J2+1))
              IF (EVALS(M-J2+1).EQ.0.0) A(J1,J2) = 0.0
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END
 
