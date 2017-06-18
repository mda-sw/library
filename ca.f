C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C  Carry out a CORRESPONDENCE ANALYSIS.                                      C
C                                                                            C
C                                                                            C
C  To call:   CALL CA(N,M,DATA,IPRINT,A1,W1,W2,A2,R,C,IERR)    where         C
C                                                                            C
C                                                                            C
C  N, M  : integer dimensions of ...                                         C
C  DATA  : input data.                                                       C
C          On output, DATA contains in first 7 columns the projections of    C
C          the row-points on the first 7 factors.                            C
C  IPRINT: print options.                                                    C
C          = 1: no printed output - arrays/vectors, only, contain items      C
C              calculated.                                                   C
C          = 2: printed output, in addition, of eigenvalues and eigenvectors.C
C          = 3: full printing of items calculated.                           C
C  A1    : matrix to be diagonalized, of dimensions M * M.                   C
C          On output, A1 contains in the first 7 columns the projections of  C
C          the column-points on the first 7 factors.                         C
C  W1,W2 : real vectors of dimension M (see called routines for use).        C
C          On output, W1 contains the eigenvalues (in increasing order of    C
C          magnitude).                                                       C
C  A2    : real array of dimensions M * M (see called routines for use).     C
C  R, C  : vectors of length N and M (respectively) containing, on output,   C
C          the masses of the row and column points.                          C
C  IERR  : error indicator (normally zero).                                  C
C                                                                            C
C                                                                            C
C  Inputs here are N, M, DATA, METHOD, IPRINT (and IERR).                    C
C  Output information is contained in DATA, A1, W1, R and C.                 C
C  All printed outputs are carried out in easily recognizable subroutines    C
C  called from the first subroutine following.                               C
C                                                                            C
C                                                                            C
C  F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.          C
C                                                                            C 
C  History                                                                   C
C                                                                            C
C  Altered subr. names (prec. by "C") to avoid potential pbs. in linking     C
C  with subr.s of same names in other programs.             F.M., Aug. 1990  C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE CA(N,M,DATA,IPRINT,A,W,FV1,Z,R,C,IERR)
        REAL    DATA(N,M), A(M,M), W(M), FV1(M), Z(M,M), R(N), C(M)
C
C          Form row sums and total
C
        TOT = 0.0
        DO 200 I = 1, N
           R(I) = 0.0
           DO 100 J = 1, M
              TOT = TOT + DATA(I,J)
              R(I) = R(I) + DATA(I,J)
  100      CONTINUE
  200   CONTINUE
C
C          Form column sums and then means
C
        DO 400 J = 1, M
           C(J) = 0.0
           DO 300 I = 1, N
              C(J) = C(J) + DATA(I,J)
  300      CONTINUE
           IF (C(J).GT.0.0) GOTO 350
              WRITE (6,*) ' ALL VALUES IN COL.',J,' ARE 0: PLS. REMOVE'
              IERR = 2
              RETURN
  350      CONTINUE
           C(J) = C(J)/TOT
  400   CONTINUE
C
C          Form row means and make data into frequencies
C
        DO 600 I = 1, N
           IF (R(I).GT.0.0) GOTO 450
              WRITE (6,*) ' ALL VALUES IN ROW',I,' ARE 0: PLS. REMOVE'
              IERR = 2
              RETURN
  450      CONTINUE
           R(I) = R(I)/TOT
           DO 500 J = 1, M
              DATA(I,J) = DATA(I,J)/TOT
  500      CONTINUE
  600   CONTINUE
C
C          Form matrix to be diagonalized
C
        DO 900 J1 = 1, M
           DO 800 J2 = 1, M
              A(J1,J2) = 0.0
              DO 700 I = 1, N
                 A(J1,J2) = A(J1,J2) + DATA(I,J1)*DATA(I,J2)/
     X                                (R(I)*SQRT(C(J1)*C(J2)))
  700         CONTINUE
  800      CONTINUE
  900   CONTINUE
C
C          Carry out eigenreduction.
C
        M2 = M
        CALL CTRED2(M,M2,A,W,FV1,Z)
        CALL CTQL2(M,M2,W,FV1,Z,IERR)
        IF (IERR.NE.0) GOTO 9000
C
C          Output eigenvalues and eigenvectors.
C
        IF (IPRINT.GT.1) CALL COUTEV(M,W)
C
C          Determine factors from the eigenvectors of the symmetric
C          matrix which has been diagonalized.
C
        DO 1100 J1 = 1, M
           DO 1000 J2 = 1, M
              Z(J1,J2) = Z(J1,J2)/SQRT(C(J1))
 1000      CONTINUE
 1100   CONTINUE
C
        IF (IPRINT.GT.1) CALL COUTEC(M,Z)
C
C          Determine projections and contributions and output them.
C
        CALL CPROJX(N,M,DATA,Z,FV1,R)
        IF (IPRINT.EQ.3) CALL COUTPX(N,M,DATA)
        CALL CPROJY(M,W,A,Z,FV1,C)
        IF (IPRINT.EQ.3) CALL COUTPY(M,A)
        IF (IPRINT.EQ.3) CALL COUTCX(N,M,DATA,R)
        IF (IPRINT.EQ.3) CALL COUTCY(M,A,C)
C
 9000   RETURN  
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C 
C Reduce a real, symmetric matrix to a symmetric, tridiagonal                C
C matrix.                                                                    C 
C                                                                            C
C To call:    CALL CTRED2(NM,N,A,D,E,Z)    where                             C
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
C F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.           C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE CTRED2(NM,N,A,D,E,Z)
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
C To call:    CALL CTQL2(NM,N,D,E,Z,IERR)    where                           C
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
C F. Murtagh, ST-ECFESA/ESO, Garching-bei-Muenchen, January 1986.            C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE CTQL2(NM,N,D,E,Z,IERR)
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
 1000   IERR = L
 1001   RETURN
        END  
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C  Output array.                                                             C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTMT(N,M,ARRAY)
C
        DIMENSION ARRAY(N,M)
C
        DO 100 K1 = 1, N
           WRITE (6,1000) (ARRAY(K1,K2),K2=1,M)
  100   CONTINUE
C
 1000   FORMAT(10(2X,F8.4))
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C  Output eigenvalues in order of decreasing value.                          C
C  The first eigenvalue in Correspondence Analysis is trivially              C
C  of value 1, and so is ignored.                                            C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTEV(NVALS,VALS)
        DIMENSION       VALS(NVALS)
C
        TOT = 0.0
        DO 100 K = 1, NVALS-1
           TOT = TOT + VALS(K)
  100   CONTINUE
C
        WRITE (6,1000)
        CUM = 0.0
        K = NVALS 
        WRITE (6,1010)
        WRITE (6,1020)
  200   CONTINUE
        K = K - 1
        CUM = CUM + VALS(K)
        VPC = VALS(K) * 100.0 / TOT
        VCPC = CUM * 100.0 / TOT
        WRITE (6,1030) VALS(K),VPC,VCPC
        IF (K.GT.1) GOTO 200
C
        RETURN
 1000   FORMAT(1H0,'EIGENVALUES FOLLOW.',/)
 1010   FORMAT
     X  (' Eigenvalues       As Percentages    Cumul. Percentages')
 1020   FORMAT
     X  (' -----------       --------------    ------------------')
 1030   FORMAT(F10.4,9X,F10.4,10X,F10.4)
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C  Output FIRST SEVEN eigenvectors associated with eigenvalues               C
C  in descending order.                                                      C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTEC(NDIM,VECS)
        DIMENSION       VECS(NDIM,NDIM)
C
        NUM = MIN0(NDIM,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K1 = 1, NDIM
        WRITE (6,1030) K1,(VECS(K1,NDIM-K2+1),K2=1,NUM)
  100   CONTINUE
C
        RETURN
 1000   FORMAT(1H0,'FACTORS FOLLOW.',/)
 1010   FORMAT('  VBLE.   EV-1    EV-2    EV-3    EV-4    EV-5    EV-6 
     X   EV-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------  
     X------')
 1030   FORMAT(I5,2X,7F8.4)
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C  Output projections of row-points on first (up to) 7 factors.              C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTPX(N,M,PRJN)
        REAL    PRJN(N,M)
C
        NUM = MIN0(M,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, N
           WRITE (6,1030) K,(PRJN(K,J),J=1,NUM)
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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C   Output projections of column points on first 7 factors.                  C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTPY(M,PRJNS)
        REAL    PRJNS(M,M)
C
        NUM = MIN0(M,7)
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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C   Form projections of row-points on first 7 factors.                       C
C   (The first - trivial - factor is ignored.)                               C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE CPROJX(N,M,DATA,EVEC,VEC,R)
        REAL    DATA(N,M), EVEC(M,M), VEC(M), R(N)
C
        NUM = MIN0(M,7)
        DO 300 K = 1, N
           DO 50 L = 1, M
              VEC(L) = DATA(K,L)
   50      CONTINUE
           DO 200 I = 1, NUM
              DATA(K,I) = 0.0
              DO 100 J = 1, M
                 DATA(K,I) = DATA(K,I) + VEC(J) * EVEC(J,M-I)
  100         CONTINUE
              IF (R(K).GT.0.0) DATA(K,I) = DATA(K,I)/R(K)
              IF (R(K).EQ.0.0) DATA(K,I) = 0.0
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C   Form projections of column points on first 7 factors.                    C
C   (The first - trivial - factor is ignored.)                               C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE CPROJY(M,EVALS,A,Z,VEC,C)
        REAL    EVALS(M), A(M,M), Z(M,M), VEC(M), C(M)
C
        NUM = MIN0(M,7)
        DO 300 J1 = 1, M
           DO 50 L = 1, M
              VEC(L) = A(J1,L)
   50      CONTINUE
           DO 200 J2 = 1, NUM
              A(J1,J2) = 0.0
              DO 100 J3 = 1, M
                 A(J1,J2) = A(J1,J2) + VEC(J3)*Z(J3,M-J2)*
     X                                   SQRT(C(J3))
  100         CONTINUE
              IF ((EVALS(M-J2).GT.0.0).AND.(C(J1).GT.0.0))
     X                  A(J1,J2) = A(J1,J2)/SQRT(EVALS(M-J2)*C(J1))
              IF ((EVALS(M-J2).EQ.0.0).OR.(C(J1).EQ.0.0)) A(J1,J2)=0.0
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C   Output contributions of row-points to first 7 factors.                   C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTCX(N,M,PRJN,R)
        REAL    PRJN(N,M), R(N)
C
        NUM = MIN0(M,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, N
           WRITE (6,1030) K,((R(K)*PRJN(K,J)**2),J=1,NUM)
  100   CONTINUE
C
 1000   FORMAT(1H0,'ABS. CONTRIBUTIONS OF OBJECTS TO FACTORS FOLLOW.',/)
 1010   FORMAT(' OBJECT  CNTR-1  CNTR-2  CNTR-3  CNTR-4  CNTR-5  CNTR-6
     X  CNTR-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------
     X  ------')
 1030   FORMAT(I5,2X,7F8.4)
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                            C
C   Output contributions of column-points to first 7 factors.                C
C                                                                            C
C----------------------------------------------------------------------------C
        SUBROUTINE COUTCY(M,PRJNS,C)
        REAL    PRJNS(M,M), C(M)
C
        NUM = MIN0(M,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, M
           WRITE (6,1030) K,(C(K)*PRJNS(K,J)**2,J=1,NUM)
  100   CONTINUE
C
 1000   FORMAT
     X    (1H0,'ABS. CONTRIBUTIONS OF VARIABLES TO FACTORS FOLLOW.',/)
 1010   FORMAT('  VBLE.  CNTR-1  CNTR-2  CNTR-3  CNTR-4  CNTR-5  CNTR-6
     X  CNTR-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------
     X  ------')
 1030   FORMAT(I5,2X,7F8.4)
        RETURN
        END
