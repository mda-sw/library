C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Carry out a LINEAR DISCRIMINANT ANALYSIS, assigning ungrouped items  C
C  to the closest group centre, using the Mahalanobis distance.         C
C                                                                       C
C                                                                       C
C  To call:   CALL LDA(N,M,DATA,GP,IPRINT,MEAN,MGP,TOTAL,DIFFVC,        C
C                   W1,W2,NOG,IW1,IW2,IERR)     where                  C
C                                                                       C
C                                                                       C
C  N, M  : integer dimensions of ...                                    C
C  DATA  : input data (real).                                           C
C          On output the first column of DATA contains projections of   C
C          all N items on the line connecting the two group means.      C
C          Zero is the boundary point.                                  C
C  GP    : Integer vector of length N giving group assignments.  An     C
C          unassigned item has group 0.  Otherwise groups 1 and 2 will  C
C          be looked for, and other values here are not acceptable.     C
C  IPRINT: integer; print options (= 3: full; otherwise none).          C
C  MEAN  : real vector of length M (number of attributes or variables). C
C  MGP   : real array of dimensions 2 by M.                             C
C  TOTAL : real array of dimensions M by M; on output contains inverse  C
C          of total variance/covariance matrix.                         C
C  DIFFVC: real vector of length M.                                    C
C  W1, W2: real vectors of length M.                                    C
C                                                                       C
C                                                                       C
C  Inputs here are N, M, DATA, GP, IPRINT (and IERR).                   C
C  The principle output information is contained in DATA.               C
C  IERR = 1 means that more than two groups have been specified; IERR   C
C  = 2 means that the total variance-covariance matrx was singular.     C
C                                                                       C
C  Note: we require N > M > 2, to prevent the seeking of the inverse    C
C        of a singular matrix.                                          C
C                                                                       C
C                                                                       C
C  F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.     C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE LDA(N,M,DATA,GP,IPRINT,MEAN,MGP,TOTAL,DIFFVC,W1,
     X                                          W2,NOG,IW1,IW2,IERR)
        REAL    DATA(N,M), TOTAL(M,M), MEAN(M), MGP(2,M), DIFFVC(M)
        REAL    W1(M), W2(M)
        INTEGER GP(N), NOG(2), IW1(M), IW2(M)
C
C          Form global mean.
C
        IERR = 0
        NEFF = 0
        DO 50 I = 1, N
           IF (GP(I).NE.0) NEFF = NEFF + 1
           IF (GP(I).LE.2) GOTO 40
              IERR = 1
              GOTO 9000
   40   CONTINUE
   50   CONTINUE
C
        DO 200 J = 1, M
           MEAN(J) = 0.0
           DO 100 I = 1, N
              IF (GP(I).NE.0) MEAN(J) = MEAN(J) + DATA(I,J)
  100      CONTINUE
           MEAN(J) = MEAN(J)/FLOAT(NEFF)
  200   CONTINUE
C
C          Form (total) variance-covariance matrix.
C
        DO 500 J1 = 1, M
           DO 400 J2 = 1, M
              TOTAL(J1,J2) = 0.0
              DO 300 I = 1, N
                 IF (GP(I).NE.0) TOTAL(J1,J2) = TOTAL(J1,J2) + 
     X               (DATA(I,J1)-MEAN(J1))*(DATA(I,J2)-MEAN(J2))
  300         CONTINUE
              TOTAL(J1,J2) = TOTAL(J1,J2)/FLOAT(NEFF)
  400      CONTINUE
  500   CONTINUE
C
C          Form group means.
C
        DO 700 J = 1, M
           DO 600 K = 1, 2
              MGP(K,J) = 0.0
  600      CONTINUE
  700   CONTINUE
C
        DO 900 I = 1, N
           G = GP(I)
           IF (G.EQ.0) GOTO 900
           NOG(G) = NOG(G) + 1
           DO 800 J = 1, M
              MGP(G,J) = MGP(G,J) + DATA(I,J)
  800      CONTINUE
  900   CONTINUE
C
        DO 1100 K = 1, 2
           DO 1000 J = 1, M
              MGP(K,J) = MGP(K,J)/NOG(K)
 1000      CONTINUE
 1100   CONTINUE
C
        IMAT = 1
        IF (IPRINT.EQ.3) CALL LOUTMT(IMAT,M,TOTAL)
C
C          Invert variance-covariance matrix.
C
        CALL LMTINV(M,TOTAL,D,IW1,IW2)
	write(6,*) ' d=',d
        IF (ABS(D).GT.0.00001) GOTO 1150
           IERR = 2
           GOTO 9000
 1150   CONTINUE
        IF (IPRINT.EQ.3) CALL LOUTMT(IMAT,M,TOTAL)
C
C          Form difference vector of group mean vectors.
C
        DO 1200 J = 1, M
           DIFFVC(J) = MGP(1,J) - MGP(2,J)     
           MEAN(J) = (MGP(1,J) + MGP(2,J))/2.
 1200   CONTINUE
C
C          Determine projections and output them.
C
        CALL LPROJX(N,M,DATA,MEAN,W1,TOTAL,DIFFVC)
        IF (IPRINT.EQ.3) CALL LOUTPX(N,M,DATA)
C
C
C       
 9000   CONTINUE
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output a matrix.                                                     C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE LOUTMT(IMAT,M,ARRAY)
C
C          Output array.
C
        DIMENSION ARRAY(M,M)
C
        IF (IMAT.EQ.1)WRITE (6,900) 
        IF (IMAT.EQ.2)WRITE (6,901)
        DO 100 K1 = 1, M
           WRITE (6,1000) (ARRAY(K1,K2),K2=1,M)
  100   CONTINUE
C
  900   FORMAT(' VARIANCE/COVARIANCE MATRIX FOLLOWS.',/)
  901   FORMAT(' INVERSE VAR/COV MATRIX FOLLOWS.',/)
 1000   FORMAT(5(2X,F18.4))
        RETURN
        END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                      C
C  Invert a symmetric matrix and calculate its determinant.            C
C                                                                      C
C                                                                      C
C  To call:      CALL LMTINV(M,ARRAY,DET,W1,W2)   where                C
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
        SUBROUTINE LMTINV(M,ARRAY,DET,IK,JK)
        REAL    ARRAY(M,M)
	INTEGER IK(M), JK(M)
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Output projections of row points.                                    C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE LOUTPX(N,M,PRJN)
        REAL    PRJN(N,M)
C
C
        NUM = 1
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, N
           WRITE (6,1030) K,(PRJN(K,J),J=1,NUM)
  100   CONTINUE
C
 1000   FORMAT(1H0,'PROJECTIONS OF ROW-POINTS FOLLOW.',/)
 1010   FORMAT(' OBJECT   PROJN')
 1020   FORMAT(' ------  ------')
 1030   FORMAT(I5,2X,F8.4)
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                       C
C  Form projections of row-points on factors.                           C
C                                                                       C
C-----------------------------------------------------------------------C
        SUBROUTINE LPROJX(N,M,DATA,MEAN,VEC,TOTINV,DIFF)
        REAL    DATA(N,M), MEAN(M), VEC(M), TOTINV(M,M), DIFF(M)
C
        NUM = 1
        DO 300 K = 1, N
           DO 50 L = 1, M
              VEC(L) = DATA(K,L)
   50      CONTINUE
           DO 200 I = 1, NUM
              DATA(K,I) = 0.0
              DO 100 J1 = 1, M
                 DO 75 J2 = 1, M
                    DATA(K,I) = DATA(K,I) + (VEC(J1)-MEAN(J1))*
     X                          TOTINV(J1,J2)*DIFF(J2)
   75            CONTINUE
  100         CONTINUE
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END
