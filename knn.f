	DIMENSION TRAIN(500,11), TEST(500,11), ICL(500), PRCL(500),
     X            ICL00(500), ICL0(500)
	DIMENSION KLIST(15),DK(15),KPOP(15)
	OPEN(UNIT=21,STATUS='OLD',FILE='ngc.dat')
	OPEN(UNIT=22,STATUS='OLD',FILE='ngc2.dat')
C
        N = 500
	N3 = 500
	M  = 11
C
	DO 40 I = 1, N
	   READ(21,120) ICL00(I),(TRAIN(I,J),J=1,M)
   40   CONTINUE
C
        N1 = 0
        N2 = 0
        DO 60 I = 1, N
           IF (ICL00(I).EQ.1) N1 = N1 + 1
           IF (ICL00(I).EQ.2) N2 = N2 + 1
   60   CONTINUE 
        WRITE (6, 680) N1, N2
C
	DO I = 1, N3
	READ(22,120) ICL0(I),(TEST(I,J),J=1,M)
	ENDDO
  120   FORMAT(I3,11f8.5)
C
        do k = 1, 15
           do iopt = 1, 2
C
        IF (IOPT.EQ.1) THEN
C          Case of equal priors
           P1 = 0.5
           P2 = 0.5
           WRITE (6, 660) K
        ELSE IF (IOPT.EQ.2) THEN
           P1 = FLOAT(N1)/FLOAT(N)
           P2 = FLOAT(N2)/FLOAT(N)
           WRITE (6, 670) K
        ENDIF
C
	CALL KNN(TRAIN,N1,N2,N,TEST,N3,M,K,KLIST,DK,KPOP,ICL,PRCL,
     X           ICL00,P1,P2)
C
C       Assess
        IA = 0
        IB = 0
        IC = 0
        ID = 0
        IX = 0
        DO 180 I = 1, N3
           IF (ICL(I).EQ.1.AND.ICL0(I).EQ.1) IA = IA + 1
           IF (ICL(I).EQ.1.AND.ICL0(I).EQ.2) IB = IB + 1
           IF (ICL(I).EQ.2.AND.ICL0(I).EQ.1) IC = IC + 1
           IF (ICL(I).EQ.2.AND.ICL0(I).EQ.2) ID = ID + 1
           IF (ICL(I).EQ.0) IX = IX + 1
  180   CONTINUE
        WRITE (6, 710) 
        WRITE (6, 720)
        WRITE (6, 730)
        WRITE (6, 740) IA, IB
        WRITE (6, 750) IC, ID
        IF (IX.GT.0) WRITE (6, 760) IX
        enddo
        enddo
C
  680   FORMAT(' Training set: cardinalities of 2 groups are',2I5,/,/)
  660   FORMAT(/,/,'       Results of K-NNs.  K = ',I2,
     X                      '.  Priors: equal.',
     X  /,'       ------------------------------------------')
  670   FORMAT(/,/,'       Results of K-NNs.  K = ',I2,
     X                      '.  Priors: proportional.'/,
     X  '       -------------------------------------------------')
  710   FORMAT('                               Teacher',/)
  720   FORMAT('                     label     1     2')
  730   FORMAT('                     -----  ----  ----')
  740   FORMAT(' Classifier            1:   ',I4,2X,I4)
  750   FORMAT('                       2:   ',I4,2X,I4)
  760   FORMAT(1H ,'Number of unclassified observations: ',I4)
C
C
	END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  Carry out a K-NN DISCRIMINANT ANALYSIS,
C  with two groups defined in a training set, and with 
C  assignment of members of a test set. 
C  
C  Parameters:
C  
C  TRAIN(N,M)     training set, where first N1 rows relate to the
C                 first group, and the next N2 rows to the second
C                 group.  Must have N1 + N2 = N. 
C  TEST(N3,M)     test set;    
C  K              number of nearest neighbours to consider;
C  KLIST(K), DK(K), KPOP(K)   are used for storing the K NNs,
C                 their distances to the object under consider-
C                 ation, and the group to which the NNs belong.
C  ICLASS, PRCL   class assignment, and probability of this assgnmt.
C                 Returned values.  Resp. integer and real of dim. N3.
C  ICL            Class assgntms. of training set.
C  P1, P2         Prior probabilities for groups 1 and 2
C 
C  F. Murtagh, ESA/ESO/STECF, Garching.  February 1986.
C  
C  HISTORY
C  Updated to return ICLASS and PRCL,                  F.M., July 1991        
C
C-----------------------------------------------------------------
        SUBROUTINE KNN(TRAIN,N1,N2,N,TEST,N3,M,K,KLIST,DK,KPOP,
     X                 ICLASS,PRCL,ICL,P1,P2)
        DIMENSION      TRAIN(N,M), TEST(N3,M)
        DIMENSION      KLIST(K), DK(K), KPOP(K)
        DIMENSION      ICL(N), ICLASS(N3), PRCL(N3)
C
        DO 90 I   = 1, N3
           DO 10 IX  = 1, K
              KLIST(IX) = 0
              DK(IX) = 1.E+15
              KPOP(IX) = 0
   10      CONTINUE
           CALL DIST(TRAIN,I,TEST,N,N3,M,K,KLIST,DK,KPOP,ICL)
           NUM1 = 0
           NUM2 = 0
           DO 80 IX = 1, K
              IF (KPOP(IX).EQ.1) NUM1 = NUM1 + 1
              IF (KPOP(IX).EQ.2) NUM2 = NUM2 + 1
   80      CONTINUE
C          (Error check:)
           IF ((NUM1+NUM2).EQ.K) GOTO 85
              STOP
   85      CONTINUE
C
C          Posteriors
           DENOM = FLOAT(NUM1)*P1 + FLOAT(NUM2)*P2
           POST1 = FLOAT(NUM1)*P1/DENOM
           POST2 = FLOAT(NUM2)*P2/DENOM
           IF (POST1.GT.POST2) THEN
              ICLASS(I) = 1
              PRCL(I) = POST1
           ENDIF
           IF (POST2.GT.POST1) THEN
              ICLASS(I) = 2
              PRCL(I) = POST2
           ENDIF
           IF (POST1.EQ.POST2) THEN
              ICLASS(I) = 0
              PRCL(I) = POST1
           ENDIF
   90   CONTINUE
C
        RETURN
        END
C-------------------------------------------------------------------------
        SUBROUTINE DIST(TRAIN,J,TEST,N,N3,M,
     X                   K,KLIST,DK,KPOP,ICL00)
        DIMENSION TRAIN(N,M), TEST(N3,M), ICL00(N)
        DIMENSION KLIST(K), DK(K), KPOP(K)
C
        D = 1000000.0
        DO 80 ITRAIN = 1, N
           D = 0.0
           DO 10 ICOLS = 1, M
              D = D + (TRAIN(ITRAIN,ICOLS)-TEST(J,ICOLS))**2
   10      CONTINUE
C
           DO 50 ILOC = 1, K
              IF (D.LT.DK(ILOC)) THEN
C Insert at locn. ILOC and shift right in the 3 length-k lists we're maintaining
                 DO 40 IILOC = K, ILOC+1, -1
C                   Protective measure:
                    IF (IILOC.LE.K) THEN
                       DK(IILOC)    = DK(IILOC-1)
                       KLIST(IILOC) = KLIST(IILOC-1)
                       KPOP(IILOC)  = KPOP(IILOC-1)  
                    ENDIF
C                   Have now freed up space at locn. ILOC
   40            CONTINUE
                 DK(ILOC)    = D
                 KLIST(ILOC) = ITRAIN
                 KPOP(ILOC)  = ICL00(ITRAIN)
                 GOTO 60
              ENDIF
   50      CONTINUE
   60      CONTINUE
   80   CONTINUE
C               
        RETURN
        END
