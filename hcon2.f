C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                C
C  HIERARCHICAL CLUSTERING using Minimum Variance Criterion,     C
C  using the O(N**2) time Nearest Neighbour Chain algorithm.     C
C                                                                C
C  Parameters:                                                   C
C                                                                C
C  DATA(N,M) :      input data,                                  C
C  IA(N), IB(N), CRIT(N) : sequence of agglomerands and          C
C                   values returned (only locations 1 to N-1     C
C                   are of interest),                            C
C  MEMBR(N), DISS(N), ICHAIN(N) : used in the routines to store  C
C                   cluster cardinalities, nearest neighbour     C    
C                   dissimilarities, and the NN-chain.           C
C  FLAG(N) :        (boolean) used to indicate agglomerable      C
C                   objects and clusters.                        C
C                                                                C
C  Reference: Murtagh, Multidimensional Clustering Algorithms,   C
C             Physica-Verlag, 1985.                              C
C                                                                C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.           C
C                                                                C
C HISTORY                                                        C
C                                                                C
C Bounds bug fix, Oct. 1990, F. Murtagh.                         C
C On line 98, alter stmt. "IDUM=ICHAIN(LEN+1)" to                C
C                         "IF (LEN.LT.N) IDUM=ICHAIN(LEN+1)"     C
C The flag in the following stmt. took care of this eventuality  C
C and the VAX/VMS compiler, without the /CHECK=(BOUNDS) option   C
C did not care about the fact that ICHAIN went out of bounds.    C 
C----------------------------------------------------------------C
      SUBROUTINE HCON2(N,M,DATA,IA,IB,CRIT,MEMBR,DISS,ICHAIN,FLAG)
      REAL    MEMBR(N), DATA(N,M), DISS(N), CRIT(N)
      INTEGER ICHAIN(N), IA(N), IB(N)
      REAL INF
      LOGICAL FLAG(N)
      DATA INF/1.E+25/
C     EQUIVALENCE (ICHAIN(1),IA(1)),(DISS(1),CRIT(1))
C
      DO 150 I=1,N
         MEMBR(I)=1
         FLAG(I)=.TRUE.
  150 CONTINUE
      NCL=N
      I1=1
C
C  Start the NN-chain:
C
  200 LEN=N
      ICHAIN(LEN)=I1
      DISS(LEN)=INF
C
C  Determine NN of object I1
C
  300 FLAG(I1)=.FALSE.
C
C  Turn off FLAG so that 0 diss. of I1 with self not obtained.
C
      D=DISS(LEN)
      IF (LEN.LT.N) I2=ICHAIN(LEN+1)
C
C  For identical diss.'s, above ensures that RNN will be found.
C
      CALL DETNN(DATA,FLAG,MEMBR,N,M,I1,I2,D)
      FLAG(I1)=.TRUE.
C
C  If LEN = 1 place obj. I2 as second obj. in NN-chain.
C
      IF (LEN.LT.N) GOTO 350
      LEN=LEN-1
      IF (LEN.LT.N-NCL) GOTO 700
      ICHAIN(LEN)=I2
      DISS(LEN)=D
      GOTO 500
C
C  If LEN < N distinguish between having RNN & continuing NN-chain.
C
  350 CONTINUE
      IF (I2.NE.ICHAIN(LEN+1)) GOTO 400
C
C  Have RNN.
C
      NCL=NCL-1
      CALL AGGLOM(I1,I2,D,DATA,MEMBR,FLAG,IA,IB,CRIT,NCL,N,M)
      LEN=LEN+2
      GOTO 500
  400 CONTINUE
C
C  Grow extra link on NN-chain.
C
      IDUM=ICHAIN(LEN+1)
      FLAG(IDUM)=.FALSE.
      LEN=LEN-1
      IF (LEN.LE.N-NCL) GOTO 700
      ICHAIN(LEN)=I2
      DISS(LEN)=D
      GOTO 500
C
C  Select obj. for continuing to grow (or restarting) NN-chain.
C
  500 CONTINUE
      IF (NCL.EQ.1) GOTO 600
      IF (LEN.EQ.N+1) GOTO 550
      I1=ICHAIN(LEN)
      FLAG(I1)=.TRUE.
      IF (LEN.LT.N) IDUM=ICHAIN(LEN+1)
      IF (LEN.LT.N) FLAG(IDUM)=.TRUE.
C
C  Reestablish agglomerability of objects in NN-chain.
C
      GOTO 300
  550 CALL NEXT(FLAG,I1,N)
      GOTO 200
C
  600 CONTINUE
      RETURN
  700 WRITE(6,750)
  750 FORMAT(' ERROR IN NN-CHAIN ROUTINE - INSUFFICIENT CHAIN SPACE'/)
      STOP
      END

      SUBROUTINE DETNN(DATA,FLAG,MEM,N,M,I1,I2,D)
C
C  Determine a nearest neighbour.
C
      REAL DATA(N,M),MEM(N)
      LOGICAL FLAG(N)
C
      DO 200 I=1,N
         IF (.NOT.FLAG(I)) GOTO 200
         DISS=0.
         DO 100 J=1,M
  100    DISS=DISS+(DATA(I1,J)-DATA(I,J))*(DATA(I1,J)-DATA(I,J))
         DISS=DISS*MEM(I)*MEM(I1)/(MEM(I1)+MEM(I))
         IF (DISS.GE.D) GOTO 200
            D=DISS
            I2=I
  200 CONTINUE
C
      RETURN
      END

      SUBROUTINE AGGLOM(I1,I2,D,DATA,MEM,FLAG,IA,IB,CRIT,NCL,N,M)
C     HISTORY
C     Bounds bug fix, Oct. 1990, F. Murtagh.
C     Insert line 25, "IF (I.EQ.0) GOTO 140".  Without /CHECK=(BOUNDS) option
C     the VAX-VMS compiler did not complain. 
C
C  Carry out an agglomeration.
C
      REAL MEM(N),DATA(N,M),CRIT(N)
      INTEGER IA(N),IB(N)
      LOGICAL FLAG(N)
      INTEGER O1,O2,LB,UB
C
C
      O1=MIN0(I1,I2)
      O2=MAX0(I1,I2)
      DO 100 J=1,M
         DATA(O1,J)=( MEM(O1)*DATA(O1,J)+MEM(O2)*DATA(O2,J) )
     X            / (MEM(O1)+MEM(O2))
         DATA(O2,J)=DATA(O1,J)
  100 CONTINUE
      NAGGL=N-NCL
      MEM(O1)=MEM(O1)+MEM(O2)
      FLAG(O2)=.FALSE.
C
C  Keep sorted list of criterion values: find 1st where new crit. fits.
C
      I=NAGGL-1
      IF (I.EQ.0) GOTO 140
  120 IF (D.GE.CRIT(I)) GOTO 140
      I=I-1
      IF (I.GE.1) GOTO 120
C
C  Arriving here must mean that D > all crit. values found so far.
C
      I=0
  140 CONTINUE
C
C  Now, shift rightwards from I+1 to AGGL-1 to make room for new crit.
C
      LB=I+1
      UB=NAGGL-1
      IF (LB.GT.UB) GOTO 180
      J=UB
  160 J1=J+1
      IA(J1)=IA(J)
      IB(J1)=IB(J)
      CRIT(J1)=CRIT(J)
      J=J-1
      IF (J.GE.LB) GOTO 160
  180 CONTINUE
      IA(LB)=O1
      IB(LB)=O2
      CRIT(LB)=D
C
      RETURN
      END

      SUBROUTINE NEXT(FLAG,I1,N)
C
C  Determine next agglomerable object/cluster.
C
      LOGICAL FLAG(N)
C
      NXT=I1+1
      IF (NXT.GT.N) GOTO 150
      DO 100 I=NXT,N
         IF (FLAG(I)) GOTO 500
  100 CONTINUE
  150 DO 200 I=1,I1
         IF (FLAG(I)) GOTO 500
  200 CONTINUE
C
      STOP
C
  500 I1=I
C
      RETURN
      END
