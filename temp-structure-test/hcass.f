C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, derive the assignments into clusters for the C
C  top LEV-1 levels of the hierarchy.                           C
C  Prepare also the required data for representing the          C
C  dendrogram of this top part of the hierarchy.                C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB, CRIT: vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  LEV:          number of clusters in largest partition.       C
C  HVALS:        vector of dim. LEV, used internally only.      C
C  ICLASS:       array of cluster assignments; dim. N by LEV.   C
C  IORDER, CRITVAL, HEIGHT: vectors describing the dendrogram,  C
C                all of dim. LEV.                               C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.          C
C                                                               C
C HISTORY                                                       C
C                                                               C
C Bounds bug fix, Oct. 1990, F. Murtagh.                        C
C Inserted line "IF (LOC.GT.LEV) GOTO 58" on line 48.  This was C
C occassioned by incorrect termination of this loop when I      C
C reached its (lower) extremity, i.e. N-LEV.  Without the       C
C /CHECK=(BOUNDS) option on VAX/VMS compilation, this inserted  C
C statement was not necessary.                                  C
C---------------------------------------------------------------C
      SUBROUTINE HCASS(N,IA,IB,CRIT,LEV,ICLASS,HVALS,IORDER,
     X        CRITVAL,HEIGHT)
      INTEGER IA(N),IB(N),ICLASS(N,LEV),HVALS(LEV),IORDER(LEV),
     X        HEIGHT(LEV)
      REAL CRIT(N),CRITVAL(LEV)
C
C  Pick out the clusters which the N objects belong to,
C  at levels N-2, N-3, ... N-LEV+1 of the hierarchy.
C  The clusters are identified by the lowest seq. no. of
C  their members.
C  There are 2, 3, ... LEV clusters, respectively, for the
C  above levels of the hierarchy.
C
      HVALS(1)=1
      HVALS(2)=IB(N-1)
      LOC=3
      DO 59 I=N-2,N-LEV,-1
         DO 52 J=1,LOC-1
            IF (IA(I).EQ.HVALS(J)) GOTO 54
  52     CONTINUE
         HVALS(LOC)=IA(I)
         LOC=LOC+1
  54     CONTINUE
         DO 56 J=1,LOC-1
            IF (IB(I).EQ.HVALS(J)) GOTO 58
  56     CONTINUE
         IF (LOC.GT.LEV) GOTO 58
         HVALS(LOC)=IB(I)
         LOC=LOC+1
  58     CONTINUE
  59  CONTINUE
C
      DO 400 LEVEL=N-LEV,N-2
         DO 200 I=1,N
            ICL=I
            DO 100 ILEV=1,LEVEL
  100       IF (IB(ILEV).EQ.ICL) ICL=IA(ILEV)
            NCL=N-LEVEL
            ICLASS(I,NCL-1)=ICL
  200    CONTINUE
  400  CONTINUE
C
      DO 120 I=1,N
      DO 120 J=1,LEV-1
      DO 110 K=2,LEV
      IF (ICLASS(I,J).NE.HVALS(K)) GOTO 110
         ICLASS(I,J)=K
         GOTO 120
  110 CONTINUE
  120 CONTINUE
C
      WRITE (6,450)
  450 FORMAT(4X,' SEQ NOS 2CL 3CL 4CL 5CL 6CL 7CL 8CL 9CL')
      WRITE (6,470)
  470 FORMAT(4X,' ------- --- --- --- --- --- --- --- --- ----')
      DO 500 I=1,N
      WRITE (6,600) I,(ICLASS(I,J),J=1,8) 
  600 FORMAT(I11,8I4)                    
  500 CONTINUE
C
C  Determine an ordering of the LEV clusters (at level LEV-1)
C  for later representation of the dendrogram.
C  These are stored in IORDER.
C  Determine the associated ordering of the criterion values
C  for the vertical lines in the dendrogram.
C  The ordinal values of these criterion values may be used in
C  preference, and these are stored in HEIGHT.
C  Finally, note that the LEV clusters are renamed so that they
C  have seq. nos. 1 to LEV.
C
      IORDER(1)=IA(N-1)
      IORDER(2)=IB(N-1)
      CRITVAL(1)=0.0
      CRITVAL(2)=CRIT(N-1)
      HEIGHT(1)=LEV
      HEIGHT(2)=LEV-1
      LOC=2
      DO 700 I=N-2,N-LEV+1,-1
         DO 650 J=1,LOC
            IF (IA(I).EQ.IORDER(J)) THEN
C              Shift rightwards and insert IB(I) beside IORDER(J):
               DO 630 K=LOC+1,J+1,-1
                  IORDER(K)=IORDER(K-1)
                  CRITVAL(K)=CRITVAL(K-1)
                  HEIGHT(K)=HEIGHT(K-1)
  630          CONTINUE
               IORDER(J+1)=IB(I)
                CRITVAL(J+1)=CRIT(I)
                HEIGHT(J+1)=I-(N-LEV)
               LOC=LOC+1
            ENDIF
  650   CONTINUE
  700 CONTINUE
      DO 705 I=1,LEV
         DO 703 J=1,LEV
            IF (HVALS(I).EQ.IORDER(J)) THEN
               IORDER(J)=I
               GOTO 705
            ENDIF
  703    CONTINUE
  705 CONTINUE
C
      RETURN
      END
