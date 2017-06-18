        SUBROUTINE DSLK(N,P,X,A,B,C,H,PH,NGR,LIMH,MAXCH,
     +    NDIMGR,LIM,IND,HEAP,PHEAP,MEM,DMEM,ILINK,JLINK,DLINK,CLUS,
     +    CLUSN,CLSEND,LGRDCD,CONST,CK,AMINP)
C--------------------------------------------------------------------
C         SINGLE LINK - MST PROGRAM USING DELTA METHOD.
C         F. JAMES ROHLF, FEB. 1977
C--------------------------------------------------------------------
C       INPUT PARAMETERS
C         P = No. of variables
C         N = No. of points
C         NDIMGR = maximum no. of dimensions to be used for the grids.
C         NGR = maximum no. of grid systems that will be used
C                    (this is a function of ndimgr and lgrdcd)
C         LIMH = length of 'hash' table.
C         X(P,N)  = input data matrix, p variables, n points
C         LGRDCD = 1 for Rabin grid system, = 2 for Yuval grid system
C       OUTPUT RESULTS
C         ILINK() = First member of each edge
C         JLINK() = second member of each edge
C         DLINK() = squared length of each edge
C                   (above three arrays are needed only if a mst is to
C                    be computed)
C         A() = list of objects in dendrogram order. During execution 
C               A(i) = -1 if object belongs to v, >0 if A(i) is a nea-
C               rest neighbor of object i.
C         B() = list of squared clustering levels for dendrogram. 
C               during execution, B(I) = distance from item I to item
C               A(I).
C       SCRATCH ARRAYS:
C         C() = list of items sorted by current cluster membership.
C         CLUS() = pointers to first member of each cluster in mem().
C         CLSEND() = pointer to last member of each cluster.
C         CLUSN() = size of each cluster.
C         MEM() = list of members of each  cluster. MEM(I) is the member
C                 that follows I. MEM(J).LT.0  if J is the last member of
C                 the cluster -MEM(J)
C         DMEM() = squared clustering levels (to be used with MEM() array
C         H() = work vector used by LAT1, LAT2, and MDIS2
C            length = LIMH
C         PH(,) = "   " dimension is NGR by N
C         IND() = work vector, length = NGR
C         LIM() = "  "    , length = NDIMGR
C         CK() = work vector to check for redundant tests, length=N
C         MAXCH() is a work vector used by SCALE,LAT1, and LAT2
C         HEAP() = a 'heap' of length NHEAP and ordered by the magnitude of
C                  B(HEAP(I))
C         PHEAP() = pointers to where an item can be found in the heap
C         note : C and HEAP can occupy the same storage locations.
C---------------------------------------------------------------------------
C
          INTEGER          P,NDIMGR,TB
          REAL*4           X(P,N),B(N),DMAX/9.99E30/,MAXCH(NDIMGR),
     +                     DMEM(N),DLINK(N),AMINP(P)
          INTEGER          A(N),C(N),R,IND(NGR),LIM(NDIMGR),H(LIMH),
     +                     PH(NGR,N),MEM(N),END,LSIZE,RSIZE,SIZE,SAVE,
     +                     PSAVE,UPTR,USIZE,USIZEP,ILINK(N),JLINK(N),
     +                     HEAP(N),PHEAP(N),CLUS(N),CLUSN(N),CLSEND(N),
     +                     CK(N),CURR
          COMMON/STOR/LASTFR,LIMLOT,LOT(3000)
C
C---------------------------------------------------------------------------
C
          FN = N
          FP = P
C
C                SCALE DATA AND FIND MAXIMUM FOR EACH VARIABLE
C
          CALL SCALE(X,N,P,MAXCH,VOL,AMINP)
          NEDGE = 0
          LASTFR = 1
          DO 120 I=1,N
             MEM(I) = -I
             CLUS(I) = I
             CLUSN(I) = 1
             CLSEND(I) = I
             C(I) = I
             PHEAP(I) = 0
             A(I) = 0
             B(I) = DMAX
             CK(I) = 0
  120     CONTINUE
          CURR = 1
          NCLUS = N
C
C                 INITIAL ESTIMATE OF DELTA
C
          CALL MDIS0(X,P,N,DMIN)
C
C                 CONSTRAIN DEL2 SO IT CANNOT BE TOO SMALL
C
          DEL2 = AMAX1(DMIN,CONST*(VOL/FN)**(2./FP))
          USIZE = 1
C
C                 GO TO SETUP GRIDS.
C
          GO TO 300
C
C                 WHILE NCLUS.GT.1 ...
C
 200      CONTINUE
C
C                 RESTORE MEMORY
C
          LASTFR = 1
C
C                 SETUP A TO POINT TO THE START OF EACH FRAGMENT AS LISTED IN C
C
          L = 1
          A(1) = L
          NEXT = CLUS(1)
          SIZE = CLUSN(1)
          DO 210  J=1,SIZE
             C(L) = NEXT
             L = L + 1
             NEXT = MEM(NEXT)
 210      CONTINUE
          DO 220  I = 2,NCLUS
             A(I) = L
             NEXT = CLUS(I)
             SIZE = CLUSN(I)
             DO 215  J = 1,SIZE
                C(L) = NEXT
                L = L + 1
                NEXT = MEM(NEXT)
 215         CONTINUE
 220      CONTINUE
C
C               MIN DIST IS RETURNED IN DMIN
C
 240      CALL MDIS1(X,P,N,C,DEL2,A,CLUSN,NCLUS)
 241      CONTINUE
C
C               RE-ETABLISH A AND B
C
          USIZE = CLUSN(1)
          DO 255  J=1,USIZE
             JJ = C(J)
             A(JJ) = -1
 255      CONTINUE
          USIZEP = USIZE + 1
          DO 265  J = USIZEP,N
             JJ = C(J)
             PHEAP(JJ) = 0
             A(JJ) = 0
             B(JJ) = DMAX
 265      CONTINUE
C
C               MAKE A SET OF GRIDS
C
 300      CONTINUE
          DEL = SQRT(DEL2*1.00005)
C
C               USE RABIN OR YUVAL GRID SYSTEM
C
          GO TO (320,330),LGRDCD
  320     CALL LAT1(N,C,P,X,DEL,H,LIMH,NGR,PH,MAXCH,NDIMGR,LIM,IND,
     +              USIZE,KMAX)
          GO TO 350
  330     CALL LAT2(N,C,P,X,DEL,H,LIMH,NGR,PH,MAXCH,NDIMGR,LIM,IND,
     +              USIZE,KMAX)
 350      CONTINUE
C
C               INIT HEAP MEMORY IN 'HEAP'
C
          NHEAP = 0
C
C               START WITH CLUSTER No 1
C
          IIC = 1
          UPTR = 1
C
C               WHILE THERE EXIST EDGES .LE. DEL, CONSIDER PIVOTAL CLUSTERS
C               FROM IIC = 1 TO NCLUS-1
C
C               FIND MIN DIST BETWEEN CLUS(IIC) AND CLUS (>IIC) USING CLASSES
C               IN H
C
 400      CONTINUE
C
C               CLOSE PAIRS OF ITEMS SAVED IN HEAP, IGNORE D**2'S > DEL2
C
          CALL MDIS2(X,P,N,A,B,HEAP,NHEAP,PHEAP,PH,H,NGR,LIMH,UPTR,MEM,
     +               USIZE,DEL2,KMAX,CK,CURR)
C
C                SEARCH FOR THE SMALLEST EDGE IN THE LIST OF SHORT EDGES.
C
           IF(NHEAP.LT.0) GO TO 600
           NEW = HEAP(1)
C
C                 DELETE THIS EDGE FROM HEAP
C
           CALL DRHEAP(HEAP,PHEAP,B,NHEAP)
C
C                 II IS A MEMBER OF CLUSTER IIC ON THE LEFT,
C                 NEW IS A MEMBER OF A CLUSTER TO ITS RIGHT.
C
           II = A(NEW)
           DMIN = B(NEW)
           NEDGE = NEDGE + 1
C
C                 ADD AN EDGE TO THE MST (NOT NEEDED FOR SINGLE LINK)
C
           ILINK(NEDGE) = II
           JLINK(NEDGE) = NEW
           DLINK(NEDGE) = DMIN
C
C                 UPDATE POINTERS
C                 FIND END OF RIGHT CLUSTER.
C
           NEXT = NEW
           NJ = 0
 510       A(NEXT) = -1
             NJ = NJ + 1
             IF(PHEAP(NEXT).GT.0) CALL DLHEAP(NEXT,HEAP,PHEAP,B,NHEAP)
             NEXT = MEM(NEXT)
           IF (NEXT.GT.0) GO TO 510
C
C                 END FOUND, GET SIZE
C
           JJC = -NEXT
           RSIZE = CLUSN(JJC)
           LSIZE = CLUSN(IIC)
           JCLUS = CLUS(JJC)
           SIZE = LSIZE + RSIZE
           CLUSN(IIC) = SIZE
C
C                 LINK THE TWO CLUSTERS TOGETHER (SINGLE LINK)
C
           MEM(CLSEND(IIC)) = JCLUS
           DMEM(CLSEND(IIC)) = DMIN
           CLSEND(IIC) = CLSEND(JJC)
           MEM(CLSEND(IIC)) = -IIC
C-------------------------------------
C
C                 CHECK IF DONE
C 
           IF (NEDGE.NE.N-1) GO TO 530
C
C                 UNTANGLE TREE MATRIX
C
              IPT =  CLUS(1)
              NM1 = N - 1
              DO 520  I = 1,NM1
                 A(I) = IPT
                 B(I) = DMEM(IPT)
                 IPT = MEM(IPT)
  520         CONTINUE
              A(N) = IPT
C
C                EXIT
C
              RETURN
C
C----------------------------------------------------------------------
C
C                CHECK IF THERE ARE ANY MORE MEMBERS OF JJC IN THE HEAP
C
  530      IF (RSIZE.EQ.NJ) GO TO 540
               NNJ = RSIZE - NJ
               NEXT = JCLUS
               DO 535  J = 1,NNJ
                  A(NEXT) = -1
                  IF (PHEAP(NEXT).GT.0) CALL DLHEAP(NEXT,HEAP,PHEAP,B,
     +                                              NHEAP)
                  NEXT = MEM(NEXT)
 535           CONTINUE
 540       CONTINUE
C
C               DELETE CLUSTER JJC, MOVE CLUSTER NCLUS UP TO TAKE ITS PLACE.
C
           IF(JJC.EQ.NCLUS) GO TO 550
             CLUS(JJC) = CLUS(NCLUS)
             CLSEND(JJC) = CLSEND(NCLUS)
             CLUSN(JJC) = CLUSN(NCLUS)
             MEM(CLSEND(JJC)) = -JJC
 550       NCLUS = NCLUS - 1
C
C               CHECK IF AT END OF CURRENT STAGE.
C
           IF(IIC.GE.NCLUS) GO TO 200
C
C               NO MORE CLUSTERS TO THE RIGHT TO BE CHECKED
C
           UPTR = JCLUS
           USIZE = RSIZE
           GO TO 400
C
C               ADVANCE TO A NEW PIVOTAL CLUSTER
C
 600       CONTINUE
           IIC = IIC + 1
           IF ( IIC.EQ.NCLUS) GO TO 200
C
C               STORE CLUSTER IIC IN SET U
C
           UPTR = CLUS(IIC)
           USIZE = CLUSN(IIC)
C
C               MARK MEMBERS OF CLUSTER IIC AS BELONGING TO SET V
C
           NEXT = UPTR
           DO 620  I = 1,USIZE
              A(NEXT) = -1
              NEXT = MEM(NEXT)
 620       CONTINUE
        GO TO 400
        END

        SUBROUTINE SCALE(X,N,P,MAX,VOL,AMINP)
C-------------------------------------------------------------------------
C               SCALE DATA MATRIX FOR USE BY DELTA METHODS
C               F. JAMES ROHLF
C-------------------------------------------------------------------------
C      INPUT
C          N = NUMBER OF OBJECTS
C          P = NO. OF VARIABLES
C          X = DATA MATRIX
C      OUTPUT
C          MATRIX WILL HAVE A MINIMUM OF ZERO FOR EACH VARIABLE
C          MAX = RANGE OF EACH VARIABLE
C          VOL = PRODUCT OF RANGES
C-------------------------------------------------------------------------
C           
        INTEGER      P,N,I,J
        REAL         X(P,N),MAX(P),VOL,AMAX,AMIN
	REAL         AMINP(P)
C
C-------------------------------------------------------------------------
C          SCALE DATA
C
        VOL = 1.
        DO I=1,P
           AMIN = X(I,1)
           AMAX = X(I,1)
           DO J=2,N
              IF (X(I,J).LT.AMIN)  AMIN = X(I,J)
              IF (X(I,J).GT.AMAX)  AMAX = X(I,J)
           END DO
	   AMINP(I) = AMIN
           DO J=1,N
              X(I,J) = X(I,J) - AMIN
           END DO
           MAX(I) = AMAX - AMIN
           VOL = VOL * MAX(I)
        END DO
        RETURN
        END

        SUBROUTINE MDIS0(X,P,N,DMIN)
C---------------------------------------------------------------------------
C                 FIND THE MAXIMA OF N DISTANCES SAMPLED AT RANDOM
C                 INITIAL ESTIMATE OF DELTA
C                 F. JAMES ROHLF
C---------------------------------------------------------------------------
C
        INTEGER           P,N
        REAL              X(P,N)
        INTEGER           MULT/1220703123/,ISEED/1643242315/
	EXTERNAL	  IOVHAND
C---------------------------------------------------------------------------
C
ccc	CALL LIB$ESTABLISH(IOVHAND)
        DMIN = 99. E30
        DO I=1,N
           ISEED = IABS(ISEED * MULT) + 1
           II = MOD(ISEED,N) + 1
 100       ISEED = IABS(ISEED * MULT)
           JJ = MOD(ISEED,N) + 1
           IF (II.EQ.JJ)  GO TO 100
           SUM = 0.0
           DO K=1,P
              SUM = SUM + (X(K,II) - X(K,JJ)) * (X(K,II) - X(K,JJ))
              IF (SUM.GT.DMIN)  GO TO 170
           END DO
           IF (SUM.GT.0.0)  DMIN = SUM
 170       CONTINUE
        END DO
ccc	CALL LIB$REVERT()
        RETURN 
        END

        SUBROUTINE MDIS1(X,P,N,C,DMIN,PC,CLUSN,NCLUS)
C----------------------------------------------------------------------------
C                  FIND MINIMUM DISTANCE BETWEEN CLUSTERS USING RANDOM SAMPLE
C                  F. JAMES ROHLF
C----------------------------------------------------------------------------
C
        INTEGER       P,C(N),PC(NCLUS),CLUSN(NCLUS)
        INTEGER       MULT/1220703123/,ISEED/1643242315/
        REAL          X(P,N)
        DATA          DMAX/99.E30/
C----------------------------------------------------------------------------
C
	EXTERNAL	  IOVHAND
C---------------------------------------------------------------------------
C
ccc	CALL LIB$ESTABLISH(IOVHAND)
        DMIN = DMAX
        DO 175  I=1,N
C
C                SELECT A PAIR OF CLUSTERS
C
 120       IC = MOD(I,NCLUS) + 1
           JC = MOD(IC,NCLUS) + 1
C
C                SELECT A MEMBER AT 'RANDOM' FROM EACH
C
           IPT = PC(IC)
           IF (CLUSN(IC).EQ.1)  GO TO 130
           ISEED = IABS (ISEED * MULT) + 1
           IPT = IPT + MOD(ISEED,CLUSN(IC))
 130       JPT = PC(JC)
           IF (CLUSN(JC).EQ.1)  GO TO 135
           ISEED = IABS(ISEED * MULT)
           JPT = JPT + MOD(ISEED,CLUSN(JC))
 135       II = C(IPT)
           JJ = C(JPT)
           SUM = 0.0
           DO 150  K=1,P
              SUM = SUM + (X(K,II) - X(K,JJ)) * (X(K,II) - X(K,JJ))
              IF (SUM.GT.DMIN)  GO TO 170
 150       CONTINUE
           DMIN = SUM
 170       CONTINUE
 175    CONTINUE
ccc	CALL LIB$REVERT()
        RETURN
        END

        SUBROUTINE MDIS2(X,P,N,A,B,HEAP,NHEAP,PHEAP,PH,H,LIMK,
     +                   LIMH,ISTART,MEM,NLEFT,DCRIT,KMAX,CK,CURR)
C---------------------------------------------------------------------------
C                  FIND NONMEMBERS OF V THAT HAVE DISTANCES < DCRIT TO
C                  MEMBERS OF V.
C                  SAVE ITEMS WITH SMALL DISTANCES ON LIST HEAP
C                  USES A HEAP STORE LIST OF SHORTEST EDGES.
C                  F. JAMES ROHLF
C---------------------------------------------------------------------------
C
        INTEGER       P,A(N),AITEM
        REAL          X(P,N),B(N)
C        INTEGER       HEAP(NHEAP),PHEAP(N)
        INTEGER       HEAP(1),PHEAP(N)
        INTEGER       MEM(N)
        INTEGER       PH(LIMK,N),H(LIMH)
        INTEGER       BACK
        INTEGER       CK(N),CURR
        COMMON/STOR/LASTFR,LIMLOT,LOT(3000)
C---------------------------------------------------------------------------
C
        II = ISTART
        DO 440  I=1,NLEFT
           CURR = CURR + 1
           DO 420  J=1,KMAX
              IH = PH(J,II)
              IF (IH.EQ.0)  GO TO 420
              IPTH = H(IH)
              BACK = 0
              IPT = IPTH
 100          IF (IPT.EQ.0) GO TO 420
 115          ITEM = LOT(IPT)
C
C                 IF ITEM IS A MEMBER OF V, DELETE IT
C
              AITEM = A(ITEM)
              IF(AITEM)  116,126,125
C
C                 DELETE NODE
C
 116            IPT = LOT(IPT+1)
C
C                 FIX UP POINTERS.
C
                IF(BACK.EQ.0) GO TO 117
                LOT(BACK+1) = IPT
                GO TO 100
 117            H(IH) = IPT
                GO TO 100
 125         IF (AITEM.EQ.II)  GO TO 170
C
C                COMPUTE DISTANCE.
C
 126         CONTINUE
             IF (CK(ITEM).EQ.CURR)  GO TO 170
             CK(ITEM) = CURR
             SUM = 0.0
             DO 150  K=1,P
                SUM = SUM + (X(K,II) - X(K,ITEM)) * (X(K,II)-X(K,ITEM))
                IF (SUM.GT.DCRIT)  GO TO 170
                IF (SUM.GE.B(ITEM))  GO TO 170
 150         CONTINUE
C
C               WE HAVE JUST FOUND A NEW LOWEST VALUE.
C               IS IT THE FIRST TIME FOR THIS ITEM ?
C
             IF (AITEM.GT.0)  GO TO 165
C
C               YES, ADD IT TO THE HEAP
C
             NHEAP = NHEAP + 1
             HEAP(NHEAP) = ITEM
             PHEAP(ITEM) = NHEAP
 165         B(ITEM) = SUM
             A(ITEM) = II
C
C               MOVE ITEM IN HEAP TO REFLECT ITS NEW (SMALLER) VALUE
C
             JHPT = PHEAP(ITEM)
 168         IHPT = JHPT / 2
             IF(IHPT.LT.1)  GO TO 170
             IHEAP = HEAP(IHPT)
             IF(B(IHEAP).LE.B(ITEM))  GO TO 170
C
C               ELSE SWAP
C
             HEAP(IHPT) = ITEM
             HEAP(JHPT) = IHEAP
             PHEAP(ITEM) = IHPT
             PHEAP(IHEAP) = JHPT
             JHPT = IHPT
             GO TO 168
 170      BACK = IPT
          IPT = LOT(IPT+1)
        GO TO 100
 420    CONTINUE
        II = MEM(II)
 440    CONTINUE
        RETURN
        END

        SUBROUTINE LAT1(N,C,P,X,DEL,H,LIMH,LIMK,PH,MAXCH,PLAT,LIM,IND,
     +                  NT,KMAX)
C---------------------------------------------------------------------------
C                  COMPUTE PROJECTIONS ONTO A LATTICE WITH GRID SIZE = 2*DEL
C                  (RABIN 1976,1977 METHOD)
C                  F. JAMES ROHLF
C---------------------------------------------------------------------------
C
        INTEGER       P,PLAT,C(N)
        REAL          X(P,N),MAXCH(PLAT)
        INTEGER       HASH/1220703123/
        INTEGER       H(LIMH),PH(LIMK,N)
        INTEGER       LIM(PLAT),IND(LIMK)
        COMMON/STOR/LASTFR,LIMLOT,LOT(3000)
C--------------------------------------------------------------------------
C
        R2DEL = 0.5 / DEL
        NONE = 0
        MXCELL = 1
        DO 110 J=1,PLAT
           IF (MAXCH(J) * R2DEL.GT.1.0)  GO TO 105
           LIM(J) = 1
           NONE = NONE + 1
           GO TO 110
 105       LIM(J) = (MAXCH(J) + DEL) * R2DEL + 1.
           MXCELL = MXCELL * LIM(J)
 110    CONTINUE
C
C                      CHECK IF ALL POINTS ARE IN THE SAME CELL
C
        IF (NONE.EQ.PLAT)  GO TO 600
        MXCELL = MIN0(MXCELL * 2 ** (PLAT-NONE),LIMH)
        DO 120 I=1,MXCELL
 120       H(I) = 0
C
C                      PROCESS _V FIRST
C
        NTP1 = NT + 1
        DO 390 I=NTP1,N
           II = C(I)
           KIND=1
           IND(1) = 0
           DO 250 J=1,PLAT
              LJ = LIM(J)
              IF (LJ.EQ.1)  GO TO 250
              IND1 = X(J,II) * R2DEL
              IND2 = (X(J,II) + DEL) * R2DEL
              DO 230 K=1,KIND
                 IND(KIND+K) = (IND(K) * LJ + IND2) * 2 + 1
                 IND(K) = (IND(K) * LJ + IND1) * 2
 230          CONTINUE
              KIND = KIND * 2
 250      CONTINUE
          DO 380 K=1,KIND
             IH=IND(K) + 1
             IF (IH.GT.MXCELL)  IH = MOD(IABS(IH * HASH),MXCELL) + 1
             PH(K,II) = IH
C
C                      LIST ITEMS IN EACH CLASS
C
329          IPT = H(IH)
C
C                      PUSH DOWN LIST OF ITEMS
C
             IF (LASTFR.GT.LIMLOT)  STOP 991
             LOT(LASTFR+1) = IPT
             H(IH) = LASTFR
             LOT(LASTFR) = II
             LASTFR = LASTFR + 2
 380      CONTINUE
 390    CONTINUE
C
C                     NOW DO MEMBERS OF V
C
        DO 490 I=1,NT
           II = C(I)
           KIND = 1
           IND(1) = 0
           DO 450 J=1,PLAT
              LJ = LIM(J)
              IF(LJ.EQ.1) GO TO 450
              IND1 = X(J,II) * R2DEL
              IND2 = (X(J,II) + DEL) * R2DEL
              DO 430 K=1,KIND
                 IND(KIND+K) = (IND(K) * LJ + IND2) * 2 + 1
                 IND(K) = (IND(K) * LJ + IND1) * 2
 430          CONTINUE
              KIND = KIND * 2
 450       CONTINUE
           DO 480 K= 1,KIND
              IH = IND(K) + 1
              IF (IH.GT.MXCELL)  IH = MOD(IABS(IH*HASH),MXCELL) + 1
              PH(K,II) = 0
              IF (H(IH).NE.0) PH(K,II) = IH
 480       CONTINUE
 490     CONTINUE
         KMAX = KIND
         RETURN
C
C                        ALL POINTS IN A SINGLE CELL
C
 600     NTP1 = NT + 1
         IPT = LASTFR
         DO 620 I = NTP1,N
            II = C(I)
            PH(1,II) = 1
            LOT(LASTFR+1) = LASTFR - 2
            LOT(LASTFR) = II
C
C                        NO DANGER OF MEMORY OVERFLOW IN THIS CASE
C
            LASTFR = LASTFR + 2
 620     CONTINUE
         LOT(IPT+1) = 0
         H(1) = LASTFR  - 2
         DO 640 I = 1,NT
            II = C(I)
            PH(1,II) = 1
 640     CONTINUE
         KMAX = 1
         RETURN
         END

        SUBROUTINE LAT2(N,C,P,X,DEL,H,LIMH,LIMK,PH,MAXCH,PLAT,LIM,IND,
     +                  NT,KMAX)
C-------------------------------------------------------------------------
C                  COMPUTE PROJECTIONS ONTO A LATTICE WITH GRID SIZE = 
C                  (P+1) * DEL
C                  METHOD OF YUVAL (1975,1976)
C                  LIMK = P + 1
C                  F. JAMES ROHLF
C-------------------------------------------------------------------------
C
        INTEGER        P,PLAT,C(N),PLP1,PLP2
        REAL           X(P,N),MAXCH(PLAT)
        INTEGER        HASH/1220703123/
        INTEGER        H(LIMH),PH(LIMK,N)
        INTEGER        LIM(PLAT),IND(LIMK)
        COMMON/STOR/LASTFR,LIMLOT,LOT(3000)
C--------------------------------------------------------------------------
C                  CHECK IF WE CAN REDUCE THE NUMBER OF DIMENSIONS
C
        R2DEL = 0.5 / DEL
        NONE = 0
        DO 80 J=1,PLAT
           LIM(J) = 0
           IF (MAXCH(J) * R2DEL.GT.1.0)  GO TO 80
           NONE = NONE + 1
           LIM(J) = 1
  80    CONTINUE
        IF (NONE.EQ.PLAT)  GO TO 600
        FPL = PLAT - NONE
        FPLP1 = FPL + 1.0
        PLP1 = PLAT - NONE + 1
        PLP2 = PLAT - NONE + 2
        RPDEL = 1./(FPLP1 * DEL)
        NONE2 = 0
        MXCELL = 1
        DO 110 J=1,PLAT
           IF(LIM(J).EQ.1)  GO TO 110
           IF ( MAXCH(J) * RPDEL.GT.1.0)  GO TO 105
           LIM(J) = 1
           NONE2 = NONE2 + 1
           GO TO 110
 105       LIM(J) = (MAXCH(J) + FPL * DEL) * RPDEL + 1.
           MXCELL = MXCELL * LIM(J)
 110    CONTINUE
        IF(NONE + NONE2.EQ.PLAT)  GO TO  600
        MXCELL = MIN0(MXCELL * (PLAT - NONE - NONE2 + 1),LIMH) !!????
        DO 120 I=1,MXCELL
 120       H(I) = 0
C
C                   PROCESS NOT V FIRST
C 
           NTP1 = NT + 1
           DO 390 I=NTP1,N
              II = C(I)
              DO 250 J=1,PLP1
                 FJ = J - 1
                 IND2 = 0
                 DO 230 IP=1,PLAT
                    IF(LIM(IP).EQ.1) GO TO 230
                    IND1 = (X(IP,II) + FJ * DEL) *RPDEL
                    IND2 = IND2 * LIM(IP) + IND1
 230             CONTINUE
                 IH = IND2 * PLP2 + J
                 IF (IH.GT.MXCELL)  IH = MOD(IABS(IH *HASH),MXCELL) + 1
                 PH(J,II) = IH
C
C                   LIST ITEMS IN EACH CLASS
C
                 IPT = H(IH)
                 IF (LASTFR.GT.LIMLOT)  STOP 991
                 LOT(LASTFR+1) = IPT
                 H(IH) = LASTFR
                 LOT(LASTFR) = II
                 LASTFR = LASTFR + 2
 250          CONTINUE
 390       CONTINUE
C
C                     NOW DO MEMBERS OF V
C
           DO 490 I=1,NT
              II = C(I)
              DO 430  J=1,PLP1
                 FJ = J - 1
                 IND2 = 0
                 DO 450 IP=1,PLAT
                    IF(LIM(IP).EQ.1)  GO TO 450
                    IND1 = (X(IP,II) + FJ * DEL) * RPDEL
                    IND2 = IND2 * LIM(IP) + IND1
 450             CONTINUE
                 IH = IND2 * PLP2 + J
                 IF(IH.GT.MXCELL)  IH = MOD(IABS(IH * HASH),MXCELL) + 1
                 PH(J,II) = 0
                 IF (H(IH).NE.0)  PH(J,II) = IH
 430          CONTINUE
 490       CONTINUE
           KMAX = PLP1
           RETURN
C
C                        ALL POINTS IN A SINGLE CELL
C
 600       NTP1 = NT + 1
           IPT = LASTFR
           DO 620 I=NTP1,N
              II = C(I)
              PH(1,II) = 1
              LOT(LASTFR+1) = LASTFR - 2
              LOT(LASTFR) = II
C
C                        NO DANGER OF MEMORY OVERFLOW IN THIS CASE
C
              LASTFR = LASTFR + 2
 620       CONTINUE
           LOT(IPT+1) = 0
           H(1) = LASTFR  - 2
           DO 640 I = 1,NT
              II = C(I)
              PH(1,II) = 1
 640       CONTINUE
           KMAX = 1
           RETURN
           END

           SUBROUTINE DLHEAP(ITEM,HEAP,PHEAP,VALUE,LAST)
C----------------------------------------------------------------------------
C                     DELETE ITEM FROM HEAP
C                     F. JAMES ROHLF
C----------------------------------------------------------------------------
C                     ARGUMENTS :
C                            HEAP = A 'HEAP'
C                            ITEM = ITEM TO BE DELETED FROM THE HEAP
C                            PHEAP = POINTER TO THE LOCATION OF OBJECT IN HEAP
C                            VALUE = MAGNITUDE OF OBJECT I
C                            LAST = CURRENT LENGTH OF THE HEAP
C-----------------------------------------------------------------------------
C
           INTEGER       HEAP(LAST),PHEAP(3)
           REAL          VALUE(3)
C-----------------------------------------------------------------------------
           LOC = PHEAP(ITEM)
           PHEAP(ITEM) = 0
           ITEML = HEAP(LAST)
           LAST = LAST - 1
           IF (LOC.GT.LAST)  RETURN
           HEAP(LOC) = ITEML
           PHEAP(ITEML) = LOC
C
C                   CHECK WICH WAY WE HAVE TO MOVE THE REPLACEMENT ITEM
C
           IF(VALUE(ITEML) - VALUE(ITEM))  100,200,300
C
C                   MOVE IT UP
C
 100       J = LOC
           JHEAP = HEAP(J)
 120          I = J / 2
              IF(I.LT.1)  RETURN
              IHEAP = HEAP(I)
              IF (VALUE(IHEAP).LE.VALUE(JHEAP))  RETURN
C
C                   SWAP
C
              HEAP(I) = JHEAP
              HEAP(J) = IHEAP
              PHEAP(JHEAP)  = I
              PHEAP(IHEAP)  = J
              J = I
           GO TO 120
C
C                   NO CHANGE
C
 200       RETURN
C
C                   MOVE IT DOWN
C
 300       I = LOC
           IHEAP = HEAP(I)
 320          J = 2 * I
              JHEAP = HEAP(J)
              IF (J - LAST) 325,330,200
 325          IF (VALUE(JHEAP).LE.VALUE(HEAP(J+1)))  GO TO 330
              J = J + 1
              JHEAP = HEAP(J)
 330          IF(VALUE(IHEAP).LE.VALUE(JHEAP))  RETURN
C
C                   SWAP
C
 350          HEAP(I) = JHEAP
              HEAP(J) = IHEAP
              PHEAP(JHEAP) = I
              PHEAP(IHEAP) = J
              I = J
           GO TO 320
           END

           SUBROUTINE DRHEAP(HEAP,PHEAP,VALUE,LAST)
C-----------------------------------------------------------------------------
C                     DELETE ROOT OF HEAP
C                     F. JAMES ROHLF
C-----------------------------------------------------------------------------
C                     ARGUMENTS :
C                            HEAP = A 'HEAP'
C                            PHEAP = POINTER TO LOCATION OF OBJECT IN HEAP
C                            VALUE = MAGNITUDE OF OBJECT I
C                            LAST = CURRENT LENGTH OF THE HEAP
C-----------------------------------------------------------------------------
           INTEGER       HEAP(1), PHEAP(3)
           REAL          VALUE(3)
C-----------------------------------------------------------------------------
           PHEAP(HEAP(1)) = 0
           ITEML = HEAP(LAST)
           LAST = LAST - 1
           IF (LAST.EQ.0)   RETURN
C
C                     PUT LAST ITEM AT THE ROOT
C
           HEAP(1) = ITEML
           PHEAP(ITEML) = 1
C
C                     NOW SHUFFLE HEAP TO FIX THINGS UP
C
           I = 1
           IHEAP = HEAP(I)
 120       CONTINUE
              J = 2 * I
              JHEAP = HEAP(J)
              IF(J - LAST)  125,130,900
 125          IF(VALUE(JHEAP).LE.VALUE(HEAP(J + 1)))  GO TO 130
              J = J + 1
              JHEAP = HEAP(J)
 130          IF (VALUE(IHEAP).LE.VALUE(JHEAP))  RETURN
C
C                     ELSE   SWAP
C
 150          HEAP(I) = JHEAP
              HEAP(J) = IHEAP
              PHEAP(JHEAP) = I
              PHEAP(IHEAP) = J
              I = J
           GO TO 120
 900       RETURN
           END

ccc	INTEGER*4 FUNCTION IOVHAND(SIGNAL,MECH)
C
C	INTERRUPT HANDLER FOR SS*INTOVF CONDITION
C
ccc	IMPLICIT INTEGER*4	(A-Z)
ccc	INTEGER*4		SIGNAL(*),MECH(5),ERR_ARRAY(9)
ccc	EXTERNAL		SS$_INTOVF, SS$_RESIGNAL
ccc	EXTERNAL		SS$_CONTINUE
C
ccc	IF (SIGNAL(2) .NE. %LOC(SS$_INTOVF)) THEN
ccc	  IOVHAND = %LOC(SS$_RESIGNAL)
ccc	ELSE
ccc	  IOVHAND = %LOC(SS$_CONTINUE)
ccc	ENDIF
ccc	RETURN
ccc	END
