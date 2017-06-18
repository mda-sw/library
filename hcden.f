C+++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                 C
C  Construct a DENDROGRAM of the top 8 levels of  C
C  a HIERARCHIC CLUSTERING.                       C
C                                                 C
C  Parameters:                                    C
C                                                 C
C  IORDER, HEIGHT, CRITVAL: vectors of length LEV C
C          defining the dendrogram.               C
C          These are: the ordering of objects     C
C          along the bottom of the dendrogram     C
C          (IORDER); the height of the vertical   C
C          above each object, in ordinal values   C
C          (HEIGHT); and in real values (CRITVAL).C
C                                                 C
C  NOTE: these vectors MUST have been set up with C
C        LEV = 9 in the prior call to routine     C
C        HCASS.
C                                                 C
C  F. Murtagh, ESA/ESO/STECF, Garching, Feb. 1986.C
C                                                 C 
C-------------------------------------------------C
      SUBROUTINE HCDEN(LEV,IORDER,HEIGHT,CRITVAL)
      CHARACTER*80 LINE
      INTEGER IORDER(LEV),HEIGHT(LEV)
      REAL CRITVAL(LEV)
      INTEGER OUT(27,27)
      INTEGER UP,ACROSS,BLANK
      DATA UP,ACROSS,BLANK/'|','-',' '/
C
C
      DO I=1,27
        DO J=1,27
          OUT(I,J)=BLANK
        ENDDO
      ENDDO
C
C
      DO I=3,27,3
         I2=I/3
C
         J2=28-3*HEIGHT(I2)
         DO J=27,J2,-1
            OUT(J,I)=UP
         ENDDO
C
         DO K=I,3,-1
            I3=INT((K+2)/3)
            IF ( (28-HEIGHT(I3)*3).LT.J2) GOTO 100
            OUT(J2,K)=ACROSS
         ENDDO
  100    CONTINUE
C
      ENDDO
C
C
      IC=3
      DO I=1,27
      IF (I.EQ.IC+1) THEN
                   IDUM=IC/3
                   IDUM=9-IDUM
                   DO L=1,9
                      IF (HEIGHT(L).EQ.IDUM) GOTO 190
                   ENDDO
  190              IDUM=L
                   WRITE(6,200) CRITVAL(IDUM),(OUT(I,J),J=1,27)
                   IC=IC+3
                   ELSE
                   LINE = ' '
                   WRITE(6,210) (OUT(I,J),J=1,27)
      ENDIF
  200 FORMAT(1H ,8X,F12.2,4X,27A1)
  210 FORMAT(1H ,24X,27A1)
      ENDDO
      WRITE(6,250)
      WRITE(6,220)(IORDER(J),J=1,9)
      WRITE(6,250)
  220 FORMAT(1H ,24X,9I3)
      WRITE(6,230)
  230 FORMAT(1H ,13X,'CRITERION        CLUSTERS 1 TO 9')
      WRITE(6,240)
  240 FORMAT(1H ,13X,'VALUES.      (TOP 8 LEVELS OF HIERARCHY).')
  250 FORMAT(/)
C
C
      RETURN
      END
