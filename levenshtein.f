C Determine Weighted Levenshtein (or Edit) Distance.
	CHARACTER*1 A(6), B(6)
	INTEGER*4   D(7,7), DIR(7,7)
	DATA A/'M','U','L','L','E','R'/
	DATA B/'H','U','L','E','R','T'/
	N=6    ! len(A)
	M=6    ! len(B)
C----------------------------------------------------------------------
C ALGORITHM.  Ref: T. Kohonen, "Self-Organization & Associative Memory",
C Springer-Verlag, 1988, p. 66. (Here, notation slightly simplified.)
C D(0,0)=0
C for i=1 to len(A) do D(i,0)=D(i-1,0)+r(Ai)     ! r(Ai)=1
C for j=1 to len(B) do D(0,j)=D(0,j-1)+q(Bj)     ! q(Bj)=1
C for i=1 to len(A) do
C     for j=1 to len(B) do
C         m1=D(i-1,j-1)+p(Ai,Bj)  ! p=1 if Ai=Bj; else 0
C         m2=D(i,j-1)+q(Bj)
C         m3=D(i-1,j)+r(Ai)
C         D(i,j)=min(m1,m2,m3)
C     end
C WLD=D(len(A),len(B))
C-----------------------------------------------------------------------
	D(1,1)=0
	DO I=1,N
	   D(I+1,1)=D(I,1)+1
	ENDDO
	DO J=1,M
           D(1,J+1)=D(1,J)+1
	ENDDO
	DO I=1,N
	   DO J=1,M
	      IP=1
	      IF (A(I).EQ.B(J)) IP=0
	      M1=D(I,J)+IP
	      M2=D(I+1,J)+1
	      M3=D(I,J+1)+1
	      D(I+1,J+1)=MIN(M1,M2,M3)
C  NEW: det. directions followed; in case of ties, favour "diag" direction.
	      IF (M2.LE.M1.AND.M2.LE.M3) DIR(I+1,J+1)=2 ! "rightwards"
              IF (M3.LE.M1.AND.M3.LE.M2) DIR(I+1,J+1)=3 ! "downwards"
	      IF (M1.LE.M2.AND.M1.LE.M3) DIR(I+1,J+1)=1 ! "diag. right/down"
	   ENDDO
	ENDDO
	WRITE(6,10) (A(I),I=1,N)
   10	FORMAT('  String-1: ',70A1) 
	WRITE(6,20) (B(J),J=1,M)
   20	FORMAT('  String-2: ',70A1) 
	WRITE(6,*) ' Weighted Levenshtein Distance:',D(N+1,M+1)
	WRITE(6,*) ' Best match:'
C Trace-back:
	I=N+1
	J=M+1
  100   CONTINUE
	IF (I.LT.2.OR.J.LT.2) GOTO 200
	WRITE(6,*) A(I-1),B(J-1)
	IF (DIR(I,J).EQ.1) THEN
	   I=I-1
	   J=J-1
	   GOTO 100
	ENDIF
	IF (DIR(I,J).EQ.2) THEN
	   J=J-1
	   GOTO 100
	ENDIF
	IF (DIR(I,J).EQ.3) THEN
	   I=I-1
	   GOTO 100
	ENDIF
C Here, DIR(I,J)=0; we have finished the trace-back.
  200   CONTINUE
	END
