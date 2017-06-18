      REAL DATA(5000,4),CRIT(5000),MEMBR(5000)
      INTEGER IA(5000),IB(5000)
      DIMENSION NN(5000),DISNN(5000)
c      REAL D(100000000)
      REAL D(5000*4999/2)
      LOGICAL FLAG(5000)
      integer*2 len

c      INTEGER ICLASS(100,4),HVALS(4)
c      REAL CRITVAL(9)
c      INTEGER IORDER(9),HEIGHT(9)
C IN ABOVE, 18=N, 16=M, 9=LEV, 153=N(N-1)/2.
C
C
c	OPEN(UNIT=21,STATUS='OLD',FILE='spectr.dat')
C
C
      N = 5000
      M = 4
      iseed = 3373
      DO I=1,N
        do j = 1, m
           data(i,j) = ran(iseed)
        enddo
      ENDDO
C
C
      LEN = (N*(N-1))/2
      IOPT=3
      CALL HC(N,M,LEN,IOPT,DATA,IA,IB,CRIT,MEMBR,NN,DISNN,FLAG,D)
c      do i = 1, n
c         write(6,*), ia(i), ib(i), crit(i)
c      enddo
C
C
c      LEV = 9
c      CALL HCASS(N,IA,IB,CRIT,LEV,ICLASS,HVALS,IORDER,CRITVAL,HEIGHT)
C
C
c      CALL HCDEN(LEV,IORDER,HEIGHT,CRITVAL)
C
C
      END

