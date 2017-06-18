	REAL DATA(150,4), g(150,3)
	REAL MEAN(4), sd(4), maxx(4)
C
	OPEN(UNIT=21,STATUS='OLD',FILE='iris.dat')
C
C	   Get input data.
C
         
	DO I = 1, 150
	   READ(21,100)(DATA(I,J),J=1,4)
           g(i,1) = 0
           g(i,2) = 0
           g(i,3) = 0
           if (i.le.50) g(i,1) = 1
           if (i.gt.50.and.i.le.100) g(i,2) = 1
           if (i.gt.100) g(i,3) = 1
  100	   FORMAT(4f4.1)
        ENDDO
        do j = 1, 4
           mean(j) = 0.
           maxx(j) = 0
           do i = 1, 150
              mean(j) = mean(j) + data(i,j)
              if (data(i,j).gt.maxx(j)) maxx(j)=data(i,j)
           enddo
           mean(j) = mean(j)/150.
         enddo
         write(6,*) ' Means: ', (mean(k),k=1,4)
         do j = 1, 4
            sd(j) = 0.
            do i = 1, 150
               sd(j) = sd(j) + (data(i,j)-mean(j))**2
            enddo
            sd(j) = sd(j)/sqrt(150.)
         enddo
         write(6,*) ' SDs: ', (sd(k),k=1,4)
         do i = 1, 150
            do j = 1, 4
c               data(i,j) = (data(i,j)-mean(j))/sd(j)
                data(i,j)=data(i,j)/maxx(j)
            enddo
          enddo
         do i = 1, 150
          write(20,140) (data(i,j),j=1,4),(g(i,k),k=1,3)
         enddo
  140    format(7f10.2)
c
	END
