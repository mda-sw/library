c     PT. MATCHING, GIVEN ROTN. ANGLE & XLN. VECT. V. 1.
c     F. MURTAGH, ST-ECF, GARCHING/MUNICH, NOV. 1991. (C)
      dimension xdata(20000), ydata(20000)
      character infil1*20, infil2*20
c
      write (6,*)
     .  ' ------------------------------------------------------------'
      write (6,*) ' PT. MATCHING, GIVEN ROTN. ANGLE & XLN. VECT., V.1.'
      write (6,*) ' Pt. matching; require input file names (ascii);'
      write (6,*) ' Free format: int.#, int.#, coordx, coordy, magn.'
      write (6,*)
     .  ' ------------------------------------------------------------'
      write (6,21)
      read  (5,30) infil1
      write (6,22)
      read  (5,30) infil2
      open  (unit=21,status='old',file=infil1)
      open  (unit=22,status='old',file=infil2)
 21   format(' Input file 1, first point set:  ',$)
 22   format(' Input file 2, second point set: ',$)
 30   format(a20)
      write (6,*)
     .  ' ------------------------------------------------------------'
      write (6,31)
      read  (5,*) xang
      write (6,32) 
      read  (5,*) c1, c2
 31   format(' Input rotation angle (for second pt. set):   ',$)
 32   format(' Input translation vector (1st->2nd pt. set): ',$)
      write (6,*)
     .  ' ------------------------------------------------------------'
      write (6,*) ' Output to unit 27.'
      write (6,*)
c
c---- Will offset and rotate xdata to align it with ydata. First deg to rad.
      xang = 360.0 - xang
      xang = (2.0*3.1415926)*xang/360.0
      cxang = cos(xang)
      sxang = sin(xang)
      ndim = 2
      iii = 0
c     List x or first list: set of pts.
      do i = 1, 5000
         read(21,*,end=50) is,in,x1,x2,xmag
         iii = iii + 1
c        Offset xdata to bring it into line with ydata; then rotate it.
         x1 = x1 - c1
         x2 = x2 - c2
         xdata(2*i-1) = cxang*x1+sxang*x2 
         xdata(2*i)  = -sxang*x1+cxang*x2 
      enddo
      write(6,*) ' Insuff. storage for input data,',
     .           ' increase in program.  Aborting.'
      stop
 50   continue
      n1 = iii - 1
      write (6,*) ' First list: # pts. (total):           ',n1
c
      iii = 0
c     List y or second list: set of pts.
      do i = 1, 5000
         read(22,*,end=60) is,in,y1,y2,xmag
         iii = iii + 1
         ydata(2*i-1) = y1
         ydata(2*i)   = y2
      enddo
      write(6,*) ' Insuff. storage for input data,',
     .           ' increase in program.  Aborting.'
      stop
 60   continue
      n2 = iii - 1
      write (6,*) ' Second list: # pts. (total):          ',n2
      write (6,*)
     .  ' ------------------------------------------------------------'
c
c---- Now carry out best match (nearest neighbor match).
      numlt1 = 0
      do i = 1, n1
c        det. closest match among pts. in 'ydata'
         mtch = 0
         dmtch = 1.e+30
         do j = 1, n2
            d = 0.0
            do k = 1, 2
               d = d + (xdata(k+(i-1)*2)-ydata(k+(j-1)*2))**2
            enddo
            if (d.lt.dmtch) then 
               dmtch = d
               mtch  = j
            endif
         enddo
         sqdmtch = sqrt(dmtch)
         if (sqdmtch.lt.1.0) numlt1 = numlt1 + 1
         write (27,200) i,mtch,sqdmtch
 200     format(i8,' in first list, --> ',i8,' in second list; dist: ',
     .                f8.5)
       enddo
       write (6,*) ' Rough measure of acceptable correspondences.'
       write (6,*) 
     .  ' Number of correspondences at distance < 1.0: ',numlt1
c
       write (6,*)
     .  ' ------------------------------------------------------------'
       end
