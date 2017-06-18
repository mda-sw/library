      dimension xdata(10000), vux(360000)
      dimension ydata(10000), vuy(360000)
c---- 'alist'    pt. list (following magnitude cut-off),
c     Note:      max # pts. (following magnitude cut-off) is 1000 at present,
c                max # pts. in total is 5000 at present.
c     Note:      in inputting arrays, the following mapping is used:
c                (i,j) --> (j+(i-1)*m)  (m = col. dimensionality)
c-----------------------------------------------------------------------------
      character  infil1*20, infil2*20, resp*1
c
      write (6,*) 
     .  ' ------------------------------------------------------------'
      write (6,*) ' POINT MATCHING, FEATURE-BASED, MAGN.-CUTOFF, V. 1.'
      write (6,*) ' Pt. matching; require input file names (ascii);'
      write (6,*) ' Free format: int.#, int.#, coordx, coordy, magn.'
      write (6,*) 
     .    ' Specify magn. limit leaving approx. 100 pts. in both lists.'
      write (6,*) 
     .  ' ------------------------------------------------------------'
c---- Get user to specify input file names, and magnitude cut-off.
c     Careful! - note expected input file format below - modify if necessary.
      write (6,21) 
      read  (5,30) infil1
      write (6,22) 
      read  (5,30) infil2
      open  (unit=21,status='old',file=infil1)
      open  (unit=22,status='old',file=infil2)
      write (6,23) 
      read  (5,*) xmglim
      write (6,*) 
     .  ' ------------------------------------------------------------'
 21   format(' Input file 1, first point set:  ',$)
 22   format(' Input file 2, second point set: ',$)
 23   format(' Input magnitude cut-off (will use pts. below this): ',$)
 30   format(a20)
c------------------------------------------------------------------------------
      ndim = 2
      ii = 1
      iii = 0
c     List x or first list: set of pts.
      do i = 1, 5000
         read(21,*,end=50) is,in,(xdata(j+(ii-1)*ndim),j=1,2),xmag
         iii = iii + 1
         if (xmag.le.xmglim) ii = ii + 1
      enddo
      write(6,*) ' Insuff. storage for input data,',
     .           ' increase in program.  Aborting.'
      stop
 50   continue
      n1 = ii - 1
      nall1 = iii
      write (6,*) ' First list: # pts. (total):           ',nall1
      write (6,*) ' First list: # pts. under mag. limit:  ',n1
c
      ii = 1
      iii = 0
c     List y or second list: set of pts.
      do i = 1, 5000
         read(22,*,end=60) is,in,(ydata(j+(ii-1)*ndim),j=1,2),xmag
         iii = iii + 1
         if (xmag.le.xmglim) ii = ii + 1
      enddo
      write (6,*) ' Insuff. storage for input data,',
     .           ' increase in program.  Aborting.'
      stop
 60   continue
c
      n2 = ii - 1
      nall2 = iii
      write (6,*) ' Second list: # pts. (total):          ',nall2
      write (6,*) ' Second list: # pts. under mag. limit: ',n2
c
c------------------------------------------------------------------------------
c---- Pause: check with the user - is he/she happy with this no. of pts.?
      write (6,54) 
      read (5,55) resp
 54   format
     .(' Are you happy with this magnitude cut-off limit? Proceed? ',$)
 55   format(a1)
      if (resp.eq.'Y'.or.resp.eq.'y') goto 58
         write (6,*) ' Okay - run again with different magnitude limit.'
         stop
 58   continue
      write (6,*) 
     .  ' ------------------------------------------------------------'
      write (6,*) 
     .' To unit 26: Pt., matched pt., confidence, rotation-angle.'
      write (6,*) 
     .' To unit 6:  Consensus rotn.; and median high-conf. transln.'
      write (6,*) 
     .  ' ------------------------------------------------------------'
c------------------------------------------------------------------------------
c---- Det. features or 360-degree profiles of all pts. in both lists.
      call profil(xdata,n1,vux)
      call profil(ydata,n2,vuy)
c
c------------------------------------------------------------------------------
c----  Prepare for matching on the basis of these features.
c----  Some fundamental constants, using by matching routine.
c      Dim. of feature space (360 one-degree steps)
       m = 360
c      Lr. and upper lts. for checking all rotns. of second pt. set: 0 and 360.
       anglo = 0.0
       anghi = 0.0
c      limit to suff. high (1, 2, ... lconf) conf. factors for det.'ing xln.
       lconf = 3
c      % limit which must be reached for consensus on rotn.
       alim = 30.0
       write (6,*) ' Input: rotation - lower and upper angle limits;'
       write (6,*) '        suff. high conf. limit value;'
       write (6,*) '        % limit for consensus on rotation.'
       write (6,*) '        (Defaults: 0  0  3  30)'
       write (6,*) '        Give four zeros to avail of these values.'
       write (6,59) 
 59    format('        Give four diff. values to override them: ',$)
       read  (5,*) v1, v2, v3, v4
       if (v1.ne.0) anglo = v1
       if (v2.ne.0) anghi = v2
       if (v3.ne.0) lconf = v3
       if (v4.ne.0) alim  = v4
       write (6,*) ' Okay, you have chosen the following:'
       write (6,*) '        Lower, upper rotation limits : ',anglo,anghi
       write (6,*) '        Conf. lim., rotn. consn. lim.: ',lconf,alim
       write (6,*) 
     .  ' ------------------------------------------------------------'
c----  Matching, using feature/profiles of all pts. in the two lists.
       call xmatch(xdata,ydata,vux,vuy,
     .             n1,n2,m,anglo,anghi,lconf,alim,
     .             xang,cxang,sxang,x,y)
c
      write (6,*) 
     .  ' ------------------------------------------------------------'
      end
c------------------------------------------------------------------------------
      subroutine profil(alist,n,vux)
      dimension alist(10000)
      dimension vecang(1000), vecd(1000), iveca(1000)
      dimension vux(360000)
c---- 'alist'    pt. list (following magnitude cut-off),
c     'vecang'   angles subtended at pt.,
c     'vecd'     distancces to other pts.,
c     'iveca'    indexes of angles subtended.
c     Note:      max # pts. (following magnitude cut-off) is 1000 at present,
c                max # pts. in total is 5000 at present.
c     Not interpolated over entire range, instead allows drop-to-zeros
c------------------------------------------------------------------------------
c
      ndim = 2
c
      vecdmax = 0.0
c---- For each pt.,
      do i = 1, n
c        Det. (i) distances to other pts. ('vecd'); and angles (in degrees).
         do j = 1, n
            if (i.eq.j) goto 100
            vecd(j) = (alist(1+(i-1)*ndim)-alist(1+(j-1)*ndim))**2+
     .                (alist(2+(i-1)*ndim)-alist(2+(j-1)*ndim))**2
            vecd(j) = sqrt(vecd(j))
            if (vecd(j).gt.vecdmax) vecdmax = vecd(j)
            if (vecd(j).ge.0.0001) then 
                    vecang(j) = 
     .               (alist(1+(j-1)*ndim)-alist(1+(i-1)*ndim))/vecd(j)
                    vecang(j) = acos(vecang(j))
                    vecang(j) = vecang(j)*180.0/3.1415926
                    if (alist(2+(j-1)*ndim).lt.alist(2+(i-1)*ndim)) 
     .               vecang(j) = 360.0 - vecang(j)
            else
                    vecang(j) = 0.0
            endif
 100        continue
         enddo
c
         do k = 1, n
            iveca(k) = k
         enddo
         vecd(i) = 0.0
         vecang(i) = 0.0
c        Convert distances to proximities; take something more that max dist
c        to prevent zero-valued proximities; normalize to 1 to near-0.
         do k = 1, n
            if (k.ne.i) vecd(k) = (1.1*vecdmax - vecd(k))/(1.1*vecdmax)
         enddo
         vecd(i) = 0.0
c        Sort the angles in-situ; sorted indexes in 'iveca'.
         call trixy(n,vecang,iveca)
c
c----    Now for each one-degree step:
         do ideg = 1, 360
c
c           Det. 'vecang' values above and below the given degree:
            do k = 2, n
               if (vecang(k).gt.float(ideg)) goto 200
            enddo
 200        ilow = k-1
            do k = 1, n
               if (vecang(k).gt.float(ideg)) goto 300
            enddo
 300        ihi = k
c
c           Get info for linear interpoln.  Allow for wrap-around (mod 360):
            if (ilow.eq.1) then
               xlo = float(ideg) + 360-vecang(n)
               vallo = vecd(iveca(n))
c              next 2 assgnts as usual
               xhi = vecang(ihi) - float(ideg)
               valhi = vecd(iveca(ihi))
            else if (ilow.eq.n) then
               xhi = 360 - float(ideg) + vecang(2)
               valhi = vecd(iveca(2))
c              next 2 assgnts as usual
               xlo = float(ideg)-vecang(ilow)
               vallo = vecd(iveca(ilow))
            else
               xlo = float(ideg)-vecang(ilow)
               vallo = vecd(iveca(ilow))
               xhi = vecang(ihi) - float(ideg)
               valhi = vecd(iveca(ihi))
            endif
c
c           Linear interpoln.
            vux(ideg+(i-1)*360) = (vallo*xhi + valhi*xlo)/(xlo+xhi)
                if (ilow.eq.1.or.ilow.eq.n) vux(indx+(i-1)*360)=0.0
         enddo
c
      enddo
c
      end
c-----------------------------------------------------------------------------
      SUBROUTINE TRIXY (N,X,Y)                                             
      DIMENSION X(N),Y(N)       
      DO 1 J=2,N      
      JJ=J-1     
      XJ=X(J)                                                                  
      DO 2 K=1,JJ   
      IF(XJ.LT.X(K)) GOTO 3                                                    
    2 CONTINUE
      GOTO 1                                                                   
    3 TR1=Y(J)
      LL=J-K
      DO  4 L=1,LL
      JML=J-L 
      X(JML+1)=X(JML)
    4 Y(JML+1)=Y(JML)
      X(K)=XJ
      Y(K)=TR1
    1 CONTINUE
      RETURN    
      END                                                                      
c-----------------------------------------------------------------------
       subroutine xmatch(xdata,ydata,xlist,ylist,
     .             n1,n2,m,anglo,anghi,lconf,alim,
     .             xang,cxang,sxang,x,y)
       dimension  xlist(360000), ylist(360000) 
       dimension  xdata(2000),   ydata(2000)
c      Orig. data: xdata, ydata (max. 1000 pts. in either list, in 2-d)
c      Feature data: xlist, ylist (max. 1000x360) in either list.
       dimension mnum(1000), valm(1000), ang(1000), 
     .           copyang(1000), mval(1000)
c----  # pts. in first pt. set: <= 1000
c      'mnum'  best matching pt. in second pt. set,
c      'valm'  dist. to best matching pt. in second pt. set,
c      'ang'   angle of rotn. for best match with this pt. in second pt. set,
c      'copyang' copy of these angles, for in-situ sorting,
c      'mval'  conf. coeff. (1 to 10) corresp. to 'valm' values.
c
       idim = 2
c
c----     For each pt. in first pt. set:
          do i = 1, n1
             xmindis = 1.e+30
c----        For each pt. in second pt. set:
             do i2 = 1, n2
c----           For each poss. angle (1-degree increments) in allowed range:
                do irot = anglo, anghi
                   d = 0.0
                   do k = 1, m
                      kk = k + irot
                      if (kk.gt.360) kk = kk - 360
c----                 Dist., given rotn.
                      d = d + (xlist(k+(i-1)*m)-ylist(kk+(i2-1)*m))**2
                   enddo
                   if (d.lt.xmindis) then
                      minx = i2
                      xmindis = d
                      angl = irot
                   endif
                enddo
              enddo
              mnum(i) = minx
              valm(i) = xmindis
              ang(i) = angl
c----         For each 'i' in 1..n1 first pt. set, have: best match in 
c----         second pt. set, 'mnum(i)', at dist. 'valm(i)', assuming angle
c----         of rotation 'ang(i)'.  
           enddo
c
c----      Det. max and min dist. values (for approx. normalizn.)
           valmmn = 1.e+30
           valmmx = -1.e+30
           do i = 1, n1
              if (valm(i).gt.valmmx) valmmx = valm(i)
              if (valm(i).lt.valmmn) valmmn = valm(i)
           enddo
c
c----      Approx. normalizn. of distances and calculn. of conf. factors.
           do i = 1, n1
c             map dist. to matched pt. (i.e. 'valmn(i)') onto 1.0..10.0 scale:
              x = 9.0*(valm(i)-valmmn)/(valmmx-valmmn)+1.0
c             ...and then onto 1..10 scale:
              mval(i) = x
              if (mval(i).eq.1) nhiconf = nhiconf+1
              write(26,740) i,mnum(i),mval(i),ang(i)
 740          format(i8,i8,i8,f10.3)
           enddo
c
c------------------------------------------------------------------------------
c----      Det. angle of rotation with most consensus
c
c----      First, take copy of set of angles (since we will sort in-situ)
           do i = 1, n1
              copyang(i) = ang(i)
           enddo
c
c----      Det. median, by first sorting in-situ
           call triage(n1,copyang)
           amed = copyang(n1/2)
c
c----      What is freq. of cases with this angle?
           mfreq = 0
           do i = 1, n1
              if (copyang(i).eq.amed) mfreq = mfreq + 1
           enddo              
c
c----      Does this indicate suff. consensus? If so, o/p rotn. If not, stop.
           percent = float(mfreq)*100.0/float(n1)
           write(6,*) ' Angle ',amed,
     .                ' degrees rotation, in % of cases:',percent
           if (percent.lt.alim) then 
              write(6,*) ' Aborting.'
              stop
           endif
c
c----      Conv. 'amed' from degrees to radians
           xang = (2.0*3.1415926)*amed/360.0
           cxang = cos(xang)
           sxang = sin(xang)
c
c------------------------------------------------------------------------------
c          Now det. mapping from 'xdata' to 'ydata', i.e. 1st and 2nd lists.
c
c          Robustify... first det. list of differences.  Use 'valm' to store.
c          Take coords. separately.
c
c          First coord.
           ngood = 0
           do i = 1, n1
              if (mval(i).le.lconf) then
                 ngood = ngood + 1
                 y1 = cxang*ydata(1+(mnum(i)-1)*idim) + 
     .                sxang*ydata(2+(mnum(i)-1)*idim)
                 valm(ngood) = xdata(1+(i-1)*idim) - y1
              endif
           enddo
           call triage(ngood,valm)
c          median translation in x
           x = valm(ngood/2)
c
c          Second coord.
           ngood = 0
           do i = 1, n1
              if (mval(i).le.lconf) then
                 ngood = ngood + 1
                 y2 = -sxang*ydata(1+(mnum(i)-1)*idim) +
     .                 cxang*ydata(2+(mnum(i)-1)*idim)
                 valm(ngood) = xdata(2+(i-1)*idim) - y2
              endif
           enddo
           call triage(ngood,valm)
c          median translation in y
           y = valm(ngood/2)
c
           write(6,710) x,y
 710       format(' Transln. vector following rotn.:', 2f10.4)
c
c----   End of routine.
        end
c
C-----------------------------------------------------------------------------
      SUBROUTINE TRIAGE (N,X)
C     IN-SITU SORTING OF VECTOR 'X'
      DIMENSION X(N)                                                
      DO 1 J=2,N                                                               
      JJ=J-1                                                                   
      XJ=X(J)                                                                  
      DO 2 K=1,JJ                                                              
      IF(XJ.LT.X(K)) GOTO 3                                                    
    2 CONTINUE                                                                 
      GOTO 1                                                                   
    3 CONTINUE                                                        
      LL=J-K                                                                   
      DO  4 L=1,LL                                                             
      JML=J-L                                                                  
      X(JML+1)=X(JML)                                                          
    4 CONTINUE
      X(K)=XJ                                                                  
    1 CONTINUE                                                                 
      RETURN                                                                   
      END                                                                      
C------------------------------------------------------------------------------
