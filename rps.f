c RPS model in May-Leonard formalism
c incorporating predation,reproduction,natural death
c and without incorporating mobility (nonspatial RPS model)
c There are vacant sites

c The program is averaged over different set of
c random number generation

c Here each time step is defined by N**2 number of
c monte carlo steps

c The program prints the spatial distribution of the
c lattice at multiple time steps

c date 03/04/2020

        implicit double precision (a-h,o-z)
        parameter(N=500,mcstep=N*N,nav=10,itmax=8.0e4)
        integer lattice(N,N),ip(N),im(N),num(0:3),step,count
        double precision d(3),r(3),p(3),frac(0:3),rho(0:3,0:itmax)


	step=1

        itp1 = itmax/4
        itp2 = 2*itmax/4
        itp3 = 3*itmax/4
        itp4 = itmax


        open(10,file='density-time-astable')
        open(20,file='fixed-point-astable')


        open(456, file='lattice0')
        open(457,file='lattice1')
        open(458, file='lattice2')
        open(459, file='lattice3')
        open(460, file='lattice4')
        
        open(461,file='latticeFinal.txt')

	open(333, file='initials')
	open(334, file='initials2')

	open(11, file='final-density')
	
	open(2545, file='extictionprob')






        do ii = 0,3
           do jj = 0,itmax
              rho(ii,jj)=0.d0
           enddo
        enddo

*********************************************************
*      Assigning probability
*********************************************************

       diffusion = 4.0e-3
       eps = 2*N*N*diffusion

       death = 2.0d0
       rep = 7.0d0
       pred = 7.0d0

       !Normalization

       tot = eps+death+rep+pred
       hop_norm = eps/tot
       death_norm = death/tot
       rep_norm = rep/tot
       pred_norm = pred/tot


       !Unnormalization!
      ! hop_norm = 0.7d0
      ! death_norm = 0.0d0
      ! rep_norm = 0.03d0
      ! pred_norm = 0.03d0

       d(1) = death_norm
       d(2) = death_norm
       d(3) = death_norm

       r(1) = rep_norm
       r(2) = rep_norm
       r(3) = rep_norm

       p(1) = pred_norm
       p(2) = pred_norm
       p(3) = pred_norm



       hop = hop_norm

       !write(*,*)hop_norm, death_norm, pred_norm, rep_norm




*********************************************************
*      Calculation of fixed points
*********************************************************
       q = (1.d0 + (d(1)/p(3)) + (d(2)/p(1)) + (d(3)/p(2)))
     ^     /(1.d0 + (r(1)/p(3)) + (r(2)/p(1)) + (r(3)/p(2)))

       afix = ((q*r(2))-(d(2)))/p(1)
       bfix = ((q*r(3))-(d(3)))/p(2)
       cfix = ((q*r(1))-(d(1)))/p(3)



       write(20,*) '# d(1), r(1), p(1)',d(1),r(1),p(1)
       write(20,*) '# d(2), r(2), p(2)',d(2),r(2),p(2)
       write(20,*) '# d(3), r(3), p(3)',d(3),r(3),p(3)
       write(20,*) '# afix, bfix, cfix',afix,bfix,cfix



*********************************************************
* Averaging with different iseed
*********************************************************
        m = 548512
        
        count=0
        
        do iim = 1,nav
          m = m + 1
          iseed = m


*********************************************************
*      INITIALIZING LATTICE
*********************************************************
        do ii = 0,3
           num(ii) = 0
        enddo

        do i=1,N
            do j=1,N
               x=random(iseed)
               if(x.lt.0.25d0) then
                 lattice(i,j)=0
                 num(0)=num(0)+1
               elseif(x.gt.0.25d0.and.x.lt.0.5d0) then
                 lattice(i,j)=1
                 num(1)=num(1)+1
               elseif(x.gt.0.5d0.and.x.lt.0.75d0) then
                 lattice(i,j)=2
                 num(2)=num(2)+1
               else
                 lattice(i,j)=3
                 num(3)=num(3)+1
               endif
            enddo
         enddo

         do ii = 1,N
            do jj = 1,N
		write(456,*)lattice(ii,jj)
            enddo
         enddo
*********************************************************
*      Periodic Boundary Condition
*********************************************************
        do i=1,N
           ip(i)=i+1
           im(i)=i-1
        enddo
        ip(N)=1
        im(1)=N
*********************************************************
*      STARTING MONTE CARLO
*********************************************************

        do ii = 0,3
           frac(ii) = dfloat(num(ii))/dfloat(N**2)
           rho(ii,0) = rho(ii,0) + frac(ii)
        enddo



        do it = 1, itmax, step
           do k=1,mcstep
666           xn=random(iseed)
              xn=(xn*(N)) + 1.d0
              ino=int(xn)

              xn=random(iseed)
              xn=(xn*(N)) + 1.d0
              jno=int(xn)

              if(ino.gt.N.or.jno.gt.N) goto 666

                 ispini=lattice(ino,jno)

                 xr= random(iseed)
                 if(random(iseed).lt.0.5) then
                   if(xr.lt.0.5)then
                     inext=ip(ino)
                   else
                     inext=im(ino)
                   endif
                   jnext=jno
                 else
                   if(xr.lt.0.5) then
                     jnext=ip(jno)
                   else
                     jnext=im(jno)
                   endif
                   inext=ino
                 endif

                 if(inext.gt.N.or.jnext.gt.N) goto 666


                 ispinj=lattice(inext,jnext)

c  The code for interactions begins here.


                aa = hop+pred_norm
       	        ab = aa + rep_norm

                rand = random(iseed)

                if(rand.le.hop) goto 121
                if(rand.gt.hop.and.rand.le.aa) goto 122
                if(rand.gt.aa.and.rand.le.ab) goto 123
                if(rand.gt.ab.and.rand.le.1.d0) goto 124


121             if(ispini.ne.ispinj) then
                  tempo = ispini
                  ispini=ispinj
                  ispinj=tempo
                end if
                goto 111


122             if(ispini.eq.1.and.ispinj.eq.2) ispinj=0
                if(ispini.eq.1.and.ispinj.eq.3) ispini=0
                if(ispini.eq.2.and.ispinj.eq.1) ispini=0
                if(ispini.eq.2.and.ispinj.eq.3) ispinj=0
                if(ispini.eq.3.and.ispinj.eq.1) ispinj=0
                if(ispini.eq.3.and.ispinj.eq.2) ispini=0
                goto 111

123             if(ispini.ne.0.and.ispinj.eq.0) ispinj=ispini
                goto 111

124             if(ispini.ne.0) ispini=0


111             lattice(ino,jno)=ispini
                lattice(inext,jnext)=ispinj

           enddo ! ends mcstep

           do ii = 0,3
              num(ii) = 0
           enddo
           do ii=1,N
              do jj=1,N
                 inum = lattice(ii,jj)
                 num(inum) = num(inum) + 1
              enddo
           enddo

           do ii = 0,3
              frac(ii) = dfloat(num(ii))/dfloat(N**2)
              rho(ii,it) = rho(ii,it) + frac(ii)
           enddo

*********************************************************
* Printing the lattice at selected time-steps
*********************************************************
           if(it.eq.itp1) then
             do ii = 1,N
                do jj = 1,N
                   write(457,*)lattice(ii,jj)
                enddo
             enddo
           endif

           if(it.eq.itp2) then
             do ii = 1,N
                do jj = 1,N
                   write(458,*)lattice(ii,jj)
                enddo
             enddo
           endif

           if(it.eq.itp3) then
             do ii = 1,N
                do jj = 1,N
	          write(459,*)lattice(ii,jj)
                enddo
             enddo
           endif

	if(it.eq.itp4) then
             do ii = 1,N
                do jj = 1,N
	          write(460,*)lattice(ii,jj)
                enddo
             enddo
           endif

        !Writing the final lattice in a different way.
        !This is for the correlation code.
        
         if(it.eq.itp4) then
             do ii = 1,N
                do jj = 1,N
	          write(461,*)ii,jj,lattice(ii,jj)
                enddo
             enddo
           endif

*********************************************************

        enddo ! ends it
        
        if(frac(1).lt.0.001)then
         count=count+1
         
        else if(frac(2).lt.0.001)then
         count=count+1

        else if(frac(3).lt.0.001)then
         count=count+1
        end if

        write(*,*) 'completed step',iim

        enddo ! ends iim
        
        write(2545,*)diffusion,count,nav,count/dfloat(nav)
        

        write(10,*) '#time 	fracN1		fracN2		fracN3'
        do iit = 1,itmax,step
c           if(mod(iit,100000).eq.0)
           write(10,*) iit,rho(1,iit)/nav,rho(2,iit)/nav,rho(3,iit)/nav
        enddo

        write(333,*)'#	N	itp1	itp2	itp3	d	r	p	hop'
        do i=1,3
        	write(333,*)N,itp1,itp2,itp3,itp4,d(i),r(i),p(i),hop
        end do
        
        write(334,*)'#  N   D  death  rep  pred'
        write(334,*)N,diffusion,death,rep,pred


        write(11,*)iseed, rho(1,itmax),rho(2,itmax),rho(3,itmax)




        close(10)
        close(20)
        close(456)
        close(457)
        close(458)
        close(459)
        close(460)
        close(461)
        close(11)
        close(333)
        close(334)
        close(2545)






      end

*********************************************************
c      Random number generating function
*********************************************************
        double precision function random(iseed)
        iseed=iseed*1566083941
        if(iseed.lt.0)iseed=iseed+2147483647+1
        iseed=iseed*1566083941
        if(iseed.lt.0)iseed=iseed+2147483647+1
        random=iseed*4.6566128752458D-10
        return
        end
*********************************************************

