!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   older_test.f90
! Author: Elaine Petro
!
! Created on September 12, 2019, 8:16 PM
!

! PIC code for computing charged particle trajectories
! Paulo Lozano, MIT-SPL

	program colloid_pic
	
	real A,It,dt,t,Dgz,Dgxy,r,th,deltaxy,tf
	integer np,j,k,l,jd,m,pc,pcc,s
! define the number of grid points, maximum number of particles and q/m peak distribution.
! An odd number of X-Y grids makes the best approximation since the jet charge
! is naturally accumulated close to the axis where it should be
! npeaks is the number of different q/m families
	integer, parameter::X=65,Y=65,Z=24,npmax=300000,mpeaks=500,npeaks=2,tpeaks=2000
! define all the program local variables
	real, dimension(npmax)::posx,posy,posz,velx,vely,velz,Exps,Eyps,Ezps,q_m,t_frag
        real, dimension(npmax)::colli,cx,cy
	real, dimension(X,Y,Z)::un,rho,Ex,Ey,Ez
	real, dimension(mpeaks,npeaks)::qms,dist,distc
	real, dimension(npeaks)::Amp,qms0,dwidth,V,Vin,vo,is,q,vi
        integer, dimension(npeaks)::npnew
        real, dimension(tpeaks)::t_range,t_dist,t_distc
        real t_amp,t_dwidth,tau,t_ulim
        integer, dimension(3)::f_ind
	real resint,errorp,rqm,tsec,rth,rph,vion,pi
        real posxn,posyn,poszn,zo,af
	integer startt(8),endt(8),elapst(3)
        character(8)::dte
        character(10)::tm
        character(5)::zn
        character(len=200)::memo
        
	common posx,posy,posz,velx,vely,velz,un,rho,Ex,Ey,Ez,Exps,Eyps,Ezps,q_m
	
        pi = 3.141592654
        jd = 0 ! counter for detection
        m=0 ! counter used to define zones for tecplot file
        pc=1 ! counter used to skip poisson solver
        
        write(*,'(a)',advance='no')'MEMO: '
	read(*,'(a)')memo
                
 	call date_and_time(dte,tm,zn,startt) ! start the timer	
	
	Dgz = 2.E-4  ! size of the grid elements, in meters
	Dgxy = 1.E-5 
	deltaxy = 15.E-9  ! Maximum deflection from the axial line during jet break-up
        dt = 5.E-9   ! simulation time step
        tf = 3.5E-6   !here goes the simulation time in seconds, usually 1.5 microsec
        af = 2 ! multiples of Dgz in which particles accelerate from Vi to V

!	Assuming that droplets are charged to the Rayleigh limit, their charge is given
!       by q = (q/m)*m, where (q/m) is the average specific charge, given by
!       (q/m) = 6*sqrt(g*eps0)/rho/r^1.5, with g=surface tension, rho=density,
!       eps0=8.854e-12 and r is the average droplet radius, and m is the average
!       mass, which is m = rho*4/3*pi*r^3. The charge is then q = A/(q/m), where 
!       A = 48*pi*g*eps0/rho. This is used to compute the current. MKS units.
	A = 48*pi*0.0583*8.854E-12/1133.4    ! Values for Formamide
!       A = 48*pi*0.0583*8.854E-12/1460    ! Values for Busek's liquid

!specify the q/m distribution according to the number of peaks defined above, the emission
!potential for each species
        qms0 = (/3.8e4,2.5e4/)   ! peak position (q/m)
        dwidth = (/1.,1./)     ! peak width (including potential distribution) from TOF data
        is = (/263.,263./)     ! current in nA for each species
        V = (/1618.,1618./)    ! Accelerating potentials (Busek's)        
        Vin = (/300.,300./)    ! Potentials at injection point
!	qms0 = (/310.,750./)    ! peak position (q/m)
!	dwidth = (/40.,80./)     ! peak width (including potential distribution) from TOF data
!	is = (/200.,120./)     ! current in nA for each species
!       V = (/1000.,1000./)  ! Accelerating potentials (Busek's)
!       qms0 = (/700.,70000./)    ! peak position (q/m)
!	dwidth = (/1320.,130./)     ! peak width (including potential distribution) from TOF data
!	is = (/263.,62./)     ! current in nA for each species
!	V = (/1618.,1418./)  ! Accelerating potentials (Busek's)

! adjust the amplitudes to obtain the correct current magnitude ratios
        do j=1,npeaks
          Amp(j) = is(j)/dwidth(j)
        enddo     

	open(1,file='dist.dat') ! save the distribution to a file
        do k=1,npeaks
	  do j=1,mpeaks
            llim=qms0(k)-3*dwidth(k)
            ulim=qms0(k)+3*dwidth(k)
            if(llim.le.1.)llim=1
	    qms(j,k) = llim + (ulim-llim)*(j-1)/(mpeaks-1)
            dist(j,k) = exp(-((qms(j,k)-qms0(k))/dwidth(k))**2)
	  enddo
	  call cumsum(mpeaks,Amp(k)*dist(:,k),distc(:,k))  ! computes the cummulative distribution function
	  do j=1,mpeaks
	    write(1,fmt=*)qms(j,k),Amp(k)*dist(j,k),distc(j,k)
	  enddo
        enddo
	close(1)
        
        !Make fragmentation time distribution
	open(2,file='t_dist.dat') ! save the distribution to a file
    
    !t_amp = 10           !amplitude of time gaussian curve
    !t_dwidth = 2E-10    !set spread of fragmentation times
    !tau = 4.5E-8         !set mean fragmentation time (1.5E-6 for 300nA in positive firing mode)
    !t_ulim = 9.E-8      !set upper time limit
    
    t_amp = 6           !amplitude of time gaussian curve
    t_dwidth = 2e6    !set spread of fragmentation times
    tau = 9.6e-5         !set mean fragmentation time (1.5E-6 for 300nA in positive firing mode)
    t_ulim = 1e-6      !set upper time limit
    
    
    do j=1,tpeaks
        t_range(j) = (t_ulim)*(j-1)/(tpeaks-1)        !compute time range for fragmentation
    enddo
    
    t_dist = 0.
    !t_dist = t_dist + t_amp*exp(-((t_range-tau)/t_dwidth)**2)     !compute Gaussian curves for fragmentation times
    t_dist = exp(-t_range/tau)
    
    call cumsum(tpeaks,t_dist,t_distc)  !compute cumulative time sum
    
    do j=1,tpeaks
	  write(2,fmt=*)t_range(j),t_dist(j),t_distc(j)    !print curves to 
    enddo

    close(2)	

! if only the distribution is wanted, stop execution
        if(memo.eq.'dist')then
          stop
        endif

! vi is the mean injection speed and vo is the final acceleration speed for each species
        do j=1,npeaks
	  vi(j) = sqrt(2.*qms0(j)*Vin(j))
          vo(j) = sqrt(2.*qms0(j)*V(j))
        enddo

! start with no particles in the computational domain
	np = 0
	posx = 0.
	posy = 0.
	posz = 0.
	velx = 0.
	vely = 0.
	velz = 0.
	q_m = 0.
	un = 0.
        cx = 0.
        cy = 0.
        colli = 0.
! or start from a previous run if only charge collection is wanted
        if(memo.eq.'detect')then
          open(7,file='auxnp.dat')
          read(7,fmt=*)np
          close(7)
	  open(2,file='vels.dat')
	  do j=1,np
	    read(2,fmt=*)velx(j),vely(j),velz(j)
	  enddo
	  close(2)
	  open(3,file='pos.dat')
	  do j=1,np
	    read(3,fmt=*)posx(j),posy(j),posz(j)
	  enddo
	  close(3)
	  open(4,file='q_m.dat')
	  do j=1,np
	    read(4,fmt=*)q_m(j)
	  enddo
	  close(4)
        endif

! inject a particle every tp seconds, such that it travels an axial distance Dp, which
! is given by the current carried by each species
        do j=1,npeaks
          q(j) = A/qms0(j)
          if(qms0(j).ge.20000)q(j)=1.6E-19
! if the time step for particle movement is dt, then each time step we inject npnew particles
	  npnew(j) = int(dt*is(j)*1E-9/q(j))
	  print*,'Number of particles injected per cycle =',npnew(j),' (for vo =',vo(j)/1000.,' km/s)'
	  print*,' '
        enddo
        print*,'Total number of particles injected per cycle =',sum(npnew)
	
! define the variables for each particle: position, velocity and specific charge
        t = 0.
        pcc=pc
        if(memo(1:3).eq.'rec')then
          open(1,file='superpos2.plt') ! save data for TecPlot
          write(1,fmt=*)'Title="colloid pic"'
          write(1,fmt=*)'variables="x" "y" "z"'
        endif
	do while(t<tf)  ! MAIN LOOP
           if(np+sum(npnew).gt.npmax)then
             print*,'Number of particles is greater than maximum allowed'
             goto 100     
           endif
           do k=1,npeaks
	     do j = 1,npnew(k)
	   	r = ran3(j+k)*deltaxy
		th = ran3(j+1+k)*2*pi
 		posx(np+j) = r*cos(th)
		posy(np+j) = r*sin(th)
		posz(np+j) = dt*vi(k) - (j-1)*dt*vi(k)/(npnew(k)-1)
! watch out for singularities here (last particle interfering with first particle of the next round?)
		rqm = ran3(j+2)
 		call interpolate(mpeaks,distc(:,k),qms(:,k),rqm,q_m(np+j)) ! compute q/m from the cummulative dist. fun.
! make a distinction here if the particle is a droplet or is an ion
                if(q_m(np+j).gt.20000)then
                  rth = ran3(j+3+k)
                  rth = rth*(5.*pi/180) ! taking 5 deg the initial maximum ion angle (see SOC model)
                  rph = ran3(j+4+k)
                  rph = rph*2*pi
                  vion = sqrt(2.*q_m(np+j)*Vin(k))
		  velz(np+j) = vion*cos(rth)
		  velx(np+j) = vion*sin(rth)*cos(rph) ! initial radial velocities
		  vely(np+j) = vion*sin(rth)*sin(rph)
                else
                  velz(np+j) = sqrt(2.*q_m(np+j)*Vin(k))
		  velx(np+j) = 0. ! initial radial velocities
		  vely(np+j) = 0.
                endif
             enddo
	     np = np + npnew(k)
	   enddo
	
           if(pc.eq.pcc)then
! call charge_density to interpolate the charge density at each grid point
	     call charge_density(A,np,npmax,X,Y,Z,Dgz,Dgxy,posx,posy,posz,q_m,rho)	   
! call SOR_poisson to solve Poisson's equation, returns the potential at each grid point
             call SOR_poisson(X,Y,Z,Dgz,Dgxy,rho,un)
! calculate the electric field E = -grad(phi)
	     call gradphi(X,Y,Z,Dgz,Dgxy,un,Ex,Ey,Ez)
             pc=0
           endif
           pc=pc+1
! call extrapol_field to extrapolate the electric field to each particle position
	   call extrapol_field(np,npmax,X,Y,Z,Dgz,Dgxy,posx,posy,posz,Ex,Ey,Ez,Exps,Eyps,Ezps)

! Now I can move the particles and repeat the process.
           s=2 ! s=1 is a droplet, s=2 for an ion. ONLY 2 FAMILES IN THIS ORDER!!!
	   do j = 1,np
             if(posz(j).le.Dgz*af)then
               if(q_m(j).ge.20000)s=2
               Ezps(j) = Ezps(j) + (V(s)-Vin(s))/Dgz/af
             endif
	     velx(j) = velx(j) + q_m(j)*dt*Exps(j)
	     vely(j) = vely(j) + q_m(j)*dt*Eyps(j)
	     velz(j) = velz(j) + q_m(j)*dt*Ezps(j)
	     posxn = posx(j) + velx(j)*dt
	     posyn = posy(j) + vely(j)*dt
	     poszn = posz(j) + velz(j)*dt
! If the current collector is turned on, then start counting charged particles
             if(memo.eq.'detect')then
               zo = Dgz*Z ! detector plane
               if(posz(j).lt.zo.and.poszn.ge.zo)then
                 jd = jd + 1
                 cx(jd) = ((posxn-posx(j))*zo+poszn*posx(j)-posz(j)*posxn)/(poszn-posz(j))
                 cy(jd) = ((posyn-posy(j))*zo+poszn*posy(j)-posz(j)*posyn)/(poszn-posz(j))
                 if(q_m(j).ge.20000)colli(jd) = 1.6e-19 ! is this an ion?
                 if(q_m(j).lt.20000)colli(jd) = A/q_m(j)**2 ! or a droplet? THIS IS BY MASS
               endif
             endif
             posx(j) = posxn
             posy(j) = posyn
             posz(j) = poszn
! Check whether particles exit the domain, if they do... ahem... just kill them, no hard
! feelings since these particles may not contribute too much to the total space charge.
 	     if(abs(posx(j)).gt.Dgxy*(X-1)/2.or.abs(posy(j)).gt.Dgxy*(Y-1)/2.or.&
 	    &posz(j).lt.0.or.posz(j).gt.Dgz*Z)then
                    do k = j,np
                        posx(k) = posx(k+1)
                        posy(k) = posy(k+1)
                        posz(k) = posz(k+1)
                        velx(k) = velx(k+1)
                        vely(k) = vely(k+1)
                        velz(k) = velz(k+1)
                        q_m(k) = q_m(k+1)
                        Exps(k) = Exps(k+1)
                        Eyps(k) = Eyps(k+1)
                        Ezps(k) = Ezps(k+1)
                    enddo
                    np = np - 1
 	     endif
             
!             !This part breaks dimers into monomers if their time has expired :)
!             if((t.ge.tfrag(j)).and.(tfrag(j).gt.0)) then
!                    !pass dimer information to monomer 
!                    posx(f_ind(1)+1) = posx(j);
!                    posy(f_ind(1)+1) = posy(j);
!                    posz(f_ind(1)+1) = posz(j); 
!
!                    velx(f_ind(1)+1) = velx(j,2); 
!                    vely(f_ind(1)+1) = vely(j,2);
!                    velz(f_ind(1)+1,1) = velz(j,2); 
!            
!                    !assign monomer charge and update counter
!                    q_m(f_ind(1)+1,1) = qms0(1);
!                    f_ind(1) = f_ind(1) + 1;    !update monomer index
!            
!                    !pass dimer information to neutral
!                    posx(f_ind(3)+1,3) = posx(j,2);
!                    posy(f_ind(3)+1,3) = posy(j,2);
!                    posz(f_ind(3)+1,3) = posz(j,2); 
!
!                    velx(f_ind(3)+1,3) = velx(j,2); 
!                    vely(f_ind(3)+1,3) = vely(j,2);
!                    velz(f_ind(3)+1,3) = velz(j,2);         
!            
!                    !assign neutral charge and update counter
!                    q_m(f_ind(3)+1,3) = 0;
!                    f_ind(3) = f_ind(3) + 1;    !update neutral index
!                  
!                    !destroy dimer and replace with the last dimer
!                    posx(j,2) = posx(f_ind(2),2);
!                    posy(j,2) = posy(f_ind(2),2);
!                    posz(j,2) = posz(f_ind(2),2); 
!
!                    velx(j,2) = velx(f_ind(2),2); 
!                    vely(j,2) = vely(f_ind(2),2);
!                    velz(j,2) = velz(f_ind(2),2);
!            
!                    q_m(j,2) = q_m(f_ind(2),2);
!                    tfrag(j,2) = tfrag(f_ind(2),2);
!            
!                    !Zero out final dimer that was transferred
!                    posx(f_ind(2),2) = 0;
!                    posy(f_ind(2),2) = 0;
!                    posz(f_ind(2),2) = 0; 
!
!                    velx(f_ind(2),2) = 0; 
!                    vely(f_ind(2),2) = 0;
!                    velz(f_ind(2),2) = 0;
!            
!                    q_m(f_ind(2),2) = 0;
!                    tfrag(f_ind(2),2) = 0;
!            
!                    f_ind(2) = f_ind(2) - 1;    !update dimer index
!                endif
           
 	   enddo
        
 	   t = t + dt
 	   print*,'time = ',t
           print*,'np = ',np
           
           if(memo(1:3).eq.'rec')then
             if(t.ge.5*(m+1)*dt-dt/10.and.t.le.5*(m+1)*dt+dt/10)then
               m=m+1
               write(1,fmt=*)'ZONE T="TIME=',t*1E6,'" I=',np,' F=POINT'
               do j=1,np
	         write(1,fmt=*)posx(j),posy(j),posz(j)
	       enddo
             endif
           endif
	enddo
	close(1)
        
    100	print*,'np = ',np

        if(memo.eq.'detect')then
	  open(2,file='coll.dat')
	  do j=1,jd
	    write(2,fmt=*)cx(j),cy(j),colli(j)
	  enddo
          stop
        endif

	open(2,file='vels.dat')
	do j=1,np
	  write(2,fmt=*)velx(j),vely(j),velz(j)
	enddo
	close(2)
	open(3,file='pos.dat')
	do j=1,np
	  write(3,fmt=*)posx(j),posy(j),posz(j)
	enddo
	close(3)
	open(4,file='q_m.dat')
	do j=1,np
	  write(4,fmt=*)q_m(j)
	enddo
	close(4)
	
	call charge_density(A,np,npmax,X,Y,Z,Dgz,Dgxy,posx,posy,posz,q_m,rho)
	open(5,file='rho.dat')
	write(5,fmt=*)X
	write(5,fmt=*)Y
	write(5,fmt=*)Z
	do j=1,X
	  do k=1,Y
	    do l=1,Z
	      write(5,fmt=*)rho(j,k,l)
	    enddo
	  enddo
	enddo
	close(5)
	
 	call date_and_time(dte,tm,zn,endt) ! stop the timer
	tsec = 3600*(endt(5)-startt(5)) + 60*(endt(6)-startt(6)) + (endt(7)-startt(7))
	elapst(1) = int(tsec/3600.)
	elapst(2) = int((tsec/3600.-elapst(1))*60)
	elapst(3) = nint(((tsec/3600.-elapst(1))*60-elapst(2))*60)
	
!       The total current is given by It = dQ/dt, and dt = dZ/vo, dZ = Dgz*Z and dQ = np*q
!	A is given above, and resint is the average (q/m)
!	It = vo*A*np/resint/Dgz/Z
	
	open(6,file='stats.dat')
	write(6,fmt=*)'Results from running the Colloid PIC code'
	write(6,fmt=*)' '
        write(6,fmt=*)'MEMO: ',memo
	write(6,fmt=*)'Started:  ',startt(2),'/',startt(3),'/',startt(1),'  at  ',startt(5),':',startt(6),':',startt(7)
	write(6,fmt=*)'Finished: ',endt(2),'/',endt(3),'/',endt(1),'  at  ',endt(5),':',endt(6),':',endt(7)
	write(6,fmt=*)' '
	if(elapst(1).ge.0)then
	 write(6,fmt=*)'Elapsed time for calculation:',elapst(1),'hrs,',elapst(2),'min,',elapst(3),'sec'
	 write(6,fmt=*)' '
	endif
	write(6,fmt=*)'Accelerating voltage =',V,' V'
	write(6,fmt=*)'Dgz =',Dgz,' m'
	write(6,fmt=*)'Dgxy =',Dgxy,' m'
	write(6,fmt=*)'Domain size in Z =',Dgz*Z*1.E3,' mm'
	write(6,fmt=*)'Domain size in XY =',Dgxy*X*1.E3,' mm'
	write(6,fmt=*)'Initial deflection =',deltaxy,' m'
	write(6,fmt=*)'Initial axial mean velocity =',vo,' m/s'
!	write(6,fmt=*)'Distance between droplets =',Dp*1.E9,' nm'
	write(6,fmt=*)'Time step =',dt,' sec'
	write(6,fmt=*)'Injected particles per cycle =',npnew
	write(6,fmt=*)'Simulation time =',t,' sec'
	write(6,fmt=*)'Total number of particles =',np
!	write(6,fmt=*)'Current =',It*1.E9,' nA'
	close(6)

        open(7,file='auxnp.dat')
        write(7,fmt=*)np
        close(7)

	end program colloid_pic
	
! 	**************************  Subroutines  ************************

! 	Charge density calculation for a 3D rectangular grid
	subroutine charge_density(A,np,npmax,X,Y,Z,Dgz,Dgxy,posx,posy,posz,q_m,rho)
	
	integer X,Y,Z,np,npmax,n,Xg,Yg,Zg
	real, dimension(npmax)::posx,posy,posz,q_m
	real, dimension(X,Y,Z)::rho
	real q,A,Dgz,Dgxy,Vol,dx,dy,dz,dxm,dym,dzm
	real V1,V2,V3,V4,V5,V6,V7,V8

	Vol = Dgz*Dgxy**2. ! Vol is the cell volume

! 	First, consider charge density = 0 for all grid points
	rho = 0.

! 	Now check the position of each particle and assign its interpolated
! 	charge density value to its neighboring grid points

	do n = 1,np
! 	determine the grid in which this particle is located
		Xg = int(posx(n)/Dgxy + (X-1)/2.) + 1
		Yg = int(posy(n)/Dgxy + (Y-1)/2.) + 1
		Zg = int(posz(n)/Dgz) + 1
	
		if(posx(n).ge.(X-1)*Dgxy/2.)Xg = X-1
		if(posy(n).ge.(Y-1)*Dgxy/2.)Yg = Y-1
		if(posz(n).ge.(Z-1)*Dgz)Zg = Z-1
	
! 		if(Xg.lt.1.or.Xg.gt.X-1)then
! 		   print*,'Xg is out: ',Xg
! 		   print*,posx(n)
! 		   stop
! 		endif
! 		if(Yg.lt.1.or.Yg.gt.Y-1)then
! 		   print*,'Yg is out: ',Yg
! 		   stop
! 		endif
! 		if(Zg.lt.1.or.Zg.gt.Z-1)then
! 		   print*,'Zg is out: ',Zg
! 		   stop
! 		endif

		dx = posx(n) + Dgxy*(X/2. - Xg) + Dgxy/2.
		dxm = Dgxy-dx
		dy = posy(n) + Dgxy*(Y/2. - Yg) + Dgxy/2.
		dym = Dgxy-dy
		dz = posz(n) + Dgz*(1-Zg)
		dzm = Dgz-dz
! 		Use Rayleigh criterion to separate the charge from the mass
                if(q_m(n).gt.20000)then
                  q = 1.6e-19
                else
                  q = A/q_m(n)
                endif

		V1 = dxm*dym*dzm
		rho(Xg,Yg,Zg) = rho(Xg,Yg,Zg) + q*V1/Vol**2.
	
		V2 = dxm*dy*dzm
		rho(Xg,Yg+1,Zg) = rho(Xg,Yg+1,Zg) + q*V2/Vol**2.
	
		V3 = dx*dym*dzm
		rho(Xg+1,Yg,Zg) = rho(Xg+1,Yg,Zg) + q*V3/Vol**2.
	
		V4 = dx*dy*dzm
		rho(Xg+1,Yg+1,Zg) = rho(Xg+1,Yg+1,Zg) + q*V4/Vol**2.
	
		V5 = dxm*dym*dz
		rho(Xg,Yg,Zg+1) = rho(Xg,Yg,Zg+1) + q*V5/Vol**2.
	
		V6 = dxm*dy*dz
		rho(Xg,Yg+1,Zg+1) = rho(Xg,Yg+1,Zg+1) + q*V6/Vol**2.
	
		V7 = dx*dym*dz
		rho(Xg+1,Yg,Zg+1) = rho(Xg+1,Yg,Zg+1) + q*V7/Vol**2.
	
		V8 = dx*dy*dz
		rho(Xg+1,Yg+1,Zg+1) = rho(Xg+1,Yg+1,Zg+1) + q*V8/Vol**2.


	enddo
	return
	end subroutine
	
! 	SOR for solving Poisson's equation on a 3D rectangular grid
	subroutine SOR_poisson(X,Y,Z,Dgz,Dgxy,rho,un)
	
	integer X,Y,Z,j,k,l,p
	real, dimension(X,Y,Z)::uo,un,rho
	real epsi,Dgz,Dgxy,w,coef,err
	epsi = 8.854E-12

! 	Define relaxation parameter
	w = 1.775

! 	Define arrays
	uo = un
	err = 0.0
	
! 	Dirichlet boundary conditions
 	un(1,:,:) = 0.0
   	un(X,:,:) = 0.0
 	un(:,1,:) = 0.0
 	un(:,Y,:) = 0.0
	un(:,:,1) = 0.0
	un(:,:,Z) = 0.0	
	
! 	Start the iterations
	p=1
	do
 
! 	Neumann boundary Conditions
! 	   un(:,:,Z) = un(:,:,Z-1)
! 	   un(1,:,:) = un(2,:,:)
!   	   un(X,:,:) = un(X-1,:,:)
! 	   un(:,1,:) = un(:,2,:)
! 	   un(:,Y,:) = un(:,Y-1,:)    
   
! 	Solve Poisson's equation
   	   do j=2,X-1
	     do k=2,Y-1
	       do l=2,Z-1
	         coef = 2./Dgz**2. + 4./Dgxy**2.
		 un(j,k,l) = (w/coef)*((1/Dgxy**2.)*(uo(j+1,k,l) + &
		&un(j-1,k,l) + uo(j,k+1,l) + un(j,k-1,l)) + (1/Dgz**2.)*(uo(j,k,l+1) + &
		&un(j,k,l-1))) + (1.-w)*uo(j,k,l) + w*rho(j,k,l)/epsi/coef	
 	        err = err + (un(j,k,l) - uo(j,k,l))**2.
	       enddo
	     enddo
   	   enddo
	   
	   uo = un
   	   err = sqrt(err)
	   if(err.lt.0.01)exit
!	   if(err.lt.0.001)exit
	   err=0.
	   p=p+1
	   
	enddo

        print*,'Converged with ',p,' iterations and u max is',maxval(un)	
 !	print*,'Error =',err,' with ',p,' iterations'
 	print*,' '
	
	return
	end subroutine
	
! 	Compute the the electric field from the potential distribution
	subroutine gradphi(X,Y,Z,Dgz,Dgxy,un,Ex,Ey,Ez)
	
	integer X,Y,Z,j,k,l,p
	real, dimension(X,Y,Z)::un,Ex,Ey,Ez
	real, allocatable::fun(:),dfun(:)
	real Dgz,Dgxy

	Ex = 0.
	Ey = 0.
	Ez = 0.

	allocate(fun(X),dfun(X))
	fun = 0.
	dfun = 0.
	do k=1,Y
	  do l=1,Z
	    fun = un(:,k,l)
	    dfun(1) = (3.*(fun(2)-fun(1))-(fun(3)-fun(2)))/(2.*Dgxy)
	    do p = 2,X-1
  	      dfun(p) = (fun(p+1)-fun(p-1))/(2*Dgxy)
	    enddo
	    dfun(X) = (3*(fun(X)-fun(X-1))-(fun(X-1)-fun(X-2)))/(2.*Dgxy)
	    Ex(:,k,l) = dfun
	  enddo
	enddo
	deallocate(fun,dfun)
	Ex = -Ex
	
	allocate(fun(Y),dfun(Y))
	fun = 0.
	dfun = 0.
	do j=1,X
	  do l=1,Z
	    fun = un(j,:,l)
	    dfun(1) = (3.*(fun(2)-fun(1))-(fun(3)-fun(2)))/(2.*Dgxy)
	    do p = 2,Y-1
  	      dfun(p) = (fun(p+1)-fun(p-1))/(2*Dgxy)
	    enddo
	    dfun(Y) = (3.*(fun(Y)-fun(Y-1))-(fun(Y-1)-fun(Y-2)))/(2.*Dgxy)
	    Ey(j,:,l) = dfun
	  enddo
	enddo
	deallocate(fun,dfun)
	Ey = -Ey

	allocate(fun(Z),dfun(Z))
	fun = 0.
	dfun = 0.
	do j=1,X
	  do k=1,Y
	    fun = un(j,k,:)
	    dfun(1) = (3.*(fun(2)-fun(1))-(fun(3)-fun(2)))/(2.*Dgz)
	    do p = 2,Z-1
  	      dfun(p) = (fun(p+1)-fun(p-1))/(2*Dgz)
	    enddo
	    dfun(Z) = (3.*(fun(Z)-fun(Z-1))-(fun(Z-1)-fun(Z-2)))/(2.*Dgz)
	    Ez(j,k,:) = dfun
	  enddo
	enddo
	Ez = -Ez
	
	return
	end subroutine
	
! 	Electric field extrapolation from a 3D rectangular grid
	subroutine extrapol_field(np,npmax,X,Y,Z,Dgz,Dgxy,posx,posy,posz,Ex,Ey,Ez,Exps,Eyps,Ezps)
	
	integer X,Y,Z,np,npmax,n,Xg,Yg,Zg
	real, dimension(npmax)::posx,posy,posz,Exps,Eyps,Ezps
	real, dimension(X,Y,Z)::Ex,Ey,Ez
	real Dgz,Dgxy,Vol,dx,dy,dz,dxm,dym,dzm
	real V1,V2,V3,V4,V5,V6,V7,V8

	Vol = Dgz*Dgxy**2. ! Vol is the cell volume
	
	Exps(1:np) = 0.
	Eyps(1:np) = 0.
	Ezps(1:np) = 0.
	
	do n = 1,np
! 	determine the grid in which this particle is located
	  Xg = int(posx(n)/Dgxy + (X-1)/2.) + 1
	  Yg = int(posy(n)/Dgxy + (Y-1)/2.) + 1
	  Zg = int(posz(n)/Dgz) + 1
	
	  if(posx(n).ge.(X-1)*Dgxy/2.)Xg = X-1
	  if(posy(n).ge.(Y-1)*Dgxy/2.)Yg = Y-1
	  if(posz(n).ge.(Z-1)*Dgz)Zg = Z-1
	
	  dx = posx(n) + Dgxy*(X/2. - Xg) + Dgxy/2.
	  dxm = Dgxy-dx
	  dy = posy(n) + Dgxy*(Y/2. - Yg) + Dgxy/2.
	  dym = Dgxy-dy
	  dz = posz(n) + Dgz*(1-Zg)
	  dzm = Dgz-dz
		
	  V1 = dxm*dym*dzm
	  V2 = dxm*dy*dzm
	  V3 = dx*dym*dzm
	  V4 = dx*dy*dzm
	  V5 = dxm*dym*dz
	  V6 = dxm*dy*dz
	  V7 = dx*dym*dz
	  V8 = dx*dy*dz
	  
	  Exps(n) = (Ex(Xg,Yg,Zg)*V1 + Ex(Xg,Yg+1,Zg)*V2 + Ex(Xg+1,Yg,Zg)*V3 +&
	 &Ex(Xg+1,Yg+1,Zg)*V4 + Ex(Xg,Yg,Zg+1)*V5 + Ex(Xg,Yg+1,Zg+1)*V6 +&
	 &Ex(Xg+1,Yg,Zg+1)*V7 + Ex(Xg+1,Yg+1,Zg+1)*V8)/Vol
	
	  Eyps(n) = (Ey(Xg,Yg,Zg)*V1 + Ey(Xg,Yg+1,Zg)*V2 + Ey(Xg+1,Yg,Zg)*V3 +&
	 &Ey(Xg+1,Yg+1,Zg)*V4 + Ey(Xg,Yg,Zg+1)*V5 + Ey(Xg,Yg+1,Zg+1)*V6 +&
	 &Ey(Xg+1,Yg,Zg+1)*V7 + Ey(Xg+1,Yg+1,Zg+1)*V8)/Vol
	
	  Ezps(n) = (Ez(Xg,Yg,Zg)*V1 + Ez(Xg,Yg+1,Zg)*V2 + Ez(Xg+1,Yg,Zg)*V3 +&
	 &Ez(Xg+1,Yg+1,Zg)*V4 + Ez(Xg,Yg,Zg+1)*V5 + Ez(Xg,Yg+1,Zg+1)*V6 +&
	 &Ez(Xg+1,Yg,Zg+1)*V7 + Ez(Xg+1,Yg+1,Zg+1)*V8)/Vol

	enddo
	
	return
	end subroutine

! this is a standard routine from numerical recipies that gives random numbers
	function ran3(idum)
	
	integer idum
	integer mbig,mseed,mz
	real ran3,fac
	parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
	integer i,iff,ii,inext,inextp,k
	integer mj,mk,ma(55)
	save iff,inext,inextp,ma
	data iff /0/
	if(idum.lt.0.or.iff.eq.0)then
	  iff=1
	  mj=mseed-iabs(idum)
	  mj=mod(mj,mbig)
	  ma(55)=mj
	  mk=1
	  do i=1,54
	    ii=mod(21*i,55)
	    ma(ii)=mk
	    mk=mj-mk
	    if(mk.lt.mz)mk=mk+mbig
	    mj=ma(ii)
	  enddo
	  do k=1,4
	    do i=1,55
	      ma(i)=ma(i)-ma(1+mod(i+30,55))
	      if(ma(i).lt.mz)ma(i)=ma(i)+mbig
	    enddo
	  enddo
	  inext=0
	  inextp=31
	  idum=1
	endif
	inext=inext+1
	if(inext.eq.56)inext=1
	inextp=inextp+1
	if(inextp.eq.56)inextp=1
	mj=ma(inext)-ma(inextp)
	if(mj.lt.mz)mj=mj+mbig
	ma(inext)=mj
	ran3=mj*fac
	return
	end
	
! Obtain the cummulative distribution function
	subroutine cumsum(mpeaks,dist,distc)
	
	integer j,k,mpeaks
	real, dimension(mpeaks)::dist,distc
        
	distc = 0.
	distc(1)=dist(1)
	do j=2,mpeaks
	  do k = 1,j
	    distc(j) = distc(j-1) + dist(k)
	  enddo
	enddo

	distc = distc/maxval(distc)

	return
	end subroutine	
	
! Integrate a function dist with respect to qms
	subroutine integ(mpeaks,qms,dist,resint)
	
	integer j,mpeaks
	real, dimension(mpeaks)::dist,qms
	real resint,hw
	
	hw = (maxval(qms)-minval(qms))/(mpeaks-1)
	resint = 0;
! 	resint = (hw/2.)*(dist(1)+dist(mpeaks))
	do j=1,mpeaks-1
 	  resint = resint + hw*(dist(j)+dist(j+1))/2.
! 	  resint = resint + hw*dist(j)	  
	enddo
	
	return
	end subroutine
	
! Interpolation routine
	subroutine interpolate(mpeaks,distc,qms,rqm,qmg)
	
	integer mpeaks,j,jl,ju,jm
	real, dimension(mpeaks)::distc,qms
	real rqm,qmg,f
	
	jl=0
	ju=mpeaks+1
    10  if(ju-jl.gt.1)then
    	  jm=(ju+jl)/2
	  if((distc(mpeaks).ge.distc(1)).eqv.(rqm.ge.distc(jm)))then
	    jl=jm
	  else
	    ju=jm
	  endif
	goto 10
	endif
	if(rqm.eq.distc(1))then
	  j=1
	else if(rqm.eq.distc(mpeaks))then
	  j=mpeaks-1
	else
	  j=jl
	endif
	
	if(distc(j+1)-distc(j).eq.0.)then
	  f=0.
	else
	  f=(rqm-distc(j))/(distc(j+1)-distc(j))
	endif
	qmg = qms(j)+f*(qms(j+1)-qms(j))
	
	return
	end subroutine


