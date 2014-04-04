Module Initial_Conditions
	Use ProblemSize
	Use Fields
	Use Parallel_Framework
	Use Fourier_Transform
	Use Legendre_Transforms, Only : Legendre_Transform 
	Use SendReceive
	Use Chebyshev_Polynomials, Only : Cheby_To_Spectral
	Use Checkpointing, Only : read_checkpoint
	Use Controls
	Implicit None
	Integer :: init_type = 1
	Integer :: magnetic_init_type = 1
	Integer :: init_tag = 8989
	Integer :: restart_iter = -1
	Real*8 :: pi = 3.1415926535897932384626433832795028841972d+0
	Real*8 :: temp_amp = 1.0d0, temp_w = 0.3d0
	Logical :: custom_t
	Character*120 :: custom_t_file
	Namelist /Initial_Conditions_Namelist/ init_type, temp_amp, temp_w, custom_t, custom_t_file, restart_iter, &
			magnetic_init_type
Contains
	
	Subroutine Initialize_Fields()
		Implicit None
		Logical :: dbtrans, dbconfig
		! When coming out of this routine, the RHS of the equation set should contain the field values.
		! This setup is consistent with the program having just completed a time step

		! wsp%p1b should contain the Adams-Bashforth terms from the previous (not current)
		! values of WPST.  This means that wsp%p1b should be zero if this init is from scratch.
		! If this init is from restart, wsp%p1b should contain the adams bashforth terms output
		! as part of the checkpoint.

		! Check control variables to see if we need want static or buffers
		dbtrans = .not. static_transpose
		dbconfig = .not. static_config


		Call wsp%init(field_count = wsfcount, config = 'p1b', &
			dynamic_transpose =dbtrans, dynamic_config = dbconfig, &
			piggyback = test_reduce, padding = pad_alltoall)		
		Call wsp%construct('p1b')	! We will always start in p1b - should do wsp%set_config('p1b')
		wsp%p1b(:,:,:,:) = 0.0d0	! All fields are zero initially

		! Allocate the Equation Set RHS
		! Set it to zero initially	
		! The equation set RHS's stays allocated throughout - it is effectively how we save the AB terms.
		Call Allocate_RHS(zero_rhs=.true.)

		If (init_type .eq. -1) Then
			Call restart_from_checkpoint(restart_iter)
		Endif
	
		If (init_type .eq. 1) Then
			call benchmark_init_hydro()
		Endif
		If (init_Type .eq. 2) Then
			call diffusion_init_hydro()
		Endif
		If (init_Type .eq. 3) Then
			call temperature_blob()
		Endif
        If (init_Type .eq. 4) Then
            call convective_init()
        Endif
        If (init_Type .eq. 5) Then
            call random_init()
        Endif

		If (magnetism) Then
			If (magnetic_init_type .eq. 1) Then
				call benchmark_insulating_init()
			Endif
		Endif
		! Fields are now initialized and loaded into the RHS. 
		! We are ready to enter the main loop

	End Subroutine Initialize_Fields

	Subroutine Restart_From_Checkpoint(iteration)
		Implicit None
		Integer, Intent(In) :: iteration
		type(SphericalBuffer) :: tempfield
		Integer :: fcount(3,2)
		fcount(:,:) = 4

		Call tempfield%init(field_count = fcount, config = 'p1a')
		Call tempfield%construct('p1a')

		wsp%p1b(:,:,:,:) = 0.0d0
		tempfield%p1a(:,:,:,:) = 0.0d0

		Call Read_Checkpoint(tempfield%p1a,wsp%p1b,iteration)

		Call Set_All_RHS(tempfield%p1a)
		Call tempfield%deconstruct('p1a')

	End Subroutine Restart_From_Checkpoint

	Subroutine Benchmark_Init_Hydro()
		Implicit None
		Real*8, Allocatable :: rfunc1(:), rfunc2(:)
		Real*8 :: x
		Integer :: i, r, l, m, mp
		Integer :: fcount(3,2)
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 1

		Allocate(rfunc1(my_r%min: my_r%max))
		Allocate(rfunc2(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = 2.0d0*radius(r)-r_inner-r_outer
			rfunc1(r) = 0.2d0*(1.0d0-3.0d0*x*x+3.0d0*x**4-x**6)
			rfunc2(r) = r_outer*r_inner/radius(r)-r_inner
		Enddo
		!write(6,*)'rf max ', maxval(rfunc1)

		! We put our temporary field in spectral space
		Call tempfield%init(field_count = fcount, config = 's2b')		
		Call tempfield%construct('s2b')		

		! Set the ell = 0 temperature and the real part of Y44		
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			tempfield%s2b(mp)%data(:,:) = 0.0d0			
			Do l = m, l_max
				if ( (l .eq. 4) .and. (m .eq. 4) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc1(r)
					Enddo
				endif

				if ( (l .eq. 0) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc2(r)*sqrt(4.0d0*pi)
					Enddo
				endif
			Enddo
		Enddo
		DeAllocate(rfunc1,rfunc2)

		Call tempfield%reform() ! goes to p1b
		If (chebyshev) Then
			! we need to load the chebyshev coefficients, and not the physical representation into the RHS
			Call tempfield%construct('p1a')
			Call Cheby_To_Spectral(tempfield%p1b,tempfield%p1a)
			tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
			Call tempfield%deconstruct('p1a')
		Endif
		! Set temperature.  Leave the other fields alone
		Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

		Call tempfield%deconstruct('p1b')
	End Subroutine Benchmark_Init_Hydro
	Subroutine Diffusion_Init_Hydro()
		Implicit None
		Real*8, Allocatable :: rfunc1(:), rfunc2(:)
		Real*8 :: x
		Integer :: i, r, l, m, mp
		Integer :: fcount(3,2)
		
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 1

		Allocate(rfunc1(my_r%min: my_r%max))
		Allocate(rfunc2(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = (radius(r)-r_inner)/(r_outer-r_inner)
			!rfunc1(r) = 0.2d0*(1.0d0-3.0d0*x*x+3.0d0*x**4-x**6)
			rfunc2(r) = 0.5d0*(1.0d0-cos(2.0d0*pi*x))
		Enddo


		! We put our temporary field in spectral space
		Call tempfield%init(field_count = fcount, config = 's2b')		
		Call tempfield%construct('s2b')		

		! Set the ell = 0 temperature and the real part of Y44		
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			tempfield%s2b(mp)%data(:,:) = 0.0d0			
			Do l = m, l_max
				if ( (l .eq. 4) .and. (m .eq. 4) ) Then
					Do r = my_r%min, my_r%max
						!tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc2(r)
					Enddo
				endif

				if ( (l .eq. 0) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc2(r)
					Enddo
				endif
			Enddo
		Enddo
		DeAllocate(rfunc1,rfunc2)

		Call tempfield%reform() ! goes to p1b
		If (chebyshev) Then
			! we need to load the chebyshev coefficients, and not the physical representation into the RHS
			Call tempfield%construct('p1a')
			Call Cheby_To_Spectral(tempfield%p1b,tempfield%p1a)
			tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
			Call tempfield%deconstruct('p1a')
		Endif
		! Set temperature.  Leave the other fields alone
		Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

		Call tempfield%deconstruct('p1b')
	End Subroutine Diffusion_Init_Hydro

	Subroutine Temperature_Blob()
		Implicit None
		Real*8, Allocatable :: rfunc(:)
		Real*8 :: x
		Integer :: r,t,p
		Integer :: fcount(3,2)
		Real*8 :: tcen,w,theta, phi
		
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 1


		Allocate(rfunc(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = (radius(r)-r_inner)/(r_outer-r_inner)
			rfunc(r) = 0.5d0*(1.0d0-cos(2.0d0*pi*x))*temp_amp
		Enddo


		! We put our temporary field in Physical Space initially
		Call tempfield%init(field_count = fcount, config = 'p3b')		
		Call tempfield%construct('p3b')		
		!write(6,*)'constructed'

		w = temp_w ! 0.3d0 by default
		tcen = pi/2.0d0
		Do t = my_theta%min, my_theta%max
			theta = acos(costheta(t))
			Do r = my_r%min, my_r%max
				Do p = 1, n_phi
					phi = 2.0d0*pi/(n_phi)*p
					tempfield%p3b(p,r,t,1)= rfunc(r)*exp(-((phi-tcen)**2)/w)*exp(-((theta-tcen)**2)/w)
				Enddo
			Enddo
		Enddo

		



		DeAllocate(rfunc)

		Call fft_to_spectral(tempfield%p3b, rsc = .true.)

		Call tempfield%reform()	! Move to p2b
		Call tempfield%construct('s2b')
		Call Legendre_Transform(tempfield%p2b,tempfield%s2b)	

		Call tempfield%deconstruct('p8b')
		tempfield%config = 's2b'
		Call tempfield%reform() ! goes to p1b

		! Set temperature.  Leave the other fields alone
		Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

		Call tempfield%deconstruct('p1b')
	End Subroutine Temperature_Blob

	Subroutine Convective_Init()
		Implicit None
		Real*8, Allocatable :: rfunc(:), rfunc1(:), rfunc2(:)
		Real*8 :: x
		Integer :: r,t,p, l,m,mp
		Integer :: fcount(3,2)
		Real*8 :: tcen,w,theta, phi
		
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 1


		Allocate(rfunc(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = (radius(r)-r_inner)/(r_outer-r_inner)
			rfunc(r) = 0.5d0*(1.0d0-cos(2.0d0*pi*x))*temp_amp
		Enddo

		Allocate(rfunc1(my_r%min: my_r%max))
		Allocate(rfunc2(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = 2.0d0*radius(r)-r_inner-r_outer
			rfunc1(r) = 0.2d0*(1.0d0-3.0d0*x*x+3.0d0*x**4-x**6)
			rfunc2(r) = r_outer*r_inner/radius(r)-r_inner
		Enddo


		! We put our temporary field in Physical Space initially
		Call tempfield%init(field_count = fcount, config = 'p3b')		
		Call tempfield%construct('p3b')		
		!write(6,*)'constructed'

		w = temp_w ! 0.3d0 by default
		tcen = pi/2.0d0
		Do t = my_theta%min, my_theta%max
			theta = acos(costheta(t))
			Do r = my_r%min, my_r%max
				Do p = 1, n_phi
					phi = 2.0d0*pi/(n_phi)*p
					tempfield%p3b(p,r,t,1)= 0.0d0*rfunc(r)*exp(-((phi-tcen)**2)/w)*exp(-((theta-tcen)**2)/w)
				Enddo
			Enddo
		Enddo

		



		DeAllocate(rfunc)

		Call fft_to_spectral(tempfield%p3b, rsc = .true.)

		Call tempfield%reform()	! Move to p2b
		Call tempfield%construct('s2b')
		Call Legendre_Transform(tempfield%p2b,tempfield%s2b)	

		! Set the ell = 0 temperature and the real part of Y44		
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)			
			Do l = m, l_max
				if ( (l .eq. 4) .and. (m .eq. 4) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = tempfield%s2b(mp)%data(l,r-my_r%min+1) + &
                            temp_amp*rfunc1(r)
					Enddo
				endif

				if ( (l .eq. 0) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = tempfield%s2b(mp)%data(l,r-my_r%min+1)+ &
                            & rfunc2(r)*sqrt(4.0d0*pi)
					Enddo
				endif
			Enddo
		Enddo
		DeAllocate(rfunc1,rfunc2)


		Call tempfield%deconstruct('p8b')
		tempfield%config = 's2b'
		Call tempfield%reform() ! goes to p1b

		! Set temperature.  Leave the other fields alone
		Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

		Call tempfield%deconstruct('p1b')
	End Subroutine Convective_Init

	Subroutine Random_Init()
		Implicit None
		Integer :: ncombinations, i, m, r, seed(1), mp,n, l, ind1, ind2
		Integer :: mode_count, my_mode_start, my_mode_end, fcount(3,2)
		Real*8, Allocatable :: rand(:), rfunc1(:), rfunc2(:), lpow(:)
		Real*8 :: amp, phase, lmid, alpha,x
		Real*8, Allocatable :: new_temp(:)
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 1

		!////////////////////////
		Allocate(rfunc1(my_r%min: my_r%max))
		Allocate(rfunc2(my_r%min: my_r%max))

		Do r = my_r%min, my_r%max
			x = 2.0d0*radius(r)-r_inner-r_outer
			rfunc1(r) = 0.2d0*(1.0d0-3.0d0*x*x+3.0d0*x**4-x**6)
			rfunc2(r) = r_outer*r_inner/radius(r)-r_inner
		Enddo
		If (custom_t) Then
			! Processors owning m = 0 (and thus ell = 0)
			! Should read a custom temperature file and replace rfunc2
			Do mp = 1, my_mp%max		
				m = m_values(mp)
				If (m .eq. 0) Then
					Allocate(new_temp(1:N_R))
					Call Load_Radial_Profile(custom_t_file,new_temp)
					rfunc2(my_r%min:my_r%max) = new_temp(my_r%min:my_r%max)
					DeAllocate(new_temp)
				Endif
			Enddo
		Endif


		! We put our temporary field in spectral space
		Call tempfield%init(field_count = fcount, config = 's2b')		
		Call tempfield%construct('s2b')		


		!///////////////////////
		ncombinations = 0
		Do m = 0, l_max
			ncombinations = ncombinations+ (l_max-m+1)
		Enddo
			!Set up the random phases and amplitudes
			Allocate(rand(1:ncombinations*2))
		If (my_rank .eq. 0) Then
			Call system_clock(seed(1))		
			Call random_seed()
			Call random_number(rand)


			Do i = 1, ncombinations
				rand(i) = 2*temp_amp*(rand(i)-0.5d0)		! first half of rand contains the amplitude
			Enddo
			! We leave the second half alone (contains phases)

			! Send rand
			Do n = 1, ncpu -1
					Call send(rand, dest = n,tag=init_tag, grp=pfi%gcomm)
			Enddo
		Else
			! receive rand
				Call receive(rand, source= 0,tag=init_tag,grp = pfi%gcomm)
		Endif	

		! Everyone establishes their range of random phases		
		mode_count = 0
		Do mp = 1, my_mp%max		
			if (mp .eq. my_mp%min) then
				my_mode_start = mode_count+1
			endif
			m = m_values(mp)
			mode_count = mode_count + (l_max-m+1)
			if (mp .eq. my_mp%max) then
				my_mode_end = mode_count
			endif
		Enddo

		Allocate(lpow(0:l_max))
		lmid = l_max/2.0d0
		alpha = lmid/3.0d0
		Do l = 0, l_max
				lpow(l) = exp(- ((l-lmid)/alpha )**2)
		Enddo


		ind1 = my_mode_start
		ind2 = ind1+ncombinations
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)			
			Do l = m, l_max
				tempfield%s2b(mp)%data(l,:) = 0.0d0
				if ( (l .eq. 0) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc2(r)*sqrt(4.0d0*pi)
					Enddo
				endif
					amp = rand(ind1)*lpow(l)
					phase = rand(ind2)
					ind1 = ind1+1
					ind2 = ind2+1
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = tempfield%s2b(mp)%data(l,r-my_r%min+1) + &
                            amp*rfunc1(r)*phase
						tempfield%s2b(mp)%data(l,r-my_r%min+1+my_r%delta) = tempfield%s2b(mp)%data(l,r-my_r%min+1+my_r%delta) + &
                            amp*rfunc1(r)*(1.0d0-phase)
					Enddo
			Enddo
		Enddo
		DeAllocate(rfunc1,rfunc2, lpow)

		Call tempfield%reform() ! goes to p1b

		! Set temperature.  Leave the other fields alone
		If (chebyshev) Then
			! we need to load the chebyshev coefficients, and not the physical representation into the RHS
			Call tempfield%construct('p1a')
			Call Cheby_To_Spectral(tempfield%p1b,tempfield%p1a)
			tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
			Call tempfield%deconstruct('p1a')
		Endif
		Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

		Call tempfield%deconstruct('p1b')

	End Subroutine Random_Init


	!////////////////////////////////////
	!  Magnetic Initialization
	Subroutine Benchmark_Insulating_Init()
		Implicit None
		Real*8, Allocatable :: rfunc1(:), rfunc2(:)
		Real*8 :: x, nrm1, nrm2
		Integer :: i, r, l, m, mp, roff
		Integer :: fcount(3,2)
		type(SphericalBuffer) :: tempfield
		fcount(:,:) = 2

		Allocate(rfunc1(my_r%min: my_r%max))
		Allocate(rfunc2(my_r%min: my_r%max))

		nrm2 = sqrt(pi/5.0d0)*(4.0d0/3.0d0)*5.0d0
		nrm1 = (5.0d0/8.0d0)*sqrt(pi/3.0d0)
		Do r = my_r%min, my_r%max

			rfunc1(r) = 8.0d0*r_outer-6.0d0*radius(r)		! functional form for c_1_0
			rfunc1(r) = rfunc1(r)-2.0d0*(r_inner**4)/(radius(r)**3)
			rfunc1(r) = rfunc1(r)*nrm1*radius(r)**2

			rfunc2(r) = sin(pi*(radius(r)-r_inner))*radius(r)*nrm2  ! function form for a_2_0
		Enddo
		!write(6,*)'rf max ', maxval(rfunc1)

		! We put our temporary field in spectral space
		Call tempfield%init(field_count = fcount, config = 's2b')		
		Call tempfield%construct('s2b')		
		roff= 2*my_r%delta
		! Set the ell = 0 temperature and the real part of Y44		
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			tempfield%s2b(mp)%data(:,:) = 0.0d0			
			Do l = m, l_max
				if ( (l .eq. 1) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1) = rfunc1(r)
					Enddo
				endif

				if ( (l .eq. 2) .and. (m .eq. 0) ) Then
					Do r = my_r%min, my_r%max
						tempfield%s2b(mp)%data(l,r-my_r%min+1+roff) = rfunc2(r)
					Enddo
				endif
			Enddo
		Enddo
		DeAllocate(rfunc1,rfunc2)

		Call tempfield%reform() ! goes to p1b
		If (chebyshev) Then
			! we need to load the chebyshev coefficients, and not the physical representation into the RHS
			Call tempfield%construct('p1a')
			Call Cheby_To_Spectral(tempfield%p1b,tempfield%p1a)
			tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
			Call tempfield%deconstruct('p1a')
		Endif

		Call Set_RHS(ceq,tempfield%p1b(:,:,:,1))
		Call Set_RHS(aeq,tempfield%p1b(:,:,:,2))

		Call tempfield%deconstruct('p1b')
	End Subroutine benchmark_insulating_init

	Subroutine Load_Radial_Profile(profile_file,profile_out)
		Implicit None
		Character*120, Intent(In) :: profile_file
		Real*8, Intent(InOut) :: profile_out(1:)
		Real*8, Allocatable :: radius_in(:), profile_in(:), spy2(:)
		Real*8 :: min_r_in, max_r_in, min_p_in, max_p_in
		Real*8 :: splx, sply
		Integer :: nr_in, r
		! Reads in a 1-D radial profile of some quantity,
		! interpolates it (using cubic splines) to the current grid, 
		! and stores it in profile_out.

      Open(unit=892, file = custom_t_file, form = 'unformatted', status = 'old')
      Read(892)nr_in
      Allocate(profile_in(1:nr_in))
      Allocate(radius_in(1:nr_in))
      Read(892)(radius_in(r),r=1,nr_in)
      Read(892)(profile_in(r),r = 1, nr_in)
      Close(892)

     !--------------------------------------------------------------
      ! Interpolate onto our grid (assume custom_radius and custom_entropy are backwards
      min_r_in = radius_in(nr_in)
      max_r_in = radius_in(1)
      max_p_in = profile_in(nr_in)
      min_p_in = profile_in(1)

      Allocate(spy2(1:nr_in))
      spy2(1:nr_in) = 0.0d0
      profile_out(1:N_R) = 0.0d0
      Call Spline(radius_in, profile_in, nr_in, 2.0D30, 2.0D30, spy2)
      Do r = 1, N_R
         If ( (radius(r) .le. max_r_in) .and. (radius(r) .ge. min_r_in) ) Then
            splx = radius(r)
            Call Splint(radius_in, profile_in,spy2,nr_in, splx, sply)
            profile_out(r) = sply
         Endif
      Enddo

      ! Take care of any out of bounds radial values
      Do r = 1, N_R
         If (radius(r) .ge. max_r_in) Then
            profile_out(r) = min_p_in
         Endif
         If (radius(r) .le. min_r_in) Then
            profile_out(r) = max_p_in
         Endif
      Enddo

		DeAllocate(radius_in, profile_in, spy2)
	End Subroutine Load_Radial_Profile

	!/////////////////////////////////////////////////////
	! Numerical Recipes Routines for Spline Interpolation
	Subroutine Spline(x,y,n,yp1,ypn,y2)
		! From Numerical Recipes in Fortran
		Integer:: n, NMAX
		Real(8) :: yp1, ypn, x(n), y(n), y2(n)
		PARAMETER (NMAX = 10000)
		Integer :: i, k
		Real(8) :: p, qn, sig, un, u(NMAX)

		If (yp1 .gt. 0.99D30) Then
			y2(1) = 0.0D0
			u(1) = 0.0D0
		Else
			y2(1) = -0.5D0
			u(1) = ( 3.0D0 / ( x(2)-x(1) ) ) * ( (y(2)-y(1))/(x(2)-x(1)) -yp1 )
		Endif

		Do i = 2, n-1
			sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
			p = sig*y2(i-1)+2.0D0
			y2(i) = (sig-1.0D0)/p
			u(i) = (6.0D0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) &
				& /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
		Enddo

		If (ypn .gt. 0.99D30) Then
			qn = 0.0D0
			un = 0.0D0
		Else
			qn = 0.5D0
			un = (3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
		Endif

		y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)

		Do k = n-1, 1, -1
			y2(k) = y2(k)*y2(k+1)+u(k)
		Enddo
		Return
	End Subroutine Spline

	Subroutine Splint(xa,ya,y2a,n,x,y)
		! From Numerical Recipes in Fortan
		Integer :: n
		Real(8) x, y, xa(n), y2a(n), ya(n)
		Integer :: k, khi, klo
		Real(8) a,b,h

		klo = 1
		khi = n
1   If ( (khi-klo) .gt. 1) Then
			k = (khi+klo)/2
			If (xa(k) .lt. x) Then		! if xa is in ascending order, change lt to gt
				khi = k
			else
         	klo = k
			endif
			Goto 1
		Endif

    	h = xa(khi)-xa(klo)
    	If (h .eq. 0.0D0) pause 'bad xa input in splint'
    	a = (xa(khi)-x)/h
    	b = (x-xa(klo))/h

    	y = a*ya(klo)+ b*ya(khi)+ &
			& ( (a**3-a)*y2a(klo)+(b**3-b)*y2a(khi) )*(h**2)/6.0D0

		Return
	End Subroutine Splint
End Module Initial_Conditions
