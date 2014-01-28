Module Initial_Conditions
	Use ProblemSize
	Use Fields
	Use Parallel_Framework
	Use Fourier_Transform
	Use Legendre_Transforms, Only : Legendre_Transform 
	Use SendReceive
	Use Chebyshev_Polynomials, Only : Cheby_To_Spectral
	Implicit None
	Integer :: init_type = 1
	Integer :: init_tag = 8989
	Real*8 :: pi = 3.1415926535897932384626433832795028841972d+0
	Real*8 :: temp_amp = 1.0d0, temp_w = 0.3d0
	Namelist /Initial_Conditions_Namelist/ init_type, temp_amp, temp_w
Contains
	
	Subroutine Initialize_Fields()
		Implicit None

		! When coming out of this routine, the RHS of the equation set should contain the field values.
		! This setup is consistent with the program having just completed a time step

		! wsp%p1b should contain the Adams-Bashforth terms from the previous (not current)
		! values of WPST.  This means that wsp%p1b should be zero if this init is from scratch.
		! If this init is from restart, wsp%p1b should contain the adams bashforth terms output
		! as part of the checkpoint.

		Call wsp%init(field_count = wsfcount, config = 'p1b')		
		Call wsp%construct('p1b')	! We will always start in p1b - should do wsp%set_config('p1b')
		wsp%p1b(:,:,:,:) = 0.0d0	! All fields are zero initially

		! Allocate the Equation Set RHS
		! Set it to zero initially	
		! The equation set RHS's stays allocated throughout - it is effectively how we save the AB terms.
		Call Allocate_RHS(zero_rhs=.true.)

	
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

		! Fields are now initialized and loaded into the RHS. 
		! We are ready to enter the main loop

	End Subroutine Initialize_Fields

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
		!write(6,*)'rf max ', maxval(rfunc1)

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

End Module Initial_Conditions
