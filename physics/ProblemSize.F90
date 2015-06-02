Module ProblemSize
	Use Parallel_Framework, Only : pfi, Load_Config, Spherical
	Use Finite_Difference, Only  : Initialize_Derivatives, Rescale_Grid_FD
	Use Legendre_Polynomials, Only : Initialize_Legendre,coloc
	Use Spectral_Derivatives, Only : Initialize_Angular_Derivatives
	Use Controls, Only : Chebyshev, use_parity, multi_run_mode, run_cpus, my_path
	Use Chebyshev_Polynomials, Only : Initialize_Chebyshev, Rescale_Grid_CP
	Use Timers

	Implicit None

	!//////////////////////////////////////////////////////////////
	! Processor Configuration
	Integer :: ncpu = 1, nprow = 1, npcol =1 , npout = 1
	Integer :: my_rank      ! This is the rank within a run communicator
    Integer :: global_rank  ! This differs from my_rank only when multi-run mode is active
    Integer :: ncpu_global  ! Same as ncpu unless multi-run mode is active
    Integer :: my_row_rank, my_column_rank ! rank *within* row and rank *within* column

	!//////////////////////////////////////////////////////////////
	! Horizontal Grid Variables
	Integer              :: n_theta, n_phi
    Integer              :: l_max = -1
	Integer              :: m_max, n_l, n_m
	Logical              :: dealias = .True.
	Integer, Allocatable :: m_values(:)
	Real*8, Allocatable  :: l_l_plus1(:), over_l_l_plus1(:)
	Real*8, Allocatable  :: costheta(:), sintheta(:), cos2theta(:), sin2theta(:), cottheta(:), csctheta(:)
	Type(Load_Config)    :: my_mp,  my_theta


	!//////////////////////////////////////////////////////////////
	!  Radial Grid Variables
	Integer             :: n_r, tnr
	Integer             :: grid_type = 1
	Real*8              :: rmin, rmax, r_inner, r_outer
    Real*8              :: stretch_factor = 0.0d0
	Real*8, Allocatable :: Radius(:), R_squared(:), One_Over_R(:)
	Real*8, Allocatable :: Two_Over_R(:), OneOverRSquared(:), Delta_R(:)
	Real*8, Allocatable :: ovrsq_repeated(:),ovr_repeated(:), radial_integral_weights(:)
    Integer :: precise_bounds = -1
    Type(Load_Config)   :: my_r

    !///////////////////////////////////////////////////////////////
    ! Radial Grid Variables Related to the Finite-Element Approach
    Integer :: fensub = -1   !  Number of subdomains to divide N_R into
    Integer :: fencheby = -1 ! Number of chebyshev modes per subdomain
	Logical :: finite_element = .false.

	Namelist /ProblemSize_Namelist/ n_r,n_theta, nprow, npcol,rmin,rmax,npout, & 
            &  precise_bounds,grid_type, stretch_factor, fensub,fencheby, l_max
Contains

	Subroutine Init_ProblemSize()
		Implicit None
		Integer :: ppars(1:10)
		Integer :: tmp,r, l
		Integer, Allocatable :: m_vals(:)
		Real*8 :: ell
        Character*120 :: grid_file
        Integer :: cpu_tmp(1)
        if (precise_bounds .eq. 1) Then
		    rmin = (7.0d0)/13.0d0		! Benchmark bounds have infinite decimal places
		    rmax = (20.0d0)/13.0d0      ! We override (if desired) the inputs for accuracy
        Endif
        if ( (fensub .gt. 0) .and. (fencheby .gt. 0) ) Then
            !Overwrite n_r
            n_r =fensub*fencheby 
            finite_element = .true.
        Endif

        If (l_max .le. 0) Then


		    If (dealias) Then
			    l_max = (2*n_theta-1)/3
		    Else
			    l_max = n_theta-1
		    Endif

        Else
            !base n_theta on l_max
            If (dealias) Then
                n_theta = (3*(l_max+1))/2
            Else
                n_theta = l_max+1
            Endif
        Endif
		n_phi = 2*n_theta
		m_max = l_max
		n_l = l_max+1
		n_m = m_max+1

		Allocate(l_l_plus1(0:l_max))
		Allocate(over_l_l_plus1(0:l_max))
		over_l_l_plus1(0) = 1.0d0 ! This keeps us from having to worry about dividing by zero
		Do l = 0, l_max
			ell = l*1.0d0
			l_l_plus1(l) = ell*(ell+1)
			if (l .ne. 0) over_l_l_plus1(l) = 1.0d0/l_l_plus1(l)
		Enddo
		!///////////////////////////////////////
		!  Check the command line to see if 
		!  any processor counts were specified.
		!  If so, these counts overwrite values
		!  read in through main input 
		
		
		ncpu = nprow*npcol
		!///////////////////////////////////////
		! Initialize load balancing
		ppars(1) = Spherical
		ppars(2) = n_r
		ppars(3) = n_r
		ppars(4) = n_theta
		ppars(5) = n_l
		ppars(6) = n_phi
		ppars(7) = n_m		
		ppars(8) = ncpu
		ppars(9) = nprow
		ppars(10) = npcol
        If (multi_run_mode) Then
    		Call pfi%init(ppars,run_cpus)
        Else
            cpu_tmp(1) = ncpu
    		Call pfi%init(ppars,cpu_tmp)
        Endif

		Call Initialize_Timers()
		Call StopWatch(init_time)%startclock()
		Call StopWatch(walltime)%startclock()
		Call Map_Indices()
		my_rank = pfi%gcomm%rank
		my_row_rank = pfi%rcomm%rank
		my_column_rank = pfi%ccomm%rank
        !//////////////////////////////////////////////////
        !Provide quick notification that n_r may have changed
        if ( finite_element) Then
            if (my_rank .eq. 0) Write(6,*)"Finite Elements.  N_R is: ", n_r
        Endif

		!//////////////////////////////////////////////////
		! Intialize Legendre Transforms & Horizontal Grid
		tmp = my_mp%delta
		allocate(m_vals(1:tmp))
		m_vals(:) = m_values(my_mp%min:my_mp%max)
		
		Call Initialize_Legendre(n_theta,l_max,m_vals,use_parity)
		tmp = my_r%delta
		Call Initialize_Angular_Derivatives(m_vals,l_max,tmp)
		DeAllocate(m_vals)


		Allocate(costheta(1:n_theta),cos2theta(1:n_theta))
		Allocate(sintheta(1:n_theta),sin2theta(1:n_theta))
		Allocate(csctheta(1:n_theta), cottheta(1:n_theta))
		costheta  = coloc	! coloc computed in init_legendre
		cos2theta = costheta*costheta
		sin2theta = 1-cos2theta
		sintheta  = sqrt(sin2theta)
		csctheta = 1/sintheta
		cottheta = costheta/sintheta

		Call Initialize_Radial_Grid()
		r_inner = rmin
		r_outer = rmax
		if (pfi%gcomm%rank .eq. 0) then
            grid_file = Trim(my_path)//'grid'
			Open(101, file = grid_file, status='replace', form = 'unformatted')
			Write(101) n_r
			Write(101) (radius(r),r=1,n_r)
			Write(101) n_theta
			Write(101) (costheta(r),r=1,n_theta)
			close(101)
		endif
		!////////////////////////////////////////
		! Make some repeated radius arrays so we can 
		! easily divide the real and imaginary parts
		! of our rlm buffers
		Allocate(ovrsq_repeated(1:2*my_r%delta))
		ovrsq_repeated(1:my_r%delta) = 1.0d0/r_squared(my_r%min:my_r%max)
		ovrsq_repeated(my_r%delta+1:2*my_r%delta) = ovrsq_repeated(1:my_r%delta)

		Allocate(ovr_repeated(1:2*my_r%delta))
		ovr_repeated(1:my_r%delta) = 1.0d0/radius(my_r%min:my_r%max)
		ovr_repeated(my_r%delta+1:2*my_r%delta) = ovr_repeated(1:my_r%delta)
		tnr = 2*my_r%delta
		

	End Subroutine Init_ProblemSize

	Subroutine Map_Indices()
		Implicit None
		Allocate(m_values(1:n_m))
		my_mp    = pfi%my_3s
		my_r     = pfi%my_1p
		my_theta = pfi%my_2p
		m_values = pfi%inds_3s

	End Subroutine Map_Indices

	Subroutine Initialize_Radial_Grid()
		Implicit None
		Integer :: r, nthr,i,j 
		real*8 :: uniform_dr, arg, pi_over_N, rmn, rmx, delta, scaling
        real*8 :: delr0
        Real*8 ::	Pi  = 3.1415926535897932384626433832795028841972d0

        !////////////////////////////////////////
        ! Variables for FE approach
        Real*8, Allocatable :: xtemp(:)
        Real*8 :: dsub, offset, xmax, xmin
        Integer :: istart, iend


		nthr = pfi%nthreads
		Allocate(Delta_r(1:N_R))
		Allocate( Radius(1:N_R))
        Allocate(Radial_Integral_Weights(1:N_R))
		If (chebyshev) Then
			grid_type = 2
			Call Initialize_Chebyshev(radius,rmin,rmax,radial_integral_weights, nthr)
			Delta_r(1) = radius(1)-radius(2)
			Do r = 2, N_R
				Delta_r(r) = radius(r)-radius(r-1)
			Enddo

        Else If (finite_element) Then
            Write(6,*)"Initializing grid..."
            Allocate(xtemp(1:fencheby))
            dsub = (rmax-rmin)/DBLE(fensub)
            xmin = 0.0d0
            xmax = dsub
            
            Call Initialize_Chebyshev(xtemp,xmin,xmax,radial_integral_weights, nthr)
            offset = rmax-xtemp(1)
            istart = 1
            iend = istart+fencheby-1
            Do j = 1, fensub
                radius(istart:iend) = xtemp(1:fencheby)+offset
                offset = radius(iend)-xtemp(1)
                istart = istart+fencheby
                iend = iend+fencheby
            Enddo 
            Do j = 1, n_r
                if (MOD(j,fencheby) .eq. 0) Then
                    delta_r(j) = radius(j-1)-radius(j)
                Else
                    delta_r(j) = radius(j)-radius(j+1)
                Endif
            Enddo
            If (my_rank .eq. 0) Then
                Do j = 1, n_r
                    write(6,*)j," : ", radius(j), " , ", delta_r(j)
                    if (MOD(j,fencheby) .eq. 0) Write(6,*)"---------------------"
                Enddo
            Endif
		Else

			Select Case (grid_type)
				Case (1)	! Uniform Grid
					Radius(N_R) = rmin	! Follow ASH convention of reversed radius
					uniform_dr = 1.0d0/(N_R-1.0d0)*(rmax-rmin)
					Do r=N_R,1,-1
						Delta_r(r) = uniform_dr
						Radius(r) = (N_R-r)*uniform_dr + rmin		
					Enddo

                Case (2)  ! Chebyshev Grid

		            pi_over_N = pi/(N_r*1.0d0)
		            arg = (0.5d0)*pi_over_N
		            Do i = 1, N_R
			            radius(i) = cos(arg)
			            arg = arg+pi_over_N
		            Enddo
                    delta = rmax-rmin
                    scaling = (radius(1)-radius(n_r) )/ delta
                    radius = radius/scaling
                    rmn = minval(radius)
                    radius(:) = radius(:)-rmn+rmin                    

                Case (3) ! Stretched Grid -  high res near boundaries
                         ! Each cell is (1+stretch_factor) bigger than the one before it
                         ! n_r is assumed to be even for this to work
                    delta = rmax-rmin
                    arg = 0.0d0
                    radius(1) = 0.0d0
                    r = (n_r/2)+1

                    Do i = 0, r-2
                        arg = arg + (1.0d0+stretch_factor)**i
                    Enddo
                    delr0 = delta/arg  ! This is the grid spacing at the outer boundary

                    ! Set up the top half of the grid (1 through nr/2+1)                    
                    Do i = 1, r-1
                        arg = delr0*( (1.0+stretch_factor)**i )
                        radius(i+1) = radius(i)+arg
                    Enddo

                    ! Reflect to get the other half.
                    Do i = 1, (n_r/2)-1
                        arg = radius(r)-radius(r-i )
                        radius(r+i )=radius(r)+arg
                    Enddo

                    !Finally, rescale the grid
                    rmx = maxval(radius)
                    radius(:) = radius(:)*delta/rmx
                    radius(:) = delta-radius(:)+rmin
                    !If (my_rank .eq. 0) Then
                    !    Do i = 1, n_r
                    !       Write(6,*)'Radius : ', radius(i)
                    !    Enddo
                    !Endif

				Case Default	! Uniform Grid - Same as case 1
					Radius(N_R) = rmin
					uniform_dr = 1.0d0/(N_R-1.0d0)*(rmax-rmin)
					Do r=N_R,1,-1
						Delta_r(r) = uniform_dr
						Radius(r) = (N_R-r)*uniform_dr + rmin		
					Enddo
			End Select
		Endif

		! Compute delta_r for CFL
		do r = N_R-1, 1, -1
			Delta_r(r) = Radius(r)-Radius(r+1)
		enddo
		Delta_r(N_R) = Delta_r(N_R-1) ! Duplicated extra value


		Allocate(OneOverRSquared(1:N_R),r_squared(1:N_R),One_Over_r(1:N_R),Two_Over_r(1:N_R))
		R_squared       = Radius**2
      One_Over_R      = (1.0d0)/Radius
      Two_Over_R      = (2.0d0)/Radius
      OneOverRSquared = (1.0d0)/r_Squared


		If (.not. chebyshev) Call Initialize_Derivatives(Radius,radial_integral_weights)
	End Subroutine Initialize_Radial_Grid
    Subroutine Rescale_Grid_and_Derivatives(length_scale)
        Implicit None
        Real*8, Intent(In) :: length_scale
        Radius = radius/length_scale
        R_squared       = Radius**2
        One_Over_R      = (1.0d0)/Radius
        Two_Over_R      = (2.0d0)/Radius
        OneOverRSquared = (1.0d0)/r_Squared
        ovr_repeated = ovr_repeated*length_scale
        ovrsq_repeated = ovrsq_repeated*length_scale
        If (chebyshev) Then
            Call Rescale_Grid_CP(length_scale)
        Else
            Call Rescale_Grid_FD(length_scale)
        Endif
        ! We need to think about how to rescale the radial integration weights
    End Subroutine 
End Module ProblemSize
