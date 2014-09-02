Module ProblemSize
	Use Parallel_Framework, Only : pfi, Load_Config, Spherical
	Use Finite_Difference, Only  : Initialize_Derivatives, Rescale_Grid_FD
	Use Legendre_Polynomials, Only : Initialize_Legendre,coloc
	Use Spectral_Derivatives, Only : Initialize_Angular_Derivatives
	Use Controls, Only : Chebyshev, use_parity, read_argv
	Use Chebyshev_Polynomials, Only : Initialize_Chebyshev, Rescale_Grid_CP
	Use Timers
	Implicit None

	!//////////////////////////////////////////////////////////////
	! Processor Configuration
	Integer :: ncpu = 1, nprow = 1, npcol =1 , npout = 1
	Integer :: my_rank, my_row_rank, my_column_rank ! rank *within* row and rank *within* column
	!//////////////////////////////////////////////////////////////
	! Horizontal Grid Variables
	Integer              :: n_theta, n_phi
	Integer              :: l_max, m_max, n_l, n_m
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
	Real*8, Allocatable :: Radius(:), R_squared(:), One_Over_R(:)
	Real*8, Allocatable :: Two_Over_R(:), OneOverRSquared(:), Delta_R(:)
	Real*8, Allocatable :: ovrsq_repeated(:),ovr_repeated(:)
	Type(Load_Config)   :: my_r

	Namelist /ProblemSize_Namelist/ n_r,n_theta, nprow, npcol,rmin,rmax,npout
Contains

	Subroutine Init_ProblemSize()
		Implicit None
		Integer :: ppars(1:10)
		Integer :: tmp,r, l, i, ii
		Integer, Allocatable :: m_vals(:)
		Real*8 :: ell
		Character*10 :: arg, arg2
		!rmin = (7.0d0)/13.0d0		! benchmark cluge
		!rmax = (20.0d0)/13.0d0

		n_phi = 2*n_theta
		If (dealias) Then
			l_max = (2*n_theta-1)/3
		Else
			l_max = n_theta-1
		Endif
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
		Call pfi%init(ppars)
		Call Initialize_Timers()
		Call StopWatch(init_time)%startclock()
		Call Map_Indices()
		my_rank = pfi%gcomm%rank
		my_row_rank = pfi%rcomm%rank
		my_column_rank = pfi%ccomm%rank
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
			Open(101, file = 'grid', status='replace', form = 'unformatted')
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
		Integer :: r
		real*8 :: uniform_dr, tmp1, tmp2

		Allocate(Delta_r(1:N_R))
		Allocate( Radius(1:N_R))
		If (chebyshev) Then
			grid_type = 2
			Call Initialize_Chebyshev(radius,rmin,rmax)
			Delta_r(1) = radius(1)-radius(2)
			Do r = 2, N_R
				Delta_r(r) = radius(r)-radius(r-1)
			Enddo
		Else

			Select Case (grid_type)
				Case (1)	! Uniform Grid
					Radius(N_R) = rmin	! Follow ASH convention of reversed radius
					uniform_dr = 1.0d0/(N_R-1.0d0)*(rmax-rmin)
					Do r=N_R,1,-1
						Delta_r(r) = uniform_dr
						Radius(r) = (N_R-r)*uniform_dr + rmin		
					Enddo

 
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


		If (.not. chebyshev) Call Initialize_Derivatives(Radius)
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
	End Subroutine 
End Module ProblemSize
