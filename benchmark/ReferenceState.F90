! REFERENCE STATE MODULE
! Contains routines for initializing the reference state structure
! Reference state structure contains all information related to background
! stratification.  It DOES NOT contain transport variable (e.g. nu, kappa) information



Module ReferenceState
	Use ProblemSize
	Implicit None
	Type ReferenceInfo
		Real*8, Allocatable :: Density(:)
		Real*8, Allocatable :: dlnrho(:)
		Real*8, Allocatable :: d2lnrho(:)

		Real*8, Allocatable :: Pressure(:)

		Real*8, Allocatable :: Temperature(:)
		Real*8, Allocatable :: dlnT(:)

		Real*8, Allocatable :: Entropy(:)
		Real*8, Allocatable :: dsdr(:)

		Real*8, Allocatable :: Gravity(:)

		Real*8 :: gamma

		Real*8 :: rho_twiddle, g_twiddle, p_twiddle, s_twiddle, t_twiddle
	End Type ReferenceInfo

	Integer :: reference_type
	
	Type(ReferenceInfo) :: ref
	Real*8 :: pressure_specific_heat  ! CP (not CV)
	Real*8 :: poly_n
   Real*8 :: poly_Nrho
   Real*8 :: poly_mass
   Real*8 :: poly_rho_i
	Real*8 :: Gravitational_Constant = 6.67d-8
	Real*8 :: rho_twiddle, g_twiddle, p_twiddle, s_twiddle, t_twiddle, length_twiddle ! also in ref structure..



	Namelist /Reference_Namelist/ reference_type,poly_n, poly_Nrho, poly_mass,poly_rho_i, pressure_specific_heat
Contains

	Subroutine Initialize_Reference()
		Implicit None
		Call Allocate_Reference_State()
		If (reference_type .eq. 1) Then
			Call Constant_Reference()
		Endif

		If (reference_type .eq. 2) Then
			Call Polytropic_Reference()
		Endif

		Call Write_Reference()
	End Subroutine Initialize_Reference

	Subroutine Allocate_Reference_State
		Implicit None
		Allocate(ref%density(1:N_R))
		Allocate(ref%pressure(1:N_R))
		Allocate(ref%temperature(1:N_R))
		Allocate(ref%entropy(1:N_R))
		Allocate(ref%gravity(1:N_R))
		Allocate(ref%dlnrho(1:N_R))
		Allocate(ref%d2lnrho(1:N_R))
		Allocate(ref%dlnt(1:N_R))
		Allocate(ref%dsdr(1:N_R))
	End Subroutine Allocate_Reference_State

	Subroutine Polytropic_Reference()
      Real*8 :: zeta_0,  c0, c1, d
      Real*8 :: rho_c, P_c, T_c,denom
      Real*8 :: beta, Gas_Constant, bottom_flux, test
      Real*8, Allocatable :: zeta(:), temp(:)
		Real*8 :: One
      Real*8 :: InnerRadius, OuterRadius
      If (my_rank .eq. 0) write(6,*)'Initializing polytropic reference state.'
      ! Adiabatic, Polytropic Reference State (see, e.g., Jones et al. 2011)
      ! The following parameters are read from the input file.
      ! poly_n
      ! poly_Nrho
      ! poly_mass
      ! poly_rho_i

      ! Note that cp must also be specified.
		InnerRadius = Radius(N_r)
		OuterRadius = Radius(1)
		

		One = 1.0d0
      !-----------------------------------------------------------
      beta = InnerRadius/OuterRadius

      denom = beta * exp(poly_Nrho / poly_n) + 1.d0
      zeta_0 = (beta+1.d0)/denom

      c0 = (2.d0 * zeta_0 - beta - 1.d0) / (1.d0 - beta)

      denom = (1.d0 - beta)**2
      c1 = (1.d0+beta)*(1.d0-zeta_0)/denom

      !-----------------------------------------------------------
      ! allocate and define zeta
      ! also rho_c, T_c, P_c

      Allocate(zeta(N_R))

      d = OuterRadius - InnerRadius    

      zeta = c0 + c1 * d / Radius

      rho_c = poly_rho_i / zeta(N_R)**poly_n

      denom = (poly_n+1.d0) * d * c1
      P_c = Gravitational_Constant * poly_mass * rho_c / denom

      T_c = (poly_n+1.d0) * P_c / (Pressure_Specific_Heat * rho_c)

      !-----------------------------------------------------------
      ! Initialize reference structure 
      ref%gamma = (poly_n+1.0D0)/(poly_n)


      Gas_Constant = (ref%Gamma-one)*Pressure_Specific_Heat/ref%Gamma

      Ref%Gravity = Gravitational_Constant * poly_mass / Radius**2

      Ref%Density = rho_c * zeta**poly_n

      Ref%dlnrho = - poly_n * c1 * d / (zeta * Radius**2)
		Ref%d2lnrho = - Ref%dlnrho*(2.0d0/Radius-c1*d/zeta/Radius**2)

      Ref%Temperature = T_c * zeta
		Ref%dlnT = -(c1*d/Radius**2)/zeta

      Ref%Pressure = P_c * zeta**(poly_n+1)

      denom = P_c**(1.d0/ref%gamma)
      Ref%Entropy = Pressure_Specific_Heat * log(denom/rho_c)

      Ref%dsdr = 0.d0

      Deallocate(zeta)



	End Subroutine Polytropic_Reference

	Subroutine Write_Reference(filename)
		Implicit None
		Character*120, Optional, Intent(In) :: filename
		Character*120 :: ref_file
		Integer :: i
		if (present(filename)) then
			ref_file = filename
		else
			ref_file = 'reference'
		endif

		If (my_rank .eq. 0) Then
			Open(unit=15,file=ref_file,form='unformatted', status='replace')
			Write(15)n_r
			Write(15)(radius(i),i=1,n_r)
			Write(15)(ref%density(i),i=1,n_r)
			Write(15)(ref%dlnrho(i),i=1,n_r)
			Write(15)(ref%d2lnrho(i),i=1,n_r)
			Write(15)(ref%pressure(i),i=1,n_r)
			Write(15)(ref%temperature(i),i=1,n_r)
			Write(15)(ref%dlnT(i),i=1,n_r)
			Write(15)(ref%dsdr(i),i=1,n_r)
			Write(15)(ref%entropy(i),i=1,n_r)			
			Write(15)(ref%gravity(i),i=1,n_r)
			Close(15)
		Endif
	End Subroutine Write_Reference

	Subroutine Constant_Reference()
			ref%density = 1.0d0
			ref%dlnrho = 0.0d0
			ref%d2lnrho = 0.0d0
			ref%pressure = 1.0d0
			ref%temperature = 1.0d0
			ref%dlnT = 0.0d0
			ref%dsdr = 0.0d0
			ref%pressure = 1.0d0
			ref%gravity = 1.0d0 ! Need to revisit this part
	End Subroutine Constant_Reference

	

End Module ReferenceState
