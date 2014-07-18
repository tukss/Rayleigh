Module TransportCoefficients
	Use ProblemSize
	Use ReferenceState
	Use Controls
	Implicit None
	Real*8, Allocatable :: nu(:), kappa(:), eta(:)
	Real*8, Allocatable :: dlnu(:), dlnkappa(:), dlneta(:)

	!//////////
	! These need to be relocated to transportcoefficients.F90 shortly
	Real*8, Allocatable :: W_Diffusion_Coefs_0(:), W_Diffusion_Coefs_1(:)
	Real*8, Allocatable :: dW_Diffusion_Coefs_0(:), dW_Diffusion_Coefs_1(:), dW_Diffusion_Coefs_2(:)
	Real*8, Allocatable :: S_Diffusion_Coefs_1(:), Z_Diffusion_Coefs_0(:), Z_Diffusion_Coefs_1(:)
	Real*8, Allocatable ::  A_Diffusion_Coefs_1(:)

	Integer :: kappa_type =1, nu_type = 1, eta_type = 1
	Real*8 :: nu_top = 1.0d0, kappa_top = 1.0d0, eta_top = 1.0d0
	Real*8 :: nu_power = 0, eta_power = 0, kappa_power = 0

	Namelist /Transport_Namelist/ nu_type, kappa_type, eta_type, nu_power, kappa_power, eta_power, &
			& nu_top, kappa_top, eta_top

Contains

	Subroutine Compute_Diffusion_Coefs()
		!////////////////////////////////////////
		! W Coefficients for W Equation
		Allocate(W_Diffusion_Coefs_0(1:N_R))
		Allocate(W_Diffusion_Coefs_1(1:N_R))
		W_Diffusion_Coefs_0 = 	-nu*(4.0d0/3.0d0)*( dlnu*ref%dlnrho + ref%d2lnrho + ref%dlnrho/radius + &
			& 3.0d0*dlnu/radius )
		W_Diffusion_Coefs_1 = nu*(2.0d0*dlnu-ref%dlnrho/3.0d0)

		!/////////////////////////////////////
		! W Coefficients for dWdr equation
		Allocate(DW_Diffusion_Coefs_0(1:N_R))
		Allocate(DW_Diffusion_Coefs_1(1:N_R))
		Allocate(DW_Diffusion_Coefs_2(1:N_R))
		DW_Diffusion_Coefs_2 = dlnu-ref%dlnrho
		DW_Diffusion_Coefs_1 = ref%d2lnrho+(2.0d0)/radius*ref%dlnrho+2.0d0/radius*dlnu+dlnu*ref%dlnrho
		DW_Diffusion_Coefs_0 = 2.0d0/radius+2.0d0*ref%dlnrho/3.0d0+dlnu
		!////////////////////////////////////////
		! S Coefficients for S Equation
		Allocate(S_Diffusion_Coefs_1(1:N_R))
		S_diffusion_Coefs_1 = kappa*(dlnkappa+ref%dlnrho+ref%dlnT)
		!////////////////////////////////////////
		! Z Coefficients for the Z Equation
		Allocate(Z_Diffusion_Coefs_0(1:N_R))
		Allocate(Z_Diffusion_Coefs_1(1:N_R))
		Z_Diffusion_Coefs_0 = -nu*( 2.0d0*dlnu/radius + ref%dlnrho*dlnu + &
			& ref%d2lnrho+2.0d0*ref%dlnrho/radius)
		Z_Diffusion_Coefs_1 = nu*(dlnu-ref%dlnrho)

		!////////////////////////////////////////
		! A (vector potential) Coefficients
		If (magnetism) Then
			Allocate(A_Diffusion_Coefs_1(1:N_R))
			A_Diffusion_Coefs_1 = eta*dlneta
		Endif
	End Subroutine Compute_Diffusion_Coefs

	Subroutine Initialize_Transport_Coefficients()
		Call Allocate_Transport_Coefficients
		Call Initialize_Nu()							! Viscosity
		Call Initialize_Kappa()						! Thermal Diffusivity
		If (magnetism) Call Initialize_Eta()	! Magnetic Diffusivity
		Call Compute_Diffusion_Coefs
	End Subroutine Initialize_Transport_Coefficients


	Subroutine Allocate_Transport_Coefficients()
		Allocate(nu(1:N_r))
		Allocate(dlnu(1:N_r))
		Allocate(kappa(1:N_r))
		Allocate(dlnkappa(1:N_r))
		
		If (magnetism) Then
			Allocate(eta(1:N_R))
			Allocate(dlneta(1:N_R))
		Endif				
	End Subroutine Allocate_Transport_Coefficients

	Subroutine Initialize_Nu()
		Select Case(nu_type)
			Case(1)	! Constant nu
				nu(:) = nu_top
				dlnu(:) = 0.0d0
			Case(2)
				Call vary_with_density(nu,dlnu,nu_top, nu_power)
			!Case(3) -> Read from file
		End Select
	End Subroutine Initialize_Nu

	Subroutine Initialize_Kappa()
		Select Case(kappa_type)
			Case(1)	! Constant Kappa
				kappa(:) = kappa_top
				dlnkappa(:) = 0.0d0
			Case(2)
				Call vary_with_density(kappa,dlnkappa,kappa_top, kappa_power)
		End Select
	End Subroutine Initialize_Kappa

	Subroutine Initialize_Eta()
		Select Case(nu_type)
			Case(1)	! Constant Eta
				eta(:) = eta_top
				dlneta(:) = 0.0d0
			Case(2)
				Call vary_with_density(eta,dlneta,eta_top, eta_power)
			!Case(3) -> Read from file
		End Select
	End Subroutine Initialize_Eta


	Subroutine Vary_With_Density(coeff, dln, coeff_top, coeff_power)
		Real*8, Intent(InOut) :: coeff(:), dln(:)
		Real*8, Intent(In) :: coeff_top, coeff_power
		! Computes a transport coefficient and its logarithmic derivative
		! using a density-dependent form for the coefficient:
		!        coeff = coeff_top*(rho/rho_top)**coeff_power
		coeff = coeff_top*(ref%density/ref%density(1))**coeff_power
		dln = coeff_power*ref%dlnrho
	End Subroutine Vary_With_Density
End Module TransportCoefficients
