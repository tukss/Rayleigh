! This Module exists to manage the equation coefficients that can vary depending on
! whether a simulation is stratified vs. Boussinesq and dimensional vs. non-dimensional
! Coefficients might be constants real*8's or they might be arrays of real*8's dimensioned (1:n_r).
Module Equation_Coefficients
    Use ReferenceState
    Use TransportCoefficients
    Use ProblemSize
    Use Controls
    Use Math_Constants
    Use NonDimensionalization
    Implicit None
    Real*8 :: coriolis_term
    Real*8 :: Lorentz_Coefficient
    Real*8 :: alf_const

    Real*8, Allocatable :: ohmic_heating_coeff(:)   ! Need to adjust for nondimensional
    Real*8, Allocatable :: viscous_heating_coeff(:) ! Need to adjust for nondimensional
Contains

Subroutine Init_Equation_Coefficients
    Implicit None

    !////// Non-Magnetic Terms
    !CORIOLIS FORCE TERM (for Omega x v)
	If (rotation) Then
        If (dimensional) Then
            coriolis_term = 2.0d0*Angular_velocity
        Else
			coriolis_term = 2.0d0/ek
        Endif
	Endif

    ! Viscous Heating Coefficient
    ! Heating coeff is 2*nu/T_bar
    If (viscous_heating) Then
    	Allocate(viscous_heating_coeff(1:N_R))
        If (.not. dimensional) Then
        	viscous_heating_coeff(1:N_R) = 2.0d0 !<-------- Check
        Else
        	viscous_heating_coeff(1:N_R) = nu(1:N_R)*2.0d0/ref%temperature(1:N_R)
        Endif
    Endif


    !////// Magnetic Terms
	If (magnetism) Then
        !Lorentz Force Coefficient (for JXB)
        If (.not. dimensional) Then
            Lorentz_Coefficient = 1.0d0/(Pm*ek)
        Else
            Lorentz_Coefficient = 1.0d0
        Endif
        ! Ohmic heating coefficient (for J^2)
        If (ohmic_heating) Then
            Allocate(ohmic_heating_coeff(1:N_R))
            If (.not. dimensional) Then
                ohmic_heating_coeff = 1.0d0     ! <--------- PROBABLY WRONG - CHECK!
            Else
                ohmic_heating_coeff = four_pi*eta/ref%density/ref%temperature
            Endif
        Endif

        !Alfven Speed Constant (for CFL)
        If (.not. dimensional) Then
            alf_const = 1.0d0
        Else
            alf_const = four_pi
        Endif
	Endif

End Subroutine Init_Equation_Coefficients

End Module Equation_Coefficients
