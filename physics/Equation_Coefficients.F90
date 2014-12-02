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
    Real*8, Allocatable :: dpdr_w_term(:), pressure_dwdr_term(:) ! For nondimensionalizing pressure
Contains

Subroutine Init_Equation_Coefficients
    Implicit None
    Integer :: i
    Real*8 :: amp, grav_r_ref
    Allocate(dpdr_w_term(1:N_R))
    dpdr_w_term = ref%density/Ekman_Number
    !dpdr_w_term(:) = ref%density  ! For now, I am nondimensionalizing differently than the benchmark
    Allocate(pressure_dwdr_term(1:N_R))
    !pressure_dwdr_term = -1.0d0*ref%density  ! This keeps 1/ek out when omega = 0
    pressure_dwdr_term(1:N_R) = -1.0d0/Ekman_Number*ref%density

    !The buoyancy term
    !    ---- Normally set as part of the reference state, but if we're nondimensional,
    !    we set it to Ra times (r/r_ref)^n
    If (.not. dimensional) Then
        amp = Rayleigh_Number/Prandtl_Number
        grav_r_ref = radius(1)
        Do i = 1, N_R
            ref%gravity_term_s(i) = amp*(radius(i)/grav_r_ref)**gravity_power
        Enddo
    Endif

    !////// Non-Magnetic Terms
    !CORIOLIS FORCE TERM (for Omega x v)
	If (rotation) Then
        If (dimensional) Then
            coriolis_term = 2.0d0*Angular_velocity
        Else
			coriolis_term = 2.0d0/Ekman_Number*Prandtl_Number
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
            Lorentz_Coefficient = Prandtl_Number/(Magnetic_Prandtl_Number*Ekman_Number)
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
