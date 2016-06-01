! This Module exists to manage the equation coefficients that can vary depending on
! whether a simulation is stratified vs. Boussinesq and dimensional vs. non-dimensional
! Coefficients might be constants real*8's or they might be arrays of real*8's dimensioned (1:n_r).
Module Equation_Coefficients
    Use ReferenceState
    Use TransportCoefficients
    Use ProblemSize
    Use Controls
    Use Math_Constants
    Use BoundaryConditions
    Implicit None
    !Real*8 :: coriolis_term
    Real*8 :: Lorentz_Coefficient

    Real*8, Allocatable :: ohmic_heating_coeff(:)   ! Need to adjust for nondimensional
    Real*8, Allocatable :: viscous_heating_coeff(:) ! Need to adjust for nondimensional
    !Real*8, Allocatable :: dpdr_w_term(:), pressure_dwdr_term(:) ! For nondimensionalizing pressure
Contains

Subroutine Init_Equation_Coefficients
    Implicit None
    Integer :: i, r
    Real*8 :: amp, grav_r_ref, dr, vhint, qadd2
    Real*8 :: lum_top, lum_bottom, sfactor, diff, r2dr, qadd, dsdr_mean
    Real*8, Allocatable :: tmp1d(:), tmp1d2(:)

    !might look into moving this somewhere else, but keepit here for now.
    If (heating_type .eq. 4) Then
        lum_top = -dtdr_top*kappa(1)*four_pi*(rmax**2)
        lum_top = lum_top*ref%density(1)*ref%temperature(1)
        lum_bottom = -dtdr_bottom*kappa(N_R)*four_pi*(rmin**2)
        lum_bottom = lum_bottom*ref%density(N_R)*ref%temperature(N_R)
        ref%heating = ref%heating*(lum_top-lum_bottom)


        !Now we build s_conductive (already set to zero)
        ! This needs to be reorganized, but do it here for now
        ! First, we build the indefinite integral of ref%heating *r^2
        Allocate(tmp1d(1:N_R),tmp1d2(1:N_R))
        tmp1d(:) =0.0d0
        tmp1d2 = r_squared*ref%density*ref%temperature*kappa
        sfactor = shell_volume/(four_pi)  !I think the one_third needs to be left out here



        tmp1d = (lum_top-lum_bottom)*(1.0d0/shell_volume)*(radius**3)/3.0d0
        
        vhint = tmp1d(1)*four_pi
        tmp1d = -tmp1d      
        !tmp1d is now r^2kappa rho T * dSdr_cond modulu an offset determined by the BCs
        diff = tmp1d2(N_R)*dTdr_bottom- tmp1d(N_R)
        tmp1d = tmp1d+diff
        tmp1d = tmp1d/tmp1d2  ! tmp1d is now dSdr_cond
        If (my_rank .eq. 0) Then
            Write(6,*)'rmax, rmin: ', rmax, rmin
            Write(6,*)'dsdr bottom is: ', tmp1d(N_R)
            Write(6,*)'dsdr top is: ', tmp1d(1)
            Write(6,*)'vhint is : ', vhint
            Write(6,*)'shell_volume: ', shell_volume
            Write(6,*)'lum_top: ', lum_top
            Write(6,*)'lum_bottom: ', lum_bottom
            Write(6,*)'ltop - lbottom = ', lum_top - lum_bottom
        Endif
        !Next, build s_conductive from dsdr_cond
        s_conductive(N_R) = 0.0d0
        Do r = N_R-1, 1,-1
            dsdr_mean = half*(tmp1d(r)+tmp1d(r+1))
            dr = radius(r)-radius(r+1)
            s_conductive(r) = s_conductive(r+1)+dsdr_mean*dr
        Enddo
        
        DeAllocate(tmp1d,tmp1d2)
    Endif


    !Allocate(dpdr_w_term(1:N_R))
    !dpdr_w_term = ref%density/Ekman_Number
    !dpdr_w_term(:) = ref%density  ! For now, I am nondimensionalizing differently than the benchmark
    !Allocate(pressure_dwdr_term(1:N_R))
    !pressure_dwdr_term = -1.0d0*ref%density  ! This keeps 1/ek out when omega = 0
    !pressure_dwdr_term(1:N_R) = -1.0d0/Ekman_Number*ref%density



    !////// Non-Magnetic Terms
    !CORIOLIS FORCE TERM (for Omega x v)
	!If (rotation) Then
    !    If (dimensional) Then
    !        coriolis_term = 2.0d0*Angular_velocity
    !    Else
	!		coriolis_term = 2.0d0/Ekman_Number*Prandtl_Number
    !    Endif
    !    If (NonDimensional_Anelastic) Then
    !        coriolis_term = 2.0d0
    !    Endif
    !
	!Endif

    ! Viscous Heating Coefficient
    ! Heating coeff is 2*nu/T_bar
    If (viscous_heating) Then
    	Allocate(viscous_heating_coeff(1:N_R))
        If (.not. dimensional) Then
        	viscous_heating_coeff(1:N_R) = 2.0d0 !<-------- Check
        Else
        	viscous_heating_coeff(1:N_R) = nu(1:N_R)*2.0d0/ref%temperature(1:N_R)
        Endif
        If (NonDimensional_Anelastic) Then
            viscous_heating_coeff(1:N_R) = nu(1:N_R)*2.0d0/ref%temperature(1:N_R)* &
                                           & Dissipation_Number/Modified_Rayleigh_Number
        Endif

    Endif


    !////// Magnetic Terms
	If (magnetism) Then
        !Lorentz Force Coefficient (for JXB)
        If (.not. dimensional) Then
            ! DOUBLE CHECK THIS LATER -- This should be consistent, but CHECK
            Lorentz_Coefficient = Prandtl_Number/(Magnetic_Prandtl_Number*Ekman_Number)
        Else
            Lorentz_Coefficient = 1.0d0/four_pi
        Endif
        ! Ohmic heating coefficient (multiplies {DelxB}^2)
        If (ohmic_heating) Then
            Allocate(ohmic_heating_coeff(1:N_R))
            ohmic_heating_coeff = lorentz_coefficient*eta/ref%density/ref%temperature
        Endif

	Endif

End Subroutine Init_Equation_Coefficients

End Module Equation_Coefficients
