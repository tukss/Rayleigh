Module NonDimensionalization
	Use Controls
	Use ReferenceState
	Use TransportCoefficients
	Use ProblemSize
    Use BoundaryConditions, Only : T_top, T_Bottom
	Implicit None
	Integer :: nd_type = 1	! Non-dimensionalization Type
	Integer :: nd_index = 1	! Radial index of reference state/transport parameters that non-dimensionalization is based off of
									! Use values at top of simulated domain by default
	Logical :: Dimensional = .false.	! Code runs non-dimensionally by default
												! Can run dimensionally for debugging purposes

	Logical :: use_dimensional_inputs = .true.

	!//////////////////////////
	!Placeholder Variables to make this compile for now...
	Real*8 :: Angular_Velocity = 1.0d0
	Logical :: Fix_Entropy_Top = .true., Fix_Entropy_Bottom = .true.
	Real*8 :: Entropy_Top = 0.0d0, Entropy_Bottom = 1000.0d0


    !////////////////////////////
    ! Dimensional-like terms

	!//////////////////////////
	! Nondimensional Parameters
	Real*8 :: Ra = 1.0d0
	Real*8 :: Ek = 1.0d0
	Real*8 :: Pr = 1.0d0
	Real*8 :: Pm = 1.0d0


	Real*8 :: nd_entropy, nd_gravity

	Real*8 :: nd_length, nd_mass, nd_time, nd_rho, nd_temperature, nd_nu, nd_kappa, nd_eta, nd_pressure

	Namelist /NonDimensionalization_Namelist/ Ra, Ek, Pr, Pm, nd_type, nd_index, dimensional, nd_entropy, angular_velocity, &
			use_dimensional_inputs, Entropy_Bottom
 

Contains

Subroutine NonDimensionalize
		If ( (.not. dimensional) .and. (use_dimensional_inputs) ) Then
			!
			Select Case(nd_type)
				Case(1)
					Call Standard_ND()
				Case(2)
				! something else
			End Select
		Endif


End Subroutine NonDimensionalize

Subroutine Standard_ND
	Implicit None
	Logical :: print_reference = .true.
	! Set standard length, time, and mass first
	nd_length = radius(1)-radius(N_r)
	!write(6,*)'nd length: ', nd_length
	nd_time = nd_length*nd_length/nu(nd_index)
	nd_rho = ref%density(nd_index)
	nd_mass = nd_rho*nd_length**3
	nd_mass = poly_mass
    nd_gravity = Gravitational_Constant*nd_mass/(nd_length**2)
	nd_temperature = ref%temperature(nd_index)
	nd_nu = nu(nd_index)
	nd_kappa = kappa(nd_index)
	If (fix_entropy_bottom .and. fix_entropy_top) Then
		If (nd_entropy .lt. 0) Then
			nd_entropy = entropy_top-entropy_bottom
		Endif
		! Otherwise, entropy should have been specified (necessary for dsdr =0 bottom simulations)
	Endif
	

	If (rotation) Then
		nd_pressure = angular_velocity*nd_rho*nd_nu
	Else
		nd_pressure = nd_rho*(nd_length/nd_time)**2
        Ek = 1.0d0  ! Technically Ek is infinity if Omega = 0
                    ! But the scaling of pressure by Ek in the nondimensionalized momentum eq
                    ! changes to no scaling at all whem Omega = 0, which is the same as Ek = 1
	Endif
	nd_entropy = T_bottom-T_top ! Entropy_Bottom-Entropy_Top
    T_bottom = T_bottom/nd_entropy      ! ugly - need to clean up the notation here
    T_top = T_top/nd_entropy            ! Mix of Boussinesq (T) and Anelastic (S)

	
	!Nondimensionalize the Transport Parameters
	nu = nu/nd_nu
	
	dlnu = dlnu*nd_length
	dlnkappa = dlnkappa*nd_length
	If (magnetism) Then
		nd_eta = eta(nd_index)
		Pm = nd_nu/nd_eta
		eta = eta/nd_eta
		dlneta = dlneta*nd_length
	Endif

	Ra = Gravitational_Constant*nd_mass*nd_entropy*nd_length/nd_nu/nd_kappa/Pressure_Specific_Heat
	Pr = nd_nu/nd_kappa
    kappa = kappa/nd_kappa/Pr
	Ek = nd_nu/Angular_Velocity/nd_length/nd_length

    ref%gravity_term_s = ref%gravity_term_s*(pressure_specific_heat/nd_gravity)
    ref%gravity_term_s = ref%gravity_term_s*Ra

	If ((my_rank .eq. 0) .and. (print_reference) ) Then
		Write(6,*)'Ra:   ', Ra
		Write(6,*)'Ek:   ', Ek
		Write(6,*)'Pr:   ', Pr
		Write(6,*)'el:   ', nd_length
		Write(6,*)'t :   ', nd_time
		Write(6,*)'rhoc: ', (ref%density(N_R/2)+ref%density(N_R/2+1))*0.5d0
		Write(6,*)'rho0: ', ref%density(1)
		Write(6,*)'Ti:   ', ref%temperature(N_R)
		Write(6,*)'Tc:   ', (ref%temperature(N_R/2)+ref%temperature(N_R/2+1))*0.5d0
		Write(6,*)'T0:   ', ref%temperature(1)
		Write(6,*)'Pi:   ', ref%pressure(N_R)
		Write(6,*)'Pc:   ', (ref%pressure(N_R/2)+ref%pressure(N_R/2+1))*0.5d0
		Write(6,*)'P0:   ', ref%pressure(1)
		Write(6,*)'Rgas: ', (ref%Gamma-1.0d0)*Pressure_Specific_Heat/ref%Gamma
	Endif
	!Nondimensionlize the Reference State	
	ref%density = ref%density/nd_rho
	ref%pressure = ref%pressure/nd_pressure
	ref%temperature = ref%temperature/nd_temperature
	ref%entropy = ref%entropy/nd_entropy
	ref%dsdr = ref%dsdr/nd_entropy*nd_length

	ref%dlnrho = ref%dlnrho*nd_length
	ref%d2lnrho = ref%d2lnrho*(nd_length**2)
	ref%dlnT = ref%dlnT*nd_length

	!Nondimensionalize the Derivatives
	Call Rescale_Grid_and_Derivatives(nd_length)
End Subroutine Standard_ND

End Module NonDimensionalization
