Module Input
	Use ProblemSize,  Only : problemsize_namelist, nprow, npcol, n_r,n_theta, npout
	Use Controls,     Only : controls_namelist, max_iterations, pad_alltoall
	Use Spherical_IO, Only : output_namelist
	Use BoundaryConditions, Only : boundary_conditions_namelist
	Use Initial_Conditions, Only : initial_conditions_namelist, alt_check
	Use TestSuite, Only : test_namelist
	!Use ArgCheck, Only : CheckArgs
	Implicit None
Contains

	Subroutine Main_Input()
		Implicit None
		! First read the main input file
		Open(unit=20, file="main_input", status="old", position="rewind")
		Read(unit=20, nml=problemsize_namelist)
		Read(unit=20, nml=controls_namelist)
		Read(unit=20, nml=output_namelist)
		Read(unit=20, nml=boundary_conditions_namelist)
		Read(unit=20, nml=initial_conditions_namelist)
		Read(unit=20, nml=test_namelist)
		Close(20)

		! Check the command line to see if any arguments were passed explicitly
		Call CheckArgs()

	End Subroutine Main_Input

	Subroutine CheckArgs()
			! Checks the command line for acceptable arguments.
			! Specified values overwrite namelist inputs.
			Implicit None
			Character*10 :: arg, arg2
			Integer :: i, itemp
			i = 1
			DO
	      	CALL get_command_argument(i, arg)
	         IF (LEN_TRIM(arg) == 0) EXIT

				arg2 = TRIM(AdjustL(arg))
				If (arg .eq. '-nprow') then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) nprow
				Endif
				If (arg .eq. '-npcol') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) npcol
				Endif
				If (arg .eq. '-npout') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) npout
				Endif
				If (arg .eq. '-niter') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) max_iterations
				Endif
				If (arg .eq. '-nr') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) n_r
				Endif
				If (arg .eq. '-ntheta') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) n_theta
				Endif
				If (arg .eq. '-pata') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) itemp
					if (itemp .eq. 1) then
						pad_alltoall = .true.
					else
						pad_alltoall = .false.
					endif
				Endif
				If (arg .eq. '-altc') Then
					CALL get_command_argument(i+1, arg)
					arg2 = TRIM(AdjustL(arg))
			      Read (arg2,*) itemp
					if (itemp .eq. 1) then
						alt_check = .true.
					else
						alt_check = .false.
					endif
				Endif
	      	i = i+1
				
	      END DO
	End Subroutine CheckArgs

End Module Input
