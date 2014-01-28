Module Input
	Use ProblemSize,  Only : problemsize_namelist
	Use Controls,     Only : controls_namelist
	Use Spherical_IO, Only : output_namelist
	Use BoundaryConditions, Only : boundary_conditions_namelist
	Use Initial_Conditions, Only : initial_conditions_namelist
	Use TestSuite, Only : test_namelist
	Implicit None
Contains

	Subroutine Main_Input()
		Implicit None
		Open(unit=20, file="main_input", status="old", position="rewind")
		Read(unit=20, nml=problemsize_namelist)
		Read(unit=20, nml=controls_namelist)
		Read(unit=20, nml=output_namelist)
		Read(unit=20, nml=boundary_conditions_namelist)
		Read(unit=20, nml=initial_conditions_namelist)
		Read(unit=20, nml=test_namelist)
		Close(20)
	End Subroutine Main_Input


End Module Input
