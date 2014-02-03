Program Main
	Use Fields
	Use Initial_Conditions
	Use Parallel_Framework
	Use ProblemSize
	Use Input
	Use Physics
	Use Diagnostics
	Use TestSuite
	Use Checkpointing
	Implicit None

	Call Main_Input()
	

	If (test_mode) Then
		Call Init_ProblemSize()
		Call Test_Lib()
	Else
		Call Main_Initialization()

	  	Call Main_Loop()
	Endif
	Call Finalization()
Contains
	Subroutine Main_Initialization()
		Implicit None
		Call Init_ProblemSize()
		Call Initialize_Diagnostics()
		Call Initialize_Field_Structure()	! organization
		!Call Initialize_Benchmark_Equations()
		!Call Compute_Benchmark_Coefficients()
		!Call Set_Boundary_Conditions()
		Call Linear_Init() 
		Call Initialize_Checkpointing()
		Call Initialize_Fields()
		
	End Subroutine Main_Initialization


	Subroutine Finalization()
		Call pfi%exit()
	End Subroutine Finalization
End Program Main
