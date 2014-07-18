Program Main
	Use Fields
	Use Initial_Conditions
	Use Parallel_Framework
	Use ProblemSize
	Use Input
	!Use Physics
	Use Diagnostics
	Use TestSuite
	Use Checkpointing

	Use Linear_Terms_Sphere
	Use Drive_Sphere, Only : Main_Loop_Sphere
	Use Timers
	Implicit None

	Call Main_Input()
	

	If (test_mode) Then
		Call Init_ProblemSize()
		Call Test_Lib()
	Else
		Call Main_Initialization()

	  	!Call Main_Loop()
		Call Main_Loop_Sphere()
	Endif
	Call Finalization()
Contains
	Subroutine Main_Initialization()
		Implicit None
		Character*120 :: ndrf='reference_nd'
		Call Init_ProblemSize()
		Call Initialize_Reference()
		
		Call Initialize_Transport_Coefficients()
		Call NonDimensionalize()
		Call Write_Reference(ndrf)
		Call pfi%exit()
		STOP
		Call Initialize_Diagnostics()
		Call Initialize_Field_Structure()	! organization
		Call Full_Barrier()
		!Call Initialize_Benchmark_Equations()
		!Call Compute_Benchmark_Coefficients()
		!Call Set_Boundary_Conditions()
		Call Linear_Init() 
		Call Initialize_Checkpointing()
		Call Initialize_Fields()
		Call StopWatch(init_time)%increment() ! started in Init_Problemsize just after MPI is started up
	End Subroutine Main_Initialization


	Subroutine Finalization()
		Call pfi%exit()
	End Subroutine Finalization
End Program Main
