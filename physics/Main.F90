Program Main
    Use Controls
	Use Fields
	Use Initial_Conditions
	Use Parallel_Framework
	Use ProblemSize
	Use Input
	Use Diagnostics
	Use TestSuite
	Use Checkpointing
    Use Equation_Coefficients
	Use Linear_Terms_Sphere
	Use Drive_Sphere, Only : Main_Loop_Sphere
	Use Timers
    Use Fourier_Transform, Only : Initialize_FFTs
	Implicit None

    Call Main_MPI_Init(global_rank)    !Initialize MPI

    Call Check_Run_Mode()   !This needs to be done before ever reading main input


	Call Main_Input()
	

	If (test_mode) Then
		Call Init_ProblemSize()
		Call Test_Lib()
	Else
		Call Main_Initialization()


		Call Main_Loop_Sphere()
	Endif
	Call Finalization()
Contains
	Subroutine Main_Initialization()
		Implicit None
		Character*120 :: ndrf='reference_nd'
        Call Initialize_Controls()
        Call Set_Math_Constants()
		Call Init_ProblemSize()
        Call Initialize_Boundary_Conditions
        Call Initialize_FFts()
		Call Initialize_Reference()
		
		Call Initialize_Transport_Coefficients()
		Call NonDimensionalize()
        Call Compute_Diffusion_Coefs()
        Call Init_Equation_Coefficients()
		Call Write_Reference(ndrf)
		
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
        If (my_rank .eq. 0) Call stdout%finalize()
		Call pfi%exit()
	End Subroutine Finalization
End Program Main
