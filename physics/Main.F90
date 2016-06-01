!Rayleigh's main program
Program Main
    Use MakeDir
    Use Controls
	Use Fields
	Use Initial_Conditions
	Use Parallel_Framework
	Use ProblemSize
	Use Input
	Use Diagnostics_Interface, Only : Initialize_Diagnostics
	Use TestSuite
	Use Checkpointing
	Use Linear_Terms_Sphere
	Use Drive_Sphere, Only : Main_Loop_Sphere
	Use Timers
    Use Fourier_Transform, Only : Initialize_FFTs
    Use Benchmarking, Only : Initialize_Benchmarking, Benchmark_Input_Reset
	Implicit None
    
    Call Main_MPI_Init(global_rank)    !Initialize MPI

    Call Check_Run_Mode()   !This needs to be done before ever reading main input


	Call Main_Input()
	Call Benchmark_Input_Reset() ! Sets run parameters to benchmark parameters if benchmark_mode .ge. 0

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

        Call Initialize_Directory_Structure()

        Call Initialize_Benchmarking()

        Call Initialize_FFts()
		Call Initialize_Reference()

		Call Initialize_Transport_Coefficients()
        Call Initialize_Boundary_Conditions()


        Call Compute_Diffusion_Coefs()
		Call Write_Reference(ndrf)
		Call Initialize_Field_Structure()
		Call Initialize_Diagnostics()

		Call Full_Barrier()
		!Call Initialize_Benchmark_Equations()
		!Call Compute_Benchmark_Coefficients()
		!Call Set_Boundary_Conditions()
		Call Linear_Init() 
		Call Initialize_Checkpointing()
		Call Initialize_Fields()
		Call StopWatch(init_time)%increment() ! started in Init_Problemsize just after MPI is started up
	End Subroutine Main_Initialization

    Subroutine Initialize_Directory_Structure()
        Implicit None
        Integer :: ecode
        If (my_rank .eq. 0) Then
            Call Make_Directory(Trim(my_path)//'G_Avgs',ecode)
            Call Make_Directory(Trim(my_path)//'Shell_Avgs',ecode)
            Call Make_Directory(Trim(my_path)//'AZ_Avgs',ecode)
            Call Make_Directory(Trim(my_path)//'Shell_Slices',ecode)
            Call Make_Directory(Trim(my_path)//'Checkpoints',ecode)
            Call Make_Directory(Trim(my_path)//'Timings',ecode)
            Call Make_Directory(Trim(my_path)//'Spherical_3D',ecode)
            Call Make_Directory(Trim(my_path)//'Shell_Spectra',ecode)
            Call Make_Directory(Trim(my_path)//'Benchmark_Reports',ecode)
        Endif
    End Subroutine Initialize_Directory_Structure

	Subroutine Finalization()
        If (.not. test_mode) Then
         If (my_rank .eq. 0) Call stdout%finalize()
        Endif
		Call pfi%exit()
	End Subroutine Finalization
End Program Main
