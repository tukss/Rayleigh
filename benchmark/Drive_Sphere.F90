Module Drive_Sphere
	Use ClockInfo
	Use Hybrid_Space_Sphere, Only : rlm_spacea, rlm_spaceb, hybrid_init
	Use Physical_Space_Sphere, Only : physical_space, coriolis_term
	Use Spectral_Space_Sphere, Only : post_solve, post_solve_cheby, advancetime, ctemp
	Use Checkpointing
	Use Controls
	Use Timers
	Use Fields
	Use Run_Parameters
	Use NonDimensionalization
Contains

	Subroutine Initialize_TimeStepping(iter)
		Implicit None
		Integer, Intent(In) :: iter
		Real*8 :: dr, tdiff

		dr = radius(1)-radius(2)  ! assume uniform grid for now
		tdiff = dr*dr				  ! Viscous diffusion time across one grid point		
		! This may seem stringent, but it is in-line with our max time step in ASH
		!max_time_step = tdiff *10.0d0 
		!max_time_step = max_time_step*4.0d0
		!min_time_step = max_time_step*1.0d-4

		If (iter .eq. 1) Then
			new_deltat   = max_time_step
			deltat   = 0.0d0
			old_deltat   = 0.0d0
			simulation_time = 0.0d0
		Else
			! We have restarted from a checkpoint
			! Change new_deltat and deltat appropriately
			new_deltat = checkpoint_newdt
			deltat = checkpoint_dt
			old_deltat = 0.0d0
			simulation_time = 0.0d0 !////////////// CHANGE THIS to WORK WITH CHECKPOINTING
		Endif
		new_timestep = .true.

	End Subroutine Initialize_TimeStepping

	Subroutine Main_Loop_Sphere()
		Implicit None
		Integer ::  last_iteration, first_iteration,i
		Real*8  :: captured_time			
		! We enter the main loop assuming that the solve has just been performed
		! and that the equation set structure's RHS contains our primary fields with 
		! radial dimension in-processor.
		! Care needs to be taken at init to ensure fields (W,Z,P,T) are stored
		! in the RHS (they are copied out upon entry into the loop).

		!Call Initialize_Timers()

		first_iteration = 1+checkpoint_iter ! checkpoint_iter is 0 by default
		last_iteration = first_iteration + max_iterations-1
		Call Initialize_TimeStepping(first_iteration)
		If (chebyshev .or. magnetism) Then
			! work structure for post_solve_cheby
			Call ctemp%init(field_count = wsfcount, config = 'p1b')
		Endif
		If (rotation) Then
            If (dimensional) Then
                coriolis_term = 2.0d0*Angular_velocity
            Else
    			coriolis_term = 2.0d0/ek
            Endif
		Endif
		If (magnetism) Then
			ovPm = 1.0d0/Pm
			ovPmek = 1.0d0/(Pm*ek)
		Endif

		Call Hybrid_Init()
		Call StopWatch(loop_time)%StartClock()
		Do iteration = first_iteration, last_iteration
			If (chebyshev) Then
				Call Post_Solve_Cheby()
			Else
				Call Post_Solve()	! Linear Solve configuration
			Endif

			If (my_rank .eq. 0) Write(6,*)'On iteration : ', iteration, deltat

			Call rlm_spacea()
			Call Physical_Space()
			Call rlm_spaceb()
			Call AdvanceTime()

			If (Mod(iteration,check_frequency) .eq. 0) Then
                Call Write_Checkpoint(wsp%p1b,iteration, deltat,new_deltat)                    
					 Call StopWatch(cwrite_time)%StartClock()
                !Call Write_Checkpoint_Alt(wsp%p1b,iteration, deltat,new_deltat)
					 Call StopWatch(cwrite_time)%Increment()
            Endif
		Enddo
		Call StopWatch(loop_time)%Increment()
		if (my_rank .eq. 0) Then
			Write(6,*)'//////////////////////////////////////////////'
			Write(6,*)' Elapsed time: ', StopWatch(loop_time)%elapsed
			Write(6,*)'  Column time: ', StopWatch(ctranspose_time)%elapsed
			Write(6,*)'     Row time: ', StopWatch(rtranspose_time)%elapsed
			Write(6,*)'Legendre time: ', StopWatch(legendre_time)%elapsed
			Write(6,*)'     FFT time: ', StopWatch(fft_time)%elapsed
			Write(6,*)'   Solve time: ', StopWatch(solve_time)%elapsed
			Write(6,*)'    rlma time: ', StopWatch(rlma_time)%elapsed
			Write(6,*)'    rlmb time: ', StopWatch(rlmb_time)%elapsed
			Write(6,*)'  pspace time: ', StopWatch(pspace_time)%elapsed
			Write(6,*)'  psolve time: ', StopWatch(psolve_time)%elapsed
			Write(6,*)'    dphi time: ', StopWatch(dphi_time)%elapsed
			captured_time = 0.0d0
			Do i = 2, 11
				captured_time = captured_time + StopWatch(i)%elapsed
			Enddo
			Write(6,*)'captured time: ', captured_time
			Write(6,*)'iter/sec     : ', max_iterations/StopWatch(loop_time)%elapsed
			Write(6,*)'//////////////////////////////////////////////'
			Write(6,*)'         sub times        '
			Write(6,*)'      nl time: ', StopWatch(nl_time)%elapsed
			Write(6,*)'    sdiv time: ', StopWatch(sdiv_time)%elapsed
			Write(6,*)'      ts time: ', StopWatch(ts_time)%elapsed
			Write(6,*)'      ar time: ', StopWatch(ar_time)%elapsed
		Endif
		Call Finalize_Timing(n_r,l_max,max_iterations)
	End Subroutine Main_Loop_Sphere

End Module Drive_Sphere
