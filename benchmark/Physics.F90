#define RHSP wsp%p3b
#define FIELDSP wsp%p3a
#define DO_IDX Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define END_DO enddo; enddo; enddo
#define IDX k,r,t
Module Physics
	Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
	Use Parallel_Framework
	Use Controls
	Use ProblemSize
	Use Fourier_Transform
	Use Legendre_Polynomials, Only : p_lm_array, pi
	Use Legendre_Transforms, Only : Legendre_Transform 
	Use Spectral_Derivatives
	Use Finite_Difference, Only : d_by_dx
	Use Linear_Solve
	Use Fields
	Use Diagnostics
	Use BoundaryConditions
	Use General_MPI, Only : global_max
	Use Chebyshev_Polynomials, Only : cheby_to_spectral, cheby_from_spectral, d_by_dr_cp
	Use Checkpointing, Only : Write_Checkpoint, checkpoint_iter, checkpoint_dt, checkpoint_newdt
	Use Timers
	Implicit None

	Real*8 :: two_over_ek, ovPm, ovPmek

	Real*8 :: old_ab_factor = 1.0d0, new_ab_factor = 1.0d0
	Real*8 :: simulation_time

	Type(rmcontainer), Allocatable :: ftemp1(:), ftemp2(:)

	Logical :: new_timestep = .true.
	Real*8  :: new_deltat, deltat, old_deltat
	Real*8  :: cflmax = 0.1d0, cflmin = 0.05d0, min_dt_change = 0.1d0
	Real*8  :: max_time_step = 5.0d-4
	Real*8  :: min_time_step = 1.0d-13


	Integer :: iteration
	Character*8 :: t_ofmt = '(ES12.5)'	! For formatted timestep output
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
		Else
			! We have restarted from a checkpoint
			! Change new_deltat and deltat appropriately
			new_deltat = checkpoint_newdt
			deltat = checkpoint_dt
			old_deltat = 0.0d0
		Endif
		new_timestep = .true.

	End Subroutine Initialize_TimeStepping
	Subroutine Main_Loop()
		Implicit None
		Integer ::  last_iteration, first_iteration,i
		Real*8  :: captured_time			
		! We enter the main loop assuming that the solve has just been performed
		! and that the equation set structure's RHS contains our primary fields with 
		! radial dimension in-processor.
		! Care needs to be taken at init to ensure fields (W,Z,P,T) are stored
		! in the RHS (they are copied out upon entry into the loop).

		Call Initialize_Timers()

		first_iteration = 1+checkpoint_iter ! checkpoint_iter is 0 by default
		last_iteration = first_iteration + max_iterations-1
		Call Initialize_TimeStepping(first_iteration)
		If (chebyshev) Then
			! work structure for post_solve_cheby
			Call ctemp%init(field_count = wsfcount, config = 'p1b')
		Endif
		If (rotation) Then
			two_over_ek = 2.0d0/ek
		Endif
		If (magnetism) Then
			ovPm = 1.0d0/Pm
			ovPmek = 1.0d0/(Pm*ek)
		Endif


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
	End Subroutine Main_Loop
	Subroutine Post_Solve()	
		Implicit None
		Integer :: m, i
		Character*12 :: tstring, otstring
		! wsp%p1b is assumed to be allocated
		Call StopWatch(psolve_time)%startclock()
		Call wsp%construct('p1a')
		wsp%config = 'p1a'

		old_deltat = deltat
		If (new_timestep) Then
			!old_deltat = deltat
			deltat = new_deltat
			new_timestep = .false.
			If (my_rank .eq. 0) Then
				Write(otstring,t_ofmt)old_deltat
				Write(tstring,t_ofmt)deltat
				Write(6,*)'Timestep has changed from '//Trim(otstring)//' to '//Trim(tstring)//'.'
			Endif

			Call Reset_Linear_Equations()
		Endif
		if (iteration .eq. 1) then
				!Euler Step
				new_ab_factor = deltat
				old_ab_factor = 0.0d0
		else
				new_ab_factor = 0.5d0*deltat*(2 + deltat/old_deltat)
				old_ab_factor = -0.5d0*deltat**2/old_deltat
		endif
		!Write(6,*)'new, old', new_ab_factor, old_ab_factor
		wsp%p1b = wsp%p1b*old_ab_factor	
		
		!Copy each variable out of the RHS into the top part of the buffer

		Call Get_All_RHS(wsp%p1a)

	
		! Now take radial derivatives.  We can automate this further later.


		!///////////////////////////////////////////////////////////////////////////	
		!Load the W derivatives into the appropriate RHS's
		Call d_by_dx(wvar,d3wdr3,wsp%p1a,3)		! d3wdr3 will be overwritten by dwdr shortly
		Call Add_Derivative(peq,wvar,3,wsp%p1b,wsp%p1a,d3wdr3)

		Call d_by_dx(wvar,dwdr   ,wsp%p1a,1)		
		Call d_by_dx(wvar,d2wdr2 ,wsp%p1a,2)    

		Call Add_Derivative(peq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(peq,wvar,1,wsp%p1b,wsp%p1a,dwdr)

		Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(weq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		Call d_by_dx(dwdr,d2wdr2,wsp%p1a,1)	! cluge like in ASH  --- seems unnecessary though.  take out once all else works

		!//////////////////////////////
		!  P Terms
		Call d_by_dx(pvar,dpdr,wsp%p1a,1)	! dpdr will be overwritten by d2tdr2 shortly
		Call Add_Derivative(weq,pvar,1,wsp%p1b,wsp%p1a,dpdr)

		Call Add_Derivative(peq,pvar,0,wsp%p1b,wsp%p1a,pvar)


		!///////////////////////////////
		! T Terms
		Call d_by_dx(tvar,d2tdr2,wsp%p1a,2)	! d2tdr2 will be overwritten by dtdr shortly
		Call Add_Derivative(teq,tvar,2,wsp%p1b,wsp%p1a,d2tdr2)

		Call d_by_dx(tvar,dtdr,wsp%p1a,1)
		Call Add_Derivative(teq,tvar,1,wsp%p1b,wsp%p1a,dtdr)	

		Call Add_Derivative(teq,tvar,0, wsp%p1b,wsp%p1a,tvar)	
		Call Add_Derivative(weq,tvar,0, wsp%p1b,wsp%p1a,tvar)	! gravity

		! Convert temperature to temperature/r (will take derivatives of this for advection)
		
		Do m = 1, my_num_lm
			Do i = 1, 2
				wsp%p1a(:,i,m,tvar) = wsp%p1a(:,i,m,tvar)*one_over_r(:)
			Enddo
		Enddo


		!///////////////////////////////
		!  Z Terms
		Call d_by_dx(zvar,d2zdr2,wsp%p1a,2)		! 2nd derivative will be overwritten with dzdr
		Call Add_Derivative(zeq,zvar,2,wsp%p1b,wsp%p1a,d2zdr2)	

		Call d_by_dx(zvar,dzdr,wsp%p1a,1)	

		Call Add_Derivative(zeq,zvar,0,wsp%p1b,wsp%p1a,zvar)	


		If (magnetism) Then
			!//////////////
			! A-terms (Toroidal magnetic field)
			Call d_by_dx(avar,d2adr2,wsp%p1a,2)		! 2nd derivative will be overwritten with dadr
			Call Add_Derivative(aeq,avar,2,wsp%p1b,wsp%p1a,d2adr2)
			Call d_by_dx(avar,dadr,wsp%p1a,1)	

			Call Add_Derivative(aeq,avar,0,wsp%p1b,wsp%p1a,avar)

			!///////////////////
			! C-terms (Poloidal magnetic field)
			Call d_by_dx(cvar,d2cdr2,wsp%p1a,2)		
			Call Add_Derivative(ceq,cvar,2,wsp%p1b,wsp%p1a,d2cdr2)
			Call d_by_dx(cvar,dcdr,wsp%p1a,1)	
			Call Add_Derivative(ceq,cvar,0,wsp%p1b,wsp%p1a,cvar)

		Endif

		!Load the old ab array into the RHS
		Call Set_All_RHS(wsp%p1b)	! RHS now holds old_AB+CN factors

		Call wsp%deconstruct('p1b')
		Call StopWatch(psolve_time)%increment()

		Call StopWatch(ctranspose_time)%startclock()
		Call wsp%reform()	! move from p1a to s2a
		Call StopWatch(ctranspose_time)%increment()
	End Subroutine Post_Solve

	Subroutine Post_Solve_Cheby()	
		Implicit None
		Integer :: m, i
		Character*12 :: tstring, otstring
		! Eventually I will merge the two post-solve routines
		! For now, keeping them separate

		! wsp%p1b is assumed to be allocated
		Call StopWatch(psolve_time)%startclock()
		Call wsp%construct('p1a')
		wsp%config = 'p1a'

		old_deltat = deltat
		If (new_timestep) Then
			!old_deltat = deltat
			deltat = new_deltat
			new_timestep = .false.
			If (my_rank .eq. 0) Then
				Write(otstring,t_ofmt)old_deltat
				Write(tstring,t_ofmt)deltat
				Write(6,*)'Timestep has changed from '//Trim(otstring)//' to '//Trim(tstring)//'.'
			Endif

			Call Reset_Linear_Equations()
		Endif
		if (iteration .eq. 1) then
				!Euler Step
				new_ab_factor = deltat
				old_ab_factor = 0.0d0
		else
				new_ab_factor = 0.5d0*deltat*(2 + deltat/old_deltat)
				old_ab_factor = -0.5d0*deltat**2/old_deltat
		endif

		wsp%p1b = wsp%p1b*old_ab_factor	
		
		!Copy each variable out of the RHS into the top part of the buffer
		! These variables are in spectral space radially
		Call Get_All_RHS(wsp%p1a)
		! This is terribly inefficient, but I just want to test the stability of Chebyshev vs. FD for not..
		! We'll create a new buffer.  ctemp
		! Store all the permanent derivatives there - in c space
		ctemp%nf1a = 4
		ctemp%nf1b = 4
		Call ctemp%construct('p1a')
		! W..
		Call d_by_dr_cp(wvar,d3wdr3,wsp%p1a,3)
		ctemp%p1a(:,:,:,1) = wsp%p1a(:,:,:,d3wdr3)
		Call d_by_dr_cp(wvar,dwdr   ,wsp%p1a,1)		
		Call d_by_dr_cp(wvar,d2wdr2 ,wsp%p1a,2)
		! P....n
		Call d_by_dr_cp(pvar,dpdr,wsp%p1a,1)
		ctemp%p1a(:,:,:,2) = wsp%p1a(:,:,:,dpdr)
		! T
		Call d_by_dr_cp(tvar,d2tdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,3) = wsp%p1a(:,:,:,d2tdr2)
		Call d_by_dr_cp(tvar,dtdr,wsp%p1a,1)
		! Z..
		Call d_by_dr_cp(zvar,d2zdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,4) = wsp%p1a(:,:,:,d2zdr2)
		Call d_by_dr_cp(zvar,dzdr,wsp%p1a,1)

		!//////////////////////////////////////////////////////////////////////////
		! Now everything we need is in the wsp or ctemp buffer
		! The ctemp terms are those terms that do not leave this configuration
		! transform them now & add them to appropriate equations
		Call ctemp%construct('p1b')
		ctemp%p1a((2*N_r)/3+1:N_r,:,:,:) = 0.0d0	! de-alias
		Call Cheby_From_Spectral(ctemp%p1a,ctemp%p1b)
		Call Add_Derivative(peq,wvar,3,wsp%p1b,ctemp%p1b,1)
		Call Add_Derivative(weq,pvar,1,wsp%p1b,ctemp%p1b,2)
		Call Add_Derivative(teq,tvar,2,wsp%p1b,ctemp%p1b,3)
		Call Add_Derivative(zeq,zvar,2,wsp%p1b,ctemp%p1b,4)
		Call ctemp%deconstruct('p1a')
		Call ctemp%deconstruct('p1b')

		!//////////////////////////////////////////////
		!  Next, we reconstruct ctemp%p1a and copy wsp%p1a into it
		ctemp%nf1a = wsp%nf1a
		Call ctemp%construct('p1a')
		ctemp%p1a(:,:,:,:) = wsp%p1a(:,:,:,:)
		ctemp%p1a((2*N_r)/3+1:N_R,:,:,:) = 0.0d0 !De-Alias
		wsp%p1a(:,:,:,:) = 0.0d0	! Shouldn't need to do this, but just to be sure
		Call Cheby_From_Spectral(ctemp%p1a,wsp%p1a)
		Call ctemp%deconstruct('p1a')

		!/////////////////////////////////////////////////////////////////
		!  The rest of the code can remain unchanged
		!/////////////////////////////////////////////////////////////////	
		!Load the W derivatives into the appropriate RHS's  

		Call Add_Derivative(peq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(peq,wvar,1,wsp%p1b,wsp%p1a,dwdr)

		Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(weq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		

		!//////////////////////////////
		!  P Terms
		
		Call Add_Derivative(peq,pvar,0,wsp%p1b,wsp%p1a,pvar)


		!///////////////////////////////
		! T Terms
	
		Call Add_Derivative(teq,tvar,1,wsp%p1b,wsp%p1a,dtdr)	

		Call Add_Derivative(teq,tvar,0, wsp%p1b,wsp%p1a,tvar)	
		Call Add_Derivative(weq,tvar,0, wsp%p1b,wsp%p1a,tvar)	! gravity

		! Convert temperature to temperature/r (will take derivatives of this for advection)
		Do m = 1, my_num_lm
			Do i = 1, 2
				wsp%p1a(:,i,m,tvar) = wsp%p1a(:,i,m,tvar)/radius(:)
			Enddo
		Enddo


		!///////////////////////////////
		!  Z Terms


		Call Add_Derivative(zeq,zvar,0,wsp%p1b,wsp%p1a,zvar)	

		!Load the old ab array into the RHS
		Call Set_All_RHS(wsp%p1b)	! RHS now holds old_AB+CN factors

		Call wsp%deconstruct('p1b')
		Call StopWatch(psolve_time)%increment()

		Call StopWatch(ctranspose_time)%startclock()
		Call wsp%reform()	! move from p1a to s2a
		Call StopWatch(ctranspose_time)%increment()
	End Subroutine Post_Solve_Cheby


	!/////////////////////////////
	!  This routine calculates terms that involve
	!      theta derivatives and loads them into
	!		 the buffer.
	Subroutine rlm_spacea()
		Implicit None
		Integer :: mp
		Call StopWatch(rlma_time)%startclock()

				! Zero out l_max mode
		Do mp = my_mp%min, my_mp%max
			wsp%s2a(mp)%data(l_max,:) = 0.0d0
		Enddo

		
		! Allocate two work arrays
		Call Allocate_rlm_Field(ftemp1)
		Call Allocate_rlm_Field(ftemp2)

		Call Velocity_Components()	
		Call Velocity_Derivatives()
		Call d_by_dtheta(wsp%s2a,tvar,dtdt)

		If (magnetism) Call compute_BandJ()

		Call DeAllocate_rlm_Field(ftemp1)
		Call DeAllocate_rlm_Field(ftemp2)

		! Zero out l_max mode
		Do mp = my_mp%min, my_mp%max
			wsp%s2a(mp)%data(l_max,:) = 0.0d0
		Enddo

		Call StopWatch(rlma_time)%increment()

		!Legendre Transform and transpose the buffer
		Call wsp%construct('p2a')
		Call StopWatch(legendre_time)%startclock()
		Call Legendre_Transform(wsp%s2a,wsp%p2a)
		Call StopWatch(legendre_time)%increment()
		Call wsp%deconstruct('s2a')
		wsp%config = 'p2a'	

		Call StopWatch(rtranspose_time)%startclock()
		Call wsp%reform()	! We are now in p3a
		Call StopWatch(rtranspose_time)%increment()		

	End Subroutine rlm_spacea

	Subroutine Find_MyMinDT()
		Implicit None
		Real*8 :: ovt2, ovht2, ovrt2, maxt2
		Real*8 ::  maxt
		Integer :: r
		Call StopWatch(ts_time)%startclock()

		ovt2 = 0.0d0	! "over t squared"
		Do r = my_r%min, my_r%max
			ovht2 = Maxval(wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2) &
								*OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
			ovt2  = Max(ovt2, ovht2)
			ovrt2 = Maxval(wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)	! radial
			ovt2  = Max(ovt2,ovrt2)
		Enddo
		If (magnetism) Then
			! Check on alfven speed as well
			Do r = my_r%min, my_r%max
				ovht2 = Maxval(wsp%p3a(:,r,:,btheta)**2+wsp%p3a(:,r,:,bphi)**2) &
								*OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
				ovt2  = Max(ovt2, ovht2)
				ovrt2 = Maxval(wsp%p3a(:,r,:,br)**2)/(delta_r(r)**2)	! radial
				ovt2  = Max(ovt2,ovrt2)
			Enddo
		Endif

		Call wsp%set_mrv(ovt2)

		Call StopWatch(ts_time)%increment()
	End Subroutine Find_MyMinDT

	Subroutine Adjust_TimeStep()
		Implicit None
		Real*8 :: maxt2, maxt
		

		Call wsp%get_mrv(maxt2)
		if (maxt2 .gt. 0.0d0) Then
			maxt = 1.0d0/sqrt(maxt2)

			if (deltat .lt. maxt*cflmin) then
				! we can increase our timestep
				new_deltat = Min(cflmax*maxt,max_time_step)

			elseif (deltat .gt. (maxt*cflmax)) then
				new_deltat = cflmax*maxt 
				if (new_deltat .gt. deltat*(1.0d0-min_dt_change)) then
					! As much as possible, we would like to avoid
					! changing the timestep (slow process).  When we do change it,
					! make sure we give it a good bump.
					new_deltat = deltat*(1.0d0-min_dt_change)
				endif
			endif
		Endif
		If (new_deltat .ne. deltat) Then
			new_timestep = .true.
		Endif
		If (new_deltat .lt. min_time_step) Then
			If (my_rank .eq. 0) Write(6,*)'Time step became too small.'
			Call pfi%exit()
			Stop
		Endif

	End Subroutine Adjust_TimeStep

	Subroutine Compute_TimeStep()
		Implicit None
		Real*8 :: ovt2, ovht2, ovrt2, maxt2, dr
		Real*8 ::  maxt
		Integer :: r
		Call StopWatch(ts_time)%startclock()
		ovt2 = 0.0d0	! "over t squared"
		dr = radius(1)-radius(2)   ! Hard-coded uniform grid temporarily
		Do r = my_r%min, my_r%max
			ovht2 = Maxval(wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2)*OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
			ovt2 = Max(ovt2, ovht2)
			ovrt2 = Maxval(wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)	! radial
			ovt2 = Max(ovt2,ovrt2)
		Enddo

		If (magnetism) Then
			! Check on the Alfven speed as well
			Do r = my_r%min, my_r%max
				ovht2 = Maxval(wsp%p3a(:,r,:,btheta)**2+wsp%p3a(:,r,:,bphi)**2)*OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
				ovt2 = Max(ovt2, ovht2)
				ovrt2 = Maxval(wsp%p3a(:,r,:,br)**2)/(delta_r(r)**2)	! radial
				ovt2 = Max(ovt2,ovrt2)
			Enddo

		Endif

		Call StopWatch(ts_time)%increment()
		Call StopWatch(ar_time)%startclock()

		Call Global_Max(ovt2, maxt2, grp = pfi%gcomm)	! Wrapper to mpi all reduce

		Call StopWatch(ar_time)%increment()



		if (maxt2 .gt. 0.0d0) Then
			maxt = 1.0d0/sqrt(maxt2)
			!if (iteration .eq. 3) maxt = 0.5*min_time_step !max_time_step = max_time_step/2.0d0
			!Write(6,*)'max timestep = ', maxt, maxt2, ovt2
			if (deltat .lt. maxt*cflmin) then
				! we can increase our timestep
				new_deltat = Min(cflmax*maxt,max_time_step)

			elseif (deltat .gt. (maxt*cflmax)) then
				new_deltat = cflmax*maxt ! min(maxt,deltat) ! just to make sure we don't accidentally hit hte 'if' again
				if (new_deltat .gt. deltat*(1.0d0-min_dt_change)) then
					! As much as possible, we would like to avoid
					! changing the timestep (slow process).  When we do change it,
					! make sure we give it a good bump.
					new_deltat = deltat*(1.0d0-min_dt_change)
				endif
			endif
		Endif
		If (new_deltat .ne. deltat) Then
			new_timestep = .true.
		Endif
		If (new_deltat .lt. min_time_step) Then
			If (my_rank .eq. 0) Write(6,*)'Time step became too small.'
			Call pfi%exit()
			Stop
		Endif
	
	End Subroutine Compute_TimeStep


	Subroutine physical_space()
		Implicit None
	
		! We aren't quite in physical space yet.
		! 1st, get the phi derivatives
		Call StopWatch(dphi_time)%startclock()
		Call Phi_Derivatives()
		Call StopWatch(dphi_time)%increment()

		! Next perform the FFT
		Call StopWatch(fft_time)%startclock()
		Call fft_to_physical(wsp%p3a,rsc = .true.)
		Call StopWatch(fft_time)%increment()

		Call StopWatch(pspace_time)%startclock()
		! Convert all our terms of the form "sintheta var" to "var"
		Call StopWatch(sdiv_time)%startclock()
		Call sintheta_div(vtheta)	! sintheta vtheta to vtheta etc.
		Call sintheta_div(vphi)
		Call sintheta_div(dvtdr)
		Call sintheta_div(dvpdr)
		Call sintheta_div(dtdt)
		Call sintheta_div(dvrdt)
		Call sintheta_div(dvpdp)
		Call sintheta_div(dvtdp)

		If (magnetism) Then
			Call rsintheta_div(jtheta)
			Call rsintheta_div(jphi)
			Call rsintheta_div(Btheta)
			Call rsintheta_div(Bphi)
		Endif

		Call StopWatch(sdiv_time)%increment()

		!////////////////////////////////////////////////////////////////////////
		!This is a good spot to do some simple diagnostic output while we debug the code
		!since velocity components, Pressure, and Temperature are all 
		!in memory and in physical space at this point in time.
		Call ps_output(wsp%p3a, iteration)
		!////////////////////////////////////////////////////////////////////////

		If (test_reduce) Then
			Call Find_MyMinDT()	! Will piggyback on transposes
		Else
			Call compute_timestep()	! Will do it all at once
		Endif		
		
		! We are now ready to built the nonlinear terms
		Call wsp%construct('p3b')
		wsp%config = 'p3b'

		!................................
		!Nonlinear Advection
		Call StopWatch(nl_time)%startclock()

		Call Temperature_Advection()		
		Call Momentum_Advection_Radial()
		Call Momentum_Advection_Theta()
		Call Momentum_Advection_Phi()

		If (magnetism) Then
			Call Compute_EMF()
		Endif

		Call StopWatch(nl_time)%increment()
		!...........................

		Call wsp%deconstruct('p3a')

		Call StopWatch(pspace_time)%increment()


		Call StopWatch(fft_time)%startclock()
		Call fft_to_spectral(wsp%p3b, rsc = .true.)
		Call StopWatch(fft_time)%increment()

		

		Call StopWatch(rtranspose_time)%startclock()
		Call wsp%reform()	! Move to p2b
		Call StopWatch(rtranspose_time)%increment()
	End Subroutine Physical_Space


	Subroutine rlm_spaceb()
		Implicit None
		Integer :: m, mp, r,rmn,rmx,rind1,rind2,rmn1, rmn2, rmn3
		! Upon entry into this routine, we have the following quantities
		! Tvar : RHS for the T equation
		! Wvar : l(l+1)*RHS for the W equation
		! Pvar : r/sintheta * [u dot grad u]_theta
		! Zvar : r/sintheta * [u dot grad u]_phi

		! The RHS for T is ready to go
		! The W, Z and dWdr RHS's need a little work

		! Transform
		Call wsp%construct('s2b')

		Call StopWatch(legendre_time)%startclock()
		Call Legendre_Transform(wsp%p2b,wsp%s2b)	
		Call StopWatch(legendre_time)%increment()

		Call wsp%deconstruct('p2b')
		wsp%config = 's2b'
		Call StopWatch(rlmb_time)%startclock()

		

		! The NL RHS for W is r^2/(l(l+1)) * the NL RHS for Ur
		! We already have the r^2 taken care of.  Now for the l(l+1)
		rmn = (wvar-1)*tnr+1
		rmx = rmn+tnr-1	! striped real, then imaginary

		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				wsp%s2b(mp)%data(m:l_max,r) = wsp%s2b(mp)%data(m:l_max,r)*over_l_l_plus1(m:l_max)
			Enddo
		Enddo

		! Now for the Z RHS, formed from the radial component of the curl of u dot grad u

		Call Allocate_rlm_Field(ftemp1)
		Call Allocate_rlm_Field(ftemp2)

		Call d_by_sdtheta(wsp%s2b, zvar,ftemp1)	! need to be sure we have this indexing correct
		Call d_by_dphi(wsp%s2b,pvar,ftemp2)

		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = 1, tnr
				ftemp1(mp)%data(m:l_max,r) = ( ftemp2(mp)%data(m:l_max,r)- &
					& ftemp1(mp)%data(m:l_max,r) )*over_l_l_plus1(m:l_max)
			Enddo
		Enddo		!ftemp1 is now the Z RHS



		Call d_by_dphi(wsp%s2b,zvar,ftemp2)
		Call d_by_sdtheta(wsp%s2b,pvar,zvar)
		rmn1 = (pvar-1)*tnr
		rmn2 = (zvar-1)*tnr
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = 1, tnr
				rind1 = rmn1+r
				rind2 = rmn2+r
				wsp%s2b(mp)%data(m:l_max,rind1) = ( wsp%s2b(mp)%data(m:l_max,rind2)+ &
					& ftemp2(mp)%data(m:l_max,r) )*over_l_l_plus1(m:l_max)
			Enddo
		Enddo		! dwdr RHS (p equation) is now loaded

		rmn = (zvar-1)*tnr
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = 1, tnr
				wsp%s2b(mp)%data(m:l_max,rmn+r) = ftemp1(mp)%data(m:l_max,r)
				!BUFFER(m:l_max,rmn:rmx) = TEMP1(m:l_max,1:tnr)	-- might try setting up macros.  Easier to debug
			Enddo
		Enddo		! Z RHS is now loaded



		!The ell =0 w and p and z equations have zero RHS
		rmn1 = (pvar-1)*tnr+1
		rmn2 = (wvar-1)*tnr+1
		rmn3 = (zvar-1)*tnr+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			if (m .eq. 0) then
				wsp%s2b(mp)%data(0,rmn1:rmn1+tnr-1) = 0.0d0				
				wsp%s2b(mp)%data(0,rmn2:rmn2+tnr-1) = 0.0d0
				wsp%s2b(mp)%data(0,rmn3:rmn3+tnr-1) = 0.0d0
			endif
		Enddo


		If (magnetism) Call adjust_emf()

		Call DeAllocate_rlm_Field(ftemp1)
		Call DeAllocate_rlm_Field(ftemp2)


		! Zero out l_max mode
		Do mp = my_mp%min, my_mp%max
			wsp%s2b(mp)%data(l_max,:) = 0.0d0
		Enddo

		Call StopWatch(rlmb_time)%increment()

		Call StopWatch(ctranspose_time)%startclock()
		Call wsp%reform() ! move to the solve space
		Call StopWatch(ctranspose_time)%increment()

		If (test_reduce) Call Adjust_TimeStep()


	End Subroutine rlm_spaceb

	Subroutine Adjust_Emf()
		Implicit None
		Integer :: m, mp, r,rmn,rmx,rind1,rind2,rmn1, rmn2, rmn3, roff,rind


		! Now for the C RHS, formed from the radial component of the curl of the emf



		Call d_by_sdtheta(wsp%s2b, emfphi,ftemp1)	
		Call d_by_dphi(wsp%s2b,emftheta,ftemp2)

		rmn = (cvar-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				rind = r+roff
				wsp%s2b(mp)%data(m:l_max,r) = ( ftemp2(mp)%data(m:l_max,rind)- &
					& ftemp1(mp)%data(m:l_max,rind) )*over_l_l_plus1(m:l_max)
			Enddo
		Enddo		

		rmn = (emfphi-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				rind = r+roff
				wsp%s2b(mp)%data(m:l_max,r) = ( ftemp2(mp)%data(m:l_max,rind)+ &
					& ftemp1(mp)%data(m:l_max,rind) )*over_l_l_plus1(m:l_max)
			Enddo
		Enddo		

		! Ensure there is no ell=0 emf
		rmn1 = (emfr-1)    *tnr+1
		rmn2 = (emftheta-1)*tnr+1
		rmn3 = (emfphi-1)  *tnr+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			if (m .eq. 0) then
				wsp%s2b(mp)%data(0,rmn1:rmn1+tnr-1) = 0.0d0				
				wsp%s2b(mp)%data(0,rmn2:rmn2+tnr-1) = 0.0d0
				wsp%s2b(mp)%data(0,rmn3:rmn3+tnr-1) = 0.0d0
			endif
		Enddo
	End Subroutine Adjust_EMF

	Subroutine Temperature_Advectiono()
		Implicit None
		Integer :: t
		! Build  -u dot grad T

		wsp%p3b(:,:,:,tvar) = -wsp%p3a(:,:,:,vr)*wsp%p3a(:,:,:,dtdr)

		! dtdt is really (1/r) d temp d theta
		wsp%p3b(:,:,:,tvar) = wsp%p3b(:,:,:,tvar) - &
									 wsp%p3a(:,:,:,dtdt)*wsp%p3a(:,:,:,vtheta)

		!dtdp is really (1/r) d temp d phi
		Do t = my_theta%min, my_theta%max
			wsp%p3b(:,:,t,tvar) = wsp%p3b(:,:,t,tvar) - &
					& wsp%p3a(:,:,t,vphi)*wsp%p3a(:,:,t,dtdp)/sintheta(t)
		Enddo									
	End Subroutine Temperature_Advectiono

	Subroutine Temperature_Advection()
		Integer :: t,r,k

		!$OMP PARALLEL DO PRIVATE(t,r,k)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
				wsp%p3b(k,r,t,tvar) = -wsp%p3a(k,r,t,vr)*wsp%p3a(k,r,t,dtdr) &
									 - wsp%p3a(k,r,t,dtdt)*wsp%p3a(k,r,t,vtheta) &
									 - wsp%p3a(k,r,t,vphi)*wsp%p3a(k,r,t,dtdp)*csctheta(t)
				Enddo
			Enddo
		Enddo				
		!$OMP END PARALLEL DO

	End Subroutine Temperature_Advection

	Subroutine Momentum_Advection_Radialo()
		Implicit None
		Integer :: t,r

		! Build -radius^2 [u dot grad u]_r

		wsp%p3b(:,:,:,wvar) = -wsp%p3a(:,:,:,wvar)*wsp%p3a(:,:,:,dvrdr)	! vr dvrdr
	
		Do t = my_theta%min, my_theta%max	! vtheta/r dvrdtheta
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar) - &
						& wsp%p3a(:,r,t,vtheta)*wsp%p3a(:,r,t,dvrdt)/radius(r)
			Enddo
		Enddo	


		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar) - &	! vphi/sintheta/r dvrdphi
						& wsp%p3a(:,r,t,vphi)*wsp%p3a(:,r,t,dvrdp)/sintheta(t)/radius(r)
			Enddo
		Enddo	

		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar) + &	! centripetal force
						& (wsp%p3a(:,r,t,vphi)**2+wsp%p3a(:,r,t,vtheta)**2)/radius(r)
			Enddo
		Enddo	

		! Add Coriolis Terms if so desired
		If (rotation) Then
			! [- 2 z_hat cross u ]_r = 2 sintheta u_phi
			Do t = my_theta%min, my_theta%max
				Do r = my_r%min, my_r%max			
					wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar) + &
						& 2.0d0*sintheta(t)*wsp%p3a(:,r,t,vphi)/ek
				Enddo
			Enddo
		Endif


		! Multiply by radius squared so that we have NL RHS U_r * r^2 (= NL RHS W * l(l+1)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar)*r_squared(r)
			Enddo
		Enddo	
	End Subroutine Momentum_Advection_Radialo


	Subroutine Momentum_Advection_Radial()
		Implicit None
		Integer :: t,r,k

		! Build -radius^2 [u dot grad u]_r

		!$OMP PARALLEL DO PRIVATE(t,r,k)

		DO_IDX
			RHSP(IDX,wvar) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dvrdr)*r_squared(r) &
				- FIELDSP(IDX,vtheta) * ( FIELDSP(IDX,dvrdt)-FIELDSP(IDX,vtheta) )*radius(r)    &
				- FIELDSP(IDX,vphi)*(FIELDSP(IDX,dvrdp)*csctheta(t)-FIELDSP(IDX,vphi) )*radius(r)  
		END_DO
	
		!$OMP END PARALLEL DO

		If (magnetism) Then
			! Add r_squared [JxB]_r
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,wvar)= RHSP(IDX,wvar) +r_squared(r)*ovPmEk* &
					(FIELDSP(IDX,jtheta)*FIELDSP(IDX,bphi)-FIELDSP(IDX,jphi)*FIELDSP(IDX,btheta))
			END_DO
			!$OMP END PARALLEL DO
		Endif

		! Add Coriolis Terms if so desired
		If (rotation) Then
		!	! [- 2 z_hat cross u ]_r = 2 sintheta u_phi
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX			
				RHSP(IDX,wvar) = RHSP(IDX,wvar) + &
					& two_over_ek*sintheta(t)*FIELDSP(IDX,vphi)*R_squared(r)
			END_DO
			!$OMP END PARALLEL DO
		Endif
	
	End Subroutine Momentum_Advection_Radial

	Subroutine Compute_EMF()
		Implicit None
		Integer :: t,r,k

		! Build the emf

		!$OMP PARALLEL DO PRIVATE(t,r,k)

		DO_IDX	
			RHSP(IDX,emfr) = &
				  FIELDSP(IDX,vtheta) *  FIELDSP(IDX,bphi)  &
				- FIELDSP(IDX,vphi)     *  FIELDSP(IDX,btheta) 
		END_DO
	
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO PRIVATE(t,r,k)

		DO_IDX	
			RHSP(IDX,emftheta) = &
				- FIELDSP(IDX,vr) *  FIELDSP(IDX,bphi)  &
				+ FIELDSP(IDX,vphi)   *  FIELDSP(IDX,br)
		END_DO
	
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO PRIVATE(t,r,k)

		DO_IDX	
			RHSP(IDX,emfphi) = &
				  FIELDSP(IDX,vr)     *  FIELDSP(IDX,btheta)  &
				- FIELDSP(IDX,vtheta) *  FIELDSP(IDX,br)
		END_DO
	
		!$OMP END PARALLEL DO

		! We need to divide by r/sintheta before taking the derivatives in the next space
		!$OMP PARALLEL DO PRIVATE(t,r,k)

		DO_IDX	
			RHSP(IDX,emfphi) = RHSP(IDX,emfphi)*csctheta(t)*radius(r)
			RHSP(IDX,emftheta) = RHSP(IDX,emftheta)*csctheta(t)*radius(r)

		END_DO
	
		!$OMP END PARALLEL DO
	
	End Subroutine Compute_EMF


	Subroutine Momentum_Advection_Thetao()
		Implicit None
		Integer :: t, r
		! Build (radius/sintheta)[u dot grad u]_theta

		! First add all the terms that get multiplied by u_theta
		wsp%p3b(:,:,:,pvar) = wsp%p3a(:,:,:,dvrdr)	! d vr dr
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar)+ &	
						& + wsp%p3a(:,r,t,dvpdp) /sintheta(t)/radius(r) & ! vphi/sintheta/r dvrdphi
						& + wsp%p3a(:,r,t,vtheta)*cottheta(t)/radius(r) & !vtheta cot(theta)/r
						& + wsp%p3a(:,r,t,vr)/radius(r)					 		!ur/r
			Enddo
		Enddo	
		! multiply by -u_theta
		! NLRHS(:,:,:,vtheta) = NLRHS(:,:,:,vtheta)*FIELDS(:,:,:,vtheta)	Later it will look more like this
		wsp%p3b(:,:,:,pvar) = -wsp%p3b(:,:,:,pvar)*wsp%p3a(:,:,:,vtheta)


		wsp%p3b(:,:,:,pvar) = wsp%p3b(:,:,:,pvar)+wsp%p3a(:,:,:,vr)*wsp%p3a(:,:,:,dvtdr) ! vr dvthetadr

		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar) + &	! vphi/sintheta/r dvtheta dphi
						& wsp%p3a(:,r,t,vphi)*wsp%p3a(:,r,t,dvtdp)/sintheta(t)/radius(r)
			Enddo
		Enddo			

		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar) - &
						& wsp%p3a(:,r,t,vphi)*wsp%p3a(:,r,t,vphi)*cottheta(t)/radius(r)  ! vphi^2 cot(theta)/r
			Enddo
		Enddo	

		If (rotation) Then
			! Add - the coriolis term (part of -RHS of theta)
			! [2 z_hat cross u]_theta = -2 costheta u_phi
			Do t = my_theta%min, my_theta%max
				Do r = my_r%min, my_r%max
					wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar) - &
						& 2.0d0*costheta(t)*wsp%p3a(:,r,t,vphi)/ek
				Enddo
			Enddo
		Endif


		! At this point, we have [u dot grad u]_theta
		! Multiply by radius/sintheta so that we have r[u dot grad u]_theta/sintheta (getting ready for Z and dWdr RHS building)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar)*radius(r)/sintheta(t)
			Enddo
		Enddo	



	End Subroutine Momentum_Advection_Thetao

	Subroutine Momentum_Advection_Theta()
		Implicit None
		Integer :: t, r,k
		! Build (radius/sintheta)[u dot grad u]_theta

		! First add all the terms that get multiplied by u_theta
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,pvar) = wsp%p3a(IDX,dvrdr)       &	
				 + ( wsp%p3a(IDX,dvpdp)*csctheta(t)    & ! vphi/sintheta/r dvrdphi		!check this comment...
				 +   wsp%p3a(IDX,vtheta)*cottheta(t)   & !vtheta cot(theta)/r
				 +   wsp%p3a(IDX,vr)  ) *one_over_r(r)					 		!ur/r
		END_DO
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,pvar) = -RHSP(IDX,pvar)*wsp%p3a(IDX,vtheta) & ! multiply by -u_theta
				+ wsp%p3a(IDX,vr  )*wsp%p3a(IDX,dvtdr)					     & ! vr dvthetadr
				+ wsp%p3a(IDX,vphi)*( wsp%p3a(IDX,dvtdp)*csctheta(t) & ! vphi/sintheta/r dvtheta dphi
				- wsp%p3a(IDX,vphi )*cottheta(t) )*one_over_r(r)    ! vphi^2 cot(theta)/r

		END_DO
		!$OMP END PARALLEL DO

		If (magnetism) Then
			! Add -[JxB]_theta
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,pvar)= RHSP(IDX,pvar) &
					- ovPmEk*(FIELDSP(IDX,jphi)*FIELDSP(IDX,br)-FIELDSP(IDX,jr)*FIELDSP(IDX,bphi))
			END_DO
			!$OMP END PARALLEL DO
		Endif

		If (rotation) Then
			! Add - the coriolis term (part of -RHS of theta)
			! [2 z_hat cross u]_theta = -2 costheta u_phi

			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,pvar) = RHSP(IDX,pvar)- two_over_ek*costheta(t)*FIELDSP(IDX,vphi)
			END_DO
			!$OMP END PARALLEL DO
		Endif


		! At this point, we have [u dot grad u]_theta
		! Multiply by radius/sintheta so that we have r[u dot grad u]_theta/sintheta (getting ready for Z and dWdr RHS building)
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,pvar) = RHSP(IDX,pvar)*radius(r)*csctheta(t)
		END_DO
		!$OMP END PARALLEL DO



	End Subroutine Momentum_Advection_Theta



	Subroutine Momentum_Advection_Phio()
		Implicit None
		Integer :: t, r
		! Build (radius/sintheta)[u dot grad u]_phi

		! terms multiplied by u_theta
		!wsp%p3b(:,:,:,zvar) = wsp%p3a(:,:,:,vtheta)*(wsp%p3a(:,:,:,zvar)+wsp%p3a(:,:,:,dvtdp))

		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,zvar) = wsp%p3a(:,r,t,vtheta)*(wsp%p3a(:,r,t,zvar)+wsp%p3a(:,r,t,dvtdp)/sintheta(t)/radius(r))
			Enddo
		Enddo
		

		wsp%p3b(:,:,:,zvar) = wsp%p3b(:,:,:,zvar)+wsp%p3a(:,:,:,vr)*wsp%p3a(:,:,:,dvpdr)	! radial advection

		! terms multiplied by u_phi
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,zvar) = wsp%p3b(:,r,t,zvar) + wsp%p3a(:,r,t,vphi)* &
					& ( wsp%p3a(:,r,t,dvpdp)/sintheta(t) + wsp%p3a(:,r,t,vr))/radius(r)	
			Enddo
		Enddo

		If (rotation) Then
			! Add - Coriolis term (we are building -RHS of vphi)
			Do t = my_theta%min, my_theta%max
				Do r = my_r%min, my_r%max
					wsp%p3b(:,r,t,zvar) = wsp%p3b(:,r,t,zvar) + &
						& + 2.0d0*costheta(t)*wsp%p3a(:,r,t,vtheta)/ek + &
						& + 2.0d0*sintheta(t)*wsp%p3a(:,r,t,vr)/ek
				Enddo
			Enddo
		Endif

		! At this point, we have [u dot grad u]_phi
		! Multiply by radius/sintheta so that we have r[u dot grad u]_phi/sintheta (getting ready for Z and dWdr RHS building)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,zvar) = wsp%p3b(:,r,t,zvar)*radius(r)/sintheta(t)
			Enddo
		Enddo	

	End Subroutine Momentum_Advection_Phio

	Subroutine Momentum_Advection_Phi()
		Implicit None
		Integer :: t, r, k
		! Build (radius/sintheta)[u dot grad u]_phi

		! terms multiplied by u_theta
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,zvar) = FIELDSP(IDX,vtheta)*(FIELDSP(IDX,zvar)  & ! terms multiplied by u_theta
									+FIELDSP(IDX,dvtdp)*csctheta(t)*one_over_r(r)) &
				+FIELDSP(IDX,vr)*FIELDSP(IDX,dvpdr)	& ! radial advection
				+ FIELDSP(IDX,vphi) & ! terms multiplied by u_phi
				* ( FIELDSP(IDX,dvpdp)*csctheta(t) + FIELDSP(IDX,vr))*one_over_r(r)
		END_DO
		!$OMP END PARALLEL DO

		If (magnetism) Then
			! Add -[JxB]_phi
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,zvar)= RHSP(IDX,zvar) - &
					ovPmEk*(FIELDSP(IDX,jr)*FIELDSP(IDX,btheta)-FIELDSP(IDX,jtheta)*FIELDSP(IDX,br))
			END_DO
			!$OMP END PARALLEL DO
		Endif


		If (rotation) Then
			! Add - Coriolis term (we are building -RHS of vphi)
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,zvar) = RHSP(IDX,zvar)  					  &
					 + two_over_ek*costheta(t)*FIELDSP(IDX,vtheta) &
					 + two_over_ek*sintheta(t)*FIELDSP(IDX,vr)
			END_DO
			!OMP END PARALLEL DO
		Endif

		! At this point, we have [u dot grad u]_phi
		! Multiply by radius/sintheta so that we have r[u dot grad u]_phi/sintheta (getting ready for Z and dWdr RHS building)
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,zvar) = RHSP(IDX,zvar)*radius(r)*csctheta(t)
		END_DO
		!OMP END PARALLEL DO
	End Subroutine Momentum_Advection_Phi

	Subroutine Phi_Derivatives()
		Implicit None
		Integer :: r,t,k
		
		DO_IDX
			FIELDSP(IDX,pvar) = FIELDSP(IDX,tvar)! cluge to keep t
		END_DO


		Call d_by_dphi(wsp%p3a,vr,dvrdp)
		Call d_by_dphi(wsp%p3a,vtheta,dvtdp)
		Call d_by_dphi(wsp%p3a,vphi,dvpdp)
		Call d_by_dphi(wsp%p3a,tvar,dtdp)
	End Subroutine Phi_Derivatives


	Subroutine sintheta_div(ind)
		! Divide by sintheta
		Implicit None
		Integer, Intent(In) :: ind
		Integer :: t,r,k
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			FIELDSP(IDX,ind) = FIELDSP(IDX,ind)*csctheta(t)	
		END_DO
		!$OMP END PARALLEL DO
	End Subroutine sintheta_div

	Subroutine rsintheta_div(ind)
		Implicit None
		!divide by rsintheta
		Integer, Intent(In) :: ind
		Integer :: t,r,k
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			FIELDSP(IDX,ind) = FIELDSP(IDX,ind)*csctheta(t)*one_over_r(r)	
		END_DO
		!$OMP END PARALLEL DO
	End Subroutine rsintheta_div

	Subroutine Velocity_Components()
		Implicit None
		Integer :: l, m, mp, rmn,rmx, r, rind

		
		!Compute the velocity vield



		!vr	overwrites w	
		rmn = (vr-1)*tnr+1
		rmx = rmn+tnr-1	! striped real, then imaginary
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r-rmn+1
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,r)*ovrsq_repeated(rind)
			Enddo
		Enddo

		!We compute sintheta v_theta
		Call d_by_dtheta(wsp%s2a,dwdr,ftemp1)	 
		Call d_by_dphi(wsp%s2a,zvar,	ftemp2)	  			

		rmn = (vtheta-1)*tnr+1
		rmx = rmn+tnr-1 

		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)+ftemp2(mp)%data(:,:)
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)
			enddo
		Enddo

		!Now sintheta v_phi
		Call   d_by_dphi(wsp%s2a,dwdr,	ftemp1) 
		Call d_by_dtheta(wsp%s2a,zvar,ftemp2)
		rmn = (vphi-1)*tnr+1
		rmx = rmn+tnr-1 
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)-ftemp2(mp)%data(:,:)
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)
			enddo
		Enddo

	End Subroutine Velocity_Components


	Subroutine Velocity_Derivatives()
		Implicit None
		Integer :: r,rind,rind2,rmn,rmn1,rmx,rmx1, rmn2, rmx2
		Integer :: l, m, mp
		!/////////////////////////////////
		!sintheta dv theta dr 
		Call d_by_dtheta(wsp%s2a,d2wdr2,ftemp1)	! Store sintheta dwdtheta there for now.  We're going to use it a bit anyway.
		Call d_by_dphi(wsp%s2a,dzdr,	ftemp2)	   ! Will overwrite this with dTdtheta shortly			

		rmn1 = (vtheta-1)*tnr+1
		rmx1 = rmn1+tnr-1 
		rmn = (dvtdr-1)*tnr+1
		rmx = rmn+tnr-1

		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)+ftemp2(mp)%data(:,:)-wsp%s2a(mp)%data(:,rmn1:rmx1)
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)
			enddo
		Enddo

		!/////////////////////////////////
		!sinphi dv phi dr 
		Call d_by_dphi(wsp%s2a,d2wdr2,ftemp1)	! Store sintheta dwdtheta there for now.  We're going to use it a bit anyway.
		Call d_by_dtheta(wsp%s2a,dzdr,	ftemp2)	   ! Will overwrite this with dTdtheta shortly			

		rmn1 = (vphi-1)*tnr+1
		rmx1 = rmn1+tnr-1 
		rmn = (dvpdr-1)*tnr+1
		rmx = rmn+tnr-1	

		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)-ftemp2(mp)%data(:,:)-wsp%s2a(mp)%data(:,rmn1:rmx1)
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)
			enddo
		Enddo


		!/////////////////////////////////////////
		!dvrdr	overwrites dwdr	
		rmn = (dvrdr-1)*tnr
		rmn2 = (vr-1)*tnr
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = 1, tnr
			rind = r+rmn
			Do l = m, l_max
				wsp%s2a(mp)%data(l,rind) = l_l_plus1(l)*wsp%s2a(mp)%data(l,rind)*ovrsq_repeated(r)	! real part
			Enddo
			Enddo
			Do r = 1,tnr
				rind = r+rmn
				rind2 = r+rmn2
				wsp%s2a(mp)%data(:,rind) = wsp%s2a(mp)%data(:,rind)-2.0d0*wsp%s2a(mp)%data(:,rind2)*ovr_repeated(r)
			Enddo
		Enddo

		Call d_by_dtheta(wsp%s2a,vr,dvrdt)

		
		! Convert Z to -H_Laplacian Z
		rmn = (zvar-1)*tnr+1
		rmx = rmn+tnr-1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn,rmx
				rind = r -rmn+1
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,r)*ovrsq_repeated(rind)
			Enddo
		Enddo
	End Subroutine Velocity_Derivatives

	Subroutine Compute_BandJ()
		Implicit None
		Integer :: l, m, mp, rmn,rmx, r, rind, rmn2,rmx2, roff, rind2, roff2


		!/////////////// BR /////////////////////		
		!First convert C to Br  !! Br overwrites C
		rmn = (br-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r +roff !-rmn+1
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,r)*ovrsq_repeated(rind)
			Enddo
		Enddo        

		!////////////////// JR ///////////////////////////
		!Compute Jr (Jr does not overwrite any existing fields)
		rmn = (jr-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		rmn2 = (avar-1)*tnr+1       
		roff2 = -rmn+1+rmn2
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r +roff 
            rind2 = r+roff2
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,rind2)*ovrsq_repeated(rind)
			Enddo
		Enddo     

        !Convert d2cdr2 to d2cdr2-Br (br = cl(l+1)/r^2
        rmn = (d2cdr2-1)*tnr+1
        rmx = rmn+tnr-1
        rmn2 = (cvar-1)*tnr+1     
        roff = -rmn+1+rmn2
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = wsp%s2a(mp)%data(m:l_max,r)-wsp%s2a(mp)%data(m:l_max,rind)
			Enddo
		Enddo             

        ! Free up the dAdr space -- get its two angular derivatives
		Call d_by_dtheta(wsp%s2a,dadr,ftemp1)	 
		Call d_by_dphi(  wsp%s2a,dadr,ftemp2)

        !////////// J _PHI //////////////////////////
        ! overwrite d_a_dr with d_d_phi(d_a_dr)
        rmn = (dadr-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = ftemp2(mp)%data(m:l_max,rind)
			Enddo
		Enddo  

        !overwrite ftemp2 with d_d_theta (d2cdr2-br)
        Call d_by_dtheta(  wsp%s2a,d2cdr2,ftemp2)

        ! Add this term to d_d_phi(d_a_dr) to build rsintheta J_phi (overwrite dadr)
        rmn = (jphi-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = wsp%s2a(mp)%data(m:l_max,r)+ftemp2(mp)%data(m:l_max,rind)
			Enddo
        Enddo

        !/////////////J Theta ///////////////////////
        Call d_by_dphi(  wsp%s2a,d2cdr2,ftemp2)       !get phi derivative of d2cdr2-Br

        ! Combine with ftemp1 to build J_theta (overwrites d2cdr2)
        rmn = (jtheta-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = ftemp1(mp)%data(m:l_max,rind)-ftemp2(mp)%data(m:l_max,rind)
			Enddo
        Enddo


        

        !////////////B Theta
        ! Free up the A space -- get its two angular derivatives
		Call d_by_dtheta(wsp%s2a,avar,ftemp1)	 
		Call d_by_dphi(  wsp%s2a,avar,ftemp2)


        ! overwrite A with dA_d_phi
        rmn = (avar-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = ftemp2(mp)%data(m:l_max,rind)
			Enddo
		Enddo  

        !overwrite ftemp2 with d_d_theta (dcdr)
        Call d_by_dtheta(  wsp%s2a,dcdr,ftemp2)

        ! Add this term to dA_d_phi to build rsintheta B_theta
        rmn = (avar-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = wsp%s2a(mp)%data(m:l_max,r)+ftemp2(mp)%data(m:l_max,rind)
			Enddo
        Enddo

        !///////////// Bphi
        Call d_by_dphi(  wsp%s2a,dcdr,ftemp2)       !get phi derivative of dcdr

        ! Combine with ftemp1 to build rsintheta B_phi
        rmn = (dcdr-1)*tnr+1
        rmx = rmn+tnr-1      
        roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r+ roff 
				wsp%s2a(mp)%data(m:l_max,r) = ftemp2(mp)%data(m:l_max,rind)-ftemp1(mp)%data(m:l_max,rind)
			Enddo
        Enddo



	End Subroutine Compute_BandJ

	Subroutine Finalize_EMF()
		Implicit None
		Integer m, i
		! we need to take one last radial derivative and combine terms
		ctemp%nf1a = 1
		Call ctemp%construct('p1a')
		Do m = 1, my_num_lm
			Do i = 1, 2
				ctemp%p1a(:,i,m,1) = wsp%p1a(:,i,m,avar)
			Enddo
		Enddo
		Call d_by_dx(emfphi,avar,wsp%p1a,1)
		Do m = 1, my_num_lm
			Do i = 1, 2
				wsp%p1a(:,i,m,avar) = ctemp%p1a(:,i,m,1) + wsp%p1a(:,i,m,avar)
			Enddo
		Enddo
		Call ctemp%deconstruct('p1a')
	End Subroutine Finalize_EMF

	Subroutine AdvanceTime
		Implicit None
		! wsp will be in 'p1b' config
		! p1b contains the new adams bashforth term
		!Call print_max_spec2(pvar)

		if (.not. nonlinear) then
			wsp%p1b(:,:,:,:) = 0.0d0
		endif
		if (magnetism) then
			Call Finalize_EMF()
		endif
		Call Add_to_All_RHS(wsp%p1b,new_ab_factor)
		Call Fix_Boundary_Conditions()
		Call StopWatch(solve_time)%startclock()
		Call Implicit_Solve()
		Call StopWatch(solve_time)%increment()
		simulation_time = simulation_time+deltat
		! The righthand side of the equation set structure
		! Now contains the updated fields.
	End Subroutine AdvanceTime




	Subroutine Linear_Init()
		Implicit None

		Call Initialize_Benchmark_Equations()
		
	End Subroutine Linear_Init

	Subroutine Reset_Linear_Equations()
		Implicit None
		Real*8 :: rhs_factor, lhs_factor
		Call Reset_Equation_Coefficients()	! Zero out all coefficients
		lhs_factor = -deltat*alpha_implicit	! Crank Nicolson scheme - alpha = 0.5
		rhs_factor = deltat*(1.0d0-alpha_implicit)
		Call Set_Time_Factors(lhs_factor,rhs_factor)
		Call Compute_Benchmark_Coefficients()
		Call Set_Boundary_Conditions()
		Call LU_Decompose_Matrices()	! Last step after all matrices have been loaded
	End Subroutine Reset_Linear_Equations


	Subroutine Initialize_Benchmark_Equations
		Implicit None
		Integer :: neq, nvar,lp, l, nlinks
		Integer, Allocatable :: eq_links(:), var_links(:)
		If (magnetism) Then
			neq  = 6
			nvar = 6
		Else
			neq  = 4
			nvar = 4
		Endif
		If (chebyshev) Call Use_Chebyshev()	! Turns chebyshev mode to "on" for the linear solve
		Call Initialize_Equation_Set(neq,nvar,N_R,my_nl_lm, my_nm_lm)

		Do lp = 1, my_nl_lm
			l = my_lm_lval(lp)
			If (l .eq.0) Then

				nlinks = 3
				Allocate(eq_links(1:3))
				Allocate(var_links(1:3))

				eq_links(1) = weq
				eq_links(2) = peq
				eq_links(3) = teq

				var_links(1) = wvar
				var_links(2) = pvar
				var_links(3) = tvar

				! Not optimal (right now if variables are linked for one mode, they are linked for all).
				Call link_equations(eq_links, var_links,nlinks,lp)

				! ell = 0 pressure equation is hydrostatic balance
				Call Initialize_Equation_Coefficients(weq, wvar, 0,lp)		! identity matrix for W
				Call Initialize_Equation_Coefficients(peq, pvar, 1,lp)
				Call Initialize_Equation_Coefficients(peq, tvar, 0,lp)
				Call Initialize_Equation_Coefficients(teq, tvar, 2,lp)

				DeAllocate(eq_links)
				DeAllocate(var_links)
			Else
			! W equation		
			Call Initialize_Equation_Coefficients(weq,wvar,2,lp)
			Call Initialize_Equation_Coefficients(weq,pvar,1,lp)
			Call Initialize_Equation_Coefficients(weq,tvar,0,lp)


			! P equation
			Call Initialize_Equation_Coefficients(peq,wvar, 3,lp) 
			Call Initialize_Equation_Coefficients(peq,pvar, 0,lp)

			! T equation
			Call Initialize_Equation_Coefficients(teq,tvar, 2,lp)

			! Z equation
			Call Initialize_Equation_Coefficients(zeq,zvar, 2,lp) 

			If (magnetism) Then
				Call Initialize_Equation_Coefficients(ceq,cvar, 2,lp)
				Call Initialize_Equation_Coefficients(aeq,avar, 2,lp)
			Endif

			nlinks = 3
			Allocate(eq_links(1:3))
			Allocate(var_links(1:3))

			eq_links(1) = weq
			eq_links(2) = peq
			eq_links(3) = teq

			var_links(1) = wvar
			var_links(2) = pvar
			var_links(3) = tvar
			Call link_equations(eq_links, var_links,nlinks,lp)



			DeAllocate(eq_links)
			DeAllocate(var_links)
			Endif
		Enddo
		Call Finalize_Equations()	

		!============================================
		!  Next decide which derivatives need to be saved.
		!  This can be either for transposition later (i.e. nonlinear terms),
		!  or for checkpointing and output (such at with P).
		!  Saved derivatives will remain in memory until explicitly deallocated
		!  by the user.  The Implicit Module simply places them into memory.  It 
		!  does not deallocate them.

		Call Set_Deriv_Save(wvar,2)
		Call Set_Deriv_Save(pvar,0)
		Call Set_Deriv_Save(tvar,1)
		Call Set_Deriv_Save(zvar,1)
		If (magnetism) Then
			Call Set_Deriv_Save(avar,1)
			Call Set_Deriv_Save(cvar,2)
		Endif
	End Subroutine Initialize_Benchmark_Equations

	Subroutine Compute_Benchmark_Coefficients()
		Implicit None

		Real*8, Allocatable :: H_Laplacian(:), amp(:)
		Integer :: l, lp
		!rmin_norm

		Allocate(amp(1:N_R))
		Allocate(H_Laplacian(1:N_R))
		Do lp = 1, my_nl_lm

			l = my_lm_lval(lp)		

			H_Laplacian = - l_l_plus1(l) * OneOverRSquared
			If (l .eq. 0) Then
				!====================================================
				!			Temperature Equation
				! T only
				amp = 1.0d0
				Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)	! Time independent part

				amp = 2.0d0/radius/Pr
				Call add_implicit_term(teq,tvar, 1, amp,lp)
				amp = 1.0d0/Pr
				Call add_implicit_term(teq,tvar, 2, amp,lp)


				!=======================================
				!   Hydrostatic balance
				! 	 This part of the equation is static (i.e. not time-evolving)

				! 	 t
				amp = -Ra/ek! *(radius/r_outer)
				Call add_implicit_term(peq, tvar, 0, amp,lp, static = .true.)			! Gravity	--- Need LHS_Only Flag

				amp = 1.0d0
				Call add_implicit_term(peq, pvar, 1, amp,lp, static = .true.)			! dPdr     --- Here too


				! Pressure equation is replaced with a statement that the ell=0 W is zero
				amp = 1.0d0
				Call add_implicit_term(weq, wvar, 0, amp,lp, static = .true.)			! --- any maybe here

				!  If band solve, do redefinition here
			Else

				!==================================================
				!				Radial Momentum Equation
				
				! Temperature
				amp = -(Ra/Ek)*( (radius/r_outer)**gpower )
				amp = amp/H_Laplacian
				Call add_implicit_term(weq, tvar, 0, amp,lp)			! Gravity

				! Pressure
				amp = 1.0d0/(Ek*H_Laplacian)		! dPdr
				Call add_implicit_term(weq,pvar, 1, amp,lp)


				! W
				amp = 1.0d0
				Call add_implicit_term(weq,wvar, 0, amp,lp,static = .true.)	! This term does not a get a dt factor

				amp = H_Laplacian		! Diffusion
				Call add_implicit_term(weq,wvar, 0, amp,lp)
				amp = 1.0d0
				Call add_implicit_term(weq,wvar, 2, amp,lp)

				!==================================================
				!				Pressure (dWdr) Equation
				
				! Pressure
				amp = -(1.0d0)/Ek	
				Call add_implicit_term(peq,pvar, 0, amp,lp)

				! W
				amp = 1.0d0
				Call add_implicit_term(peq,wvar, 1, amp,lp, static = .true.)	! Time independent term
				amp =-H_Laplacian*2.0d0/radius	
				Call add_implicit_term(peq,wvar, 0, amp,lp)
				amp = H_Laplacian
				Call add_implicit_term(peq,wvar, 1, amp,lp)
				amp = 1.0d0
				Call add_implicit_term(peq,wvar, 3, amp,lp)


				!====================================================
				!			Temperature Equation

				! T 
				amp = 1.0d0
				Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)		! Time independent term
				amp = H_Laplacian/Pr		! Diffusion
				Call add_implicit_term(teq,tvar, 0, amp,lp)
				amp = 2.0d0/radius/Pr
				Call add_implicit_term(teq,tvar, 1, amp,lp)
				amp = 1.0d0/Pr
				Call add_implicit_term(teq,tvar, 2, amp,lp)


				!=====================================================
				!	Z Equation
				amp = 1.0d0
				Call add_implicit_term(zeq,zvar, 0, amp,lp, static = .true.)	! Time-independent piece
				amp = H_Laplacian
				Call add_implicit_term(zeq,zvar, 0, amp,lp)				
				amp = 1.0d0
				Call add_implicit_term(zeq,zvar, 2, amp,lp)				

				If (magnetism) Then
					!=========================================
					!  Btor Equation
					amp = 1.0d0
					Call add_implicit_term(aeq,avar, 0, amp,lp, static = .true.)	! Time-independent piece
					amp = H_Laplacian/Pm
					Call add_implicit_term(aeq,avar, 0, amp,lp)					
					amp = 1.0d0/Pm
					Call add_implicit_term(aeq,avar, 2, amp,lp)	

					!=========================================
					!  Bpol Equation
					amp = 1.0d0
					Call add_implicit_term(ceq,cvar, 0, amp,lp, static = .true.)	! Time-independent piece
					amp = H_Laplacian/Pm
					Call add_implicit_term(ceq,cvar, 0, amp,lp)
					amp = 1.0d0/Pm
					Call add_implicit_term(ceq,cvar, 2, amp,lp)					
				Endif
			
				! If band solve, do the redefinition of the matrix here

			Endif
		Enddo
		DeAllocate(amp)
		DeAllocate(H_Laplacian)
	End Subroutine Compute_Benchmark_Coefficients

	Subroutine Set_Boundary_Conditions
		Implicit None
		Real*8 :: samp,one
		Integer :: l, r,lp,  dorder
		one = 1.0d0

		Do lp = 1, my_nl_lm
			l = my_lm_lval(lp)

			If (l .eq. 0) Then
				Call Clear_Row(peq,lp,1)			! Pressure only has one boundary condition
				Call Clear_Row(teq,lp,1)
				Call Clear_Row(teq,lp,N_R)

				! Temperature Boundary Conditions (T fixed bottom and top)
				r = 1
				Call Load_BC(lp,r,teq,tvar,one,0)	!upper boundary
				r = N_R
				Call Load_BC(lp,r,teq,tvar,one,0)	! lower boundary


				! The ell=0 pressure is really a diagnostic of the system.
				! It doesn't drive anything.  The simplist boundary condition
				! is to enforce a pressure node at the top.
				r = 1	
				Call Load_BC(lp,r,peq,pvar,one,0)


			Else

				!*******************************************************
				!		Clear the boundary rows
				Call Clear_Row(weq,lp,1)
				Call Clear_Row(weq,lp,N_R)
				Call Clear_Row(peq,lp,1)
				Call Clear_Row(peq,lp,N_R)			
				Call Clear_Row(teq,lp,1)
				Call Clear_Row(teq,lp,N_R)
				Call Clear_Row(zeq,lp,1)
				Call Clear_Row(zeq,lp,N_R)


				!*******************************************************
				! Entropy Boundary Conditions

				! Temperature Boundary Conditions (T fixed bottom and top)
				r = 1
				Call Load_BC(lp,r,teq,tvar,one,0)	!upper boundary
				r = N_R
				Call Load_BC(lp,r,teq,tvar,one,0)	! lower boundary
				



				!************************************************************
				! Velocity Boundary Conditions
		
				! Impenetrable top and bottom
				! W vanishes at the boundaries
				r = 1
				Call Load_BC(lp,r,weq,wvar,one,0)
				r = N_R
				Call Load_BC(lp,r,weq,wvar,one,0)

		
				! No Slip Top and Bottom
				! Z and dWdr vanish at the boundaries
				r = 1
				Call Load_BC(lp,r,peq,wvar,one,1)
				Call Load_BC(lp,r,zeq,zvar,one,0)
				r = N_R
				Call Load_BC(lp,r,peq,wvar,one,1)
				Call Load_BC(lp,r,zeq,zvar,one,0)


				!*******************************************************
				!		Magnetic Boundary Conditions

				If (Magnetism) Then
					!  Clear the boundary rows
					Call Clear_Row(ceq,lp,1)
					Call Clear_Row(ceq,lp,N_R)
					Call Clear_Row(aeq,lp,1)
					Call Clear_Row(aeq,lp,N_R)


					! Match to a potential field at top and bottom
					! Btor = 0 at top and bottom
					r = 1
					Call Load_BC(lp,r,aeq,avar,one,0)
					r = N_R
					Call Load_BC(lp,r,aeq,avar,one,0)

					! dBpol/dr+ell*Bpol/r = 0 at outer boundary
					r = 1
					Call Load_BC(lp,r,ceq,cvar,one,1)
					samp = my_lm_lval(lp)*one_over_r(r)
					Call Load_BC(lp,r,ceq,cvar,samp,0)

					! dBpol/dr-ell(ell+1)*Bpol/r = 0 at inner boundary
					r = N_R
					Call Load_BC(lp,r,ceq,cvar,one,1)	
					samp = - l*(l+1)*One_Over_R(r)
					Call Load_BC(lp,r,ceq,cvar,samp,0)	

				Endif	! Magnetism


			Endif ! l = 0 or not
		Enddo


	End Subroutine Set_Boundary_Conditions

	Subroutine Fix_Boundary_Conditions()
		Implicit None
		Integer :: l, indx, ii,lp
      ! start applying the boundary and continuity conditions by setting 
      ! the appropriate right hand sides.

		! This is ugly, and I have no idea how to make this pretty.
		! Might zero these by default and then call a routine for exceptions
		! such as fixed entropy top.
      ii = 2*N_R

      
		indx = 1

		Do lp = 1, my_nl_lm
			n_m = my_nm_lm(lp)-1	! really n_m, but the indexing below is from the old implicit solv
			l = my_lm_lval(lp)
 
			If (l /= 0) Then

				equation_set(1,weq)%RHS(1+2*N_R  ,:,indx:indx+n_m)    = zero
				equation_set(1,weq)%RHS(N_R+2*N_R,:,indx:indx+n_m)    = zero
	
            equation_set(1,weq)%RHS(1    ,:,indx:indx+n_m) = Zero
            equation_set(1,weq)%RHS(N_R  ,:,indx:indx+n_m) = Zero


				equation_set(1,weq)%RHS(1+N_R,:,indx:indx+n_m) = Zero
				equation_set(1,weq)%RHS(2*N_R,:,indx:indx+n_m) = Zero

            equation_set(1,zeq)%RHS(1  ,:,indx:indx+n_m)    = Zero
            equation_set(1,zeq)%RHS(N_R,:,indx:indx+n_m)    = Zero



				If (Magnetism) Then
					equation_set(1,ceq)%RHS(1,:,indx:indx+n_m) = Zero
					equation_set(1,ceq)%RHS(N_R,:,indx:indx+n_m) = Zero
 
					equation_set(1,aeq)%RHS(1,:,indx:indx+n_m) = Zero
					equation_set(1,aeq)%RHS(N_R,:,indx:indx+n_m) = Zero

          
				Endif

			Else

				equation_set(1,weq)%RHS(:,2,indx)   = Zero ! no imaginary part for any ell=0 equations
				equation_set(1,zeq)%RHS(:,:,indx)   = Zero	! no ell = 0 z_equation

				equation_set(1,weq)%rhs(1:N_R,:,indx) = zero				! ell =0 W is zero
				equation_set(1,weq)%rhs(N_R+1,1,indx) = zero	! Pressure node



				!Top temperature (in spectral space, but BC's specified in physical space
				!    so multiply by sqrt4pi)
				equation_set(1,weq)%RHS(1+ii,1,indx)   = T_Top*sqrt(4.0D0*Pi)
				!Bottom Temperature
				equation_set(1,weq)%RHS(N_R+ii,1,indx) = T_Bottom*sqrt(4.0D0*Pi)

            
         Endif



         indx = indx + n_m + 1
      Enddo
    End Subroutine Fix_Boundary_Conditions


	Subroutine Allocate_rlm_Field(arr)
		Implicit None
		Type(rmcontainer), Intent(InOut), Allocatable :: arr(:)
		Integer :: mp,m


		Allocate(arr(my_mp%min:my_mp%max))
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Allocate(arr(mp)%data(m:l_max,1:tnr))
		Enddo
	End Subroutine Allocate_rlm_Field

	Subroutine DeAllocate_rlm_Field(arr)
		Implicit None
		Type(rmcontainer), Intent(InOut), Allocatable :: arr(:)
		Integer :: mp
		Do mp = my_mp%min, my_mp%max
			DeAllocate(arr(mp)%data)
		Enddo
		DeAllocate(arr)
	End Subroutine DeAllocate_rlm_Field


	Subroutine print_max_rlm(lcheck,mcheck,varind)
		! Simple debugging routine
		Implicit None
		Integer, Intent(In) :: mcheck, lcheck,varind
		Integer :: m, l, mp
		Real*8 :: mxv
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			If (m .eq. mcheck) Then
				rmin = (varind-1)*2*n_r+1
				rmax = rmin+2*n_r-1
				if (wsp%config .eq. 's2a') then
					mxv = maxval(abs(wsp%s2a(mp)%data(lcheck,rmin:rmax)))
					write(6,*)'maxval s2a: ', mxv, lcheck, mcheck,varind
				else
					mxv = maxval(wsp%s2b(mp)%data(lcheck,rmin:rmax))
					write(6,*)'maxval s2b: ', mxv, lcheck, mcheck,varind
				endif				
			Endif
		Enddo
	End Subroutine print_max_rlm



	Subroutine print_max_spec(lcheck,mcheck,varind)
		! Simple debugging routine
		! Prints malval of field in spectral space config at lcheck mcheck mode
		Implicit None
		Integer, Intent(In) :: mcheck, lcheck,varind
		Integer :: m, l, mp, lp
		Real*8 :: mxv

		Do lp = my_lm_min, my_lm_max
			l = l_lm_values(lp)
			m = m_lm_values(lp)
			If ( (l .eq. lcheck) .and. (m .eq. mcheck) ) Then
				mxv = maxval(abs(wsp%p1a(:,:,lp,varind)))
				write(6,*)'maxval p1a: ', mxv, lcheck, mcheck,varind
			Endif
		Enddo
	End Subroutine print_max_spec

	Subroutine zero_mode_spec(lcheck,mcheck,varind)
		! Simple debugging routine
		! Prints malval of field in spectral space config at lcheck mcheck mode
		Implicit None
		Integer, Intent(In) :: mcheck, lcheck,varind
		Integer :: m, l, mp, lp
		Real*8 :: mxv

		Do lp = my_lm_min, my_lm_max
			l = l_lm_values(lp)
			m = m_lm_values(lp)
			If ( (l .eq. lcheck) .and. (m .eq. mcheck) ) Then
				wsp%p1a(:,:,lp,varind) = 0.0d0
				write(6,*)'zeroed: ', lcheck, mcheck,varind
			Endif
		Enddo
	End Subroutine zero_mode_spec


	Subroutine print_max_spec2(varind)
		! Simple debugging routine
		! Prints max of field in spectral space at all ell's and m's
		Implicit None
		Integer, Intent(In) :: varind
		Real*8 :: mxv

		if (wsp%config .eq. 'p1a') Then
			mxv = maxval(abs(wsp%p1a(:,:,:,varind)))
			write(6,*)'maxval p1a: ', mxv, varind
		else
			mxv = maxval(abs(wsp%p1b(:,:,:,varind)))
			write(6,*)'maxval p1b: ', mxv, varind			
		endif

	End Subroutine print_max_spec2




	Subroutine find_power_rlm(varind)
		! Simple debugging routine
		Implicit None
		Integer, Intent(In) :: varind
		Integer :: m, l, mp
		Real*8 :: mxv
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			rmin = (varind-1)*2*n_r+1
			rmax = rmin+2*n_r-1
			Do l = m, l_max
				if (wsp%config .eq. 's2a') then
					mxv = maxval(abs(wsp%s2a(mp)%data(l,:)))
					if (mxv .gt. 0.0d0) then
						write(6,*)'s2a power at: ', mxv, m, l, varind
					endif
				else
					!mxv = maxval(wsp%s2b(mp)%data(lcheck,rmin:rmax))
					!write(6,*)'maxval s2b: ', mxv, lcheck, mcheck,varind
				endif				
			Enddo
		Enddo
	End Subroutine find_power_rlm

End Module Physics
