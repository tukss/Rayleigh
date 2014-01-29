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
	Implicit None

	Real*8 :: old_ab_factor = 1.0d0, new_ab_factor = 1.0d0
	Real*8 :: simulation_time

	Type(rmcontainer), Allocatable :: ftemp1(:), ftemp2(:)

	Logical :: new_timestep = .true., already_at_max
	Real*8  :: new_deltat, deltat, old_deltat
	Real*8  :: cflmax = 0.1d0, cflmin = 0.05d0, min_dt_change = 0.1d0
	Real*8  :: max_time_step = 5.0d-4
	Real*8  :: min_time_step = 1.0d-13
	Integer :: iteration
	Character*8 :: t_ofmt = '(ES12.5)'	! For formatted timestep output
Contains

	Subroutine Initialize_TimeStepping()
		Implicit None
		Real*8 :: dr, tdiff

		dr = radius(1)-radius(2)  ! assume uniform grid for now
		tdiff = dr*dr				  ! Viscous diffusion time across one grid point		
		! This may seem stringent, but it is in-line with our max time step in ASH
		!max_time_step = tdiff *10.0d0 
		!max_time_step = max_time_step*4.0d0
		!min_time_step = max_time_step*1.0d-4


		new_deltat   = max_time_step
		    deltat   = 0.0d0
		old_deltat   = 0.0d0
		new_timestep = .true.
		already_at_max = .true.	! Used to keep from changing timestep unecessarily
		! If checkpointing, change new_deltat and deltat appropriately
	End Subroutine Initialize_TimeStepping
	Subroutine Main_Loop()
		Implicit None
		Integer ::  last_iteration, first_iteration
			
		! We enter the main loop assuming that the solve has just been performed
		! and that the equation set structure's RHS contains our primary fields with 
		! radial dimension in-processor.
		! Care needs to be taken at init to ensure fields (W,Z,P,T) are stored
		! in the RHS (they are copied out upon entry into the loop).



		first_iteration = 1  ! hard-coded for now
		last_iteration = first_iteration + max_iterations-1
		Call Initialize_TimeStepping()
		If (chebyshev) Then
			! work structure for post_solve_cheby
			Call ctemp%init(field_count = wsfcount, config = 'p1b')
		Endif
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
			
		Enddo
	End Subroutine Main_Loop
	Subroutine Post_Solve()	
		Implicit None
		Integer :: m, i
		Character*12 :: tstring, otstring
		! wsp%p1b is assumed to be allocated
		
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
				wsp%p1a(:,i,m,tvar) = wsp%p1a(:,i,m,tvar)/radius(:)
			Enddo
		Enddo


		!///////////////////////////////
		!  Z Terms
		Call d_by_dx(zvar,d2zdr2,wsp%p1a,2)		! 2nd derivative will be overwritten with dzdr
		Call Add_Derivative(zeq,zvar,2,wsp%p1b,wsp%p1a,d2zdr2)	

		Call d_by_dx(zvar,dzdr,wsp%p1a,1)	

		Call Add_Derivative(zeq,zvar,0,wsp%p1b,wsp%p1a,zvar)	

		!Load the old ab array into the RHS
		Call Set_All_RHS(wsp%p1b)	! RHS now holds old_AB+CN factors

		Call wsp%deconstruct('p1b')

		Call wsp%reform()	! move from p1a to s2a
		
	End Subroutine Post_Solve

	Subroutine Post_Solve_Cheby()	
		Implicit None
		Integer :: m, i
		Character*12 :: tstring, otstring
		! Eventually I will merge the two post-solve routines
		! For now, keeping them separate

		! wsp%p1b is assumed to be allocated
		
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

		Call wsp%reform()	! move from p1a to s2a
		
	End Subroutine Post_Solve_Cheby


	!/////////////////////////////
	!  This routine calculates terms that involve
	!      theta derivatives and loads them into
	!		 the buffer.
	Subroutine rlm_spacea()
		Implicit None
		Integer :: mp
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

		! If magnetism compute B and J

		Call DeAllocate_rlm_Field(ftemp1)
		Call DeAllocate_rlm_Field(ftemp2)

		! Zero out l_max mode
		Do mp = my_mp%min, my_mp%max
			wsp%s2a(mp)%data(l_max,:) = 0.0d0
		Enddo


		!Legendre Transform and transpose the buffer
		Call wsp%construct('p2a')
		Call Legendre_Transform(wsp%s2a,wsp%p2a)
		Call wsp%deconstruct('s2a')
		wsp%config = 'p2a'	
		Call wsp%reform()	! We are now in p3a
		

	End Subroutine rlm_spacea




	Subroutine Compute_TimeStep()
		Implicit None
		Real*8 :: ovt2, ovht2, ovrt2, maxt2, dr
		Real*8 ::  maxt
		Integer :: r

		ovt2 = 0.0d0	! "over t squared"
		dr = radius(1)-radius(2)   ! Hard-coded uniform grid temporarily
		Do r = my_r%min, my_r%max
			ovht2 = Maxval(wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2)*OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
			ovt2 = Max(ovt2, ovht2)
			ovrt2 = Maxval(wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)	! radial
			ovt2 = Max(ovt2,ovrt2)
		Enddo

		Call Global_Max(ovt2, maxt2, grp = pfi%gcomm)	! Wrapper to mpi all reduce
		
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

		Call Phi_Derivatives()
		! Next perform the FFT
		Call fft_to_physical(wsp%p3a,rsc = .true.)

		! Convert all our terms of the form "sintheta var" to "var"
		Call sintheta_div(vtheta)	! sintheta vtheta to vtheta etc.
		Call sintheta_div(vphi)
		Call sintheta_div(dvtdr)
		Call sintheta_div(dvpdr)
		Call sintheta_div(dtdt)
		Call sintheta_div(dvrdt)


		Call sintheta_div(dvpdp)
		Call sintheta_div(dvtdp)
		

		!////////////////////////////////////////////////////////////////////////
		!This is a good spot to do some simple diagnostic output while we debug the code
		!since velocity components, Pressure, and Temperature are all 
		!in memory and in physical space at this point in time.
		Call ps_output(wsp%p3a, iteration)
		!////////////////////////////////////////////////////////////////////////

		Call compute_timestep()
		

		! We are now ready to built the nonlinear terms
		Call wsp%construct('p3b')
		wsp%config = 'p3b'

		Call Temperature_Advection()
		Call Momentum_Advection_Radial()
		Call Momentum_Advection_Theta()
		Call Momentum_Advection_Phi()

		Call wsp%deconstruct('p3a')
		Call fft_to_spectral(wsp%p3b, rsc = .true.)
		Call wsp%reform()	! Move to p2b

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
		Call Legendre_Transform(wsp%p2b,wsp%s2b)	

		Call wsp%deconstruct('p2b')
		wsp%config = 's2b'


		

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

		Call DeAllocate_rlm_Field(ftemp1)
		Call DeAllocate_rlm_Field(ftemp2)

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
		! Zero out l_max mode
		Do mp = my_mp%min, my_mp%max
			wsp%s2b(mp)%data(l_max,:) = 0.0d0
		Enddo

		Call wsp%reform() ! move to the solve space

	End Subroutine rlm_spaceb

	Subroutine Temperature_Advection()
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
	End Subroutine Temperature_Advection

	Subroutine Momentum_Advection_Radial()
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

		! Multiply by radius squared so that we have NL RHS U_r * r^2 (= NL RHS W * l(l+1)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,wvar) = wsp%p3b(:,r,t,wvar)*r_squared(r)
			Enddo
		Enddo	
	End Subroutine Momentum_Advection_Radial

	Subroutine Momentum_Advection_Theta()
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

		! At this point, we have [u dot grad u]_theta
		! Multiply by radius/sintheta so that we have r[u dot grad u]_theta/sintheta (getting ready for Z and dWdr RHS building)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,pvar) = wsp%p3b(:,r,t,pvar)*radius(r)/sintheta(t)
			Enddo
		Enddo	



	End Subroutine Momentum_Advection_Theta

	Subroutine Momentum_Advection_Phi()
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

		! At this point, we have [u dot grad u]_phi
		! Multiply by radius/sintheta so that we have r[u dot grad u]_phi/sintheta (getting ready for Z and dWdr RHS building)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				wsp%p3b(:,r,t,zvar) = wsp%p3b(:,r,t,zvar)*radius(r)/sintheta(t)
			Enddo
		Enddo	

	End Subroutine Momentum_Advection_Phi

	Subroutine Phi_Derivatives()
		Implicit None
		wsp%p3a(:,:,:,pvar) = wsp%p3a(:,:,:,tvar) ! cluge to keep t
		Call d_by_dphi(wsp%p3a,vr,dvrdp)
		Call d_by_dphi(wsp%p3a,vtheta,dvtdp)
		Call d_by_dphi(wsp%p3a,vphi,dvpdp)
		Call d_by_dphi(wsp%p3a,tvar,dtdp)
	End Subroutine Phi_Derivatives

	Subroutine sintheta_div(ind)
		Implicit None
		Integer, Intent(In) :: ind
		Integer :: t
		Do t = my_theta%min, my_theta%max
			wsp%p3a(:,:,t,ind) = wsp%p3a(:,:,t,ind)/sintheta(t)	! may want to do x csctheta (not sure)
		Enddo

	End Subroutine sintheta_div

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




	Subroutine AdvanceTime
		Implicit None
		! wsp will be in 'p1b' config
		! p1b contains the new adams bashforth term
		!Call print_max_spec2(pvar)

		if (.not. nonlinear) then
			wsp%p1b(:,:,:,:) = 0.0d0
		endif
		Call Add_to_All_RHS(wsp%p1b,new_ab_factor)
		Call Fix_Boundary_Conditions()
		Call Implicit_Solve()
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
		neq = 4
		nvar = 4

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
				amp = -Ra! *(radius/r_outer)
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
				amp = -(Ra/Ek)*radius/r_outer/H_Laplacian
				amp = -(Ra/Ek) ! *(radius/r_outer)**(-2)
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
					amp = H_Laplacian/Pm
					Call add_implicit_term(bteq,btvar, 0, amp,lp)					
					amp = 1.0d0/Pm
					Call add_implicit_term(bteq,btvar, 2, amp,lp)	

					!=========================================
					!  Bpol Equation
					amp = H_Laplacian/Pm
					Call add_implicit_term(bteq,btvar, 0, amp,lp)
					amp = 1.0d0/Pm
					Call add_implicit_term(bteq,btvar, 2, amp,lp)					
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
					Call Clear_Row(bpeq,lp,1)
					Call Clear_Row(bpeq,lp,N_R)
					Call Clear_Row(bteq,lp,1)
					Call Clear_Row(bteq,lp,N_R)


					! Match to a potential field at top and bottom
					! Btor = 0 at top and bottom
					r = 1
					Call Load_BC(lp,r,bteq,btvar,one,0)
					r = N_R
					Call Load_BC(lp,r,bteq,btvar,one,0)

					! dBpol/dr+ell*Bpol/r = 0 at outer boundary
					r = 1
					Call Load_BC(lp,r,bpeq,bpvar,one,1)
					samp = my_lm_lval(lp)*one_over_r(r)
					Call Load_BC(lp,r,bpeq,bpvar,samp,0)

					! dBpol/dr-ell(ell+1)*Bpol/r = 0 at inner boundary
					r = N_R
					Call Load_BC(lp,r,bpeq,bpvar,one,1)	
					samp = - l*(l+1)*One_Over_R(r)
					Call Load_BC(lp,r,bpeq,bpvar,samp,0)	

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
					equation_set(1,bpeq)%RHS(1,:,indx:indx+n_m) = Zero
					equation_set(1,bpeq)%RHS(N_R,:,indx:indx+n_m) = Zero
 
					equation_set(1,bteq)%RHS(1,:,indx:indx+n_m) = Zero
					equation_set(1,bteq)%RHS(N_R,:,indx:indx+n_m) = Zero

          
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
