Module Sphere_Spectral_Space
	Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
	Use Parallel_Framework
	Use Controls
	Use ProblemSize
	Use Finite_Difference, Only : d_by_dx
	Use Linear_Solve
	Use Fields
	Use BoundaryConditions
	Use Chebyshev_Polynomials, Only : cheby_to_spectral, cheby_from_spectral, d_by_dr_cp, cheby_from_spectralFE, &
            d_by_dr_cpFE, cheby_to_spectralFE
	Use ClockInfo
	Use Timers
	Use Sphere_Linear_Terms
	Implicit None
	Type(SphericalBuffer) :: ctemp ! workspace
Contains

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
			deltat = new_deltat
			new_timestep = .false.
			If (my_rank .eq. 0) Then
				Write(otstring,t_ofmt)old_deltat
				Write(tstring,t_ofmt)deltat
				Call stdout%print(' Timestep has changed from '//Trim(otstring)//' to '//Trim(tstring)//'.')
                Call stdout%partial_flush()  ! Make SURE that a changing timestep is recorded ...
                                             ! ... even at the expense of additional file I/O for redirected stdout
			Endif
			Call StopWatch(seteq_time)%startclock()
			Call Reset_Linear_Equations()
			Call StopWatch(seteq_time)%increment()
		Endif
		if (euler_step) then
				!Euler Step
				new_ab_factor = deltat
				old_ab_factor = 0.0d0
                euler_step = .false.
		else
				new_ab_factor = 0.5d0*deltat*(2 + deltat/old_deltat)
				old_ab_factor = -0.5d0*deltat**2/old_deltat
		endif
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
        Call Add_Derivative(peq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)
	    Call Add_Derivative(weq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
		Call Add_Derivative(weq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		If (deriv_cluge) Then
			Call d_by_dx(dwdr,d2wdr2,wsp%p1a,1)	! cluge like in ASH  --- seems unnecessary though.  take out once all else works
		Endif
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
        Call Add_Derivative(zeq,zvar,1,wsp%p1b,wsp%p1a,dzdr)

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
        If (output_iteration) Then
            !Convert p/rho to p
            ! We already took d/dr(p/rho), so we'll fix that later
		    Do m = 1, my_num_lm
			    Do i = 1, 2
				    wsp%p1a(:,i,m,pvar) = wsp%p1a(:,i,m,pvar)*ref%density(:)
			    Enddo
		    Enddo
            Call wsp%reform(nextra_recv = nicknum) ! The s2a buffer needs to be larger than normal
        Else
    		Call wsp%reform()	! move from p1a to s2a
        Endif
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
			deltat = new_deltat
			new_timestep = .false.
			If (my_rank .eq. 0) Then
				Write(otstring,t_ofmt)old_deltat
				Write(tstring,t_ofmt)deltat
				Call stdout%print(' Timestep has changed from '//Trim(otstring)//' to '//Trim(tstring)//'.')
                Call stdout%partial_flush()  ! Make SURE that a changing timestep is recorded ...
                                             ! ... even at the expense of additional file I/O for redirected stdout

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
		!!!DDDD Write(6,*)'I am getting the new rhs: ', my_rank
		Call Get_All_RHS(wsp%p1a)
		wsp%p1a((2*N_r)/3+1:N_r,:,:,:) = 0.0d0	! de-alias


		! This is terribly inefficient, but I just want to test the stability of Chebyshev vs. FD for not..
		! We'll create a new buffer.  ctemp
		! Store all the permanent derivatives there - in c space
		ctemp%nf1a = 4
		ctemp%nf1b = 4
		If (magnetism) then
			ctemp%nf1a = 5
			ctemp%nf1b = 5
		Endif
		Call ctemp%construct('p1a')
		! W..
		Call d_by_dr_cp(wvar,d3wdr3,wsp%p1a,3)
		ctemp%p1a(:,:,:,1) = wsp%p1a(:,:,:,d3wdr3)

		Call d_by_dr_cp(wvar,dwdr   ,wsp%p1a,1)		
		Call d_by_dr_cp(wvar,d2wdr2 ,wsp%p1a,2)
		! P....n
		Call d_by_dr_cp(pvar,dpdr1,wsp%p1a,1)
		ctemp%p1a(:,:,:,2) = wsp%p1a(:,:,:,dpdr1)
		! T
		Call d_by_dr_cp(tvar,d2tdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,3) = wsp%p1a(:,:,:,d2tdr2)
		Call d_by_dr_cp(tvar,dtdr,wsp%p1a,1)
		! Z..
		Call d_by_dr_cp(zvar,d2zdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,4) = wsp%p1a(:,:,:,d2zdr2)
		Call d_by_dr_cp(zvar,dzdr,wsp%p1a,1)

		! Magnetism
		If (magnetism) Then
			Call d_by_dr_cp(avar,d2adr2,wsp%p1a,2)	
			ctemp%p1a(:,:,:,5) = wsp%p1a(:,:,:,d2adr2)
			Call d_by_dr_cp(avar,dadr  ,wsp%p1a,1)
			Call d_by_dr_cp(cvar,dcdr  ,wsp%p1a,1)
			Call d_by_dr_cp(cvar,d2cdr2,wsp%p1a,2)
		Endif

		!//////////////////////////////////////////////////////////////////////////
		! Now everything we need is in the wsp or ctemp buffer
		! The ctemp terms are those terms that do not leave this configuration
		! transform them now & add them to appropriate equations
		Call ctemp%construct('p1b')
		ctemp%p1a((2*N_r)/3+1:N_r,:,:,:) = 0.0d0	! de-alias

		Call Cheby_From_Spectral(ctemp%p1a,ctemp%p1b)
        If (output_iteration) Then
            ! Grab dpdr
            Call cobuffer%construct('p1a')
            cobuffer%p1a(:,:,:,dpdr_cb) = ctemp%p1b(:,:,:,2)
        Endif


		Call Add_Derivative(peq,wvar,3,wsp%p1b,ctemp%p1b,1)
		Call Add_Derivative(weq,pvar,1,wsp%p1b,ctemp%p1b,2)
		Call Add_Derivative(teq,tvar,2,wsp%p1b,ctemp%p1b,3)
		Call Add_Derivative(zeq,zvar,2,wsp%p1b,ctemp%p1b,4)
		If (magnetism) Then
			Call Add_Derivative(aeq,avar,2,wsp%p1b,ctemp%p1b,5)
		Endif
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
		Call Add_Derivative(peq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(weq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
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
		Call Add_Derivative(zeq,zvar,1,wsp%p1b,wsp%p1a,dzdr)

		!///////////////////////////////////////
		!  Magnetic Terms
		If (magnetism) Then
			!//////////////
			! A-terms (Toroidal magnetic field)
		
			Call Add_Derivative(aeq,avar,0,wsp%p1b,wsp%p1a,avar)

			!///////////////////
			! C-terms (Poloidal magnetic field)
	
			Call Add_Derivative(ceq,cvar,2,wsp%p1b,wsp%p1a,d2cdr2)

			Call Add_Derivative(ceq,cvar,0,wsp%p1b,wsp%p1a,cvar)

		Endif


		!Load the old ab array into the RHS
		Call Set_All_RHS(wsp%p1b)	! RHS now holds old_AB+CN factors


		Call wsp%deconstruct('p1b')
		Call StopWatch(psolve_time)%increment()

		Call StopWatch(ctranspose_time)%startclock()


    	

        If (output_iteration) Then
            !Convert p/rho to p
            ! We already took d/dr(p/rho), so we'll fix that later
		    Do m = 1, my_num_lm
			    Do i = 1, 2
				    wsp%p1a(:,i,m,pvar) = wsp%p1a(:,i,m,pvar)*ref%density(:)
			    Enddo
		    Enddo
            Call cobuffer%reform()
        Endif
        Call wsp%reform()	! move from p1a to s2a
		Call StopWatch(ctranspose_time)%increment()

	End Subroutine Post_Solve_Cheby


	Subroutine Post_Solve_FE()	
		Implicit None
		Integer :: m, i,j, jstart, jend
		Character*12 :: tstring, otstring
		! Eventually I will merge the two post-solve routines
		! For now, keeping them separate

		! wsp%p1b is assumed to be allocated
		Call StopWatch(psolve_time)%startclock()
		Call wsp%construct('p1a')
		wsp%config = 'p1a'

		old_deltat = deltat
		If (new_timestep) Then
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
		!!!DDDD Write(6,*)'I am getting the new rhs: ', my_rank
		Call Get_All_RHS(wsp%p1a)
        
        !Write(6,*)'MAXVAL:  ', maxval(wsp%p1a)
        !Write(6,*)'MINVAL:  ', minval(wsp%p1a)
		! de-alias  each subdomain

        jstart = (2*fencheby)/3  ! might play with de-aliasing type...
        jend = fencheby
        Do j = 1, fensub
            wsp%p1a(jstart:jend,:,:,:) = 0.0d0
            jstart = jstart+fencheby
            jend = jend+fencheby
        Enddo
		! We'll create a new buffer.  ctemp
		! Store all the permanent derivatives there - in c space
		ctemp%nf1a = 4
		ctemp%nf1b = 4
		If (magnetism) then
			ctemp%nf1a = 5
			ctemp%nf1b = 5
		Endif
		Call ctemp%construct('p1a')
		! W..
		Call d_by_dr_cpFE(wvar,d3wdr3,wsp%p1a,3)
		ctemp%p1a(:,:,:,1) = wsp%p1a(:,:,:,d3wdr3)

		Call d_by_dr_cpFE(wvar,dwdr   ,wsp%p1a,1)		
		Call d_by_dr_cpFE(wvar,d2wdr2 ,wsp%p1a,2)
		! P....n
		Call d_by_dr_cpFE(pvar,dpdr,wsp%p1a,1)
		ctemp%p1a(:,:,:,2) = wsp%p1a(:,:,:,dpdr)
		! T
		Call d_by_dr_cpFE(tvar,d2tdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,3) = wsp%p1a(:,:,:,d2tdr2)
		Call d_by_dr_cpFE(tvar,dtdr,wsp%p1a,1)
		! Z..
		Call d_by_dr_cpFE(zvar,d2zdr2,wsp%p1a,2)
		ctemp%p1a(:,:,:,4) = wsp%p1a(:,:,:,d2zdr2)
		Call d_by_dr_cpFE(zvar,dzdr,wsp%p1a,1)

		! Magnetism
		If (magnetism) Then
			Call d_by_dr_cpFE(avar,d2adr2,wsp%p1a,2)	
			ctemp%p1a(:,:,:,5) = wsp%p1a(:,:,:,d2adr2)
			Call d_by_dr_cpFE(avar,dadr  ,wsp%p1a,1)
			Call d_by_dr_cpFE(cvar,dcdr  ,wsp%p1a,1)
			Call d_by_dr_cpFE(cvar,d2cdr2,wsp%p1a,2)
            !wsp%p1a(:,:,:,d2cdr2) = 0.0d0 ! DEBUG
		Endif

		!//////////////////////////////////////////////////////////////////////////
		! Now everything we need is in the wsp or ctemp buffer
		! The ctemp terms are those terms that do not leave this configuration
		! transform them now & add them to appropriate equations
		Call ctemp%construct('p1b')
		
		! de-alias  each subdomain
        jstart = (2*fencheby)/3  ! might play with de-aliasing type...
        jend = fencheby
        Do j = 1, fensub
            ctemp%p1a(jstart:jend,:,:,:) = 0.0d0
        Enddo

        !Do m = my_lm_min, my_lm_max
        !    If (l_lm_values(m) .eq. 0) Then
        !        Write(6,*)'tvar 0 (A):   ', wsp%p1a(:,:,m,tvar)
        !    Endif
        !Enddo
		Call Cheby_From_SpectralFE(ctemp%p1a,ctemp%p1b)


		Call Add_Derivative(peq,wvar,3,wsp%p1b,ctemp%p1b,1)
		Call Add_Derivative(weq,pvar,1,wsp%p1b,ctemp%p1b,2)
		Call Add_Derivative(teq,tvar,2,wsp%p1b,ctemp%p1b,3)
		Call Add_Derivative(zeq,zvar,2,wsp%p1b,ctemp%p1b,4)
		If (magnetism) Then
			Call Add_Derivative(aeq,avar,2,wsp%p1b,ctemp%p1b,5)
		Endif
		Call ctemp%deconstruct('p1a')
		Call ctemp%deconstruct('p1b')

		!//////////////////////////////////////////////
		!  Next, we reconstruct ctemp%p1a and copy wsp%p1a into it
		ctemp%nf1a = wsp%nf1a
		Call ctemp%construct('p1a')
		ctemp%p1a(:,:,:,:) = wsp%p1a(:,:,:,:)

		! de-alias  each subdomain
        jstart = (2*fencheby)/3  ! might play with de-aliasing type...
        jend = fencheby
        Do j = 1, fensub
            ctemp%p1a(jstart:jend,:,:,:) = 0.0d0
        Enddo

		wsp%p1a(:,:,:,:) = 0.0d0	! Shouldn't need to do this, but just to be sure
		Call Cheby_From_SpectralFE(ctemp%p1a,wsp%p1a)


        !Write(6,*)'MAXVAL2:  ', maxval(wsp%p1a)
        !Write(6,*)'MINVAL2:  ', minval(wsp%p1a)
        !Do m = my_lm_min, my_lm_max
        !    If (l_lm_values(m) .eq. 0) Then
        !        Write(6,*)'tvar 0 (Post):   ', wsp%p1a(:,:,m,tvar)
        !    Endif
        !Enddo


		Call ctemp%deconstruct('p1a')

		!/////////////////////////////////////////////////////////////////
		!  The rest of the code can remain unchanged
		!/////////////////////////////////////////////////////////////////	
		!Load the W derivatives into the appropriate RHS's  

		Call Add_Derivative(peq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(peq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
		Call Add_Derivative(peq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

		Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)	
		Call Add_Derivative(weq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
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
		Call Add_Derivative(zeq,zvar,1,wsp%p1b,wsp%p1a,dzdr)

		!///////////////////////////////////////
		!  Magnetic Terms
		If (magnetism) Then
			!//////////////
			! A-terms (Toroidal magnetic field)
		
	

			Call Add_Derivative(aeq,avar,0,wsp%p1b,wsp%p1a,avar)

			!///////////////////
			! C-terms (Poloidal magnetic field)
	
			Call Add_Derivative(ceq,cvar,2,wsp%p1b,wsp%p1a,d2cdr2)

			Call Add_Derivative(ceq,cvar,0,wsp%p1b,wsp%p1a,cvar)

		Endif


		!Load the old ab array into the RHS
        !wsp%p1a = 0.0d0
        !wsp%p1b = 0.0d0
		Call Set_All_RHS(wsp%p1b)	! RHS now holds old_AB+CN factors


		Call wsp%deconstruct('p1b')
		Call StopWatch(psolve_time)%increment()

		Call StopWatch(ctranspose_time)%startclock()

        If (output_iteration) Then
            Call wsp%reform(nextra_recv = nicknum) ! The s2a buffer needs to be larger than normal
        Else
    		Call wsp%reform()	! move from p1a to s2a
        Endif

		Call StopWatch(ctranspose_time)%increment()
	End Subroutine Post_Solve_FE


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
		Call Enforce_Boundary_Conditions()
		Call StopWatch(solve_time)%startclock()
		Call Implicit_Solve()
		Call StopWatch(solve_time)%increment()
		simulation_time = simulation_time+deltat
		! The righthand side of the equation set structure
		! Now contains the updated fields.
	End Subroutine AdvanceTime
	Subroutine Finalize_EMF()
		Implicit None
		Integer m, i,j,jstart,jend
		! we need to take one last radial derivative and combine terms

		If (chebyshev .or. finite_element) Then
			! Again, terribly inefficient, but we are looking to check the MHD right now.
			! Will optimize this later.
			ctemp%nf1a = 2
			ctemp%nf1b = 2
			Call ctemp%construct('p1a')
			Call ctemp%construct('p1b')
			ctemp%p1a(:,:,:,:) = 0.0d0
			Do m = 1, my_num_lm
				Do i = 1, 2
					ctemp%p1a(:,i,m,1) = wsp%p1b(:,i,m,emfphi)
				Enddo
			Enddo

            If (finite_element) Then
    			Call Cheby_To_SpectralFE(ctemp%p1a,ctemp%p1b)
    			Call d_by_dr_cpFE(1,2,ctemp%p1b,1)


		        ! de-alias  each subdomain
                jstart = (2*fencheby)/3  ! might play with de-aliasing type...
                jend = fencheby
                Do j = 1, fensub
                    ctemp%p1b(jstart:jend,:,:,:) = 0.0d0
                    jstart = jstart+fencheby
                    jend = jend+fencheby
                Enddo
    			Call Cheby_From_SpectralFE(ctemp%p1b,ctemp%p1a)
                !Write(6,*)'Checking: ', maxval(ctemp%p1a), minval(ctemp%p1a)
            Else
    			Call Cheby_To_Spectral(ctemp%p1a,ctemp%p1b)
    			Call d_by_dr_cp(1,2,ctemp%p1b,1)
    			ctemp%p1b((2*N_r)/3+1:N_r,:,:,:) = 0.0d0
    			Call Cheby_From_Spectral(ctemp%p1b,ctemp%p1a)
            Endif


			Do m = 1, my_num_lm
				Do i = 1, 2
					wsp%p1b(:,i,m,avar) = wsp%p1b(:,i,m,avar) + ctemp%p1a(:,i,m,2)
				Enddo
			Enddo
			Call ctemp%deconstruct('p1a')
			Call ctemp%deconstruct('p1b')
		Else
			ctemp%nf1a = 1
			Call ctemp%construct('p1a')
			ctemp%p1a(:,:,:,:) = 0.0d0
			Do m = 1, my_num_lm
				Do i = 1, 2
					ctemp%p1a(:,i,m,1) = wsp%p1b(:,i,m,avar)
				Enddo
			Enddo
			Call d_by_dx(emfphi,avar,wsp%p1b,1)
		
			Do m = 1, my_num_lm
				Do i = 1, 2
					wsp%p1b(:,i,m,avar) = ctemp%p1a(:,i,m,1) + wsp%p1b(:,i,m,avar)
				Enddo
			Enddo
			Call ctemp%deconstruct('p1a')
		Endif

	End Subroutine Finalize_EMF

End Module Sphere_Spectral_Space
