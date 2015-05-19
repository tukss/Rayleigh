Module Linear_Terms_Sphere
	Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
	Use Controls
	Use ProblemSize
	Use Linear_Solve
	Use Fields
	Use BoundaryConditions
	Use Timers
	Use ClockInfo
	!Use Legendre_Polynomials, Only : pi
	Use ReferenceState
	Use TransportCoefficients
	Use NonDimensionalization
    Use Equation_Coefficients
    !Use Math_Constants, Only : pi
    Implicit None
    Real*8, Allocatable :: Lconservation_weights(:)
	
Contains	
	Subroutine Linear_Init()
		Implicit None
        Real*8 :: amp, T,arg
        Integer :: n, r
		Call Initialize_Benchmark_Equations()
		If (strict_L_conservation) Then
            Allocate(Lconservation_weights(1:N_R))
            Lconservation_weights(1:N_R) = 0.0d0
            If (chebyshev) Then
                            
                amp = Pi / (N_R*1.0d0) 
                do n = 1, N_R
                    do r = 1, N_R
                        arg = (n-1.d0) * (r-1.d0+0.5d0) * amp
                        T = Cos(arg) 
                        Lconservation_weights(n) = Lconservation_weights(n) + radial_integral_weights(r) * T
                    enddo
                enddo
                Lconservation_weights(1) = Lconservation_weights(1)*0.5d0
                Lconservation_weights(N_R) = Lconservation_weights(N_r)*0.5d0
                Lconservation_weights( (2*N_R)/3+1: ) = 0.0d0  ! De-Alias here for now
            Else
                Lconservation_weights(1:N_R) = radial_integral_weights(1:N_R)
            Endif
        Endif
	End Subroutine Linear_Init

	Subroutine Reset_Linear_Equations()
		Implicit None
		Real*8 :: rhs_factor, lhs_factor
		Call Reset_Equation_Coefficients()	! Zero out all coefficients
		lhs_factor = -deltat*alpha_implicit	! Crank Nicolson scheme - alpha = 0.5
		rhs_factor = deltat*(1.0d0-alpha_implicit)
		Call Set_Time_Factors(lhs_factor,rhs_factor)
		Call Compute_Benchmark_Coefficients()
		!Call Set_Boundary_Conditions()
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
        If (finite_element) Call Use_Finite_Elements()
		Call Initialize_Equation_Set(neq,nvar,N_R,my_nl_lm, my_nm_lm,2)

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
		If (bandsolve) Call Use_BandSolve()
		!============================================
		!  Next decide which derivatives need to be saved.
		!  This can be either for transposition later (i.e. nonlinear terms),
		!  or for checkpointing and output (such at with P).
		!  Saved derivatives will remain in memory until explicitly deallocated
		!  by the user.  The Implicit Module simply places them into memory.  It 
		!  does not deallocate them.

		!Call Set_Deriv_Save(wvar,2)
		!Call Set_Deriv_Save(pvar,0)
		!Call Set_Deriv_Save(tvar,1)
		!Call Set_Deriv_Save(zvar,1)
		!If (magnetism) Then
		!	Call Set_Deriv_Save(avar,1)
		!	Call Set_Deriv_Save(cvar,2)
		!Endif
	End Subroutine Initialize_Benchmark_Equations

	Subroutine Compute_Benchmark_Coefficients()
		Implicit None

		Real*8, Allocatable :: H_Laplacian(:), amp(:)
		Integer :: l, lp
		!rmin_norm

		Allocate(amp(1:N_R))
		Allocate(H_Laplacian(1:N_R))
		Do lp = 1, my_nl_lm
            If (bandsolve) Call DeAllocate_LHS(lp)
            Call Allocate_LHS(lp)
			l = my_lm_lval(lp)		

			H_Laplacian = - l_l_plus1(l) * OneOverRSquared
			If (l .eq. 0) Then
				!====================================================
				!			Temperature Equation
				! T only
				amp = 1.0d0
				Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)	! Time independent part

				!amp = 2.0d0/radius/Pr
                amp = 2.0d0/radius*kappa
				Call add_implicit_term(teq,tvar, 1, amp,lp)
				!amp = 1.0d0/Pr
                amp = 1.0d0*kappa
				Call add_implicit_term(teq,tvar, 2, amp,lp)

				! Kappa,rho, T variation in radius
				amp = S_Diffusion_Coefs_1
				Call add_implicit_term(teq,tvar,1,amp,lp)

				!=======================================
				!   Hydrostatic balance
				! 	 This part of the equation is static (i.e. not time-evolving)

				! 	 t

                amp = -ref%gravity_term_s
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

                amp = -ref%gravity_term_s/H_Laplacian
				Call add_implicit_term(weq, tvar, 0, amp,lp)			! Gravity

				! Pressure
				!amp = 1.0d0/(Ek*H_Laplacian)*ref%density		! dPdr
                amp = dpdr_W_term/H_Laplacian
				Call add_implicit_term(weq,pvar, 1, amp,lp)


				! W
				amp = 1.0d0
				Call add_implicit_term(weq,wvar, 0, amp,lp,static = .true.)	! This term does not a get a dt factor

				!amp = H_Laplacian		! Diffusion
                amp = H_Laplacian*nu
				Call add_implicit_term(weq,wvar, 0, amp,lp)
				!amp = 1.0d0
                amp = nu
				Call add_implicit_term(weq,wvar, 2, amp,lp)

				! These two diffusion bits are different 
				! depending on variation of rho and nu
				amp = W_Diffusion_Coefs_0		
				Call add_implicit_term(weq,wvar, 0, amp,lp)
				amp = W_Diffusion_Coefs_1
				Call add_implicit_term(weq,wvar, 1, amp,lp)

				!==================================================
				!				Pressure (dWdr) Equation
				
				! Pressure
				!amp = -(1.0d0)/Ek*ref%density	
                amp = pressure_dwdr_term
				Call add_implicit_term(peq,pvar, 0, amp,lp)

				! W
				amp = 1.0d0
				Call add_implicit_term(peq,wvar, 1, amp,lp, static = .true.)	! Time independent term
				!amp =-H_Laplacian*2.0d0/radius	
				amp =-nu*H_Laplacian*2.0d0/radius
				Call add_implicit_term(peq,wvar, 0, amp,lp)
                !amp = H_Laplacian
				amp = H_Laplacian*nu
				Call add_implicit_term(peq,wvar, 1, amp,lp)
				!amp = 1.0d0
                amp = nu
				Call add_implicit_term(peq,wvar, 3, amp,lp)


				! Again, these two bits depend on variation of rho and nu
				amp = dW_Diffusion_Coefs_0*H_Laplacian
				Call add_implicit_term(peq,wvar, 0, amp,lp)
				amp = dW_Diffusion_Coefs_1
				Call add_implicit_term(peq,wvar, 1, amp,lp)				
				amp = dW_Diffusion_Coefs_2	
				Call add_implicit_term(peq,wvar, 2, amp,lp)
				!====================================================
				!			Temperature Equation

				! T 
				amp = 1.0d0
				Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)		! Time independent term

				!amp = H_Laplacian/Pr		! Diffusion
                amp = H_Laplacian*kappa
				Call add_implicit_term(teq,tvar, 0, amp,lp)
				!amp = 2.0d0/radius/Pr
				amp = 2.0d0/radius*kappa
				Call add_implicit_term(teq,tvar, 1, amp,lp)

                ! amp = 1.0d0/Pr
				amp = kappa
				Call add_implicit_term(teq,tvar, 2, amp,lp)

				! Kappa,rho, T variation in radius
				amp = S_Diffusion_Coefs_1          !/Pr
				Call add_implicit_term(teq,tvar,1,amp,lp)

                !Reference State Advection (only do this if reference state is non-adiabatic)
                If (advect_reference_state) Then
                    amp = -H_Laplacian/ref%density*ref%dsdr
                    Call add_implicit_term(teq,wvar,0,amp,lp)
                Endif
				
				!=====================================================
				!	Z Equation
				amp = 1.0d0
				Call add_implicit_term(zeq,zvar, 0, amp,lp, static = .true.)	! Time-independent piece

				!amp = H_Laplacian
                amp = H_Laplacian*nu
				Call add_implicit_term(zeq,zvar, 0, amp,lp)				
				!amp = 1.0d0
                amp = nu
				Call add_implicit_term(zeq,zvar, 2, amp,lp)				

				! Variation of rho and nu
				amp = Z_Diffusion_Coefs_0
				Call add_implicit_term(zeq,zvar, 0, amp,lp)
				amp = Z_Diffusion_Coefs_1
				Call add_implicit_term(zeq,zvar, 1, amp,lp)
				If (magnetism) Then
					!=========================================
					!  Btor Equation
					amp = 1.0d0
					Call add_implicit_term(aeq,avar, 0, amp,lp, static = .true.)	! Time-independent piece
					amp = H_Laplacian*eta
					Call add_implicit_term(aeq,avar, 0, amp,lp)					
					amp = 1.0d0*eta
					Call add_implicit_term(aeq,avar, 2, amp,lp)	

					! Eta variation in radius
					amp = A_Diffusion_Coefs_1
					Call add_implicit_term(aeq,avar,1,amp,lp)

					!=========================================
					!  Bpol Equation
					amp = 1.0d0
					Call add_implicit_term(ceq,cvar, 0, amp,lp, static = .true.)	! Time-independent piece
					amp = H_Laplacian*eta
					Call add_implicit_term(ceq,cvar, 0, amp,lp)
					amp = 1.0d0*eta
					Call add_implicit_term(ceq,cvar, 2, amp,lp)					
				Endif
			
				! If band solve, do the redefinition of the matrix here

			Endif
            Call Set_Boundary_Conditions(lp)
            If (bandsolve) Then
                Call Band_Arrange(weq,lp)
				Call Band_Arrange(zeq,lp)
                If (magnetism) Then
				    Call Band_Arrange(aeq,lp)
					Call Band_Arrange(ceq,lp)
                Endif
            Endif
		Enddo
		DeAllocate(amp)
		DeAllocate(H_Laplacian)
	End Subroutine Compute_Benchmark_Coefficients

	Subroutine Set_Boundary_Conditions(mode_ind)
        ! Modified version of set_boundary_conditions
        ! Designed to work with more memory friendly logic
        ! Sets boundary condition of indicated l-value (index lp)
        ! only.  Does not loop over lp.
		Implicit None
		Real*8 :: samp,one
        Integer, Intent(In) :: mode_ind
		Integer :: l, r,lp
		one = 1.0d0
        lp = mode_ind
		!Do lp = 1, my_nl_lm
			l = my_lm_lval(lp)

			If (l .eq. 0) Then
				Call Clear_Row(peq,lp,1)			! Pressure only has one boundary condition
				Call Clear_Row(teq,lp,1)
				Call Clear_Row(teq,lp,N_R)
                If (finite_element) Then
                    Call FEContinuity(peq,lp,pvar,fencheby,0) ! Pressure is continuous
                    Call FEContinuity(peq,lp,pvar,1,1)   ! dPdr is continuous (for ell = 0)

                    Call FEContinuity(teq,lp,tvar,fencheby,0) ! T/S is continuous
                    Call FEContinuity(teq,lp,tvar,1,1)   ! T' / S' is continuous (for ell = 0)

                Endif

				!*******************************************************
				! Entropy Boundary Conditions

				! Temperature Boundary Conditions (T fixed bottom and top)
				r = 1
                If (fix_tvar_top) Then
    				Call Load_BC(lp,r,teq,tvar,one,0)	!upper boundary
                Endif
                If (fix_dtdr_top) Then
                    Call Load_BC(lp,r,teq,tvar,one,1)	
                Endif

				r = N_R
                If (fix_tvar_bottom) Then
    				Call Load_BC(lp,r,teq,tvar,one,0)	! lower boundary
                Endif
                If (fix_dtdr_bottom) Then
                    Call Load_BC(lp,r,teq,tvar,one,1)
                Endif	    


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
                If (finite_element) Then
                    Call FEContinuity(peq,lp,pvar,fencheby,0)   ! Pressure is continuous
                    Call FEContinuity(peq,lp,wvar,1,2)          ! W'' is continuous 

                    Call FEContinuity(weq,lp,wvar,fencheby,0)   ! W is continuous
                    Call FEContinuity(weq,lp,wvar,1,1)          ! W' is continuous 

                    Call FEContinuity(teq,lp,tvar,fencheby,0)   ! T/S is continuous
                    Call FEContinuity(teq,lp,tvar,1,1)          ! T' / S' is continuous (for ell = 0)


                    Call FEContinuity(zeq,lp,zvar,fencheby,0)   ! Z is continuous
                    Call FEContinuity(zeq,lp,zvar,1,1)          ! Z' is continuous 
                Endif

				!*******************************************************
				! Entropy Boundary Conditions

				! Temperature Boundary Conditions (T fixed bottom and top)
				r = 1
                If (fix_tvar_top) Then
                    If (.not. fix_divrfc_top) Then
        				Call Load_BC(lp,r,teq,tvar,one,0)	!upper boundary
                    Endif
                Endif
                If (fix_dtdr_top) Then
                    Call Load_BC(lp,r,teq,tvar,one,1)	
                Endif

				r = N_R
                If (fix_tvar_bottom) Then
    				Call Load_BC(lp,r,teq,tvar,one,0)	! lower boundary
                Endif
                If (fix_dtdr_bottom) Then
                    Call Load_BC(lp,r,teq,tvar,one,1)
                Endif	    			

                !////////////////////////////////////////
                ! These are four different boundary conditions that are similar, though
                ! slightly different, in nature.  The idea is to try some boundary conditions
                ! that allow entropy and it's derivatives to vary on the boundary
                ! Either Del dot Grad S, Del dot F_conductive, or Del_r dot Grad S or F_conductive is zero
                
                If (fix_divrt_top) Then
                    r = 1
                    !d2sdr2 + 2/r dsdr = 0 at r = r_top
                    Call Load_BC(lp,r,teq,tvar,one,2)
                    samp = 2.0d0/radius(r)
                    Call Load_BC(lp,r,teq,tvar,samp,1)
                Endif

                If (fix_divt_top) Then
                    r = 1
                    !d2sdr2 + 2/r dsdr -l(l+1)/r^2= 0 at r = r_top
                    Call Load_BC(lp,r,teq,tvar,one,2)
                    samp = 2.0d0/radius(r)
                    Call Load_BC(lp,r,teq,tvar,samp,1)
                    samp = -l*1.0d0*(l*1.0d0+1)/radius(r)**2
                    Call Load_BC(lp,r,teq,tvar,samp,0)
                Endif


                If (fix_divrFc_top) Then
                    r = 1
                    !d2sdr2 + 2/r dsdr = 0 at r = r_top
                    Call Load_BC(lp,r,teq,tvar,one,2)
                    samp = 2.0d0/radius(r)+(dlnkappa(r)+ref%dlnrho(r)+ref%dlnT(r))
                    Call Load_BC(lp,r,teq,tvar,samp,1)
                Endif
                If (fix_divFc_top) Then
                    r = 1
                    if ( l .gt. 1) Then
                        !d2sdr2 + 2/r dsdr -l(l+1)/r^2= 0 at r = r_top
                        Call Load_BC(lp,r,teq,tvar,one,2)
                        samp = 2.0d0/radius(r)+(dlnkappa(r)+ref%dlnrho(r)+ref%dlnT(r))
                        Call Load_BC(lp,r,teq,tvar,samp,1)
                        samp = -l*1.0d0*(l*1.0d0+1)/radius(r)**2
                        Call Load_BC(lp,r,teq,tvar,samp,0)
                    Else
                        ! This boundary condition generally works well, but it can develop an ell=1
                        ! convection mode that causes undesired (and probably unphysical) behavior
                        samp = 1.0d0
                        Call Load_BC(lp,r,teq,tvar,samp,0)
                    Endif
                Endif


				!************************************************************
				! Velocity Boundary Conditions
		
				! Impenetrable top and bottom
				! W vanishes at the boundaries
				r = 1
				Call Load_BC(lp,r,weq,wvar,one,0)
				r = N_R
				Call Load_BC(lp,r,weq,wvar,one,0)

		
                If (no_slip_boundaries) Then		
				! No Slip Top and Bottom
				! Z and dWdr vanish at the boundaries
                    r = 1

                    Call Load_BC(lp,r,zeq,zvar,one,0)
				    Call Load_BC(lp,r,peq,wvar,one,1)
				    r = N_R
				    Call Load_BC(lp,r,peq,wvar,one,1)
				    Call Load_BC(lp,r,zeq,zvar,one,0)
                Else
                    ! stress-free boundaries
                    r = 1
                    samp = -(2.0d0/radius(r)+ref%dlnrho(r))
                    Call Load_BC(lp,r,peq,wvar,one,2)
                    Call Load_BC(lp,r,peq,wvar,samp,1)


                    Call Load_BC(lp,r,zeq,zvar,one,1)
                    Call Load_BC(lp,r,zeq,zvar,samp,0)


                    r = N_R
                    samp = -(2.0d0/radius(r)+ref%dlnrho(r))
                    Call Load_BC(lp,r,peq,wvar,one,2)
                    Call Load_BC(lp,r,peq,wvar,samp,1)
                    Call Load_BC(lp,r,zeq,zvar,one,1)
                    Call Load_BC(lp,r,zeq,zvar,samp,0)
                Endif
                If ((l .eq. 1) .and. (strict_L_Conservation) ) then
                   !    write(6,*)'Conserving Angular Momentum'
                    Call Clear_Row(zeq,lp,1)
    			    Call Load_BC(lp,1,zeq,zvar,one,0,integral = Lconservation_weights)
                Endif


				!*******************************************************
				!		Magnetic Boundary Conditions

				If (Magnetism) Then
					!  Clear the boundary rows
					Call Clear_Row(ceq,lp,1)
					Call Clear_Row(ceq,lp,N_R)
					Call Clear_Row(aeq,lp,1)
					Call Clear_Row(aeq,lp,N_R)

                    If (finite_element) Then
                        Call FEContinuity(aeq,lp,avar,fencheby,0)   ! A is continuous
                        Call FEContinuity(aeq,lp,avar,1,1)          ! A' is continuous

                        Call FEContinuity(ceq,lp,cvar,fencheby,0)   ! C is continuous
                        Call FEContinuity(ceq,lp,cvar,1,1)          ! C' is continuous
 
                    Endif

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

					! dBpol/dr-(ell+1)*Bpol/r = 0 at inner boundary
					r = N_R
					Call Load_BC(lp,r,ceq,cvar,one,1)	
					samp = - (l+1)*One_Over_R(r)
					Call Load_BC(lp,r,ceq,cvar,samp,0)	


                    If (fix_poloidalfield_top) Then
					    Call Clear_Row(ceq,lp,1)
					    Call Clear_Row(aeq,lp,1)
                        r = 1
					    Call Load_BC(lp,r,ceq,cvar,one,0)
                        Call Load_BC(lp,r,aeq,avar,one,0)
                    Endif
                    If (fix_poloidalfield_bottom) Then
					    Call Clear_Row(aeq,lp,N_R)
					    Call Clear_Row(ceq,lp,N_R)
                        r = N_R
					    Call Load_BC(lp,r,ceq,cvar,one,0)
					    Call Load_BC(lp,r,aeq,avar,one,0)
                    Endif

				Endif	! Magnetism


			Endif ! l = 0 or not
		!Enddo


	End Subroutine Set_Boundary_Conditions



	Subroutine Fix_Boundary_Conditions()
		Implicit None
		Integer :: l, indx, ii,lp, j, k
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

                    If (fix_poloidalfield_top) Then
                        If (l .eq. 1) Then
                        Do k = indx, indx+n_m
                            If (m_lm_values(k) .eq. 0) Then
                                equation_set(1,ceq)%RHS(1,1,k) = C10_top
                                equation_set(1,ceq)%RHS(1,2,k) = 0.0d0
                            Endif
                            If (m_lm_values(k) .eq. 1) Then
                                equation_set(1,ceq)%RHS(1,1,k) = C11_top
                                equation_set(1,ceq)%RHS(1,2,k) = C1m1_top
                            Endif
                        Enddo
                        Endif
                    Endif

                    If (fix_poloidalfield_bottom) Then
                        If (l .eq. 1 ) Then
                        Do k = indx, indx+n_m
                            If (m_lm_values(k) .eq. 0) Then
                                equation_set(1,ceq)%RHS(N_R,1,k) = C10_bottom
                                equation_set(1,ceq)%RHS(N_R,2,k) = 0.0d0
                            Endif
                            If (m_lm_values(k) .eq. 1) Then
                                equation_set(1,ceq)%RHS(N_R,1,k) = C11_bottom
                                equation_set(1,ceq)%RHS(N_R,2,k) = C1m1_bottom
                            Endif
                        Enddo
                        Endif
                    Endif

          
				Endif

			Else

				equation_set(1,weq)%RHS(:,2,indx)   = Zero ! no imaginary part for any ell=0 equations
				equation_set(1,zeq)%RHS(:,:,indx)   = Zero	! no ell = 0 z_equation
				if (magnetism) then
					equation_set(1,aeq)%RHS(:,:,indx)   = Zero
					equation_set(1,ceq)%RHS(:,:,indx)   = Zero
				endif
				equation_set(1,weq)%rhs(1:N_R,:,indx) = zero				! ell =0 W is zero
				equation_set(1,weq)%rhs(N_R+1,1,indx) = zero	! Pressure node


                If (fix_tvar_top) Then
				!Top temperature (in spectral space, but BC's specified in physical space
				!    so multiply by sqrt4pi)
				equation_set(1,weq)%RHS(1+ii,1,indx)   = T_Top*sqrt(4.0D0*Pi)
                Endif
                If (fix_dtdr_top) Then
				!Top temperature (in spectral space, but BC's specified in physical space
				!    so multiply by sqrt4pi)
				equation_set(1,weq)%RHS(1+ii,1,indx)   = dTdr_Top*sqrt(4.0D0*Pi)
                Endif

                If (fix_tvar_bottom) Then
				!Bottom Temperature
				equation_set(1,weq)%RHS(N_R+ii,1,indx) = T_Bottom*sqrt(4.0D0*Pi)
                Endif
            
                If (fix_dtdr_bottom) Then
				!Bottom Temperature
				equation_set(1,weq)%RHS(N_R+ii,1,indx) = dTdr_Bottom*sqrt(4.0D0*Pi)
                Endif

         Endif

         If (finite_element) Then
            ! Apply continuity conditions across the subdomains
            
            Do j = fencheby, N_R-1, fencheby
                equation_set(1,weq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
                equation_set(1,weq)%RHS(N_R  +j, : , indx:indx+n_m) = 0.0d0
                equation_set(1,weq)%RHS(2*N_R+j, : , indx:indx+n_m) = 0.0d0
                if (l .ne. 0) equation_set(1,zeq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
            Enddo
            If (Magnetism) Then
                if (l .ne. 0) then
                Do j = fencheby, N_R-1, fencheby
                    equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                    equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                Enddo
                endif
            Endif

            Do j = fencheby+1, N_R-1, fencheby
                equation_set(1,weq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                equation_set(1,weq)%RHS(N_R+j  ,: , indx:indx+n_m) = 0.0d0
                equation_set(1,weq)%RHS(2*N_R+j,: , indx:indx+n_m) = 0.0d0
                if (l .ne. 0) equation_set(1,zeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
            Enddo
            If (Magnetism) Then
                if (l .ne. 0) then
                Do j = fencheby+1, N_R-1, fencheby
                    equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                    equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                Enddo
                endif
            Endif


         Endif

         indx = indx + n_m + 1
      Enddo
    End Subroutine Fix_Boundary_Conditions
End Module Linear_Terms_Sphere
