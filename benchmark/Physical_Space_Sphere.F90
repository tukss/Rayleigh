#define RHSP wsp%p3b
#define DO_IDX Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define END_DO enddo; enddo; enddo
#define IDX k,r,t
#define FIELDSP wsp%p3a
Module Physical_Space_Sphere
	Use Parallel_Framework
	Use Controls
	Use ProblemSize
	Use Fourier_Transform
	Use Spectral_Derivatives
	Use Fields
	Use Diagnostics
	Use General_MPI, Only : global_max
	Use Timers
	Use Run_Parameters
	Use ClockInfo
	Use ReferenceState
	Use TransportCoefficients
	Use NonDimensionalization
	Implicit None

Contains
	Subroutine physical_space()
		Implicit None
		!write(6,*)'MAX: ', maxval(wsp%p3a(:,:,:,jphi))
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

		Call Compute_dvtheta_by_dtheta()
		Call Compute_dvphi_by_dtheta()

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
		Call ps_output(wsp%p3a, iteration,simulation_time)
		!////////////////////////////////////////////////////////////////////////


		Call Find_MyMinDT()	! Piggyback CFL communication on transposes

		
		! We are now ready to build the nonlinear terms
		Call wsp%construct('p3b')
		wsp%config = 'p3b'

		!................................
		!Nonlinear Advection
		Call StopWatch(nl_time)%startclock()

		Call Temperature_Advection()	
		Call Viscous_Heating()
		Call Momentum_Advection_Radial()
		Call Momentum_Advection_Theta()
		Call Momentum_Advection_Phi()

		If (magnetism) Then
			Call Ohmic_Heating()
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

	Subroutine Compute_dvtheta_by_dtheta()
		Implicit None
		Integer :: t, r,k

		DO_IDX
			wsp%p3a(IDX,dvtdt) = -wsp%p3a(IDX,vr)*(radius(r)*ref%dlnrho(r)-2.0d0) &
										- radius(r)*wsp%p3a(IDX,dvrdr) &
										- wsp%p3a(IDX,vtheta)*cottheta(t) &
										- wsp%p3a(IDX,dvpdp)*csctheta(t)
		END_DO

	End Subroutine Compute_dvtheta_by_dtheta

	Subroutine Compute_dvphi_by_dtheta()
		Implicit None
		Integer :: t, r,k
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			wsp%p3a(IDX,dvpdt) = radius(r)*(wsp%p3a(IDX,zvar)+wsp%p3a(IDX,dvtdp) &
										-wsp%p3a(IDX,vphi)*cottheta(t) )
		END_DO
		!$OMP END PARALLEL DO
	End Subroutine Compute_dvphi_by_dtheta

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

	Subroutine Viscous_Heating()
		Implicit None
		Integer :: t,r,k
		Real*8 :: tmp, tmp2
		Real*8, Allocatable :: htemp(:,:,:), heating_coef(:)

		Allocate(htemp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))

		! Need to optimize these loops later, but for now, let's write this in
		! easily debuggable way.

		!Contributions from E_rr, E_theta_theta 	& E_phi_phi

		!$OMP PARALLEL DO PRIVATE(t,r,k,tmp,tmp2)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
					tmp = (wsp%p3a(k,r,t,dvpdp)*csctheta(t) +wsp%p3a(k,r,t,vr) &
							+wsp%p3a(k,r,t,vtheta)*cottheta(t))*one_over_r(r)	!e_phi_phi
					tmp2 = (wsp%p3a(k,r,t,dvtdt)+wsp%p3a(k,r,t,vr))*one_over_r(r) ! e_theta_theta
					htemp(k,r,t) = wsp%p3a(k,r,t,dvrdr)*wsp%p3a(k,r,t,dvrdr)+tmp*tmp +tmp2*tmp2
					
				Enddo
			Enddo
		Enddo				
		!$OMP END PARALLEL DO

		!E_r_phi
		!$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
					tmp = (wsp%p3a(IDX,dvrdp)*csctheta(t)*one_over_r(r) &
							+wsp%p3a(IDX,dvpdr) &
							-wsp%p3a(IDX,vphi)*one_over_r(r) )*0.5d0		! e_r_phi
					
					htemp(IDX) = htemp(IDX)+tmp*tmp*2.0d0
					
				Enddo
			Enddo
		Enddo			
		!$OMP END PARALLEL DO


		!E_r_theta
		!$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
					tmp = (wsp%p3a(IDX,dvrdt)*one_over_r(r) &
							+wsp%p3a(IDX,dvtdr) &
							-wsp%p3a(IDX,vtheta)*one_over_r(r) )*0.5d0		! e_r_theta
					
					htemp(IDX) = htemp(IDX)+tmp*tmp*2.0d0
					
				Enddo
			Enddo
		Enddo			
		!$OMP END PARALLEL DO


		!E_phi_theta
		!$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
					tmp = (wsp%p3a(IDX,dvpdt) &
							+wsp%p3a(IDX,dvtdp)*csctheta(t) &
							-wsp%p3a(IDX,vphi)*cottheta(t) )*0.5d0*one_over_r(r)		! e_phi_theta
					
					htemp(IDX) = htemp(IDX)+tmp*tmp*2.0d0
					
				Enddo
			Enddo
		Enddo			
		!$OMP END PARALLEL DO


		! Allocate heating_coeff
		! Heating coeff is 2*nu/T_bar
		Allocate(heating_coef(1:N_R))
		heating_coef(1:N_R) = nu(1:N_R)*2.0d0/ref%temperature(1:N_R)

		!$OMP PARALLEL DO PRIVATE(t,r,k)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi

					wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar)+heating_coef(r)*htemp(k,r,t)
					
				Enddo
			Enddo
		Enddo				
		!$OMP END PARALLEL DO

		DeAllocate(htemp)
		DeAllocate(heating_coef)

	End Subroutine Viscous_Heating


	Subroutine Ohmic_Heating()
		Implicit None
		Integer :: t,r,k

		!We need a prefactor here for nondimensionalization

		!$OMP PARALLEL DO PRIVATE(t,r,k)
		Do t = my_theta%min, my_theta%max
			Do r = my_r%min, my_r%max
				Do k =1, n_phi
				wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar) &
									 + wsp%p3a(k,r,t,jr)*wsp%p3a(k,r,t,jr) &
									 + wsp%p3a(k,r,t,jtheta)*wsp%p3a(k,r,t,jtheta) &
									 + wsp%p3a(k,r,t,jphi)*wsp%p3a(k,r,t,jphi)
				Enddo
			Enddo
		Enddo				
		!$OMP END PARALLEL DO

	End Subroutine Ohmic_Heating

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


		! Multiply advection/coriolis pieces by rho
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,wvar) = RHSP(IDX,pvar)*ref%density(r)
		END_DO
		!OMP END PARALLEL DO	


		If (magnetism .and. lorentz_forces) Then
			! Add r_squared [JxB]_r
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,wvar)= RHSP(IDX,wvar) +r_squared(r)*ovPmEk* &
					(FIELDSP(IDX,jtheta)*FIELDSP(IDX,bphi)-FIELDSP(IDX,jphi)*FIELDSP(IDX,btheta))
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

		If (rotation) Then
			! Add - the coriolis term (part of -RHS of theta)
			! [2 z_hat cross u]_theta = -2 costheta u_phi

			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,pvar) = RHSP(IDX,pvar)- two_over_ek*costheta(t)*FIELDSP(IDX,vphi)
			END_DO
			!$OMP END PARALLEL DO
		Endif

		! Multiply advection/coriolis pieces by rho
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,pvar) = RHSP(IDX,pvar)*ref%density(r)
		END_DO
		!OMP END PARALLEL DO	

		If (magnetism .and. lorentz_forces) Then
			! Add -[JxB]_theta
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,pvar)= RHSP(IDX,pvar) &
					- ovPmEk*(FIELDSP(IDX,jphi)*FIELDSP(IDX,br)-FIELDSP(IDX,jr)*FIELDSP(IDX,bphi))
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

		! Multiply advection/coriolis pieces by rho
		!$OMP PARALLEL DO PRIVATE(t,r,k)
		DO_IDX
			RHSP(IDX,zvar) = RHSP(IDX,zvar)*ref%density(r)
		END_DO
		!OMP END PARALLEL DO		

		If (magnetism .and. lorentz_forces) Then
			! Add -[JxB]_phi
			!$OMP PARALLEL DO PRIVATE(t,r,k)
			DO_IDX
				RHSP(IDX,zvar)= RHSP(IDX,zvar) - &
					ovPmEk*(FIELDSP(IDX,jr)*FIELDSP(IDX,btheta)-FIELDSP(IDX,jtheta)*FIELDSP(IDX,br))
			END_DO
			!$OMP END PARALLEL DO
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

End Module Physical_Space_Sphere
