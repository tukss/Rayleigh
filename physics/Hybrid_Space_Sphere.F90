Module Hybrid_Space_Sphere
	Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
	Use Parallel_Framework
	Use Controls
	Use ProblemSize
	Use Legendre_Polynomials, Only : p_lm_array
	Use Legendre_Transforms, Only : Legendre_Transform 
	Use Spectral_Derivatives
	Use Fields
	Use Timers
	Use ClockInfo
	Use ReferenceState
    Use Equation_Coefficients
	Implicit None
	Real*8, Allocatable :: rho_rep(:), dlnrho_rep(:)

	Type(rmcontainer), Allocatable :: ftemp1(:), ftemp2(:),ftemp3(:)
Contains
!///// TODO
!		change /rho_rep to *one_over_rhorep

!/////////////////////////////
	!  This routine calculates terms that involve
	!      theta derivatives and loads them into
	!		 the buffer.
	Subroutine Hybrid_Init()
		! This routine will be unnecessary once the data storage
		! in this configuration is reconfigured
		Allocate(rho_rep(1:2*my_r%delta))
		rho_rep(1:my_r%delta)              = ref%density(my_r%min:my_r%max)
		rho_rep(my_r%delta+1:2*my_r%delta) = ref%density(my_r%min:my_r%max)

		Allocate(dlnrho_rep(1:2*my_r%delta))
		dlnrho_rep(1:my_r%delta)              = ref%dlnrho(my_r%min:my_r%max)
		dlnrho_rep(my_r%delta+1:2*my_r%delta) = ref%dlnrho(my_r%min:my_r%max)


	End Subroutine Hybrid_Init

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

		Call Adjust_TimeStep()


	End Subroutine rlm_spaceb

	Subroutine Velocity_Components()
		Implicit None
		Integer ::  m, mp, rmn,rmx, r, rind

		
		!Compute the velocity vield



		!vr	overwrites w	
		rmn = (vr-1)*tnr+1
		rmx = rmn+tnr-1	! striped real, then imaginary
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r-rmn+1
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,r)*ovrsq_repeated(rind)/ &
					rho_rep(rind)
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
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)/rho_rep(rind)
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
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)/rho_rep(rind)
			enddo
		Enddo

	End Subroutine Velocity_Components


	Subroutine Velocity_Derivatives()
		Implicit None
		Integer :: r,rind,rind2,rmn,rmn1,rmx,rmx1, rmn2
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
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)+ftemp2(mp)%data(:,:) !-wsp%s2a(mp)%data(:,rmn1:rmx1)  !PROBLEM
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)/rho_rep(rind)
			enddo
		Enddo

			!.... Small correction for density variation  :  - u_theta*dlnrhodr (added -u_theta/r as well here)
            ! Notice that there is a -u_theta/r term above.  These should be combined
            ! for efficiency later
		rmn1 = (vtheta-1)*tnr+1
		rmx1 = rmn1+tnr-1 
		rmn = (dvtdr-1)*tnr+1
		rmx = rmn+tnr-1

		Do mp = my_mp%min, my_mp%max
			do r = rmn, rmx
				rind = r-rmn+1
				rind2 = r-rmn+rmn1			
				wsp%s2a(mp)%data(:,r) = wsp%s2a(mp)%data(:,r)- &
                    & wsp%s2a(mp)%data(:,rind2)*(dlnrho_rep(rind)+ovr_repeated(rind))
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
			ftemp1(mp)%data(:,:) = ftemp1(mp)%data(:,:)-ftemp2(mp)%data(:,:)
			do r = rmn, rmx
				rind = r-rmn+1			
				wsp%s2a(mp)%data(:,r) = ftemp1(mp)%data(:,rind)*ovr_repeated(rind)/rho_rep(rind)
			enddo
		Enddo

			!.... Small correction for density variation  :  - u_phi*dlnrhodr
            ! .... moved -u_phi/r here as well
		rmn1 = (vphi-1)*tnr+1
		rmx1 = rmn1+tnr-1 
		rmn = (dvpdr-1)*tnr+1
		rmx = rmn+tnr-1

		Do mp = my_mp%min, my_mp%max
			do r = rmn, rmx
				rind = r-rmn+1
				rind2 = r-rmn+rmn1			
				wsp%s2a(mp)%data(:,r) = wsp%s2a(mp)%data(:,r)- &
                        &  wsp%s2a(mp)%data(:,rind2)*(dlnrho_rep(rind)+ovr_repeated(rind))
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
				wsp%s2a(mp)%data(l,rind) = l_l_plus1(l)*wsp%s2a(mp)%data(l,rind)*ovrsq_repeated(r)/rho_rep(r)
			Enddo
			Enddo
			Do r = 1,tnr
				rind = r+rmn
				rind2 = r+rmn2
				wsp%s2a(mp)%data(:,rind) = wsp%s2a(mp)%data(:,rind)-2.0d0*wsp%s2a(mp)%data(:,rind2)*ovr_repeated(r)
			Enddo
		Enddo

			!.... Small correction for density variation  :  - u_r*dlnrhodr
		rmn1 = (vr-1)*tnr+1
		rmx1 = rmn1+tnr-1 
		rmn = (dvrdr-1)*tnr+1
		rmx = rmn+tnr-1

		Do mp = my_mp%min, my_mp%max
			do r = rmn, rmx
				rind = r-rmn+1
				rind2 = r-rmn+rmn1			
				wsp%s2a(mp)%data(:,r) = wsp%s2a(mp)%data(:,r)-wsp%s2a(mp)%data(:,rind2)*dlnrho_rep(rind)
			enddo
		Enddo		


		Call d_by_dtheta(wsp%s2a,vr,dvrdt)

		
		! Convert Z to ell(ell+1) Z/r^2  (i.e. omega_r)		
		rmn = (zvar-1)*tnr+1
		rmx = rmn+tnr-1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn,rmx
				rind = r -rmn+1
				wsp%s2a(mp)%data(m:l_max,r) = l_l_plus1(m:l_max)*wsp%s2a(mp)%data(m:l_max,r)*&
					& ovrsq_repeated(rind)/rho_rep(rind)
			Enddo
		Enddo
	End Subroutine Velocity_Derivatives

	Subroutine Compute_BandJ()
		Implicit None
		Integer :: m, mp, rmn,rmx, r, rind, rmn2, roff, rind2, roff2
       
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
		roff2 = -rmn+rmn2
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Do r = rmn, rmx
				rind = r +roff 
            rind2 = r+roff2
				wsp%s2a(mp)%data(m:l_max,r) = alf_const*l_l_plus1(m:l_max) &
                    *wsp%s2a(mp)%data(m:l_max,rind2)*ovrsq_repeated(rind)
			Enddo
		Enddo     

        !Convert d2cdr2 to d2cdr2-Br (br = cl(l+1)/r^2
        rmn = (d2cdr2-1)*tnr+1
        rmx = rmn+tnr-1
        rmn2 = (cvar-1)*tnr+1     
        roff = -rmn+rmn2
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
                wsp%s2a(mp)%data(m:l_max,r) = alf_const*wsp%s2a(mp)%data(m:l_max,r)
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
				wsp%s2a(mp)%data(m:l_max,r) = alf_const*(ftemp1(mp)%data(m:l_max,rind)-ftemp2(mp)%data(m:l_max,rind))
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
	Subroutine Adjust_Emf()
		Implicit None
		Integer :: m, mp, r,rmn,rmx,roff,rind

		Call d_by_sdtheta(wsp%s2b, emfphi,ftemp1)	
		Call d_by_dphi(wsp%s2b,emftheta,ftemp2)

		Call Allocate_rlm_Field(ftemp3)
		! Copy out emf_theta before we overwrite it
		rmn = (emftheta-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				rind = r+roff		
				ftemp3(mp)%data(m:l_max,rind) = wsp%s2b(mp)%data(m:l_max,r)
			Enddo
		Enddo

		! Now for the C RHS, formed from the radial component of the curl of the emf
		! cvar overwrites emftheta
		rmn = (cvar-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				rind = r+roff
				wsp%s2b(mp)%data(m:l_max,r) = ( ftemp1(mp)%data(m:l_max,rind)- &
					& ftemp2(mp)%data(m:l_max,rind) )*over_l_l_plus1(m:l_max)
			Enddo
		Enddo		

			
		Call d_by_dphi(wsp%s2b,emfphi,ftemp2)
		! Move ftemp3 (emftheta) into emfphi's old spot
		rmn = (emfphi-1)*tnr+1
		rmx = rmn+tnr-1
		roff = -rmn+1
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			do r = rmn, rmx
				rind = r+roff		
				wsp%s2b(mp)%data(m:l_max,r)=ftemp3(mp)%data(m:l_max,rind) 
			Enddo
		Enddo
		Call d_by_sdtheta(wsp%s2b, emfphi,ftemp1)

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
		Call DeAllocate_rlm_Field(ftemp3)
		! Ensure there is no ell=0 emf  -- should I do this?
		!rmn1 = (emfr-1)    *tnr+1
		!rmn2 = (emftheta-1)*tnr+1
		!rmn3 = (emfphi-1)  *tnr+1
		!Do mp = my_mp%min, my_mp%max
		!	m = m_values(mp)
		!	if (m .eq. 0) then
		!		wsp%s2b(mp)%data(0,rmn1:rmn1+tnr-1) = 0.0d0				
		!		wsp%s2b(mp)%data(0,rmn2:rmn2+tnr-1) = 0.0d0
		!		wsp%s2b(mp)%data(0,rmn3:rmn3+tnr-1) = 0.0d0
		!	endif
		!Enddo
	End Subroutine Adjust_EMF

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

	Subroutine Allocate_rlm_Field(arr)
		Implicit None
		Type(rmcontainer), Intent(InOut), Allocatable :: arr(:)
		Integer :: mp,m


		Allocate(arr(my_mp%min:my_mp%max))
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Allocate(arr(mp)%data(m:l_max,1:tnr))
			arr(mp)%data(:,:) = 0.0d0
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
End Module Hybrid_Space_Sphere
