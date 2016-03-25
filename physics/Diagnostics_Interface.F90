Module Diagnostics_Interface
    Use ProblemSize
    Use Controls
    Use Spherical_IO
    Use Fields
    Use Legendre_Polynomials, Only : gl_weights
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Use Diagnostics_Base
    Use Diagnostics_Lorentz_Forces
    Use Diagnostics_Inertial_Forces
    Implicit None


    !///////////////////////////////////
    !Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    !Real*8, Allocatable :: tmp1(:,:,:)
    !Real*8, Allocatable :: rweights(:), tweights(:)

    !//////////////////////////////////

    Integer, private :: reboot_count = 0  ! Number of times diagnostics has been rebooted during this run
    
Contains

    !//////////////////////////////////////////////////////////////////////////////////////////
    ! When entering PS_OUTPUT, the indices below may be used to reference the 4th dimension of
    ! the buffer array.  That array is dimensioned as:
    !   buffer(1:n_phi+2, my_r%min:my_r%max, my_theta%min:my_theta%max,1:nvariables)
    ! 
    ! The extra 2 in the first index is needed for the in-place FFTs.  Care should be taken
    !   to only loop over 1 to n_phi.   
    !
    ! Each index along the 4th dimension of buffer corresponds to a different variable.  
    !  These indices (left) and the variables they correspond to (right) are given below.

    ! Field variables:
    !   vr      -- radial velocity
    !   vtheta  -- theta velocity
    !   vphi    -- phi velocity
    !   tout    -- temperature or entropy (note that this is NOT tvar -- that is T/radius)
    !   pvar    -- pressure
    !   zvar    -- l(l+1)*Z/r^2  where Z is the toroidal streamfunction

    ! Radial Derivatives:
    !   dvrdr   -- d(v_r)/dr
    !   dvtdr   -- d(v_theta)/dr
    !   dvpdr   -- d(v_phi)/dr
    !   dtdr    -- d(temperature or entropy)/dr

    
    ! Theta Derivatives:
    !   dvrdt   -- d(v_r)/dtheta
    !   dvtdt   -- d(v_theta)/dtheta
    !   dvpdt   -- d(v_phi)/dtheta
    !   dtdt    -- (1/r)*d(temperature or entropy)/dtheta (<--- Note 1/r)
    !

    ! Phi Derivatives:
    !   dvrdp   --  d(v_r)/dphi
    !   dvtdp   --  d(v_theta)/dphi
    !   dvpdp   --  d(v_phi)/dphi
    !   dtdp    --  (1/r)*d(temperature or entropy)/dphi   (<--- Note 1/r)


    ! If Magnetism is On, six additional variables are present:
    !   br      -- radial magnetic field
    !   btheta  -- theta magnetic field
    !   bphi    -- phi magnetic field
    !   jr      -- radial current density
    !   jtheta  -- theta current density
    !   jphi    -- phi current density

    Subroutine PS_Output(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(In) :: current_time
        Real*8 :: mypi, over_n_phi, tmp, tmp2, tmp3, dt_by_dp, dt_by_ds, tpert

        Integer :: p,t,r, nfields, bdims(1:4), pass_num


        If (time_to_output(iteration)) Then
            Call Begin_Outputting(iteration)

            !//////////////////////////////////////////////////////////////////////////
            !We go ahead and compute two averages over the buffer:
            !   1.  The ell=0 component of each field in the buffer (requires comm)
            !   2.  The m=0   component of each field in the buffer (requires no comm)
            bdims = shape(buffer)
            nfields = bdims(4)
            Allocate(ell0_values(my_r%min:my_r%max,1:nfields))
            Allocate(m0_values(my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields))
            Call ComputeEll0(buffer,ell0_values)
            Call ComputeM0(buffer,m0_values)




            Allocate(qty(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            over_n_phi = 1.0d0/dble(n_phi)

            Do pass_num = 1, 2
            !////////////////////////
            ! All requested Shell_Average quantities are computed twice
            ! During the first pass, ell = 0 and m = 0 averages are computed
            !   (for the shell_average quantities)
            ! During the second pass, all quantities are computed, 
            !   file output is conducted, and the previously
            !   computed averages are used for moments in the shell_average output
            ! Compute_quantity returns false on the first pass for everything but shell_averages
            Call Set_Avg_Flag(pass_num)  ! This sets the averaging flag, so that all quantities or only shell averages are computed
            Call Compute_Inertial_Terms(buffer)
            If (compute_quantity(visc_flux_r)) Then
                !- v dot D |_r
                Allocate(tmp1(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
                !Radial contribution (mod rho*nu)
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = buffer(p,r,t,vr)*ref%dlnrho(r)/3.0d0+buffer(p,r,t,dvrdr)
                            qty(p,r,t) = qty(p,r,t)*buffer(p,r,t,vr)*2.0d0
                        Enddo
                    Enddo
                Enddo

                !Theta contribution (mod rho*nu)
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            tmp1(p,r,t) = (2.0d0/3.0d0)*buffer(p,r,t,vr)*ref%dlnrho(r)
                            tmp1(p,r,t) = tmp1(p,r,t)+buffer(p,r,t,dvtdr)-buffer(p,r,t,vtheta)/radius(r)
                            tmp1(p,r,t) = tmp1(p,r,t)+buffer(p,r,t,dvrdt)/radius(r)
                            tmp1(p,r,t) = tmp1(p,r,t)*buffer(p,r,t,vtheta)
                        Enddo
                    Enddo
                Enddo                


                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = qty(p,r,t)+tmp1(p,r,t)
                        Enddo
                    Enddo
                Enddo

                !phi contribution (mod rho*nu)
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            tmp1(p,r,t) = (2.0d0/3.0d0)*buffer(p,r,t,vr)*ref%dlnrho(r)
                            tmp1(p,r,t) = tmp1(p,r,t)+buffer(p,r,t,dvpdr)-buffer(p,r,t,vphi)/radius(r)
                            tmp1(p,r,t) = tmp1(p,r,t)+buffer(p,r,t,dvrdp)/radius(r)/sintheta(t)
                            tmp1(p,r,t) = tmp1(p,r,t)*buffer(p,r,t,vphi)
                        Enddo
                    Enddo
                Enddo                

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = qty(p,r,t)+tmp1(p,r,t)
                        Enddo
                    Enddo
                Enddo

                !Multiply by rho and nu
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = qty(p,r,t)*nu(r)*ref%density(r)                            
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
                DeAllocate(tmp1)
            Endif

            If (compute_quantity(v_r)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(v_theta)) Then	
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(v_phi)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(vr2)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)**2
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(vt2)) Then	
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)**2
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(vp2)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(vr3)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)**3
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(vt3)) Then	
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)**3
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(vp3)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**3
                Call Add_Quantity(qty)
            Endif	


            If (compute_quantity(rhov_r)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = buffer(p,r,t,vr)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(rhov_theta)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = buffer(p,r,t,vtheta)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif				

            If (compute_quantity(rhov_phi)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi				
                            qty(p,r,t) = buffer(p,r,t,vphi)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(temperature)) Then
                ! This is really d_by_dphi temperature/r with the current logic in Physics.F90
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = buffer(p,r,t,tout)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(pressure)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = buffer(p,r,t,pvar)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(gradt_r)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = buffer(p,r,t,dtdr)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif		

            If (compute_quantity(cond_flux_r)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = ref%density(r) &
                             & *ref%temperature(r)*kappa(r)*buffer(p,r,t,dtdr)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(zonal_ke)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        ! compute mean v_phi here
                        tmp = 0.0d0
                        Do p = 1, n_phi
                            tmp = tmp+buffer(p,r,t,vphi)
                        Enddo
                        tmp = tmp*over_n_phi
                        tmp = 0.5d0*ref%density(r)*tmp**2
                        qty(:,r,t) = tmp
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif	


            If (compute_quantity(merid_ke)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        ! compute mean v_phi here
                        tmp = 0.0d0
                        tmp2 = 0.0d0
                        Do p = 1, n_phi
                            tmp = tmp+buffer(p,r,t,vr)
                            tmp2 = tmp2+buffer(p,r,t,vtheta)
                        Enddo
                        tmp = tmp*over_n_phi
                        tmp2 = tmp2*over_n_phi
                        tmp = 0.5d0*ref%density(r)*(tmp**2+tmp2**2)                       
                        qty(:,r,t) = tmp
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(v_sq)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2
                Call Add_Quantity(qty)
            Endif	
            If (compute_quantity(radial_ke)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)**2
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = qty(p,r,t)*ref%density(r)*0.5d0
                        Enddo
                    Enddo
                Enddo  
                Call Add_Quantity(qty)
            Endif	
            If (compute_quantity(kinetic_energy)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = qty(p,r,t)*ref%density(r)*0.5d0
                        Enddo
                    Enddo
                Enddo                
                Call Add_Quantity(qty)
            Endif	


            If (compute_quantity(ke_flux_radial)) Then
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2

                qty(1:n_phi,:,:) = qty(1:n_phi,:,:)*buffer(1:n_phi,:,:,vr)

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = qty(p,r,t)*ref%density(r)*0.5d0
                        Enddo
                    Enddo
                Enddo                
                Call Add_Quantity(qty)
            Endif	




            If (compute_quantity(Enth_flux_radial)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        dt_by_ds = ref%temperature(r)/pressure_specific_heat
                        dt_by_dp = 1.0d0/pressure_specific_heat
                        Do p = 1, n_phi
                            tpert = dt_by_ds*buffer(p,r,t,tout) &
                             & + dt_by_dp*buffer(p,r,t,pvar)
                            tpert = tpert*ref%density(r) ! This is now T'*rho_bar
                            qty(p,r,t) = tpert*buffer(p,r,t,vr)*pressure_specific_heat
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif


            If (compute_quantity(thermalE_flux_radial)) Then

                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr) &
                    & *buffer(1:n_phi,:,:,tout)   
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = qty(p,r,t)*ref%density(r)*ref%temperature(r)
                        Enddo
                    Enddo
                Enddo                
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(buoyancy_work)) Then

                qty(1:n_phi,:,:) = -buffer(1:n_phi,:,:,vr) &
                 & *buffer(1:n_phi,:,:,tout)   
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = qty(p,r,t)*ref%gravity_term_s(r)
                        Enddo
                    Enddo
                Enddo                
                Call Add_Quantity(qty)
            Endif	


            If (compute_quantity(vol_heating)) Then
                If (allocated(ref%heating)) Then
                    Do t = my_theta%min, my_theta%max
                        Do r = my_r%min, my_r%max
                            Do p = 1, n_phi
                                qty(p,r,t) = ref%heating(r) &
                                 & *ref%density(r)*ref%temperature(r)
                            Enddo
                        Enddo
                    Enddo                
                Else
                    qty(:,:,:) = 0.0d0
                Endif
                Call Add_Quantity(qty)
            Endif	

            If (compute_quantity(enstrophy)) Then
                Allocate(tmp1(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            Endif

            If (compute_quantity(vort_r) .or. compute_quantity(enstrophy)) Then

                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,zvar)

                If (compute_quantity(vort_r)) Call Add_Quantity(qty)
                If (compute_quantity(enstrophy)) Then
                    tmp1 = qty**2
                Endif
            Endif

            If (compute_quantity(vort_theta) .or. compute_quantity(enstrophy)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = buffer(p,r,t,dvrdp)*csctheta(t) &
                             & - buffer(p,r,t,dvpdr) &
                             & - buffer(p,r,t,vphi)*one_over_r(r)

                        Enddo
                    Enddo
                Enddo

                If (compute_quantity(vort_theta)) Call Add_Quantity(qty)
                If (compute_quantity(enstrophy)) Then
                    tmp1 = tmp1+qty**2
                Endif

            Endif

		
            If (compute_quantity(vort_phi) .or. compute_quantity(enstrophy)) Then
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = buffer(p,r,t,vtheta)*one_over_r(r) &
                             & + buffer(p,r,t,dvtdr) &
                             & - buffer(p,r,t,dvrdt)*one_over_r(r)

                        Enddo
                    Enddo
                Enddo

                If (compute_quantity(vort_phi)) Call Add_Quantity(qty)

                If (compute_quantity(enstrophy)) Then
                    tmp1 = tmp1+qty**2
                Endif

            Endif

            If (compute_quantity(enstrophy)) Then
                Call Add_Quantity(tmp1)
                DeAllocate(tmp1)
            Endif
            
            If (compute_quantity(amom_fluct_r)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *((buffer(:,r,t,vr)-m0_values(r,t,vr))*(buffer(:,r,t,vphi) &
                            & -m0_values(r,t,vphi)))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_fluct_theta)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *((buffer(:,r,t,vtheta)-m0_values(r,t,vtheta)) & 
                            & *(buffer(:,r,t,vphi)-m0_values(r,t,vphi)))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 


            If (compute_quantity(amom_dr_r)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *(m0_values(r,t,vr)*m0_values(r,t,vphi))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_dr_theta)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *(m0_values(r,t,vtheta)*m0_values(r,t,vphi))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_mean_r)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                            & *(m0_values(r,t,vr)*Angular_Velocity)
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_mean_theta)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                            & *(m0_values(r,t,vtheta)*Angular_Velocity)
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif              
               


            !////////////////////////////////////////////////////////
            ! Diagnostics for verifying output is working
            If (compute_quantity(diagnostic1)) Then
                mypi = acos(-1.0d0)
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = sin(p*2.0d0*mypi/n_phi) &
                             & *(sintheta(t)**2)*radius(r)
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(diagnostic2)) Then
                mypi = acos(-1.0d0)
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do p = 1, n_phi
                            qty(p,r,t) = sin(p*4.0d0*mypi/n_phi) &
                             & *(sintheta(t)*costheta(t))*radius(r)**2
                        Enddo
                    Enddo
                Enddo
                Call Add_Quantity(qty)
            Endif

            !//////////////////// Magnetic Quantities
            If (magnetism) Then

                Call Compute_Lorentz_Forces(buffer)

                If (compute_quantity(B_r)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,br)
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(B_theta)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,btheta)
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(b_phi)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)
                    Call Add_Quantity(qty)
                Endif	

                If (compute_quantity(B_r2)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,br)**2
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(B_theta2)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,btheta)**2
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(b_phi2)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
                    Call Add_Quantity(qty)
                Endif	

                If (compute_quantity(J_r)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jr)
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(J_theta)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jtheta)
                    Call Add_Quantity(qty)
                Endif		

                If (compute_quantity(j_phi)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jphi)
                    Call Add_Quantity(qty)
                Endif	

                If (compute_quantity(b_sq)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
                    qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,br)**2
                    qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,btheta)**2
                    Call Add_Quantity(qty)
                Endif	

                If (compute_quantity(magnetic_energy)) Then
                    qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
                    qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,br)**2
                    qty(1:n_phi,:,:) = (qty(1:n_phi,:,:) &
                     & +buffer(1:n_phi,:,:,btheta)**2)*over_eight_pi
                    Call Add_Quantity(qty)
                Endif	

                If (compute_quantity(zonal_me)) Then
                    Do t = my_theta%min, my_theta%max
                        Do r = my_r%min, my_r%max
                            ! compute mean b_phi here
                            tmp = 0.0d0
                            Do p = 1, n_phi
                                tmp = tmp+buffer(p,r,t,bphi)
                            Enddo
                            tmp = tmp*over_n_phi
                            tmp = over_eight_pi*tmp**2
                            qty(:,r,t) = tmp
                        Enddo
                    Enddo
                    Call Add_Quantity(qty)
                Endif	


                If (compute_quantity(merid_me)) Then
                    Do t = my_theta%min, my_theta%max
                        Do r = my_r%min, my_r%max
                            tmp = 0.0d0
                            tmp2 = 0.0d0
                            Do p = 1, n_phi
                                tmp = tmp+buffer(p,r,t,br)
                                tmp2 = tmp2+buffer(p,r,t,btheta)
                            Enddo
                            tmp = tmp*over_n_phi
                            tmp2 = tmp2*over_n_phi
                            tmp = over_eight_pi*(tmp**2+tmp2**2)                       
                            qty(:,r,t) = tmp
                        Enddo
                    Enddo
                    Call Add_Quantity(qty)
                Endif	

			Endif !Magnetism
            If (pass_num .eq. 1) Call Finalize_Averages()
            Enddo !Pass_num

			DeAllocate(qty)
			Call Complete_Output(iteration, current_time)

            DeAllocate(ell0_values,m0_values)
        Endif  ! time_to_output(iteration)
    End Subroutine PS_Output

    Subroutine Read_Output_Namelist()
        Implicit None
        Character*120 :: input_file
        input_file = Trim(my_path)//'main_input'

        ! First read the main input file
        Open(unit=20, file=input_file, status="old", position="rewind")
        Read(unit=20, nml=output_namelist)
        Close(20)

    End Subroutine Read_Output_Namelist

    Subroutine Initialize_Diagnostics()
        Implicit None
        Integer :: i
        Real*8 :: delr
        
        Allocate(tweights(1:n_theta))
        tweights(:) = gl_weights(:)/2.0d0

        Allocate(rweights(1:n_r))

        If (chebyshev .or. finite_element) Then
            If (my_rank .eq. 0) Write(6,*)'Integrating using Chebyshev Quadrature'
            rweights = radial_integral_weights

        Else
            Do i = 2, n_r-1
                delr = (radius(i-1)-radius(i+1))/2.0d0
                rweights(i) = delr*radius(i)**2
            Enddo
            delr = ( radius(1)-radius(2) )/ 2.0d0
            rweights(1) = delr*radius(1)**2

            delr = (radius(n_r-1)-radius(n_r))/2.0d0
            rweights(n_r) = delr*radius(n_r)**2

            rweights = rweights/sum(rweights)

        Endif

        Call Initialize_Spherical_IO(radius,sintheta,rweights,tweights,costheta,my_path)	


        !DeAllocate(tweights)  !<---- Used to deallocate these.  We now use these for the computing the ell0 components
        !DeAllocate(rweights)
        
        !Call Set_Spherical_IO_Integration_Weights(gl_weights, r_int_weights)
    End Subroutine Initialize_Diagnostics


    Function count_digits(n) result(ndigits)
        !Counts the number of digits in the integer n
        Implicit None
        Integer, Intent(in) :: n
        Integer :: ndigits, tmp
        ndigits = 1
        tmp = abs(n)/10
            Do While( tmp .gt. 0)
                ndigits = ndigits+1
                tmp = tmp/10
            Enddo
    End Function count_digits
    Subroutine Reboot_Diagnostics(iteration,force_reboot)
        Implicit None
        ! Checks to see if a reboot file is found.  If so, reboot the diagnostics
        Integer, Intent(In) :: iteration
        Logical, Intent(In), Optional :: force_reboot

        Character*120 :: reboot_file
        Character*1 :: ndigstr
        Character*6 :: digfmt
        Character*4 :: suffix
        Logical :: reboot_now = .false.

        Integer :: ndigits

        If (present(force_reboot)) Then
            If (force_reboot) Then
                !This functionality is used in conjunction with a full reboot.
                reboot_count = 0
                reboot_now = .true.
            Endif
        Else
            If (MOD(iteration,diagnostic_reboot_interval) .eq. 0) Then

                ! Find the name of the current reboot file.
                ndigits = count_digits(reboot_count)
                Write(ndigstr,'(i1.1)') ndigits
                digfmt = '(i'//ndigstr//'.'//ndigstr//')'
                If (my_rank .eq. 0) Write(suffix,digfmt)reboot_count
                reboot_file = 'reboot_diagnostics_'//Trim(suffix)


                INQUIRE(file = reboot_file, exist = reboot_now)
                If (reboot_now) reboot_count = reboot_count+1
            Endif
        Endif


        If (reboot_now) Then     
            If (my_rank .eq. 0) Call stdout%print('Reboot file found.  Rebooting diagnostics.')   
            Call   CleanUP_Spherical_IO()
            Call   Read_Output_Namelist()
            Call Initialize_Diagnostics()
            reboot_now = .false.
        Endif
    End Subroutine Reboot_Diagnostics

End Module Diagnostics_Interface
