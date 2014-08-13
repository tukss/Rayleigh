Module Diagnostics
	Use ProblemSize
	Use Spherical_IO
	Use Fields
	Use Legendre_Polynomials, Only : gl_weights, pi
    Use ReferenceState
    Use TransportCoefficients
	Implicit None
	!/////////////////////////////////////////////////////////
	!  Quantity Codes (some care must be taken to
	Integer, Parameter, Private :: V_r = 1,   V_theta = 2, V_phi = 3
  	Integer, Parameter, Private :: Temperature = 4,    Pressure = 5
	
	Integer, Parameter, Private :: v_sq = 6, kinetic_energy = 7
    Integer, Parameter, Private :: gradt_r = 8, cond_flux_r = 9
    Integer, Parameter, Private :: zonal_ke = 10, merid_ke = 11
    Integer, Parameter, Private :: vol_heating = 12

	Integer, Parameter, Private :: rhoV_r = 13,   rhoV_theta = 14, rhoV_phi = 15
	! We have some "known" outputs as well that allow us to verify that
	! the spherical_io interface is functional
	Integer, Parameter, Private :: diagnostic1 = 99, diagnostic2 = 100

	!/////////// Magnetic Outputs.  Start at 200 to organization room for hydro
	Integer, Parameter, Private :: B_r = 201, B_theta = 202, B_phi = 203
	Integer, Parameter, Private :: J_r = 204, J_theta = 205, J_phi = 206
	Integer, Parameter, Private :: B_sq = 207, magnetic_energy=208, zonal_me = 209
    Integer, Parameter, Private :: merid_me = 210

	!///////////////////////////////////
	Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: rweights(:)
    Real*8 :: over_eight_pi
Contains

	Subroutine Initialize_Diagnostics()
		Implicit None
        Integer :: i
        Real*8 :: delr
        Real*8, Allocatable :: tweights(:)
        over_eight_pi = 1.0d0/pi
        Allocate(tweights(1:n_theta))
        tweights(:) = gl_weights(:)/2.0d0

        Allocate(rweights(1:n_r))
        Do i = 2, n_r-1
            delr = (radius(i-1)-radius(i+1))/2.0d0
            rweights(i) = delr*radius(i)**2
        Enddo
        delr = ( radius(1)-radius(2) )/ 2.0d0
        rweights(1) = delr*radius(1)**2

        delr = (radius(n_r-1)-radius(n_r))/2.0d0
        rweights(n_r) = delr*radius(n_r)**2

        rweights = rweights/sum(rweights)
		  Call Initialize_Spherical_IO(radius,sintheta,rweights,tweights,costheta)	


        DeAllocate(tweights)
		!Call Set_Spherical_IO_Integration_Weights(gl_weights, r_int_weights)
	End Subroutine Initialize_Diagnostics

	Subroutine PS_Output(buffer,iteration, current_time)
		Implicit None
		Integer, Intent(In) :: iteration
		Real*8, Intent(InOut) :: buffer(:,my_r%min:,my_theta%min:,:)
		Real*8, Intent(In) :: current_time
		Real*8 :: mypi, over_n_phi, tmp, tmp2
		Integer :: p,t,r
		
		If (mod(iteration,output_frequency) .eq. 0) Then
			Call Begin_Outputting(iteration)
			Allocate(qty(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            over_n_phi = 1.0d0/dble(n_phi)
			If (compute_q(v_r) .ne. 0) Then
				
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)
				!write(6,*)'Computing vr ', maxval(qty)
				Call Add_Quantity(v_r,qty)
			Endif		

			If (compute_q(v_theta) .ne. 0) Then
				
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)
				Call Add_Quantity(v_theta,qty)
			Endif		

			If (compute_q(v_phi) .ne. 0) Then
				
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)
				Call Add_Quantity(v_phi,qty)
			Endif	

			If (compute_q(rhov_r) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi				
				            qty(p,r,t) = buffer(p,r,t,vr)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
				Call Add_Quantity(rhov_r,qty)
			Endif		

			If (compute_q(rhov_theta) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi				
				            qty(p,r,t) = buffer(p,r,t,vtheta)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
				Call Add_Quantity(rhov_theta,qty)
			Endif				

			If (compute_q(rhov_phi) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi				
				            qty(p,r,t) = buffer(p,r,t,vphi)*ref%density(r)
                        Enddo
                    Enddo
                Enddo
				Call Add_Quantity(rhov_phi,qty)
			Endif	

			If (compute_q(temperature) .ne. 0) Then
				! This is really d_by_dphi temperature/r with the current logic in Physics.F90
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
						qty(p,r,t) = buffer(p,r,t,pvar)*radius(r)
						Enddo
					Enddo
				Enddo
				Call Add_Quantity(temperature,qty)
			Endif		

			If (compute_q(gradt_r) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
						qty(p,r,t) = buffer(p,r,t,dtdr)
						Enddo
					Enddo
				Enddo
				Call Add_Quantity(gradt_r,qty)
			Endif		

			If (compute_q(cond_flux_r) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
						qty(p,r,t) = ref%density(r)*ref%temperature(r)*kappa(r)*buffer(p,r,t,dtdr)
						Enddo
					Enddo
				Enddo
				Call Add_Quantity(cond_flux_r,qty)
			Endif	

			If (compute_q(zonal_ke) .ne. 0) Then
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
				Call Add_Quantity(zonal_ke,qty)
			Endif	

            ! Here we can do something like:
            !If (use_mean_vr .or. use_mean_vphi) then
            ! then we can individual quantities within - taking moments for example
			If (compute_q(merid_ke) .ne. 0) Then
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
				Call Add_Quantity(merid_ke,qty)
			Endif	

			If (compute_q(v_sq) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
				qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
				qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2
				Call Add_Quantity(v_sq,qty)
			Endif	

			If (compute_q(kinetic_energy) .ne. 0) Then
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
				Call Add_Quantity(kinetic_energy,qty)
			Endif	

			If (compute_q(vol_heating) .ne. 0) Then
                If (allocated(ref%heating)) Then
				    Do t = my_theta%min, my_theta%max
					    Do r = my_r%min, my_r%max
					    	Do p = 1, n_phi
					    	    qty(p,r,t) = ref%heating(r)*ref%density(r)*ref%temperature(r)
						    Enddo
					    Enddo
				    Enddo                
                Else
                    qty(:,:,:) = 0.0d0
                Endif
				Call Add_Quantity(vol_heating,qty)
			Endif	

			If (compute_q(diagnostic1) .ne. 0) Then
				mypi = acos(-1.0d0)
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
							qty(p,r,t) = sin(p*2.0d0*mypi/n_phi)*(sintheta(t)**2)*radius(r)

						Enddo
					Enddo
				Enddo
				Call Add_Quantity(diagnostic1,qty)
			Endif

			If (compute_q(diagnostic2) .ne. 0) Then
				mypi = acos(-1.0d0)
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
							qty(p,r,t) = sin(p*4.0d0*mypi/n_phi)*(sintheta(t)*costheta(t))*radius(r)**2

						Enddo
					Enddo
				Enddo
				Call Add_Quantity(diagnostic2,qty)
			Endif

			If (magnetism) Then
			!//////////////////// Magnetic Quantities
			If (compute_q(B_r) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,br)
				Call Add_Quantity(B_r,qty)
			Endif		

			If (compute_q(B_theta) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,btheta)
				Call Add_Quantity(b_theta,qty)
			Endif		

			If (compute_q(b_phi) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)
				Call Add_Quantity(b_phi,qty)
			Endif	


			If (compute_q(J_r) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jr)
				Call Add_Quantity(J_r,qty)
			Endif		

			If (compute_q(J_theta) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jtheta)
				Call Add_Quantity(j_theta,qty)
			Endif		

			If (compute_q(j_phi) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jphi)
				Call Add_Quantity(j_phi,qty)
			Endif	

			If (compute_q(b_sq) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
				qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,br)**2
				qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,btheta)**2
				Call Add_Quantity(b_sq,qty)
			Endif	

			If (compute_q(magnetic_energy) .ne. 0) Then
				qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
				qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,br)**2
				qty(1:n_phi,:,:) = (qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,btheta)**2)*over_eight_pi
				Call Add_Quantity(magnetic_energy,qty)
			Endif	

			If (compute_q(zonal_me) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
                        ! compute mean v_phi here
                        tmp = 0.0d0
						Do p = 1, n_phi
						    tmp = tmp+buffer(p,r,t,bphi)
						Enddo
                        tmp = tmp*over_n_phi
                        tmp = over_eight_pi*tmp**2
                        qty(:,r,t) = tmp
					Enddo
				Enddo
				Call Add_Quantity(zonal_me,qty)
			Endif	

            ! Here we can do something like:
            !If (use_mean_vr .or. use_mean_vphi) then
            ! then we can individual quantities within - taking moments for example
			If (compute_q(merid_me) .ne. 0) Then
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
                        ! compute mean v_phi here
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
				Call Add_Quantity(merid_me,qty)
			Endif	


			Endif
			DeAllocate(qty)
			Call Complete_Output(iteration, current_time)

		Endif
	End Subroutine PS_Output

End Module Diagnostics
