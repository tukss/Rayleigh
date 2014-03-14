Module Diagnostics
	Use ProblemSize
	Use Spherical_IO
	Use Fields
	Use Legendre_Polynomials, Only : gl_weights
	Implicit None
	!/////////////////////////////////////////////////////////
	!  Quantity Codes (some care must be taken to
	Integer, Parameter, Private :: V_r = 1,   V_theta = 2, V_phi = 3
  	Integer, Parameter, Private :: Temperature = 4,    Pressure = 5
	
	Integer, Parameter, Private :: v_sq = 6

	! We have some "known" outputs as well that allow us to verify that
	! the spherical_io interface is functional
	Integer, Parameter, Private :: diagnostic1 = 99, diagnostic2 = 100

	!/////////// Magnetic Outputs.  Start at 200 to organization room for hydro
	Integer, Parameter, Private :: B_r = 201, B_theta = 202, B_phi = 203
	Integer, Parameter, Private :: J_r = 204, J_theta = 205, J_phi = 206
	Integer, Parameter, Private :: B_sq = 207

	!///////////////////////////////////
	Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output

Contains

	Subroutine Initialize_Diagnostics()
		Implicit None
        Integer :: i
        Real*8 :: delr
        Real*8, Allocatable :: rweights(:), tweights(:)

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

        DeAllocate(rweights)
        DeAllocate(tweights)
		!Call Set_Spherical_IO_Integration_Weights(gl_weights, r_int_weights)
	End Subroutine Initialize_Diagnostics

	Subroutine PS_Output(buffer,iteration)
		Implicit None
		Integer, Intent(In) :: iteration
		Real*8, Intent(InOut) :: buffer(:,my_r%min:,my_theta%min:,:)
		Real*8 :: mypi
		Integer :: p,t,r
		
		If (mod(iteration,output_frequency) .eq. 0) Then
			Call Begin_Outputting(iteration)
			Allocate(qty(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))

			If (compute_q(v_r) .ne. 0) Then
				
				qty(:,:,:) = buffer(:,:,:,vr)
				!write(6,*)'Computing vr ', maxval(qty)
				Call Add_Quantity(v_r,qty)
			Endif		

			If (compute_q(v_theta) .ne. 0) Then
				
				qty(:,:,:) = buffer(:,:,:,vtheta)
				Call Add_Quantity(v_theta,qty)
			Endif		

			If (compute_q(v_phi) .ne. 0) Then
				
				qty(:,:,:) = buffer(:,:,:,vphi)
				Call Add_Quantity(v_phi,qty)
			Endif	

			If (compute_q(v_sq) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,vphi)**2
				qty(:,:,:) = qty(:,:,:)+buffer(:,:,:,vr)**2
				qty(:,:,:) = qty(:,:,:)+buffer(:,:,:,vtheta)**2
				Call Add_Quantity(v_sq,qty)
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
				qty(:,:,:) = buffer(:,:,:,br)
				Call Add_Quantity(B_r,qty)
			Endif		

			If (compute_q(B_theta) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,btheta)
				Call Add_Quantity(b_theta,qty)
			Endif		

			If (compute_q(b_phi) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,bphi)
				Call Add_Quantity(b_phi,qty)
			Endif	


			If (compute_q(J_r) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,jr)
				Call Add_Quantity(J_r,qty)
			Endif		

			If (compute_q(J_theta) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,jtheta)
				Call Add_Quantity(j_theta,qty)
			Endif		

			If (compute_q(j_phi) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,jphi)
				Call Add_Quantity(j_phi,qty)
			Endif	

			If (compute_q(b_sq) .ne. 0) Then
				qty(:,:,:) = buffer(:,:,:,bphi)**2
				qty(:,:,:) = qty(:,:,:)+buffer(:,:,:,br)**2
				qty(:,:,:) = qty(:,:,:)+buffer(:,:,:,btheta)**2
				Call Add_Quantity(b_sq,qty)
			Endif	

			Endif
			DeAllocate(qty)
			Call Complete_Output(iteration)

		Endif
	End Subroutine PS_Output

End Module Diagnostics
