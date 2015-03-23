Module Diagnostics
	Use ProblemSize
    Use Controls
	Use Spherical_IO
	Use Fields
	Use Legendre_Polynomials, Only : gl_weights
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
    Integer, Parameter, Private :: thermalE_flux_radial = 16, radial_ke = 17
    Integer, Parameter, Private :: ke_flux_radial = 18, enth_flux_radial = 19
    Integer, Parameter, Private :: buoyancy_work = 20
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
    !//////////////////////////////////

    Integer, private :: reboot_count = 0  ! Number of times diagnostics has been rebooted during this run
    
Contains

	Subroutine Initialize_Diagnostics()
		Implicit None
        Integer :: i
        Real*8 :: delr
        Real*8, Allocatable :: tweights(:)
        over_eight_pi = 1.0d0/(8.0d0*pi)
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
		  Call Initialize_Spherical_IO(radius,sintheta,rweights,tweights,costheta,my_path)	


        DeAllocate(tweights)
        DeAllocate(rweights)
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

	Subroutine Read_Output_Namelist()
		Implicit None
        Character*120 :: input_file
        input_file = Trim(my_path)//'main_input'

		! First read the main input file
		Open(unit=20, file=input_file, status="old", position="rewind")
		Read(unit=20, nml=output_namelist)
		Close(20)




	End Subroutine Read_Output_Namelist

	Subroutine PS_Output(buffer,iteration, current_time)
		Implicit None
		Integer, Intent(In) :: iteration
		Real*8, Intent(InOut) :: buffer(:,my_r%min:,my_theta%min:,:)
		Real*8, Intent(In) :: current_time
		Real*8 :: mypi, over_n_phi, tmp, tmp2, dt_by_dp, dt_by_ds, tpert
		Integer :: p,t,r
		
		If (time_to_output(iteration)) Then
			Call Begin_Outputting(iteration)
			Allocate(qty(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            over_n_phi = 1.0d0/dble(n_phi)
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
						qty(p,r,t) = ref%density(r)*ref%temperature(r)*kappa(r)*buffer(p,r,t,dtdr)
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

            ! Here we can do something like:
            !If (use_mean_vr .or. use_mean_vphi) then
            ! then we can individual quantities within - taking moments for example
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
                tpert = dt_by_ds*buffer(p,r,t,tout)+dt_by_dp*buffer(p,r,t,pvar)
                                        tpert = tpert*ref%density(r) ! This is now T'*rho_bar
                                        qty(p,r,t) = tpert*buffer(p,r,t,vr)*pressure_specific_heat
                                    Enddo
                            Enddo
                    Enddo
                    Call Add_Quantity(qty)
            Endif


			If (compute_quantity(thermalE_flux_radial)) Then
                
                qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)*buffer(1:n_phi,:,:,tout)   !pvar is temperature/r
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
                
                qty(1:n_phi,:,:) = -buffer(1:n_phi,:,:,vr)*buffer(1:n_phi,:,:,tout)   !pvar is temperature/r
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
					    	    qty(p,r,t) = ref%heating(r)*ref%density(r)*ref%temperature(r)
						    Enddo
					    Enddo
				    Enddo                
                Else
                    qty(:,:,:) = 0.0d0
                Endif
				Call Add_Quantity(qty)
			Endif	

			If (compute_quantity(diagnostic1)) Then
				mypi = acos(-1.0d0)
				Do t = my_theta%min, my_theta%max
					Do r = my_r%min, my_r%max
						Do p = 1, n_phi
							qty(p,r,t) = sin(p*2.0d0*mypi/n_phi)*(sintheta(t)**2)*radius(r)

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
							qty(p,r,t) = sin(p*4.0d0*mypi/n_phi)*(sintheta(t)*costheta(t))*radius(r)**2

						Enddo
					Enddo
				Enddo
				Call Add_Quantity(qty)
			Endif

			If (magnetism) Then
			!//////////////////// Magnetic Quantities
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
				qty(1:n_phi,:,:) = (qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,btheta)**2)*over_eight_pi
				Call Add_Quantity(qty)
			Endif	

			If (compute_quantity(zonal_me)) Then
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
				Call Add_Quantity(qty)
			Endif	

            ! Here we can do something like:
            !If (use_mean_vr .or. use_mean_vphi) then
            ! then we can individual quantities within - taking moments for example
			If (compute_quantity(merid_me)) Then
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
				Call Add_Quantity(qty)
			Endif	


			Endif
			DeAllocate(qty)
			Call Complete_Output(iteration, current_time)

		Endif
	End Subroutine PS_Output

End Module Diagnostics
