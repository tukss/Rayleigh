!This module allows benchmarking measurements to be carried out in-situ
!Benchmark 1: Christensen et al. Model X (hydro)
!Benchmark 2: Christensen et al. Model Y (mhd)
!Benchmark 3: Jones et al. hydro + steady
!Benchmark 4: Jones et al. mhd + steady
!Benchmark 5: Jones et al. mhd + unsteady (In Development)

!  This module is intended for checking ACCURACY of Rayleigh's results.
!  This module has NOT been programmed with EFFICIENCY in mind (much).
Module Benchmarking
    Use ProblemSize
    Use Controls
    Use Spherical_IO
    Use Fields
    Use Legendre_Polynomials, Only : gl_weights
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Implicit None

    Integer, Private :: nobs
    Integer, Private :: max_numt, numt_ind, global_count, num_int
    Integer, Private :: report_interval = 500
    Integer, Private :: drift_check_interval = 5 , integration_interval = 50
    Real*8 :: mag_factor
    Real*8, Allocatable :: time_series(:,:), time_saves(:)

    Integer :: r_one, r_two, theta_one, theta_two
    Integer :: r1_id, r2_id, theta1_id, theta2_id
    Logical :: have_r_one, have_r_two, have_theta_one, have_theta_two
    Logical :: have_r1t1, have_r1t2, have_r2t1, have_r2t2

Contains

    Subroutine Initialize_Benchmarking
        Implicit None
        Integer :: p, your_r_min, your_r_max, your_theta_min, your_theta_max
        global_count = 0
        numt_ind = 1
        If (benchmark_mode .eq. 1) Then

            max_numt = 1000
            integration_interval = 10 !0
            report_interval = 100 !000
            mag_factor = 1.0d0/(2*ekman_number*magnetic_prandtl_number)

            ! Ideally, we override namelist values with benchmark values here

        Endif

        !NOTE:  for dimensional anelastic runs, use mag_factor = over_eight_pi

        ! We need to bracket the equator and the radial center of the domain

        !n_theta should always be even
        theta_one = n_theta/2
        theta_two = theta_one+1

        !n_r is always even in Chebyshev runs - assume it's even for now
        r_one = n_r/2
        r_two = r_one+1

        !///////////////////////////////////////////////////////////
        ! Figure out if I own one of the desired strips
        have_r_one = .false.
        have_r_two = .false.
        have_theta_two = .false.
        have_theta_one = .false.

    
        have_r1t1 = .false.
        have_r1t2 = .false.
        have_r2t1 = .false.
        have_r2t2 = .false.


        If ((r_one .le. my_r%max) .and. (r_one .ge. my_r%min) ) have_r_one = .true.
        If ((r_two .le. my_r%max) .and. (r_two .ge. my_r%min) ) have_r_two = .true.
        If ((theta_one .le. my_theta%max) .and. (theta_one .ge. my_theta%min) ) have_theta_one = .true.
        If ((theta_two .le. my_theta%max) .and. (theta_two .ge. my_theta%min) ) have_theta_two = .true.


        If (have_r_one) Then
            if (have_theta_one) have_r1t1 = .true.
            if (have_theta_two) have_r1t2 = .true.
        Endif

        If (have_r_two) Then
            if (have_theta_one) have_r2t1 = .true.
            if (have_theta_two) have_r2t2 = .true.
        Endif



        If (my_rank .eq. 0) Then 
            Do p = 0, npcol-1
                your_r_min = pfi%all_1p(p)%min
                your_r_max = pfi%all_1p(p)%max 
                If ((r_one .le. your_r_max) .and. (r_one .ge. your_r_min) ) r1_id = p
                If ((r_two .le. your_r_max) .and. (r_two .ge. your_r_min) ) r2_id = p
            Enddo
        Endif

        If (my_row_rank .eq. 0) Then 
            Do p = 0, nprow-1
                your_theta_min = pfi%all_2p(p)%min
                your_theta_max = pfi%all_2p(p)%max 
                If ((theta_one .le. your_theta_max) .and. & 
                    & (theta_one .ge. your_theta_min) ) theta1_id = p
                If ((theta_two .le. your_theta_max) .and. &
                    & (theta_two .ge. your_theta_min) ) theta2_id = p
            Enddo
        Endif


        If (magnetism) Then
            num_int = 6
            nobs = 2
        Else
            num_int = 3
            nobs = 3
        Endif
        Allocate(time_series(1:max_numt,1:num_int))
        Allocate(time_saves(1:max_numt))
    End Subroutine Initialize_Benchmarking

    Subroutine Benchmark_Checkup(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(:,my_r%min:,my_theta%min:,:)
        Real*8, Intent(In) :: current_time
        Real*8 :: tmp, tmp2, tmp3, time_passed, over_n_phi, shell_volume

        Integer :: i,p,t,r
        Real*8, Allocatable :: ell0_values(:,:), volume_integrals(:)	
        Real*8, Allocatable :: qty(:,:,:,:)

        If (mod(iteration,integration_interval) .eq. 0) Then
            !Write(6,*)'Integrating!'
            !First we grab the volume-integrated quantities

            !We keep a time-series of volume-integrated quantities in memory
            over_n_phi = 1.0d0/dble(n_phi)


            Allocate(qty(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:num_int))
            
            qty(:,:,:,:) = 0.0d0


            !KE
            qty(1:n_phi,:,:,1) = buffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:,1) = qty(1:n_phi,:,:,1)+buffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:,1) = qty(1:n_phi,:,:,1)+buffer(1:n_phi,:,:,vtheta)**2
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do p = 1, n_phi
                        qty(p,r,t,1) = qty(p,r,t,1)*ref%density(r)*0.5d0
                    Enddo
                Enddo
            Enddo                


            !Zonal KE        
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    ! compute mean v_phi here
                    tmp = 0.0d0
                    Do p = 1, n_phi
                        tmp = tmp+buffer(p,r,t,vphi)
                    Enddo
                    tmp = tmp*over_n_phi
                    tmp = 0.5d0*ref%density(r)*tmp**2
                    qty(:,r,t,2) = tmp
                Enddo
            Enddo

            !Meridional KE
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
                    qty(:,r,t,3) = tmp
                Enddo
            Enddo

            If (magnetism) Then

                ! ME
                qty(1:n_phi,:,:,4) = buffer(1:n_phi,:,:,bphi)**2
                qty(1:n_phi,:,:,4) = qty(1:n_phi,:,:,4)+buffer(1:n_phi,:,:,br)**2
                qty(1:n_phi,:,:,4) = (qty(1:n_phi,:,:,4) &
                    & + buffer(1:n_phi,:,:,btheta)**2)


                ! Zonal ME
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        ! compute mean b_phi here
                        tmp = 0.0d0
                        Do p = 1, n_phi
                            tmp = tmp+buffer(p,r,t,bphi)
                        Enddo
                        tmp = tmp*over_n_phi
                        tmp = tmp**2
                        qty(:,r,t,5) = tmp
                    Enddo
                Enddo

                ! Meridional ME
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
                        tmp = (tmp**2+tmp2**2)                       
                        qty(:,r,t,6) = tmp
                    Enddo
                Enddo

            Endif ! (magnetism)

            Allocate(ell0_values(my_r%min:my_r%max,1:num_int))
            Allocate(volume_integrals(1:num_int))
            Call ComputeEll0(qty,ell0_values)   ! Requires communication across row
            DeAllocate(qty)

            Call Compute_Radial_Average(ell0_values,volume_integrals) !Requires communication across column
            DeAllocate(ell0_values)


            
            

            !Now, add these into the time_series and update the counter
            If (numt_ind .gt. max_numt) numt_ind = 1
            time_series(numt_ind,1:num_int) = volume_integrals(1:num_int)
            time_saves(numt_ind) = current_time
            numt_ind = numt_ind+1
            global_count = global_count+1

            !////////////////////////////////////////////////////////
            Call Grab_Slices(buffer)


            If (Mod(iteration,report_interval) .eq. 0) Then
                ! Generate a benchmark report
                ! Re-task the volume_integrals array
                If (global_count .ge. max_numt) Then
                    Do i = 1, num_int
                        volume_integrals(i) = SUM(time_series(:,i))/max_numt 
                    Enddo
                Else
                    Do i = 1, num_int
                        volume_integrals(i) = SUM(time_series(1:global_count,i))/global_count  
                    Enddo
                Endif
                shell_volume = four_pi*one_third*(radius(1)**3-radius(N_R)**3)
                !volume_integrals = volume_integrals*shell_volume ! We want integral, not average
                ! Note that we'll need to do some readjustments here for 
                ! non-dimensional and dimensional runs.

                If (my_rank .eq. 0) Then
                    Write(6,*)'Kinetic Energy: ', volume_integrals(1)
                    Write(6,*)'Magnetic Energy: ', volume_integrals(4)*mag_factor
                Endif

            Endif

            DeAllocate(volume_integrals)
        Endif
    End Subroutine Benchmark_Checkup

    Subroutine Grab_Slices(inbuff)
        Real*8, Intent(In) :: inbuff(1:,my_r%min:,my_theta%min:,:)
        Integer :: i

        ! At the end of this routine, rank zero will have
        ! All four slices in memory 
        Allocate(slices(1:n_phi,1:nobs,4))
        Allocate(obs_inds(1:nobs))
        obs_inds(1) = vr
        obs_inds(2) = vphi
        obs_inds(3) = tvar
        if (magnetism) obs_inds(4) = btheta

        If (my_rank .eq. 0) Then
            If (.not. have_r1t1) Then
                !IReceive
            Endif
            If (.not. have_r1t2) Then
                !IReceive
            Endif
            If (.not. have_r2t1) Then
                !IReceive
            Endif
            If (.not. have_r2t2) Then
                !IReceive
            Endif

            !Wait on the receives to complete
            Allocate(avg_slice(1:n_phi,1:nobs))        
            avg_slice(:,:) = 0.0d0
            Do j = 1, 4
                Do i = 1, nobs
                    avg_slice(1:n_phi,i) = avg_slice(1:n_phi,i)+slices(1:n_phi,i,j)
                Enddo
            Enddo
            avg_slice = avg_slice*0.25d0 ! convert to an average
        Endif


        If (have_r1t1) Then
            Do i = 1, nobs
                slices(1:n_phi,i,1) = inbuff(1:n_phi,r_one,theta_one,obs_inds(i))
            Enddo

            If (my_rank .ne. 0) Then
                !Note --- check the interface for ISEND
                Call ISend(slices(1,1,1), sirq11,ncount, dest = 0, tag = r11_tag, grp = pfi%gcomm)
            Endif
        Endif

        If (have_r1t2) Then
            Do i = 1, nobs
                slices(1:n_phi,i,2) = inbuff(1:n_phi,r_one,theta_two,obs_inds(i))
            Enddo
            If (my_rank .ne. 0) Then
                ! ISend
            Endif
        Endif

        If (have_r2t1) Then
            Do i = 1, nobs
                slices(1:n_phi,i,3) = inbuff(1:n_phi,r_two,theta_one,obs_inds(i))
            Enddo
            If (my_rank .ne. 0) Then
                ! ISend
            Endif
        Endif

        If (have_r2t2) Then
            Do i = 1, nobs
                slices(1:n_phi,i,4) = inbuff(1:n_phi,r_two,theta_two,obs_inds(i))
            Enddo

            If (my_rank .ne. 0) Then
                ! ISend
            Endif
        Endif

        If (my_rank .ne. 0) Then
            If (have_r1t1) ! Wait
            If (have_r1t2) ! Wait
            If (have_r2t1) ! Wait
            If (have_r2t2) ! Wait
        Endif


    End Subroutine
End Module Benchmarking
