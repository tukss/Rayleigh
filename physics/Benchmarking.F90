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

    Integer, Private :: max_numt, numt_ind
    Integer, Private :: report_interval = 5
    Integer, Private :: drift_check_interval = 5 , integration_interval = 1
    Real*8, Allocatable :: time_series(:,:), time_saves(:)
Contains

    Subroutine Initialize_Benchmarking
        Implicit None
        If (benchmark_mode .eq. 1) Then
            max_numt = 1000
            integration_interval = 10
            report_interval = 5000


            ! Ideally, we override namelist values with benchmark values here

        Endif
    End Subroutine Initialize_Benchmarking

    Subroutine Benchmark_Checkup(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(:,my_r%min:,my_theta%min:,:)
        Real*8, Intent(In) :: current_time
        Real*8 :: tmp, tmp2, tmp3, time_passed, over_n_phi, shell_volume

        Integer :: i,p,t,r, num_int
        Real*8, Allocatable :: ell0_values(:,:), volume_integrals(:)	
        Real*8, Allocatable :: qty(:,:,:,:)



        !First we grab the volume-integrated quantities

        !We keep a time-series of volume-integrated quantities in memory
        over_n_phi = 1.0d0/dble(n_phi)
        If (magnetism) Then
            num_int = 6
        Else
            num_int = 3
        Endif

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
                & + buffer(1:n_phi,:,:,btheta)**2)*over_eight_pi


            ! Zonal ME
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    ! compute mean b_phi here
                    tmp = 0.0d0
                    Do p = 1, n_phi
                        tmp = tmp+buffer(p,r,t,bphi)
                    Enddo
                    tmp = tmp*over_n_phi
                    tmp = over_eight_pi*tmp**2
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
                    tmp = over_eight_pi*(tmp**2+tmp2**2)                       
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

        If (Mod(iteration,report_interval) .eq. 0) Then
            ! Generate a benchmark report
            ! Re-task the volume_integrals array
            Do i = 1, num_int
                volume_integrals(i) = SUM(time_series(:,i))  ! We want integral, not average
            Enddo
            shell_volume = four_pi*one_third*(radius(1)**3-radius(N_R)**3)
            volume_integrals = volume_integrals*(shell_volume/max_numt)
            ! Note that we'll need to do some readjustments here for 
            ! non-dimensional and dimensional runs.

            If (my_rank .eq. 0) Then
                Write(6,*)'Kinetic Energy: ', volume_integrals(1)

            Endif

        Endif

        DeAllocate(volume_integrals)
    End Subroutine Benchmark_Checkup

End Module Benchmarking
