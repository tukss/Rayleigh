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

    Integer, Private :: nobs, msymm
    Integer, Private :: max_numt, numt_ind, global_count, num_int
    Integer, Private :: report_interval = 500
    Integer, Private :: drift_check_interval = 5 , integration_interval = 50
    Real*8 :: mag_factor
    Real*8, Allocatable :: time_series(:,:), time_saves(:)

    Integer :: r_one, r_two, theta_one, theta_two
    Integer :: strip_owners(1:4), btags(1:4), tvals(1:4), rvals(1:4)
    Integer :: num_strips = 4  !The number of strips we intend to average over
    Logical :: have_strip(1:4) 
    Real*8, Allocatable :: strips(:,:) ! equatorial-mid-shell strips used for tracking pattern drift rate
    Real*8, Allocatable :: observations(:)
Contains

    Subroutine Initialize_Benchmarking
        Implicit None

        Integer :: p,i
        Integer :: your_row_rank, your_col_rank
        Integer :: your_r_min, your_r_max
        Integer :: your_theta_min, your_theta_max
        Logical :: have_r_one, have_r_two, have_theta_one, have_theta_two

        global_count = 0
        numt_ind = 1
        If (magnetism) Then
            num_int = 6
            nobs = 4
        Else
            num_int = 3
            nobs = 3
        Endif
     
        If (benchmark_mode .eq. 1) Then

            max_numt = 1000
            integration_interval = 5 !0
            report_interval = 50 !000
            mag_factor = 1.0d0/(2*ekman_number*magnetic_prandtl_number)
            msymm = 4
            ! Ideally, we override namelist values with benchmark values here

        Endif



        Allocate(time_series(1:max_numt,1:num_int))  ! Possible that only rank 0 needs these -- check
        Allocate(time_saves(1:max_numt))   


        !NOTE:  for dimensional anelastic runs, use mag_factor = over_eight_pi

        ! We need to bracket the equator and the radial center of the domain

        Do i = 1, 4
            btags(i)=100+i  !MPI Tags for each of the four strips
        Enddo

        !n_theta should always be even
        theta_one = n_theta/2
        theta_two = theta_one+1

        !n_r is always even in Chebyshev runs - assume it's even for now
        r_one = n_r/2
        r_two = r_one+1

        rvals(1:2) = r_one
        rvals(3:4) = r_two
        tvals(1) = theta_one
        tvals(2) = theta_two
        tvals(3) = theta_one
        tvals(4) = theta_two

        !///////////////////////////////////////////////////////////
        ! Figure out if I own one of the desired strips
        have_r_one = .false.
        have_r_two = .false.
        have_theta_two = .false.
        have_theta_one = .false.
        have_strip(1:4) = .false.
    

        If ((r_one .le. my_r%max) .and. (r_one .ge. my_r%min) ) have_r_one = .true.
        If ((r_two .le. my_r%max) .and. (r_two .ge. my_r%min) ) have_r_two = .true.
        If ((theta_one .le. my_theta%max) .and. (theta_one .ge. my_theta%min) ) have_theta_one = .true.
        If ((theta_two .le. my_theta%max) .and. (theta_two .ge. my_theta%min) ) have_theta_two = .true.


        If (have_r_one) Then
            if (have_theta_one) have_strip(1) = .true.
            if (have_theta_two) have_strip(2) = .true.
        Endif

        If (have_r_two) Then
            if (have_theta_one) have_strip(3) = .true.
            if (have_theta_two) have_strip(4) = .true.
        Endif



        If (my_rank .eq. 0) Then 
            Allocate(strips(1:n_phi,1:nobs))
            Allocate(observations(1:nobs))
            Do p = 0, (npcol*nprow)-1
		        your_row_rank = mod(p,nprow)
		        your_col_rank = p/nprow

                your_r_min = pfi%all_1p(your_col_rank)%min
                your_r_max = pfi%all_1p(your_col_rank)%max 
                If ((r_one .le. your_r_max) .and. (r_one .ge. your_r_min) ) Then
                    your_theta_min = pfi%all_2p(your_row_rank)%min
                    your_theta_max = pfi%all_2p(your_row_rank)%max 
                    If ((theta_one .le. your_theta_max) .and. &
                        &  (theta_one .ge. your_theta_min) ) Then

                        strip_owners(1) = p

                    Endif
                    If ((theta_two .le. your_theta_max) .and. &
                        & (theta_two .ge. your_theta_min) ) Then

                        strip_owners(2) = p

                    Endif
                Endif
                If ((r_two .le. your_r_max) .and. (r_two .ge. your_r_min) ) Then
                    your_theta_min = pfi%all_2p(your_row_rank)%min
                    your_theta_max = pfi%all_2p(your_row_rank)%max 
                    If ((theta_one .le. your_theta_max) .and. &
                        &  (theta_one .ge. your_theta_min) ) Then

                        strip_owners(3) = p

                    Endif
                    If ((theta_two .le. your_theta_max) .and. &
                        & (theta_two .ge. your_theta_min) ) Then

                        strip_owners(4) = p

                    Endif
                Endif


            Enddo
        Endif



    End Subroutine Initialize_Benchmarking

    Subroutine Benchmark_Checkup(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,:)
        Real*8, Intent(In) :: current_time
        Real*8 :: tmp, tmp2, tmp3, time_passed, over_n_phi, shell_volume

        Integer :: i,p,t,r
        Real*8, Allocatable :: ell0_values(:,:), volume_integrals(:)	
        Real*8, Allocatable :: qty(:,:,:,:)
        Character*120 :: strip_file = 'strips'


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
            Call Assemble_Strips(buffer)
            If (my_rank .eq. 0) Then
                Call Point_Observations()
                !Call Calculate_Drift()
            Endif

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
                    Call Write_Array(strips,strip_file)
                    Write(6,*)maxval(strips(:,1)), maxval(strips(:,2)), maxval(strips(:,3)), maxval(strips(:,4))
                    Write(6,*)observations
                Endif

            Endif

            DeAllocate(volume_integrals)
        Endif
    End Subroutine Benchmark_Checkup

    Subroutine Assemble_Strips(inbuff)
        Real*8, Intent(In) :: inbuff(1:,my_r%min:,my_theta%min:,:)
        Real*8, Allocatable :: all_strips(:,:,:)
        Integer :: i,j, indst(3), ndata
        Integer :: rirqs(1:4), sirqs(1:4)

        Integer, Allocatable :: obs_inds(:)

        ! At the end of this routine, rank zero will have
        ! All four slices in memory 
        Allocate(all_strips(1:n_phi,1:nobs,4))
        Allocate(obs_inds(1:nobs))
        obs_inds(1) = vr
        obs_inds(2) = vphi
        obs_inds(3) = tvar
        if (magnetism) obs_inds(4) = btheta

        indst(:) = 1
        ndata = n_phi*nobs !number of data points communicated during each send

        !Rank 0 Post receives as appropriate

        If (my_rank .eq. 0) Then
            Do j = 1, num_strips
                If (.not. have_strip(j)) Then
                    indst(3) = j
                    Call Ireceive(all_strips, rirqs(j), n_elements = ndata,source= strip_owners(j), &
                        &  tag=btags(j),grp = pfi%gcomm, indstart = indst)
                Endif
            Enddo
        Endif

        ! Ranks that own strips send them to rank 0 as appropriate
        Do j = 1, num_strips
            If (have_strip(j)) Then
                Do i = 1, nobs
                    all_strips(1:n_phi,i,1) = inbuff(1:n_phi,rvals(j),tvals(j),obs_inds(i))
                Enddo

                If (my_rank .ne. 0) Then
                    ! If this was rank 0, we only need to copy the strip into the strips array
                    indst(3) = j
                    Call ISend(all_strips, sirqs(j),n_elements = ndata, dest = 0, &
                        &  tag = btags(j), grp = pfi%gcomm)
                Endif                
            Endif
        Enddo

        !Next, we wait on sends and receives to complete and compute the average
        If (my_rank .ne. 0) Then
            Do j = 1, num_strips
                If (have_strip(j)) Then
                    Call IWait(sirqs(j))
                Endif
            Enddo
        Else
            Do j = 1, num_strips
                If (.not. have_strip(j)) Then
                    Call IWait(rirqs(j))
                Endif
            Enddo


                    
            strips(:,:) = 0.0d0
            Do j = 1, num_strips
                Do i = 1, nobs
                    strips(1:n_phi,i) = strips(1:n_phi,i)+all_strips(1:n_phi,i,j)
                Enddo
            Enddo
            strips = strips*(1.0d0/num_strips) ! convert to an average

        Endif
        DeAllocate(all_strips, obs_inds)

    End Subroutine Assemble_Strips

    Subroutine Point_Observations()
        Implicit None
        !Conduct point-wise observations of T, u_phi, and b_theta
        ! Observations taken where vr = 0 and dvr/dphi > 0
        Integer :: i, j, xind, im1
        Real*8 :: vrlast, vrnext, dvr, test,x, y2, y1, slope
        Real*8, Allocatable :: vsave(:,:)
        Allocate(vsave(1:msymm,1:nobs))
        xind = 1
        im1 = n_phi
        Do i = 1, n_phi
            vrnext = strips(i,1)
            vrlast = strips(im1,1)
            dvr = vrnext-vrlast
            test = vrnext*vrlast
            If ( (test .lt. 0) .and. (xind .le. msymm) ) Then
                ! We have crossed a zero point
                If (dvr .gt. 0) Then
                    !dvr/dphi is positive at this zero point
                    x = -vrlast/dvr  ! zero crossing (relative to im1 = 0)
                    Do j = 1, nobs
                        ! Linearly interpolate to find the values at the zero crossing of vr
                        ! Interpolated vr is zero by definition -- good sanity check to save
                        y1 = strips(im1,j)
                        y2 = strips(i,j)
                        slope = y2-y1
                        vsave(xind,j) = slope*x+y1
                    Enddo
                    xind = xind+1
                Endif
            Endif
            im1 = i
        Enddo

        Do i = 1, nobs 
            observations(i) = SUM(vsave(:,i))/msymm
        Enddo
        Write(6,*)xind, msymm
        DeAllocate(vsave)
    End Subroutine Point_Observations

	Subroutine Write_Array(arr,filename)
		Implicit None
		Character*120, Optional, Intent(In) :: filename
		Integer :: i,j,sig = 314, nx,ny
        Real*8, Intent(In) :: arr(1:,1:)
        nx = size(arr,1)
        ny = size(arr,2)

		Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
        Write(15)sig
		Write(15)nx
        Write(15)ny
		Write(15)((arr(i,j),i=1,nx),j = 1, ny)

		Close(15)

	End Subroutine Write_Array
End Module Benchmarking
