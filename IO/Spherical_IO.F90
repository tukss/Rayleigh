Module Spherical_IO
	Use Parallel_Framework
	Use SendReceive

	! This module contains routines for outputing spherical data as:
	! 1. Slices of sphere
	! 2. Phi-Averages over sphere (f(r,theta))
	! 3. Phi/theta Averages over sphere ((f(r))
	! 4. Volume averages over sphere

	!////////////////////////////////////////////
		Integer, Parameter :: nqmax=300, nshellmax=100
		Integer :: shellavg_values(1:nqmax)=-1, globalavg_values(1:nqmax)=-1
		Integer :: shellslice_values(1:nqmax) =-1, shellslice_levels(1:nshellmax)=-1, azavg_values(1:nqmax)=-1
		Integer :: full3d_values(1:nqmax) = -1
		Integer :: output_frequency, eng_frequency, checkpoint_frequency
		Real(8) :: checkpoint_time
      Namelist /output_namelist/shellavg_values, globalavg_values, &
			& shellslice_values, shellslice_levels, azavg_values, output_frequency, &
			& checkpoint_frequency, checkpoint_time, eng_frequency, full3d_values

		Character*120 :: run_dir
		Integer :: itmax, restart_checkpoint
		Real(8) :: tend, wtlimit
		Namelist /input_namelist/run_dir, itmax, restart_checkpoint, &
			& tend, wtlimit

		Integer, Allocatable :: compute_q(:)
		Integer, Allocatable, Private :: compute_azav(:), compute_shellav(:), compute_globav(:)
		Integer, Allocatable, Private :: compute_shell(:), shell_levels(:)
		Integer, Allocatable, Private :: step_zero(:), step_one(:), step_two(:), step_three(:), step_five(:)
		Integer, Allocatable, Private :: qvals_azav(:), qvals_shellav(:), qvals_globav(:), qvals_shell(:)
		Integer, Allocatable, Private :: qvals_3d(:)
		Integer, Private :: nq_max, nq_azav, nq_shellav, nq_globav, azav_ind, shellav_ind
		Integer, Private :: first_azav_q, last_azav_q, first_shellav_q, last_shellav_q
		Integer, Private :: first_globav_q, last_globav_q, globav_ind, nq_shell, nshell_levels
		Integer, Private :: first_shell_q, last_shell_q, shell_ind, my_nshells, nshell_r_ids
		Integer, Private :: shell_data_written, nq_3d
		Integer, Private :: az_avg_tag = 54, shell_avg_tag = 55, global_avg_tag = 56, shell_slice_tag = 57, full_3d_tag = 58
		Integer, Private, Allocatable :: my_shell_levs(:), have_shell(:), shell_r_ids(:), my_shell_ind(:)
		Integer, Private, Allocatable :: nshells_at_rid(:)
		Real*8, Private, Allocatable :: circumference(:,:), qty(:,:,:), f_of_r_theta(:,:)
		Real*8, Private, Allocatable :: azav_outputs(:,:,:), f_of_r(:), rdtheta_total(:)
		Real*8, Private, Allocatable :: shellav_outputs(:,:), globav_outputs(:), shell_slice_outputs(:,:,:,:)
		Real*8, Private :: da_total, int_vol, int_dphi, int_rsquared_dr, int_sintheta_dtheta
		Real*8, Private, Allocatable :: sintheta_dtheta(:), rsquared_dr(:)
		Character*6, Public :: i_ofmt = '(i8.8)', i_pfmt = '(i5.5)'
		integer :: output_ireq(1), output_status(1)

		!////////////////////////////////////////////////////////////////////
		Logical :: start_az_average, start_shell_average, start_shell_slice, start_global_average

        Integer :: io_node = 0

	Integer, Private :: current_iteration
	!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	! We store some redundant information just to keep the IO level as independent of the Physics level as possible
	! These variables are all private.
	Real*8, Allocatable,Private :: theta_integration_weights(:), r_integration_weights(:)
	Integer, Private :: my_theta_min, my_theta_max, my_rmin, my_rmax
	Integer, Private :: nr, nphi, ntheta
	Integer, Private :: nproc1, nproc2, myid, nproc, my_row_rank, my_column_rank, my_nr
	Real*8, Allocatable, Private :: radius(:),sintheta(:),costheta(:) 
Contains

    Subroutine Begin_Outputting(iter)
		  Implicit None
		  Integer, Intent(In) :: iter
		  current_iteration = iter
        start_az_average = .true.
        start_shell_average = .true.
        start_shell_slice = .true.
        start_global_average = .true.
        Allocate(f_of_r_theta(my_rmin:my_rmax,my_theta_min:my_theta_max))
        Allocate(f_of_r(my_rmin:my_rmax))
    End Subroutine Begin_Outputting

	Subroutine Initialize_Spherical_IO(rad_in,sintheta_in, rw_in, tw_in, costheta_in)
		Implicit None
		Integer :: i,j, k, ind, ilocal, pcount, your_rmin, your_rmax
		Integer, Allocatable :: ptemp(:), ptemp2(:), ptemp3(:)
		Real*8, Intent(In) :: rad_in(:), sintheta_in(:), rw_in(:), tw_in(:), costheta_in(:)
		! Handles output bookkeeping

		my_theta_min = pfi%my_2p%min
		my_theta_max = pfi%my_2p%max
		my_rmin = pfi%my_1p%min
		my_rmax = pfi%my_1p%max
		my_nr = pfi%my_1p%delta
        nproc  = pfi%gcomm%np
		nproc1 = pfi%ccomm%np		! processor's per column (radius split among these)
		nproc2 = pfi%rcomm%np      ! processor's per row  (theta split among these)
		myid   = pfi%gcomm%rank
		my_row_rank = pfi%rcomm%rank
		my_column_rank = pfi%ccomm%rank
		nphi = pfi%n3p
		ntheta = pfi%n2p
		nr = pfi%n1p

		Allocate(sintheta(1:ntheta))
		Allocate(costheta(1:ntheta))
		sintheta(:) = sintheta_in(:)
		costheta(:) = costheta_in(:)
		Allocate(radius(1:nr))
		radius(:) = rad_in(:)

        Allocate(r_integration_weights(1:nr))
        Allocate(theta_integration_weights(1:ntheta))

        r_integration_weights(:) = rw_in(:)
        theta_integration_weights(:) = tw_in(:)		

		nq_max = nqmax
		nq_azav = 0
		nq_shellav = 0
		nq_globav = 0
		nq_shell = 0
		nq_3d = 0
		nshell_levels = 0

		do i = 1, nq_max
			if (shellslice_values(i) .gt. 0) nq_shell = nq_shell+1
			if (azavg_values(i) .gt. 0)       nq_azav = nq_azav+1
			if (globalavg_values(i) .gt. 0) nq_globav = nq_globav+1
			if (shellavg_values(i) .gt. 0) nq_shellav = nq_shellav+1
			if (full3d_values(i) .gt. 0) nq_3d = nq_3d+1
		enddo

		do i = 1, nshellmax
			if (shellslice_levels(i) .gt. 0) then 
				nshell_levels = nshell_levels+1
			endif
		enddo
		Allocate(qvals_azav(1:nq_azav))
		Allocate(qvals_shellav(1:nq_shellav))
		Allocate(qvals_globav(1:nq_globav))
		Allocate(qvals_shell(1:nq_shell))
		Allocate(qvals_3d(1:nq_3d))
		Allocate(shell_levels(1:nshell_levels))

		do i = 1, nq_shell
			qvals_shell(i) = shellslice_values(i)
		enddo

		do i = 1, nq_azav
			qvals_azav(i) = azavg_values(i)
		enddo

		do i = 1, nq_shellav
			qvals_shellav(i) = shellavg_values(i)
		enddo

		do i = 1, nq_globav
			qvals_globav(i) = globalavg_values(i)
		enddo

		do i = 1, nq_3d
			qvals_3d(i) = full3d_values(i)
		enddo

		do i = 1, nshell_levels
			shell_levels(i) = shellslice_levels(i)
		enddo

		If (nq_shell .gt. 0) Then
			first_shell_q = qvals_shell(1)
			last_shell_q = qvals_shell(nq_shell)
		Endif
		Allocate(compute_shell(1:nq_max))
		compute_shell(1:nq_max) = 0

		If (nq_azav .gt. 0) Then
			first_azav_q = qvals_azav(1)
			last_azav_q = qvals_azav(nq_azav)
		Endif
		Allocate(compute_azav(1:nq_max))
		compute_azav(1:nq_max) = 0

		If (nq_shellav .gt. 0) Then
			first_shellav_q = qvals_shellav(1)
			last_shellav_q = qvals_shellav(nq_shellav)
		Endif
		Allocate(compute_shellav(1:nq_max))
		compute_shellav(1:nq_max) = 0

		If (nq_globav .gt. 0) Then
			first_globav_q = qvals_globav(1)
			last_globav_q = qvals_globav(nq_globav)
		Endif
		Allocate(compute_globav(1:nq_max))
		compute_globav(1:nq_max) = 0			

		Allocate(compute_q(1:nq_max))
		compute_q(:) = 0

		Allocate(step_zero(1:nq_max))
		Allocate(step_one(1:nq_max))
		Allocate(step_two(1:nq_max))
		Allocate(step_three(1:nq_max))
		Allocate(step_five(1:nq_max))


		



		
		step_zero(1:nq_max)  = 0
		step_one(1:nq_max)   = 0
		step_two(1:nq_max)   = 0
		step_three(1:nq_max) = 0
		step_five(1:nq_max)  = 0	! For 3-D Output

		Do i = 1, nq_globav
			ind = qvals_globav(i)
			compute_globav(ind) = 1
			compute_q(ind) = 1
			step_one(ind) = 1
			step_two(ind) = 1
			step_three(ind) = 1
		Enddo
		Do i = 1, nq_shellav
			ind = qvals_shellav(i)
			compute_shellav(ind) = 1
			compute_q(ind) = 1
			step_one(ind) = 1
			step_two(ind) = 1
		Enddo	
		Do i = 1, nq_azav
			ind = qvals_azav(i)
			compute_q(ind) = 1
			compute_azav(ind) = 1
			step_one(ind) = 1
		Enddo
		Do i = 1, nq_shell
			ind = qvals_shell(i)
			compute_q(ind) = 1
			compute_shell(ind) =1
			step_zero(ind) = 1
		Enddo

		Do i = 1, nq_3d
			ind = qvals_3d(i)
			step_five(ind) = 1
		Enddo

		my_nshells = 0
		Allocate(my_shell_levs(1:nshell_levels))
		Allocate(have_shell(1:nshell_levels))
		Allocate(my_shell_ind(1:nshell_levels))		
		Do j = 1, nshell_levels
			ilocal = shell_levels(j)
			If ((ilocal .ge. my_rmin) .and. (ilocal .le. my_rmax)) Then ! my processor has this radius
			   have_shell(j) = 1
			   my_nshells = my_nshells+1
			   my_shell_levs(my_nshells) = shell_levels(j)
			   my_shell_ind(my_nshells) = j
			Endif
		Enddo
		If (myid .eq. 0) Then
			!Use process zero to figure out which radial processors will
			! actually output shells.  This will make it easier to read in the input later.
			Allocate(ptemp(1:nproc1))
			Allocate(ptemp2(1:nproc1))
			Allocate(ptemp3(1:nproc1))
			ptemp(:) = 0
			ptemp2(:) = 0
			do i = 0, nproc1-1
				your_rmin = pfi%all_1p(i)%min
				your_rmax = pfi%all_1p(i)%max
				Do j = 1, nshell_levels
					ilocal = shell_levels(j)
					If ( (ilocal .ge. your_rmin) .and. (ilocal .le. your_rmax) ) Then
						ptemp(i+1) = ptemp(i+1)+1	!  radial id "i" has another shell that we want to output
					Endif	
				Enddo
			enddo
			pcount =0
			do i = 0, nproc1-1
				if (ptemp(i+1) .ge. 1) then 
					pcount = pcount+1					! pcount is the number of unique radial id's that have shells
					ptemp2(pcount) = i				! ptemp2 is for storing the radial id's of that have shells to output
					ptemp3(pcount) = ptemp(i+1)	! ptemp3 is the number of shells this radial id has
				endif
			enddo
			!   Resize the temporary arrays
			Allocate(shell_r_ids(1:pcount))
			Allocate(nshells_at_rid(1:pcount))
			shell_r_ids(:) = ptemp2(1:pcount)	! These are the radial ids of the processors that have shells
			nshells_at_rid(:) = ptemp3(1:pcount) ! How many shells this rid has
			nshell_r_ids = pcount
			DeAllocate(ptemp3)
			DeAllocate(ptemp2)
			DeAllocate(ptemp)
		Endif

	   !int_dphi = 0.0D0	
	   !do k = ks, ke
		!	int_dphi = int_dphi+dx3a(k)
		!enddo

		!Allocate(sintheta_dtheta(js:je))
		!sintheta_dtheta(:) = 0.0D0
		!int_sintheta_dtheta = 0.0D0
		!do j = js, je
		!  sintheta_dtheta(j) = g32b(j)*dx2a(j)
		!  int_sintheta_dtheta = int_sintheta_dtheta+sintheta_dtheta(j)
		!enddo

		!Allocate(rsquared_dr(is:ie))
		!rsquared_dr(:) = 0.0D0
		!int_rsquared_dr = 0.0D0
		!do i = is, ie
		!	rsquared_dr(i) = g2b(i)*g31b(i)*dx1a(i)					
		!	int_rsquared_dr = int_rsquared_dr+rsquared_dr(i)
		!enddo
	
   	!int_vol = int_rsquared_dr
		!int_vol = int_vol*int_sintheta_dtheta
		!int_vol = int_vol*int_dphi
	

   End Subroutine Initialize_Spherical_IO


	Subroutine Get_Shell_Slice(qval,qty)
		Implicit None
		Integer :: j, ilocal, shell_lev_ind
		Integer, Intent(In) :: qval
		Real*8, Intent(In) :: qty(:,:,my_theta_min:)

		If (my_nshells .gt. 0) Then


		  !If (qval .eq. first_shell_q) Then
		  If (start_shell_slice) Then
			Allocate(shell_slice_outputs(1:nphi,my_theta_min:my_theta_max,1:my_nshells,1:nq_shell))
			shell_ind = 1	
			start_shell_slice = .false.
		  Endif

		  shell_lev_ind =1
		  Do j = 1, nshell_levels
		    ilocal = shell_levels(j)-my_rmin+1
		    If (have_shell(j) .eq. 1) Then ! my processor has this radius
				!Write(6,*)'my shells: ', have_shell(j), shell_levels(j), my_rmin, my_rmax
				shell_slice_outputs(:,my_theta_min:my_theta_max,shell_lev_ind,shell_ind) = qty(:,ilocal,my_theta_min:my_theta_max)
				shell_lev_ind = shell_lev_ind +1
		    Endif
		  Enddo

		  shell_ind = shell_ind+1	! advance counter for next quantity to store (if any)
		Endif

	End Subroutine Get_Shell_Slice

	Subroutine Write_Shell_Slices(this_iter)
		Implicit None
		Real*8, Allocatable :: buff(:,:,:,:), all_shell_slices(:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,ii,qq, this_id2, this_id1, i_start, i_end, j_start, j_end, k_total
		Integer :: this_iter, n, nn, this_id,this_nshell
		Character*8 :: iterstring
		Character*120 :: shell_slice_file

		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: nelem

		! Later, we may want to generalize this, but for now just assume the process 0 handles the output
		responsible = 0
		If (myid .eq. io_node) responsible = 1


		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the shell slices from the other nodes

		
			Allocate(all_shell_slices(1:nphi,1:ntheta,nshell_levels,nq_shell))
			all_shell_slices(:,:,:,:) = 0.0d0
			current_shell = 1
			Do n = 1, nshell_r_ids  
				 this_rid = shell_r_ids(n)
				 this_nshell = nshells_at_rid(n)		
				 s_start = current_shell
				 s_end = s_start+this_nshell-1 


				 Do nn = 0, nproc2-1
					your_id = nn+this_rid*nproc2

					your_ntheta    = pfi%all_2p(nn)%delta
					your_theta_min = pfi%all_2p(nn)%min
					your_theta_max = pfi%all_2p(nn)%max

				 	Allocate(buff(1:nphi,1:your_ntheta,1:this_nshell,1:nq_shell))


					nelem = nphi*your_ntheta*this_nshell*nq_shell
				 	If (your_id .ne. io_node) then
						! Receive and copy

             		Call receive(buff, source= your_id,tag=shell_slice_tag,grp = pfi%gcomm)

						all_shell_slices(1:nphi,your_theta_min:your_theta_max,s_start:s_end,1:nq_shell) = &
							& buff(1:nphi , 1:your_ntheta, 1:this_nshell, 1:nq_shell)
					Else
						If (my_nshells .gt. 0) Then
							! Copy my shells into the main output array

							all_shell_slices(1:nphi,your_theta_min:your_theta_max,s_start:s_end,1:nq_shell) &
								& = shell_slice_outputs(:,:,:,:)
						Endif			
					Endif
					DeAllocate(buff)
 				Enddo
				current_shell = s_end+1
				!DeAllocate(buff)	! lot of allocation/de-allocation here.  Probably doesn't matter, but can optimize later.
			Enddo

      	Write(iterstring,i_ofmt) this_iter
         !shell_slice_file = trim(run_dir)//'/Shell_Slices/'//trim(iterstring)
			shell_slice_file = 'Shell_Slices/'//trim(iterstring)
	 		Open(unit=15,file=shell_slice_file,form='unformatted', status='replace')
         Write(15)nphi,ntheta,nshell_levels,nq_shell
         Write(15)(qvals_shell(i),i=1,nq_shell)
	 		Write(15)(radius(shell_levels(i)),i=1,nshell_levels)
	 		Write(15)(sintheta(j),j=1,ntheta)
			Write(15)((((all_shell_slices(k,j,i,qq),k=1,nphi),j=1,ntheta),i=1,nshell_levels),qq=1,nq_shell)
			Write(15)(costheta(j),j=1,ntheta)
			Close(15)

			DeAllocate(all_shell_slices)

		Else
			!  Non responsible nodes send their info
			If (my_nshells .gt. 0) Then
				Call send(shell_slice_outputs, dest = 0,tag=shell_slice_tag, grp = pfi%gcomm)
			Endif
		Endif
		If (my_nshells .gt. 0) Then
			DeAllocate(shell_slice_outputs)
		Endif

	End Subroutine Write_Shell_Slices


	Subroutine Add_Quantity( qval,qty)
		Implicit None
		Logical :: log_exist
		Integer :: i, this_iter, qval, oerr
		character*120 :: ologfile
		Real*8, Intent(In) :: qty(:,:,:)

		! write(6,*)'main_output'
		!ologfile = 'Raw/output.log'
		!shell_data_written = 0

		!If (myid .eq. 0) Then
		!   inquire(file='Raw/output.log',exist=log_exist)
		!   If (log_exist .eq. .true.) Then
		!      open(unit=13,file=ologfile,status='old', form='formatted', position='append')
		!   Else
		!      open(unit=13,file=ologfile,status='new', form='formatted')
		!   Endif
		!   write(13,i_ofmt)this_iter
		!   close(13)	
		!Endif

		!Call Output_Allocation()


		If (compute_q(qval) .eq. 1) Then
			If (step_zero(qval) .eq. 1) Then
				Call get_shell_slice(qval,qty)
			Endif
			
			If (step_one(qval) .eq. 1) Then 
				!!write(6,*)'in step 1'
			    Call get_azimuthal_average(this_iter, qty,qval)
				!output_ireq(1) = az_avg_tag
				!call MPI_BARRIER(MPI_COMM_WORLD,oerr)
					!write(6,*)'out step 1'
			Endif
			If (step_two(qval) .eq. 1) Then
			    Call get_shell_average(this_iter,qval)
				!output_ireq(1) = shell_avg_tag
				!call MPI_BARRIER(MPI_COMM_WORLD,oerr)
			Endif
			If (step_three(qval) .eq. 1) Then
			    Call get_global_average(this_iter,qval)
				!output_ireq(1) = global_avg_tag
				!call MPI_BARRIER(MPI_COMM_WORLD,oerr)
			Endif

			! Step 5 is separate from the others.  
			! For 3-D output, each output quantity is written as it is computed
			If (step_five(qval) .eq. 1) Then
				Write(6,*)'Entering 3D output'
				Call write_full_3d(qty,qval)
			Endif
		Endif



	!	output_ireq(1) = shell_slice_tag
	!	call MPI_WAITALL(1,output_ireq,output_status,oerr)
	
		!Call MPI_BARRIER(MPI_COMM_WORLD,oerr)


	End Subroutine	Add_Quantity

	Subroutine Complete_Output(iter)
		If (maxval(step_zero) .eq. 1) Then
			 Call Write_Shell_Slices(iter)
		Endif
		If (maxval(step_one) .eq. 1) Then
			 Call Write_Azimuthal_Average(iter)
		Endif
		If (maxval(step_two) .eq. 1) Then
			 Call Write_Shell_Average(iter)
		Endif
		If (maxval(step_three) .eq. 1) Then
			 Call Write_Global_Average(iter)
		Endif
        DeAllocate(f_of_r_theta)
        DeAllocate(f_of_r)
	End Subroutine Complete_Output

	Subroutine Get_Azimuthal_Average(this_iter,qty,qval)
		Implicit None
		Integer :: this_iter
		Integer :: r, t
		Integer, Intent(in) :: qval
		Real*8, Intent(In) :: qty(:,my_rmin:,my_theta_min:)
		If (start_az_average) Then
			Allocate(azav_outputs(my_rmin:my_rmax,my_theta_min:my_theta_max,1:nq_azav))
			azav_ind = 1			
			start_az_average = .false.
		Endif
		
		f_of_r_theta(:,:) = 0.0D0
		
		Do t = my_theta_min, my_theta_max
			Do r = my_rmin, my_rmax
				f_of_r_theta(r,t) = sum(qty(:,r,t))
			Enddo
		Enddo

		!Write(6,*)'my max: ', maxval(f_of_r_theta)
		f_of_r_theta = f_of_r_theta/dble(nphi)     ! average in phi

		If (compute_azav(qval) .eq. 1) Then
			azav_outputs(:,:,azav_ind) = f_of_r_theta
			azav_ind = azav_ind+1
		Endif
		
	END Subroutine Get_Azimuthal_Average

	Subroutine Get_Shell_Average(this_iter, qval)
		Implicit None
		Integer, Intent(In) :: this_iter, qval
        Integer :: t


		If (start_shell_average) Then
			shellav_ind = 1			
			Allocate(shellav_outputs(my_rmin:my_rmax,1:nq_shellav))			
            start_shell_average = .false.
		Endif
		f_of_r(:) = 0.0D0

        Do t = my_theta_min, my_theta_max
            f_of_r(:) = f_of_r(:) + f_of_r_theta(:,t)*theta_integration_weights(t)
        Enddo


		If (compute_shellav(qval) .eq. 1) Then
			shellav_outputs(:,shellav_ind) = f_of_r(:)
			shellav_ind = shellav_ind+1
		Endif

	END Subroutine Get_Shell_Average		

	Subroutine Get_Global_Average(this_iter, qval)
		Implicit None
		Integer, Intent(In) :: this_iter, qval
        Integer :: i
        Real*8 :: this_average

		If (start_global_average) Then
			globav_ind = 1			
			Allocate(globav_outputs(1:nq_globav))			
            start_global_average = .false.
		Endif
        this_average =0.0d0
        do i = my_rmin, my_rmax
            this_average = this_average+f_of_r(i)*r_integration_weights(i)
        enddo

		If (compute_globav(qval) .eq. 1) Then
			globav_outputs(globav_ind) = this_average
			globav_ind = globav_ind+1
		Endif

	END Subroutine Get_Global_Average

	Subroutine Write_Azimuthal_Average(this_iter)
		Implicit None
		Real*8, Allocatable :: buff(:,:,:), full_azavg(:,:,:)

		Integer :: responsible
		Integer :: i, j, k
		Integer, Intent(In) :: this_iter
        Integer :: n, nn
		Character*8 :: iterstring
		Character*120 :: azfile	
        Integer :: your_r_min, your_r_max, your_theta_min
        Integer :: your_id, your_theta_max, your_nr, your_ntheta

		responsible = 0
		If (myid .eq. io_node) responsible = 1

		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the az_averages from the other nodes


            Allocate(full_azavg(nr,ntheta,nq_azav))
 
            Do n = 0, nproc1-1  
                your_r_min = pfi%all_1p(n)%min
                your_r_max = pfi%all_1p(n)%max
                your_nr = pfi%all_1p(n)%delta                

					Do nn = 0, nproc2-1
                    your_id = nn+n*nproc2
                    your_ntheta    = pfi%all_2p(nn)%delta
                    your_theta_min = pfi%all_2p(nn)%min
                    your_theta_max = pfi%all_2p(nn)%max

                    Allocate(buff(your_r_min:your_r_max,your_theta_min:your_theta_max,1:nq_azav))
				    If (your_id .ne. io_node) then		
                        Call receive(buff, source= your_id,tag=az_avg_tag,grp = pfi%gcomm)
				    Else
					    buff(:,:,:) = azav_outputs(:,:,:)
                        DeAllocate(azav_outputs)
                    Endif
                    full_azavg(your_r_min:your_r_max,your_theta_min:your_theta_max,1:nq_azav) &
                        & = buff(:,:,:)
                    DeAllocate(buff)
                Enddo
            Enddo			
            write(iterstring,i_ofmt) this_iter
            azfile = 'AZ_Avgs/'//trim(iterstring)
            Open(unit=15,file=azfile,form='unformatted', status='replace')
            Write(15)nr, ntheta,nq_azav
            Write(15)(qvals_azav(i),i=1,nq_azav)
            Write(15)(radius(i),i=1,nr)
            Write(15)(sintheta(i),i=1,ntheta)
            Write(15)(((full_azavg(i,j,k),i=1,nr),j=1,ntheta),k=1,nq_azav)
            Close(15)
		Else
            Call send(azav_outputs, dest = 0,tag=az_avg_tag, grp=pfi%gcomm)
            DeAllocate(azav_outputs)
		Endif
    End Subroutine Write_Azimuthal_Average

	Subroutine Write_Global_Average(this_iter)
		Implicit None
		Real*8, Allocatable :: buff(:), full_avg(:)

		Integer :: responsible
		Integer :: i,n, your_id
		Integer, Intent(In) :: this_iter
		Character*8 :: iterstring
		Character*120 :: gfile	

		responsible = 0
		If (myid .eq. io_node) responsible = 1

		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the az_averages from the other nodes


            Allocate(full_avg(nq_globav))
            Allocate(buff(nq_globav))
            Do n = 0, nproc-1
				If (n .ne. io_node) then		
                    Call receive(buff, source= n,tag=az_avg_tag,grp = pfi%gcomm)
				Else
				    buff(:) = globav_outputs(:)
                    DeAllocate(globav_outputs)
                Endif
                full_avg(:) = full_avg(:)+buff(:)
            Enddo			
            DeAllocate(buff)

            write(iterstring,i_ofmt) this_iter
            gfile = 'G_Avgs/'//trim(iterstring)
            Open(unit=15,file=gfile,form='unformatted', status='replace')
            Write(15)nq_globav
            Write(15)(qvals_globav(i),i=1,nq_globav)
            Write(15)(full_avg(i),i=1,nq_globav)
            Close(15)
		Else
            Call send(globav_outputs, dest = 0,tag=az_avg_tag, grp=pfi%gcomm)
            DeAllocate(globav_outputs)
		Endif
    End Subroutine Write_Global_Average


	Subroutine Write_Shell_Average(this_iter)
		Implicit None
        Integer, Intent(In) :: this_iter
		Integer :: responsible
		Integer :: i, j, k, n, nn
		Character*8 :: iterstring
		Character*120 :: shellav_file
        Real*8, Allocatable :: full_shellavg(:,:), buff(:,:)		
        Integer :: your_r_min, your_r_max, your_nr, your_id

		responsible = 0
		If (myid .eq. io_node) responsible = 1

		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the shell_averages from the other nodes

			Allocate(full_shellavg(1:nr,nq_shellav))
			full_shellavg(:,:) = 0.0D0

			Do n = 0, nproc1-1  

                your_r_min = pfi%all_1p(n)%min
                your_r_max = pfi%all_1p(n)%max
                your_nr = pfi%all_1p(n)%delta                

                Do nn = 0, nproc2-1
                    your_id = nn+n*nproc2

	     		 
                    Allocate(buff(your_r_min:your_r_max,nq_shellav))
                    If (your_id .ne. io_node) then
                        Call receive(buff, source= your_id,tag=shell_avg_tag,grp=pfi%gcomm)
                    Else
                        buff = shellav_outputs
                        DeAllocate(shellav_outputs)
                    Endif
                    full_shellavg(your_r_min:your_r_max,:) = &
                        & full_shellavg(your_r_min:your_r_max,:) + &
                        & buff(your_r_min:your_r_max,:)
                    DeAllocate(buff)
                Enddo
			Enddo

            write(iterstring,i_ofmt) this_iter
            shellav_file = 'Shell_Avgs/'//trim(iterstring)
            open(unit=15,file=shellav_file,form='unformatted', status='replace')
            Write(15)nr, nq_shellav
            Write(15)(qvals_shellav(i),i=1,nq_shellav)
            Write(15)(radius(i),i=1,nr)
            Write(15)((full_shellavg(i,k),i=1,nr),k=1,nq_shellav)
            Close(15)

		Else
            !  Non responsible nodes send their info				 
            Call send(shellav_outputs, dest = 0,tag=shell_avg_tag, grp = pfi%gcomm)
            DeAllocate(shellav_outputs)
		Endif
	End Subroutine Write_Shell_Average

	Subroutine Write_Full_3D(qty,qtag)
		Use MPI_BASE ! Doing this here for now.  No other routine above sees MPI_Base, and I may want to keep it that way.
		Implicit None		
		Real*8, Intent(In) :: qty(:,my_rmin:,my_theta_min:)
		Integer, Intent(In) :: qtag
		Real*8, Allocatable :: my_shells(:,:,:), buff(:,:,:)
		Integer :: i, j, k
		Character*2 :: qstring
		Character*8 :: iterstring
		Character*120 :: cfile

		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: np, jind, buffsize, p
		Integer(kind=MPI_OFFSET_KIND) :: my_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
		Integer :: funit, ierr

		! qty is dimensioned 1:n_phi, my_rmin:my_rmax, my_theta_min:my_theta_max


		If (my_row_rank .eq. 0) Then
			! Everyone in the row communicates to row-rank zero.
			! Each row owns a specific rank of radii
			Allocate(my_shells(1:nphi, 1:ntheta, my_rmin:my_rmax))
			my_shells(:,:,:) = 0.0d0					

			! First each rank stripes its own data into the new array
			! Not that striping is not in the "native" order
			Do j = my_theta_min, my_theta_max
				Do i = my_rmin, my_rmax
					my_shells(:,j,i) = qty(:,i,j)
				Enddo
			Enddo
			np = pfi%rcomm%np
			Do p = 1, np-1	
				your_theta_min = pfi%all_2p(p)%min
				your_theta_max = pfi%all_2p(p)%max
				your_ntheta    = pfi%all_2p(p)%delta
				Allocate(buff(1:nphi,my_rmin:my_rmax, your_theta_min:your_theta_max))
				Call receive(buff, source= p,tag=full_3d_tag,grp = pfi%rcomm)
				Do j = your_theta_min, your_theta_max
					Do i = my_rmin, my_rmax
						my_shells(:,j,i) = buff(:,i,j)
					Enddo
				Enddo
				DeAllocate(buff)
			Enddo
			! Now do the MPI write
			np = pfi%ccomm%np
			my_disp = 0
			Do p = 1, my_column_rank
				my_disp = my_disp+pfi%all_1p(p-1)%delta
			Enddo
			my_disp = my_disp*ntheta*nphi*8	! Displacment of this rank in the MPI File
			buffsize = my_nr*nphi*ntheta		! Number of elements to write


			write(iterstring,i_ofmt) current_iteration
			write(qstring,'(i2.2)') qtag
         cfile = 'Spherical_3D/'//trim(iterstring)//'_'//qstring
			Write(6,*)cfile
			call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                   MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                   MPI_INFO_NULL, funit, ierr) 

			call MPI_FILE_SET_VIEW(funit, my_disp, MPI_DOUBLE_PRECISION, & 
                   MPI_DOUBLE_PRECISION, 'native', & 
                   MPI_INFO_NULL, ierr) 
			call MPI_FILE_WRITE(funit, my_shells, buffsize, MPI_DOUBLE_PRECISION, & 
                   mstatus, ierr) 

			call MPI_FILE_CLOSE(funit, ierr) 
			If (my_column_rank .eq. 0) Then
				! row/column 0 writes out a file with the grid, etc.
				! This file should contain everything that needs to be known for processing later
	         write(iterstring,'(i8.8)') current_iteration
            cfile = 'Spherical_3D/'//trim(iterstring)//'_'//'grid'
	         open(unit=15,file=cfile,form='unformatted', status='replace')
	         Write(15)nr
				Write(15)ntheta
				Write(15)nphi
	         Write(15)(radius(i),i=1,nr)
				Write(15)(acos(costheta(i)),i = 1, ntheta)
	         Close(15)
			Endif


		Else
			! Send an array that's indexed starting at 1.  Shouldn't be necessary, but just in case.
			!Allocate(buff(1:nphi,1:my_ntheta,1:my_nr))
			
			!Do j = my_theta_min, my_theta_max
			!	jind = j-my_theta_min+1
			!	Do i = my_rmin, my_rmax
			!		buff(:,jind,i) = qty(:,i,j)
			!	Enddo
			!Enddo
			Call send(qty, dest= 0,tag=full_3d_tag,grp = pfi%rcomm)
			!DeAllocate(buff)

		Endif	


	End Subroutine Write_Full_3D

End Module Spherical_IO
