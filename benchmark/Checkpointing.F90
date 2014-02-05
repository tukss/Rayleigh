Module Checkpointing
	Use ProblemSize
	Use Parallel_Framework
	Use Linear_Solve, Only : get_all_rhs
	Use SendReceive
	Use MPI_BASE
	! Simple Checkpointing Module
	! Uses MPI-IO to split writing of files amongst rank zero processes from each row
	Implicit None
	Type(SphericalBuffer) :: chktmp
	Integer, private :: numfields = 4 ! 6 for hydro
	Integer,private,Allocatable :: mode_count(:)
	Integer,private :: nlm_total, checkpoint_tag = 425
	Integer, Allocatable, Private :: lmstart(:)
	Integer, Private :: my_check_disp, buffsize, full_disp		! for writing checkpoints
	Integer, Private :: my_in_disp, buffsize_in, full_in_disp		! for reading checkpoints
	Character*3 :: wchar = 'W', pchar = 'P', tchar = 'T', zchar = 'Z'
	Integer :: checkpoint_iter = 0
	Real*8 :: checkpoint_dt, checkpoint_newdt
Contains

	Subroutine Initialize_Checkpointing()
		Implicit None
		Integer :: nfs(6)
		Integer :: p, np, nl, m, mp
		nfs(:) = numfields*2
		Call chktmp%init(field_count = nfs, config = 'p1a')			! This structure hangs around through the entire run

		!////////////////////////
		np = pfi%rcomm%np
		Allocate(mode_count(0:np-1))
		mode_count(:)  = 0	! This is how many total l-m combinations rank p of a row owns

		
		
		nlm_total = 0			! This is the total number of l-m combinations
		Do p = 0, np -1
			Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
				m = m_values(mp)			
				nl = l_max-m+1
				mode_count(p) = mode_count(p)+nl
				nlm_total = nlm_total+nl
			Enddo
		Enddo
		if (my_row_rank .eq. 0) Then
			Allocate(lmstart(0:l_max))
			lmstart(0) = 1
			do m = 0, l_max-1
				nl = l_max-m+1
				lmstart(m+1) = lmstart(m)+nl
			enddo
			np = pfi%ccomm%np
			my_check_disp = 0
			Do p = 1, my_column_rank
				my_check_disp = my_check_disp+pfi%all_1p(p-1)%delta
			enddo
			my_check_disp = my_check_disp*nlm_total !*2
			buffsize = nlm_total*my_r%delta
			full_disp = N_r*nlm_total
		Endif
	End Subroutine Initialize_Checkpointing



	Subroutine Write_Checkpoint(abterms,iteration,dt,new_dt)
		Implicit None
		Real*8, Intent(In) :: abterms(:,:,:,:)
		Real*8, Intent(In) :: dt, new_dt
		Integer, Intent(In) :: iteration
		Integer :: mp, m, nm,nmodes, offset,nl,p,np
		Integer :: dim2, lstart, i
		Real*8, Allocatable :: myarr(:,:), rowstrip(:,:)
		Character*8 :: iterstring
		Character*120 :: cfile
		np = pfi%rcomm%np

		Call chktmp%construct('p1a')
		chktmp%config = 'p1a'
		!Copy the RHS into chtkmp
		Call Get_All_RHS(chktmp%p1a)
		chktmp%p1a(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
		!Now we want to move from p1a to s2a (rlm space)
		Call chktmp%reform()

		! Next, each process stripes their s2a array into a true 2-D array
		dim2 = tnr*numfields*2
		Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
		offset =1 
		Do mp = my_mp%min, my_mp%max
				m = m_values(mp)
				nl = l_max-m+1
				myarr(offset:offset+nl-1,:) = chktmp%s2a(mp)%data(m:l_max,:)
				offset = offset+nl
		Enddo
		Call chktmp%deconstruct('s2a')

		! Everyone sends to then 0 process of each row, who organizes the data into one large strip for output.		


		If (my_row_rank .ne. 0) Then
			! Send myarr

            Call Send(myarr, dest = 0, tag = checkpoint_tag, grp = pfi%rcomm)
			DeAllocate(myarr)
		Else
			Allocate( rowstrip(1:nlm_total, 1:tnr*numfields*2))
			! first, copy myarr into the larger array
			offset = 1
			Do mp = my_mp%min, my_mp%max
				m = m_values(mp)
				nl = l_max-m+1
				lstart = lmstart(m)
				rowstrip(lstart:lstart+nl-1,:) = myarr(offset:offset+nl-1,:)				
				offset = offset+nl
			Enddo
			DeAllocate(myarr)
			Do p = 1, np -1
					Allocate(myarr(1:mode_count(p),1:dim2))
					! Receive
                    Call receive(myarr, source= p,tag=checkpoint_tag,grp = pfi%rcomm)
					offset = 1
					Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
						m = m_values(mp)
						nl = l_max-m+1
						lstart = lmstart(m)
						rowstrip(lstart:lstart+nl-1,:) = myarr(offset:offset+nl-1,:)				
						offset = offset+nl
					Enddo										
					DeAllocate(myarr)
			Enddo
				Call Write_Field(rowstrip,1,wchar, iteration)
				Call Write_Field(rowstrip,2,pchar, iteration)
				Call Write_Field(rowstrip,3,tchar, iteration)
				Call Write_Field(rowstrip,4,zchar, iteration)

				Call Write_Field(rowstrip,5,'WAB', iteration)
				Call Write_Field(rowstrip,6,'PAB', iteration)
				Call Write_Field(rowstrip,7,'TAB', iteration)
				Call Write_Field(rowstrip,8,'ZAB', iteration)

            DeAllocate(rowstrip)
				If (my_column_rank .eq. 0) Then
					! row/column 0 writes out a file with the grid, etc.
					! This file should contain everything that needs to be known
	         	write(iterstring,'(i8.8)') iteration
            	cfile = 'Checkpoints/'//trim(iterstring)//'_'//'grid_etc'
	            open(unit=15,file=cfile,form='unformatted', status='replace')
	            Write(15)n_r
					Write(15)grid_type
					Write(15)l_max
					Write(15)dt
					Write(15)new_dt
	            Write(15)(radius(i),i=1,N_R)
	            Close(15)

				Endif
		Endif

		
	End Subroutine Write_Checkpoint

	Subroutine Read_Checkpoint(fields, abterms,iteration)
		Implicit None
		Integer, Intent(In) :: iteration
		Real*8, Intent(InOut) :: fields(:,:,:,:), abterms(:,:,:,:)
		Integer :: n_r_old, l_max_old, grid_type_old, nr_read
		Integer :: i, ierr, nlm_total_old, m, nl,p, np, mxread
		Integer :: maxl, dim2,offset, nl_load,lstart,mp
		Integer :: old_pars(3)
		Integer, Allocatable :: lmstart_old(:)
		Real*8, Allocatable :: old_radius(:)
		Real*8, Allocatable :: rowstrip(:,:), myarr(:,:), sendarr(:,:)
		Real*8 :: dt_pars(2),dt,new_dt
		Character*8 :: iterstring
		Character*120 :: cfile
		dim2 = tnr*numfields*2
		checkpoint_iter = iteration
		Write(iterstring,'(i8.8)') iteration
		If (my_rank .eq. 0) Then
			!process zero reads all the old info and broadcasts to all other ranks
          cfile = 'Checkpoints/'//trim(iterstring)//'_'//'grid_etc'
	       open(unit=15,file=cfile,form='unformatted', status='old')
	       Read(15)n_r_old
			 Read(15)grid_type_old
			 Read(15)l_max_old
			 Read(15)dt
			 Read(15)new_dt
			 Allocate(old_radius(1:N_r_old))
	       Read(15)(old_radius(i),i=1,N_R)
	       Close(15)				
			 old_pars(1) = n_r_old
			 old_pars(2) = grid_type_old
			 old_pars(3) = l_max_old
			 dt_pars(1) = dt
			 dt_pars(2) = new_dt

			If (l_max_old .lt. l_max) Then
					Write(6,*)' '
					Write(6,*)'#####################################################################'
					Write(6,*)'# '
					Write(6,*)'#  Checkpoint horizontal resolution is lower than current resolution.'
					Write(6,*)'#  The old solution will be interpolated onto horizontal grid with '
					Write(6,*)'#  higher resolution corresponding to the new l_max.'
					Write(6,*)'#  Old l_max: ', l_max_old
					Write(6,*)'#  New l_max: ', l_max
					Write(6,*)'# '
					Write(6,*)'#####################################################################'
					Write(6,*)' '
			Endif
			If (l_max_old .gt. l_max) Then
					Write(6,*)' '
					Write(6,*)'#####################################################################'
					Write(6,*)'# '
					Write(6,*)'#  Checkpoint horizontal resolution is higher than current resolution.'
					Write(6,*)'#  The old SPH expansion will be truncated at the new l_max.'
					Write(6,*)'#  This might not be a good idea.'
					Write(6,*)'#  Old l_max: ', l_max_old
					Write(6,*)'#  New l_max: ', l_max
					Write(6,*)'# '
					Write(6,*)'#####################################################################'
					Write(6,*)' '
			Endif
		Endif
		Call MPI_Bcast(old_pars,3, MPI_INTEGER, 0, pfi%gcomm%comm, ierr)

		n_r_old = old_pars(1)
		grid_type_old = old_pars(2)
		l_max_old = old_pars(3)
		!///////// Later we only want to do this if the grid is actually different
		If (my_rank .ne. 0) Then
			Allocate(old_radius(1:n_r_old))
		Endif
		Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%gcomm%comm, ierr)
		Call MPI_Bcast(dt_pars,2, MPI_DOUBLE_PRECISION, 0, pfi%gcomm%comm, ierr)
		checkpoint_dt = dt_pars(1)
		checkpoint_newdt = dt_pars(2)

		Call chktmp%construct('s2b')
		chktmp%config = 's2b'

		! Rank zero from each row participates in the read
		If (my_row_rank .eq. 0) Then

			n_r_old = old_pars(1)
			grid_type_old = old_pars(2)
			l_max_old = old_pars(3)

			! Column zero reads in the old checkpoint no matter what
			! the old dimensions are. We'll broadcast back to all members of the row later

			! Interpolation/truncation of the old spherical harmonic basis is easy
			! Interpolation up in radius is easy, but down is more difficult to program
			! and unlikely to be used.  I will write a serial version of the input
			! for those cases where interpolation down is required.

			nlm_total_old = 0			! This is the total number of l-m combinations in the CHECKPOINT FILE
			Do m = 0, l_max_old
				nl = l_max_old-m+1
				nlm_total_old = nlm_total_old+nl
			Enddo

			np = pfi%ccomm%np
			my_in_disp = 0
			Do p = 1, my_column_rank
				my_in_disp = my_in_disp+pfi%all_1p(p-1)%delta
			enddo
			my_in_disp = my_in_disp*nlm_total_old !*2
			buffsize_in = nlm_total_old*my_r%delta
			full_in_disp = nlm_total_old*n_r_old

			mxread = n_r_old-my_r%min+1
			nr_read = min(tnr,mxread)
			if(my_r%min .gt. n_r_old) nr_read = 0

			Allocate( rowstrip(1:nlm_total_old, 1:tnr*numfields*2))	
			rowstrip(:,:) = 0

			Call Read_Field(rowstrip,1,wchar, iteration)
			Call Read_Field(rowstrip,2,pchar, iteration)
			Call Read_Field(rowstrip,3,tchar, iteration)
			Call Read_Field(rowstrip,4,zchar, iteration)

			Call Read_Field(rowstrip,5,'WAB', iteration)
			Call Read_Field(rowstrip,6,'PAB', iteration)
			Call Read_Field(rowstrip,7,'TAB', iteration)
			Call Read_Field(rowstrip,8,'ZAB', iteration)			

			! Now the head of each row owns all modes of each field at the
			! radii owned by that row.  The different modes now need to be 
			! distributed to their respect owners.
			! This is where horizontal interpolation or truncation is done (implicitly)
			! by only loading  appropriate modes into the send arrays

			!/////////////////////////////
			!  Do some book-keeping related to the old l_max
			!  This array tells us where in the first dimension of rowstrip, ell values for mode m start
			np = pfi%rcomm%np

			Allocate(lmstart_old(0:l_max_old))
			lmstart_old(0) = 1
			do m = 0, l_max_old-1
				nl = l_max_old-m+1
				lmstart_old(m+1) = lmstart_old(m)+nl
			enddo

			maxl = min(l_max,l_max_old)	! Take care to only read in modes that are common
													! to the checkpoint & and the current simulation

			! First, each row-head pulls out their own modes			
			Do mp = my_mp%min, my_mp%max
				m = m_values(mp)
				chktmp%s2b(mp)%data(:,:) = 0.0d0
				If (m .le. l_max_old) Then
					nl = maxl-m+1
					lstart = lmstart_old(m)
					chktmp%s2b(mp)%data(m:maxl,:) = rowstrip(lstart:lstart+nl-1,:)
				Endif
			Enddo


			!/////////////////////////////////////////////////////////////
			! Next send each of the other processors in the row their information
			Do p = 1, np -1
					Allocate(sendarr(1:mode_count(p),1:dim2))
					sendarr(:,:) = 0.0d0
					offset = 1
					Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
						m = m_values(mp)
						nl = l_max-m+1
						If (m .le. l_max_old) Then
							nl_load = maxl-m+1
							lstart = lmstart_old(m)
							sendarr(offset:offset+nl_load-1,:)	= rowstrip(lstart:lstart+nl_load-1,:) 			
						Endif
						offset = offset+nl
					Enddo										
               Call send(sendarr, dest= p,tag=checkpoint_tag,grp = pfi%rcomm)
					DeAllocate(sendarr)
			Enddo
			!/////////

			DeAllocate(lmstart_old)
			DeAllocate(rowstrip)

		Else			
			! Receive my modes
		
			Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
			Call receive(myarr, source= 0,tag=checkpoint_tag,grp = pfi%rcomm)
			offset =1 
			Do mp = my_mp%min, my_mp%max
				m = m_values(mp)
				nl = l_max-m+1
				chktmp%s2b(mp)%data(m:l_max,:) = myarr(offset:offset+nl-1,:)
				offset = offset+nl
			Enddo
			DeAllocate(myarr)

		Endif
		Call chktmp%reform()	! move to p1b

		! NOW, if n_r_old and grid_type_old are the same, we can copy chtkmp%p1b into abterms and 
		! fields.  Otherwise, we need to interpolate onto the current grid
		If  ((n_r_old .ne. n_r) .or. (grid_type_old .ne. grid_type) ) Then
			! Interpolate
			! We will assume the user kept the same radial domain bounds.
			! If they  have not, this will end badly.
			If (my_rank .eq. 0) Then
				Write(6,*)'Grid has changed.  Interpolating onto new grid.'
				Write(6,*)'Old grid_type:     ', grid_type_old
				Write(6,*)'Current grid_type: ', grid_type
				Write(6,*)'Old N_R:           ', n_r_old
				Write(6,*)'Current N_R:           ', n_r
			Endif
		Endif
		
		! Interpolation is complete, now we just copy into the other arrays
		fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)
		abterms(:,:,:,1:numfields) = chktmp%p1b(:,:,:,numfields+1:numfields*2)

		Call chktmp%deconstruct('p1b')
		DeAllocate(old_radius)

	End Subroutine Read_Checkpoint

	Subroutine Write_Field(arr,ind,tag,iter)
				Implicit None
				Integer, Intent(In) :: ind, iter
				Real*8, Intent(In) :: arr(1:,1:)
				Character*8 :: iterstring
				Character*3, Intent(In) :: tag
				Character*120 :: cfile

				integer ierr, i, funit , v_offset1, v_offset2
				integer(kind=MPI_OFFSET_KIND) disp1,disp2 
				Integer :: mstatus(MPI_STATUS_SIZE)
	         write(iterstring,'(i8.8)') iter
            cfile = 'Checkpoints/'//trim(iterstring)//'_'//trim(tag)


 				! We have to be careful here.  Each processor does TWO writes. 
				! The first write places the real part of the field into the file.
				! The view then changes and advances to the appropriate location of the
				! imaginary part.  This step is crucial for checkpoints to work with
				! Different processor configurations.
 				v_offset1 = (ind-1)*tnr+1
				v_offset2 = v_offset1+my_r%delta

				call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                       MPI_INFO_NULL, funit, ierr) 
				disp1 = my_check_disp*8
				disp2 = (my_check_disp+full_disp)*8

				call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, & 	! Real part
                           MPI_DOUBLE_PRECISION, 'native', & 
                           MPI_INFO_NULL, ierr) 
				call MPI_FILE_WRITE(funit, arr(1,v_offset1), buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

				call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, & 	! Imaginary part
                           MPI_DOUBLE_PRECISION, 'native', & 
                           MPI_INFO_NULL, ierr) 
				call MPI_FILE_WRITE(funit, arr(1,v_offset2), buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

				call MPI_FILE_CLOSE(funit, ierr) 
     			

                                      
                                          
	End Subroutine Write_Field


	Subroutine Read_Field(arr,ind,tag,iter)
				Implicit None
				Integer, Intent(In) :: ind, iter
				Real*8, Intent(In) :: arr(1:,1:)
				Character*8 :: iterstring
				Character*3, Intent(In) :: tag
				Character*120 :: cfile

				integer ierr, i, funit , v_offset1, v_offset2
				integer(kind=MPI_OFFSET_KIND) disp1,disp2 
				Integer :: mstatus(MPI_STATUS_SIZE)
	         write(iterstring,'(i8.8)') iter
            cfile = 'Checkpoints/'//trim(iterstring)//'_'//trim(tag)


  
 				v_offset1 = (ind-1)*tnr+1
				v_offset2 = v_offset1+my_r%delta

				call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                       MPI_MODE_RDONLY, & 
                       MPI_INFO_NULL, funit, ierr) 
				disp1 = my_in_disp*8
				disp2 = (my_in_disp+full_in_disp)*8
				call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, & 
                           MPI_DOUBLE_PRECISION, 'native', & 
                           MPI_INFO_NULL, ierr) 
				call MPI_FILE_READ(funit, arr(1,v_offset1), buffsize_in, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

				call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, & 
                           MPI_DOUBLE_PRECISION, 'native', & 
                           MPI_INFO_NULL, ierr) 
				call MPI_FILE_READ(funit, arr(1,v_offset2), buffsize_in, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

				call MPI_FILE_CLOSE(funit, ierr) 
     
				If (ind .eq. 3) then
					if (my_rank .eq. 0) then
						!Write(6,*)arr(1,v_offset:v_offset+tnr-1)
					endif
				Endif
                                      
                                          
	End Subroutine Read_Field

End Module Checkpointing
