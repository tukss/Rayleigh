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
	Integer, Private :: my_check_disp, buffsize
	Character*3 :: wchar = 'W', pchar = 'P', tchar = 'T', zchar = 'Z'
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
			my_check_disp = my_check_disp*nlm_total*2
			buffsize = nlm_total*2*my_r%delta

		Endif
	End Subroutine Initialize_Checkpointing



	Subroutine Write_Checkpoint(abterms,iteration)
		Implicit None
		Real*8, Intent(In) :: abterms(:,:,:,:)
		Integer, Intent(In) :: iteration
		Integer :: mp, m, nm,nmodes, offset,nl,p,np
		Integer :: dim2, lstart, grid_type, i
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
					grid_type = 0	! uniform grid
					if (chebyshev) grid_type = 1
	         	write(iterstring,'(i8.8)') iteration
            	cfile = 'Checkpoints/'//trim(iterstring)//'_'//'grid_etc'
	            open(unit=15,file=cfile,form='unformatted', status='replace')
	            Write(15)n_r
					Write(15)grid_type
					Write(15)l_max
	            Write(15)(radius(i),i=1,N_R)
	            Close(15)

				Endif
		Endif

		
	End Subroutine Write_Checkpoint

	Subroutine Read_Checkpoint(abterms,iteration)
		Implicit None
		Integer :: n_r_old, l_max_old, grid_type_old
		Integer :: old_pars(3)
		Real*8, Allocatable :: old_radius(:)
		! Rank zero from each row participates in the read
		If (my_row_rank .eq. 0) Then
			If (my_column_rank .eq. 0) 
				!process zero reads all the old info and broadcasts to the other row rank zeros
            	cfile = 'Checkpoints/'//trim(iterstring)//'_'//'grid_etc'
	            open(unit=15,file=cfile,form='unformatted', status='old')
	            Read(15)n_r_old
					Read(15)grid_type_old
					Read(15)l_max_old
	            Read(15)(radius(i),i=1,N_R)
	            Close(15)				
			Endif
		Endif

	End Subroutine Read_Checkpoint

	Subroutine Write_Field(arr,ind,tag,iter)
				Implicit None
				Integer, Intent(In) :: ind, iter
				Real*8, Intent(In) :: arr(1:,1:)
				Character*8 :: iterstring
				Character*3, Intent(In) :: tag
				Character*120 :: cfile

				integer ierr, i, funit , v_offset
				integer(kind=MPI_OFFSET_KIND) disp 
				Integer :: mstatus(MPI_STATUS_SIZE)
	         write(iterstring,'(i8.8)') iter
            cfile = 'Checkpoints/'//trim(iterstring)//'_'//trim(tag)


  
 				v_offset = (ind-1)*tnr+1

				call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                       MPI_INFO_NULL, funit, ierr) 
				disp = my_check_disp*8
				call MPI_FILE_SET_VIEW(funit, disp, MPI_DOUBLE_PRECISION, & 
                           MPI_DOUBLE_PRECISION, 'native', & 
                           MPI_INFO_NULL, ierr) 
				call MPI_FILE_WRITE(funit, arr(1,v_offset), buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
				call MPI_FILE_CLOSE(funit, ierr) 
     

                                      
                                          
	End Subroutine Write_Field

End Module Checkpointing
