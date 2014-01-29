Module Checkpointing
	Use ProblemSize
	Use Parallel_Framework
	Use Linear_Solve, Only : get_all_rhs
	! Simple Checkpointing Module
	! Uses MPI-IO to split writing of files amongst rank zero processes from each row
	Implicit None
	Type(SphericalBuffer) :: chktmp
	Integer, private :: numfields = 4 ! 6 for hydro
	Integer,private,Allocatable :: mode_count(:)
	Integer,private :: nlm_total
	Integer, Allocatable, Private :: lmstart(:)
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
		Endif
	End Subroutine Initialize_Checkpointing

	Subroutine Write_Spherical_Array(abterms)
		Implicit None
		Real*8, Intent(In) :: abterms(:,:,:,:)
		Integer :: mp, m, nm,nmodes, offset,nl,p,np
		Integer :: dim2, lstart
		Real*8, Allocatable :: myarr(:,:), rowstrip(:,:)

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
		Endif

		
	End Subroutine Write_Spherical_Array

End Module Checkpointing
