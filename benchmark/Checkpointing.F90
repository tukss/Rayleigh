Module Checkpointing
	Use ProblemSize
	! Simple Checkpointing Module
	! Uses MPI-IO to split writing of files amongst rank zero processes from each row
	Implicit None
	Type(SphericalBuffer) :: chktmp
	Integer, private :: numfields = 4 ! 6 for hydro
	Integer,private,Allocatable :: m_count(:)
	Integer,private :: nlm_total
	Subroutine Initialize_Checkpointing()
		Implicit None
		Integer :: nfs(6)
		Integer :: p
		nfs(:) = numfields*2
		Call chktmp%init(field_count = nfs, config = 'p1a')			! This structure hangs around through the entire run

		!////////////////////////
		np = pfi%rcomm%np
		Allocate(m_count(0:np-1))
		m_count(:) = 0
		
		
		
		nlm_total = 0
		Do p = 0, np -1
			Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
				m = m_values(mp)			
				nm = l_max-m+1
				mode_count(p) = mode_count(p)+nm
				nlm_total = nlm_total+nm
			Enddo
		Enddo

	End Subroutine Initialize_Checkpointing

	Subroutine Write_Spherical_Array(abterms)
		Implicit None
		Type(SphericalBuffer), Intent(In) :: abterms
		Integer :: mp, m, nm,nmodes, offset
		Integer :: colrank, rowrank
		Real*8, Allocatable :: myarr(:,:,:)
		chktmp%construct('p1a')
		chktmp%config = 'p1a'
		!Copy the RHS into chtkmp
		Call Get_All_RHS(chktmp%p1a)
		chktmp(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
		!Now we want to move from p1a to s2a (rlm space)
		Call chktmp%reform()

		! Next, each process stripes their s2a array into a true 2-D array
		nmodes =  0
		Do mp = my_mp%min,my_mp%max
			m = m_values(mp)
			nm = l_max-m+1
			nmodes = nmodes+nm
		Enddo
		Allocate(myarr(1:nmodes,1:tnr*numfields*2))	! numfields is either 4 (hydro) or 6 (MHD)
		offset =1 
		Do mp = my_mp%min, my_mp%max
				m = m_values(mp)
				nm = l_max-m+1
				myarr(offset:offset+nm-1,:) = chktmp(mp)%s2a%data(m:l_max,:)
				offset = offset+nm
		Enddo

		! Everyone sends to then 0 process of each row, who organizes the data into one large strip for output.		
		colrank = pfi%romm%rank
		report = .false.
		rowrank = pfi%rcomm%rank
		if (rowrank .eq. 0) report = .true.


	End Subroutine Write_Spherical_Array

End Module Checkpointing
