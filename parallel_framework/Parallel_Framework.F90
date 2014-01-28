Module Parallel_Framework
	Use MPI_LAYER
	Use Load_Balance
	Use Structures
	Implicit None
	Private

	!///////////////////////////////////////////////////////////////////////////
	!  
	Integer, Parameter, Public :: Cartesian = 1, Cylindrical = 2, Spherical = 3
	Integer, Parameter, Public :: p1 =1 ,s1 = 2, p2 = 3,s2a =4, p3a=5,p3b=6
	Public :: Load_Config
	Character*6 :: ifmt = '(i4.4)' ! Integer format for indicating processor tag numbers in output
	Type, Public :: Parallel_Interface
		Type(Load_Config) :: my_1p, my_2p, my_3p	!  like my_r, my_theta in ASH 
		Type(Load_Config) :: my_1s, my_2s, my_3s	! 	like my_mp, my_n etc.
		Type(Load_Config), Allocatable :: all_1p(:), all_2p(:)
		Type(Load_Config), Allocatable :: all_3s(:)
		Integer :: Geometry
		! Global dimensions in 'a' and 'b' configurations
		! Allows primarily for dealiasing right now
		Integer :: n1p, n1s, n2p, n2s, n3p, n3s	
		Integer :: npe, nprow, npcol, npio,npc
		Integer, Allocatable :: inds_3s(:)
		Type(communicator) :: rcomm ! row communicator
		Type(communicator) :: ccomm ! column communicator
		Type(communicator) :: gcomm ! global communicator
		!##################################################
		Type(Load_Config), Allocatable, Public :: lb_1p(:), lb_1s(:)


		Contains

		Procedure :: Init => Initialize_Parallel_Interface
		Procedure :: Exit => Finalize_Framework
		Procedure ::  Spherical_Init 
		Procedure ::  Init_Geometry
	End Type Parallel_Interface

	Type, Public :: mcontainer
		Complex*16, Allocatable :: data(:,:,:)
	End Type mcontainer
	!Type, Public :: rmcontainer
	!	Real*8, Allocatable :: data(:,:)
	!End Type rmcontainer
	Type, Public :: eqcontainer
		Real*8, Allocatable :: data(:,:,:)
	End Type eqcontainer

	Type, Public :: SphericalBuffer
		! The buffer object for buffer moving between spaces
		! The buffer can advance configurations

		Type(rmcontainer), Allocatable :: s2b(:),s2a(:)

		Real*8, Allocatable :: p2b(:,:,:),p2a(:,:,:)		! for dgemm
		Real*8, Allocatable :: p3a(:,:,:,:)  ! m/r/delta_theta/field
		Real*8, Allocatable :: p3b(:,:,:,:)  ! something
		Real*8, Allocatable :: p1b(:,:,:,:),p1a(:,:,:,:)  ! something

		Integer :: nf1a = 1
		Integer :: nf2a = 1
		Integer :: nf3a = 1
		Integer :: nf3b = 1
		Integer :: nf2b = 1
		Integer :: nf1b = 1


		Integer, Allocatable :: sdisp12(:), rdisp12(:)
		Integer, Allocatable :: scount12(:), rcount12(:)

		Integer, Allocatable :: sdisp21(:), rdisp21(:)
		Integer, Allocatable :: scount21(:), rcount21(:)

		Integer, Allocatable :: sdisp23(:), rdisp23(:)
		Integer, Allocatable :: scount23(:), rcount23(:)

		Integer, Allocatable :: sdisp32(:), rdisp32(:)
		Integer, Allocatable :: scount32(:), rcount32(:)

		Integer, Allocatable :: sdisp32v2(:), rdisp32v2(:)
		Integer, Allocatable :: scount32v2(:), rcount32v2(:)

		Character*3 :: config
		!scount12(0) => number I send to rank 0 when going FROM 1 TO 2
		!scount21(0) => number I send to rank 0 when going FROM 2 TO 1
		Contains
		Procedure :: Init  => Initialize_Spherical_Buffer
		Procedure :: construct => Allocate_Spherical_Buffer
		Procedure :: deconstruct => DeAllocate_Spherical_Buffer
		Procedure :: reform => advance_configuration
		Procedure :: transpose_1a2a
		Procedure :: transpose_2a3a  ! Move from 2a to 3a
		Procedure :: transpose_3b2b
		Procedure :: transpose_2b1b
		Procedure :: write_space
	End Type SphericalBuffer


	Type(Parallel_Interface), Public :: pfi	



Contains
	Subroutine Initialize_Spherical_Buffer(self,report, field_count,config)
		! Buffer initialization
		! Handles send/receive disp/counts
		Implicit None
		Integer :: np, p
		Integer :: report_unit = 500
		Logical, Intent(In), Optional :: report
		Integer, Intent(In), Optional :: field_count(3,2)
		Character*120 :: report_file,report_tag
		Character*10 :: gtag, rtag, ctag
		Character*3, Intent(In), Optional :: config
		Class(SphericalBuffer) :: self
		If (present(report)) Then
			Write(gtag,ifmt)pfi%gcomm%rank
			Write(rtag,ifmt)pfi%rcomm%rank
			Write(ctag,ifmt)pfi%ccomm%rank
			report_tag = 'g'//Trim(gtag)//'_r'//Trim(rtag)//'_c'//Trim(ctag)
		Endif
		If (present(config)) Then
			self%config = config
		Else
			self%config = 'p1a'	! Physical space by default, configuration 1 by default
		Endif
		self%nf1a = field_count(1,1)
		self%nf2a = field_count(2,1)
		self%nf3a = field_count(3,1)
		self%nf3b = field_count(3,2)
		self%nf2b = field_count(2,2)
		self%nf1b = field_count(1,2)
		! Going from configuration 1 to configuration 2
		! Suppose we are going from r-in processor, modes distributed
		! --- immediately following the implicit solve		
		! -- to r-distributed, delta_lm(column) in processor
		! The message size I would send to rank r is
		!scount(r,1) = my_nl_lm*delta_r(r)*nfsend(1,2)		
		! where my_nl_lm is the number of l-m modes I hold during the solve
		! and nfsend is fields sent+radial derivatives sent
		! i.e. WSZ, and radial derivatives (P + radial derivatives do not need to be transposed)

		!Configuration 2 to 3
		!Next, suppose we are going from theta in processor, r distributed, m distributed
		! --- to --- phi in-processor, r-distributed, theta-distributed
		! The message size i would send to rank r is
		!scount23(r,1) = delta_r(my_row_rank)*delta_theta(column(r))*my_nm(local)*nfsend(2,3)

		np = pfi%rcomm%np
		Allocate(self%scount23(0:np-1))
		Allocate(self%rcount23(0:np-1))
		Allocate(self%sdisp23(0:np-1))
		Allocate(self%rdisp23(0:np-1))
		Do p = 0, pfi%rcomm%np-1
			self%scount23(p) = (pfi%my_1p%delta) * (pfi%all_2p(p)%delta) * (pfi%my_3s%delta) * (self%nf2a)
			self%rcount23(p) = (pfi%my_1p%delta) * (pfi%my_2p%delta) * (pfi%all_3s(p)%delta) * (self%nf2a)
		Enddo
		self%sdisp23(0) = 0
		self%rdisp23(0) = 0
		If (np .gt. 1) Then
			Do p = 1, np -1	
				self%sdisp23(p) = self%sdisp23(p-1)+self%scount23(p-1)
				self%rdisp23(p) = self%rdisp23(p-1)+self%rcount23(p-1)
			Enddo
		Endif
		If (present(report)) Then
			report_file = 'parallel_framework/reports/buff_init/sr23_'//Trim(report_tag)
			Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
			Write(report_unit)np
			Write(report_unit)(self%scount23(p),p=0,np-1)
			Write(report_unit)(self%sdisp23(p),p=0,np-1)
			Write(report_unit)(self%rcount23(p),p=0,np-1)
			Write(report_unit)(self%rdisp23(p),p=0,np-1)
			Close(report_unit)
		Endif

		!////////////////////////////////////////////////////////////
		! ------------ Now the reverse-----------------
		!Configuration 3 to 2
		! if nf2a and nf3b were the same, rcount32 would equal scount23, and scount32 would equal rcount23
		Allocate(self%scount32(0:np-1))
		Allocate(self%rcount32(0:np-1))
		Allocate(self%sdisp32(0:np-1))
		Allocate(self%rdisp32(0:np-1))
		Do p = 0, pfi%rcomm%np-1
			self%rcount32(p) = (pfi%my_1p%delta) * (pfi%all_2p(p)%delta) * (pfi%my_3s%delta) * (self%nf3b)
			self%scount32(p) = (pfi%my_1p%delta) * (pfi%my_2p%delta) * (pfi%all_3s(p)%delta) * (self%nf3b)
		Enddo
		self%sdisp32(0) = 0
		self%rdisp32(0) = 0
		If (np .gt. 1) Then
			Do p = 1, np -1	
				self%sdisp32(p) = self%sdisp32(p-1)+self%scount32(p-1)
				self%rdisp32(p) = self%rdisp32(p-1)+self%rcount32(p-1)
			Enddo
		Endif

		!///////////////////////////////////////////////////////
		!  ----- rcount32 and scount32 were set up as though
		!			the buffer was complex
		!        use these instead if it is real
		Allocate(self%scount32v2(0:np-1))
		Allocate(self%rcount32v2(0:np-1))
		Allocate(self%sdisp32v2(0:np-1))
		Allocate(self%rdisp32v2(0:np-1))

		self%scount32v2 = 2*self%scount32
		self%rcount32v2 = 2*self%rcount32
		self%sdisp32v2  = 2*self%sdisp32
		self%rdisp32v2  = 2*self%rdisp32

		If (present(report)) Then
			If (report == .true.) Then
			report_file = 'parallel_framework/reports/buff_init/sr32_'//Trim(report_tag)
			Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
			Write(report_unit)np
			Write(report_unit)(self%scount32(p),p=0,np-1)
			Write(report_unit)(self%sdisp32(p),p=0,np-1)
			Write(report_unit)(self%rcount32(p),p=0,np-1)
			Write(report_unit)(self%rdisp32(p),p=0,np-1)
			Close(report_unit)
			Endif
		Endif






		! Suppose we are going from 3 to 2 (reverse of above)
		!  The message size I would send to rank r is
		! scount(r,2) =  delta_r(my_row_rank)*delta_theta(column(my_rank))*nm(r)*nfsend(3,2)

		! Configuration 1 to 2
		np = pfi%ccomm%np
		Allocate(self%scount12(0:np-1))
		Allocate(self%rcount12(0:np-1))
		Allocate(self%sdisp12(0:np-1))
		Allocate(self%rdisp12(0:np-1))
		Do p = 0, pfi%ccomm%np-1
			self%scount12(p) = (pfi%all_1p(p)%delta) * (my_num_lm)  * (self%nf1a)*2
			self%rcount12(p) = (pfi%my_1p%delta)     * (num_lm(p)) * (self%nf1a)*2
		Enddo
		self%sdisp12(0) = 0
		self%rdisp12(0) = 0
		If (np .gt. 1) Then
			Do p = 1, np -1	
				self%sdisp12(p) = self%sdisp12(p-1)+self%scount12(p-1)
				self%rdisp12(p) = self%rdisp12(p-1)+self%rcount12(p-1)
			Enddo
		Endif


		! Configuration 2 to 1
		!  Suppose we are going from 2 to 1 (reverse direction)
		! The message size I would send to rank r is
		!scount(r,1) = nl_lm(r)*delta_r(my_row_rank)*nfsend(2,1)
		np = pfi%ccomm%np
		Allocate(self%scount21(0:np-1))
		Allocate(self%rcount21(0:np-1))
		Allocate(self%sdisp21(0:np-1))
		Allocate(self%rdisp21(0:np-1))
		Do p = 0, np-1
			self%rcount21(p) = (pfi%all_1p(p)%delta) * (my_num_lm)  * (self%nf2b)*2	! 2 is because we split the complex into two real pieces
			self%scount21(p) = (pfi%my_1p%delta)     * (num_lm(p)) * (self%nf2b)*2
		Enddo
		self%sdisp21(0) = 0
		self%rdisp21(0) = 0
		If (np .gt. 1) Then
			Do p = 1, np -1	
				self%sdisp21(p) = self%sdisp21(p-1)+self%scount21(p-1)
				self%rdisp21(p) = self%rdisp21(p-1)+self%rcount21(p-1)
			Enddo
		Endif



	End Subroutine Initialize_Spherical_Buffer

	Subroutine DeAllocate_Spherical_Buffer(self,config)
		Class(SphericalBuffer) :: self
		Character*3, Intent(In) :: config
		Integer :: mn1, mn2, mn3, mn4
		Integer :: mx1,mx2,mx3,mx4,i
		Select Case(config)
			Case('p1a')
				!Write(6,*)'De-allocating p1a'
				If (allocated(self%p1a)) Then
					DeAllocate(self%p1a)
				Else
					Write(6,*)'p1a does not appear to be allocated'
				Endif
			Case('p1b')
				!Write(6,*)'De-allocating p1b'
				If (allocated(self%p1b)) Then
					DeAllocate(self%p1b)
				Else
					Write(6,*)'p1b does not appear to be allocated'
				Endif
			Case('p2a')
				!Write(6,*)'De-allocating p2a'
				If (allocated(self%p2a)) Then
					DeAllocate(self%p2a)
				Else
					Write(6,*)'p2a does not appear to be allocated'
				Endif
			Case('p2b')
				!Write(6,*)'De-allocating p2b'
				If (allocated(self%p2b)) Then
					DeAllocate(self%p2b)
				Else
					Write(6,*)'p2b does not appear to be allocated'
				Endif

			Case('p3a')
				!Write(6,*)'De-allocating p3a'
				If (allocated(self%p3a)) Then
					DeAllocate(self%p3a)
				Else
					Write(6,*)'p3a does not appear to be allocated'
				Endif
			Case('p3b')
				!Write(6,*)'De-allocating p3b'
				If (allocated(self%p3b)) Then
					DeAllocate(self%p3b)
				Else
					Write(6,*)'p3b does not appear to be allocated'
				Endif
			Case('s2a')
				! Appropriate for a triangular truncation 
				If (allocated(self%s2a)) Then
					!Write(6,*)'De-allocating s2a'
					mn1 = pfi%my_3s%min
					mx1 = pfi%my_3s%max
					Do i = mn1, mx1
						If (Allocated(self%s2a(i)%data)) Then	
							DeAllocate(self%s2a(i)%data)
						Endif
					Enddo
					DeAllocate(self%s2a)
				Else
					Write(6,*)'s2a does not appear to be allocated.'
				Endif

			Case('s2b')
				! Appropriate for a triangular truncation 
				If (allocated(self%s2b)) Then
					!Write(6,*)'De-allocating s2b'
					mn1 = pfi%my_3s%min
					mx1 = pfi%my_3s%max
					Do i = mn1, mx1
						If (Allocated(self%s2b(i)%data)) Then	
							DeAllocate(self%s2b(i)%data)
						Endif
					Enddo
					DeAllocate(self%s2b)
				Else
					Write(6,*)'s2b does not appear to be allocated.'
				Endif

		End Select
	End Subroutine DeAllocate_Spherical_Buffer
	Subroutine Allocate_Spherical_Buffer(self,config)
		Class(SphericalBuffer) :: self
		Character*3, Intent(In) :: config
		Integer :: mn1, mn2, mn3, mn4
		Integer :: mx1,mx2,mx3,mx4
		Integer :: i
		Select Case(config)
			Case('p1a')
				mn1 = 1
				mx1 = pfi%n1p
				mn2 = 1
				mx2 = 2
				mn3 = 1
				mx3 = my_num_lm
				mn4 = 1
				mx4 = self%nf1a
				Allocate(self%p1a(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
			Case('p1b')
				mn1 = 1
				mx1 = pfi%n1p
				mn2 = 1
				mx2 = 2
				mn3 = 1
				mx3 = my_num_lm
				mn4 = 1
				mx4 = self%nf1b
				Allocate(self%p1b(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))


			Case('p2a')
				 
				mn1 = 1
				mx1 = pfi%n2p
				mn2 = pfi%my_1p%min
				mx2 = pfi%my_1p%max
				mn4 = pfi%my_3s%min
				mx4 = pfi%my_3s%max

				mx3 = self%nf2a*2*(mx2-mn2+1)
				Allocate(self%p2a(mn1:mx1, 1:mx3, mn4:mx4))
				
			Case('p2b')
				! 
				mn1 = 1
				mx1 = pfi%n2p
				mn2 = pfi%my_1p%min
				mx2 = pfi%my_1p%max
				mn4 = pfi%my_3s%min
				mx4 = pfi%my_3s%max

				mx3 = self%nf2b*2*(mx2-mn2+1)
				Allocate(self%p2b(mn1:mx1, 1:mx3, mn4:mx4))		
			Case('p3a')
				mn1 = 1
				mx1 = pfi%n3p+2
				mn2 = pfi%my_1p%min
				mx2 = pfi%my_1p%max
				mn3 = pfi%my_2p%min
				mx3 = pfi%my_2p%max
				mn4 = 1
				mx4 = self%nf3a
				! might think of calling this p3a rather than rdata 3a
				Allocate(self%p3a(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
			Case('p3b')
				mn1 = 1
				mx1 = pfi%n3p+2 ! necessary for an in place transform
				mn2 = pfi%my_1p%min
				mx2 = pfi%my_1p%max
				mn3 = pfi%my_2p%min
				mx3 = pfi%my_2p%max
				mn4 = 1
				mx4 = self%nf3b
				!Write(6,*)'Allocating p3b'
				Allocate(self%p3b(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))


			Case('s2a')
				! We use a real array here instead of a complex array 
				! Fields are striped rmin-rmax real then rmax+1-2*rmax imaginary in second index
				! Appropriate for a triangular truncation 
				! --- partly why spherical and cartesian buffers will be separate
				mn1 = pfi%my_3s%min
				mx1 = pfi%my_3s%max
				mx4  = self%nf2a
				mn3 = pfi%my_1p%min
				mx3 = pfi%my_1p%max
				
				Allocate(self%s2a(mn1:mx1))
				mx2 = maxval(pfi%inds_3s)	! l_max = m_max

				!/////////////////////////////////


				mx3 = self%nf2a*2*pfi%my_1p%delta

				!////////////////////////////////

				Do i = mn1, mx1
					mn2 = pfi%inds_3s(i)		!l_min = m
					Allocate(self%s2a(i)%data(mn2:mx2,1:mx3))
				Enddo
			Case('s2b')
				! This is just s2b, but we are no longer complex, and 
				! have fields, imaginary/real, and radius all in second index
				! Appropriate for a triangular truncation 
				! --- partly why spherical and cartesian buffers will be separate
				mn1 = pfi%my_3s%min
				mx1 = pfi%my_3s%max
				mx4  = self%nf2b
				mn3 = pfi%my_1p%min
				mx3 = pfi%my_1p%max
				
				Allocate(self%s2b(mn1:mx1))
				mx2 = maxval(pfi%inds_3s)	! l_max = m_max

				!/////////////////////////////////


				mx3 = self%nf2b*2*pfi%my_1p%delta

				!////////////////////////////////

				Do i = mn1, mx1
					mn2 = pfi%inds_3s(i)		!l_min = m
					Allocate(self%s2b(i)%data(mn2:mx2,1:mx3))
				Enddo


		End Select

	End Subroutine Allocate_Spherical_Buffer

	Subroutine Advance_Configuration(self)
		Class(SphericalBuffer) :: self		
		Select Case(self%config)
			Case ('p2a')
				Call self%transpose_2a3a()
			Case ('p3b')
				Call self%transpose_3b2b()
			Case ('s2b')
				Call self%transpose_2b1b()
			Case ('p1a')
				Call self%transpose_1a2a
		End Select
	End Subroutine Advance_Configuration	

	Subroutine Write_Space(self,extra_tag)
		Class(SphericalBuffer) :: self
		character*120, Optional, Intent(In) :: extra_tag
		Character*120 :: report_file,report_tag
		Character*10 :: gtag, rtag, ctag
		Integer :: report_unit = 500
		Integer :: i,j,k,f
		Integer :: imin,imax,jmin,jmax,kmin,kmax,nf
		Write(gtag,ifmt)pfi%gcomm%rank
		Write(rtag,ifmt)pfi%rcomm%rank
		Write(ctag,ifmt)pfi%ccomm%rank
		report_tag = 'g'//Trim(gtag)//'_r'//Trim(rtag)//'_c'//Trim(ctag)
		Select Case(self%config)

			Case ('p3a')		
				imin = 1
				imax = pfi%n3p
				jmin = pfi%my_1p%min
				jmax = pfi%my_1p%max
				nf = self%nf3a
				kmin = pfi%my_2p%min
				kmax = pfi%my_2p%max
				if (present(extra_tag)) Then

					report_file = 'parallel_framework/reports/workspace/p3a_'//Trim(report_tag)//'_'//Trim(extra_tag)
				else		
					report_file = 'parallel_framework/reports/workspace/p3a_'//report_tag
				endif
				Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
				Write(report_unit)pfi%n3p
				Write(report_unit)pfi%my_1p%delta
				Write(report_unit)pfi%my_2p%delta
				Write(report_unit)nf
				Write(report_unit)((((self%p3a(i,j,k,f),i=imin,imax),j=jmin,jmax),k = kmin,kmax),f = 1, nf)
				Write(report_unit)pfi%my_1p%min
				Write(report_unit)pfi%my_1p%max
				Write(report_unit)pfi%my_2p%min
				Write(report_unit)pfi%my_2p%max
				Close(report_unit)
			Case ('p3b')		
				imin = 1
				imax = pfi%n3p
				jmin = pfi%my_1p%min
				jmax = pfi%my_1p%max
				nf = self%nf3a
				kmin = pfi%my_2p%min
				kmax = pfi%my_2p%max
				if (present(extra_tag)) Then

					report_file = 'parallel_framework/reports/workspace/p3b_'//Trim(report_tag)//'_'//Trim(extra_tag)
				else		
					report_file = 'parallel_framework/reports/workspace/p3b_'//report_tag
				endif
				Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
				Write(report_unit)pfi%n3p
				Write(report_unit)pfi%my_1p%delta
				Write(report_unit)pfi%my_2p%delta
				Write(report_unit)nf
				Write(report_unit)((((self%p3b(i,j,k,f),i=imin,imax),j=jmin,jmax),k = kmin,kmax),f = 1, nf)
				Write(report_unit)pfi%my_1p%min
				Write(report_unit)pfi%my_1p%max
				Write(report_unit)pfi%my_2p%min
				Write(report_unit)pfi%my_2p%max
				Close(report_unit)
		End Select
	End Subroutine Write_Space




	Subroutine Transpose_2a3a(self)
		Class(SphericalBuffer) :: self
		Real*8, Allocatable :: send_buff(:), recv_buff(:)
		Integer :: send_size,np, recv_size
		Integer :: imin, imax, jmin, jmax, kmin,kmax,ii,nf
		Integer :: i,f,j,p,k,k_ind,delf,delj
		! This is where we we move from theta, delta_r, delta_m 
		!  to m, delta_r, delta_theta
		send_size = sum(self%scount23)*2
		recv_size = sum(self%rcount23)*2
		Allocate(send_buff(1:send_size))

		!--- Not sure if this is good or bad, but copy out the bounds of the loop for now
		nf = self%nf2a
		kmin = pfi%my_3s%min
		kmax = pfi%my_3s%max
		!jmin = pfi%my_1p%min
		!jmax = pfi%my_1p%max

		jmin = 1
		jmax = pfi%my_1p%delta

		np = pfi%rcomm%np
		! Possibly better ways to stripe, but for now, we will stripe each processors data all at once
		! This means the the send buffer is accessed "naturally," but that the p2a array
		! gets jumped around in
		delj = pfi%my_1p%delta

		ii = 1
		Do p = 0, np -1 
			imin = pfi%all_2p(p)%min
			imax = pfi%all_2p(p)%max	
			Do k = kmin, kmax
			Do f = 1,nf
			delf = (f-1)*delj*2
			Do j = jmin,jmax
			Do i = imin,imax

				send_buff(ii) = self%p2a(i,j+delf,k)
				send_buff(ii+1) = self%p2a(i,j+delf+delj,k)
				ii = ii+2
			Enddo
			Enddo
			Enddo
			Enddo
		Enddo

		Call self%deconstruct('p2a')
		Allocate(recv_buff(1:recv_size))
		recv_buff(:) = 0.0d0
		!----- This is where alltoall will be called
		Call Standard_Transpose(send_buff, recv_buff, self%scount23*2, self%sdisp23*2, self%rcount23*2, self%rdisp23*2, pfi%rcomm)
		!--------------------------------------------------
		DeAllocate(send_buff)
		Call self%construct('p3a')
		self%p3a(:,:,:,:) = 0.0d0	! This is important because we are going to take an fft later (De-aliasing is implicit here because
		! we only stripe in data of the de-aliased m's, but we need to make sure the higher m's are zero!
		!Stripe from the receive buff
		!Here we access the receive buffer in the 'natural' order
		imin = pfi%my_2p%min
		imax = pfi%my_2p%max	
		jmin = pfi%my_1p%min
		jmax = pfi%my_1p%max
		ii = 1
		Do p = 0, np -1
			kmin = pfi%all_3s(p)%min
			kmax = pfi%all_3s(p)%max	
			Do k = kmin, kmax
			Do f = 1,nf
			Do j = jmin,jmax
			Do i = imin,imax

				k_ind = pfi%inds_3s(k)*2+1  ! (real) m=0 stored in p3b(1,:,:,:) (img in p3b(2,:,:,:))

				self%p3a(k_ind,j,i,f)=recv_buff(ii) 
				self%p3a(k_ind+1,j,i,f)=recv_buff(ii+1) 			

				ii = ii+2

			Enddo
			Enddo
			Enddo
			Enddo
		Enddo
		self%config = 'p3a'
		DeAllocate(recv_buff)


	End Subroutine Transpose_2a3a

	Subroutine Transpose_3b2b(self)		
		! Version 2
		! This assumes that s3b and p2b are actually real arrays
		! p3b was fft'd in place
		Class(SphericalBuffer) :: self
!		Complex*16, Allocatable :: send_buff(:), recv_buff(:)
		Real*8, Allocatable :: send_buff(:), recv_buff(:)
		Integer :: send_size,np, recv_size
		Integer :: imin, imax, jmin, jmax, kmin,kmax,ii,nf
		Integer :: i,f,j,p,k,k_ind
		Integer :: delf, delj
		! This is where we we move from theta, delta_r, delta_m 
		!  to m, delta_r, delta_theta
		send_size = sum(self%scount32v2)
		recv_size = sum(self%rcount32v2)
		Allocate(send_buff(1:send_size))
		!write(6,*)'executing new transpose'
		!--- Not sure if this is good or bad, but copy out the bounds of the loop for now
		nf = self%nf3b

		jmin = pfi%my_1p%min
		jmax = pfi%my_1p%max

		np = pfi%rcomm%np


		!///////////////////
		!  Again stripe in the natural order of the send buffer
		imin = pfi%my_2p%min
		imax = pfi%my_2p%max	
		ii = 1
		Do p = 0, np -1
			kmin = pfi%all_3s(p)%min 
			kmax = pfi%all_3s(p)%max	
			Do k = kmin, kmax
			Do f = 1,nf
			Do j = jmin,jmax
			Do i = imin,imax		! interleave real and imaginary parts
				k_ind = pfi%inds_3s(k)*2+1  ! (real) m=0 stored in p3b(1,:,:,:) (img in p3b(2,:,:,:))

				send_buff(ii) = self%p3b(k_ind,j,i,f)! real part
				send_buff(ii+1) = self%p3b(k_ind+1,j,i,f)! complex part				

				ii = ii+2
			Enddo
			Enddo
			Enddo
			Enddo
		Enddo

		!/////////////////////////////////////
		Call self%deconstruct('p3b')
		Allocate(recv_buff(1:recv_size))
		recv_buff(:) = 0.0d0
		!----- This is where alltoall will be called
		Call Standard_Transpose(send_buff, recv_buff, self%scount32v2, self%sdisp32v2, self%rcount32v2, self%rdisp32v2, pfi%rcomm)
		!--------------------------------------------------
		DeAllocate(send_buff)
	
		Call self%construct('p2b')		! p2a and p2b can share the same buffer space... maybe just call this p2...

		delj = pfi%my_1p%delta
		

		!///////////////////////////////////////
		! Read from the receive buffer in its natural order
		kmin = pfi%my_3s%min
		kmax = pfi%my_3s%max

		jmin = 1
		jmax = pfi%my_1p%delta

		ii = 1
		Do p = 0, np -1 
			imin = pfi%all_2p(p)%min
			imax = pfi%all_2p(p)%max	
			Do k = kmin, kmax
			Do f = 1,nf
				delf = (f-1)*delj*2
			Do j = jmin,jmax
			Do i = imin,imax
				!p2b needs to be reshaped
				!   self%p2b(i,j*2*f,k) -- since dgemm needs a 2D array

				self%p2b(i,j+delf ,k)   = recv_buff(ii)  ! real
				self%p2b(i,j+delf+delj,k) = recv_buff(ii+1)	! complex
				ii = ii+2 ! added +2 
			Enddo
			Enddo
			Enddo
			Enddo
		Enddo



		self%config = 'p2b'
		DeAllocate(recv_buff)
	End Subroutine Transpose_3b2b



    Subroutine Transpose_2b1b(self)
      ! Go from Explicit_Part(IDX) configuration to the new Implicit%RHS configuration
      ! Communication is now done entirely within the radial group (processors that 
      ! share a common set of ell-m values).  
      Implicit None

      Integer :: r,l,m, mp, lp, indx, rank, rrank, r_min, r_max, dr, mm, arank, cnt,i
		Integer :: n1, n, nfields, offset, delta_r, rmin, rmax, n_lm_local, np,p


		Real*8, Allocatable :: send_buff(:),recv_buff(:)
      Integer :: send_size,recv_size, lmax, tnr, send_offset

		Class(SphericalBuffer) :: self

		n1 = pfi%n1p

      !Allocate(send(my_r%min:my_r%max, 2, n_fields, N_lm_local))
		nfields = self%nf2b

		!/////
		! First loading, should be relatively straight forward
		! We have nfields (nf2b)  and the data is dimensioned x%(l,r,field_num)

      ! Load the send array in the ordering of the l-m values

		delta_r = pfi%my_1p%delta
		!rmin = pfi%my_1p%min
		!rmax = pfi%my_1p%max
		!modification for new configuration
		rmin = 1
		rmax = 2*delta_r
		tnr = 2*delta_r
		lmax = maxval(pfi%inds_3s)
		n_lm_local = 0
		Do i = pfi%my_3s%min, pfi%my_3s%max
			n_lm_local = n_lm_local+(lmax-pfi%inds_3s(i)+1)
		Enddo
		send_size = (pfi%my_1p%delta)*2*(nfields)*(n_lm_local)
		Allocate(send_buff(1:send_size))
		send_offset = 0
      Do i = 1, lm_count
         mp = mp_lm_values(i)
         l = l_lm_values(i)
			offset = 0
			Do n = 1, nfields

         	Do r = rmin,rmax
					send_buff(send_offset+r) = self%s2b(mp)%data(l,offset+r)
         	End Do
				send_offset = send_offset+tnr
				offset = offset+tnr
			Enddo
      Enddo



		Call self%deconstruct('s2b')

		recv_size = (pfi%n1p)*2*(nfields)*num_lm(pfi%ccomm%rank)
		Allocate(recv_buff(1:recv_size))


		Call Standard_Transpose(send_buff, recv_buff, self%scount21, self%sdisp21, self%rcount21, self%rdisp21, pfi%ccomm)
		DeAllocate(send_buff)

		Call self%construct('p1b')
		! Now, the receive striping needs a little mapping
		! WPS are coupled, but Z,Btor, and Bpol are not
		! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
		! We may want to modify this later on to mesh with linear equation structure
      indx = 1
		np = pfi%ccomm%np
      Do p = 0, np - 1
         !r_min = 1 + p*n1/np
			r_min = pfi%all_1p(p)%min
			r_max = pfi%all_1p(p)%max
         !r_max = (p+1)*n1/np
         dr = r_max - r_min
         ! Each processor in the radial group will have given me the same number of 
         ! l-m combos in the same (correct) order
         cnt = 1
         Do lp = 1, my_num_lm
				Do n = 1, nfields
					self%p1b(r_min:r_max,1,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
					self%p1b(r_min:r_max,2,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
				Enddo
				cnt = cnt+1
			Enddo
      End Do

		self%config='p1b'
      Deallocate(recv_buff)
		
    End Subroutine Transpose_2b1b



    Subroutine Transpose_1a2a(self)
      ! Go from implicit configuration (1 physical) to configuration 2 (spectral)
      Implicit None

      Integer :: r,l,m, mp, lp, indx, rank, rrank, r_min, r_max, dr, mm, arank, cnt,i
		Integer :: n1, n, nfields, offset, delta_r, rmin, rmax, n_lm_local, np,p
		Integer :: recv_offset, tnr

		Real*8, Allocatable :: send_buff(:),recv_buff(:)
      Integer :: send_size,recv_size, lmax

		Class(SphericalBuffer) :: self


		! This goes at the beginning
		
		nfields = self%nf1a
		send_size = (pfi%n1p)*2*(nfields)*num_lm(pfi%ccomm%rank)
		Allocate(send_buff(1:send_size))


		! Now, the receive striping needs a little mapping
		! WPS are coupled, but Z,Btor, and Bpol are not
		! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
		! We may want to modify this later on to mesh with linear equation structure
      indx = 1
		np = pfi%ccomm%np
      Do p = 0, np - 1
         !r_min = 1 + p*n1/np
			r_min = pfi%all_1p(p)%min
			r_max = pfi%all_1p(p)%max
         !r_max = (p+1)*n1/np
         dr = r_max - r_min
         ! Each processor in the radial group will have given me the same number of 
         ! l-m combos in the same (correct) order
         cnt = 1
         Do lp = 1, my_num_lm
				Do n = 1, nfields
					!self%p1b(r_min:r_max,1,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
					!self%p1b(r_min:r_max,2,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
					send_buff(indx:indx+dr) = self%p1a(r_min:r_max,1,cnt,n) ; indx = indx + dr + 1
					send_buff(indx:indx+dr) = self%p1a(r_min:r_max,2,cnt,n) ; indx = indx + dr + 1
				Enddo
				cnt = cnt+1
			Enddo
      End Do
		Call self%deconstruct('p1a')



		!////


		!/////
		! First loading, should be relatively straight forward
		! We have nfields (nf2b)  and the data is dimensioned x%(l,r,field_num)

      ! Load the send array in the ordering of the l-m values
      indx = 1

		delta_r = pfi%my_1p%delta
		rmin = pfi%my_1p%min
		rmax = pfi%my_1p%max
		lmax = maxval(pfi%inds_3s)
		n_lm_local = 0
		Do i = pfi%my_3s%min, pfi%my_3s%max
			n_lm_local = n_lm_local+(lmax-pfi%inds_3s(i)+1)
		Enddo
		recv_size = (pfi%my_1p%delta)*2*(nfields)*(n_lm_local)
		Allocate(recv_buff(1:recv_size))

		!Write(6,*)'maxval sendbuff : ', maxval(abs(send_buff))
		Call Standard_Transpose(send_buff, recv_buff, self%scount12, self%sdisp12, self%rcount12, self%rdisp12, pfi%ccomm)
		DeAllocate(send_buff)

		Call self%construct('s2a')

		rmin = 1
		rmax = delta_r*2
		recv_offset = 0
		tnr = 2*delta_r
      Do i = 1, lm_count
         mp = mp_lm_values(i)
         l = l_lm_values(i)
			offset = 0
			Do n = 1, nfields
         Do r = rmin,rmax 
            !self%s2a(mp)%data(l,r,n) = Cmplx(recv_buff(offset+indx),recv_buff(offset+indx+delta_r)) !  I am keeping this here as a reminder  ,Precision)
				!indx = indx+1
				self%s2a(mp)%data(l,offset+r) = recv_buff(recv_offset+r)
         End Do
			offset = offset+tnr
			recv_offset = recv_offset+tnr
			Enddo

         indx = indx+1
      Enddo
		!////
		!Write(6,*)'maxval recv : ', maxval(abs(recv_buff))
		self%config='s2a'
      Deallocate(recv_buff)
		
    End Subroutine Transpose_1a2a



	Subroutine Initialize_Parallel_Interface(self, pars)
		Integer, Intent(In) :: pars(1:)
		Integer :: pcheck, error
		Class(Parallel_Interface) :: self	
		self%geometry = pars(1)
		self%n1p = pars(2)
		self%n1s = pars(3)
		self%n2p = pars(4)
		self%n2s = pars(5)
		self%n3p = pars(6)
		self%n3s = pars(7)
		pcheck = size(pars,1)
		If (pcheck .gt. 7) Then
			self%npe = pars(8)
			self%nprow = pars(9)
			self%npcol = pars(10)
			self%npio = pars(8)-pars(9)*pars(10)
			self%npc = pars(10)*pars(9)
		Else
			! get the process grid info from the environment
		Endif
		self%gcomm = init_main_group(error)
		if (self%gcomm%np .ne. self%npe) Then
			Write(6,*)'Error np does not agree with number of processes.'
			Write(6,*)'gcomm np : ', self%gcomm%np
			Write(6,*)'specified np : ', self%npe
		Endif

		Call rowcolsplit(self%gcomm,self%rcomm,self%ccomm,self%nprow,self%npcol,error)
		Call self%Init_Geometry()
		!Write(6,*)'Hello ', self%gcomm%rank, self%rcomm%rank, self%ccomm%rank
		!Write(6,*)'Hello: ', self%ccomm%rank, self%ccomm%rank, self%my_1p%min, self%my_1p%max, self%my_1p%delta
		!Write(6,*)'Hello: ', self%rcomm%rank,self%ccomm%rank, self%my_2p%min, self%my_2p%max, self%my_2p%delta		

	End Subroutine Initialize_Parallel_Interface

	Subroutine Init_Geometry(self)
		Class(Parallel_Interface) :: self		
		If (self%geometry == Spherical) Then
			Call self%Spherical_Init()
		Endif
	End Subroutine Init_Geometry


	Subroutine Spherical_Init(self)
		Integer :: r, unit1
		Class(Parallel_Interface) :: self	

		!Distribute radii amongst members of a column 	
		Allocate(self%all_1p(0:self%npcol-1))
		Call Standard_Balance(self%all_1p,self%n1p,self%ccomm)
		self%my_1p = self%all_1p(self%ccomm%rank)

		!Distribute theta amongst members of a row 	
		Allocate(self%all_2p(0:self%nprow-1))
		Call Standard_Balance(self%all_2p,self%n2p,self%rcomm)
		self%my_2p = self%all_2p(self%rcomm%rank)

		! Determine number of m-values per process
		
		Allocate(self%all_3s(0:self%nprow-1))
		Call Standard_Balance(self%all_3s,self%n3s,self%rcomm)
		self%my_3s = self%all_3s(self%rcomm%rank)
		!Write(6,*), self%my_3s%min, self%my_3s%max,self%rcomm%rank, self%my_3s%delta
		! We assume a high-low pairing of m-values.
		! This means that we need to set up an index array
		! for the m-values that tells how many each processor has.		
		Allocate(self%inds_3s(1:self%n3s))
		Call m_balance(self%all_3s, self%inds_3s, self%rcomm)
		Call LM_Load_Balance(pfi%my_3s,self%inds_3s,self%ccomm)
		If (self%gcomm%rank .eq. -100) Then
						unit1 = 10
			         Open(unit1, file = 'verification/m_check', status='replace', form = 'unformatted')
						Write(unit1) self%n3s
						Write(unit1) (self%inds_3s(r),r=1,self%n3s)
			
						Close(unit1)
		Endif
	End Subroutine Spherical_Init



	Subroutine Finalize_Framework(self)
		Class(Parallel_Interface) :: self
		Integer :: error
		Call Exit_Comm_Lib(error)	
		
	End Subroutine Finalize_Framework

End Module Parallel_Framework
