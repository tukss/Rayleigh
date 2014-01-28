Module Field_Manipulation
	Implicit None

	!============================================================================
	!  		The Field Object
	!			Can contain a field or a series of fields
	!			Primary methods are transpose and transform
	Type, Public :: Field
		Integer :: Geometry =1 			! 1 for Cartesian, 2 for Cylindrical, 3 for Spherical		
		Integer :: Config_Order(1:3) 
		Integer :: Configuration = 1	! Changes as the field is transposed
		Logical :: Physical = .true.	! Determines if we are in real space
		Logical :: Spectral = .false.	! or spectral space
		Integer*8 :: fftw_plan1, fftw_plan2, fftw_plan3	! FFT plans for each configuration space
		Real*8, Pointer, dimension(:,:,:) :: rdata
		Complex*16, Pointer, dimension(:,:,:) :: cdata
		!----------- Size of field (total - not distributed)
		Integer, Target :: nx1, nx2, nx3, nx1s, nx2s, nx3s

		!---------- Indices corresponding to each dimension
		Integer, Target :: x1min, x1max, x1smin, x1smax
		Integer, Target :: x2min, x2max, x2smin, x2smax
		Integer, Target :: x3min, x3max, x3smin, x3smax

		Contains

		!Procedure, Public :: transpose => Transpose_Field
		Procedure, Public :: advance_config_serial => Permute_Forward
		!Procedure, Public :: reverse_config_serial => Permute_Reverse
		Procedure, Public :: print_dims => print_rdims
	End Type Field

	Type, Public :: Spherical_Field
		
		
	End Type Spherical_Field
	

Contains
	Subroutine print_rdims(self)
		Class(Field) :: self
		Real*8, Pointer, dimension(:,:,:) :: temp
		nullify(temp)
		temp => self%rdata
		Write(6,*)size(temp,1), size(temp,2), size(temp,3)
	End Subroutine print_rdims

	Subroutine Permute_Indices(field_rpointer)
		! Simple transpose
		! Permutes (x,y,z) ordering to (y,z,x) ordering
		Real*8, Pointer, dimension(:,:,:) :: field_rpointer
		Real*8, Pointer, dimension(:,:,:) :: field_copy
		Integer :: ilb, iub, jlb, jub, klb, kub
		Integer :: i,j,k,ni, nj, nk
		ilb = Lbound(field_rpointer,1)
		iub = Ubound(field_rpointer,1)
		jlb = Lbound(field_rpointer,2)
		jub = Ubound(field_rpointer,2)
		klb = Lbound(field_rpointer,3)
		kub = Ubound(field_rpointer,3)


		Nullify(field_copy)
		field_copy => field_rpointer
		nullify(field_rpointer)
		Allocate(field_rpointer(jlb:jub,klb:kub,ilb:iub))		
		Do i = ilb, iub
			Do k = klb, kub

				Do j = jlb, jub
					field_rpointer(j,k,i) = field_copy(i,j,k)
				Enddo
			Enddo
		Enddo
		DeAllocate(field_copy)
		Nullify(field_copy)
	End Subroutine Permute_Indices

	Subroutine Permute_Forward(self)
		! Simple transpose
		! Permutes (x,y,z) ordering to (y,z,x) ordering
		Class(Field) :: self
		Real*8, Pointer, dimension(:,:,:) :: field_rpointer
		Real*8, Pointer, dimension(:,:,:) :: field_copy
		Integer :: ilb, iub, jlb, jub, klb, kub
		Integer :: i,j,k,ni, nj, nk
		nullify(field_rpointer)
		field_rpointer => self%rdata
		ilb = Lbound(field_rpointer,1)
		iub = Ubound(field_rpointer,1)
		jlb = Lbound(field_rpointer,2)
		jub = Ubound(field_rpointer,2)
		klb = Lbound(field_rpointer,3)
		kub = Ubound(field_rpointer,3)


		Nullify(field_copy)
		field_copy => self%rdata
		nullify(self%rdata)
		Allocate(self%rdata(jlb:jub,klb:kub,ilb:iub))		
		Do i = ilb, iub
			Do k = klb, kub

				Do j = jlb, jub
					self%rdata(j,k,i) = field_copy(i,j,k)
				Enddo
			Enddo
		Enddo
		DeAllocate(field_copy)
		Nullify(field_copy)
	End Subroutine Permute_Forward

	Subroutine Transform(self)
		Implicit None
		Class(Field) :: self 
		
		Call dfftw_plan_many_dft_r2c (planrc,1,info%length,Size(x,Dim=2)*Size(x,Dim=3),&
			& x,0,1,2*(info%length/2+1),&
			& x,0,1,info%length/2+1,FFTW_ESTIMATE)
    	Call dfftw_execute(planrc)
    	Call dfftw_destroy_plan (planrc)

	End Subroutine Transform(self)

	Subroutine TransformNew(self)
		! Forget the FFT for now - let's just set up the basic logic to change sizes

		If (self%configuration .eq 3) Then ! we're in physical space
			nx3s = nx3/2
			Allocate(self%cdata(1:nx3s,1:nx2,1:nx1))
			! Call the transform to
			DeAllocate(self%rdata)
		Endif
		If (self%configuration .eq. 2) Then	! We are in the hybrid space
				field_copy => self%cdata
				nullify(self%rdata)
				Allocate(self%cdata(1:nx2s,1:nx3s,1:nx1))		
		Endif		
		! This was scratchwork  N.F. Wednesday, March 20, 2013
	End Subroutine

	Subroutine Transpose_Scratch()
		Integer :: i, j, k
		Integer :: nm, delta_r, delta_theta
		Integer :: chunk_size
		! Suppose we are striped phi, r, theta (x3,delta_x1,delta_x2)
		! and we want to go to (x2, delta_x3, delta_x1)

		chunk_size = delta_r*delta_theta
		!======= This is like ASH
		

		! Striping method 1
		! Loop over the (send) buffer in natural order.  Place elements in the receive buffer
		! in an unnatural order
		! this takes (m, delta_r, delta_theta) to (delta_r,delta_theta, m)
		do k = 1, delta_theta
			do j = 1, delta_r
				do i = 1, nm
					offset = (i-1)*chunk_size + (k-1)*delta_r+j					
					send_buffer(offset) = field_buffer(i,j,k)
				enddo
			enddo
		enddo

		!Alternative method
		! take (m,delta_r,delta_theta) to (delta_theta,delta_r,m)
		do k = 1, delta_theta
			do j = 1, delta_r
				do i = 1, nm
					offset = (i-1)*chunk_size + (j-1)*delta_theta+k					
					send_buffer(offset) = field_buffer(i,j,k)
				enddo
			enddo
		enddo

		! The two examples above are really only appropriate when processor 0 has m1,m2,m3 etc
		! They are not so great for when we do high/low m pairing.
		! This requires a simple modification
		! We have a low-high ordered list of m values that are ordered by how they are distributed
		! amongst the processors.  For instance, suppose processor 0 has m1,m8,m2,m7 and processor
		! 1 has m3,m6,m4,m5
		! then we make an array named ind
		! with ind(1) = 1, ind(8) = 2, ind(2) = 3, etc
		do k = 1, delta_theta
			do j = 1, delta_r
				do i = 1, nm
					offset = (ind(i)-1)*chunk_size + (j-1)*delta_theta+k					
					send_buffer(offset) = field_buffer(i,j,k)
				enddo
			enddo
		enddo		

		!========================================================================================================
		!  Now we've done the alltoallv
		!  The receive buffer looks like:  (data from p0, data from p1, data from p2)
		!  This is striped as either:   ( [delta_r,delta_theta,m_local]_0, [delta_r,delta_theta,m_local]_1, etc)
		!                               or
		!										  ( [delta_theta,delta_r,m_local]_0, [delta_theta,delta_r,m_local]_1, etc)
		!  Note that when the number of ell-values varies for each m, we want m to run the slowest, so that we can group
		!		all the data for 1 m value together (making the legendre transform more efficient)
		!  That is why m is the slowest index, and not r
		!  For cartesian, this shouldn't matter since the fft operation doesn't have similar contraints.
		!    --- though we'll probably keep this format for cartesian at first just to make the rest of the code jive.
		!
		!  When we leave the transpose, we want a configuration that looks like (theta, delta_r, m_local)   ***m_local means delta_m
		!  At this stage, it appears that the second option is preferable - so assume the data is striped [delta_theta,delta_r, delta_m]_i
		do p = 1, np  ! loop over all the processors (who each send a different set of thetas
			do m = 1, nm_local
				do r = 1, delta_r
					field_buffer(tmin(p):tmax(p),r,m) = receive_buffer(ind:ind+tmax(p)-1)		! Field buffer will have been reallocated/resized
					ind = ind+dtheta(p)																			! It might even have a different structure alltogether (like ASH)
				enddo
		enddo			
		!
		!==================================================================================

		!========================================
		!  Finally, there is one caveat.  We would like our buffer to contain multiple fields.
		!  The send buffer really should have been striped as [delta_theta, delta_r * nfields, delta_m]
		!  Let us assume that the typical field buffer layout is given by
		!  (nphi,delta_r, delta_theta*nfields)  in physical space.
		!  This will allow us to construct efficient points vr, vp for instance that point to pieces of this buffer
		!  e.g. vr => buff, vp => buff.  vrvp = vr*vp would be efficient
		! this is a minor modification to the striping

		!  Maybe actually, field should just be a 4-D buffer
		chunk_size = delta_r*delta_theta*nfields		! modify chunk size
		small_chunk = delta_r*delta_theta
		do nn = 1, nfields
		do k = 1, delta_theta
			do j = 1, delta_r
				do i = 1, nm
					offset = (ind(i)-1)*chunk_size + (j-1)*delta_theta+k	+ (nn-1)*small_chunk			
					send_buffer(offset) = field_buffer(i,j,k,n)
				enddo
			enddo
		enddo			
		enddo

		!  When we unload this, we just need to assume delta_r is really delta_rnf = delta_r*nfields
		do p = 1, np  ! loop over all the processors (who each send a different set of thetas
			do m = 1, nm_local
				do r = 1, delta_rnf
					field_buffer(tmin(p):tmax(p),r,m) = receive_buffer(ind:ind+tmax(p)-1)		! Field buffer will have been reallocated/resized
					ind = ind+dtheta(p)																			! It might even have a different structure alltogether (like ASH)
				enddo
		enddo	
		! The field buffer in this configuration will be dimensioned (theta, delta_r*nfields, nm_local)

	End Subroutine Transpose_Scratch

	Subroutine Reverse_Angular_Transpose()
		! Now that I've written the forward (physical to spectral) transpose of theta, let's code the reverse transpose
		! Reverse is pretty easy, we use the same loop, but switch the order

		! We want to take (theta,delta_rnf, delta_m) and pack it as ( [delta_theta0,delta_rnf,delta_m], [delta_theta1,delta_rnf,delta_m] )

		do p = 1, np  ! loop over all the processors 
			do m = 1, nm_local
				do r = 1, delta_rnf
					receive_buffer(ind:ind+tmax(p)-1) = field_buffer(tmin(p):tmax(p),r,m)
					ind = ind+dtheta(p)
				enddo
		enddo	

		! Call the alltoall

		!  The receive buffer is striped [(delta_theta,delta_rnf,delta_m_0),(delta_theta,delta_rnf,delta_m1)....]
		!  We want to change it to [m,delta_r, delta_theta*nfields] (our buffer config before transforming m => phi)

		! If we've done this properly, we should be able to just reverse the logic from the other loop
		! --- NEED TO VERIFY THIS
		!-------------------------------
		do k = 1, delta_theta
			do j = 1, delta_r
				do i = 1, nm
					offset = (ind(i)-1)*chunk_size + (j-1)*delta_theta+k					
					field_buffer(i,j,k) = send_buffer(offset)
				enddo
			enddo
		enddo		


	End Subroutine Reverse_Angular_Tranpose

End Module Field_Manipulation
