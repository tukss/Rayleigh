    Subroutine Transpose_2b1b(self)
      ! Go from Explicit_Part(IDX) configuration to the new Implicit%RHS configuration
      ! Communication is now done entirely within the radial group (processors that 
      ! share a common set of ell-m values).  
      Implicit None

      Integer :: r,l, mp, lp, indx, r_min, r_max, dr,  cnt,i, imi
		Integer :: n1, n, nfields, offset, delta_r, rmin, rmax, np,p

		Real*8, Allocatable :: send_buff(:),recv_buff(:)
      Integer :: tnr, send_offset

        ! piggyback information
        Integer :: inext, pcurrent

		Class(SphericalBuffer) :: self

		n1 = pfi%n1p
        np = pfi%ccomm%np
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

		Allocate(send_buff(1:self%send_size21))


    !//  pcurrent = 0
    !//  inext = num_lm(0)
	!//	send_offset = 0
    !//  Do i = 1, lm_count
    !//     mp = mp_lm_values(i)
    !//     l = l_lm_values(i)
	!//		offset = 0 

	!//		Do n = 1, nfields

    !//     	Do r = rmin,rmax
	!//				send_buff(send_offset+r) = self%s2b(mp)%data(l,offset+r)
    !//     	End Do
	!//			send_offset = send_offset+tnr
	!//			offset = offset+tnr
	!//		Enddo
	!//		If (i .eq. inext) Then
	!//        If (self%do_piggy) Then
    !//        ! load the piggyback value into the next buffer slot.
    !//            send_buff(send_offset+1) = self%mrv
    !//            send_offset = send_offset+1            
	!//        	Endif
	!//			pcurrent = pcurrent+1
    !//        inext = self%istop(pcurrent)
	!//			if (i .lt. lm_count) send_offset = self%sdisp21(pcurrent)
	!//		Endif
    !//  Enddo

    !///////////////////////////////////////
    ! New way (for new layout)
        r_min = pfi%my_1p%min
        r_max = pfi%my_1p%max
        pcurrent = 0
        inext = num_lm(0)
        send_offset = 0
        Do i = 1, lm_count
            mp = mp_lm_values(i)
            l = l_lm_values(i)
            offset = 0 

            Do n = 1, nfields
              Do imi = 1, 2
                Do r = rmin,rmax
                  send_buff(send_offset+r) = self%s2b(mp)%data(l,r,imi,n)
                End Do
                send_offset = send_offset+delta_r
              Enddo
            Enddo
            If (i .eq. inext) Then
              If (self%do_piggy) Then
                ! load the piggyback value into the next buffer slot.
                send_buff(send_offset+1) = self%mrv
                send_offset = send_offset+1            
              Endif
              pcurrent = pcurrent+1
              inext = self%istop(pcurrent)
              If (i .lt. lm_count) send_offset = self%sdisp21(pcurrent)
            Endif
        Enddo

    !///////////////////////////////////////


		Call self%deconstruct('s2b')


		Allocate(recv_buff(1:self%recv_size21))


		Call Standard_Transpose(send_buff, recv_buff, self%scount21, self%sdisp21, self%rcount21, &
			self%rdisp21, pfi%ccomm, self%pad_buffer)
		DeAllocate(send_buff)

		Call self%construct('p1b')
		! Now, the receive striping needs a little mapping
		! WPS are coupled, but Z,Btor, and Bpol are not
		! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
		! We may want to modify this later on to mesh with linear equation structure
      
		
        Do p = 0, np - 1
            indx = self%rdisp21(p)+1

            r_min = pfi%all_1p(p)%min
            r_max = pfi%all_1p(p)%max

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
            if (self%do_piggy) then
                self%mrv = max(self%mrv,recv_buff(indx))
                indx = indx+1
            endif
        End Do

		self%config='p1b'
      Deallocate(recv_buff)
		
    End Subroutine Transpose_2b1b

    Subroutine Transpose_1a2a(self)
      ! Go from implicit configuration (1 physical) to configuration 2 (spectral)
      Implicit None

      Integer :: r,l, mp, lp, indx, r_min, r_max, dr, cnt,i
		Integer :: n, nfields, offset, delta_r, rmin, rmax, np,p
		Integer :: recv_offset, tnr

		Real*8, Allocatable :: send_buff(:),recv_buff(:)
      Integer :: inext, pcurrent

		Class(SphericalBuffer) :: self


		! This goes at the beginning
		
		nfields = self%nf1a

		Allocate(send_buff(1:self%send_size12))


		! Now, the receive striping needs a little mapping
		! WPS are coupled, but Z,Btor, and Bpol are not
		! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
		! We may want to modify this later on to mesh with linear equation structure

		np = pfi%ccomm%np
      Do p = 0, np - 1
			indx = self%sdisp12(p)+1
			r_min = pfi%all_1p(p)%min
			r_max = pfi%all_1p(p)%max
         dr = r_max - r_min
         ! Each processor in the radial group will have given me the same number of 
         ! l-m combos in the same (correct) order
         cnt = 1
         Do lp = 1, my_num_lm
				Do n = 1, nfields
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
		! We have nfields (nf1a)  and the data is dimensioned x%(l,r,field_num)

      ! Load the send array in the ordering of the l-m values
      indx = 1

		delta_r = pfi%my_1p%delta
		rmin = pfi%my_1p%min
		rmax = pfi%my_1p%max

		Allocate(recv_buff(1:self%recv_size12))


		Call Standard_Transpose(send_buff, recv_buff, self%scount12, self%sdisp12, &
			 self%rcount12, self%rdisp12, pfi%ccomm, self%pad_buffer)
		DeAllocate(send_buff)

		Call self%construct('s2a')

		rmin = 1
		rmax = delta_r*2
		recv_offset = 0
		tnr = 2*delta_r
		pcurrent = 0
		inext = num_lm(0)
		recv_offset = 0
        !Old way
        !//Do i = 1, lm_count
            !//mp = mp_lm_values(i)
            !//l = l_lm_values(i)

            !//offset = 0
            !//Do n = 1, nfields
            !//Do r = rmin,rmax 

            !//self%s2a(mp)%data(l,offset+r) = recv_buff(recv_offset+r)
            !//End Do
            !//offset = offset+tnr
            !//recv_offset = recv_offset+tnr
            !//Enddo
            !//If (i .eq. inext) Then
            !//			pcurrent = pcurrent+1
            !//      inext = self%istop(pcurrent)
            !//			if (i .lt. lm_count) recv_offset = self%rdisp12(pcurrent)
            !//		Endif
        !//Enddo
        !///////////////////////////////
      ! New way for new data_layout
        r_min = pfi%my_1p%min
        r_max = pfi%my_1p%max
        Do i = 1, lm_count
            mp = mp_lm_values(i)
            l = l_lm_values(i)


            Do n = 1, nfields
              Do imi = 1, 2
                Do r = rmin,rmax 
	              self%s2a(mp)%data(l,r,imi,n) = recv_buff(recv_offset+r)
                EndDo
                recv_offset = recv_offset+delta_r
              Enddo
            Enddo
            If (i .eq. inext) Then
	            pcurrent = pcurrent+1
                inext = self%istop(pcurrent)
	            if (i .lt. lm_count) recv_offset = self%rdisp12(pcurrent)
            Endif
        Enddo


		!////
		!Write(6,*)'maxval recv : ', maxval(abs(recv_buff))
		self%config='s2a'
      Deallocate(recv_buff)
		
    End Subroutine Transpose_1a2a


	Subroutine Allocate_Spherical_Buffer(self,config)
		Class(SphericalBuffer) :: self
		Character*3, Intent(In) :: config
		Integer :: mn1, mn2, mn3, mn4
		Integer :: mx1,mx2,mx3,mx4
		Integer :: i
		Select Case(config)
			Case('p1a')
				If (.not. Allocated(self%p1a)) Then
					mn1 = 1
					mx1 = pfi%n1p
					mn2 = 1
					mx2 = 2
					mn3 = 1
					mx3 = my_num_lm
					mn4 = 1
					mx4 = self%nf1a
					Allocate(self%p1a(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
				Endif
			Case('p1b')
				If (.not. Allocated(self%p1b)) Then
					mn1 = 1
					mx1 = pfi%n1p
					mn2 = 1
					mx2 = 2
					mn3 = 1
					mx3 = my_num_lm
					mn4 = 1
					mx4 = self%nf1b
					Allocate(self%p1b(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
				Endif

			Case('p2a')
				If (.not. Allocated(self%p2a)) Then				 
					mn1 = 1
					mx1 = pfi%n2p
					mn2 = pfi%my_1p%min
					mx2 = pfi%my_1p%max
					mn4 = pfi%my_3s%min
					mx4 = pfi%my_3s%max

					mx3 = self%nf2a*2*(mx2-mn2+1)
					Allocate(self%p2a(mn1:mx1, 1:mx3, mn4:mx4))
				Endif				
			Case('p2b')
				If (.not. Allocated(self%p2b)) Then 
					mn1 = 1
					mx1 = pfi%n2p
					mn2 = pfi%my_1p%min
					mx2 = pfi%my_1p%max
					mn4 = pfi%my_3s%min
					mx4 = pfi%my_3s%max

					mx3 = self%nf2b*2*(mx2-mn2+1)
					Allocate(self%p2b(mn1:mx1, 1:mx3, mn4:mx4))		
				Endif
			Case('p3a')
				If (.not. Allocated(self%p3a)) Then
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
				Endif
			Case('p3b')
				If (.not. Allocated(self%p3b)) Then
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
				Endif

			Case('s2a')
				If (.not. Allocated(self%s2a)) Then
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



					Do i = mn1, mx1
						mn2 = pfi%inds_3s(i)		!l_min = m
                        Allocate(self%s2a(i)%data(mn2:mx2,mn3:mx3,1:2,1:mx4)
					Enddo
				Endif
			Case('s2b')
				If (.not. Allocated(self%s2b)) Then
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

					Do i = mn1, mx1
						mn2 = pfi%inds_3s(i)		!l_min = m
						Allocate(self%s2b(i)%data(mn2:mx2,mn3:mx3,1:2,1:mx4))
					Enddo

				Endif
		End Select

	End Subroutine Allocate_Spherical_Buffer
