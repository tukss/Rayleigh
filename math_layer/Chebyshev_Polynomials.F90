Module Chebyshev_Polynomials
	! Module for computing Chebyshev Polynomial Arrays and the associated Derivative Arrays
	Implicit None
	Integer, Private :: N_max, N_even, N_odd, n_x
	Real*8, Allocatable, Private :: x(:)  ! The colocation points
	Real*8, Allocatable :: cheby(:,:)		! cheby(r,k) is chebyshev polynomial of degree k-1 at radius r
	Real*8, Allocatable :: cheby_even(:,:), cheby_odd(:,:) ! even and odd chebyshev arrays

	Real*8, Allocatable :: dcheby(:,:,:)	! The Chebyshev Derivative Arrays
	Real*8, private ::	Pi  = 3.1415926535897932384626433832795028841972d+0
	Real*8, private :: pi_over_N
	Logical :: DeAlias = .true.
	Logical :: Parity = .true.
	Logical :: initialized = .false.
	Real*8, Private :: scaling ! x runs from -0.5 to 0.5 by default
	Interface Cheby_To_Spectral
		Module Procedure To_Spectral_1D, To_Spectral_2D, To_Spectral_3D, To_Spectral_4D
	End Interface

	Interface Cheby_From_Spectral
		Module Procedure From_Spectral_1D, From_Spectral_2D, From_Spectral_3D, From_Spectral_4D
	End Interface

	Interface d_by_dr_cp
		Module Procedure Cheby_Deriv_Buffer_4D
	End Interface

Contains

	Subroutine Initialize_Chebyshev(grid, xmin,xmax, integration_weights)
		Implicit None
		Real*8, Intent(InOut) :: grid(:)
		Real*8, Intent(In), Optional :: xmin, xmax
		Real*8 ::delta, gmin, tmp, xx
        Real*8, Intent(InOut) :: integration_weights(1:)
        Integer :: r
		If (.not. initialized) Then
		    N_max = size(grid)
		    Call gen_colocation_points()
		    grid(:) = x(:)
		    Call Gen_Tn()
		    Call Gen_Tn_Deriv_Arrays(3)
		    If (present(xmin)) Then
			    delta = xmax-xmin
			    scaling = (x(1)-x(N_max))/delta
			    dcheby(:,:,1) = dcheby(:,:,1)*scaling
			    dcheby(:,:,2) = dcheby(:,:,2)*(scaling**2)
			    dcheby(:,:,3) = dcheby(:,:,3)*(scaling**3)
			    grid(:) = grid(:)/scaling
			    gmin = grid(N_max)
			    grid(:) = grid(:)-gmin+xmin
		    Endif
		    initialized = .true.
		Else
			grid(:) = x(:)
		Endif


        integration_weights(1:n_max) = 0.0d0


		tmp = 1.5d0*Pi * (grid(1)-grid(N_max)) / &
			& ( (grid(1)**3 - grid(N_max)**3) * N_max )
		Do r=1,N_max
			xx = (2.0d0*grid(r)-grid(N_max)-grid(1))/(grid(1)-grid(N_max))
			integration_weights(r) = grid(r)**2 * tmp * sqrt(1.0d0-xx*xx)
		Enddo


	End Subroutine Initialize_Chebyshev

	Subroutine Rescale_Grid_CP(length_scale)
		Implicit None
		Real*8, Intent(In) :: length_scale
		! Following initialization, we can rescale the chebyshev arrays if we choose
		! This is useful when nondimensionalizing after the reference state has been set up
		! (which typically requires a radial grid to have been established)
		
		dcheby(:,:,1) = dcheby(:,:,1)*length_scale
		dcheby(:,:,2) = dcheby(:,:,2)*(length_scale**2)
		dcheby(:,:,3) = dcheby(:,:,3)*(length_scale**3)

	End Subroutine Rescale_Grid_CP

	Subroutine Gen_Colocation_Points()
		Implicit None
		Integer :: i
		Real*8 :: arg
		Allocate(x(1:N_max))
		pi_over_N = pi/(N_max*1.0d0)
		arg = (0.5d0)*pi_over_N
		Do i = 1, N_Max
			x(i) = cos(arg)
			arg = arg+pi_over_N
		Enddo
	End Subroutine Gen_Colocation_Points

	Subroutine Gen_Tn()
		Implicit None
		Integer :: i, k, r
		Real*8 :: acx, arg
		Allocate(cheby(1:N_max,1:N_max))
		Do r = 1, N_max
			acx = pi_over_n*(r-1+0.5d0)
			Do i = 1, N_max
				k = i -1
				arg = k*acx
				cheby(r,i) = cos(arg)
			Enddo
		Enddo
		If (parity) Then
			n_odd = N_max/2
			n_even = n_odd+mod(N_max,2)
			n_x = n_even
			Allocate(cheby_even(1:N_x,1:N_even))
			Allocate(cheby_odd(1:N_x,1:N_odd))
			Do i = 1, n_even
				cheby_even(1:N_x,i) = cheby(1:N_x,2*i-1)
			Enddo
			Do i = 1,n_odd
				cheby_odd(1:N_x,i) = cheby(1:N_x,2*i)
			Enddo
			If (n_x .ne. n_odd) Then
				! We actually have an x = 0 point
				! We must be careful not to double count the power here
				! when exploiting parity to speed up the transforms.
				! No such adjustment need be made for the regular cheby array
				cheby_even(n_x,:) = 0.5d0*cheby_even(n_x,:)
				cheby_odd(n_x,:) = 0.5d0*cheby_odd(n_x,:)
				! There is really no need to modify cheby_odd, but this way it is consistently stored.
			Endif
		Endif
	End Subroutine Gen_Tn

	Subroutine Gen_Tn_Deriv_Arrays(dmax)
		Implicit None
		Integer, Intent(In) :: dmax
		Integer :: i, k,n,d 
		Real*8, Allocatable :: alpha(:,:)

		! sum_n (alpha_kn  c_n) = c'_k
		Allocate(alpha(1:N_max,1:N_max))
        alpha(:,:) = 0.0d0
		alpha(N_max,:) = 0.0d0
		alpha(N_max-1,N_max) = 2.0d0*(N_max-1)
		Do k = N_max-2, 1, -1
			alpha(k,k+1) = 2.0d0*k
			alpha(k,:) = alpha(k,:)+alpha(k+2,:)
		Enddo	 
		Allocate(dcheby(1:N_max,1:N_max,0:dmax))
		dcheby(:,:,:) = 0.0d0
		dcheby(:,:,0) = cheby(:,:)
		dcheby(:,1,0) = dcheby(:,1,0) - 0.5d0	! This accounts for the -1/2c_0 in going from_spectral
		If (dmax .ge. 1) Then
			Do d = 1, dmax
			Do n = 1, N_Max
			Do k = 1, N_max
				Do i = 1, N_max
					dcheby(k,n,d) = dcheby(k,n,d) + dcheby(k,i,d-1)*alpha(i,n)
				Enddo
			Enddo
			Enddo
			Enddo
		Endif
		DeAllocate(alpha)
	End Subroutine Gen_Tn_Deriv_Arrays

	Subroutine To_Spectral_1D(f_in,c_out)
		Implicit None
		Real*8, Intent(In) :: f_in(:)
		Real*8, Intent(InOut) :: c_out(:)
		Real*8 :: alpha, beta
		Real*8, Allocatable :: f_even(:), f_odd(:), c_temp(:)
		Integer :: i
		alpha = 2.0d0/n_max
		beta = 0.0d0
		
		If (parity) Then
			Allocate(c_temp(1:n_even ))
			Allocate(f_even(1:n_x ))
			Allocate( f_odd(1:n_x  ))
			Do i = 1, N_x
				f_even(i) = f_in(i)+f_in(N_max-i+1)
				f_odd(i)  = f_in(i)-f_in(N_max-i+1)
			Enddo
			CALL DGEMM('T','N',N_even,1,N_x, alpha, cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
			Do i = 1, n_even
				c_out(2*i-1) = c_temp(i)
			Enddo
			If (n_even .ne. n_odd) Then
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd))
			Endif
			Call DGEMM('T','N',N_odd,1,N_x, alpha, cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
			Do i = 1, N_odd
				c_out(2*i) = c_temp(i)
			Enddo
			DeAllocate(f_even,f_odd,c_temp)
			
		Else
			CALL DGEMM('T','N',N_max,1,N_Max, alpha, cheby, N_max,f_in , N_Max, beta,c_out,N_max)
		Endif
	End Subroutine To_Spectral_1D

	Subroutine To_Spectral_2D(f_in,c_out)
		Implicit None
		Real*8, Intent(In) :: f_in(:,:)
		Real*8, Intent(InOut) :: c_out(:,:)
		Real*8 :: alpha, beta
		Real*8, Allocatable :: f_even(:,:), f_odd(:,:), c_temp(:,:)
		Integer :: i, j, n2, dims(2)
		alpha = 2.0d0/n_max
		beta = 0.0d0
		dims = shape(f_in)
		n2 = dims(2)
		If (parity) Then
			Allocate(c_temp(1:n_even, 1:n2 ))
			Allocate(f_even(1:n_x   , 1:n2 ))
			Allocate( f_odd(1:n_x   , 1:n2 ))
			Do j = 1, n2
			Do i = 1, N_x
				f_even(i,j) = f_in(i,j)+f_in(N_max-i+1,j)
				f_odd(i,j)  = f_in(i,j)-f_in(N_max-i+1,j)
			Enddo
			Enddo
			CALL DGEMM('T','N',N_even,n2,N_x, alpha, cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
			Do j = 1, n2
			Do i = 1, n_even
				c_out(2*i-1,j) = c_temp(i,j)
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2))
			Endif
			Call DGEMM('T','N',N_odd,n2,N_x, alpha, cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
			Do j = 1, n2
			Do i = 1, N_odd
				c_out(2*i,j) = c_temp(i,j)
			Enddo
			Enddo
			DeAllocate(f_even,f_odd,c_temp)
			
		Else
			CALL DGEMM('T','N',N_max,n2,N_Max, alpha, cheby, N_max,f_in , N_Max, beta,c_out,N_max)
		Endif
	End Subroutine To_Spectral_2D

	Subroutine To_Spectral_3D(f_in,c_out)
		Implicit None
		! Play a sneaky trick on dgemm and see if it sticks without complaining
		Real*8, Intent(In) :: f_in(:,:,:)
		Real*8, Intent(InOut) :: c_out(:,:,:)
		Real*8 :: alpha, beta
		Real*8, Allocatable :: f_even(:,:,:), f_odd(:,:,:), c_temp(:,:,:)
		Integer :: i, j, k, n2, n3, dims(3)
		alpha = 2.0d0/n_max
		beta = 0.0d0
		dims = shape(f_in)
		n2 = dims(2)
		n3 = dims(3)
		If (parity) Then
			Allocate(c_temp(1:n_even, 1:n2, 1:n3 ))
			Allocate(f_even(1:n_x   , 1:n2, 1:n3 ))
			Allocate( f_odd(1:n_x   , 1:n2, 1:n3 ))
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, N_x
				f_even(i,j,k) = f_in(i,j,k)+f_in(N_max-i+1,j,k)
				f_odd(i,j,k)  = f_in(i,j,k)-f_in(N_max-i+1,j,k)
			Enddo
			Enddo
			Enddo
			CALL DGEMM('T','N',N_even,n2*n3,N_x, alpha, cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				c_out(2*i-1,j,k) = c_temp(i,j,k)
			Enddo
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2,1:n3))
			Endif
			Call DGEMM('T','N',N_odd,n2*n3,N_x, alpha, cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, N_odd
				c_out(2*i,j,k) = c_temp(i,j,k)
			Enddo
			Enddo
			Enddo
			DeAllocate(f_even,f_odd,c_temp)
			
		Else
			CALL DGEMM('T','N',N_max,n2*n3,N_Max, alpha, cheby, N_max,f_in , N_Max, beta,c_out,N_max)
		Endif
	End Subroutine To_Spectral_3D

	Subroutine To_Spectral_4D(f_in,c_out)
		Implicit None
		! Play a sneaky trick on dgemm and see if it sticks without complaining
		Real*8, Intent(In) :: f_in(:,:,:,:)
		Real*8, Intent(InOut) :: c_out(:,:,:,:)
		Real*8 :: alpha, beta
		Real*8, Allocatable :: f_even(:,:,:,:), f_odd(:,:,:,:), c_temp(:,:,:,:)
		Integer :: i, j, k, kk, n2, n3, n4, dims(4)
		alpha = 2.0d0/n_max
		beta = 0.0d0
		dims = shape(f_in)
		n2 = dims(2)
		n3 = dims(3)
		n4 = dims(4)
		If (parity) Then
			Allocate(c_temp(1:n_even, 1:n2, 1:n3, 1:n4 ))
			Allocate(f_even(1:n_x   , 1:n2, 1:n3, 1:n4 ))
			Allocate( f_odd(1:n_x   , 1:n2, 1:n3, 1:n4 ))
			Do kk= 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, N_x
				f_even(i,j,k, kk) = f_in(i,j,k,kk)+f_in(N_max-i+1,j,k,kk)
				f_odd(i,j,k, kk)  = f_in(i,j,k,kk)-f_in(N_max-i+1,j,k,kk)
			Enddo
			Enddo
			Enddo
			Enddo
			CALL DGEMM('T','N',N_even,n2*n3*n4,N_x, alpha, cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
			Do kk = 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				c_out(2*i-1,j,k,kk) = c_temp(i,j,k,kk)
			Enddo
			Enddo
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2,1:n3,1:n4))
			Endif
			Call DGEMM('T','N',N_odd,n2*n3*n4,N_x, alpha, cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
			Do kk =1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, N_odd
				c_out(2*i,j,k,kk) = c_temp(i,j,k,kk)
			Enddo
			Enddo
			Enddo
			Enddo
			DeAllocate(f_even,f_odd,c_temp)
			
		Else
			CALL DGEMM('T','N',N_max,n2*n3*n4,N_Max, alpha, cheby, N_max,f_in , N_Max, beta,c_out,N_max)
		Endif
	End Subroutine To_Spectral_4D

	Subroutine From_Spectral_1D(c_in,f_out)
		Implicit None
		Real*8, Intent(In) :: c_in(:)
		Real*8, Intent(InOut) :: f_out(:)
		Real*8, Allocatable :: c_temp(:), f_temp(:)
		Real*8 :: alpha, beta
		Integer :: i
		alpha = 1.0d0
		beta = 0.0d0



		If (parity) Then
			Allocate(c_temp(1:n_even))
			Allocate(f_temp(1:n_x))
			Do i = 1, n_even
				c_temp(i) = c_in(2*i-1) 
			Enddo
			CALL DGEMM('N','N',N_x,1,N_even, alpha, cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
			Do i = 1, n_even
				f_out(i) = f_temp(i)
				f_out(N_max-i+1) = f_temp(i)
			Enddo
			If (n_even .ne. n_odd) Then
				f_out(n_x) = f_out(n_x)*2.0d0
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd))
			Endif
			Do i = 1, n_odd
				c_temp(i) = c_in(2*i) 
			Enddo
			CALL DGEMM('N','N',N_x,1,N_odd, alpha, cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)
			Do i = 1, n_odd
				f_out(i) = f_out(i) + f_temp(i)
				f_out(N_max-i+1) = f_out(N_max-i+1)-f_temp(i)
			Enddo
			If (n_even .ne. n_odd) Then
				f_out(n_x) = f_out(n_x) + f_temp(n_x)*2.0d0
			Endif
			DeAllocate(c_temp, f_temp)
		Else
			CALL DGEMM('N','N',N_max,1,N_Max, alpha, cheby, N_max,c_in , N_Max, beta,f_out,N_max)

		Endif
		f_out(:) = f_out(:) -c_in(1)/2.0d0
	End Subroutine From_Spectral_1D

	Subroutine From_Spectral_2D(c_in,f_out)
		Implicit None
		Real*8, Intent(In) :: c_in(:,:)
		Real*8, Intent(InOut) :: f_out(:,:)
		Real*8, Allocatable :: c_temp(:,:), f_temp(:,:)
		Real*8 :: alpha, beta
		Integer :: i, j, n2, dims(2)
		alpha = 1.0d0
		beta = 0.0d0

		dims = shape(c_in)
		n2 = dims(2)

		If (parity) Then
			Allocate(c_temp(1:n_even,1:n2))
			Allocate(f_temp(1:n_x,1:n2))
			Do j = 1, n2
			Do i = 1, n_even
				c_temp(i,j) = c_in(2*i-1,j) 
			Enddo
			enddo

			CALL DGEMM('N','N',N_x,n2,N_even, alpha, cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
			Do j = 1, n2
			Do i = 1, n_even
				f_out(i,j) = f_temp(i,j)
				f_out(N_max-i+1,j) = f_temp(i,j)
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				Do j = 1, n2
					f_out(n_x,j) = f_out(n_x,j)*2.0d0
				Enddo
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2))
			Endif
			Do j = 1, n2
			Do i = 1, n_odd
				c_temp(i,j) = c_in(2*i,j) 
			Enddo
			Enddo
			CALL DGEMM('N','N',N_x,n2,N_odd, alpha, cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

			Do j = 1, n2
			Do i = 1, n_odd
				f_out(i,j) = f_out(i,j) + f_temp(i,j)
				f_out(N_max-i+1,j) = f_out(N_max-i+1,j)-f_temp(i,j)
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				f_out(n_x,:) = f_out(n_x,:) + f_temp(n_x,:)*2.0d0
			Endif
			DeAllocate(c_temp, f_temp)
		Else
			CALL DGEMM('N','N',N_max,n2,N_Max, alpha, cheby, N_max,c_in , N_Max, beta,f_out,N_max)

		Endif
		Do j = 1, n2
			f_out(:,j) = f_out(:,j) -c_in(1,j)/2.0d0
		Enddo
	End Subroutine From_Spectral_2D

	Subroutine From_Spectral_3D(c_in,f_out)
		Implicit None
		Real*8, Intent(In) :: c_in(:,:,:)
		Real*8, Intent(InOut) :: f_out(:,:,:)
		Real*8, Allocatable :: c_temp(:,:,:), f_temp(:,:,:)
		Real*8 :: alpha, beta
		Integer :: i, j, k, n2, n3, dims(3)
		alpha = 1.0d0
		beta = 0.0d0

		dims = shape(c_in)
		n2 = dims(2)
		n3 = dims(3)
		If (parity) Then
			Allocate(c_temp(1:n_even,1:n2, 1:n3))
			Allocate(f_temp(1:n_x   ,1:n2, 1:n3))
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				c_temp(i,j,k) = c_in(2*i-1,j,k) 
			Enddo
			Enddo
			Enddo

			CALL DGEMM('N','N',N_x,n2*n3,N_even, alpha, cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				f_out(i,j,k) = f_temp(i,j,k)
				f_out(N_max-i+1,j,k) = f_temp(i,j,k)
			Enddo
			Enddo
			Enddo

			If (n_even .ne. n_odd) Then
				Do k = 1, n3
				Do j = 1, n2
					f_out(n_x,j,k) = f_out(n_x,j,k)*2.0d0
				Enddo
				Enddo
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2,1:n3))
			Endif
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_odd
				c_temp(i,j,k) = c_in(2*i,j,k) 
			Enddo
			Enddo
			Enddo
			CALL DGEMM('N','N',N_x,n2*n3,N_odd, alpha, cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_odd
				f_out(i,j,k) = f_out(i,j,k) + f_temp(i,j,k)
				f_out(N_max-i+1,j,k) = f_out(N_max-i+1,j,k)-f_temp(i,j,k)
			Enddo
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				f_out(n_x,:,:) = f_out(n_x,:,:) + f_temp(n_x,:,:)*2.0d0
			Endif
			DeAllocate(c_temp, f_temp)
		Else
			CALL DGEMM('N','N',N_max,n2*n3,N_Max, alpha, cheby, N_max,c_in , N_Max, beta,f_out,N_max)

		Endif
		Do k = 1, n3
		Do j = 1, n2
			f_out(:,j,k) = f_out(:,j,k) -c_in(1,j,k)/2.0d0
		Enddo
		Enddo
	End Subroutine From_Spectral_3D

	Subroutine From_Spectral_4D(c_in,f_out)
		Implicit None
		Real*8, Intent(In) :: c_in(:,:,:,:)
		Real*8, Intent(InOut) :: f_out(:,:,:,:)
		Real*8, Allocatable :: c_temp(:,:,:,:), f_temp(:,:,:,:)
		Real*8 :: alpha, beta
		Integer :: i, j, k, kk, n2, n3, n4, dims(4)
		alpha = 1.0d0
		beta = 0.0d0

		dims = shape(c_in)
		n2 = dims(2)
		n3 = dims(3)
		n4 = dims(4)
		If (parity) Then
			Allocate(c_temp(1:n_even,1:n2, 1:n3, 1:n4))
			Allocate(f_temp(1:n_x   ,1:n2, 1:n3, 1:n4))
			Do kk = 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				c_temp(i,j,k,kk) = c_in(2*i-1,j,k,kk) 
			Enddo
			Enddo
			Enddo
			Enddo

			CALL DGEMM('N','N',N_x,n2*n3*n4,N_even, alpha, cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)

			Do kk = 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_even
				f_out(i,j,k,kk) = f_temp(i,j,k,kk)
				f_out(N_max-i+1,j,k,kk) = f_temp(i,j,k,kk)
			Enddo
			Enddo
			Enddo
			Enddo

			If (n_even .ne. n_odd) Then
				Do kk = 1, n4
				Do k = 1, n3
				Do j = 1, n2
					f_out(n_x,j,k,kk) = f_out(n_x,j,k,kk)*2.0d0
				Enddo
				Enddo
				Enddo
				DeAllocate(c_temp)
				Allocate(c_temp(1:n_odd,1:n2,1:n3,1:n4))
			Endif
			Do kk = 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_odd
				c_temp(i,j,k,kk) = c_in(2*i,j,k,kk) 
			Enddo
			Enddo
			Enddo
			Enddo
			CALL DGEMM('N','N',N_x,n2*n3*n4,N_odd, alpha, cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

			Do kk = 1, n4
			Do k = 1, n3
			Do j = 1, n2
			Do i = 1, n_odd
				f_out(i,j,k,kk) = f_out(i,j,k,kk) + f_temp(i,j,k,kk)
				f_out(N_max-i+1,j,k,kk) = f_out(N_max-i+1,j,k,kk)-f_temp(i,j,k,kk)
			Enddo
			Enddo
			Enddo
			Enddo
			If (n_even .ne. n_odd) Then
				f_out(n_x,:,:,:) = f_out(n_x,:,:,:) + f_temp(n_x,:,:,:)*2.0d0
			Endif
			DeAllocate(c_temp, f_temp)
		Else
			CALL DGEMM('N','N',N_max,n2*n3*n4,N_Max, alpha, cheby, N_max,c_in , N_Max, beta,f_out,N_max)

		Endif
		Do kk = 1, n4
		Do k = 1, n3
		Do j = 1, n2
			f_out(:,j,k,kk) = f_out(:,j,k,kk) -c_in(1,j,k,kk)/2.0d0
		Enddo
		Enddo
		Enddo
	End Subroutine From_Spectral_4D

	Subroutine Cheby_Deriv_Buffer_4D(ind,dind,buffer,dorder)
		Implicit None
		Real*8,  Intent(InOut) :: buffer(0:,1:,1:,1:)	! Makes it easier to reconcile with my IDL code
		Integer, Intent(In)    :: ind, dind, dorder
		Real*8, Allocatable :: dbuffer(:,:)
		Integer :: dims(4), n,n2,n3, i,j,k, order
		dims = shape(buffer)
		n = dims(1)
		n2 = dims(2)
		n3 = dims(3)
		If (ind .ne. dind) Then
			Do k = 1, n3
				Do j = 1, n2
					buffer(n-1,j,k,dind) = 0.0d0
					buffer(n-2,j,k,dind) = 2.0d0*(n-1)*buffer(n-1,j,k,ind)*scaling
					Do i = n-3,0, -1
						buffer(i,j,k,dind) = buffer(i+2,j,k,dind)+2.0d0*(i+1)*buffer(i+1,j,k,ind)*scaling
					Enddo
				Enddo
			Enddo
			If (dorder .gt. 1) Then 
				Allocate(dbuffer(0:n-1,1:dorder))
				Do k = 1, n3
					Do j = 1, n2
						dbuffer(:,1) = buffer(:,j,k,dind)
						Do order = 2, dorder
							dbuffer(n-1,order) = 0.0d0
							dbuffer(n-2,order) = 2.0d0*(n-1)*dbuffer(n-1,order-1)*scaling
							Do i = n -3, 0, -1
								dbuffer(i,order) = dbuffer(i+2,order)+2.0d0*(i+1)*dbuffer(i+1,order-1)*scaling						
							Enddo
						Enddo
						buffer(:,j,k,dind) = dbuffer(:,dorder)
					Enddo
				Enddo

				DeAllocate(dbuffer)
			Endif
		Else
			! In-place
		Endif
		buffer(:,:,:,dind) = buffer(:,:,:,dind)
	End Subroutine Cheby_Deriv_Buffer_4D	


	!////////////////////////////////////////////////////////////////////
	! Matrix Row-Loading routines for use with the Implicit Solve
	!=============================
	! Implicit.F row-loading routines  -- this should maybe be moved out later, but OK for now
	!Subroutine Load_Interior_Rows(row,col,amp,dorder,mpointer)
	!	Integer, Intent(In) :: row, col, dorder ! ,rr,rp,n,np
	!	real*8, Intent(In) :: amp(:)
	!	real*8, Pointer, Dimension(:,:), Intent(In) :: mpointer
	!	
	!			Call Load_Interior_RowsFD(row,col,amp,dorder,mpointer)
	!
	!
	!	
	!End Subroutine Load_Interior_Rows
	!Subroutine Load_Single_Row(r,row,col,amp,dorder,mpointer, clear_row, amp_arr, boundary)
	!	Integer, Intent(In) :: r,row, col, dorder ! ,rr,rp,n,np
	!	real*8, Intent(In) :: amp
	!	real*8, Intent(In), Optional :: amp_arr(:)	! lets us load every column of the row
	!	real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer
	!	Logical, Intent(In), Optional :: clear_row, boundary

	
	!			Call Load_Single_RowFD(r,row,col,amp,dorder,mpointer, clear_row, amp_arr, boundary)


		
	!End Subroutine Load_Single_Row
	Subroutine Load_Interior_Rows_Cheby(row,col,amp,dorder,mpointer)
		Integer, Intent(In) :: row, col, dorder 
		Integer :: r, n
		real*8, Intent(In) :: amp(:)
		real*8, Pointer, Dimension(:,:), Intent(In) :: mpointer
		Do r = 1, N_max
			Do n = 1, N_max
				mpointer(row+r,col+n) = mpointer(row+r,col+n)+amp(r)*dcheby(r,n,dorder)
			Enddo
		Enddo
	End Subroutine Load_Interior_Rows_Cheby


	Subroutine Load_Single_Row_Cheby(r,row,col,amp,dorder,mpointer, clear_row, boundary)
		Integer, Intent(In) :: r,row, col, dorder
		Integer :: n
		real*8, Intent(In) :: amp
		real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer
		Logical, Intent(In), Optional :: clear_row, boundary
		Logical :: bjunk
		If (present(clear_row)) Then
			! clear everything in this row
			mpointer(r+row,:) = 0.0d0
		Endif
		If (present(boundary)) Then
			! Do nothing at the moment
			bjunk = boundary	! This is in development, but placeholder to avoid Intel compiler warnings (unused vars)
		Endif
		Do n = 1, (2*N_max)/3 	! De-Alias at boundaries (single rows are really just for boundaries)
			mpointer(row+r,col+n) = mpointer(row+r,col+n)+amp*dcheby(r,n,dorder)
		Enddo

	End Subroutine Load_Single_Row_Cheby

End Module Chebyshev_Polynomials
