Module Grid
	Implicit None
	! ================================
	!  Grid dimensions
	Integer, Public :: nx1 = 1
	Integer, Public :: nx2 = 1
	Integer, Public :: nx3 = 1
	Real*8, Public  :: x1min = 0.0d0, x1max = 0.0d0
	Real*8, Public  :: x2min = 0.0d0, x2max = 0.0d0
	Real*8, Public  :: x3min = 0.0d0, x3max = 0.0d0
	Real*8, Allocatable :: x1(:), x2(:), x3(:)
	Real*8 :: dx1, dx2, dx3

	Logical :: x_periodic = .false.		! Force this direction to span 0 to 2Pi
	Logical :: y_periodic = .false.
	Logical :: z_periodic = .false.

	!=================================
	!	Renaming grid dimensions (for convenience)
	Integer, Pointer, Public :: nx, ny, nz, nr, ntheta,nphi
	Real*8,  Pointer, Public :: x_min, x_max, y_min, y_max, z_min, z_max
	Real*8,  Pointer, Public :: r_min, r_max, theta_min, theta_max, phi_min, phi_max
	Real*8,  Pointer, Dimension(:) :: xgrid(:), ygrid(:), zgrid(:), radius(:), theta(:), phi(:)


	Namelist /Grid_Namelist/ nx1, x1min, x1max,nx2,x2min,x2max,nx3,x3min,x3max, &
		& x_periodic, y_periodic, z_periodic

Contains
	Subroutine Initialize_Grid()
		Implicit None
		Integer :: i 
		Real*8 :: Pi

		Allocate(x1(1:nx1))
		Allocate(x2(1:nx2))
		Allocate(x3(1:nx3))

		If (x_periodic) Then 
			Pi = acos(-1.0d0)
			x1min = 0.0d0
			x1max = Pi*2.0d0
		Endif
		If (y_periodic) Then 
			Pi = acos(-1.0d0)
			x2min = 0.0d0
			x2max = Pi*2.0d0
		Endif
		If (z_periodic) Then 
			Pi = acos(-1.0d0)
			x3min = 0.0d0
			x3max = Pi*2.0d0
		Endif

		x1(1) = x1min
		x2(1) = x2min
		x3(1) = x3min

		If (nx1 .gt. 1) Then

				
			dx1 = (x1max-x1min)/real(nx1-1)



			Do i = 2, nx1
				x1(i) = x1(i-1)+dx1
			Enddo
		Endif

		If (nx2 .gt. 1) Then

		
			dx2 = (x2max-x2min)/real(nx2-1)



			Do i = 2, nx2
				x2(i) = x2(i-1)+dx2
			Enddo
		Endif

		If (nx3 .gt. 1) Then

		
			dx3 = (x3max-x3min)/real(nx3-1)

			Do i = 2, nx3
				x3(i) = x3(i-1)+dx3
			Enddo
		Endif

	End Subroutine Initialize_Grid

End Module Grid
