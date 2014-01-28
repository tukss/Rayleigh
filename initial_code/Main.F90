Program Main
	Use Input
	Use Grid
	Use Output
	Use Field_Manipulation
	Implicit None
	include 'fftw3.f'
	Character*120 :: filename
	Real*8, Allocatable :: arr1d(:), temp1(:)
	Real*8, Allocatable :: arr2d(:,:)
	Real*8, Allocatable :: arr3d(:,:,:)
	Real*8, Pointer, Dimension(:,:,:) :: temp
	Complex*16, Allocatable :: the_fft(:)
	Integer :: nd1,nd2,nd3,i,j,k
	Integer :: nkx	
	Integer*8 :: plan1
	Type(Field) :: field1
	Call Main_Input()	
	Call Initialize_Grid()
	
	Allocate(temp1(1:nx1))
	Allocate(arr1d(1:nx1))
	arr1d(:) = 0.0d0
	nkx = 5
	Allocate(the_fft(1:nx1/2+1))
	Do i = 0, nkx
		temp1(:) = sin(i*x1)		
		arr1d = arr1d+temp1
	Enddo

	filename = 'function'
	call Write_Array(arr1d,filename)

	call dfftw_plan_dft_r2c_1d(plan1,nx1,arr1d,the_fft,FFTW_ESTIMATE)
	call dfftw_execute_dft_r2c(plan1, arr1d, the_fft)
	call dfftw_destroy_plan(plan1)

	filename = 'fft'
	call Write_Array(the_fft,filename)

	nullify(field1%rdata)
	Allocate(field1%rdata(1:nx1,1:nx2,1:nx3))


	Call Field1%print_dims()
	!Call Field1%advance_config()
	Write(6,*)'Advancing configurations'
	Do i = 1, 9 
		Call Field1%advance_config_serial()
		Call Field1%print_dims()
	Enddo
End Program Main
