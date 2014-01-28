Module Output
	Implicit None
	Interface Write_Array
		Module Procedure Write_DArray_1D , Write_DArray_2D , Write_DArray_3D
		Module Procedure Write_DCArray_1D
	End Interface

Contains

	Subroutine Write_DArray_1D(array, filename, funit)
		Real*8, Intent(In) :: array(:)
		Character(len=*), Intent(In) :: filename
		Integer, Intent(In), Optional :: funit
		Integer :: ndim, file_unit
		Integer :: i, ni, ilb, iub

		If (present(funit)) Then
			file_unit = funit	
		else
			file_unit = 39
		Endif	
		ilb = lbound(array,1)
		iub = ubound(array,1)
		ni = iub-ilb+1

		ndim = 1		
	   open(file_unit,FILE=filename,form='unformatted', status ='replace')
		write(file_unit) ndim
		Write(file_unit) ni
		write(file_unit) (array(i),i=ilb,iub)
	   close(file_unit)
	End Subroutine Write_DArray_1D

	Subroutine Write_DCArray_1D(array, filename, funit)
		Complex*16, Intent(In) :: array(:)
		Character(len=*), Intent(In) :: filename
		Integer, Intent(In), Optional :: funit
		Integer :: ndim, file_unit
		Integer :: i, ni, ilb, iub

		If (present(funit)) Then
			file_unit = funit	
		else
			file_unit = 39
		Endif	
		ilb = lbound(array,1)
		iub = ubound(array,1)
		ni = iub-ilb+1

		ndim = 1		
	   open(file_unit,FILE=filename,form='unformatted', status ='replace')
		write(file_unit) ndim
		Write(file_unit) ni
		write(file_unit) (array(i),i=ilb,iub)
	   close(file_unit)
	End Subroutine Write_DCArray_1D

	Subroutine Write_DArray_2D(array, filename, funit)
		Real*8, Intent(In) :: array(:,:)
		Character(len=*), Intent(In) :: filename
		Integer, Intent(In), Optional :: funit
		Integer :: ndim, file_unit
		Integer :: i, ni, ilb, iub
		Integer :: j, nj, jlb, jub

		If (present(funit)) Then
			file_unit = funit	
		else
			file_unit = 39
		Endif	
		ilb = lbound(array,1)
		iub = ubound(array,1)
		ni = iub-ilb+1
		jlb = lbound(array,2)
		jub = ubound(array,2)
		nj = jub-jlb+1

		ndim = 2		
	   open(file_unit,FILE=filename,form='unformatted', status ='replace')
		write(file_unit) ndim
		Write(file_unit) ni
		Write(file_unit) nj
		write(file_unit) ((array(i,j),i=ilb,iub),j=jlb,jub)
	   close(file_unit)
	End Subroutine Write_DArray_2D

	Subroutine Write_DArray_3D(array, filename, funit)
		Real*8, Intent(In) :: array(:,:,:)
		Character(len=*), Intent(In) :: filename
		Integer, Intent(In), Optional :: funit
		Integer :: ndim, file_unit
		Integer :: i, ni, ilb, iub
		Integer :: j, nj, jlb, jub
		Integer :: k, nk, klb, kub
		If (present(funit)) Then
			file_unit = funit	
		else
			file_unit = 39
		Endif	
		ilb = lbound(array,1)
		iub = ubound(array,1)
		ni = iub-ilb+1
		jlb = lbound(array,2)
		jub = ubound(array,2)
		nj = jub-jlb+1
		klb = lbound(array,3)
		kub = ubound(array,3)
		nk = kub-klb+1

		ndim = 3		
	   open(file_unit,FILE=filename,form='unformatted', status ='replace')
		write(file_unit) ndim
		Write(file_unit) ni
		Write(file_unit) nj
		Write(file_unit) nk
		write(file_unit) (((array(i,j,k),i=ilb,iub),j=jlb,jub),k=klb,kub)
	   close(file_unit)
	End Subroutine Write_DArray_3D

End Module Output
