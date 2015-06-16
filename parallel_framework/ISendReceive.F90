Module ISendReceive
	Use MPI_BASE

	Interface ISend
		Module Procedure D_ISend_1D, D_ISend_2D, D_ISend_3D, D_ISend_4D
        Module Procedure D_ISend_5D
		Module Procedure Z_ISend_1D, Z_ISend_2D, Z_ISend_3D
	End Interface 

 	Interface IReceive
		Module Procedure D_IReceive_1D, D_IReceive_2D, D_IReceive_3D, D_IReceive_4D
        Module Procedure D_IReceive_5D
		Module Procedure Z_IReceive_1D, Z_IReceive_2D, Z_IReceive_3D
	End Interface 

  Integer :: mpi_err

Contains

	Subroutine IWait(irq)
		Implicit None
		Integer :: irq, status(MPI_STATUS_SIZE), mpi_err
		Call MPI_WAIT(irq,status,mpi_err)
	End Subroutine IWait

	Subroutine IWaitAll(n,irq)
		Integer :: irq(:)
		Integer, Intent(In) :: n
		Integer :: mpi_err
		Integer, Allocatable :: istat(:,:)
		Allocate(istat(MPI_STATUS_SIZE,1:n))
		Call MPI_WAITALL(n,irq,istat,mpi_err)
		DeAllocate(istat)
	End Subroutine IWaitAll

	Subroutine D_ISend_1D(x, irq,n_elements, istart,dest, tag, grp)
		Real(8), Intent(in)  :: x(:)

		Integer, Optional :: dest, n_elements, tag,istart
		Type(communicator), optional :: grp
		Integer :: p, n, comm2, tag2, irq,ione
	
		If (Present(n_elements)) Then
			n = n_elements
		Else
			n = Size(x)
		End If

		If (Present(dest)) Then
			p = dest
		Else
			p = 0
		End If

		If (Present(grp)) Then
			comm2 = grp%comm
		Else
			comm2 = MPI_COMM_WORLD
		End If

		If (Present(tag)) Then
			tag2 = tag
		Else
			tag2 = p
		End If

		If (Present(istart)) Then
			ione = istart
		Else
			ione = 1
		End if
    
		Call mpi_isend(x(ione), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)
    
	End Subroutine D_ISend_1D
	Subroutine D_ISend_4D(x, irq,n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(1:,1:,1:,1:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:4)
	 Integer :: istart, kstart, jstart,lstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2,irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
 	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
		lstart = indstart(4)
	Else
		istart = 1
		jstart = 1
		kstart = 1
		lstart = 1
	Endif   
    Call mpi_isend(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)
    !write(6,*)'zs ', p
	End Subroutine D_ISend_4D

	Subroutine D_IReceive_4D(x, irq,n_elements, source, tag, grp,indstart)
		Real*8, Intent(out)  :: x(1:,1:,1:,1:)

    Integer, Optional :: source, n_elements, tag,indstart(1:4)
	 Integer :: istart,jstart,kstart,lstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
		lstart = indstart(4)
	Else
		istart = 1
		jstart = 1
		kstart = 1
		lstart = 1
	Endif
    
    Call mpi_irecv(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
     !   write(6,*)'zr ', p

	End Subroutine D_IReceive_4D

	Subroutine D_ISend_5D(x, irq,n_elements, dest, tag, grp, indstart)
        Real*8, Intent(in)  :: x(1:,1:,1:,1:,1:)

        Integer, Optional :: dest, n_elements, tag,indstart(1:5)
	     Integer :: istart, kstart, jstart,lstart, mstart
        Type(communicator), optional :: grp
        Integer :: p, n, comm2, tag2,irq

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(dest)) Then
           p = dest
        Else
           p = 0
        End If

        If (Present(grp)) Then
           comm2 = grp%comm
        Else
           comm2 = MPI_COMM_WORLD
        End If

        If (Present(tag)) Then
           tag2 = tag
        Else
           tag2 = p
        End If
     	If (Present(indstart)) Then
		    istart = indstart(1)
		    jstart = indstart(2)
		    kstart = indstart(3)
		    lstart = indstart(4)
            mstart = indstart(5)
	    Else
		    istart = 1
		    jstart = 1
		    kstart = 1
		    lstart = 1
            mstart = 1
	    Endif   
        Call mpi_isend(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_5D

    Subroutine D_IReceive_5D(x, irq,n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(1:,1:,1:,1:,1:)

        Integer, Optional :: source, n_elements, tag,indstart(1:5)
	     Integer :: istart,jstart,kstart,lstart,mstart
        Type(communicator), optional :: grp
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), irq

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
        Else
           p = MPI_ANY_SOURCE
        End If

        If (Present(grp)) Then
           comm2 = grp%comm
        Else
           comm2 = MPI_COMM_WORLD
        End If

        If (Present(tag)) Then
           tag2 = tag
        Else
           tag2 = MPI_ANY_TAG
        End If

	    If (Present(indstart)) Then
		    istart = indstart(1)
		    jstart = indstart(2)
		    kstart = indstart(3)
		    lstart = indstart(4)
            mstart = indstart(5)
	    Else
		    istart = 1
		    jstart = 1
		    kstart = 1
		    lstart = 1
            mstart = 1
	    Endif
        
        Call mpi_irecv(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

        If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
        If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)


	End Subroutine D_IReceive_5D


	Subroutine Z_ISend_1D(x, irq,n_elements, istart,dest, tag, grp)
		Complex*16, Intent(in)  :: x(:)

		Integer, Optional :: dest, n_elements, tag,istart
		Type(communicator), optional :: grp
		Integer :: p, n, comm2, tag2, irq,ione
	
		If (Present(n_elements)) Then
			n = n_elements
		Else
			n = Size(x)
		End If

		If (Present(dest)) Then
			p = dest
		Else
			p = 0
		End If

		If (Present(grp)) Then
			comm2 = grp%comm
		Else
			comm2 = MPI_COMM_WORLD
		End If

		If (Present(tag)) Then
			tag2 = tag
		Else
			tag2 = p
		End If

		If (Present(istart)) Then
			ione = istart
		Else
			ione = 1
		End if
    
		Call mpi_isend(x(ione), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)
    
	End Subroutine Z_ISend_1D

Subroutine D_ISend_2D(x, irq,n_elements, dest, tag, grp,indstart)
    Real(8), Intent(in)  :: x(1:,1:)

    Integer, Optional :: dest, n_elements, tag, indstart(2)
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq
	 Integer :: ione, jone
    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If
	 
    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
    If (Present(indstart)) Then
			ione = indstart(1)
			jone = indstart(2)
	 Else
			ione = 1
			jone = 1
    Endif
    Call mpi_isend(x(ione,jone), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)
    
  End Subroutine D_ISend_2D

Subroutine D_ISend_3D(x, irq,n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(1:,1:,1:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:3)
	 Integer :: istart, kstart, jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2,irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
 	If (Present(indstart)) Then

		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
	Else
		istart = 1
		jstart = 1
		kstart = 1
	Endif   
    Call mpi_isend(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)
    !write(6,*)'zs ', p
  End Subroutine D_ISend_3D

Subroutine Z_ISend_2D(x, irq,n_elements, dest, tag, grp)
    Complex*16, Intent(in)  :: x(:,:)

    Integer, Optional :: dest, n_elements, tag
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
    
    Call mpi_isend(x(1,1), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)
    
  End Subroutine Z_ISend_2D

Subroutine Z_ISend_3D(x, irq,n_elements, dest, tag, grp, indstart)
    Complex*16, Intent(in)  :: x(:,:,:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:3)
	 Integer :: istart, kstart, jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2,irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
 	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
	Else
		istart = 1
		jstart = 1
		kstart = 1
	Endif   
    Call mpi_isend(x(istart,jstart,kstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)
    !write(6,*)'zs ', p
  End Subroutine Z_ISend_3D

	Subroutine D_IReceive_1D(x, irq, n_elements, istart, source, tag, grp)
		Real(8), Intent(out)  :: x(:)

    Integer, Optional :: source, n_elements, tag, istart
	 Integer ::  ione
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq, status(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If


	 If (Present(istart)) Then
		ione = istart
	 Else
		ione = 1
	 Endif
    
    Call mpi_irecv(x(ione), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
    !write(6,*)'zr ', p
  End Subroutine D_IReceive_1D

	Subroutine Z_IReceive_1D(x, irq, n_elements, istart, source, tag, grp)
		Complex*16, Intent(out)  :: x(:)

    Integer, Optional :: source, n_elements, tag, istart
	 Integer ::  ione
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq, status(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If


	 If (Present(istart)) Then
		ione = istart
	 Else
		ione = 1
	 Endif
    
    Call mpi_irecv(x(ione), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)

  End Subroutine Z_IReceive_1D

	Subroutine D_IReceive_2D(x, irq, n_elements, source, tag, grp, indstart)
		Real(8), Intent(out)  :: x(1:,1:)
    Integer, Optional :: source, n_elements, tag, indstart(1:2)
	 Integer :: ione, jone
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq, status(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

	 If (present(indstart)) Then
		ione = indstart(1)
		jone = indstart(2)
	 Else
		ione = 1
		jone = 1
	 Endif
	
    
    Call mpi_irecv(x(ione,jone), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
  End Subroutine D_IReceive_2D

Subroutine D_IReceive_3D(x, irq,n_elements, source, tag, grp,indstart)
    Real*8, Intent(out)  :: x(1:,1:,1:)

    Integer, Optional :: source, n_elements, tag,indstart(1:3)
	 Integer :: istart,jstart,kstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
	Else
		istart = 1
		jstart = 1
		kstart = 1
	Endif
    
    Call mpi_irecv(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
     !   write(6,*)'zr ', p
  End Subroutine D_IReceive_3D

	Subroutine Z_IReceive_2D(x, irq, n_elements, source, tag, grp, indstart)
		Complex*16, Intent(out)  :: x(:,:)

    Integer, Optional :: source, n_elements, tag, indstart(1:2)
	 Integer :: istart, jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, irq, status(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
	Else
		istart = 1
		jstart = 1
	Endif
    
    Call mpi_irecv(x(istart,jstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
    !write(6,*)'zr ', p
  End Subroutine Z_IReceive_2D

Subroutine Z_IReceive_3D(x, irq,n_elements, source, tag, grp,indstart)
    Complex*16, Intent(out)  :: x(:,:,:)

    Integer, Optional :: source, n_elements, tag,indstart(1:3)
	 Integer :: istart,jstart,kstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), irq

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm
    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

	If (Present(indstart)) Then
		istart = indstart(1)
		jstart = indstart(2)
		kstart = indstart(3)
	Else
		istart = 1
		jstart = 1
		kstart = 1
	Endif
    
    Call mpi_irecv(x(istart,jstart,kstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    If (p == MPI_ANY_SOURCE) source = status(MPI_SOURCE)
    If (tag2 == MPI_ANY_TAG) tag = status(MPI_TAG)
     !   write(6,*)'zr ', p
  End Subroutine Z_IReceive_3D

End Module ISendReceive
