Module All_To_All

  Use MPI_BASE

  Implicit None

  Private
  Public :: Standard_Transpose





	Interface Standard_Transpose
		Module Procedure Z_Transpose_v_1D, D_Transpose_v_1D, D_Transpose_choose_1D
	End Interface



Contains

Subroutine Z_Transpose_v_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp)
    Complex*16 :: send_buf(:), recv_buf(:)

    Integer, Intent(in) :: send_count(:), send_displ(:)
    Integer, Intent(in) :: recv_count(:), recv_displ(:)
    Type(communicator), Intent(in) :: grp
    Integer :: MPI_err

   call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_COMPLEX, recv_buf, &
         & recv_count, recv_displ, MPI_DOUBLE_COMPLEX, grp%comm, MPI_err)
  End Subroutine Z_Transpose_v_1D

Subroutine D_Transpose_v_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp)
    Real*8 :: send_buf(:), recv_buf(:)

    Integer, Intent(in) :: send_count(:), send_displ(:)
    Integer, Intent(in) :: recv_count(:), recv_displ(:)
    Type(communicator), Intent(in) :: grp
    Integer :: MPI_err

   call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_PRECISION, recv_buf, &
         & recv_count, recv_displ, MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
End Subroutine D_Transpose_v_1D


Subroutine D_Transpose_choose_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp, normal)
    Real*8 :: send_buf(:), recv_buf(:)

    Integer, Intent(in) :: send_count(:), send_displ(:)
    Integer, Intent(in) :: recv_count(:), recv_displ(:)
    Type(communicator), Intent(in) :: grp
    Integer :: MPI_err
	 Logical, Intent(in) :: normal
	 If (normal) Then
		Write(6,*)'Normal alltoall'
   	call MPI_ALLTOALL(send_buf, send_count(1), MPI_DOUBLE_PRECISION, recv_buf, &
         & recv_count(1), MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
	 Else

   call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_PRECISION, recv_buf, &
         & recv_count, recv_displ, MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
	 Endif
End Subroutine D_Transpose_choose_1D

End Module All_To_All
