Module General_MPI
	Use MPI_BASE
  Implicit None
Contains

Subroutine Global_Max(sendbuf, recvbuf, grp)
  Real*8 :: sendbuf, recvbuf
  Type(communicator), Optional :: grp
  Integer :: icount,  comm
  Integer :: MPI_err

  icount = 1
  
  If (Present(grp)) Then
     comm = grp%comm
  Else
     comm = MPI_COMM_WORLD
  End If


  Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_MAX, comm, MPI_err)
  
End Subroutine Global_Max

Subroutine Global_IMax(sendbuf, recvbuf, grp)
  Integer*4 :: sendbuf, recvbuf
  Type(communicator), Optional :: grp
  Integer :: icount,  comm
  Integer :: MPI_err

  icount = 1
  
  If (Present(grp)) Then
     comm = grp%comm
  Else
     comm = MPI_COMM_WORLD
  End If


  Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_INTEGER, MPI_MAX, comm, MPI_err)
  
End Subroutine Global_IMax

Subroutine DSUM1D(sendbuf, recvbuf, grp, dest)
  Real*8 :: sendbuf(:), recvbuf(:)
  Type(communicator), Optional :: grp
  Integer, Intent(In), Optional :: dest
  Integer :: icount,  comm
  Integer :: MPI_err, root
  If (present(dest)) then
    root = dest
  Else
    root = 0
  Endif

  icount = size(sendbuf)
  
  If (Present(grp)) Then
     comm = grp%comm
  Else
     comm = MPI_COMM_WORLD
  End If


  Call MPI_REDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, MPI_err)
  
End Subroutine DSUM1D

Subroutine DSUM2D(sendbuf, recvbuf, grp, dest)
  Real*8 :: sendbuf(:,:), recvbuf(:,:)
  Type(communicator), Optional :: grp
  Integer, Intent(In), Optional :: dest
  Integer :: icount,  comm
  Integer :: MPI_err, root
  If (present(dest)) then
    root = dest
  Else
    root = 0
  Endif

  icount = size(sendbuf)
  
  If (Present(grp)) Then
     comm = grp%comm
  Else
     comm = MPI_COMM_WORLD
  End If


  Call MPI_REDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, MPI_err)
  
End Subroutine DSUM2D

End Module General_MPI
