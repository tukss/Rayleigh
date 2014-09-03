Module MPI_LAYER
	Use MPI_BASE  
	Use All_to_All
	Implicit None
	!Public :: Standard_Transpose
	



Contains


	Function Init_Main_Group( err) result(grp)
		Type(communicator) :: grp
		Integer, Intent(out) :: err
 		Call mpi_init(err)
		grp%comm = mpi_comm_world
    	Call mpi_comm_size(grp%comm, grp%np, err)
		Call mpi_comm_rank(grp%comm, grp%rank, err)

	End Function Init_Main_Group

	Subroutine RowColSplit(grp,rgrp,cgrp,nprow, err)
		! Take one group and split it using a row/column decomposition
		Type(communicator), Intent(InOut) :: grp, rgrp, cgrp
		Integer, Intent(out) :: err
		Integer, Intent(In) ::  nprow
		Integer :: row_rank, col_rank

		row_rank = mod(grp%rank,nprow)
		col_rank = grp%rank/nprow

		Call mpi_comm_split(grp%comm, col_rank, grp%rank, rgrp%comm, err)
		Call mpi_comm_split(grp%comm, row_rank, grp%rank, cgrp%comm, err)

		
    	Call mpi_comm_size(rgrp%comm, rgrp%np, err)
		Call mpi_comm_rank(rgrp%comm, rgrp%rank, err)

    	Call mpi_comm_size(cgrp%comm, cgrp%np, err)
		Call mpi_comm_rank(cgrp%comm, cgrp%rank, err)
		
		If (cgrp%rank .ne. col_rank) Write(6,*)'Error - ', cgrp%rank, col_rank
		If (rgrp%rank .ne. row_rank) Write(6,*)'Error - ', rgrp%rank, row_rank
	End Subroutine RowColSplit
	
	Subroutine Exit_Comm_Lib(err)
		Integer, Intent(out) :: err

		Call mpi_finalize(err)

	End Subroutine Exit_Comm_Lib

    Subroutine Barrier(grp)
        Implicit None
        Type(communicator) :: grp
        Integer :: ierr
        Call MPI_Barrier(grp%comm,ierr)
    End Subroutine Barrier

End Module MPI_LAYER
