Module Communication_Library

	Implicit None
	Integer, Public :: nproc, myid
	include 'mpif.h'

Contains

	Subroutine Initialize_Communication
		Implicit None
		Integer :: err,i
	   Call mpi_init(err)	
    	Call mpi_comm_size(MPI_COMM_WORLD, nproc, err)
    	Call mpi_comm_rank(MPI_COMM_WORLD, myid, err)
	End Subroutine Initialize_Communication

	Subroutine Finalize_Communication
		Implicit None
		Integer :: err	 
		Call mpi_finalize(err)
	End Subroutine Finalize_Communication

	! Use this to kill all processes somewhat cleanly
	SUBROUTINE Kill_All (proc)
		IMPLICIT NONE
		INTEGER :: ierr, proc
		WRITE(*,'(A,I0,A)') "Proc ",proc," called Kill_All. Aborting MPI.."
		CALL MPI_Abort(MPI_COMM_WORLD, ierr)
	END SUBROUTINE Kill_All

End Module Communication_Library
