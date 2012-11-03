
! Attempt at a timing module
! May not be the best, but it doesn't really matter for routine tracking jobs

MODULE Timing
	USE Communication_Library
	IMPLICIT NONE

	INTEGER :: currtime, timesalloc
	! 'flag' for each timestamp
	INTEGER, ALLOCATABLE :: timeflag(:), timeflag2(:)
	! timestamps
	DOUBLE PRECISION, ALLOCATABLE :: timearr(:), timearr2(:)

CONTAINS

	SUBROUTINE InitTimer (myid, num)
		INTEGER :: myid, num
		timesalloc = num
		currtime = 1
		ALLOCATE(timearr(timesalloc))
		ALLOCATE(timeflag(timesalloc))
		timeflag(:) = -1
	END SUBROUTINE InitTimer

	SUBROUTINE KillTimer ()
		DEALLOCATE(timearr)
		DEALLOCATE(timeflag)
	END SUBROUTINE KillTimer

	! Add a single timestamp with flag
	SUBROUTINE AddTime (flag)
		INTEGER :: flag
		timearr(currtime) = MPI_Wtime()
		timeflag(currtime) = flag
		currtime = currtime + 1
		IF (currtime .GT. timesalloc) &
			CALL ResizeArray()
	END SUBROUTINE AddTime

	! Dynamic array sizing courtesy programming 101
	SUBROUTINE ResizeArray ()
		ALLOCATE(timearr2(timesalloc))
		ALLOCATE(timeflag2(timesalloc))
		timearr2(:) = timearr(:)
		timeflag2(:) = timeflag(:)
		DEALLOCATE(timearr)
		DEALLOCATE(timeflag)
		ALLOCATE(timearr(timesalloc*2))
		ALLOCATE(timeflag(timesalloc*2))
		timeflag(:) = -1
		timearr(1:timesalloc) = timearr2(:)
		timeflag(1:timesalloc) = timeflag2(:)
		DEALLOCATE(timearr2)
		DEALLOCATE(timeflag2)
		timesalloc = timesalloc * 2
	END SUBROUTINE ResizeArray

	! Save timing array to file
	! have processors take turns appending to file
	SUBROUTINE SaveTimer (myid, numproc, fname)
		INTEGER :: myid, numproc, ii, ij, ierr
		CHARACTER(LEN=*) :: fname
	
		DO ii=0,nproc-1
			IF (myid .EQ. ii) THEN
				! append to file
				OPEN(4,FILE=fname,ACCESS='APPEND')
					DO ij=1,currtime-1
						WRITE(4,'(I6,I6,I6,F18.5)') myid, ij, timeflag(ij),&
							timearr(ij)
					ENDDO
				CLOSE(4)
			ENDIF
			CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
		ENDDO
	END SUBROUTINE

END MODULE Timing
