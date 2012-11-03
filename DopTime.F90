
! the purpose of this is to enable scripts to extract the time info from
!  dopplergrams

! read a dopplergram and output time info
PROGRAM DopTime
	IMPLICIT NONE

	INTEGER :: argcount, stat, blocksize
	CHARACTER(LEN=200) :: fname, record, tobs, trec
	REAL :: cmlon
	INTEGER :: crot

	! Check number of arguments
	argcount = IARGC()
	IF (argcount .lt. 1) THEN
		PRINT*, "Usage: ./DopTime filename"
		PRINT*, "Output:"
		PRINT*, "  T_OBS"
		PRINT*, "  T_REC"
		PRINT*, "  CAR_ROT"
		PRINT*, "  CRLN_OBS"
		PRINT*, "  Carrington Time"
		STOP
	ENDIF

	! Get the filename
	CALL GETARG(1,fname)

	! Open FITS file
	stat = 0
	CALL FTOPEN(30,fname,0,blocksize,stat)

	IF (stat .NE. 0) THEN
		PRINT*, "Error reading dopplergram."
		STOP
	ENDIF

	! read keys
	CALL FTGKYE(30,"CRLN_OBS",cmlon,record,stat)
	CALL FTGKYJ(30,"CAR_ROT", crot, record,stat)
	CALL FTGKYS(30, "T_OBS", tobs, record, stat)
	CALL FTGKYS(30, "T_REC", trec, record, stat)

	! Close FITS file
	CALL FTCLOS(30,stat)

	! Output keys
	WRITE(*,'(A)') TRIM(tobs)
	WRITE(*,'(A)') TRIM(trec)
	WRITE(*,'(I0)') crot
	WRITE(*,'(F12.7)') cmlon
	WRITE(*,'(F13.7)') FLOAT(crot)+1D0-(cmlon/360D0)

END PROGRAM
