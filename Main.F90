
! CU Tracking program

! Basic code flow:

! Make list of all tiles to track
! Make list of all relevant dopplergrams
! For each tilesize
! {
!   While there are untracked tiles
!   {
!     Let each processor track some number of tiles that fits in memory
!     Let each processor save the recently tracked tiles
!     Increment the number of completed tiles
!   }
! }

! Loading dopplergrams happens in Track.F90 and is done only 
!  when the dopplergrams are needed for tracking.
!  (dopplergram loading is cheap)

PROGRAM FRACK

	USE ParseInput
	USE Communication_Library
	USE Track
	USE Grid
	USE File_Management
	USE Timing
	IMPLICIT NONE

	INTEGER :: ierr, ii, ij, it, tilecount, currsize
	INTEGER :: currtile, ix, iy, completed, tilesperproc, starttile
	INTEGER :: endtile, maxtiles, num

	! Init MPI and timer
		! Communication_Library.F90
	CALL Initialize_Communication()
		! Timing.F90
	CALL InitTimer(myid,10)
	CALL AddTime(1)

	! Set defaults
		! ParseInput.F90
	CALL SetDefaults()
	! Sort out command line arguments
	CALL ReadCommandLine()
	CALL AddTime(2)

	! Output some stats to stdout
	CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
	IF (verbose .AND. myid .EQ. 0) &
		WRITE(*,'(I3,A)') nproc," processors checking in."

	IF (verbose .AND. myid .EQ. 0) CALL PrintDetails()


	! Count number of tiles needed for each tilesize
	CALL CountTiles(myid, tilecount, dotilesize, numtiles, lonrn, latrn, apode, verbose)
	! Allocate info foreach tile
	ALLOCATE(tilesize(tilecount))
	ALLOCATE(lon(tilecount))
	ALLOCATE(lat(tilecount))
	ALLOCATE(delta_rot(tilecount))
	
	! Go back through and load details of each tile
	currtile = 1
	DO ii=6,1,-1
		IF (dotilesize(ii)) THEN
			currsize = 2**(ii-1)
			tilesize(currtile:currtile+numtiles(ii)-1) = currsize
				! Grid.F90
			CALL GenerateGrid(currsize,lonrn,latrn,clon,clat,apode,&
				lon(currtile:currtile+numtiles(ii)-1),lat(currtile:currtile+numtiles(ii)-1))
				! Grid.F90
			CALL SetTrackingRate(delta_rot,a0,a2,a4,lat,numtiles(ii))
			currtile = currtile + numtiles(ii)
		ENDIF
	ENDDO
	CALL AddTime(3)

	! Let proc0 read the master dop list
	CALL AddTime(4)
	IF (myid .EQ. 0) THEN
		IF (verbose) WRITE(*,'(A)') "Reading Master List.."
			! File_Management.F90
!		CALL Make_Dopplergram_List(myid, crot, cmlon, nsteps, masterlist, &
!		dopfname, doptime, dopinterp, dopfname_ends, doptime_ends)
		CALL Load_Dopplergram_List(myid, nsteps, masterlist, dopfname, doptime, dopinterp)
	ENDIF
	CALL AddTime(5)

	! Communicate master dop list to all procs
	IF (myid .EQ. 0 .AND. verbose) &
		WRITE(*,'(A)') "Broadcasting dopplergram list to all procs.."
	CALL MPI_BCAST(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	IF (myid .NE. 0) THEN
		ALLOCATE(doptime(nsteps))
		ALLOCATE(dopfname(nsteps))
		ALLOCATE(dopinterp(nsteps))
	ENDIF
	ALLOCATE(dop_p_angle(nsteps))
	ALLOCATE(dop_cen_lon(nsteps))
	ALLOCATE(dop_cen_lat(nsteps))
	ALLOCATE(dop_cen_xpix(nsteps))
	ALLOCATE(dop_cen_ypix(nsteps))
	ALLOCATE(dop_r_sun_pix(nsteps))
	CALL MPI_BCAST(doptime, nsteps, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!	CALL MPI_BCAST(doptime_ends, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	CALL MPI_BCAST(dopfname, nsteps*200, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
!	CALL MPI_BCAST(dopfname_ends, 400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
	CALL MPI_BCAST(dopinterp, nsteps, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
	CALL AddTime(6)

	! Can't interpolate first or last dopplergram
	IF (dopinterp(1) .OR. dopinterp(nsteps)) THEN
		WRITE(*,'(A)') "FATAL ERROR: Can't interpolate first or last dopplergram :("
		STOP
	ENDIF

	! Main loop through tile size
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (verbose .AND. myid .EQ. 0) PRINT*, "Entering Main Loop"
	currtile = 1
	DO ii=6,1,-1
		CALL AddTime(1000+ii)
		IF (dotilesize(ii)) THEN
			completed = 0 ! completed of this tilesize
			currsize = 2**(ii-1)

			IF (verbose .AND. myid .EQ. 0) WRITE(*,'(A,I0,A,I0,A)') &
					"Tracking tilesize ",currsize," (",numtiles(ii)," tiles)"

			! fill up each processor and run tracking
			! until no tiles are left
			maxtiles = FLOOR(MIN(memlimit,memlimittotal/nproc) &
				/memusage(tilesize(currtile), nsteps))
			
			DO WHILE (completed .LT. numtiles(ii))
				
				! fill up each processor with tiles
				num = MIN(maxtiles*nproc, numtiles(ii)-completed)
				
				! call main tracking
					! Track.F90
				CALL TrackTiles(myid,nproc,num,currsize,nsteps,completed,&
					dopfname,doptime,dopfname_ends,doptime_ends,dopinterp,loaddops,verbose)
				completed = completed + num
				currtile = currtile + num

				CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
			
			ENDDO
		ENDIF
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL AddTime(7)

	! Free ALL the things
	DEALLOCATE(tilesize)
	DEALLOCATE(lon)
	DEALLOCATE(lat)
	DEALLOCATE(delta_rot)
	DEALLOCATE(doptime)
	DEALLOCATE(dopfname)
	DEALLOCATE(dopinterp)
	DEALLOCATE(dop_p_angle)
	DEALLOCATE(dop_cen_lon)
	DEALLOCATE(dop_cen_lat)
	DEALLOCATE(dop_cen_xpix)
	DEALLOCATE(dop_cen_ypix)
	DEALLOCATE(dop_r_sun_pix)

	! End everything
	IF (verbose .AND. myid .EQ. 0) PRINT*, "Main loop done, Finalizing MPI.."
	CALL AddTime(8)
	CALL SaveTimer(myid, nproc, "timing")
	CALL KillTimer()
	CALL Finalize_Communication()

END PROGRAM
