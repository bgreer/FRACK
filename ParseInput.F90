
! Contains global variables, tracking settings, etc.

MODULE ParseInput

	IMPLICIT NONE

	LOGICAL :: verbose
	INTEGER :: nsteps ! number of timesteps
	INTEGER :: trackrate ! 0=carrington, 1=snodgrass, 2=custom
	INTEGER :: crot, cmlon ! defines the central time
	LOGICAL :: dotilesize(6) ! dotilesize(i) is for 2^(i-1) degrees
	REAL :: lonrn, latrn, clon, clat ! ranges and center lat/lons
	REAL :: memlimittotal ! max GB of memory to allocate for tiles
	CHARACTER(LEN=200) :: masterlist, background
	INTEGER :: backgroundset
	CHARACTER(LEN=400) :: outdir
	INTEGER :: loaddops ! number of dopplergrams to keep in memory
	REAL :: apode ! apodization something
	REAL :: a0,a2,a4 ! tracking rate
	LOGICAL :: dotiming ! record timing info
	REAL :: densepack ! fraction to overlap by

CONTAINS

	SUBROUTINE Usage()
		IMPLICIT NONE
		PRINT*, "Minimal Usage: ./frack"
		PRINT*, "Other Options:"
		PRINT*, " -v      Verbose Mode"
		PRINT*, " -r [#]  Tracking rate (0=car, 1=snod, 2=custom)"
		PRINT*, " -ml [file]  Master list of dopplergrams"
		PRINT*, " -bk [file]  Background fits file"
		PRINT*, " -ts[32,16,8,4,2,1]  Do tilesize [#]"
		PRINT*, " -clon [#], -clat [#]  Central lon/lat"
		PRINT*, " -lonrn [#], -latrn [#]  Lon/lat ranges"
		PRINT*, " -memlimit [#]  Total memory limit in GB"
		PRINT*, " -loaddops [#]  Number of dopplergrams to load at a time"
		PRINT*, " -time  Record timing info for performance analysis"
		PRINT*, " -outdir [dir]  Directory ti save output in"
	END SUBROUTINE Usage

	! Note the lack of IMPLICIT NONE
	! Poor coding? Maybe.
	SUBROUTINE ReadCommandLine()
		INTEGER :: ii, argcount, currarg, tempint, tempint2
		REAL :: tempreal
		CHARACTER(LEN=200) :: strbuffer

		argcount = IARGC()

		! Read command line arguments
		DO ii=1,argcount
			CALL getarg(ii,strbuffer)

			IF (strbuffer .EQ. "--help" .OR. strbuffer .EQ. "-h") THEN
				CALL Usage()
				STOP
			ENDIF

			! verbose mode
			IF (strbuffer .EQ. "-v") THEN
				verbose = .TRUE.

			! parse time
			ELSEIF (strbuffer .EQ. "-t") THEN
				CALL getarg(ii+1,strbuffer)
				IF (INDEX(strbuffer,":") .EQ. 5) THEN ! assume crot has 4 digits
					READ(strbuffer(1:4),*) tempint
					READ(strbuffer(6:LEN_TRIM(strbuffer)),*) tempint2
					IF (tempint .GT. 2000 .AND. tempint .LT. 3000 .AND. &
							tempint2 .GE. 0 .AND. tempint2 .LE. 360) THEN
						crot = tempint
						cmlon = tempint2
					ELSE
						PRINT*, "Invalid Time"
					ENDIF
				ELSE
					PRINT*, "Invalid Time. Example:  2096:180"
				ENDIF

			! number of time steps
			ELSEIF (strbuffer .EQ. "-l") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempint
				IF (tempint .GE. 3 .AND. tempint .LE. 15000) THEN
					nsteps = tempint
				ELSE
					PRINT*, "Invalid number of steps, using default."
				ENDIF
			
			! tracking rate
			ELSEIF (strbuffer .EQ. "-r") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) trackrate
				IF (trackrate .EQ. 0) THEN ! carrington
					a0 = 0D0
					a2 = 0D0
					a4 = 0D0
				ELSE IF (trackrate .EQ. 1) THEN ! snodgrass
					a0 = -0.02893D0
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 2) THEN ! custom
					a0 = -0.04330D0 ! snodgrass + 10 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 3) THEN ! custom
					a0 = -0.06485D0 ! snodgrass + 25 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 4) THEN ! custom
					a0 = -0.10077D0 ! snodgrass + 50 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 5) THEN ! custom
					a0 = -0.17261D0 ! snodgrass + 100 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 6) THEN ! custom
					a0 = -0.38813D0 ! snodgrass + 250 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 7) THEN ! custom
					a0 = 0.33027D0 ! snodgrass - 250 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 8) THEN ! custom
					a0 = -0.747321D0 ! snodgrass + 500 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE IF (trackrate .EQ. 9) THEN ! custom
					a0 = 0.689461D0 ! snodgrass - 500 m/s
					a2 = -0.341D0
					a4 = -0.5037D0
				ELSE
					PRINT*, "Invalid tracking rate, defaulting to snodgrass"
					a0 = -0.02893D0
					a2 = -0.341D0
					a4 = -0.5037D0
				ENDIF
			
			! dopplergram master list
			ELSEIF (strbuffer .EQ. "-ml") THEN
				masterlist = ""
				CALL getarg(ii+1,masterlist)

			! background file
			ELSEIF (strbuffer .EQ. "-bk") THEN
				background = ""
				CALL getarg(ii+1,background)
				backgroundset = 1

			! tile save directory
			ELSEIF (strbuffer .EQ. "-outdir") THEN
				outdir = ""
				CALL getarg(ii+1,outdir)

			! do tilesizes
			ELSEIF (strbuffer .EQ. "-ts32") THEN
				dotilesize(6) = .TRUE.
			ELSEIF (strbuffer .EQ. "-ts16") THEN
				dotilesize(5) = .TRUE.
			ELSEIF (strbuffer .EQ. "-ts8") THEN
				dotilesize(4) = .TRUE.
			ELSEIF (strbuffer .EQ. "-ts4") THEN
				dotilesize(3) = .TRUE.
			ELSEIF (strbuffer .EQ. "-ts2") THEN
				dotilesize(2) = .TRUE.
			ELSEIF (strbuffer .EQ. "-ts1") THEN
				dotilesize(1) = .TRUE.
			
			! central lat/lon
			ELSEIF (strbuffer .EQ. "-clon") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				clon = tempreal
			ELSEIF (strbuffer .EQ. "-clat") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				clat = tempreal

			! lat/lon ranges
			ELSEIF (strbuffer .EQ. "-lonrn") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				lonrn = tempreal
			ELSEIF (strbuffer .EQ. "-latrn") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				latrn = tempreal

			! memory limit
			ELSEIF (strbuffer .EQ. "-memlimit") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				memlimittotal = tempreal

			! number of dops to load at a time
			ELSEIF (strbuffer .EQ. "-loaddops") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempint
				loaddops = tempint
			
			! record timing info
			ELSEIF (strbuffer .EQ. "-time") THEN
				dotiming = .TRUE.

			! densepack
			ELSEIF (strbuffer .EQ. "-densepack") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) tempreal
				IF (tempreal .GT. 0D0 .AND. tempreal .LT. 1D0) THEN
					densepack = tempreal
				ENDIF

			ENDIF
		ENDDO
	END SUBROUTINE ReadCommandLine



	! Set the default values of various run parameters
	! It is best to call this before reading in any conflicting information..
	SUBROUTINE SetDefaults()
		verbose = .FALSE.
		nsteps = 2048
		trackrate = 1
		a0 = -0.02893D0
		a2 = -0.341D0
		a4 = -0.5037D0
		lonrn = 0D0
		latrn = 0D0
		clon = 120D0
		clat = 0D0
		dotilesize(1) = .FALSE.
		dotilesize(2) = .FALSE.
		dotilesize(3) = .FALSE.
		dotilesize(4) = .FALSE.
		dotilesize(5) = .FALSE.
		dotilesize(6) = .FALSE.
		memlimittotal = 8.0
		masterlist = "doplist"
		loaddops = 8
		apode = 0.9375D0
		dotiming = .FALSE.
		densepack = 0.5D0
		background = ""
		backgroundset = 0
	END SUBROUTINE SetDefaults

	! Print the run details to stdout
	SUBROUTINE PrintDetails()
		INTEGER :: ii
		PRINT*, "Run Details:"
!		WRITE(*,'(A,I4)') "  Carrington Rotation = ",crot
!		WRITE(*,'(A,I3)') "  Central Meridian Longitude = ", cmlon
		SELECT CASE (trackrate)
			CASE (0)
				PRINT*, "  Tracking at Carrington Rate"
			CASE (1)
				PRINT*, "  Tracking at Snodgrass Rate"
			CASE (2)
				PRINT*, "  Tracking at Custom Rate"
		END SELECT
		WRITE(*,'(A,A)') "  Master List File: ", TRIM(masterlist)
		WRITE(*,'(A$)') "  Tracking Tilesizes: "
		DO ii=1,6
			IF (dotilesize(ii)) THEN
				WRITE(*,'(I0,A$)'), 2**(ii-1), " "
			ENDIF
		ENDDO
		WRITE(*,'(A)') ""
		PRINT*, "  Central Lon/Lat: ", clon, clat
		PRINT*, "  Lon/Lat Ranges: ", lonrn, latrn
		WRITE(*,'(A)') ""
	END SUBROUTINE PrintDetails


END MODULE
