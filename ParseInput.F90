
! Contains global variables, tracking settings, etc.

MODULE ParseInput

	IMPLICIT NONE

	LOGICAL :: verbose
	INTEGER :: nsteps ! number of timesteps
	INTEGER :: trackrate ! 0=carrington, 1=snodgrass, 2=custom
	INTEGER :: crot, cmlon ! defines the central time
	LOGICAL :: dotilesize(6) ! dotilesize(i) is for 2^(i-1) degrees
	REAL :: lonrn, latrn, clon, clat ! ranges and center lat/lons
	REAL :: memlimit, memlimittotal ! max GB of memory to allocate for tiles
	CHARACTER(LEN=200) :: masterlist
	INTEGER :: loaddops ! number of dopplergrams to keep in memory
	REAL :: apode ! apodization something
	REAL :: a0,a2,a4 ! tracking rate
	REAL :: a_carr ! Carrington rotation rate

CONTAINS

	SUBROUTINE Usage()
		IMPLICIT NONE
		PRINT*, "Minimal Usage: ./frack"
		PRINT*, "Other Options:"
		PRINT*, " -v      Verbose Mode"
		PRINT*, " -r [#]  Tracking rate (0=car, 1=snod, 2=custom)"
		PRINT*, " -ml [file]  Master list of dopplergrams"
		PRINT*, " -and others"
	END SUBROUTINE Usage

	! Note the lack of IMPLICIT NONE
	! Poor coding? Maybe.
	! TODO: add more command-line options
	SUBROUTINE ReadCommandLine()
		INTEGER :: ii, argcount, currarg, tempint, tempint2
		CHARACTER(LEN=200) :: strbuffer

		argcount = IARGC()

		! Read command line arguments
		DO ii=1,argcount
			CALL getarg(ii,strbuffer)

			IF (strbuffer .EQ. "--help") THEN
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
			
			! tracking rate TODO insert correct coefs
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
					a0 = 0D0
					a2 = 0D0
					a4 = 0D0
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
		a_carr = 2.6662237D0
		lonrn = 0D0
		latrn = 0D0
		clon = 120D0
		clat = 0D0
		dotilesize(1) = .FALSE.
		dotilesize(2) = .FALSE.
		dotilesize(3) = .FALSE.
		dotilesize(4) = .FALSE.
		dotilesize(5) = .TRUE.
		dotilesize(6) = .FALSE.
		memlimit = 4.0
		memlimittotal = 8.0
		masterlist = "doplist"
		loaddops = 8
		apode = 0.9375D0
	END SUBROUTINE SetDefaults

	! Print the run details to stdout
	SUBROUTINE PrintDetails()
		PRINT*, "Run Details:"
		WRITE(*,'(A,I4)') "  Carrington Rotation = ",crot
		WRITE(*,'(A,I3)') "  Central Meridian Longitude = ", cmlon
		SELECT CASE (trackrate)
			CASE (0)
				PRINT*, "  Tracking at Carrington Rate"
			CASE (1)
				PRINT*, "  Tracking at Snodgrass Rate"
			CASE (2)
				PRINT*, "  Tracking at Custom Rate"
		END SELECT
		WRITE(*,'(A,A)') "  Master List File: ", TRIM(masterlist)
	END SUBROUTINE PrintDetails


END MODULE
