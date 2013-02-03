Module File_Management
	USE Communication_Library
	Implicit None

	! Data for dopplergrams
	LOGICAL, ALLOCATABLE :: dopinterp(:)
	CHARACTER(LEN=200), ALLOCATABLE :: dopfname(:)
	INTEGER, ALLOCATABLE :: doptime(:)
	REAL, ALLOCATABLE :: dop_cen_lon(:), dop_cen_lat(:)
	REAL, ALLOCATABLE :: dop_cen_xpix(:), dop_cen_ypix(:)
	REAL, ALLOCATABLE :: dop_p_angle(:), dop_r_sun_pix(:)

	! needed in case of interpolation on ends
	CHARACTER(LEN=200) :: dopfname_ends(2)
	INTEGER :: doptime_ends(2)

	INTEGER, DIMENSION(6) :: numtiles ! for each tile size

	! Data for tiles
	INTEGER, ALLOCATABLE :: tilesize(:)
	REAL, ALLOCATABLE :: lon(:), lat(:)
	REAL, ALLOCATABLE :: delta_rot(:)

	REAL :: backframe(4096,4096)
Contains

	SUBROUTINE Read_Background(myid, fname, arr)
		IMPLICIT NONE
		CHARACTER(LEN=200) :: fname
		REAL :: arr(4096,4096)
		INTEGER :: nv, g, stat, blocksize, myid, nfnd
		INTEGER :: naxes(2)
		LOGICAL :: anyf

		stat = 0
		! Open input FITS file
		CALL FTOPEN(40,fname,0,blocksize,stat)
		IF (stat .NE. 0) THEN
			WRITE(*,'(A,A)') "Error reading background file: ",fname
			CALL Kill_All(myid)
		ENDIF
		! Get dimensions of data cube
		CALL FTGKNJ(40,'NAXIS',1,2,naxes,nfnd,stat)

		nv = -999
		g = 1
		CALL FTG2DE(40,g,nv,naxes(1),naxes(1),naxes(2),arr,anyf,stat)
		CALL FTCLOS(40,stat)
	END SUBROUTINE Read_Background

	! Reads a single dopplergram into memory
	! INPUT: fname, myid, verbose
	! OUTPUT: arr(4096,4096), everything else
	SUBROUTINE Read_Dopplergram(fname, arr, dop_cen_lon, dop_cen_lat, &
		dop_cen_xpix, dop_cen_ypix, dop_p_angle, dop_r_sun_pix, myid, verbose)
		IMPLICIT NONE
		CHARACTER(LEN=200) :: fname, record
		REAL :: arr(4096,4096)
		INTEGER :: stat, blocksize, nfnd, nv, g, offset
		INTEGER :: naxes(2), ii, myid
		LOGICAL :: anyf, verbose
		REAL :: dop_cen_lon, dop_cen_lat
		REAL :: dop_cen_xpix, dop_cen_ypix
		REAL :: dop_p_angle, dop_r_sun_pix
		REAL :: corr1, corr2, dsun_obs, RSUNM, temp, degrad
		RSUNM = 6.96D8
		degrad = 180D0/3.14159265D0

		stat = 0
		! Open input FITS file
		CALL FTOPEN(30,fname,0,blocksize,stat)
		IF (stat .NE. 0) THEN
			WRITE(*,'(A,A)') "Error reading file ",fname
			CALL Kill_All(myid)
		ELSE
			IF (verbose) WRITE(*,'(A,I0,A,A)') "Proc",myid," reading file ",TRIM(fname)
		ENDIF
		! Get dimensions of data cube
		CALL FTGKNJ(30,'NAXIS',1,2,naxes,nfnd,stat)

		nv = -999
		g = 1
		offset = 1
		CALL FTG2DE(30,g,nv,naxes(1),naxes(1),naxes(2),arr,anyf,stat)

		! load important keys
		CALL FTGKYE(30,"CRLN_OBS",dop_cen_lon,record,stat) ! is this right?
		CALL FTGKYE(30,"CRLT_OBS",dop_cen_lat,record,stat)
		CALL FTGKYE(30,"CRPIX1",dop_cen_xpix,record,stat)
		CALL FTGKYE(30,"CRPIX2",dop_cen_ypix,record,stat)
		CALL FTGKYE(30,"CROTA2",dop_p_angle,record,stat)
		CALL FTGKYE(30,"RSUN_OBS",dop_r_sun_pix,record,stat)
		CALL FTGKYE(30,"DSUN_OBS",dsun_obs,record,stat)
		CALL FTGKYE(30,"CDELT1",corr1,record,stat) ! plate scale (arcsec per pixel, x-direction)
		CALL FTGKYE(30,"CDELT2",corr2,record,stat)

		temp = dop_r_sun_pix / corr1
		
		dop_r_sun_pix = asin(RSUNM/dsun_obs) ! r_sun in radians
		dop_r_sun_pix = dop_r_sun_pix * degrad ! convert to degrees
		dop_r_sun_pix = dop_r_sun_pix * 3600D0 ! convert to arcseconds
		dop_r_sun_pix = dop_r_sun_pix / corr1 ! convert to pixels

		CALL FTCLOS(30,stat)

	END SUBROUTINE Read_Dopplergram

	! Given a CROT:CMLON time and nsteps, make a list of
	!  of dopplergram filenames. Flag missing ones for interpolation
	! INPUT: myid, crot, cmlon, nsteps
	! OUTPUT: masterlist, dop*
	SUBROUTINE Make_Dopplergram_List(myid, crot, cmlon, nsteps, masterlist, &
	dopfname, doptime, dopinterp, dopfname_ends, doptime_ends)
		IMPLICIT NONE

		INTEGER :: crot, cmlon, nsteps, stat, intbuffer, myid, ii
		CHARACTER(LEN=200) :: dopfname(nsteps), dopfname_ends(2)
		INTEGER :: doptime(nsteps), doptime_ends(2)
		INTEGER :: numentries, starttime, tdiff, closest
		CHARACTER(LEN=*) :: masterlist
		CHARACTER(LEN=200) :: strbuffer
		LOGICAL :: dopinterp(nsteps)

		! the master list should be a 2 column ascii text file
		! each line is a possible dopplergram
		! col1 is the filename on our system (or a 0 for missing file)
		! col2 is some integer identifier for time (unix time?)
		! assume list is sorted! for now?

		! decide the time for the first dopplergram
		! find middle time, then subtract 45*nsteps/2 seconds
		starttime = 1040 ! fix this

		! read master list and find closest to start
		closest = 1
		tdiff = 1000
		OPEN(2, FILE=masterlist)
		stat = 0
		numentries = 0
		DO WHILE (stat .EQ. 0)
			READ(2,*,IOSTAT=stat) strbuffer, intbuffer
			numentries = numentries + 1
			IF (stat .EQ. 0) THEN
				IF (ABS(intbuffer - starttime) .LT. tdiff) THEN
					closest = numentries
					tdiff = ABS(intbuffer-starttime)
				ENDIF
			ENDIF
		ENDDO
		CLOSE(2)
		numentries = numentries - 1

		! error checking
		IF (closest-1 .LT. 1 .OR. closest+nsteps .GT. numentries) THEN
			WRITE(*,'(A,I0)') "Too few dopplergrams in master list"
			CALL Kill_All(myid)
		ENDIF

		! currently, 'closest' is the index of the first dop in the
		!  master list, and 'tdiff' is how far off this

		! go back into the file and load information into memory
		OPEN(2, FILE=masterlist)
		! seek to closest-1
		DO ii=1,closest-2
			READ(2,*) strbuffer, intbuffer
		ENDDO
		! read closest-1 into _end(1)
		READ(2,*) strbuffer, intbuffer
		dopfname_ends(1) = TRIM(strbuffer)
		doptime_ends(1) = intbuffer
		! read nsteps into normal list
		DO ii=1,nsteps
			READ(2,*) strbuffer, intbuffer
			dopfname(ii) = TRIM(strbuffer)
			doptime(ii) = intbuffer
			IF (dopfname(ii) .EQ. '0') THEN
				dopinterp(ii) = .TRUE.
			ELSE
				dopinterp(ii) = .FALSE.
			ENDIF
		ENDDO
		! read end+1 into _end(2)
		READ(2,*) strbuffer, intbuffer
		dopfname_ends(2) = TRIM(strbuffer)
		doptime_ends(2) = intbuffer
		CLOSE(2)

		! if we need to inteprolate ends, make sure we can
		IF (dopinterp(1) .AND. dopfname_ends(1) .EQ. '0') THEN
			WRITE(*,'(A)') "Need more dopplergrams at start."
			CALL Kill_All(myid)
		ENDIF
		IF (dopinterp(nsteps) .AND. dopfname_ends(2) .EQ. '0') THEN
			WRITE(*,'(A)') "Need more dopplergrams at end."
			CALL Kill_All(myid)
		ENDIF
	END SUBROUTINE Make_Dopplergram_List

	! loads a list of dopplergrams, assume you are tracking through all of them
	! overwrites nsteps
	SUBROUTINE Load_Dopplergram_List(myid, nsteps, masterlist, dopfname, doptime, dopinterp)
		INTEGER :: myid, nsteps, ii, intbuffer, stat
		CHARACTER(LEN=*) :: masterlist
		CHARACTER(LEN=200), ALLOCATABLE :: dopfname(:)
		CHARACTER(LEN=200) :: strbuffer
		INTEGER, ALLOCATABLE :: doptime(:)
		LOGICAL, ALLOCATABLE :: dopinterp(:)

		! count number of lines
		OPEN(2, FILE=masterlist)
		nsteps = 0
		stat = 0
		DO WHILE (stat .EQ. 0)
			READ(2,*,IOSTAT=stat) strbuffer
			IF (stat .EQ. 0) nsteps = nsteps + 1
		ENDDO
		CLOSE(2)
		! Allocate space
		ALLOCATE(dopfname(nsteps))
		ALLOCATE(doptime(nsteps))
		ALLOCATE(dopinterp(nsteps))
		! read nsteps into normal list
		OPEN(2, FILE=masterlist)
		DO ii=1,nsteps
			READ(2,*) strbuffer
			dopfname(ii) = TRIM(strbuffer)
			doptime(ii) = 0
			IF (dopfname(ii) .EQ. '0') THEN
				dopinterp(ii) = .TRUE.
			ELSE
				dopinterp(ii) = .FALSE.
			ENDIF
		ENDDO
		CLOSE(2)
	END SUBROUTINE Load_Dopplergram_List

	! Save a single tile to disk
	! INPUT: arr, dim*, filename, vebrose
	! OUTPUT: (none)
	SUBROUTINE Save_Tile(arr, dim1, dim2, dim3, filename, verbose)
		IMPLICIT NONE

		CHARACTER(LEN=*) :: filename
		LOGICAL :: verbose, simple, extend
		INTEGER :: blocksize, stat, bitpix, offset, g
		INTEGER :: ii, ij, dim1, dim2, dim3
		INTEGER, DIMENSION(3) :: naxes
		REAL, DIMENSION(:,:,:) :: arr
		REAL :: dnu, dk, mapscale
	
		! Open file
		stat = 0
		CALL FTINIT(20,filename,blocksize,stat)
		! Set file attributes?
		simple=.TRUE.
		bitpix=-32
		extend=.TRUE.
		naxes(3) = dim3
		naxes(2) = dim2
		naxes(1) = dim1
		CALL FTPHPR(20,simple,bitpix,3,naxes,0,1,extend,stat)

		! make some new header entries
		mapscale = 1D0/24D0
		dnu = 1e6/(45.*dim3)
		dk = 360./(mapscale*dim1*696.0)
		! these are probably right..
		CALL FTPKYE(20,"DELTA_NU",dnu,6,"",stat)
		CALL FTPKYE(20,"DELTA_K",dk,6,"",stat) ! this probably isnt needed
		CALL FTPKYE(20,"MAPSCALE",mapscale,6,"",stat)
!		CALL FTCLOS(30,stat)
		g = 1
		offset = 1
!		CALL FTP3DE(20,g,naxes(1),naxes(2),naxes(1),naxes(2),naxes(3),arr,stat)
		DO ii=1,naxes(3) ! TODO: speed this up?
			DO ij=1,naxes(2)
				CALL FTPPRE(20,g,offset,naxes(1),arr(:,ij,ii),stat)
				offset = offset + naxes(1)
			ENDDO
		ENDDO
		CALL FTCLOS(20,stat)
		
	END SUBROUTINE

End Module File_Management
