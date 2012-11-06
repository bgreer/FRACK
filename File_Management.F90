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
Contains

	! Reads a single dopplergram into memory
	! INPUT: fname, myid, verbose
	! OUTPUT: arr(4096,4096), everything else
	SUBROUTINE Read_Dopplergram(fname, arr, dop_cen_lon, dop_cen_lat, &
		dop_cen_xpix, dop_cen_ypix, dop_p_angle, dop_r_sun_pix, myid, verbose)
		IMPLICIT NONE
		CHARACTER(LEN=200) :: fname, record
		REAL :: arr(4096,4096)
		REAL, ALLOCATABLE :: buff(:)
		INTEGER :: stat, blocksize, nfnd, nv, g, offset
		INTEGER :: naxes(2), ii, myid
		LOGICAL :: anyf, verbose
		REAL :: dop_cen_lon, dop_cen_lat
		REAL :: dop_cen_xpix, dop_cen_ypix
		REAL :: dop_p_angle, dop_r_sun_pix

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

		ALLOCATE(buff(naxes(1)))
		nv = -999
		g = 1
		offset = 1
		buff = 0
		DO ii=1,naxes(2)
			CALL FTGPVE(30,g,offset,naxes(1),nv,buff,anyf,stat)
			arr(ii,:) = buff
			offset = offset + naxes(1)
		ENDDO
		DEALLOCATE(buff)

		! load important keys
		CALL FTGKYE(30,"CRLN_OBS",dop_cen_lon,record,stat) ! is this right?
		CALL FTGKYE(30,"CRLT_OBS",dop_cen_lat,record,stat)
		CALL FTGKYE(30,"CRPIX1",dop_cen_xpix,record,stat)
		CALL FTGKYE(30,"CRPIX2",dop_cen_ypix,record,stat)
		CALL FTGKYE(30,"CROTA2",dop_p_angle,record,stat)
		CALL FTGKYE(30,"RSUN_OBS",dop_r_sun_pix,record,stat)

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

		! Copy old header info TODO: decide what is needed
!		CALL FTRDEF(20,stat)
!		CALL FTGHSP(30,nkeys,nmore,stat)
!		DO ii=7,nkeys
!		   CALL FTGREC(30,ii,record,stat)
!		   CALL FTPREC(20,record,stat)
!		ENDDO
!		CALL FTGKYE(30,"MAPSCALE",mapscale,record,stat)
!		dk = 360./(mapscale*npix*696.0)
!		PRINT*, " Adding header keys:"
!		PRINT*, "   DELTA_NU=",dnu
!		PRINT*, "   DELTA_K=",dk
!		CALL FTPKYE(20,"DELTA_NU",dnu,6,"",stat)
!		CALL FTPKYE(20,"DELTA_K",dk,6,"",stat)
!		CALL FTCLOS(30,stat)
		g = 1
		offset = 1
		DO ii=1,naxes(3) ! TODO: speed this up?
			DO ij=1,naxes(2)
				CALL FTPPRE(20,g,offset,naxes(1),arr(:,ij,ii),stat)
				offset = offset + naxes(1)
			ENDDO
		ENDDO
		CALL FTCLOS(20,stat)
		
	END SUBROUTINE

End Module File_Management
