
! The meat of the code

MODULE Track
	USE File_Management
	USE Timing
	USE Interpolation
	IMPLICIT NONE

CONTAINS

	! INPUT: everything
	! OUTPUT: nothing
	SUBROUTINE TrackTiles(myid,nproc,num,tilesize,nsteps,startindex,&
	dopfname, doptime, dopfname_ends, doptime_ends, dopinterp,loaddops,verbose)
		IMPLICIT NONE
		LOGICAL :: verbose
		INTEGER :: myid, nproc, num, tilesize, startindex,nsteps,ierr
		REAL, ALLOCATABLE :: tdata(:,:,:,:)
		INTEGER :: starttile, endtile, ix, numtiles, ii,ij, stepsdone
		CHARACTER(LEN=200) :: outfile, buff
		CHARACTER(LEN=200) :: dopfname(nsteps), dopfname_ends(2)
		INTEGER :: doptime(nsteps), doptime_ends(2)
		LOGICAL :: dopinterp(nsteps)
		INTEGER :: loaddops, numdops, dopstart, dopend
		REAL, ALLOCATABLE :: dopdata(:,:,:), dopdata_ends(:,:,:), tempslice(:,:)
		REAL :: tempdop(4096,4096)
		INTEGER :: sendcnts(nproc), sdispls(nproc)
		INTEGER :: recvcnts(nproc), rdispls(nproc)
		INTEGER :: ti, di, ii2
		REAL :: currlon, currlat, delta_time

		! some local definitions
		! let each proc decide what tiles to handle
		starttile = startindex+FLOOR(num*1.0*myid/nproc)+1
		endtile = startindex+FLOOR(num*1.0*(myid+1)/nproc)
		numtiles = endtile-starttile+1

		! tilesize in pixels
		ix = tilesize*24
		ALLOCATE(tempslice(ix,ix))

		! output how much memory will be used for tiles
		IF (verbose) WRITE(*,'(A,I0,A,I0,A,F6.3,A)') "Proc",myid," allocating " &
			,(endtile-starttile+1)," tiles, (",memusage(tilesize, nsteps)*numtiles,"GB)"

		! only allocate if we have tiles to track
		IF (endtile .GE. starttile) THEN
			ALLOCATE(tdata(ix,ix,nsteps,numtiles))
		ENDIF

		! allocate space for dopplergrams
		ALLOCATE(dopdata(4096,4096,loaddops)) ! allow for interp at ends
		dopdata(:,:,:) = 0.0
		! even if i have no tiles to track, still help load
		!  dops for other processors; sharing is caring.
		! loop through dopplergrams in chunks to track all tiles
		stepsdone = 0
		DO WHILE (stepsdone .LT. nsteps)
			! load some dopplergrams in parallel!
			CALL AddTime(2000)
			numdops = MIN(nsteps-stepsdone, loaddops)
			dopstart = FLOOR(numdops*1.0*myid/nproc)+1
			dopend = FLOOR(numdops*1.0*(myid+1)/nproc)
			dopstart = 1 ! remove this when the alltoallv works
			dopend = numdops ! and this
			DO ii=dopstart,dopend
				IF (.NOT.dopinterp(stepsdone+ii)) THEN
					ii2 = stepsdone+ii
					CALL Read_Dopplergram(dopfname(ii2),&
						tempdop,dop_cen_lon(ii2),dop_cen_lat(ii2),&
						dop_cen_xpix(ii2),dop_cen_ypix(ii2),dop_p_angle(ii2),&
						dop_r_sun_pix(ii2),myid,verbose)
					dopdata(:,:,ii) = tempdop(:,:)
				ENDIF
			ENDDO
			CALL AddTime(2001)
			CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

			! do an alltoallv to communicate dops
!			sendcnts(:) = INT((dopend-dopstart+1)*4096*4096)
!			sdispls(:) = INT((dopstart-1)*4096*4096) ! dont add 1?
!			DO ii=0,nproc-1
!				recvcnts(ii+1) = INT(FLOOR(numdops*1.0*(ii+2)/nproc)-&
!					FLOOR(numdops*1.0*(ii+1)/nproc))*4096*4096
!				rdispls(ii+1) = (FLOOR(numdops*1.0*(ii)/nproc))*4096*4096
!			ENDDO
			! TODO: fix alltoall, make nick do this
!			CALL MPI_Alltoallv(dopdata(1,1,1),sendcnts,sdispls,MPI_REAL,&
!				dopdata(1,1,1),recvcnts,rdispls,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL AddTime(2002)
			! do an alltoallv to communicate all accessory dop info
			sendcnts(:) = dopend-dopstart+1
			sdispls(:) = dopstart
			DO ii=0,nproc-1
				recvcnts(ii+1) = FLOOR(numdops*1.0*(ii+1)/nproc)-&
					FLOOR(numdops*1.0*ii/nproc)
				rdispls(ii+1) = FLOOR(numdops*1.0*(ii)/nproc)+1
			ENDDO
			! TODO: fix alltoall, make nick do this
!			CALL MPI_Alltoallv(dop_cen_xpix,sendcnts,sdispls,MPI_REAL,&
!			dop_cen_xpix,recvcnts,rdispls,MPI_REAL,MPI_COMM_WORLD,ierr)

			! currently have:
			!  tiles 'starttile' to 'endtile'
			!  dops 'stepsdone' to 'stepsdone+numdops'

			! track
			WRITE(*,'(A,I0,A,I0,A,I0,A)') "proc",myid," tracking ",numtiles,&
				" tiles through ",numdops," dops"

			! loop through each dopplergram
			DO ii=1,numdops
				di = ii+stepsdone
				IF (.NOT.dopinterp(di)) THEN
					delta_time = 0D0
					! loop through each tile
					DO ij=1,numtiles
						! project dop ii into tile ij
						ti = ij + starttile - 1
						currlon = lon(ti)+delta_rot(ti)*delta_time
						currlat = lat(ti)
						CALL AddTime(3000)
						CALL Projection(currlon,currlat,ix,dopdata(:,:,ii),tempslice,&
							dop_cen_xpix(di),dop_cen_ypix(di),dop_cen_lon(di),&
							dop_cen_lat(di),dop_p_angle(di),dop_r_sun_pix(di))
						CALL AddTime(3001)
						! copy result to tile
						tdata(:,:,di,ij) = tempslice
					ENDDO
				ENDIF
			ENDDO
			CALL AddTime(2003)

			stepsdone = stepsdone + numdops
		ENDDO
		
		CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

		! go back and interpolate missing frames
		DO ii=2,nsteps-1
			! check for interpolation
			IF (dopinterp(ii) .AND. .NOT. dopinterp(ii-1) .AND. &
				.NOT. dopinterp(ii-1)) THEN
				PRINT*, "interpolating step ", ii
			ENDIF
		ENDDO
		
		DEALLOCATE(tempslice)
		DEALLOCATE(dopdata)
		! save and deallocate all tiles on this proc
		IF (endtile .GE. starttile) THEN
			IF (verbose) WRITE(*,'(A,I0,A,I0,A)') "Proc ",myid," writing ",numtiles," tiles to disk.."
			DO ij=1,numtiles
				ti = ij + starttile - 1
				! make an appropriate filename TODO this
				CALL AddTime(4000)
				! need to say tile size, coordinates, carrington time
				WRITE(outfile,'(A5,I0,A1,SP,F0.4,A1,F0.4,A5)') "tile_",tilesize,"_",lon(ti),'_',lat(ti),".fits"
				IF (verbose) WRITE(*,'(A,I0,A,A)') "Proc ",myid, " saving tile ",TRIM(outfile)
				CALL Save_Tile(tdata(:,:,:,ij), ix, ix, nsteps, TRIM(outfile), verbose)
				CALL AddTime(4001)
			ENDDO
			DEALLOCATE(tdata)
		ENDIF
	END SUBROUTINE TrackTiles

	! PROJECTION
	! given a central lat/lon, tilesize, and dopplergram
	!  produce a projected square cut
	! pixels = size of tile in pixels (16 deg = 384 pix)
	! x_disk_center,y_disk_center = pixel coords of center of sun on dop image
	! p_angle = p-angle, including telescope rotation (prob around 180deg)
	! lon_disk_center,lat_disk_center = lon/lat coord of current disk center
	SUBROUTINE Projection (lon_tile_center, lat_tile_center, pixels, dop, res, x_disk_center, &
			y_disk_center,lon_disk_center, lat_disk_center, p_angle, &
			r_sun_pix)
		IMPLICIT NONE
		REAL :: lon_tile_center, lat_tile_center, dop(4096,4096)
		INTEGER :: pixels, ii, ij
		REAL :: res(pixels,pixels), val
		REAL :: xx, yy, xr, yr, x_disk_center, y_disk_center, delta_lon
		REAL :: p_angle, cos_p_angle, sin_p_angle
		REAL :: cos_lat_disk_center, sin_lat_disk_center
		REAL :: plon, plat, r, r_sun_pix, r_sun_norm, lon_disk_center, lat_disk_center
		REAL :: cos_lat_lon, cos_cang, sin_lat, cos_lat, cos_lat_tile_center, sin_lat_tile_center
		REAL :: sinlat,sinlon, sinr, cosr, x, y, sinp, cosp, coslat

		! things to possibly load?
		REAL :: mapscale = 1D0/24D0 ! degrees per pixel
		REAL :: cos_asd = 0.99998914D0 ! angular size of sun in sky
		REAL :: sin_asd = 0.004660D0
		REAL :: rad_per_deg = 0.01745329251994D0 ! convert deg to rad
		REAL :: PI = 3.1415926D0

		cos_p_angle = cos(p_angle*rad_per_deg)
		sin_p_angle = sin(p_angle*rad_per_deg)
		
		r_sun_norm = r_sun_pix * 2D0 / 4096

		PRINT*, cos_p_angle, sin_p_angle

		! pre-computations for efficiency
		! the b-angle is hidden in the lat_disk_center
		cos_lat_disk_center = cos(lat_disk_center*rad_per_deg)
		sin_lat_disk_center = sin(lat_disk_center*rad_per_deg)

		cos_lat_tile_center = cos(lat_tile_center*rad_per_deg)
		sin_lat_tile_center = sin(lat_tile_center*rad_per_deg)

		! loop over each pixel in tile
		DO ii=1,pixels
			DO ij=1,pixels
				! find pos of pixel relative to tile center in radians
				x = (ii-pixels/2.0+0.5)*mapscale*rad_per_deg
				y = (ij-pixels/2.0+0.5)*mapscale*rad_per_deg
				! polar coords (r,p) of pixel on tile
				r = sqrt(x*x + y*y)
				sinp = y / r
				cosp = x / r
				sinr = sin(r)
				cosr = cos(r)
				!
				sinlat = sin_lat_tile_center * cosr + cos_lat_tile_center * sinr * sinp
				plat = asin(sinlat)
				coslat = cos(plat)
				IF (coslat .EQ. 0.0) THEN
					sinlon = 0.0
				ELSE
					!sinlon = sinr*cosp / coslat
					sinlon = MIN(MAX(sinr*cosp / coslat,-1D0),1D0)
				ENDIF
				plon = asin(sinlon)
				! "Correction suggested by Dick Shine" - mtrack
				IF (cosr .LT. (sinlat * sin_lat_tile_center)) plon = PI-plon
				! rotate lon to tile center
				plon = plon + lon_tile_center*rad_per_deg
		
				! INTERMISSION
				! plon is the absolute heliographic angle for the current pixel
				! plat is the relative heliographic angle to tile center
				
				sin_lat = sin(plat)
				cos_lat = cos(plat)
				
				! starting projection math
				cos_lat_lon = cos_lat * cos(plon-lon_disk_center*rad_per_deg)
				cos_cang = sin_lat * sin_lat_disk_center + cos_lat_disk_center * cos_lat_lon

				! check if off disk
				IF (cos_cang < 0.0) THEN
					val = 0D0 ! missingval
				ELSE
					! r = scale factor on disk?
					r = r_sun_norm * cos_asd / (1.0 - cos_cang * sin_asd)
					xr = r * cos_lat * sin(plon-lon_disk_center*rad_per_deg)
					yr = r * (sin_lat*cos_lat_disk_center - sin_lat_disk_center*cos_lat_lon)
					! rotation matrix (camera and p-angle)
					xx = xr*cos_p_angle - yr*sin_p_angle
					yy = xr*sin_p_angle + yr*cos_p_angle
					! scale to pixels
					xx = xx * r_sun_pix
					yy = yy * r_sun_pix
					! offset from center coords
					xx = xx + x_disk_center
					yy = yy + y_disk_center
					! xx,yy in pixel coords on dopplergram
					CALL INTERP(xx,yy,dop,4096,val)
				ENDIF
				! put the interpolated (or missingval) value in the result
				res(ii,ij) = val
			ENDDO
		ENDDO

	END SUBROUTINE Projection

	! give an estimate of memory usage for a single tile
	REAL FUNCTION memusage(ts,ns)
		IMPLICIT NONE
		INTEGER :: ts, ns
		memusage = (4D0*FLOAT(ns)*FLOAT(ts*24)**2D0)/(1024D0*1024D0*1024D0)
	END FUNCTION

END MODULE Track
