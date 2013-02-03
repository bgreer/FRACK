
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
	dopfname, doptime, dopfname_ends, doptime_ends, dopinterp,loaddops,outdir,backframe,verbose)
		IMPLICIT NONE
		LOGICAL :: verbose
		INTEGER :: myid, nproc, num, tilesize, startindex,nsteps,ierr
		REAL, ALLOCATABLE :: tdata(:,:,:,:)
		INTEGER :: starttile, endtile, ix, numtiles, ii,ij, stepsdone
		CHARACTER(LEN=200) :: outfile, buff
		CHARACTER(LEN=400) :: outdir
		CHARACTER(LEN=200) :: dopfname(nsteps), dopfname_ends(2)
		INTEGER :: doptime(nsteps), doptime_ends(2)
		LOGICAL :: dopinterp(nsteps)
		INTEGER :: loaddops, numdops, dopstart, dopend
		REAL, ALLOCATABLE :: dopdata(:,:,:), dopdata_ends(:,:,:), tempslice(:,:)
		REAL :: tempdop(4096,4096), backframe(4096,4096)
		INTEGER :: sendcnts(nproc), sdispls(nproc), otherstart, otherend
		INTEGER :: recvcnts(nproc), rdispls(nproc)
		INTEGER :: ti, di, ii2, interp1, interp2
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
!			dopstart = 1 ! remove this when the alltoallv works
!			dopend = numdops ! and this
			DO ii=dopstart,dopend
				IF (.NOT.dopinterp(stepsdone+ii)) THEN
					ii2 = stepsdone+ii
					CALL Read_Dopplergram(dopfname(ii2),&
						tempdop,dop_cen_lon(ii2),dop_cen_lat(ii2),&
						dop_cen_xpix(ii2),dop_cen_ypix(ii2),dop_p_angle(ii2),&
						dop_r_sun_pix(ii2),myid,verbose)
					dopdata(:,:,ii) = tempdop(:,:)
					! subtract background here
					dopdata(:,:,ii) = dopdata(:,:,ii) - backframe
				ENDIF
			ENDDO
			CALL AddTime(2001)
			CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

			! do an alltoallv to communicate dops
			IF (dopend .GE. dopstart) THEN
				sendcnts(:) = (dopend-dopstart+1)
				sdispls(:) = (dopstart-1)
			ELSE
				sendcnts(:) = 0
				sdispls(:) = 0
			ENDIF
			DO ii=0,nproc-1
				otherstart = FLOOR(numdops*1.0*ii/nproc)+1
				otherend = FLOOR(numdops*1.0*(ii+1)/nproc)
				recvcnts(ii+1) = (otherend-otherstart+1)
				rdispls(ii+1) = (otherstart-1)
				IF (otherend .LT. otherstart) THEN
					recvcnts(ii+1) = 0
					rdispls(ii+1) = 0
				ENDIF
			ENDDO

			! Share dopplergram data as well as header info
			CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
			IF (verbose .AND. myid .EQ. 0) WRITE(*,'(A)') "Entering Alltoallv.."
			CALL MPI_Alltoallv(dopdata,sendcnts*4096*4096,sdispls*4096*4096,MPI_REAL,&
				dopdata,recvcnts*4096*4096,rdispls*4096*4096,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_cen_lon,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_cen_lon,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_cen_lat,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_cen_lat,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_cen_xpix,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_cen_xpix,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_cen_ypix,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_cen_ypix,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_p_angle,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_p_angle,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			CALL MPI_Alltoallv(dop_r_sun_pix,sendcnts,sdispls+stepsdone,MPI_REAL,&
				dop_r_sun_pix,recvcnts,rdispls+stepsdone,MPI_REAL,MPI_COMM_WORLD,ierr)
			IF (verbose .AND. myid .EQ. 0) WRITE(*,'(A)') "Finished Alltoallv."
!			PRINT*, myid, dopdata(2000,2000,:)

			CALL AddTime(2002)

			CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

			! currently have:
			!  tiles 'starttile' to 'endtile'
			!  dops 'stepsdone' to 'stepsdone+numdops'

			! track
			WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "Proc",myid," tracking ",numtiles,&
				" tiles through dops ",(stepsdone)," to ",(numdops+stepsdone),&
				" : ", numdops," total"

			! loop through each dopplergram
			DO ii=1,numdops
				di = ii+stepsdone
				IF (.NOT.dopinterp(di)) THEN
					delta_time = 45D0*(di-0.5*nsteps)
					! loop through each tile
					DO ij=1,numtiles
!						PRINT*, myid, "projecting ", ii, ij
						! project dop ii into tile ij
						ti = ij + starttile - 1
						currlon = lon(ti)+delta_rot(ti)*delta_time
						currlat = lat(ti)
						CALL AddTime(3000)
!						PRINT*, myid, "calling Projection..", dop_cen_xpix(di)
						CALL Projection(currlon,currlat,ix,dopdata(:,:,ii),tempslice,&
							dop_cen_xpix(di),dop_cen_ypix(di),dop_cen_lon(di),&
							dop_cen_lat(di),dop_p_angle(di),dop_r_sun_pix(di))
						CALL AddTime(3001)
!						PRINT*, myid, "done projecting, copying result.."
						! copy result to tile
						tdata(:,:,di,ij) = tempslice
!						PRINT*, myid, "done copying"
					ENDDO
				ENDIF
			ENDDO
!			PRINT*, myid, "done projecting, barrier"
			CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
!			PRINT*, myid, "done"
			CALL AddTime(2003)

			stepsdone = stepsdone + numdops
		ENDDO
		
!		PRINT*, myid, "reaches barrier"
		CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

		! go back and interpolate missing frames
		DO ii=2,nsteps-1
			! check for interpolation
			IF (dopinterp(ii)) THEN
				! find nearest non-interpolated neighbors
				interp1 = ii-1
				DO WHILE (interp1 .GT. 1 .AND. dopinterp(interp1))
					interp1 = interp1 - 1
				ENDDO
				interp2 = ii+1
				DO WHILE (interp2 .LT. nsteps .AND. dopinterp(interp2))
					interp2 = interp2 + 1
				ENDDO
				PRINT*, "interpolating step ", ii, " from ", interp1, " and ", interp2
				! make each pixel a linear combination
				! loop through each local tile
				DO ij=1,numtiles
					tdata(:,:,ii,ij) = ((tdata(:,:,interp2,ij)-tdata(:,:,interp1,ij))*(ii-interp1)) &
						/ FLOAT(interp2 - interp1) + tdata(:,:,interp1,ij)
				ENDDO
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
				WRITE(outfile,'(A,A5,I0,A1,SP,F0.4,A1,F0.4,A5)') TRIM(outdir),"/tile_",tilesize,"_",lon(ti),'_',lat(ti),".fits"
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
		REAL :: sinlat,sinlon, sinr, cosr, x, y, sinp, cosp, coslat, vobs

		! things to possibly load?
		REAL :: mapscale = 1D0/24D0 ! degrees per pixel
		REAL :: cos_asd = 0.99998914D0 ! angular size of sun in sky
		REAL :: sin_asd = 0.004660D0
		REAL :: rad_per_deg = 0.01745329251994D0 ! convert deg to rad
		REAL :: PI = 3.1415926D0

		cos_p_angle = cos(p_angle*rad_per_deg)
		sin_p_angle = sin(p_angle*rad_per_deg)
		
		r_sun_norm = r_sun_pix * 2D0 / 4096D0

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
					r = r_sun_norm * cos_asd / (1.D0 - cos_cang * sin_asd)
					xr = r * cos_lat * sin(plon-lon_disk_center*rad_per_deg)
					yr = r * (sin_lat*cos_lat_disk_center - sin_lat_disk_center*cos_lat_lon)
					! rotation matrix (camera and p-angle)
					xx = xr*cos_p_angle - yr*sin_p_angle
					yy = xr*sin_p_angle + yr*cos_p_angle
					! scale to pixels
!					xx = xx * r_sun_pix
!					yy = yy * r_sun_pix
					xx = xx * 2048D0
					yy = yy * 2048D0
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

		! compute LOS correction
		CALL CorrectLOS(vobs)
		! apply
		res(:,:) = res(:,:) - vobs

	END SUBROUTINE Projection

	SUBROUTINE CorrectLOS(vobs)
		IMPLICIT NONE
		REAL :: vobs, sig, chi
		REAL :: orbvr, orbvw, orbvn
		REAL :: x, y, tanang_r

		x = coslat * sin(clon)
		y = sin (mapclat[m]) * coslatc - sinlatc * coslat * cos(clon)
		chi = atan2(x, y)
		sig = atan(sqrt(x*x + y*y) * tanang_r)

		vobs = orbvr * cos(sig)
		vobs = vobs - orbvw * sin(sig) * sin(chi)
		vobs = vobs - orbvn * sin(sig) * cos(chi)

!  double chi, clon, coslat, sig, x, y;
!  double coslatc = cos (latc), sinlatc = sin (latc);
!  int m, n, mset = 0;

!  double tanang_r = tan (semidiam);
!				       /*  No correction for foreshortening  */
!  for (m = 0; m < mapct; m++) {
!    coslat = cos (mapclat[m]);
!    clon = mapclon[m] - lonc;
!    x = coslat * sin (clon);
!    y = sin (mapclat[m]) * coslatc - sinlatc * coslat * cos (clon);
!  }

	END SUBROUTINE CorrectLOS

	! give an estimate of memory usage for a single tile
	REAL FUNCTION memusage(ts,ns)
		IMPLICIT NONE
		INTEGER :: ts, ns
		memusage = (4D0*FLOAT(ns)*FLOAT(ts*24)**2D0)/(1024D0*1024D0*1024D0)
	END FUNCTION

END MODULE Track
