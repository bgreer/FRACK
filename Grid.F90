
MODULE Grid
	IMPLICIT NONE

CONTAINS

	! Read a grid file, each line is a tile, cols are:
	!  [(int)tilesize]  [(real)lon]  [(real)lat]
	SUBROUTINE ReadGrid ()
		IMPLICIT NONE
		! TODO implement reading grid file, someday
	END SUBROUTINE ReadGrid

	! Generate central tile coords for all tiles in a region
	! INPUT: ts, lonrn, latrn, clon, clat, apode
	! OUTPUT: lon(), lat()  assumed to be allocated with correct size
	! Only one kind of grid supported, dense-pack
	SUBROUTINE GenerateGrid (ts, lonrn, latrn, clon, clat, apode, lon, lat)
		IMPLICIT NONE

		INTEGER :: ts, numx, numy, ii, ij, currtile
		REAL :: lonrn, latrn, clon, clat, apode, space
		REAL :: lonmin, latmin
		REAL :: lon(:), lat(:)

		space = ts*0.5*apode

		numx = INT(lonrn/space+1)
		numy = INT(latrn/space+1)

		lonmin = clon - 0.5*lonrn
		latmin = clat - 0.5*latrn

		currtile = 1
		DO ii=1,numx
			DO ij=1,numy
				lon(currtile) = lonmin+(ii-1)*space
				! Wrap longitude around
				IF (lon(currtile) .GT. 360D0) &
					lon(currtile) = lon(currtile) - 360D0
				IF (lon(currtile) .LT. 0D0) &
					lon(currtile) = lon(currtile) + 360D0
				lat(currtile) = latmin+(ij-1)*space
				currtile = currtile + 1
			ENDDO
		ENDDO

	END SUBROUTINE GenerateGrid

	! Set delta_rot for each tile
	! INPUT: a0,a2,a4, lat(num), num
	! OUTPUT: delta_rot(num)
	SUBROUTINE SetTrackingRate (delta_rot, a0, a2, a4, lat, num)
		IMPLICIT NONE
		INTEGER :: num, ii
		REAL :: delta_rot(num), lat(num)
		REAL :: sinlat, a0, a2, a4

		DO ii=1,num
			sinlat = sin(lat(ii))
			sinlat = sinlat * sinlat
			! units of radians per second
			delta_rot(ii) = (a0 + sinlat * (a2 + a4*sinlat))*1E-6
		ENDDO
	END SUBROUTINE SetTrackingRate

	! Given a lat/lon range and tilesize, give number of tiles
	INTEGER FUNCTION GetTileCount (ts, lonrn, latrn, apode)
		IMPLICIT NONE
		INTEGER :: ts
		REAL :: lonrn, latrn, space, apode
		space = ts*0.5*apode
		GetTileCount = INT(lonrn/space+1)*INT(latrn/space+1)
	END FUNCTION
END MODULE
