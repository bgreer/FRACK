
! Just the cubic convolution interpolation routine

MODULE Interpolation
	IMPLICIT NONE
CONTAINS

	! ASSUME: x,y in pixel coords
	!  so x=[0,dimx]
	!  also, arr is square
	! IN: x,y,arr,dimx
	! OUT: val
	SUBROUTINE INTERP(x,y,arr,dimx,val)
		IMPLICIT NONE
		INTEGER :: ii, ij, dimx
		REAL :: x, y, val, weight, temp
		REAL :: arr(dimx,dimx)
		! Coordinates of surrounding points
		INTEGER, DIMENSION(4) :: ix, iy
		! Imaginary data points at lat=y
		REAL, DIMENSION(4) :: rx
		
		! For bicubic convolution interpolation, need nearest 2 points
		!  in each direction. Find them in O(1) time.
		ix(2) = MAX(FLOOR(x), 1)
		iy(2) = MAX(FLOOR(y), 1)
		ix(3) = MIN(ix(2)+1, dimx)
		iy(3) = MIN(iy(2)+1, dimx)
		ix(1) = MAX(ix(2)-1,1)
		iy(1) = MAX(iy(2)-1,1)
		ix(4) = MIN(ix(2)+2,dimx)
		iy(4) = MIN(iy(2)+2,dimx)

		! Interpolate at current y-val
		rx = 0.0
		DO ii=1,4 ! x-index
			weight = 0.0
			DO ij=1,4 ! sum over y
				temp = KERNEL(iy(ij)-y)
				weight = weight + temp
				rx(ii) = rx(ii) + temp*arr(ix(ii),iy(ij))
			ENDDO
			rx(ii) = rx(ii) / weight
		ENDDO

		! Interpolate at current x-val
		val = 0.0
		weight = 0.0
		DO ii=1,4
			temp = KERNEL(ix(ii)-x)
			weight = weight + temp
			val = val + temp*rx(ii)
		ENDDO

		val = val / weight
		RETURN
	END SUBROUTINE INTERP

	REAL FUNCTION KERNEL (x)
		IMPLICIT NONE
		REAL :: x, res
		x = ABS(x)
		IF (x .le. 1.0) THEN
			res = 1.0 + x*x*(-2.5 + 1.5*x)
		ELSE IF (x .lt. 2.0) THEN
			res = 2.0 + x*(-4.0 + x*(2.5 - 0.5*x))
		ELSE
			res = 0.0
		ENDIF
		KERNEL = res
	END FUNCTION KERNEL
END MODULE Interpolation
