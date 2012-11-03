FITSLOC=/shared/cfitsio
FITSLIB=$(FITSLOC)/lib
FITSINC=$(FITSLOC)/include

F90  = /usr/local/mpich2-1.0.7/bin/mpif90 -L$(FITSLIB) -lcfitsio

all : clean frack DopTime

frack : 
	$(F90) -c ParseInput.F90
	$(F90) -c Communication_Library.F90
	$(F90) -c File_Management.F90
	$(F90) -c Interpolation.F90
	$(F90) -c Timing.F90
	$(F90) -c Track.F90
	$(F90) -c Grid.F90
	$(F90) Main.F90 ParseInput.F90 Communication_Library.F90 File_Management.F90 Track.F90 Grid.F90 Timing.F90 Interpolation.F90 -o frack


DopTime : 
	$(F90) DopTime.F90 -o DopTime

clean : cleandop
	rm -f frack *.o *.mod

cleandop :
	rm DopTime
