
OBJECTS=arrays.o io_cubes.o  initcalc.o rotate.o

####################################################################
## --with-cfitsio=/opt/local/cfitsio/cfitsio-3.350#


all: main 

main:    $(OBJECTS)  
	 $(FC) $(FPPFLAGS)  -o cubeconv_azimuth.exe main.F90  $(OBJECTS) $(INCLUDE) $(FITSLIB) $(NETCDFLIB)   
	mv cubeconv_azimuth.exe  ../bin 

arrays.o : arrays.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) arrays.F90 

io_cubes.o : io_cubes.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) io_cubes.F90 $(INCLUDE) $(FITSLIB) 

rotate.o : rotate.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) rotate.F90

initcalc.o : initcalc.f90
	$(FC) $(FFLAGS) initcalc.f90
 
# ========================================================================
clean:
	rm -f *~ *.o *.mod main~  
