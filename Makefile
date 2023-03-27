#FC=ifort 
FC= mpif90    
PROG=testcube

#ifort preprocessor flags
FPPFLAGS = "-DMPIF -DVELO" 
FFLAGS = "-c -traceback -no-wrap-margin  -heap-arrays -check bounds" 
#  -stand f90  -assume realloc_lhs  -check all  -traceback   -fstack-protector  -assume protect_parens"  #-O2 


# NETCDF library routines
# NETCDF library routines
NFDIR="/mpcdf/soft/SLE_12/packages/skylake/netcdf-mpi/intel_19.1.1-19.1.1-impi_2019.7-2019.7.217/4.4.1"
INCLUDE="-I${NFDIR}/include"
NETCDFLIB="-L../netcdf -lnet  -L${NFDIR}/lib -lnetcdff "


####################################################################

SUBDIRS = netcdf  src
 
  

all: testmain 

testmain:
	 @for i in $(SUBDIRS) ; do \
                cd $$i ; \
                $(MAKE)         \
                        FC=$(FC) \
                        EXEC=$(PROG) \
                        FFLAGS=$(FFLAGS) \
                        FPPFLAGS=$(FPPFLAGS) \
                        FFTLIB=$(FFTLIB) \
			FITSLIB=$(FITSLIB)\
                        INCLUDE=$(INCLUDE) \
                        NETCDFLIB=$(NETCDFLIB) \
                        STATIC=$(STATIC) \
                        all ; cd .. ;\
        done



clean:
	@for i in $(SUBDIRS); do \
                (cd $$i; \
                $(MAKE) clean ); \
        done
	rm -f *~

