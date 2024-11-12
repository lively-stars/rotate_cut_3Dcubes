FC= mpiifort
PROG=testcube


#ifort preprocessor flags
FPPFLAGS = "-DMPIF -DVELO " 
FFLAGS = "-c -traceback -no-wrap-margin  -heap-arrays -check bounds" 
#  -stand f90  -assume realloc_lhs  -check all  -traceback   -fstack-protector  -assume protect_parens"  #-O2 


# NETCDF library routines
# NETCDF library routines
NFDIR="/mpcdf/soft/SLE_15/packages/skylake/netcdf-mpi/intel_19.1.3-19.1.3-impi_2019.9-2019.9.304/4.4.1/"
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

