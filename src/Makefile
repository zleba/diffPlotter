###########################################################################################
#
#  simple Makefile for fastNLO for Diffraction
#
###########################################################################################

# please set your fastNLO prefix path here, or use fnlo-tk-config
#FastNLOInstallDir=$(fnlo-tk-config --prefix)
FastNLOInstallDir=../fastnloBuild/

ROOT_LIBS=$(shell root-config --libs)
ROOT_INCLUDE=$(shell root-config --cflags)

LHA_LIBS=$(shell lhapdf-config --libdir)
LHA_INCLUDE=$(shell lhapdf-config --cflags)

ARMA_DIR=../arma/install

DPDF_LIBS=../dpdfLib/


# flags and variables
OBJS 	= fnlo-tk-h1diffpdf.o 
LDFLAGS	= -lm -lpthread  -lgfortran  -lstdc++   # `gtk-config --cflags` `gtk-config --libs` 
CC 	= g++
CFLAGS 	= -Wall -g -O0 -I$(FastNLOInstallDir)/include 
EXEC	= fastNLO


fastCalc:  obj/diffCalc.o obj/fastDiffLib.o obj/fastNLODiffAlphas.o obj/diffPlotter.o
	$(CC) -g  -fopenmp $^  $(LDFLAGS) \
	 -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(ROOT_LIBS)    \
	-L$(DPDF_LIBS) -lDPDFset \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-L$(LHA_LIBS)   -lLHAPDF \
	-Wl,-rpath  $(LHA_LIBS)  \
	-Wl,-rpath  $(DPDF_LIBS)  \
	-Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@





fastNLO: $(OBJS) fastNLODiffAlphas.h
	$(CC) -g  $(OBJS) $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(LHA_INCLUDE)    -I$(DPDF_LIBS)  \
	-L$(LHA_LIBS) -lLHAPDF \
	-L$(DPDF_LIBS) -lDPDFset \
	-Wl,-rpath $(LHA_LIBS) \
    -Wl,-rpath $(DPDF_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@


testFast: testFast.o fastNLODiffAlphas.h
	$(CC) -g -fopenmp  testFast.o $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(LHA_INCLUDE)    -I$(DPDF_LIBS)  \
	-L$(LHA_LIBS) -lLHAPDF \
	-L$(DPDF_LIBS) -lDPDFset \
	-Wl,-rpath $(LHA_LIBS) \
    -Wl,-rpath $(DPDF_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@




plotter: plotter.o diffpdferr.f 
	$(CC)   $^  $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(ROOT_LIBS)    \
	-Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@

# ../data/Measurement/fullErrors/readcov.o
#	-I../../../Helenka/dPDFs/install/include/    -I../../../Helenka/dPDFs/dpdfLib  
plotterNew: plotterNew.o  ../data/Measurement/fullErrors/readcov.o
	$(CC)  -fopenmp $^  $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include   -I$(DPDF_LIBS) -I$(LHA_INCLUDE)  -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(ROOT_LIBS)    \
	-L$(ARMA_DIR)/lib64 -larmadillo    \
	-L$(DPDF_LIBS) -lDPDFset \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-L$(LHA_LIBS)   -lLHAPDF \
	-Wl,-rpath  $(LHA_LIBS)  \
	-Wl,-rpath  $(DPDF_LIBS)  \
	-Wl,-rpath $(FastNLOInstallDir)/lib \
	-Wl,-rpath ../arma/install/lib64 \
	-Wl,-rpath /opt/anaconda/3/pkgs/mkl-2019.1-144/lib/    \
	-o $@




clean:
	rm -f obj/*.o fastDiff



clean-objs:
	rm -f *.o *.a *.so *Dict* 2>/dev/null


obj/%.o: src/%.cc 
	$(CC) -g -fopenmp   $(CFLAGS) -I$(DPDF_LIBS) -Iinc  -c $(ROOT_INCLUDE) $(LHA_INCLUDE) -I../PlottingHelper/      -o $@ $^   


all: fnlo-tk-h1diffpdf clean-objs


