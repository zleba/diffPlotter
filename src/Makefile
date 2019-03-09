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

DPDF_LIBS=../dpdfLib/


# flags and variables
OBJS 	= fnlo-tk-h1diffpdf.o 
LDFLAGS	= -lm -lpthread  -lgfortran  -lstdc++   # `gtk-config --cflags` `gtk-config --libs` 
CC 	= g++
CFLAGS 	= -Wall -g -O2 -I$(FastNLOInstallDir)/include 
EXEC	= fastNLO





fastNLO: $(OBJS) fastNLODiffAlphas.h
	$(CC) -g  $(OBJS) $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(LHA_INCLUDE)    -I$(DPDF_LIBS)  \
	-L$(LHA_LIBS) -lLHAPDF \
	-L$(DPDF_LIBS) -lDPDFset \
	-Wl,-rpath $(LHA_LIBS) \
    -Wl,-rpath $(DPDF_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $(EXEC)

plotter: plotter.o diffpdferr.f 
	$(CC)   $^  $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	$(ROOT_LIBS)    \
	-Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@

plotterNew: plotterNew.o diffpdferr.f fastNLODiffAlphas.h ../data/Measurement/fullErrors/readcov.o
	$(CC)   $^  $(LDFLAGS) \
	-I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-I../../../Helenka/dPDFs/install/include/    -I../../../Helenka/dPDFs/dpdfLib  \
	$(ROOT_LIBS)    \
	-larmadillo    \
	-L../../../Helenka/dPDFs/dpdfLib -lDPDFset \
	h1fitJets/h12007jetsdpdf_pdfMy.o h1fitJets/h12007jetsdpdf_charm.o h1fitJets/h12007jetsdpdf_gluon.o \
	h1fitJets/h12007jetsdpdf_prcoeff.o h1fitJets/h12007jetsdpdf_singlet.o zeuspdf/lib/libzdpdf09.a \
	MRW2006/fint.o MRW2006/mrw2006.o \
	-L/home/radek/Dropbox/patrick/undolding/PlottingHelper -lPlottingHelper -Wl,-rpath,/home/radek/Dropbox/patrick/undolding/PlottingHelper   \
	-L../../../Helenka/dPDFs/install/lib   -lLHAPDF \
	-Wl,-rpath  ../../../Helenka/dPDFs/install/lib  \
	-Wl,-rpath  ../../../Helenka/dPDFs/dpdfLib  \
	-Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@




clean: clean-objs
	rm $(EXEC) 2>/dev/null


clean-objs:
	rm -f *.o *.a *.so *Dict* 2>/dev/null


%.o: %.cc 
	$(CC)  $(CFLAGS)  -c $(ROOT_INCLUDE) $(LHA_INCLUDE) -I$(DPDF_LIBS)     -o $@ $^   


all: fnlo-tk-h1diffpdf clean-objs


