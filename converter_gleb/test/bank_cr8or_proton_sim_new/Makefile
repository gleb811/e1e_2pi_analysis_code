PROG     = bank_cr8or
EXE      = $(PROG)_$(OSNAME)
#EXE      = $(HOME)/bin/$(OSNAME)/$(PROG)
SOURCES  = $(wildcard *.[c] *.[f])
INC      = -I/group/clas/builds/64bit/STABLE/build/include  -I/apps/root/PRO/root/include
OBJ      = $(shell echo "$(SOURCES)" | sed -e "s/c\ /o\ /g" -e "s/f\ /o\ /g") 
#OBJ      = BOSREAD.o
CERNLIBS = -lpid -ltag -llac -lseb -lst -lclasutil -lsc -lc_bos_io -ldc -lec -lcc -ltrk -ldc -lc_bos_io -lsc -lmapmanager -lfputil -lfpack -lrecutl -lonline_dummy -lc_cern -lclasutil -lbos -lfpack -lbankdefs -L/u/home/gleb/cern/2005/lib -lpacklib -lkernlib -lnsl -lgfortran -lmathlib -lpacklib -lkernlib -lpawlib
ROOTLIBS = -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lz -pthread -lm -ldl -rdynamic


include include.mk

CXX           = g++ -Wno-write-strings -Wno-pragmas
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINCLUDE  := -I$(shell root-config --incdir)

#include include.mk

%.o:    %.cxx 
	$(COMPILER) $(INC) $(GPPOPT) -c $< -o $@

all:    $(EXE)

$(EXE):   $(OBJ) 
	$(CXX) $(INC) $(GPPOPT) -o $(EXE) $(OBJ) $(LIBS)

clean:
	rm -f  *.o core

delete:
	rm -f $(EXE)

