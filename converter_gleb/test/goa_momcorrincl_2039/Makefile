PROG     = analyser_goa
EXE      = $(PROG)_$(OSNAME)
#EXE      = $(HOME)/bin/$(OSNAME)/$(PROG)
SOURCES  = $(wildcard *.[c] *.[f])
INC      = $(wildcard *.[inc])
OBJ      = $(shell echo "$(SOURCES)" | sed -e "s/c\ /o\ /g" -e "s/f\ /o\ /g") 

CERNLIBS = -lmathlib -lphtools  -lgeant321 -lpacklib -lkernlib -lpawlib
CLASLIBS = -lc_bos_io -lrecutl  -lfputil -lbos -lfpack -lcc -lmapmanager

include include.mk
#include include.mk

%.o:    %.f $(INC)
	$(COMPILER) $(F77OPT) -c $< -o $@

all:    $(EXE)

$(EXE): $(OBJ) $(INC)
	$(COMPILER) $(F77OPT) -o $(EXE) $(OBJ) $(LIBS)

clean:
	rm -f  *.o core

delete:
	rm -f $(EXE)

