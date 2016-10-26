#
# if DEBUG then make exe in current directory
#
ifdef DEBUG
     EXE    = $(PROG)_$(OSNAME)_debug
endif




  COMPILER = g++
  GPPOPT   = -O2 -Wall -fPIC  -D_compile_g11pcor  -pthread


LIBS     = -L/usr/lib64 -L/home/fedotov/GLEB -L/u/home/gleb/lib/LinuxRHFC8 $(CERNLIBS)  -L/apps/root/PRO/root/lib $(ROOTLIBS)

