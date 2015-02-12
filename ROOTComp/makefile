ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs) -lMinuit -lGeom -lTMVA -lMLP -lEG

#COMPILER
CFLAGS		= -O -Wall -fPIC -g3 -std=c++11
CC		= g++
 
all:	pi0d_select

pi0d_select.o: pi0d_Selection.C ../userinc/funLib.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc

pi0d_select: pi0d_select.o ../obj/funLib.o
	$(CC) -o $@ $< $(CFLAGS) ../obj/funLib.o $(ROOTLIBS) -L../obj -lexportClasses

#Clean
.PHONY: clean cleanall

cleanall:	clean

clean: 
	rm -f pi0d_select
	rm -f *.o
