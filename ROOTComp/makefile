ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs) -lMinuit -lGeom -lTMVA -lMLP -lEG

#COMPILER
CFLAGS		= -O -Wall -fPIC -g3 -std=c++11
CC		= g++
 
all:	pi0d_select

#Executable
pi0d_select: pi0d_select.o compactimport.o 
	$(CC) -o $@ $^ ../usersrc/funLib.c $(CFLAGS) -I../userinc -I$(ROOTINCDIR) $(ROOTLIBS) -L$(BOOST)/lib -L../obj -lexportClasses -lboost_program_options-gcc47-mt-1_53 -DOUTSIDECOMPACT=1

#objects
pi0d_select.o: pi0d_Selection.C ../userinc/funLib.h ../userinc/exportClasses.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc -I$(BOOST)/include/boost-1_53

compactimport.o: CompactImport.cpp ../userinc/exportClasses.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc

# other
funLib.o: ../usersrc/funLib.c ../userinc/funLib.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc -DOUTSIDECOMPACT=1

#Clean
.PHONY: clean cleanall

cleanall:	clean

clean: 
	rm -f pi0d_select
	rm -f *.o
