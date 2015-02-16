ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs) -lMinuit -lGeom -lTMVA -lMLP -lEG

#COMPILER
CFLAGS		= -O -Wall -fPIC -g3 -std=c++11
CC		= g++
 
all:	pi0d_select

pi0d_select.o: pi0d_Selection.C ../userinc/funLib.h ../userinc/exportClasses.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc

pi0d_select: pi0d_select.o
	$(CC) -o $@ $< ../usersrc/funLib.c $(CFLAGS) -I../userinc -I$(ROOTINCDIR) $(ROOTLIBS) -L../obj -lexportClasses -DOUTSIDECOMPACT=1

funLib.o: ../usersrc/funLib.c ../userinc/funLib.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR) -I../userinc -DOUTSIDECOMPACT=1

#Clean
.PHONY: clean cleanall

cleanall:	clean

clean: 
	rm -f pi0d_select
	rm -f *.o
