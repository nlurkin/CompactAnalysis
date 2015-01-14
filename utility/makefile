ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs) -lMinuit -lGeom -lTMVA -lMLP -lEG

#COMPILER
CFLAGS		= -O -Wall -fPIC -g3 
CC		= g++
 
all:	combine fit prepareWeights make_bins

combine.o: combine_batch.C pi0DalitzHeader.h
	@$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR)

combine: combine.o style.o root
	@$(CC) -o $@ $< style.o $(ROOTLIBS) -L../obj -lmystructs

make_bins.o: make_bins.C pi0DalitzHeader.h
	@$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR)

make_bins: make_bins.o style.o root
	@$(CC) -o $@ $< style.o $(ROOTLIBS) -L../obj -lmystructs

fit.o: fit_batch.C pi0DalitzHeader.h
	@$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR)

fit: fit.o style.o  root
	@$(CC) -o $@ $< style.o $(ROOTLIBS) -L../obj -lmystructs

prepareWeights.o: prepareWeights.C
	@$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR)

prepareWeights: prepareWeights.o root
	@$(CC) -o $@ $< $(ROOTLIBS)  -L../obj -lmystructs

style.o: style.C
	@$(CC) -o $@ -c $< $(CFLAGS) -I$(ROOTINCDIR)

root:../obj/libmystructs.so

../obj/libmystructs.so:../usersrc/mystructs.c ../obj/Dict_mystructs.cpp
	g++ -W -Wall -fPIC -g3 -shared $(ROOTCFLAGS) $(ROOTLIBS) -I$(ROOTINCDIR) ../usersrc/mystructs.c ../obj/Dict_mystructs.cpp -o ../obj/libmystructs.so -I../userinc -D__CINT_NICO__


../obj/Dict_mystructs.cpp: ../userinc/mystructs.h ../userinc/LinkDef_mystructs.h
	rootcint -v -f ../obj/Dict_mystructs.cpp -c ../userinc/mystructs.h ../userinc/LinkDef_mystructs.h

#Clean
.PHONY: clean cleanall

cleanall:	clean

clean: 
	rm -f fit
	rm -f combine
	rm -f *.o
