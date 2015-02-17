/* File to define variables to be used in usersrc routines */
#include <stdio.h>
#include <vector>
#include "mystructs.h"
#include "TTree.h"
#include "exportClasses.h"

//Output
FILE * fprt, *fprt2;
TTree *outTree, *headerTree;

//CPD Cells things
float CPDpos_leftDownCorner[256][2];      	// CPDpos_leftDownCorner[256CPDs][x,y]
float CELLpos_leftDownCorner[256][64][2];  	// CELLpos_leftDownCorner[256CPDs][64Cells][x,y]
int   CPDindex, CELLindex;
float CPDlength, CELLlength;
float EopCorr[256][64];           			// correction for each cell
bool eopLoaded;

//Global variables
int iEvent;
int channel;
std::map<std::string,std::string> opts; //parsed string options
bool optDebug;
bool exportAllEvents;
bool noOutput;
int outputMod;
cutsValues cutsDefinition;
int periodKeep;

//Event filtering
vector<eventID> badEventsList;

bool mcBranched;

ROOTRawEvent rawEvent;
ROOTCorrectedEvent corrEvent;
ROOTBurst rootBurst;
ROOTFileHeader rootFileHeader;
NGeom rootGeom;
ROOTMCEvent rootMC;

