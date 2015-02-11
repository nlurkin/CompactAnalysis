/* File to define variables to be used in usersrc routines */
#include <stdio.h>
#include <vector>
#include "TVector3.h"
#include "mystructs.h"
#include "TTree.h"
#include "exportClasses.h"

//Output
FILE * fprt, *fprt2;
TTree *outTree;

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
bool mcOnly;
bool dataOnly;
bool optDebug;
int ffWeightType;
bool noOutput;
int outputMod;
cutsValues cutsDefinition;

//Kaon
TVector3 kaonMomentum;
double kaonP;
int beamCharge;
bool noTestCharge;
bool pbWall;
double alpha;
int period;
int periodKeep;

//Event filtering
vector<eventID> badEventsList;

//containers for corrected tracks and clusters
//std::vector<CorrectedTrack*> vtrack;
//std::vector<CorrectedCluster*> vCluster;
//std::vector<CorrectedCluster*> closeClusters;

//Container for selected corrected tracks and clusters (exported to TTree)
//std::vector<CorrectedTrack*> goodTracks;
//std::vector<CorrectedCluster*> assocClusters;
bool cutsWord[19];

//pi0dEvent fullEvent;

double Mpi0 = 0.1349766;
double Mpic = 0.139570;
//double Me = 0.000510998928;
double Me = 0.00051099891;

ROOTRawEvent rawEvent;
ROOTCorrectedEvent corrEvent;
ROOTBurst rootBurst;



