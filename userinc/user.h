#include <vector>
#include "TVector3.h"
#include <sstream>
#include "mystructs.h"
#include "TTree.h"
#include <iostream>
#include "ana3pi.h"
#include "TH1D.h"
#include "exportClasses.h"
using namespace std;

#ifndef USER_HH
#define USER_HH

/*
 * Channels defines
 */
#define KE2 1
#define PI0DALITZ 2
#define NONE 3

/*
 * Output
 */
extern FILE *fprt, *fprt2;		// Pointers to output files (passing events, failing events)
extern TTree *outTree;			// Output TTree

/*
 * CPD Cells things
 */
extern float CPDpos_leftDownCorner[256][2];       	// CPDpos_leftDownCorner[256CPDs][x,y]
extern float CELLpos_leftDownCorner[256][64][2];  	// CELLpos_leftDownCorner[256CPDs][64Cells][x,y]
extern int   CPDindex, CELLindex;					//
extern float CPDlength, CELLlength;					//
extern float EopCorr[256][64];           			// correction for each cell
extern bool eopLoaded;								// Did we already load the eop data for the current run?

/*
 * Global variables
 */
extern int iEvent;						// Current event index
extern char gString[200];				// Global string passed via CLI
extern map<string,string> opts;			// Parsed options map
extern int channel;						// Selected channel
extern bool mcOnly;						// For MC only actions
extern bool dataOnly;					// For Data only actions
extern bool optDebug;					// Activate debugging
extern int ffWeightType;				// Type of form factor weight (MC only)
										// 0:FF=1 , 1:FF=x, 2:FF=x^2
extern bool noOutput;					// Do not create output txt files
extern int outputMod;					// Output events index every outputMod
extern cutsValues cutsDefinition;		// Cuts values

/*
 * Kaon
 */
extern TVector3 kaonMomentum;			// Kaon unit momentum
extern double kaonP;					// Kaon momentum magnitude
extern int beamCharge;					// Kaon charge
extern bool noTestCharge;				// Beam charge not defined
extern abcog_params_t abcog_params;		// Definition of abcog_params for eclipse indexer
extern bool pbWall;						// Is the Pb wall present
extern double alpha;					// alpha correction to mc kaon weight
extern int period;						// Period number
extern int periodKeep;					// Period number to keep
/*
 * Event filtering
 */
extern vector<eventID> badEventsList;	// List of events extracted

/*
 * Containers for corrected tracks and clusters
 */
//extern vector<CorrectedTrack*> vtrack;
//extern vector<CorrectedCluster*> vCluster;
//extern vector<CorrectedCluster*> closeClusters;

/*
 * Container for selected corrected tracks and clusters (exported to TTree)
 */
//extern vector<CorrectedTrack*> goodTracks;
//extern vector<CorrectedCluster*> assocClusters;
extern bool cutsWord[19];

//extern pi0dEvent fullEvent;

extern double Mpi0;
extern double Mpic;
extern double Me;

extern ROOTRawEvent rawEvent;
extern ROOTCorrectedEvent corrEvent;
extern ROOTBurst rootBurst;


/*
 * Function declarations
 */
//compact
extern "C" int 	LKr_acc			(int,float,float,float);
void 			GetCpdCellIndex	(double pos_x, double pos_y, int *cpd_index, int *cell_index);
float 			correctedEP		(superCmpEvent* sevt, trak t, float &eOverP);
float 			correctClusterE	(NPhysicsCluster *c);
void 			loadEOPData		(superBurst *sbur);
void 			defineBeamCharge(superBurst *sbur);

//TVector related
string 		printVector3		(TVector3 v);
template<typename T>
string printSTLvector(vector<T> v){
	ostringstream ss;
	typename vector<T>::iterator it;

	ss << "Vector values :" << endl;
	for(it=v.begin(); it!=v.end();it++){
		ss << "\t" << (*it);
		ss << endl;
	}
	ss << endl;
	return ss.str();
}
string printSTLvector(vector<TVector3> v);

TH1D* gH(TString name);

double		distance2D			(TVector3 v1, TVector3 v2);
void 		propagateBefore		(float &x, float &y, float &z, float zplane, trak t);
TVector3 	propagateBefore		(float zplane, NPhysicsTrack t);
TVector3 	propagateCorrBefore	(float zplane, NPhysicsTrack t);
void 		propagateAfter		(float &x, float &y, float &z, float zplane, trak t);
TVector3 	propagateAfter		(float zplane, NPhysicsTrack t);
TVector3 	propagate			(float zplane, NPhysicsTrack t);
TVector3 	propagate			(float zplane, TVector3 pos, TVector3 p);

double 	missMass2	(double m1, double m2, TVector3 p1, TVector3 p2);
double 	missMass2Man(double m1, double m2, float kp, float kdxdz, float kdydz, float ep, float edxdz, float edydz);
double	invMass2	(vector<double> mass, vector<TVector3> p);

bool 	isFilteredEvent	(int nrun, int nburst, int timestamp);

//Corrections
NPhysicsTrack	correctTrack	(superCmpEvent *sevt, trak t);
void			CreateTracks	(superCmpEvent *sevt);
NPhysicsCluster	correctCluster	(cluster c);
void		 	CreateClusters	(superCmpEvent *sevt);

//Channel selects
int nico_ke2Select(superBurst *sbur,superCmpEvent *sevt);
int nico_pi0DalitzSelect();

//Channel analysis
int nico_pi0DalitzAna(superBurst *sbur,superCmpEvent *sevt);
int nico_pi0DalitzFilter(superBurst *sbur,superCmpEvent *sevt);



#endif
