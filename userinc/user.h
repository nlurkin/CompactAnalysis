#include "mystructs.h"
#include "TTree.h"
#include "ana3pi.h"
#include "exportClasses.h"

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
extern TTree *outTree, *headerTree;			// Output TTree

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
extern std::map<std::string,std::string> opts;			// Parsed options map
extern int channel;						// Selected channel
extern bool optDebug;					// Activate debugging
extern bool exportAllEvents;			// Export all events. If false, export only events passing the cuts
extern bool noOutput;					// Do not create output txt files
extern int outputMod;					// Output events index every outputMod
extern cutsValues cutsDefinition;		// Cuts values
extern int periodKeep;					// Period number to keep

/*
 * Kaon
 */
extern abcog_params_t abcog_params;		// Definition of abcog_params for eclipse indexer

/*
 * Event filtering
 */
extern std::vector<eventID> badEventsList;	// List of events extracted

extern bool mcBranched;

extern ROOTRawEvent rawEvent;
extern ROOTCorrectedEvent corrEvent;
extern ROOTBurst rootBurst;
extern ROOTFileHeader rootFileHeader;
extern NGeom rootGeom;
extern ROOTMCEvent rootMC;

/*
 * Function declarations
 */
//compact
extern "C" int 	LKr_acc			(int,float,float,float);

//Channel selects
int nico_ke2Select(superBurst *sbur,superCmpEvent *sevt);
int nico_pi0DalitzSelect();

//Channel analysis
int nico_pi0DalitzAna(superBurst *sbur,superCmpEvent *sevt);
int nico_pi0DalitzFilter(superBurst *sbur,superCmpEvent *sevt);



#endif
