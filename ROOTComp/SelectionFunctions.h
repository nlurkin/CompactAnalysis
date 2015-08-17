/*
 * SelectionFunctions.h
 *
 *  Created on: Aug 13, 2015
 *      Author: nlurkin
 */

#ifndef SELECTIONFUNCTIONS_H_
#define SELECTIONFUNCTIONS_H_

#include <TH2D.h>

// Local includes
#include "CompactIO.h"
#include "OptionsParser.h"

extern CompactIO io;
extern OptionsParser options;

///### TTree objects
extern ROOTRawEvent &rawEvent;
extern ROOTCorrectedEvent &corrEvent;
extern ROOTBurst &rootBurst;
extern ROOTFileHeader &rootFileHeader;
extern NGeom &rootGeom;
extern ROOTMCEvent &rootMC;
extern ROOTPhysicsEvent &rootPhysics;
extern ROOTFileHeader &outputFileHeader;

// ### Database objects
extern NAbcog_params abcog_params;

// ### K2pi kmu3 selectors
//extern double Mx;
//extern int xPDGId;
//extern NRecoParticle *xRootPhysics;

// ### Histograms
extern TH1D meeTrue;
extern TH1D mkTrue;
extern TH2D meeexTrue;
extern TH2D meekTrue;

extern TH1D meeFalse;
extern TH1D mkFalse;
extern TH2D meeexFalse;
extern TH2D meekFalse;

extern TH1D meeTotal;
extern TH1D mkTotal;
extern TH2D meeexTotal;
extern TH2D meekTotal;

extern TH1D meeDiffTrue;
extern TH1D mkDiffTrue;

extern TH1D meeDiffFalse;
extern TH1D mkDiffFalse ;

extern TH1D nxCandidates;
extern TH1D nxCandidatesNew;

extern TH1D nBeamSign;

extern TH2D xTruexFalse;
extern TH2D xTruexFalseMany;
extern TH2D xTruexFalseNo;
extern TH1D xMCNoID;
extern TH1D xMCManyID;

extern TH2D xTruexMCMany;
extern TH2D xTruexMCNo;

extern TH2D combi2Dpi;
extern TH2D combi2Dk;
extern TH2D combi2Dk_exclu;

extern TH1D eopx;
extern TH1D eope;

extern TH1D eope_goodx;
extern TH1D eope_badx;
extern TH1D eopx_bade;
extern TH1D eopx_goode;

extern TH1D eoplowest;
extern TH1D eopsecond;
extern TH1D eophighest;

extern int em;
extern int ep;
extern int xPart;
extern bool flBad;
extern double xTrue;
extern double xFalse;

//extern struct alt_pid_res pid_res;

//extern bool isK2PiType;
//extern bool isKMu3Type;

#define PRINTVAR(v) #v << "= " << v << " "

int pid(int &xCandidate, TLorentzVector &gamma, OptionsParser::ESelectionType t);
int pi0d_tracksAcceptance();
int pi0d_trackCombinationVeto_loose();
int pi0d_trackCombinationVeto_tight(NRecoParticle &xParticle);
int pi0d_goodClusters_loose();
int pi0d_goodClusters_tight(NRecoParticle &xParticle, ROOTPhysicsEvent &event);
bool pi0d_failCut(int i);
bool pi0d_failCutInc(int i, bool assoc, bool good, bool bad, struct alt_pid_res *pid_res);
void pi0d_passSelection();
vector<float> sortEOP();
int pi0d_identifyPiEOP(int &piCandidate, bool &badElectron);
bool associateMCTracks(struct alt_pid_res &pid_res);

void savePlots();

#endif /* SELECTIONFUNCTIONS_H_ */
