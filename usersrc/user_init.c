/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <sstream>
#include <iostream>
#include <map>
#include <iomanip>
#include "funLib.h"
#include "compactLib.h"
using namespace std;

int nico_ke2Init(){
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	//goodTracks.clear();
	//assocClusters.clear();
	outTree = new TTree("event", "Event");
	//outTree->Branch("goodTrack", "std::vector<CorrectedTrack>", &goodTracks, 64000, 1);
	//outTree->Branch("assocCluster", "std::vector<CorrectedCluster>", &assocClusters, 64000, 1);

	new TH1D("Cuts", "Failed Cuts", 20, 0.5, 20.5);
	new TH1D("CutsKe2", "Failed Ke2 Cuts", 20, 0.5, 20.5);

	//Histograms for common cuts
	new TH1D("timeDchOffset", "timeDchOffset", 100, 0, 100);
	new TH1D("deadCellDist", "deadCellDist", 100, 0, 15);
	new TH2D("DCHAcceptance1", "DCHAcceptance1", 400, -200, 200, 400, -200, 200);
	new TH2D("DCHAcceptance4", "DCHAcceptance4", 400, -200, 200, 400, -200, 200);
	new TH1D("RDCH1", "RDCH1", 120, 0, 120);
	new TH1D("RDCH4", "RDCH4", 140, 0, 140);
	new TH2D("excludedZone", "excludedZone", 300, -300, 300, 600, -300, 300);
	new TH2D("hotCellDistT", "hotCellDistT", 400, -200, 200, 400, -200, 200);
	new TH2D("hotCellDistC", "hotCellDistC", 400, -200, 200, 400, -200, 200);
	new TH2D("lkrAcceptance", "lkrAcceptance", 400, -200, 200, 400, -200, 200);
	new TH1D("cda", "cda", 200, 0, 20);

	new TH1D("trackMomentum", "trackMomentum", 500, 0, 100);
	new TH1D("trackQuality", "trackQuality", 200, 0, 2);

	//Histogram for ke2 Cuts
	new TH1D("mmass2", "mmass2", 400, -1, 1);
	new TH1D("clusterStatus", "clusterStatus", 10, 0, 10);
	new TH1D("trckClDist", "trckClDist", 100, 0, 10);
	new TH1D("trckClTime", "trckClTime", 300, 0, 30);
	new TH1D("eop", "eop", 100, 0, 2);

	//Final histos
	new TH1D("fmmass2", "fmmass2", 400, -1, 1);

	return 0;
}

int nico_pi0DalitzInit(){
	rawEvent.clear();
	corrEvent.clear();
	outTree = new TTree("event", "Event");
	headerTree = new TTree("header", "Header");
	//outTree->Branch("goodTrack", "std::vector<CorrectedTrack*>", &goodTracks, 64000, 1);
	//outTree->Branch("assocCluster", "std::vector<CorrectedCluster*>", &assocClusters, 64000, 1);
	//outTree->Branch("cutsWord", &cutsWord, "cutsWord[19]/O");

	outTree->Branch("rawBurst" ,"ROOTBurst", &rootBurst);
	outTree->Branch("rawEvent" ,"ROOTRawEvent", &rawEvent);
	outTree->Branch("corrEvent" ,"ROOTCorrectedEvent", &corrEvent);
	headerTree->Branch("geom" ,"NGeom", &rootGeom);
	headerTree->Branch("header" ,"ROOTFileHeader", &rootFileHeader);

	//Vertex
	new TH1D("vertexN", "vertexN", 10, 0, 10); 			///
	new TH1D("vertexZ", "vertexZ", 1300, -3000, 10000); ///
	new TH1D("vertexChi2", "vertexChi2", 1000, 0, 100); 	///
	new TH1D("vertexQ", "vertexQ", 10, -5, 5); 			///

	//Cluster
	new TH1D("clusterN", "clusterN", 10, 0, 10);					///
	new TH1D("clusterX", "clusterX", 600, -150, 150);			///
	new TH1D("clusterY", "clusterY", 600, -150, 150);			///
	new TH1D("clusterRadius", "clusterRadius", 600, 0, 150);		///
	new TH1D("clusterDPi", "clusterDPi", 500, 0, 250);			///
	new TH1D("clusterDep", "clusterDep", 500, 0, 250);			///
	new TH1D("clusterDem", "clusterDem", 500, 0, 250);			///
	new TH1D("clusterDUnep", "clusterDUnep", 500, 0, 250);		///
	new TH1D("clusterDUnem", "clusterDUnem", 500, 0, 250);		///
	new TH1D("clusterVertexTime", "clusterVertexTime", 200, 0, 200);///
	new TH1D("clusterGoodN", "clusterGoodN", 10, 0, 10);			///

	//Tracks
	new TH1D("trackN", "trackN", 10, 0, 10); 						///
	new TH1D("trackDCHTime", "trackDCHTime", 100, 0, 100);			///
	new TH1D("trackTrackTime", "trackTrackTime", 100, 0, 50); 		///
	new TH1D("trackVertexTime", "trackVertexTime", 100, 0, 100); 	///
	new TH1D("trackDCH1Radius", "trackDCH1Radius", 300, 0, 150); 	///
	new TH1D("trackDCH2Radius", "trackDCH2Radius", 300, 0, 150); 	///
	new TH1D("trackDCH4Radius", "trackDCH4Radius", 300, 0, 150); 	///
	new TH1D("trackLKrRadius", "trackLKrRadius", 400, 0, 200); 	///
	new TH1D("trackEOP", "trackEOP", 150, 0, 1.5);					///
	new TH1D("trackP", "trackP", 800, 0, 80); 						///
	new TH1D("track_ij_DCH1", "track_ij_DCH1", 400, 0, 200);		///
	new TH1D("track_ee_LKr", "track_ee_LKr", 500, 0, 250);		///
	new TH1D("track_epi_LKr", "track_epi_LKr", 500, 0, 250);		///

	//Gamma cluster
	new TH1D("gammaX", "gammaX", 600, -150, 150);				///
	new TH1D("gammaY", "gammaY", 600, -150, 150);				///
	new TH1D("gammaRadius", "gammaRadius", 300, 0, 150);			///
	new TH1D("gammaEnergy", "gammaEnergy", 800, 0, 80);				///
	new TH1D("gammaDeadCell", "gammaDeadCell", 60, 0, 30);			///
	new TH1D("gammaDCH1Radius", "gammaDCH1Radius", 300, 0, 150);	///

	//Kinematic before cuts
	new TH1D("bkinKP", "bkinKP", 100, 0, 80);			///
	new TH1D("bkinPiP", "bkinKPi", 100, 0, 80);			///
	new TH1D("bkinepP", "bkinKep", 100, 0, 80);			///
	new TH1D("bkinemP", "bkinKem", 100, 0, 80);			///
	new TH1D("bkinGammaP", "bkinKgamma", 100, 0, 80);	///
	new TH1D("bkinMee", "bkinMee", 160, 0, 0.16);			///
	new TH1D("bkinx", "bkinx", 1000, 0, 1);			///
	new TH1D("bkinMeeg", "bkinMeeg", 500, 0, 0.5);		///
	new TH1D("bkinMeegDiff", "bkinMeegDiff", 400, -0.2, 0.2);		///
	new TH1D("bkinMeegpi", "bkinMeegpi", 1000, 0, 0.6);	///
	new TH1D("bkinPTot", "bkinPTot", 300, 65, 80);		///
	new TH1D("bkinPt2", "bkinPt2", 1000, 0, 0.01);			///

	//Kinematic after cuts
	new TH1D("akinKP", "akinKP", 100, 0, 80);			///
	new TH1D("akinPiP", "akinKPi", 100, 0, 80);			///
	new TH1D("akinepP", "akinKep", 100, 0, 80);			///
	new TH1D("akinemP", "akinKem", 100, 0, 80);			///
	new TH1D("akinGammaP", "akinKgamma", 100, 0, 80);	///
	new TH1D("akinMee", "akinMee", 160, 0, 0.16);			///
	new TH1D("akinx", "akinx", 1000, 0, 1);			///
	new TH1D("akinMeeg", "akinMeeg", 500, 0, 0.5);		///
	new TH1D("akinMeegDiff", "akinMeegDiff", 400, -0.2, 0.2);		///
	new TH1D("akinMeegpi", "akinMeegpi", 1000, 0, 0.6);	///
	new TH1D("akinPTot", "akinPTot", 300, 65, 80);		///
	new TH1D("akinPt2", "akinPt2", 1000, 0, 0.01);			///



	//Selection histo
	new TH1D("Cuts", "Failed Cuts", 20, 0.5, 20.5);

	//Ana histo
	new TH1D("sel_NVertices", "Number of vertices", 10, 0, 10);
	new TH1D("sel_NTracks", "Number of tracks", 10, 0, 10);

	new TH1D("sel_DCH1Rad", "Extrapolated tracks radius on DCH1", 150, 0, 150);

	new TH1D("peegMass", "peegMass", 1000, 0.4, 0.6);
	new TH1D("eegMass", "eegMass", 1000, 0, 0.2);
	new TH1D("eeMass", "eeMass", 1000, 0, 0.2);

	new TH1D("xDistrib", "xDistrib", 1000, 0, 1);
	return 0;
}

int user_init() {
	selectOptions(gString);

	if(channel==KE2) nico_ke2Init();
	if(channel==PI0DALITZ) nico_pi0DalitzInit();
	/*----------- End of user C code -----------*/
	return 0;
}
