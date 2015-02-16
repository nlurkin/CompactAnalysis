#define OUTSIDECOMPACT

#include "exportClasses.h"
#include "mystructs.h"
#include "funLib.h"
#include <iostream>
using namespace std;

#include <TDirectory.h>
#include <TChain.h>
#include <cmath>
#include <TFile.h>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <boost/program_options.hpp>

///### Objects
ROOTRawEvent rawEvent;
ROOTCorrectedEvent corrEvent;
ROOTBurst rootBurst;
ROOTFileHeader rootFileHeader;
NGeom rootGeom;
ROOTMCEvent rootMC;
ROOTPhysicsEvent rootPhysics;
ROOTFileHeader outputFileHeader;

//### Ptr for TTree
namespace root_ptr{
	ROOTRawEvent *rawEvent_ptr = &rawEvent;
	ROOTCorrectedEvent *corrEvent_ptr = &corrEvent;
	ROOTBurst *rootBurst_ptr = &rootBurst;
	ROOTFileHeader *rootFileHeader_ptr = &rootFileHeader;
	NGeom *Geom = &rootGeom;
	ROOTMCEvent *rootMC_ptr = &rootMC;
}

NAbcog_params abcog_params;

//### Global Options
namespace options{
	bool optDebug;
	bool doOutput;
	int outputModulo;
	int periodKeep;
	bool exportAllEvents;
}

namespace parameters{
	cutsValues cutsDefinition;
	vector<eventID> badEventsList;
}

bool cutsWord[19];
bool mcBranched = false;

//### IO variables
namespace io_ptr{
	FILE* fprt, *fprt2;
	TTree *outTree, *outHeaderTree;
	TFile *outFile;

	TChain *inTree;
	TChain *headerTree;
}

ostream& operator<<(ostream &s, TVector3 v){
	s.precision(7);
	s << std::fixed;
	s << "( " << v.X() << " , " << v.Y() << " , " << v.Z() << " )";
	return s;
}

ostream& operator<<(ostream &s, TLorentzVector v){
	s.precision(7);
	s << std::fixed;
	s << "( " << v.X() << " , " << v.Y() << " , " << v.Z() << " , " << v.E() << " ) mass=" << v.M() << endl;
	return s;
}

bool readFile(TString fName, bool isList){

	if(isList){
		if(fName.Contains(".root")){
			cerr << "ROOT file cannot be a list" << endl;
			return false;
		}
		TString inputFileName;
		ifstream inputList(fName.Data());
		while(inputFileName.ReadLine(inputList)){
			if(inputFileName.Contains("/castor/") && !inputFileName.Contains("root://castorpublic.cern.ch//")){
								TString svcClass = getenv("STAGE_SVCCLASS");
								if(svcClass=="") svcClass="na62";
					inputFileName = "root://castorpublic.cern.ch//"+inputFileName+"?svcClass="+svcClass;
			}
			if(inputFileName.Contains("/eos/") && !inputFileName.Contains("root://eosna62.cern.ch//")){
				inputFileName = "root://eosna62.cern.ch//"+inputFileName;
			}
			io_ptr::inTree->AddFile(inputFileName);
			io_ptr::headerTree->AddFile(inputFileName);
		}
	}
	else{
		io_ptr::inTree->AddFile(fName);
		io_ptr::headerTree->AddFile(fName);
	}

	io_ptr::inTree->SetBranchAddress("rawBurst", &root_ptr::rootBurst_ptr);
	io_ptr::inTree->SetBranchAddress("rawEvent", &root_ptr::rawEvent_ptr);
	io_ptr::inTree->SetBranchAddress("corrEvent", &root_ptr::corrEvent_ptr);
	io_ptr::headerTree->SetBranchAddress("header", &root_ptr::rootFileHeader_ptr);
	io_ptr::inTree->SetBranchAddress("geom", &root_ptr::Geom);
	if(io_ptr::inTree->GetListOfBranches()->Contains("mc")){
		io_ptr::inTree->SetBranchAddress("mc", &root_ptr::rootMC_ptr);
		mcBranched = true;
	}

	return true;
}

bool openOutput(){
	io_ptr::outFile = gFile;
	io_ptr::outTree = new TTree("event", "Event");
	io_ptr::outHeaderTree = new TTree("header", "Header");

	io_ptr::outTree->Branch("rawBurst" ,"ROOTBurst", &rootBurst);
	io_ptr::outTree->Branch("rawEvent" ,"ROOTRawEvent", &rawEvent);
	io_ptr::outTree->Branch("corrEvent" ,"ROOTCorrectedEvent", &corrEvent);
	io_ptr::outTree->Branch("geom" ,"NGeom", &rootGeom);
	if(mcBranched) io_ptr::outTree->Branch("mc" ,"ROOTMCEvent", &rootMC);

	io_ptr::outTree->Branch("pi0dEvent" ,"ROOTPhysicsEvent", &rootPhysics);
	io_ptr::outHeaderTree->Branch("header" ,"ROOTFileHeader", &outputFileHeader);

	return true;
}

int pi0d_tracksAcceptance(){
	bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(root_ptr::Geom->Dch[0].PosChamber.x,root_ptr::Geom->Dch[0].PosChamber.y,root_ptr::Geom->Dch[0].PosChamber.z);
	TVector3 dch2(root_ptr::Geom->Dch[1].PosChamber.x,root_ptr::Geom->Dch[1].PosChamber.y,root_ptr::Geom->Dch[1].PosChamber.z);
	TVector3 dch4(root_ptr::Geom->Dch[3].PosChamber.x,root_ptr::Geom->Dch[3].PosChamber.y,root_ptr::Geom->Dch[3].PosChamber.z);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		propPos = propagateAfter(root_ptr::Geom->Lkr.z, t);
		lkrAcceptance = t.lkr_acc;
		if(options::optDebug) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(rootBurst.pbWall){
			if(options::optDebug) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}


		propPos = propagateCorrBefore(root_ptr::Geom->Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		if(options::optDebug) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(root_ptr::Geom->Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		if(options::optDebug) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(root_ptr::Geom->Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		if(options::optDebug) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;
	}

	return badTrack;
}

int pi0d_trackCombinationVeto(){
	int ntracks = 3;

	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			int trackID1 = corrEvent.goodTracks[i];
			int trackID2 = corrEvent.goodTracks[j];
			NPhysicsTrack t1 = corrEvent.pTrack[trackID1];
			NPhysicsTrack t2 = corrEvent.pTrack[trackID2];
			bad = false;
			if(options::optDebug) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(root_ptr::Geom->Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(root_ptr::Geom->Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(options::optDebug) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(rootBurst.isData){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				if(options::optDebug) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(root_ptr::Geom->Lkr.z, t1);
			propPos2 = propagateAfter(root_ptr::Geom->Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(trackID1==rootPhysics.pic.parentTrack || trackID2==rootPhysics.pic.parentTrack){
				if(options::optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				if(options::optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
				if(RLKr<=20) bad = true;
			}

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_identifyPi(int &piCandidate, bool &badElectron){
	double eop;

	int piCandidatesNb = 0;
	badElectron = false;

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		int goodTrackID = corrEvent.goodTracks[i];
		if(options::optDebug) cout << "\tTrying track :\t\t" << i <<  endl;

		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		if(options::optDebug) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
		if(eop<0.85){
			piCandidatesNb++;
			piCandidate = i;
		}
		else if(eop>1.15) badElectron = true;
	}

	return piCandidatesNb;
}

int pi0d_goodClusters(){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) conditions = 6;
	else conditions=5;

	if(options::optDebug) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options::optDebug) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options::optDebug) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(options::optDebug) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(options::optDebug) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack pi = corrEvent.pTrack[rootPhysics.pic.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[rootPhysics.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[rootPhysics.em.parentTrack];

		// separation from pi impact point >30cm
		propPos = propagateAfter(root_ptr::Geom->Lkr.z, pi);
		distance = distance2D(propPos, c.position);
		if(options::optDebug) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(root_ptr::Geom->Lkr.z, ep);
		distance = distance2D(propPos, c.position);
		if(options::optDebug) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(root_ptr::Geom->Lkr.z, em);
		distance = distance2D(propPos, c.position);
		if(options::optDebug) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		propPos = propagate(root_ptr::Geom->Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options::optDebug) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(root_ptr::Geom->Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options::optDebug) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(options::optDebug) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(options::optDebug) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			rootPhysics.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

int pi0d_failCut(int i){
	cutsWord[i] = false;
	if(options::optDebug) cout << "Event is not passing selection" << endl;
	if(options::doOutput) fprintf(io_ptr::fprt, "%i %i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, i);
	return 0;
}
void pi0d_passSelection(){
	if(options::optDebug) cout << "Event is passing selection" << endl;
	if(options::doOutput) fprintf(io_ptr::fprt2, "%i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp);
}

int nico_pi0DalitzSelect(){
	bool badAcceptance;

	int badCombis=0;
	int piTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int piCandNb;
	bool badElectron;

	int goodClusters;
	//int goodClusterID=-1;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	double pt;

	vector<double> vMass;
	vector<TVector3> vP;

	if(corrEvent.failedCond>=0) {
		outputFileHeader.NFailedEvents--;
		pi0d_failCut(corrEvent.failedCond);
		return -1;
	}
	if(options::optDebug) cout << endl;

	// 6) Track DCH time
	if(options::optDebug) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(rootBurst.isData){
		if(options::optDebug) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=parameters::cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
		if(options::optDebug) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=parameters::cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
		if(options::optDebug) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=parameters::cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
	}
	else{
		if(options::optDebug) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(options::optDebug) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(options::optDebug) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==parameters::cutsDefinition.boolBadTrack) {pi0d_failCut(7); return -1;}

	// 9) Identify candidates
	if(options::optDebug) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pi0d_identifyPi(piTrack, badElectron);
	if(options::optDebug) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=parameters::cutsDefinition.numPiCandidates) {pi0d_failCut(9); return -1;}

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if((int)i==piTrack) continue;

		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	//Create physics event from tracks
	rootPhysics.em.parentTrack = corrEvent.goodTracks[emTrack];
	rootPhysics.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.em.parentVertex = corrEvent.goodVertexID;
	rootPhysics.ep.parentTrack = corrEvent.goodTracks[epTrack];
	rootPhysics.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.ep.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pic.parentTrack = corrEvent.goodTracks[piTrack];
	rootPhysics.pic.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pic.parentVertex = corrEvent.goodVertexID;

	rootPhysics.em.P.SetVectM(corrEvent.pTrack[rootPhysics.em.parentTrack].momentum*corrEvent.pTrack[rootPhysics.em.parentTrack].p, Me);
	rootPhysics.ep.P.SetVectM(corrEvent.pTrack[rootPhysics.ep.parentTrack].momentum*corrEvent.pTrack[rootPhysics.ep.parentTrack].p, Me);
	rootPhysics.pic.P.SetVectM(corrEvent.pTrack[rootPhysics.pic.parentTrack].momentum*corrEvent.pTrack[rootPhysics.pic.parentTrack].p, Mpic);
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		rootPhysics.pic.pdgID = 211;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		rootPhysics.pic.pdgID = -211;
	}

	// 8) Track combination veto
	if(options::optDebug) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto();
	if(options::optDebug) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=parameters::cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return -1;}

	// 10) Bad electron cluster
	if(options::optDebug) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options::optDebug) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==parameters::cutsDefinition.boolBadECandidates) {pi0d_failCut(10); return -1;}

	// 12) Exactly 1 good LKr cluster
	if(options::optDebug) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters();
	if(options::optDebug) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=parameters::cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return -1;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event
	rootPhysics.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.gamma.parentVertex = corrEvent.goodVertexID;
	rootPhysics.gamma.P.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rootPhysics.gamma.vertex).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);
	rootPhysics.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pi0.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pi0.P = rootPhysics.em.P + rootPhysics.ep.P + rootPhysics.gamma.P;
	rootPhysics.kaon.P = rootPhysics.pic.P + rootPhysics.pi0.P;

	// 13) Photon candidate in LKr acceptance
	if(options::optDebug) cout << "~~~~ Cut 13 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options::optDebug) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=parameters::cutsDefinition.lkrAcceptance) {pi0d_failCut(13); return -1;}

	// 14) E_gamma>3GeV
	if(options::optDebug) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(options::optDebug) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << rootPhysics.gamma.P.E() << "\t <= 3 : rejected" << endl;
	if(rootPhysics.gamma.P.E()<=parameters::cutsDefinition.minGammaEnergy) {pi0d_failCut(14); return -1;}

	// 15) D_deadcell>2cm
	if(options::optDebug) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(options::optDebug) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=parameters::cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return -1;}

	if(options::optDebug) cout << "d_deadcell(pi) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.pic.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.pic.parentTrack].trackID].dDeadCell<=parameters::cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return -1;}

	if(options::optDebug) cout << "d_deadcell(ep) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.ep.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.ep.parentTrack].trackID].dDeadCell<=parameters::cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return -1;}

	if(options::optDebug) cout << "d_deadcell(em) :\t\t\t" << rawEvent.track[corrEvent.pTrack[rootPhysics.em.parentTrack].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[rootPhysics.em.parentTrack].trackID].dDeadCell<=parameters::cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return -1;}

	// 16) Photon DCH1 intercept >13cm
	if(options::optDebug) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(root_ptr::Geom->Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, rootPhysics.gamma.P.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options::optDebug) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=parameters::cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16); return -1;}

	//Start Kinematic cuts
	//Pt
	pt = rootPhysics.kaon.P.Perp2(corrEvent.kaonMomentum);

	//Mee
	rootPhysics.mee = (rootPhysics.em.P + rootPhysics.ep.P).M();
	rootPhysics.x = pow(rootPhysics.mee/Mpi0, 2.);
	rootPhysics.y = 2.*(rootPhysics.em.P.E() - rootPhysics.ep.P.E())/(Mpi0*(1-rootPhysics.x));

	// 11) Tracks momenta
	if(options::optDebug) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(options::optDebug) cout << "p_pi :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(options::optDebug) cout << "p_e+ :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(options::optDebug) cout << "p_e- :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p<=parameters::cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p>=parameters::cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}
	if(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p<=parameters::cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p>=parameters::cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}
	if(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p<=parameters::cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p>=parameters::cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}

	// 17) Total momentum 70<p<78
	if(options::optDebug) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options::optDebug) cout << "p_pieeg :\t\t\t" << rootPhysics.kaon.P.Vect().Mag() << "\t <70 || >78 : rejected" << endl;
	if(rootPhysics.kaon.P.Vect().Mag()<parameters::cutsDefinition.minTotalMomentum || rootPhysics.kaon.P.Vect().Mag()>parameters::cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return -1;}

	// 18) Transverse momentum^2 < 5E-4
	if(options::optDebug) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(options::optDebug) cout << "P_t^2 :\t\t" << pt << "\t >= " << parameters::cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=parameters::cutsDefinition.maxPt) {pi0d_failCut(18); return -1;}

	// 19) |M_eeg - M_pi0|<8 MeV
	if(options::optDebug) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options::optDebug) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options::optDebug) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=parameters::cutsDefinition.maxPi0MassDiff) {pi0d_failCut(19); return -1;}

	// 20) 0.475 < M_pieeg < 0.510
	if(options::optDebug) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(options::optDebug) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
	if(rootPhysics.kaon.P.M()<parameters::cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>parameters::cutsDefinition.maxKaonMassDiff) {pi0d_failCut(20); return -1;}

	pi0d_passSelection();
	return 0;
}

bool newEvent(int i, int &nevt){
	rootPhysics.clear();

	if(options::periodKeep!= 0 && rootBurst.period!=options::periodKeep) return false;
	if(i==0) cout << "First event: ";
	if(i % options::outputModulo == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r";
	if(parameters::badEventsList.size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, parameters::badEventsList)) return false;
	}
	if(i==0) cout << endl;

	++nevt;
	abcog_params = rootBurst.abcog_params;
	return nico_pi0DalitzSelect()==0;
}

int main(int argc, char **argv){
	namespace po = boost::program_options;

	//Options variables
	int nevt = -1;
	bool fileList = false;
	string fileName;
	string optString;
	string prefix;
	string cutsFile;
	string filterFile;

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("nevt,n", po::value<int>(), "max number of events")
		("file,i", po::value<string>(), "input file name")
		("list,l", po::value<string>(), "list of input files")
		("prefix,p", po::value<string>(), "prefix for output files")
		("debug,d", po::value<bool>(), "Activate verbose debugging")
		("filter,f", po::value<string>(), "Filter file")
		("dooutput", po::value<bool>(), "Activate output text files")
		("period", po::value<int>(), "Keep only events from specified period")
		("mod,m", po::value<int>(), "Event number printing modulo")
		("cuts,c", po::value<string>(), "Cuts file")
		("eall,e", po::value<bool>(), "Export all events, even failed")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	std::cout << ">>>> Received parameters:" << std::endl;
	for(const auto &it : vm){
		cout << "\t" << it.first << " = ";
		auto& value = it.second.value();
		if (auto v = boost::any_cast<int>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<bool>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<std::string>(&value))
			std::cout << *v;
		else if (auto v = boost::any_cast<TString>(&value))
			std::cout << *v;
		else
			std::cout << "**cast error**";
		cout << endl;
	}

	if (vm.count("nevt")) nevt = vm["nevt"].as<int>();

	/// String options
	if (vm.count("prefix")) prefix = vm["prefix"].as<string>();
	if (vm.count("debug")) options::optDebug = vm["debug"].as<bool>();
	if (vm.count("period")) options::periodKeep = vm["period"].as<int>();
	if (vm.count("dooutput")) options::doOutput = vm["dooutput"].as<bool>();
	if (vm.count("mod")) options::outputModulo = vm["mod"].as<int>();
	if (vm.count("cuts")) cutsFile = vm["cuts"].as<string>();
	if (vm.count("eall")) options::exportAllEvents = vm["eall"].as<bool>();

	/// Input (one of them mandatory)
	if (vm.count("list") and vm.count("file")){
		cerr << "Cannot specify both an input file and an input list at the same time" << endl;
		return -1;
	}
	else if(!vm.count("list") and !vm.count("file")){
		cerr << "Must specify either an input file or an input list" << endl;
		return -1;
	}
	else{
		if (vm.count("list")){
			fileList = true;
			fileName = vm["list"].as<string>();
		}
		if (vm.count("file")) fileName = vm["file"].as<string>();
	}

	//selectOptions(optString);
	common_init(prefix, filterFile, parameters::badEventsList, options::doOutput, io_ptr::fprt, io_ptr::fprt2);

	parseCutsValues("", parameters::cutsDefinition);
	printCuts(parameters::cutsDefinition);


	std::cout << "Output events every " << options::outputModulo << " events" << std::endl;
	if(!cutsFile.empty()) std::cout << "Using cuts at: " << cutsFile << std::endl;
	std::cout << "Debugging activated: " << (options::optDebug==true ? "Yes" : "No") << std::endl;
	if(options::periodKeep==0) std::cout << "Keeping period: All" << std::endl;
	else std::cout << "Keeping period: " << options::periodKeep << std::endl;
	if(options::doOutput) std::cout << "Text output files requested" << std::endl;
	if(options::exportAllEvents) std::cout << "Export all events requested" << std::endl;
	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cout << std::endl << std::endl;

	io_ptr::inTree = new TChain("event");
	io_ptr::headerTree = new TChain("header");

	if(!readFile(fileName, fileList)) return -1;
	openOutput();
	TString oldFile = "";
	int currFile = -1;
	int nevent = 0;


	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io_ptr::inTree->GetEntries() << endl;
	for(int i=0; i<io_ptr::inTree->GetEntries() && (nevt<0 || nevent<nevt ); ++i){
		io_ptr::inTree->GetEntry(i);
		if(oldFile!=io_ptr::inTree->GetCurrentFile()->GetName()){
			oldFile = io_ptr::inTree->GetCurrentFile()->GetName();
			++currFile;
			io_ptr::headerTree->GetEntry(currFile);
			outputFileHeader.NProcessedEvents += rootFileHeader.NProcessedEvents;
			outputFileHeader.NFailedEvents += rootFileHeader.NFailedEvents;
		}
		if(newEvent(i, nevent)){
			io_ptr::outTree->Fill();
			outputFileHeader.NPassedEvents++;
		}
		else{
			outputFileHeader.NFailedEvents++;
			if(options::exportAllEvents) io_ptr::outTree->Fill();
		}
	}
	cout << endl;

	cout << outputFileHeader.NProcessedEvents << endl;
	cout << outputFileHeader.NFailedEvents << endl;
	cout << outputFileHeader.NPassedEvents << endl;

	io_ptr::outTree->Write();
	io_ptr::outHeaderTree->Fill();
	io_ptr::outHeaderTree->Write();

	io_ptr::outFile->Close();
	return 0;
}
