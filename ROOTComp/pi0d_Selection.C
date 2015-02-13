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
ROOTRawEvent *rawEvent_ptr = &rawEvent;
ROOTCorrectedEvent *corrEvent_ptr = &corrEvent;
ROOTBurst *rootBurst_ptr = &rootBurst;
ROOTFileHeader *rootFileHeader_ptr = &rootFileHeader;
NGeom *Geom = &rootGeom;
ROOTMCEvent *rootMC_ptr = &rootMC;

NAbcog_params abcog_params;

//### Global Options
bool optDebug;
bool noOutput;
cutsValues cutsDefinition;
map<string,string> opts;
vector<eventID> badEventsList;
int channel;
int outputMod;
int periodKeep;
bool exportAllEvents;

bool cutsWord[19];
bool mcBranched = false;

//### IO variables
FILE* fprt, *fprt2;
TTree *outTree, *outHeaderTree;
TFile *outFile;

TChain *inTree;
TChain *headerTree;

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
			inTree->AddFile(inputFileName);
			headerTree->AddFile(inputFileName);
		}
	}
	else{
		inTree->AddFile(fName);
		headerTree->AddFile(fName);
	}

	inTree->SetBranchAddress("rawBurst", &rootBurst_ptr);
	inTree->SetBranchAddress("rawEvent", &rawEvent_ptr);
	inTree->SetBranchAddress("corrEvent", &corrEvent_ptr);
	headerTree->SetBranchAddress("header", &rootFileHeader_ptr);
	inTree->SetBranchAddress("geom", &Geom);
	if(inTree->GetListOfBranches()->Contains("mc")){
		inTree->SetBranchAddress("mc", &rootMC_ptr);
		mcBranched = true;
	}

	return true;
}

bool openOutput(){
	outFile = gFile;
	outTree = new TTree("event", "Event");
	outHeaderTree = new TTree("header", "Header");

	outTree->Branch("pi0dBurst" ,"ROOTBurst", &rootBurst);
	outTree->Branch("rawEvent" ,"ROOTRawEvent", &rawEvent);
	outTree->Branch("corrEvent" ,"ROOTCorrectedEvent", &corrEvent);
	outTree->Branch("geom" ,"NGeom", &rootGeom);
	if(mcBranched) outTree->Branch("mc" ,"ROOTMCEvent", &rootMC);

	outTree->Branch("pi0dEvent" ,"ROOTPhysicsEvent", &rootPhysics);
	outHeaderTree->Branch("header" ,"ROOTFileHeader", &outputFileHeader);

	return true;
}

int pi0d_tracksAcceptance(){
	bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(Geom->Dch[0].PosChamber.x,Geom->Dch[0].PosChamber.y,Geom->Dch[0].PosChamber.z);
	TVector3 dch2(Geom->Dch[1].PosChamber.x,Geom->Dch[1].PosChamber.y,Geom->Dch[1].PosChamber.z);
	TVector3 dch4(Geom->Dch[3].PosChamber.x,Geom->Dch[3].PosChamber.y,Geom->Dch[3].PosChamber.z);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		propPos = propagateAfter(Geom->Lkr.z, t);
		lkrAcceptance = t.lkr_acc;
		if(optDebug) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(rootBurst.pbWall){
			if(optDebug) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}


		propPos = propagateCorrBefore(Geom->Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		if(optDebug) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(Geom->Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		if(optDebug) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(Geom->Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		if(optDebug) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
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
			if(optDebug) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(Geom->Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(Geom->Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(optDebug) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(rootBurst.isData){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				if(optDebug) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(Geom->Lkr.z, t1);
			propPos2 = propagateAfter(Geom->Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(trackID1==rootPhysics.pic.parentTrack || trackID2==rootPhysics.pic.parentTrack){
				if(optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				if(optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
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
		if(optDebug) cout << "\tTrying track :\t\t" << i <<  endl;

		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		if(optDebug) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
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

	if(optDebug) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(optDebug) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(optDebug) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(optDebug) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(optDebug) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack pi = corrEvent.pTrack[rootPhysics.pic.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[rootPhysics.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[rootPhysics.em.parentTrack];

		// separation from pi impact point >30cm
		propPos = propagateAfter(Geom->Lkr.z, pi);
		distance = distance2D(propPos, c.position);
		if(optDebug) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(Geom->Lkr.z, ep);
		distance = distance2D(propPos, c.position);
		if(optDebug) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(Geom->Lkr.z, em);
		distance = distance2D(propPos, c.position);
		if(optDebug) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		propPos = propagate(Geom->Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(optDebug) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(Geom->Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(optDebug) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(optDebug) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(optDebug) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			rootPhysics.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

int pi0d_failCut(int i){
	cutsWord[i] = false;
	if(optDebug) cout << "Event is not passing selection" << endl;
	if(!noOutput) fprintf(fprt, "%i %i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, i);
	return 0;
}
void pi0d_passSelection(){
	if(optDebug) cout << "Event is passing selection" << endl;
	if(!noOutput) fprintf(fprt2, "%i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp);
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
	if(optDebug) cout << endl;

	// 6) Track DCH time
	if(optDebug) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(rootBurst.isData){
		if(optDebug) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
		if(optDebug) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
		if(optDebug) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) {pi0d_failCut(6); return -1;}
	}
	else{
		if(optDebug) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(optDebug) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(optDebug) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==cutsDefinition.boolBadTrack) {pi0d_failCut(7); return -1;}

	// 9) Identify candidates
	if(optDebug) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pi0d_identifyPi(piTrack, badElectron);
	if(optDebug) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=cutsDefinition.numPiCandidates) {pi0d_failCut(9); return -1;}

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
	rootPhysics.ep.parentTrack = corrEvent.goodTracks[epTrack];
	rootPhysics.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pic.parentTrack = corrEvent.goodTracks[piTrack];
	rootPhysics.pic.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;

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
	if(optDebug) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto();
	if(optDebug) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return -1;}

	// 10) Bad electron cluster
	if(optDebug) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(optDebug) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==cutsDefinition.boolBadECandidates) {pi0d_failCut(10); return -1;}

	// 12) Exactly 1 good LKr cluster
	if(optDebug) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters();
	if(optDebug) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return -1;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event
	rootPhysics.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.gamma.P.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rootPhysics.gamma.vertex).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);
	rootPhysics.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pi0.P = rootPhysics.em.P + rootPhysics.ep.P + rootPhysics.gamma.P;
	rootPhysics.kaon.P = rootPhysics.pic.P + rootPhysics.pi0.P;

	// 13) Photon candidate in LKr acceptance
	if(optDebug) cout << "~~~~ Cut 13 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(optDebug) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=cutsDefinition.lkrAcceptance) {pi0d_failCut(13); return -1;}

	// 14) E_gamma>3GeV
	if(optDebug) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(optDebug) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << rootPhysics.gamma.P.E() << "\t <= 3 : rejected" << endl;
	if(rootPhysics.gamma.P.E()<=cutsDefinition.minGammaEnergy) {pi0d_failCut(14); return -1;}

	// 15) D_deadcell>2cm
	if(optDebug) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(optDebug) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return -1;}

	// 16) Photon DCH1 intercept >13cm
	if(optDebug) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(Geom->Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, rootPhysics.gamma.P.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(optDebug) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16); return -1;}

	//Start Kinematic cuts
	//Pt
	pt = rootPhysics.kaon.P.Perp2(corrEvent.kaonMomentum);

	//Mee
	rootPhysics.mee = (rootPhysics.em.P + rootPhysics.ep.P).M();
	rootPhysics.x = pow(rootPhysics.mee/Mpi0, 2.);
	rootPhysics.y = 2.*(rootPhysics.em.P.E() - rootPhysics.ep.P.E())/(Mpi0*(1-rootPhysics.x));

	// 11) Tracks momenta
	if(optDebug) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(optDebug) cout << "p_pi :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(optDebug) cout << "p_e+ :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(optDebug) cout << "p_e- :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p>=cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}
	if(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p>=cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}
	if(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p>=cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return -1;}

	// 17) Total momentum 70<p<78
	if(optDebug) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(optDebug) cout << "p_pieeg :\t\t\t" << rootPhysics.kaon.P.Vect().Mag() << "\t <70 || >78 : rejected" << endl;
	if(rootPhysics.kaon.P.Vect().Mag()<cutsDefinition.minTotalMomentum || rootPhysics.kaon.P.Vect().Mag()>cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return -1;}

	// 18) Transverse momentum^2 < 5E-4
	if(optDebug) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(optDebug) cout << "P_t^2 :\t\t" << pt << "\t >= " << cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=cutsDefinition.maxPt) {pi0d_failCut(18); return -1;}

	// 19) |M_eeg - M_pi0|<8 MeV
	if(optDebug) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(optDebug) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(optDebug) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=cutsDefinition.maxPi0MassDiff) {pi0d_failCut(19); return -1;}

	// 20) 0.475 < M_pieeg < 0.510
	if(optDebug) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(optDebug) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
	if(rootPhysics.kaon.P.M()<cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>cutsDefinition.maxKaonMassDiff) {pi0d_failCut(20); return -1;}

	pi0d_passSelection();
	return 0;
}

bool newEvent(int i, int &nevt){
	rootPhysics.clear();

	if(periodKeep!= 0 && rootBurst.period!=periodKeep) return false;
	if(i==0) cout << "First event: ";
	if(i % outputMod == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r";
	if(opts.count("filter")>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp)) return false;
	}
	if(i==0) cout << endl;

	++nevt;
	abcog_params = rootBurst.abcog_params;
	return nico_pi0DalitzSelect()==0;
}

int main(int argc, char **argv){
	int opt;
	string optString;
	int nevt = -1;
	bool fileList = false;
	TString fileName;

	while ((opt = getopt(argc, argv, "s:n:l:i:")) != -1) {
		switch (opt) {
		case 's':
			optString = optarg;
			break;
		case 'n':
			nevt = TString(optarg).Atoi();
			break;
		case 'l':
			fileName = optarg;
			fileList = true;
			break;
		case 'i':
			fileName = optarg;
			break;
		}
	}

	selectOptions(optString);

	parseCutsValues("");
	printCuts();
	inTree = new TChain("event");
	headerTree = new TChain("header");

	if(!readFile(fileName, fileList)) return -1;
	openOutput();
	TString oldFile = "";
	int currFile = -1;
	int nevent = 0;


	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << inTree->GetEntries() << endl;
	for(int i=0; i<inTree->GetEntries() && (nevt<0 || nevent<nevt ); ++i){
		inTree->GetEntry(i);
		if(oldFile!=inTree->GetCurrentFile()->GetName()){
			oldFile = inTree->GetCurrentFile()->GetName();
			++currFile;
			headerTree->GetEntry(currFile);
			outputFileHeader.NProcessedEvents += rootFileHeader.NProcessedEvents;
			outputFileHeader.NFailedEvents += rootFileHeader.NFailedEvents;
		}
		if(newEvent(i, nevent)){
			outTree->Fill();
			outputFileHeader.NPassedEvents++;
		}
		else{
			outputFileHeader.NFailedEvents++;
			if(exportAllEvents) outTree->Fill();
		}
	}
	cout << endl;

	cout << outputFileHeader.NProcessedEvents << endl;
	cout << outputFileHeader.NFailedEvents << endl;
	cout << outputFileHeader.NPassedEvents << endl;

	outTree->Write();
	outHeaderTree->Fill();
	outHeaderTree->Write();

	outFile->Close();
	return 0;
}
