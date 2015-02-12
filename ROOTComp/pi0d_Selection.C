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

#define FULLRUN false
#if FULLRUN
#define FAILCUT(i) pi0d_failCut(i);
#else
#define FAILCUT(i) {pi0d_failCut(i); return -1;}
#endif

#define DEBUG_ALL false

#define DEBUG_1 DEBUG_ALL || false
#define DEBUG_2 DEBUG_ALL || false
#define DEBUG_3 DEBUG_ALL || false
#define DEBUG_4 DEBUG_ALL || false
#define DEBUG_5 DEBUG_ALL || false
#define DEBUG_6 DEBUG_ALL || false
#define DEBUG_7 DEBUG_ALL || false
#define DEBUG_8 DEBUG_ALL || false
#define DEBUG_9 DEBUG_ALL || false
#define DEBUG_10 DEBUG_ALL || false
#define DEBUG_11 DEBUG_ALL || false
#define DEBUG_12 DEBUG_ALL || false
#define DEBUG_13 DEBUG_ALL || false
#define DEBUG_14 DEBUG_ALL || false
#define DEBUG_15 DEBUG_ALL || false
#define DEBUG_16 DEBUG_ALL || false
#define DEBUG_17 DEBUG_ALL || false
#define DEBUG_18 DEBUG_ALL || false
#define DEBUG_19 DEBUG_ALL || false
#define DEBUG_20 DEBUG_ALL || false

#define DEBUG_ DEBUG_1 || DEBUG_2 || DEBUG_3 || DEBUG_4 || DEBUG_5 || DEBUG_6 || DEBUG_7 || DEBUG_8 || DEBUG_9 || DEBUG_10 ||  DEBUG_11 || DEBUG_12 || DEBUG_13 || DEBUG_15 || DEBUG_16

///### Objects
ROOTRawEvent rawEvent;
ROOTCorrectedEvent corrEvent;
ROOTBurst rootBurst;
ROOTFileHeader rootFileHeader;
NGeom rootGeom;

//### Ptr for TTree
ROOTRawEvent *rawEvent_ptr = &rawEvent;
ROOTCorrectedEvent *corrEvent_ptr = &corrEvent;
ROOTBurst *rootBurst_ptr = &rootBurst;
ROOTFileHeader *rootFileHeader_ptr = &rootFileHeader;
NGeom *Geom = &rootGeom;

TChain *inTree;
TChain *headerTree;

NAbcog_params abcog_params;

double Mpi0 = 0.1349766;
double Mpic = 0.139570;
double Me = 0.00051099891;

//TODO TOFILL
bool optDebug;
bool pbWall;
bool dataOnly;
bool mcOnly;
bool noOutput;
FILE* fprt, *fprt2;
cutsValues cutsDefinition;
TTree *outTree;
bool cutsWord[19];
map<string,string> opts;
vector<eventID> badEventsList;
int channel;
int outputMod;
int periodKeep;
int iEvent;

/*
 * Channels defines
 */
#define KE2 1
#define PI0DALITZ 2
#define NONE 3



bool readFile(TString fName){
	inTree->AddFile(fName);
	headerTree->AddFile(fName);

	inTree->SetBranchAddress("pi0dBurst", &rootBurst_ptr);
	inTree->SetBranchAddress("rawEvent", &rawEvent_ptr);
	inTree->SetBranchAddress("corrEvent", &corrEvent_ptr);
	headerTree->SetBranchAddress("header", &rootFileHeader_ptr);
	inTree->SetBranchAddress("geom", &Geom);

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

	for(int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		propPos = propagateAfter(Geom->Lkr.z, t);
		lkrAcceptance = t.lkr_acc;
		//gH("trackLKrRadius")->Fill(distance2D(TVector3(0,0,0), propPos), corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(pbWall){
			if(DEBUG_8 || optDebug) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}


		propPos = propagateCorrBefore(Geom->Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		//gH("trackDCH1Radius")->Fill(radius, corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(Geom->Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		//gH("trackDCH2Radius")->Fill(radius, corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(Geom->Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		//gH("trackDCH4Radius")->Fill(radius, corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;
	}

	return badTrack;
}

int pi0d_trackCombinationVeto(int piTrack){
	int ntracks = 3;

	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;
	int totalCharge = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			NPhysicsTrack t1 = corrEvent.pTrack[corrEvent.goodTracks[i]];
			NPhysicsTrack t2 = corrEvent.pTrack[corrEvent.goodTracks[j]];
			bad = false;
			if(DEBUG_8 || optDebug) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(Geom->Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(Geom->Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			//gH("track_ij_DCH1")->Fill(RDCH1, corrEvent.weight);
			if(DEBUG_8 || optDebug) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(dataOnly){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				//gH("trackTrackTime")->Fill(tDiff, corrEvent.weight);
				if(DEBUG_8 || optDebug) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(Geom->Lkr.z, t1);
			propPos2 = propagateAfter(Geom->Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(i==piTrack || j==piTrack){
				//gH("track_epi_LKr")->Fill(RLKr, corrEvent.weight);
				if(DEBUG_8 || optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				//gH("track_ee_LKr")->Fill(RLKr, corrEvent.weight);
				if(DEBUG_8 || optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
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

	for(int i=0; i<corrEvent.goodTracks.size(); i++){
		int goodTrackID = corrEvent.goodTracks[i];
		if(DEBUG_9 || optDebug) cout << "\tTrying track :\t\t" << i <<  endl;

		//cout << corrEvent.pTrack[goodTrackID].E <<" " <<corrEvent.pTrack[goodTrackID].p << endl;
		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		//gH("trackEOP")->Fill(eop, corrEvent.weight);
		if(DEBUG_9 || optDebug) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
		if(eop<0.85){
			piCandidatesNb++;
			piCandidate = i;
		}
		else if(eop>1.15) badElectron = true;
	}

	return piCandidatesNb;
}

int pi0d_goodClusters(int piTrack, int epTrack, int emTrack, int ivtx, int &goodClusterID, double vertexTime){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(dataOnly) conditions = 6;
	else conditions=5;
	//if(mcOnly) conditions++;

	if(DEBUG_12 || optDebug) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(DEBUG_12 || optDebug) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	//gH("clusterN")->Fill(corrEvent.pCluster.size(), corrEvent.weight);
	for(int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		//gH("clusterX")->Fill(c.position.X(), corrEvent.weight);
		//gH("clusterY")->Fill(c.position.Y(), corrEvent.weight);
		//gH("clusterRadius")->Fill(distance2D(TVector3(0,0,0), c.position), corrEvent.weight);

		if(DEBUG_12 || optDebug) cout << "\tTrying cluster :\t" << i << endl;

		if(pbWall){
			if(DEBUG_12 || optDebug) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(DEBUG_12 || optDebug) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack pi = corrEvent.pTrack[corrEvent.goodTracks[piTrack]];
		NPhysicsTrack ep = corrEvent.pTrack[corrEvent.goodTracks[epTrack]];
		NPhysicsTrack em = corrEvent.pTrack[corrEvent.goodTracks[emTrack]];

		// separation from pi impact point >30cm
		propPos = propagateAfter(Geom->Lkr.z, pi);
		distance = distance2D(propPos, c.position);
		//gH("clusterDPi")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(Geom->Lkr.z, ep);
		//cout << printVector3(propPos) << endl;
		distance = distance2D(propPos, c.position);
		//gH("clusterDep")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(Geom->Lkr.z, em);
		distance = distance2D(propPos, c.position);
		//gH("clusterDem")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[epTrack]);
		propPos = propagate(Geom->Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		//gH("clusterDUnep")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[emTrack]);
		propPos = propagate(Geom->Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		//gH("clusterDUnem")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// |t_g - t_vtx|<10ns
		if(dataOnly){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - vertexTime);
			//gH("clusterVertexTime")->Fill(tDiff, corrEvent.weight);
			if(DEBUG_12 || optDebug) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(DEBUG_12 || optDebug) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			goodClusterID = i;
		}
	}
	return goodClusters;
}

int pi0d_failCut(int i){
	//TH1D *fails = (TH1D*)gDirectory->Get("Cuts");

	//fails->Fill(i, 1);
	cutsWord[i] = false;
	if(DEBUG_ || optDebug) cout << "Event is not passing selection" << endl;
	if(!noOutput) fprintf(fprt, "%i %i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, i);
	return 0;
}
void pi0d_passSelection(){
	if(DEBUG_ || optDebug) cout << "Event is passing selection" << endl;
	if(!noOutput) fprintf(fprt2, "%i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp);
}

int nico_pi0DalitzSelect(){
	int vtxNb = 0;
	int ivtx = -1;
	int extraTracks = 0;

	double vertexTime;
	bool badAcceptance;

	int badCombis=0;
	int piTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int piCandNb;
	bool badElectron;

	int goodClusters;
	int goodClusterID=-1;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	double pt;

	vector<double> vMass;
	vector<TVector3> vP;

	//To fill
	TVector3 kaonMomentum;
	double kaonP;
	double alpha;

	if(corrEvent.failedCond>=0) FAILCUT(corrEvent.failedCond);
	if(DEBUG_ || optDebug) cout << endl;

	ivtx = corrEvent.goodVertexID;
	vertexTime = rawEvent.vtx[ivtx].time;
	if(rawEvent.vtx[ivtx].charge==1){
		kaonMomentum = TVector3(abcog_params.pkdxdzp, abcog_params.pkdydzp, 1.).Unit();
		kaonP = abcog_params.pkp*(1+abcog_params.beta);
	}
	else{
		kaonMomentum = TVector3(abcog_params.pkdxdzm, abcog_params.pkdydzm, 1.).Unit();
		kaonP = abcog_params.pkm*(1+abcog_params.beta);
	}

	double weight = 1.;
	if(mcOnly && false){
		weight = 1 + alpha*pow(kaonP-74.,2);
	}
	corrEvent.weight = weight;

	// 6) Track DCH time
	if(DEBUG_6 || optDebug) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(dataOnly){
		if(DEBUG_6 || optDebug) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) FAILCUT(6)
		if(DEBUG_6 || optDebug) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) FAILCUT(6)
		if(DEBUG_6 || optDebug) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=cutsDefinition.maxTrackTime) FAILCUT(6)
	}
	else{
		if(DEBUG_6 || optDebug) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(DEBUG_7 || optDebug) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(DEBUG_7 || optDebug) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==cutsDefinition.boolBadTrack) FAILCUT(7)

	// 9) Identify candidates
	if(DEBUG_9 || optDebug) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pi0d_identifyPi(piTrack, badElectron);
	if(DEBUG_9 || optDebug) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=cutsDefinition.numPiCandidates) FAILCUT(9)

	for(int i=0; i<corrEvent.goodTracks.size(); i++){
		if(i==piTrack) continue;

		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	// 8) Track combination veto
	if(DEBUG_8 || optDebug) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto(piTrack);
	if(DEBUG_8 || optDebug) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=cutsDefinition.numBadTrackCombi) FAILCUT(8)


					//if(piTrack==-1 || epTrack==-1 || emTrack==-1)
		//FILCUT(9)


	//fullEvent.em = goodTracks[emTrack];
	//fullEvent.ep = goodTracks[epTrack];
	//fullEvent.pip= goodTracks[piTrack];

	// 10) Bad electron cluster
	if(DEBUG_10 || optDebug) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(DEBUG_10 || optDebug) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==cutsDefinition.boolBadECandidates) FAILCUT(10)

	// 12) Exactly 1 good LKr cluster
	if(DEBUG_12 || optDebug) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters(piTrack, epTrack, emTrack, ivtx, goodClusterID, vertexTime);
	//gH("clusterGoodN")->Fill(goodClusters, corrEvent.weight);
	if(DEBUG_12 || optDebug) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=cutsDefinition.numAddGoodCluster) FAILCUT(12)

	//fullEvent.gamma = vCluster[goodClusterID];

	if(goodClusterID==-1){
		return 0;
	}

	TVector3 pGamma = (corrEvent.pCluster[goodClusterID].position - rawEvent.vtx[rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[piTrack]].trackID].vtxID].position).Unit();
	//fullEvent.pGamma = (vCluster[goodClusterID]->position - goodTracks[piTrack]->vertex).Unit();
	//gH("gammaEnergy")->Fill(corrEvent.pCluster[goodClusterID].E, corrEvent.weight);
	//gH("gammaX")->Fill(corrEvent.pCluster[goodClusterID].position.X(), corrEvent.weight);
	//gH("gammaY")->Fill(corrEvent.pCluster[goodClusterID].position.Y(), corrEvent.weight);
	//gH("gammaRadius")->Fill(distance2D(TVector3(0,0,0), corrEvent.pCluster[goodClusterID].position), corrEvent.weight);
	//gH("gammaDeadCell")->Fill(rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].dDeadCell, corrEvent.weight);


	// 13) Photon candidate in LKr acceptance
	//TODO check it's the good way
	if(DEBUG_13 || optDebug) cout << "~~~~ Cut 13 ~~~~" << endl;
	//propPos = propagate(Geom->Lkr.z, vCluster[goodClusterID]->position, photonMomentum);
	//lkrAcceptance = LKr_acc(sbur->nrun, vCluster[goodClusterID]->position.X(), vCluster[goodClusterID]->position.Y(), 8);
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].lkr_acc;
	//lkrAcceptance = LKr_acc(sbur->nrun, propPos.X(), propPos.Y(), 8);
	if(DEBUG_13 || optDebug) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=cutsDefinition.lkrAcceptance) FAILCUT(13)

	// 14) E_gamma>3GeV
	if(DEBUG_14 || optDebug) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(DEBUG_14 || optDebug) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << corrEvent.pCluster[goodClusterID].E << "\t <= 3 : rejected" << endl;
	if(corrEvent.pCluster[goodClusterID].E<=cutsDefinition.minGammaEnergy) FAILCUT(14)

	// 15) D_deadcell>2cm
	if(DEBUG_15 || optDebug) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(DEBUG_15 || optDebug) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].dDeadCell<=cutsDefinition.minDeadCellDist) FAILCUT(15)

	//cout << printVector3(goodTracks[piTrack].vertex) << endl;
	//cout << printVector3(vCluster[goodClusterID]->position) << endl;
	//cout << printVector3(vCluster[goodClusterID]->energy*photonMomentum) << endl;

	// 16) Photon DCH1 intercept >13cm
	if(DEBUG_16 || optDebug) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(Geom->Dch[0].PosChamber.z, corrEvent.pCluster[goodClusterID].position, pGamma);
	//radius = sqrt(pow(propPos.X()-Geom->Dch[0].PosChamber.x,2) + pow(propPos.Y()-Geom->Dch[0].PosChamber.y,2));
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	//gH("gammaDCH1Radius")->Fill(radius, corrEvent.weight);
	if(DEBUG_16 || optDebug) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=cutsDefinition.minGammaDCHRadius) FAILCUT(16)

	//Start Kinematic cuts
	//Total P
	TVector3 pTotal = corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[piTrack]].momentum +
	corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[emTrack]].momentum +
	corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[epTrack]].momentum +
	corrEvent.pCluster[goodClusterID].E*pGamma;
	/*fullEvent.pTotal = goodTracks[piTrack]->pMag*goodTracks[piTrack]->momentum +
			goodTracks[emTrack]->pMag*goodTracks[emTrack]->momentum +
			goodTracks[epTrack]->pMag*goodTracks[epTrack]->momentum +
			vCluster[goodClusterID]->energy*fullEvent.pGamma;*/
	//Pt
	//pt = fullEvent.pTotal.Perp2(kaonMomentum);
	pt = pTotal.Perp2(kaonMomentum);
	//fullEvent.pt2 = pt;

	//Mee
	vMass.push_back(Me);
	vP.push_back(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[epTrack]].momentum);
	vMass.push_back(Me);
	vP.push_back(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[emTrack]].momentum);
	//fullEvent.mee = sqrt(invMass2(vMass, vP));
	//fullEvent.x = pow(fullEvent.mee/Mpi0, 2.);
	double mee = sqrt(invMass2(vMass, vP));
	double x = pow(mee/Mpi0, 2.);

	//Meeg
	vMass.push_back(0.0);
	vP.push_back(corrEvent.pCluster[goodClusterID].E*pGamma);
	//fullEvent.mPi0 = sqrt(invMass2(vMass, vP));
	double mPi0 = sqrt(invMass2(vMass, vP));

	//Mpieeg
	vMass.push_back(Mpic);
	vP.push_back(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p*corrEvent.pTrack[corrEvent.goodTracks[piTrack]].momentum);
	//fullEvent.mK = sqrt(invMass2(vMass, vP));
	double mK = sqrt(invMass2(vMass, vP));

	//gH("bkinKP")->Fill(kaonMomentum.Mag(), corrEvent.weight);
	//gH("bkinPiP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p, corrEvent.weight);
	//gH("bkinemP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p, corrEvent.weight);
	//gH("bkinepP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p, corrEvent.weight);
	//gH("bkinGammaP")->Fill((corrEvent.pCluster[goodClusterID].E*pGamma).Mag(), corrEvent.weight);

	//gH("bkinMee")->Fill(mee, corrEvent.weight);
	//gH("bkinx")->Fill(mee/Mpi0, corrEvent.weight);
	//gH("bkinMeeg")->Fill(mPi0, corrEvent.weight);
	//gH("bkinMeegDiff")->Fill(mPi0-Mpi0, corrEvent.weight);
	//gH("bkinMeegpi")->Fill(mK, corrEvent.weight);

	//gH("bkinPTot")->Fill(pTotal.Mag(), corrEvent.weight);
	//gH("bkinPt2")->Fill(pt, corrEvent.weight);

	// 11) Tracks momenta
	if(DEBUG_11 || optDebug) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(DEBUG_11 || optDebug) cout << "p_pi :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(DEBUG_11 || optDebug) cout << "p_e+ :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(DEBUG_11 || optDebug) cout << "p_e- :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p>=cutsDefinition.maxTrackMomentum) FAILCUT(11)
	if(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p>=cutsDefinition.maxTrackMomentum) FAILCUT(11)
	if(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p<=cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p>=cutsDefinition.maxTrackMomentum) FAILCUT(11)

	// 17) Total momentum 70<p<78
	if(DEBUG_17 || optDebug) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(DEBUG_17 || optDebug) cout << "p_pieeg :\t\t\t" << pTotal.Mag() << "\t <70 || >78 : rejected" << endl;
	if(pTotal.Mag()<cutsDefinition.minTotalMomentum || pTotal.Mag()>cutsDefinition.maxTotalMomentum) FAILCUT(17)

	// 18) Transverse momentum^2 < 5E-4
	if(DEBUG_18 || optDebug) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(DEBUG_18 || optDebug) cout << "P_t^2 :\t\t" << pt << "\t >= " << cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=cutsDefinition.maxPt) FAILCUT(18)

	// 19) |M_eeg - M_pi0|<8 MeV
	if(DEBUG_19 || optDebug) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(DEBUG_19 || optDebug) cout << "M_ee :\t\t" << mee << endl;
	if(DEBUG_19 || optDebug) cout << "|M_eeg - M_pi0| :\t\t" << fabs(mPi0-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(mPi0-Mpi0)>=cutsDefinition.maxPi0MassDiff) FAILCUT(19)

	// 20) 0.475 < M_pieeg < 0.510
	if(DEBUG_20 || optDebug) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(DEBUG_20 || optDebug) cout << "M_pieeg :\t\t" << mK << "\t <0.475 || >0.510: rejected" << endl;
	if(mK<cutsDefinition.minKaonMassDiff || mK>cutsDefinition.maxKaonMassDiff) FAILCUT(20)
	//gH("akinKP")->Fill(kaonMomentum.Mag(), corrEvent.weight);
	//gH("akinPiP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p, corrEvent.weight);
	//gH("akinemP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p, corrEvent.weight);
	//gH("akinepP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p, corrEvent.weight);
	//gH("akinGammaP")->Fill((corrEvent.pCluster[goodClusterID].E*pGamma).Mag(), corrEvent.weight);

	//gH("akinMee")->Fill(mee, corrEvent.weight);
	//gH("akinx")->Fill(pow(mee/Mpi0,2), corrEvent.weight);
	//gH("akinMeeg")->Fill(mPi0, corrEvent.weight);
	//gH("akinMeegDiff")->Fill(mPi0-Mpi0, corrEvent.weight);
	//gH("akinMeegpi")->Fill(mK, corrEvent.weight);

	//gH("akinPTot")->Fill(pTotal.Mag(), corrEvent.weight);
	//gH("akinPt2")->Fill(pt, corrEvent.weight);
	//outTree->FlushBaskets();
	//outTree->Fill();

	pi0d_passSelection();
	return 0;
}

int common_init(string filePrefix){
	int i, j, k;
	int l, m, n;
	int cpd, cell;

	int runNum, burstNum, timestamp;
	vector<eventID>::iterator it;


	string outRoot = "outfile.root";
	string outFile = "compact.txt";
	string outPass = "compactpass.txt";
	if(filePrefix.find('~')!=string::npos) filePrefix=filePrefix.replace(filePrefix.find('~'), 1, string("/afs/cern.ch/user/n/nlurkin"));
	if(filePrefix.length()>0){
		outRoot = filePrefix + ".root";
		outFile = filePrefix + ".txt";
		outPass = filePrefix + "pass.txt";
	}

	//Load badEvents list for debugging
	if(opts.count("filter")>0){
		cout << ">>>> Filtering events from file " << opts["filter"] << endl;
		FILE *badEvents = fopen(opts["filter"].c_str(), "r");
		if(badEvents!=NULL){
			while(fscanf(badEvents, "%i %i %i", &runNum, &burstNum, &timestamp) != EOF){
				badEventsList.push_back(eventID(runNum, burstNum, timestamp));
			}
			fclose(badEvents);
		}
		else{
			cout << "Unable to open filter file" << endl;
		}
		cout << "\t" << badEventsList.size() << " events in filter list" << endl;
		for(it=badEventsList.begin(); it!=badEventsList.end();it++){
			cout << "\t\t" << (*it).rnum << " " << (*it).bnum << " " << (*it).timestamp << endl;
		}
	}

	if(opts.count("nooutput")==0){
		noOutput = false;
		fprt=fopen(outFile.c_str(),"w");
		fprt2=fopen(outPass.c_str(),"w");
	}
	else noOutput = true;

	gFile = TFile::Open(outRoot.c_str(), "RECREATE");

	return 0;
}

void selectOptions(string s){
	map<string,string>::iterator it;

	opts = parseOptions(s);

	cout << endl << ">>>>>>>>>>>>>>>>>>>>> Initialization" << endl;
	if(opts.count("h")!=0){
		cout << ">>>> Help " << endl;
		cout << ">>>> Syntax: param=value:param=value" << endl;
		cout << ">>>> List of parameters:" << endl;
		cout << ">>>> h: This help" << endl;
		cout << ">>>> prefix: Output file names to use (without extension)" << endl;
		cout << ">>>> can: ke2 | pi0d" << endl;
		cout << ">>>> debug: Activate debugging" << endl;
		cout << ">>>> ff: Type of form factor (0=1, 1=x, 2=x^2)" << endl;
		cout << ">>>> nooutput: Don't create output txt files" << endl;
		cout << ">>>> period: keep only events from this period" << endl;
		cout << ">>>> mod: print events index every mod events" << endl;
		cout << ">>>> cuts: specify cuts file" << endl;
		cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		exit(0);
	}
	cout << ">>>> Received parameters:" << endl;
	for(it=opts.begin(); it!=opts.end(); it++){
		cout << "\t" << it->first << " = " << it->second << endl;
	}

	common_init(opts["prefix"]);
	string chanName;

	if(strcmp(opts["can"].c_str(), "ke2")==0){
		channel=KE2;
		chanName = "KE2";
	}
	else if(strcmp(opts["can"].c_str(), "pi0d")==0){
		channel=PI0DALITZ;
		chanName = "PI0Dalitz";
	}
	else if(strcmp(opts["can"].c_str(), "none")==0){
		channel=NONE;
		chanName = "None";
	}
	else if(opts["can"].length()==0){
		channel=PI0DALITZ;
		chanName = "Pi0Dalitz";
	}

	if(opts["mod"].length()!=0) outputMod = atoi(opts["mod"].c_str());
	else outputMod = 1;

	if(opts.count("debug")!=0) optDebug = true;
	else optDebug = false;

	if(opts.count("period")!=0) periodKeep = atoi(opts["period"].c_str());
	else periodKeep = 0;

	string cutsFileName;
	if(opts.count("cuts")!=0) cutsFileName = opts["cuts"];
	else cutsFileName = "";
	parseCutsValues(cutsFileName);
	printCuts();

	cout << "Starting on channel: " << chanName << endl;
	cout << "Output events every " << outputMod << " events" << endl;
	if(cutsFileName.length()>0) cout << "Using cuts at: " << cutsFileName << endl;
	cout << "Debugging activated: " << (optDebug==true ? "Yes" : "No") << endl;
	if(periodKeep==0) cout << "Keeping period: All" << endl;
	else cout << "Keeping period: " << periodKeep << endl;
	if(noOutput) cout << "No file output requested" << endl;
	cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	cout << endl << endl;
}
int main(int argc, char **argv){
	optDebug = false;
	pbWall = false;
	dataOnly = true;
	mcOnly = false;
	noOutput = false;
	outputMod = 10000;

	int opt;
	string optString;
	int nevt = -1;
	while ((opt = getopt(argc, argv, "s:n:")) != -1) {
		switch (opt) {
		case 's':
			optString = optarg;
			break;
		case 'n':
			//Input file
			nevt = TString(optarg).Atoi();
			break;
		}
	}

	selectOptions(optString);
	outTree = outTree = new TTree("event", "Event");

	parseCutsValues("");
	printCuts();
	inTree = new TChain("event");
	headerTree = new TChain("header");

	readFile("../outfile.root");

	for(int i=0; i<inTree->GetEntries() && (nevt<0 || i<nevt ); ++i){
		inTree->GetEntry(i);

		//if(periodKeep!= 0 && period!=periodKeep) return 0;
		if(i==0) cout << "First event: ";
		if(i % outputMod == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r";
		int passEvent;
		if(opts.count("filter")>0){
			if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp)) continue;
		}
		if(i==0) cout << endl;

		iEvent++;

		abcog_params = rootBurst.abcog_params;
		pbWall = rootBurst.pbWall;
		dataOnly = rootBurst.isData;
		mcOnly = rootBurst.isMC;
		nico_pi0DalitzSelect();
	}
	cout << endl;
}
