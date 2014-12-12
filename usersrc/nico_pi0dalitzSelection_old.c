/*
 * nico_pi0dalitzSelection.c
 *
 *  Created on: 18 Feb 2014
 *      Author: ncl
 */

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <math.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <iterator>
#include <bitset>
#include <iomanip>
using namespace std;

#define FULLRUN false
#if FULLRUN
#define FAILCUT(i) pi0d_failCut(sbur,sevt,i);
#else
#define FAILCUT(i) {pi0d_failCut(sbur,sevt,i); return -1;}
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

double pi0d_getVertexTime(superCmpEvent *sevt, int ivtx){
	double mean = 0;

	for(int i=0; i<sevt->vtx[ivtx].Nvtxtrack; i++){
		mean += sevt->track[vtrack[sevt->vtx[ivtx].vtxtrack[i].iTrack]->trackID].time;
	}

	return mean/(double)sevt->vtx[ivtx].Nvtxtrack;
}

int pi0d_countVtx3Tracks(superBurst *sbur, superCmpEvent *sevt, int &ivtx){
	int vtxNb = 0;
	//int totalCharge;

	TH1D *nTracks = (TH1D*)gFile->Get("NTracks");
	TH1D *nVertices = (TH1D*)gFile->Get("NVertices");
	nTracks->Reset("M");
	nVertices->Reset("M");

	nVertices->Fill(sevt->Nvtx, fullEvent.weight);

	if(DEBUG_2 || optDebug) cout << "Total Number of vertices :\t" << sevt->Nvtx << " ==1" << endl;
	if(sevt->Nvtx!=1) return sevt->Nvtx;
	for(int i=0; i<sevt->Nvtx; i++){
		if(DEBUG_2 || optDebug) cout << "\tTrying vertex :\t\t" << i << endl;
		if(DEBUG_2 || optDebug) cout << "\tNumber of tracks:\t" << sevt->vtx[i].Nvtxtrack << endl;
		nTracks->Fill(sevt->vtx[i].Nvtxtrack, fullEvent.weight);
		if(sevt->vtx[i].Nvtxtrack==3){
			if(DEBUG_2 || optDebug) cout << "\tTotal charge :\t" << sevt->vtx[i].charge << endl;
			if(sevt->vtx[i].charge==beamCharge || noTestCharge){
				vtxNb++;
				ivtx = i;
			}
		}
	}

	return vtxNb;
}

int pi0d_extraTrackVeto(superBurst *sbur, superCmpEvent *sevt, int ivtx, double vertexTime){
	CorrectedTrack *t;
	vector<CorrectedTrack*>::iterator it;

	int extra = 0;
	int extraCond = 0;

	double tDiff;

	//For CDA
	double p1[3], p2[3], v1[3], v2[3];
	float dmin;
	float vertex[3];
	float corrdxdz, corrdydz;

	int conditions;
	int vtxIndex = 0;

	if(dataOnly) conditions=5;
	else conditions=4;

	if(DEBUG_5 || optDebug) cout << "Vertex tracks are :\t\t" << sevt->vtx[ivtx].vtxtrack[0].iTrack << " " << sevt->vtx[ivtx].vtxtrack[1].iTrack << " " << sevt->vtx[ivtx].vtxtrack[2].iTrack << endl;

	for(int i=0; i<vtrack.size(); i++){
		t = vtrack[i];
		extraCond = 0;

		//Blue field correction
		//vertex[0] = sevt->vtx[ivtx].x;
		//vertex[1] = sevt->vtx[ivtx].y;
		//vertex[2] = sevt->vtx[ivtx].z;
		//blue_1trk(vertex, i+1, &corrdxdz, &corrdydz, sevt);

		//t->momentum = TVector3(corrdxdz, corrdydz, 1.).Unit();
		//t->corrdxdz = corrdxdz;
		//t->corrdydz = corrdydz;


		//Is it one of the vertex track?
		if(DEBUG_5 || optDebug) cout << "\tTrying track :\t\t" << i << "\t\t not vertex track: ++" <<  endl;

		if( (sevt->vtx[ivtx].vtxtrack[0].iTrack==i) ||
				(sevt->vtx[ivtx].vtxtrack[1].iTrack==i) ||
				(sevt->vtx[ivtx].vtxtrack[2].iTrack==i)){

			t->cda = sevt->vtx[ivtx].cda;
			t->vertex = TVector3(sevt->vtx[ivtx].x, sevt->vtx[ivtx].y, sevt->vtx[ivtx].z);


			if(sevt->vtx[ivtx].vtxtrack[0].iTrack==i) vtxIndex=0;
			else if(sevt->vtx[ivtx].vtxtrack[1].iTrack==i) vtxIndex=1;
			else if(sevt->vtx[ivtx].vtxtrack[2].iTrack==i) vtxIndex=2;
			t->corrdxdz = sevt->vtx[ivtx].vtxtrack[vtxIndex].bdxdz;
			t->corrdydz = sevt->vtx[ivtx].vtxtrack[vtxIndex].bdydz;
			t->momentum = TVector3(t->corrdxdz, t->corrdydz, 1).Unit();

			//cout << "Blue field \t" << corrdxdz << "\t " << corrdydz << endl;
			//cout << "Vertex \t" << t->corrdxdz << "\t " << t->corrdydz << endl;

			goodTracks.push_back(t);
		}
		else{
			extraCond++;
		}

		// |t-t_vtx|<20ns
		if(dataOnly){
			tDiff = fabs(sevt->track[t->trackID].time - vertexTime);
			if(DEBUG_5 || optDebug) cout << "\t\t|t-t_vtx| :\t" << tDiff << "\t\t <20: ++" << endl;
			if(tDiff<20) extraCond++;
		}

		// p<60GeV/c
		if(DEBUG_5 || optDebug) cout << "\t\tp :\t\t" << t->pMag << "\t\t <74: ++" << endl;
		if(DEBUG_5 || optDebug) cout << "\t\tq :\t\t" << sevt->track[t->trackID].q << "\t\t <74: ++" << endl;
		if(t->pMag<74) extraCond++;


		//TODO Do I correct for Blue field here
		// cda to z axis
		//wrt K axis
		/*p1[0] = abcog_params.pkxoffp;
		p1[1] = abcog_params.pkyoffp;
		p1[2] = 0;
		v1[0] = abcog_params.pkdxdzp;
		v1[1] = abcog_params.pkdydzp;
		v1[2] = 1;*/
		//wrt to z axis
		p1[0] = 0;
		p1[1] = 0;
		p1[2] = 0;
		v1[0] = 0;
		v1[1] = 0;
		v1[2] = 1;

		p2[0] = t->detBPos.X();
		p2[1] = t->detBPos.Y();
		p2[2] = Geom->DCH.bz;
		v2[0] = sevt->track[t->trackID].bdxdz;
		v2[1] = sevt->track[t->trackID].bdydz;
		v2[2] = 1;
		closestApproach(p1, p2, v1, v2, &dmin, vertex);
		// -2000<z<9000
		if(DEBUG_5 || optDebug) cout << "\t\tZ_vtx :\t\t" << vertex[2] << "\t\t >-2000 || <9000: ++" << endl;
		if(vertex[2]>-2000 && vertex[2]<9000) extraCond++;

		// cda<10cm
		if(DEBUG_5 || optDebug) cout << "\t\tCDA :\t\t" << dmin << "\t\t <10: ++" << endl;
		if(dmin<10) extraCond++;

		if(DEBUG_5 || optDebug) cout << "\tExtraCond :\t\t" << extraCond << "\t ==" << conditions << ": extra" << endl;
		if(extraCond==conditions) extra++;


		// LKr Acceptance

	}

	return extra;
}

int pi0d_tracksAcceptance(superBurst *sbur){
	bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(Geom->Dch[0].PosChamber.x,Geom->Dch[0].PosChamber.y,Geom->Dch[0].PosChamber.z);
	TVector3 dch2(Geom->Dch[1].PosChamber.x,Geom->Dch[1].PosChamber.y,Geom->Dch[1].PosChamber.z);
	TVector3 dch4(Geom->Dch[3].PosChamber.x,Geom->Dch[3].PosChamber.y,Geom->Dch[3].PosChamber.z);

	TH2D *DCH1Pos = (TH2D*)gFile->Get("DCH1Pos");
	TH1D *DCH1Rad = (TH1D*)gFile->Get("DCH1Rad");
	TH2D *LKrPos = (TH2D*)gFile->Get("LKrPos");

	DCH1Rad->Reset("M");

	for(int i=0; i<goodTracks.size(); ++i){
		propPos = propagateAfter(Geom->Lkr.z, goodTracks[i]);
		lkrAcceptance = LKr_acc(sbur->nrun, propPos.X(), propPos.Y(), 8);
		LKrPos->Fill(propPos.X(), propPos.Y(),fullEvent.weight);
		if(DEBUG_7 || optDebug) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(pbWall){
			if(DEBUG_8 || optDebug) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}


		propPos = propagateCorrBefore(Geom->Dch[0].PosChamber.z, goodTracks[i]);
		radius = distance2D(dch1, propPos);
		DCH1Pos->Fill(propPos.X(), propPos.Y(), fullEvent.weight);
		DCH1Rad->Fill(radius, fullEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 120 : rejected" << endl;
		if(radius<12 || radius>120) badTrack = true;

		propPos = propagateCorrBefore(Geom->Dch[1].PosChamber.z, goodTracks[i]);
		radius = distance2D(dch2, propPos);
		if(DEBUG_7 || optDebug) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 120 : rejected" << endl;
		if(radius<12 || radius>120) badTrack = true;

		propPos = propagateAfter(Geom->Dch[3].PosChamber.z, goodTracks[i]);
		radius = distance2D(dch4, propPos);
		if(DEBUG_7 || optDebug) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 120 : rejected" << endl;
		if(radius<12 || radius>120) badTrack = true;
	}

	return badTrack;
}

int pi0d_trackCombinationVeto(superBurst *sbur, superCmpEvent *sevt){
	int ntracks = 3;
	CorrectedTrack *t1, *t2;

	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;
	int totalCharge = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			bad = false;
			if(DEBUG_8 || optDebug) cout << "\tTrying combination :\t" << i << " " << j << endl;

			t1 = goodTracks[i];
			t2 = goodTracks[j];

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(Geom->Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(Geom->Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(DEBUG_8 || optDebug) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <1: rejected" << endl;
			if(RDCH1<=1) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(dataOnly){
				tDiff = fabs(sevt->track[t1->trackID].time - sevt->track[t2->trackID].time);
				if(DEBUG_8 || optDebug) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(Geom->Lkr.z, t1);
			propPos2 = propagateAfter(Geom->Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(DEBUG_8 || optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <15: rejected" << endl;
			if(RLKr<=15) bad = true;

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_identifyPi(superBurst *sbur, superCmpEvent *sevt, int &piCandidate, bool &badElectron){
	double eop;

	int piCandidatesNb = 0;
	badElectron = false;

	for(int i=0; i<goodTracks.size(); i++){
		if(DEBUG_9 || optDebug) cout << "\tTrying track :\t\t" << i <<  endl;

		eop = goodTracks[i]->clusterEnergy/goodTracks[i]->pMag;

		if(DEBUG_9 || optDebug) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
		if(eop<0.85){
			piCandidatesNb++;
			piCandidate = i;
		}
		else if(eop>1.15) badElectron = true;
	}

	return piCandidatesNb;
}

int pi0d_goodClusters(superBurst *sbur, superCmpEvent *sevt, int piTrack, int epTrack, int emTrack, int ivtx, int &goodClusterID, double vertexTime){
	CorrectedCluster *c;
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(dataOnly) conditions = 6;
	else conditions=5;
	//if(mcOnly) conditions++;

	if(DEBUG_12 || optDebug) cout << "\tNumber of vclusters :\t" << vCluster.size() << endl;
	if(DEBUG_12 || optDebug) cout << "\tNumber of clusters :\t" << sevt->Ncluster << endl;

	for(int i=0; i<vCluster.size(); i++){
		c = vCluster[i];
		cond = 0;

		if(DEBUG_12 || optDebug) cout << "\tTrying cluster :\t" << i << endl;

		if(pbWall){
			if(DEBUG_12 || optDebug) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c->position.Y() << " < -11.850 : reject" << endl;
			if(c->position.Y()>-33.575 && c->position.Y() < -11.850) continue;
		}
		if(DEBUG_12 || optDebug) cout << "\tEnergy :\t" << c->energy << endl;
		//cout << printVector3(c->position) << endl;
		//cout << printVector3(goodTracks[epTrack].unMomentum) << endl;
		//cout << printVector3(goodTracks[epTrack].momentum) << endl;
		//cout << printVector3(goodTracks[epTrack].detPos) << endl;
		//cout << printVector3(goodTracks[epTrack].detBPos) << endl;
		//cout << printVector3(goodTracks[epTrack].vertex) << endl;
		//cout << printVector3(goodTracks[epTrack].middlePos) << endl;
		//cout << goodTracks[epTrack].t.x << " " << goodTracks[epTrack].t.y << endl;
		//cout << goodTracks[epTrack].t.dxdz << " " << goodTracks[epTrack].t.dydz << endl;
		//double x = goodTracks[epTrack].t.dxdz*(Geom->Lkr.z-Geom->DCH.z) + goodTracks[epTrack].t.x;
		//double y = goodTracks[epTrack].t.dydz*(Geom->Lkr.z-Geom->DCH.z) + goodTracks[epTrack].t.y;
		//cout << x << " " << y << endl;
		//cout << printVector3(goodTracks[emTrack].unMomentum) << endl;
		//cout << printVector3(goodTracks[emTrack].momentum) << endl;
		// separation from pi impact point >30cm
		propPos = propagateAfter(Geom->Lkr.z, goodTracks[piTrack]);
		distance = distance2D(propPos, c->position);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 30 : ++" << endl;
		if(distance>30) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(Geom->Lkr.z, goodTracks[epTrack]);
		//cout << printVector3(propPos) << endl;
		distance = distance2D(propPos, c->position);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 10 : ++" << endl;
		if(distance>10) cond++;

		propPos = propagateAfter(Geom->Lkr.z, goodTracks[emTrack]);
		distance = distance2D(propPos, c->position);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 10 : ++" << endl;
		if(distance>10) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[epTrack]);
		propPos = propagate(Geom->Lkr.z, goodTracks[epTrack]->detBPos, goodTracks[epTrack]->unBMomentum);
		distance = distance2D(propPos, c->position);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[emTrack]);
		propPos = propagate(Geom->Lkr.z, goodTracks[emTrack]->detBPos, goodTracks[emTrack]->unBMomentum);
		distance = distance2D(propPos, c->position);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// |t_g - t_vtx|<10ns
		if(dataOnly){
			tDiff = fabs(sevt->cluster[c->clusterID].time - vertexTime);
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

int pi0d_failCut(superBurst *sbur, superCmpEvent *sevt, int i){
	TH1D *fails = (TH1D*)gDirectory->Get("Cuts");

	fails->Fill(i, 1);
	cutsWord[i] = false;
	if(DEBUG_ || optDebug) cout << "Event is not passing selection" << endl;
	if(!noOutput) fprintf(fprt, "%i %i %i %i\n", sbur->nrun, sbur->time, sevt->timeStamp, i);
	return 0;
}
void pi0d_passSelection(superBurst *sbur, superCmpEvent *sevt){
	if(DEBUG_ || optDebug) cout << "Event is passing selection" << endl;
	if(!noOutput) fprintf(fprt2, "%i %i %i\n", sbur->nrun, sbur->time, sevt->timeStamp);
}

int nico_pi0DalitzFilter(superBurst *sbur,superCmpEvent *sevt){
	int i = 0;
	double eop;
	int eTracks = 0;
	CorrectedTrack *t;

	// 1) Trigger Q2xMBX(1VTX||2VTX||1TRK)
	if((sevt->trigWord & 0x400) ==0) return 1;

	// 2) 1 vertex that is a 3-track vertex
	if(sevt->Nvtx!=1) return 1;
	if(sevt->vtx[i].Nvtxtrack!=3) return 1;

	// 9) 2 tracks with E/p > 0.5
	t = vtrack[sevt->vtx[i].vtxtrack[0].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;

	t = vtrack[sevt->vtx[i].vtxtrack[1].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;
	t = vtrack[sevt->vtx[i].vtxtrack[2].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;

	if(eTracks<2) return 1;

	return -1;
}

int nico_pi0DalitzSelect(superBurst *sbur,superCmpEvent *sevt){
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

	fullEvent.runNumber = sbur->nrun;
	fullEvent.burstNumber = sbur->time;
	fullEvent.timestamp = sevt->timeStamp;

	if(DEBUG_ || optDebug) cout << endl;

	memset(&cutsWord, true, sizeof(bool)*19);
	// 1) Trigger Q2xMBX(1VTX||2VTX||1TRK)
	int triggerMask;
	if(sbur->nrun<20209) triggerMask = (0x800 | 0x400);
	else triggerMask = 0x400;

	if(DEBUG_1 || optDebug) cout << "~~~~ Cut 1 ~~~~" << endl;
	bitset<32> triggBit(sevt->trigWord);
	bitset<32> triggMask(triggerMask);
	bitset<32> comb(sevt->trigWord & triggerMask);
	if(DEBUG_1 || optDebug) cout << "Word :\t " << triggBit << endl;
	if(DEBUG_1 || optDebug) cout << "Mask :\t " << triggMask << endl;
	if(DEBUG_1 || optDebug) cout << "comb :\t " << comb << endl;
	if(DEBUG_1 || optDebug) cout << "L2 Trigger word :\t " << hex << (sevt->trigWord & triggerMask) << dec << "\t == 0: rejected" << endl;
	if(dataOnly && ((sevt->trigWord & triggerMask) ==0)) FAILCUT(1)

	// 2) Exactly one 3-track vertex with correct charge
	if(DEBUG_2 || optDebug) cout << "~~~~ Cut 2 ~~~~" << endl;
	vtxNb = pi0d_countVtx3Tracks(sbur,sevt,ivtx);
	if(DEBUG_2 || optDebug) cout << "3-track vertex :\t\t" << vtxNb << "\t != 1 : rejected" << endl;
	if(vtxNb!=1) FAILCUT(2)

	fullEvent.zVtx = sevt->vtx[ivtx].z;

	if(sevt->vtx[ivtx].charge==1){
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
	fullEvent.weight = weight;

	vertexTime = pi0d_getVertexTime(sevt, ivtx);

	// 3) Vertex position -18m<Z_vtx<90m
	if(DEBUG_3 || optDebug) cout << "~~~~ Cut 3 ~~~~" << endl;
	if(DEBUG_3 || optDebug) cout << "vertex Z :\t\t\t" << sevt->vtx[ivtx].z << "\t <-1800 || >9000 : rejected" << endl;
	if(sevt->vtx[ivtx].z<=-1800 || sevt->vtx[ivtx].z>=9000) FAILCUT(3)

	// 4) Vertex quality chi2<25
	if(DEBUG_4 || optDebug) cout << "~~~~ Cut 4 ~~~~" << endl;
	if(DEBUG_4 || optDebug) cout << "vertex chi2 :\t\t\t" << sevt->vtx[ivtx].chi2 << "\t\t > 25: rejected" << endl;
	if(sevt->vtx[ivtx].chi2>=25) FAILCUT(4)

	// 5) No extra track
	if(DEBUG_5 || optDebug) cout << "~~~~ Cut 5 ~~~~" << endl;
	extraTracks = pi0d_extraTrackVeto(sbur,sevt,ivtx, vertexTime);
	if(DEBUG_5 || optDebug) cout << "Extra tracks :\t\t\t" << extraTracks << "\t > 0 : rejected" << endl;
	if(extraTracks>0) FAILCUT(5)

	// 6) Track DCH time
	if(DEBUG_6 || optDebug) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(dataOnly){
		if(DEBUG_6 || optDebug) cout << "|t_1| :\t\t\t\t" << fabs(sevt->track[goodTracks[0]->trackID].time - sbur->tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(sevt->track[goodTracks[0]->trackID].time - sbur->tOffst.Dch)>=25) FAILCUT(6)
				if(DEBUG_6 || optDebug) cout << "|t_2| :\t\t\t\t" << fabs(sevt->track[goodTracks[1]->trackID].time - sbur->tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(sevt->track[goodTracks[1]->trackID].time - sbur->tOffst.Dch)>=25) FAILCUT(6)
				if(DEBUG_6 || optDebug) cout << "|t_3| :\t\t\t\t" << fabs(sevt->track[goodTracks[2]->trackID].time - sbur->tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(sevt->track[goodTracks[2]->trackID].time - sbur->tOffst.Dch)>=25) FAILCUT(6)
	}
	else{
		if(DEBUG_6 || optDebug) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(DEBUG_7 || optDebug) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance(sbur);
	if(DEBUG_7 || optDebug) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==true) FAILCUT(7)

	// 8) Track combination veto
	if(DEBUG_8 || optDebug) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto(sbur,sevt);
	if(DEBUG_8 || optDebug) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=0) FAILCUT(8)

	// 9) Identify candidates
	if(DEBUG_9 || optDebug) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pi0d_identifyPi(sbur, sevt, piTrack, badElectron);
	if(DEBUG_9 || optDebug) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=1) FAILCUT(9)

	for(int i=0; i<goodTracks.size(); i++){
		if(i==piTrack) continue;

		if(sevt->track[goodTracks[i]->trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(sevt->track[goodTracks[i]->trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}
	//if(piTrack==-1 || epTrack==-1 || emTrack==-1)
		//FILCUT(9)



	fullEvent.em = goodTracks[emTrack];
	fullEvent.ep = goodTracks[epTrack];
	fullEvent.pip= goodTracks[piTrack];

	// 10) Bad electron cluster
	if(DEBUG_10 || optDebug) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(DEBUG_10 || optDebug) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==true) FAILCUT(10)

	// 11) Tracks momenta
	if(DEBUG_11 || optDebug) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(DEBUG_11 || optDebug) cout << "p_pi :\t\t\t\t" << goodTracks[piTrack]->pMag << "\t <5 || > 60 : rejected" << endl;
	if(DEBUG_11 || optDebug) cout << "p_e+ :\t\t\t\t" << goodTracks[epTrack]->pMag << "\t <5 || > 60 : rejected" << endl;
	if(DEBUG_11 || optDebug) cout << "p_e- :\t\t\t\t" << goodTracks[emTrack]->pMag << "\t <5 || > 60 : rejected" << endl;
	if(goodTracks[piTrack]->pMag<=5 || goodTracks[piTrack]->pMag>=60) FAILCUT(11)
	if(goodTracks[epTrack]->pMag<=5 || goodTracks[epTrack]->pMag>=60) FAILCUT(11)
	if(goodTracks[emTrack]->pMag<=5 || goodTracks[emTrack]->pMag>=60) FAILCUT(11)

	// 12) Exactly 1 good LKr cluster
	if(DEBUG_12 || optDebug) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters(sbur, sevt, piTrack, epTrack, emTrack, ivtx, goodClusterID, vertexTime);
	if(DEBUG_12 || optDebug) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=1) FAILCUT(12)

	fullEvent.gamma = vCluster[goodClusterID];

	if(goodClusterID==-1){ //Automatically fails all subsequent cuts
		//pi0d_failCut(sbur,sevt,12);
		//pi0d_failCut(sbur,sevt,13);
		//pi0d_failCut(sbur,sevt,14);
		//pi0d_failCut(sbur,sevt,15);
		//pi0d_failCut(sbur,sevt,16);
		//pi0d_failCut(sbur,sevt,17);
		//pi0d_failCut(sbur,sevt,18);
		//pi0d_failCut(sbur,sevt,19);
		//outTree->Fill();
		return 0;
	}
	fullEvent.pGamma = (vCluster[goodClusterID]->position - goodTracks[piTrack]->vertex).Unit();

	// 13) Photon candidate in LKr acceptance
	//TODO check it's the good way
	if(DEBUG_13 || optDebug) cout << "~~~~ Cut 13 ~~~~" << endl;
	//propPos = propagate(Geom->Lkr.z, vCluster[goodClusterID]->position, photonMomentum);
	//lkrAcceptance = LKr_acc(sbur->nrun, vCluster[goodClusterID]->position.X(), vCluster[goodClusterID]->position.Y(), 8);
	lkrAcceptance = LKr_acc(sbur->nrun, sevt->cluster[vCluster[goodClusterID]->clusterID].x, sevt->cluster[vCluster[goodClusterID]->clusterID].y, 8);
	//lkrAcceptance = LKr_acc(sbur->nrun, propPos.X(), propPos.Y(), 8);
	if(DEBUG_13 || optDebug) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=0) FAILCUT(13)

	// 14) E_gamma>3GeV
	if(DEBUG_14 || optDebug) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(DEBUG_14 || optDebug) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << vCluster[goodClusterID]->energy << "\t < 3 : rejected" << endl;
	if(vCluster[goodClusterID]->energy<=3) FAILCUT(14)

	// 15) D_deadcell>2cm
	if(DEBUG_15 || optDebug) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(DEBUG_15 || optDebug) cout << "d_deadcell :\t\t\t" << sevt->cluster[vCluster[goodClusterID]->clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(sevt->cluster[vCluster[goodClusterID]->clusterID].dDeadCell<=2) FAILCUT(15)

	//cout << printVector3(goodTracks[piTrack].vertex) << endl;
	//cout << printVector3(vCluster[goodClusterID]->position) << endl;
	//cout << printVector3(vCluster[goodClusterID]->energy*photonMomentum) << endl;

	// 16) Photon DCH1 intercept >13cm
	if(DEBUG_16 || optDebug) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(Geom->Dch[0].PosChamber.z, vCluster[goodClusterID]->position, fullEvent.pGamma);
	//radius = sqrt(pow(propPos.X()-Geom->Dch[0].PosChamber.x,2) + pow(propPos.Y()-Geom->Dch[0].PosChamber.y,2));
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(DEBUG_16 || optDebug) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=13) FAILCUT(16)

	// 17) Total momentum 70<p<78
	if(DEBUG_17 || optDebug) cout << "~~~~ Cut 17 ~~~~" << endl;

	fullEvent.pTotal = goodTracks[piTrack]->pMag*goodTracks[piTrack]->momentum +
			goodTracks[emTrack]->pMag*goodTracks[emTrack]->momentum +
			goodTracks[epTrack]->pMag*goodTracks[epTrack]->momentum +
			vCluster[goodClusterID]->energy*fullEvent.pGamma;

	if(DEBUG_17 || optDebug) cout << "p_pieeg :\t\t\t" << fullEvent.pTotal.Mag() << "\t <70 || >78 : rejected" << endl;
	if(fullEvent.pTotal.Mag()<70 || fullEvent.pTotal.Mag()>78) FAILCUT(17)

	// 18) Transverse momentum^2 < 5E-4
	if(DEBUG_18 || optDebug) cout << "~~~~ Cut 18 ~~~~" << endl;
	pt = fullEvent.pTotal.Perp2(kaonMomentum);
	//pt = totalMomentum.Perp2(TVector3(0,0,1));
	if(DEBUG_18 || optDebug) cout << "P_t^2 :\t\t" << pt << "\t >= 0.0005 : rejected" << endl;
	if(pt>=5e-4) FAILCUT(18)

	// 19) |M_eeg - M_pi0|<8 MeV
	if(DEBUG_19 || optDebug) cout << "~~~~ Cut 19 ~~~~" << endl;
	vMass.push_back(Me);
	vP.push_back(goodTracks[epTrack]->pMag*goodTracks[epTrack]->momentum);
	vMass.push_back(Me);
	vP.push_back(goodTracks[emTrack]->pMag*goodTracks[emTrack]->momentum);
	fullEvent.mee = sqrt(invMass2(vMass, vP));
	if(DEBUG_19 || optDebug) cout << "M_ee :\t\t" << fullEvent.mee << endl;
	vMass.push_back(0.0);
	vP.push_back(vCluster[goodClusterID]->energy*fullEvent.pGamma);

	fullEvent.mPi0 = sqrt(invMass2(vMass, vP));
	if(DEBUG_19 || optDebug) cout << "|M_eeg - M_pi0| :\t\t" << fabs(fullEvent.mPi0-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(fullEvent.mPi0-Mpi0)>=0.008) FAILCUT(19)

	// 20) 0.475 < M_pieeg < 0.510
	if(DEBUG_20 || optDebug) cout << "~~~~ Cut 20 ~~~~" << endl;
	vMass.push_back(Mpic);
	vP.push_back(goodTracks[piTrack]->pMag*goodTracks[piTrack]->momentum);
	fullEvent.mK = sqrt(invMass2(vMass, vP));

	if(DEBUG_20 || optDebug) cout << "M_pieeg :\t\t" << fullEvent.mK << "\t <0.475 || >0.510: rejected" << endl;
	if(fullEvent.mK<0.475 || fullEvent.mK>0.510) FAILCUT(20)

	outTree->FlushBaskets();
	outTree->Fill();

	pi0d_passSelection(sbur, sevt);
	return 0;
}
