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

double pi0d_getVertexTime(int ivtx){
	double mean = 0;

	for(int i=0; i<rawEvent.vtx[ivtx].Nvtxtrack; i++){
		mean += rawEvent.track[corrEvent.pTrack[rawEvent.vtx[ivtx].vtxtrack[i].iTrack].trackID].time;
	}

	return mean/(double)rawEvent.vtx[ivtx].Nvtxtrack;
}

int pi0d_countVtx3Tracks(int &ivtx){
	int vtxNb = 0;
	//int totalCharge;


	gH("vertexN")->Fill(rawEvent.Nvtx, corrEvent.weight);

	if(DEBUG_2 || optDebug) cout << "Total Number of vertices :\t" << rawEvent.Nvtx << " ==1" << endl;
	if(rawEvent.Nvtx!=1) return rawEvent.Nvtx;
	for(int i=0; i<rawEvent.Nvtx; i++){
		if(DEBUG_2 || optDebug) cout << "\tTrying vertex :\t\t" << i << endl;
		if(DEBUG_2 || optDebug) cout << "\tNumber of tracks:\t" << rawEvent.vtx[i].Nvtxtrack << endl;
		gH("trackN")->Fill(rawEvent.vtx[i].Nvtxtrack, corrEvent.weight);
		if(rawEvent.vtx[i].Nvtxtrack==3){
			if(DEBUG_2 || optDebug) cout << "\tTotal charge :\t" << rawEvent.vtx[i].charge << endl;
			if(rawEvent.vtx[i].charge==beamCharge || noTestCharge){
				vtxNb++;
				ivtx = i;
			}
		}
	}

	return vtxNb;
}

int pi0d_extraTrackVeto(int ivtx, double vertexTime){

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

	if(DEBUG_5 || optDebug) cout << "Vertex tracks are :\t\t" << rawEvent.vtx[ivtx].vtxtrack[0].iTrack << " " << rawEvent.vtx[ivtx].vtxtrack[1].iTrack << " " << rawEvent.vtx[ivtx].vtxtrack[2].iTrack << endl;

	for(int i=0; i<corrEvent.pTrack.size(); i++){
		NPhysicsTrack &pt = corrEvent.pTrack[i];
		NTrak &t = rawEvent.track[pt.trackID];
		extraCond = 0;

		//Is it one of the vertex track?
		if(DEBUG_5 || optDebug) cout << "\tTrying track :\t\t" << i << "\t\t not vertex track: ++" <<  endl;

		if( (rawEvent.vtx[ivtx].vtxtrack[0].iTrack==i) ||
				(rawEvent.vtx[ivtx].vtxtrack[1].iTrack==i) ||
				(rawEvent.vtx[ivtx].vtxtrack[2].iTrack==i)){

			t.vtxID = ivtx;

			if(rawEvent.vtx[ivtx].vtxtrack[0].iTrack==i) vtxIndex=0;
			else if(rawEvent.vtx[ivtx].vtxtrack[1].iTrack==i) vtxIndex=1;
			else if(rawEvent.vtx[ivtx].vtxtrack[2].iTrack==i) vtxIndex=2;

			pt.momentum = TVector3(rawEvent.vtx[ivtx].vtxtrack[vtxIndex].bdxdz, rawEvent.vtx[ivtx].vtxtrack[vtxIndex].bdydz, 1).Unit();

			corrEvent.goodTracks.push_back(i);
		}
		else{
			extraCond++;
		}

		// |t-t_vtx|<20ns
		if(dataOnly){
			gH("trackVertexTime")->Fill(rawEvent.track[pt.trackID].time - vertexTime, corrEvent.weight);
			tDiff = fabs(rawEvent.track[pt.trackID].time - vertexTime);
			if(DEBUG_5 || optDebug) cout << "\t\t|t-t_vtx| :\t" << tDiff << "\t\t <20: ++" << endl;
			if(tDiff<20) extraCond++;
		}

		// p<74GeV/c
		gH("trackP")->Fill(pt.p, corrEvent.weight);
		if(DEBUG_5 || optDebug) cout << "\t\tp :\t\t" << pt.p << "\t\t <74: ++" << endl;
		if(DEBUG_5 || optDebug) cout << "\t\tq :\t\t" << rawEvent.track[pt.trackID].q << "\t\t <74: ++" << endl;
		if(pt.p<74) extraCond++;

		//wrt to z axis
		p1[0] = 0;
		p1[1] = 0;
		p1[2] = 0;
		v1[0] = 0;
		v1[1] = 0;
		v1[2] = 1;

		p2[0] = t.bDetPos.X();
		p2[1] = t.bDetPos.Y();
		p2[2] = Geom->DCH.bz;
		v2[0] = rawEvent.track[pt.trackID].bdxdz;
		v2[1] = rawEvent.track[pt.trackID].bdydz;
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
	}
	return extra;
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
		lkrAcceptance = LKr_acc(rootBurst.nrun, propPos.X(), propPos.Y(), 8);
		gH("trackLKrRadius")->Fill(distance2D(TVector3(0,0,0), propPos), corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(pbWall){
			if(DEBUG_8 || optDebug) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}


		propPos = propagateCorrBefore(Geom->Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		gH("trackDCH1Radius")->Fill(radius, corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(Geom->Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		gH("trackDCH2Radius")->Fill(radius, corrEvent.weight);
		if(DEBUG_7 || optDebug) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(Geom->Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		gH("trackDCH4Radius")->Fill(radius, corrEvent.weight);
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
			gH("track_ij_DCH1")->Fill(RDCH1, corrEvent.weight);
			if(DEBUG_8 || optDebug) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(dataOnly){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				gH("trackTrackTime")->Fill(tDiff, corrEvent.weight);
				if(DEBUG_8 || optDebug) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			propPos1 = propagateAfter(Geom->Lkr.z, t1);
			propPos2 = propagateAfter(Geom->Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(i==piTrack || j==piTrack){
				gH("track_epi_LKr")->Fill(RLKr, corrEvent.weight);
				if(DEBUG_8 || optDebug) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				gH("track_ee_LKr")->Fill(RLKr, corrEvent.weight);
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

		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		gH("trackEOP")->Fill(eop, corrEvent.weight);
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

	gH("clusterN")->Fill(corrEvent.pCluster.size(), corrEvent.weight);
	for(int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		gH("clusterX")->Fill(c.position.X(), corrEvent.weight);
		gH("clusterY")->Fill(c.position.Y(), corrEvent.weight);
		gH("clusterRadius")->Fill(distance2D(TVector3(0,0,0), c.position), corrEvent.weight);

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
		gH("clusterDPi")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(Geom->Lkr.z, ep);
		//cout << printVector3(propPos) << endl;
		distance = distance2D(propPos, c.position);
		gH("clusterDep")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(Geom->Lkr.z, em);
		distance = distance2D(propPos, c.position);
		gH("clusterDem")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[epTrack]);
		propPos = propagate(Geom->Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		gH("clusterDUnep")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		//propPos = propagateCorrBefore(c->position.Z(), goodTracks[emTrack]);
		propPos = propagate(Geom->Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		gH("clusterDUnem")->Fill(distance, corrEvent.weight);
		if(DEBUG_12 || optDebug) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// |t_g - t_vtx|<10ns
		if(dataOnly){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - vertexTime);
			gH("clusterVertexTime")->Fill(tDiff, corrEvent.weight);
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
	TH1D *fails = (TH1D*)gDirectory->Get("Cuts");

	fails->Fill(i, 1);
	cutsWord[i] = false;
	if(DEBUG_ || optDebug) cout << "Event is not passing selection" << endl;
	if(!noOutput) fprintf(fprt, "%i %i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, i);
	return 0;
}
void pi0d_passSelection(){
	if(DEBUG_ || optDebug) cout << "Event is passing selection" << endl;
	if(!noOutput) fprintf(fprt2, "%i %i %i\n", rootBurst.nrun, rootBurst.time, rawEvent.timeStamp);
}

int nico_pi0DalitzFilter(){
/*	int i = 0;
	double eop;
	int eTracks = 0;
	CorrectedTrack *t;

	// 1) Trigger Q2xMBX(1VTX||2VTX||1TRK)
	if((corrEvent.trigWord & 0x400) ==0) return 1;

	// 2) 1 vertex that is a 3-track vertex
	if(corrEvent.Nvtx!=1) return 1;
	if(corrEvent.vtx[i].Nvtxtrack!=3) return 1;

	// 9) 2 tracks with E/p > 0.5
	t = vtrack[corrEvent.vtx[i].vtxtrack[0].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;

	t = vtrack[corrEvent.vtx[i].vtxtrack[1].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;
	t = vtrack[corrEvent.vtx[i].vtxtrack[2].iTrack];
	eop = t->clusterEnergy/t->pMag;
	if(eop>0.5) eTracks++;

	if(eTracks<2) return 1;
*/
	return -1;
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

	if(DEBUG_ || optDebug) cout << endl;

	memset(&cutsWord, true, sizeof(bool)*19);
	// 1) Trigger Q2xMBX(1VTX||2VTX||1TRK)
	int triggerMask;
	if(rootBurst.nrun<20209) triggerMask = 0x800;
	else triggerMask = cutsDefinition.triggerMask;

	if(DEBUG_1 || optDebug) cout << "~~~~ Cut 1 ~~~~" << endl;
	bitset<32> triggBit(rawEvent.trigWord);
	bitset<32> triggMask(triggerMask);
	bitset<32> comb(rawEvent.trigWord & triggerMask);
	if(DEBUG_1 || optDebug) cout << "Word :\t " << triggBit << endl;
	if(DEBUG_1 || optDebug) cout << "Mask :\t " << triggMask << endl;
	if(DEBUG_1 || optDebug) cout << "comb :\t " << comb << endl;
	if(DEBUG_1 || optDebug) cout << "L2 Trigger word :\t " << hex << (rawEvent.trigWord & triggerMask) << dec << "\t == 0: rejected" << endl;
	if(dataOnly && ((rawEvent.trigWord & triggerMask) ==0)) FAILCUT(1)

	// 2) Exactly one 3-track vertex with correct charge
	if(DEBUG_2 || optDebug) cout << "~~~~ Cut 2 ~~~~" << endl;
	vtxNb = pi0d_countVtx3Tracks(ivtx);
	if(DEBUG_2 || optDebug) cout << "3-track vertex :\t\t" << vtxNb << "\t != 1 : rejected" << endl;
	if(vtxNb!=cutsDefinition.numVertex3) FAILCUT(2)

	//corrEvent.zVtx = corrEvent.vtx[ivtx].z;

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

	vertexTime = pi0d_getVertexTime(ivtx);

	gH("vertexZ")->Fill(rawEvent.vtx[ivtx].position.Z(), corrEvent.weight);
	gH("vertexChi2")->Fill(rawEvent.vtx[ivtx].chi2, corrEvent.weight);
	gH("vertexQ")->Fill(rawEvent.vtx[ivtx].charge, corrEvent.weight);

	// 3) Vertex position -18m<Z_vtx<90m
	if(DEBUG_3 || optDebug) cout << "~~~~ Cut 3 ~~~~" << endl;
	if(DEBUG_3 || optDebug) cout << "vertex Z :\t\t\t" << rawEvent.vtx[ivtx].position.Z() << "\t <-1800 || >9000 : rejected" << endl;
	if(rawEvent.vtx[ivtx].position.Z()<=cutsDefinition.minZVertex || rawEvent.vtx[ivtx].position.Z()>=cutsDefinition.maxZVertex) FAILCUT(3)


	// 4) Vertex quality chi2<25
	if(DEBUG_4 || optDebug) cout << "~~~~ Cut 4 ~~~~" << endl;
	if(DEBUG_4 || optDebug) cout << "vertex chi2 :\t\t\t" << rawEvent.vtx[ivtx].chi2 << "\t\t > 25: rejected" << endl;
	if(rawEvent.vtx[ivtx].chi2>=cutsDefinition.maxChi2Vertex) FAILCUT(4)

	// 5) No extra track
	if(DEBUG_5 || optDebug) cout << "~~~~ Cut 5 ~~~~" << endl;
	extraTracks = pi0d_extraTrackVeto(ivtx, vertexTime);
	if(DEBUG_5 || optDebug) cout << "Extra tracks :\t\t\t" << extraTracks << "\t > 0 : rejected" << endl;
	if(extraTracks>cutsDefinition.maxExtraTracks) FAILCUT(5)

	if(dataOnly){
		gH("trackDCHTime")->Fill(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch), corrEvent.weight);
		gH("trackDCHTime")->Fill(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch), corrEvent.weight);
		gH("trackDCHTime")->Fill(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch), corrEvent.weight);
	}

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
	gH("clusterGoodN")->Fill(goodClusters, corrEvent.weight);
	if(DEBUG_12 || optDebug) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=cutsDefinition.numAddGoodCluster) FAILCUT(12)

	//fullEvent.gamma = vCluster[goodClusterID];

	if(goodClusterID==-1){
		return 0;
	}

	TVector3 pGamma = (corrEvent.pCluster[goodClusterID].position - rawEvent.vtx[rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[piTrack]].trackID].vtxID].position).Unit();
	//fullEvent.pGamma = (vCluster[goodClusterID]->position - goodTracks[piTrack]->vertex).Unit();
	gH("gammaEnergy")->Fill(corrEvent.pCluster[goodClusterID].E, corrEvent.weight);
	gH("gammaX")->Fill(corrEvent.pCluster[goodClusterID].position.X(), corrEvent.weight);
	gH("gammaY")->Fill(corrEvent.pCluster[goodClusterID].position.Y(), corrEvent.weight);
	gH("gammaRadius")->Fill(distance2D(TVector3(0,0,0), corrEvent.pCluster[goodClusterID].position), corrEvent.weight);
	gH("gammaDeadCell")->Fill(rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].dDeadCell, corrEvent.weight);


	// 13) Photon candidate in LKr acceptance
	//TODO check it's the good way
	if(DEBUG_13 || optDebug) cout << "~~~~ Cut 13 ~~~~" << endl;
	//propPos = propagate(Geom->Lkr.z, vCluster[goodClusterID]->position, photonMomentum);
	//lkrAcceptance = LKr_acc(sbur->nrun, vCluster[goodClusterID]->position.X(), vCluster[goodClusterID]->position.Y(), 8);
	lkrAcceptance = LKr_acc(rootBurst.nrun, rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].position.X(), rawEvent.cluster[corrEvent.pCluster[goodClusterID].clusterID].position.Y(), 8);
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
	gH("gammaDCH1Radius")->Fill(radius, corrEvent.weight);
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

	gH("bkinKP")->Fill(kaonMomentum.Mag(), corrEvent.weight);
	gH("bkinPiP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p, corrEvent.weight);
	gH("bkinemP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p, corrEvent.weight);
	gH("bkinepP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p, corrEvent.weight);
	gH("bkinGammaP")->Fill((corrEvent.pCluster[goodClusterID].E*pGamma).Mag(), corrEvent.weight);

	gH("bkinMee")->Fill(mee, corrEvent.weight);
	gH("bkinx")->Fill(mee/Mpi0, corrEvent.weight);
	gH("bkinMeeg")->Fill(mPi0, corrEvent.weight);
	gH("bkinMeegDiff")->Fill(mPi0-Mpi0, corrEvent.weight);
	gH("bkinMeegpi")->Fill(mK, corrEvent.weight);

	gH("bkinPTot")->Fill(pTotal.Mag(), corrEvent.weight);
	gH("bkinPt2")->Fill(pt, corrEvent.weight);

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
	gH("akinKP")->Fill(kaonMomentum.Mag(), corrEvent.weight);
	gH("akinPiP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[piTrack]].p, corrEvent.weight);
	gH("akinemP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[emTrack]].p, corrEvent.weight);
	gH("akinepP")->Fill(corrEvent.pTrack[corrEvent.goodTracks[epTrack]].p, corrEvent.weight);
	gH("akinGammaP")->Fill((corrEvent.pCluster[goodClusterID].E*pGamma).Mag(), corrEvent.weight);

	gH("akinMee")->Fill(mee, corrEvent.weight);
	gH("akinx")->Fill(pow(mee/Mpi0,2), corrEvent.weight);
	gH("akinMeeg")->Fill(mPi0, corrEvent.weight);
	gH("akinMeegDiff")->Fill(mPi0-Mpi0, corrEvent.weight);
	gH("akinMeegpi")->Fill(mK, corrEvent.weight);

	gH("akinPTot")->Fill(pTotal.Mag(), corrEvent.weight);
	gH("akinPt2")->Fill(pt, corrEvent.weight);
	outTree->FlushBaskets();
	outTree->Fill();

	pi0d_passSelection();
	return 0;
}
