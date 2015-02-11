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

	outTree->FlushBaskets();
	outTree->Fill();

	pi0d_passSelection();
	return 0;
}
