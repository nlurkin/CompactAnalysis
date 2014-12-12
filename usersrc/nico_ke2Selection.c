/*
 * ke2Selection.c
 *
 *  Created on: 18 Feb 2014
 *      Author: ncl
 */
#include <math.h>

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include <iostream>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <bitset>
using namespace std;

#define DEBUG_ALL false

#define DEBUG_ALL_COM DEBUG_ALL || false
#define DEBUG_1 DEBUG_ALL_COM || false
#define DEBUG_2 DEBUG_ALL_COM || false
#define DEBUG_3 DEBUG_ALL_COM || false
#define DEBUG_4 DEBUG_ALL_COM || false
#define DEBUG_5 DEBUG_ALL_COM || false
#define DEBUG_6 DEBUG_ALL_COM || false
#define DEBUG_7 DEBUG_ALL_COM || false
#define DEBUG_8 DEBUG_ALL_COM || false
#define DEBUG_9 DEBUG_ALL_COM || false
#define DEBUG_10 DEBUG_ALL_COM || false
#define DEBUG_11 DEBUG_ALL_COM || false
#define DEBUG_12 DEBUG_ALL_COM || false
#define DEBUG_13 DEBUG_ALL_COM || false
#define DEBUG_15 DEBUG_ALL_COM || false
#define DEBUG_16 DEBUG_ALL_COM || false

#define DEBUG_ALL_KE2 DEBUG_ALL || false
#define DEBUG_KE2_1 DEBUG_ALL_KE2 || false
#define DEBUG_KE2_2 DEBUG_ALL_KE2 || false
#define DEBUG_KE2_3 DEBUG_ALL_KE2 || false
#define DEBUG_KE2_4 DEBUG_ALL_KE2 || false
#define DEBUG_KE2_5 DEBUG_ALL_KE2 || false
#define DEBUG_KE2_6 DEBUG_ALL_KE2 || false

#define DEBUG_COM DEBUG_1 || DEBUG_2 || DEBUG_3 || DEBUG_4 || DEBUG_5 || DEBUG_6 || DEBUG_7 || DEBUG_8 || DEBUG_9 || DEBUG_10 ||  DEBUG_11 || DEBUG_12 || DEBUG_13 || DEBUG_15 || DEBUG_16
#define DEBUG_KE2 DEBUG_KE2_1 || DEBUG_KE2_2 || DEBUG_KE2_3 || DEBUG_KE2_4 || DEBUG_KE2_5 || DEBUG_KE2_6
#define DEBUG_ DEBUG_COM || DEBUG_KE2


void reassociateCluster(superCmpEvent *sevt){
	bool assoc = false;
	TVector3 propPos;
	double dist;
	int iClosest=-1;
	double dClosest=99999;

	for(int i=0; i<vCluster.size(); i++){
		propPos = propagateAfter(vCluster[i]->position.Z(), goodTracks[0]);
		if(sevt->cluster[vCluster[i]->clusterID].iTrack == goodTracks[0]->trackID) assoc = true;
		dist = distance2D(propPos, vCluster[i]->position);
		if(dist<dClosest){
			dClosest = dist;
			iClosest = i;
		}
	}

	if(!assoc && iClosest>-1){
		//No cluster associated to the goodTracks[0]s[0] found. Associate the closest.
		sevt->cluster[vCluster[iClosest]->clusterID].iTrack = goodTracks[0]->trackID;
	}
}

int ke2_vetoLKr(superCmpEvent *sevt){
	CorrectedCluster *clst;

	//Counters
	int cond = 0;
	int bad = 0;

	//track compare
	TVector3 propPos;
	double distance;

	if(DEBUG_5) cout << "Number of clusters : " << vCluster.size() << endl;
	for(int i=0; i<vCluster.size(); i++){
		cond = 0;
		if(DEBUG_5) cout << "Working on cluster  " << i << endl;
		clst = vCluster[i];

		propPos = propagateAfter(clst->position.Z(), goodTracks[0]);
		distance = distance2D(propPos, clst->position);
		//store clusters id when distance is <10cm for later
		if(distance<10.){
			closeClusters.push_back(vCluster[i]);
		}

		propPos = propagateBefore(Geom->Lkr.z+30, goodTracks[0]);
		distance = distance2D(propPos, clst->position);

		if(DEBUG_5) cout << "\t Cluster Energy :\t" << sevt->cluster[clst->clusterID].energy << "\t > 2: ++" << endl;
		if(sevt->cluster[clst->clusterID].energy>2) cond++; //!cluster energy <2GeV

		if(DEBUG_5) cout << "\t Cluster track :\t" << sevt->cluster[clst->clusterID].iTrack << "\t != " << goodTracks[0]->trackID << ": ++" << endl;
		if(sevt->cluster[clst->clusterID].iTrack != goodTracks[0]->trackID){
			cond++; //!cluster associated to track
		}
		else assocClusters.push_back(clst);

		if(DEBUG_5) cout << "\t Distance :\t\t" << distance << "\t > 6: ++" << endl;
		if(distance > 6) cond++; //!undeflected track to cluster distance R>6cm in Z plane 30cm downstream LKr face;

		if(DEBUG_5) cout << "\t Time :\t\t\t" << fabs(sevt->cluster[clst->clusterID].time - sevt->track[goodTracks[0]->trackID].time) << "\t < 12: ++" << endl;
		if(fabs(sevt->cluster[clst->clusterID].time - sevt->track[goodTracks[0]->trackID].time)<12) cond++; //!|Delta t|<12ns

		if(DEBUG_5) cout << "\t Conditions :\t\t" << cond << "\t == 4 : bad cluster" << endl;
		if(cond==4) bad++; //All conditions fulfilled
	}

	return bad;
}

int ke2_getGoodTracks(superBurst *sbur,superCmpEvent *sevt){
	//select exactly 1 good track from DCH
	CorrectedTrack *t;

	//Counters
	int bad = 0;
	int fBad;
	int goodTrackID = -1;
	int good=0;

	//Compare with other tracks
	TVector3 pos1,pos2;
	double rSquare;
	double y4a[2], y4b[2];

	//For CDA
	double p1[3], p2[3], v1[3], v2[3];
	float dmin;
	float vertex[3];
	float corrdxdz, corrdydz;

	for(int i=0; i<vtrack.size(); i++){
		//Select the track
		if(DEBUG_1) cout << "Working on track :\t" << i << endl;
		t = vtrack[i];

		//Apply corrections and create CorrectedObject
		//tCorr = correctTrack(sevt, t, p, e, eop);

		if(DEBUG_1) cout << "\tMomentum :\t" << t->pMag << "\t < 3 || > 75" << endl;
		if(t->pMag<3 || t->pMag>75) {bad++; continue;} // p<3GeV || p>75Gev with alpha, beta corrections
		fBad=0;

		pos1 = propagateBefore(Geom->Dch[0].PosChamber.z, t);

		//Loop over other tracks
		for(int j=0; j<vtrack.size(); j++){
			if(j==i) continue;
			if(vtrack[j]->pMag<3 || vtrack[j]->pMag>75) continue;
			if(DEBUG_1) cout << "\tComparing with track :\t" << j << endl;
			pos2 = propagateBefore(Geom->Dch[0].PosChamber.z, vtrack[j]);
			rSquare = pow(pos1.X() - pos2.X(), 2) + pow(pos1.Y() - pos2.Y(), 2);
			if(DEBUG_1) cout << "\t\tRadius^2 :\t" << rSquare << "\t < 0.25" << endl;
			if(rSquare<0.25){
				if(DEBUG_1) cout << "\t\tt1 Quality :\t" << sevt->track[t->trackID].quality << endl;
				if(DEBUG_1) cout << "\t\tt2 Quality :\t" << sevt->track[vtrack[j]->trackID].quality << endl;
				if(sevt->track[t->trackID].quality<sevt->track[vtrack[j]->trackID].quality) {fBad=1;break;}
				else if(sevt->track[t->trackID].quality==sevt->track[vtrack[j]->trackID].quality){
					y4a[0] = propagateAfter(Geom->Dch[3].PosChamber.z, t).Y();
					y4b[0] = propagateBefore(Geom->Dch[3].PosChamber.z, t).Y();
					y4a[1] = propagateAfter(Geom->Dch[3].PosChamber.z, vtrack[j]).Y();
					y4b[1] = propagateBefore(Geom->Dch[3].PosChamber.z, vtrack[j]).Y();

					if(DEBUG_1) cout << "\t\tt1 y4Diff :\t" << fabs(y4a[0]-y4b[0]) << endl;
					if(DEBUG_1) cout << "\t\tt2 y4Diff :\t" << fabs(y4a[1]-y4b[1]) << endl;
					if(fabs(y4a[0]-y4b[0]) > fabs(y4a[1]-y4b[1])){
						fBad=1;
						break;
					}
				}
			}
		}
		if(fBad) {bad++; continue;} // pair of tracks with R<0.5cm at DCH1 plane, this track quality was lower

		//TODO Beta corrected Kaon???
		p1[0] = abcog_params.pkxoffp;
		p1[1] = abcog_params.pkyoffp;
		p1[2] = 0;
		v1[0] = abcog_params.pkdxdzp;
		v1[1] = abcog_params.pkdydzp;
		v1[2] = 1;

		p2[0] = t->detBPos.X();
		p2[1] = t->detBPos.Y();
		p2[2] = Geom->DCH.bz;
		v2[0] = sevt->track[t->trackID].bdxdz;
		v2[1] = sevt->track[t->trackID].bdydz;
		v2[2] = 1;
		closestApproach(p1, p2, v1, v2, &dmin, vertex);

		t->V0 = TVector3(vertex);
		t->unCda = dmin;

		if(DEBUG_1) cout << "\tCDA (no BF correction):\t\t" << t->unCda << "\t > 10" << endl;
		if(t->unCda>10) {bad++; continue;} //uncorrected CDA>10cm

		if(DEBUG_1) cout << "\tVertex Z (no BF correction):\t" << t->V0.Z() << "\t < -2000 || > 9000" << endl;
		if(t->V0.Z()<-2000 || t->V0.Z()>9000) {bad++; continue;} //Z(vertex)<-20m || Z(vertex)>90m (uncorrected)

		if(DEBUG_1) cout << "\ttime :\t\t" << fabs(sevt->track[t->trackID].time - sbur->tOffst.Dch) << "\t > 62.5" << endl;
		if(fabs(sevt->track[t->trackID].time - sbur->tOffst.Dch)>62.5){ bad++; continue;} //|t|>62.5ns

		if(DEBUG_1){
			cout << "First CDA " << t->unCda << endl;
			cout << "First Z " << t->V0.Z() << endl;
		}

		blue_1trk(vertex, i+1, &corrdxdz, &corrdydz, sevt);

		t->momentum = TVector3(corrdxdz, corrdydz, 1.).Unit();
		t->corrdxdz = corrdxdz;
		t->corrdydz = corrdydz;

		if(DEBUG_1){
			cout << "Detected position before: " << printVector3(t->detBPos) << endl;
			cout << "bdxdz: " << sevt->track[t->trackID].bdxdz << endl;
			cout << "bdydz: " << sevt->track[t->trackID].bdydz << endl;
			cout << "Detected position after: " << printVector3(t->detPos) << endl;
			cout << "dxdz: " << sevt->track[t->trackID].dxdz << endl;
			cout << "dydz: " << sevt->track[t->trackID].dydz << endl;

			cout << "corrdxdz: " << corrdxdz << endl;
			cout << "corrdydz: " << corrdydz << endl;
		}

		t->middlePos = propagateBefore((vertex[2]+Geom->Dch[0].PosChamber.z)/2., t);

		if(DEBUG_1){
			cout << "Corrected middle position: " << printVector3(t->middlePos) << endl;
		}

		p2[0] = t->middlePos.X();
		p2[1] = t->middlePos.Y();
		p2[2] = t->middlePos.Z();
		v2[0] = corrdxdz;
		v2[1] = corrdydz;
		v2[2] = 1;
		closestApproach(p1, p2, v1, v2, &dmin, vertex);

		t->cda = dmin;
		t->vertex = TVector3(vertex);
		if(DEBUG_1) cout << "Second CDA " << t->cda << endl;

		goodTrackID = i;
		//goodVertex = dmin;
		//vz = vertex[2];
		if(DEBUG_1) cout << "\tTrack is good" << endl;
		good++;
	}

	//trackID = goodTracks[0]ID;
	if(goodTrackID>=0) goodTracks.push_back(vtrack[goodTrackID]);
	//vertexCDA = goodVertex;
	return good;
}

/*int makePlots(superBurst *sbur,superCmpEvent *sevt, float tp){
	TH1D *mmass = (TH1D*)gDirectory->Get("mmass");
	TH1D *ncluster = (TH1D*)gDirectory->Get("NCluster");
	TH1D *ntrack = (TH1D*)gDirectory->Get("NTrack");

	float trackE = sqrt(tp*tp + pow(0.13957018,2));
	float beamEnergy = sqrt(65*65 + pow(0.493667,2));
	//float clusterEnergy =
	//float mmiss2 = pow(beamEnergy-trackE + ,2) - pow(65 - tp,2)

	ncluster->Fill(sevt->Ncluster);
	ntrack->Fill(sevt->Ntrack);
	return 0;
}*/

int ke2_getMomentumBin(float tp){
	double min = 13;
	int bin = -1;

	if(tp < min+7) bin=0;
	else bin = int((tp-20)/5)+1;

	if(bin==10) bin=9;

	if(DEBUG_7) cout << "Momentum :\t" << tp << endl;
	if(DEBUG_7) cout << "Bin :\t\t" << bin << endl;
	return bin;
}

double ke2_getZMinBin(int bin){
	double binsVals[10] = {-1000,-1000,0,500,500,500,1000,1000,1500,2500};

	if(DEBUG_7) cout << "ZMin :\t\t" << binsVals[bin] << endl;
	return binsVals[bin];
}

bool ke2_isCloseToHotCell(){
	int cpd,cell;
	TVector3 propPos;

	TH2D *hhotCellDistT = (TH2D*)gDirectory->Get("hotCellDistT");
	TH2D *hhotCellDistC = (TH2D*)gDirectory->Get("hotCellDistC");

	propPos = propagateAfter(Geom->Lkr.z, goodTracks[0]);
	GetCpdCellIndex(propPos.X(),propPos.Y(), &cpd, &cell);
	hhotCellDistT->Fill(propPos.X(), propPos.Y());

	if(DEBUG_10) cout << "Position : \t" << printVector3(propPos) << endl;
	if(DEBUG_10) cout << "Track CPD :\t" << cpd << endl;
	if(DEBUG_10) cout << "Track Cell :\t" << cell << endl;

	if((cpd==134) && (cell==56 || cell==57)){
		if(DEBUG_10) cout << "Track close to hot cell"<< endl;
		return true;
	}

	//int i=0;
	//int clID;
	for(int i=0; i<closeClusters.size(); i++){
		if(DEBUG_10) cout << "Working with cluster :\t" << closeClusters[i]->clusterID << endl;
		if(DEBUG_10) cout << "Position : \t" << printVector3(closeClusters[i]->position) << endl;
		GetCpdCellIndex(closeClusters[i]->unPosition.X(), closeClusters[i]->unPosition.Y(), &cpd,&cell);
		if(DEBUG_10) cout << "\tCluster CPD :\t" << cpd << endl;
		if(DEBUG_10) cout << "\tCluster Cell :\t" << cell << endl;
		hhotCellDistC->Fill(closeClusters[i]->position.X(), closeClusters[i]->position.Y());
		if((cpd==134) && (cell==56 || cell==57)){
			if(DEBUG_10) cout << "Cluster close to hot cell" << endl;
			return true;
		}
		i++;
	}

	if(DEBUG_10) cout << "None close to hot cell"<< endl;
	return false;
}

bool ke2_isTrackInExcludedZone(){
	TVector3 propPos;
	TH2D *hexcludedZone = (TH2D*)gDirectory->Get("excludedZone");

	propPos = propagateBefore(-180, goodTracks[0]);
	hexcludedZone->Fill(propPos.X(), propPos.Y());

	if(DEBUG_12) cout << "Propagated position : \t" << printVector3(propPos) << endl;

	if((propPos.X()>-40.0 && propPos.X()<-20.0 && propPos.Y()>-32.5 && propPos.Y()< 32.5) ||
			(propPos.X()>-50.0 && propPos.X()<-40.0 && propPos.Y()>-25.0 && propPos.Y()< 25.0) ||
			(propPos.X()>-55.0 && propPos.X()<-50.0 && propPos.Y()>-20.0 && propPos.Y()< 15.0) ||
			(propPos.X()>-20.0 && propPos.X()< 00.0 && propPos.Y()> 27.5 && propPos.Y()< 45.0) ||
			(propPos.X()> 00.0 && propPos.X()< 12.5 && propPos.Y()> 30.0 && propPos.Y()< 47.5) ||
			(propPos.X()> 12.5 && propPos.X()< 20.0 && propPos.Y()> 30.0 && propPos.Y()< 45.0) ||
			(propPos.X()>-15.0 && propPos.X()< 00.0 && propPos.Y()>-50.0 && propPos.Y()<-30.0) ||
			(propPos.X()> 00.0 && propPos.X()< 10.0 && propPos.Y()>-45.0 && propPos.Y()<-30.0) ||
			(propPos.X()> 10.0 && propPos.X()< 15.0 && propPos.Y()>-40.0 && propPos.Y()<-30.0)){
		if(DEBUG_12) cout << "Is in the excluded zone" << endl;
		return true;
	}
	if(DEBUG_12) cout << "Is out of the excluded zone" << endl;
	return false;
}

bool ke2_isTrackInDCHAcceptance(){
	float RDCH1, RDCH4;
	TVector3 propPos;
	TH2D *hDCHAcceptance1 = (TH2D*)gDirectory->Get("DCHAcceptance1");
	TH2D *hDCHAcceptance4 = (TH2D*)gDirectory->Get("DCHAcceptance4");
	TH1D *hRDCH1 = (TH1D*)gDirectory->Get("RDCH1");
	TH1D *hRDCH4 = (TH1D*)gDirectory->Get("RDCH4");

	if(DEBUG_13) cout << "\t\tPropagate" << endl;
	propPos = propagateBefore(Geom->Dch[0].PosChamber.z, goodTracks[0]);
	RDCH1 = sqrt(propPos.X()*propPos.X() + propPos.Y()*propPos.Y());
	if(DEBUG_13) cout << "Track at DCH1 :" << printVector3(propPos) << endl;
	if(DEBUG_13) cout << "RDCH1 :" << RDCH1 << endl;
	hRDCH1->Fill(RDCH1);
	hDCHAcceptance1->Fill(propPos.X(), propPos.Y());
	propPos = propagateAfter(Geom->Dch[3].PosChamber.z, goodTracks[0]);
	RDCH4 = sqrt(propPos.X()*propPos.X() + propPos.Y()*propPos.Y());
	if(DEBUG_13) cout << "Track at DCH4 :" << printVector3(propPos) << endl;
	if(DEBUG_13) cout << "RDCH4 :" << RDCH4 << endl;
	hRDCH4->Fill(RDCH4);
	hDCHAcceptance4->Fill(propPos.X(), propPos.Y());

	if((RDCH1<12 || RDCH1>115) || (RDCH4<14 || RDCH4>115)) return false;
	return true;
}

int ke2_failCut(superBurst *sbur, superCmpEvent *sevt, int i){
	TH1D *fails = (TH1D*)gDirectory->Get("Cuts");

	fails->Fill(i, 1);

	if(DEBUG_COM) cout << "Event is not passing selection" << endl;
	fprintf(fprt, "%i %i %i %i\n", sbur->nrun, sbur->time, sevt->timeStamp, i);
	return 0;
}

int ke2_failKe2Cut(superBurst *sbur, superCmpEvent *sevt, int i){
	TH1D *fails = (TH1D*)gDirectory->Get("CutsKe2");

	fails->Fill(i, 1);

	if(DEBUG_KE2) cout << "Event is not passing selection" << endl;
	fprintf(fprt, "%i %i %i ke2:%i\n", sbur->nrun, sbur->time, sevt->timeStamp, i);
	return 1;
}

double ke2_getMMinBin(int bin){
	double binsVals[10] = {16, 14, 13, 13, 13, 13, 13, 14, 15, 16};
	return binsVals[bin]*1.E-3;
}
double ke2_getMMaxBin(int bin){
	double binsVals[10] = {13, 11, 10, 10, 10, 10, 10, 11, 12, 12};
	return binsVals[bin]*1.E-3;
}

int ke2_ke2Selection(superBurst *sbur, superCmpEvent *sevt, double mmass2){
	int pbin;
	double minBin, maxBin;

	TVector3 propPos;
	double radius;
	double radius0;
	double trckClTime;
	double eopMin;
	double eop;

	//Histograms
	TH1D *hmmass2= (TH1D*)gDirectory->Get("mmass2");
	TH1D *hclusterStatus= (TH1D*)gDirectory->Get("clusterStatus");
	TH1D *htrckClDist= (TH1D*)gDirectory->Get("trckClDist");
	TH1D *htrckClTime= (TH1D*)gDirectory->Get("trckClTime");
	TH1D *heop= (TH1D*)gDirectory->Get("eop");

	//1) Ke2 L2 trigger word bit ON [for info: TW bit #10, LKR_MB bit = PU chan 6 bit 1, 1TRKL bit = PU chan 4 bit 9, PU since run 20154] [Data only];
	if(DEBUG_KE2_1) cout << "~~~~ Cut Ke2 1 ~~~~" << endl;
	bitset<32> triggBit(sevt->trigWord);
	bitset<32> triggMask(0x400);
	bitset<32> comb(sevt->trigWord & 0x400);
	if(DEBUG_KE2_1) cout << "Word :\t " << triggBit << endl;
	if(DEBUG_KE2_1) cout << "Mask :\t " << triggMask << endl;
	if(DEBUG_KE2_1) cout << "comb :\t " << comb << endl;
	if(DEBUG_KE2_1) cout << "L2 Trigger word :\t " << hex << (sevt->trigWord & 0x400) << dec << "\t == 0: rejected" << endl;
	if((sevt->trigWord & 0x400) ==0) return ke2_failKe2Cut(sbur, sevt, 1);

	//What if no cluster?
	if(assocClusters.size()>0){
		//3) Associated LKr cluster quality: status<4.   [reject clusters containing cells with ADC saturation]
		if(DEBUG_KE2_3) cout << "~~~~ Cut Ke2 3 ~~~~" << endl;
		if(DEBUG_KE2_3) cout << "Cluster status :\t " << sevt->cluster[assocClusters[0]->clusterID].status << "\t >= 4: rejected" << endl;
		hclusterStatus->Fill(sevt->cluster[assocClusters[0]->clusterID].status);
		if(sevt->cluster[assocClusters[0]->clusterID].status>=4) return ke2_failKe2Cut(sbur, sevt, 3);


		//4) Track-cluster distance: R<R0 (in the track-cluster matching Z plane defined above).
		//	R0=5cm (for p<25 GeV/c), R0=1.5cm (for p≥25 GeV/c).
		if(DEBUG_KE2_4) cout << "~~~~ Cut Ke2 4 ~~~~" << endl;
		propPos = propagateAfter(assocClusters[0]->position.Z(),goodTracks[0]);
		radius = distance2D(propPos, assocClusters[0]->position);

		if(goodTracks[0]->pMag<25){
			radius0 = 5;
			eopMin = 0.90;
		}
		else{
			radius0 = 1.5;
			eopMin = 0.95;
		}
		if(DEBUG_KE2_4) cout << "R(t,cl) :\t " << radius << "\t > " << radius0 << ": rejected" << endl;
		htrckClDist->Fill(radius);
		if(radius>radius0) return ke2_failKe2Cut(sbur, sevt, 4);


		//5) Track-cluster timing: |Δt|<12ns [Data only].
		if(DEBUG_KE2_5) cout << "~~~~ Cut Ke2 5 ~~~~" << endl;
		trckClTime = fabs(sevt->track[goodTracks[0]->trackID].time-sevt->cluster[assocClusters[0]->clusterID].time);
		if(DEBUG_KE2_5) cout << "DeltaT :\t " << trckClTime << "\t > 12: rejected" << endl;
		htrckClTime->Fill(trckClTime);
		if(trckClTime>12) return ke2_failKe2Cut(sbur, sevt, 5);
	}

	if(DEBUG_KE2_6) cout << "~~~~ Cut Ke2 6 ~~~~" << endl;
	eop = goodTracks[0]->clusterEnergy / goodTracks[0]->pMag;
	if(DEBUG_KE2_6) cout << "E/p :\t " << eop << "\t < " << eopMin << " || > 1.10: rejected" << endl;
	//6) (E/p)min < E/p < 1.10. The lower limit: (E/p)min = 0.90 (for p<25 GeV/c), (E/p)min = 0.95 (for p≥25 GeV/c).
	heop->Fill(eop);
	if(eop<eopMin || eop>1.10) return ke2_failKe2Cut(sbur, sevt, 6);

	//2) −Mmin ×10−3 < Mmiss2(e) < Mmax ×10−3 (GeV/c2)2.
	//	Mmin in track momentum bins: {16, 14, 13, 13, 13, 13, 13, 14, 15, 16}.
	//	Mmax in track momentum bins (no Pb wall): {13, 11, 10, 10, 10, 10, 10, 11, 12, 12}.
	//	Mmax in track momentum bins (Pb wall): {10, 10, 10, 10, 10, 10, 10, 11, 11, 11}.
	if(DEBUG_KE2_2) cout << "~~~~ Cut Ke2 2 ~~~~" << endl;
	pbin = ke2_getMomentumBin(goodTracks[0]->pMag);
	minBin = -ke2_getMMinBin(pbin);
	maxBin = ke2_getMMaxBin(pbin);
	if(DEBUG_KE2_2) cout << "mmass2 :\t " << mmass2 << "\t < " << minBin << " || > " << maxBin << ": rejected" << endl;
	hmmass2->Fill(mmass2);
	if(mmass2<minBin || mmass2>maxBin) return ke2_failKe2Cut(sbur, sevt, 2);


	return 0;
}

void ke2_passSelection(superBurst *sbur, superCmpEvent *sevt){
	if(DEBUG_) cout << "Event is passing selection" << endl;
	fprintf(fprt2, "%i %i %i\n", sbur->nrun, sbur->time, sevt->timeStamp);
}

int nico_ke2Select(superBurst *sbur,superCmpEvent *sevt) {
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	int goodTracksNb;
	int badClusters;
	float minZ;

	//Acceptance
	TVector3 lkrTrackPos;
	int lkrAcceptance;
	bool excludedZone;
	bool dchAcceptance;

	double MK = abcog_params.mkp;
	double Me = 0.00051099891;

	//Histograms
	TH1D *htimeDchOffset = (TH1D*)gDirectory->Get("timeDchOffset");
	TH1D *hdeadCellDist = (TH1D*)gDirectory->Get("deadCellDist");
	TH1D *hcda = (TH1D*)gDirectory->Get("cda");
	TH2D *hlkrAcceptance = (TH2D*)gDirectory->Get("lkrAcceptance");
	TH1D *hvertexZ = (TH1D*)gDirectory->Get("vertexZ");
	TH1D *htrackMomentum = (TH1D*)gDirectory->Get("trackMomentum");
	TH1D *htrackQuality = (TH1D*)gDirectory->Get("trackQuality");

	if(DEBUG_) cout << endl;

	// 2) Common selection
	if(DEBUG_1) cout << "~~~~ Cut 1 ~~~~" << endl;
	goodTracksNb = ke2_getGoodTracks(sbur, sevt);
	if(DEBUG_1) cout << "Good tracks :\t\t" << goodTracksNb << "\t != 1 : rejected" << endl;
	if(goodTracksNb!=1) return ke2_failCut(sbur, sevt, 1); //!(Exactly one good track) DCH Veto
	if(DEBUG_1) cout << "Good track ID :\t\t" << goodTracks[0]->trackID << endl;
	goodTracks[0]->energy = sqrt(Me*Me + goodTracks[0]->pMag*goodTracks[0]->pMag);

	if(DEBUG_2) cout << "~~~~ Cut 2 ~~~~" << endl;
	if(DEBUG_2) cout << "Track charge :\t" << sevt->track[goodTracks[0]->trackID].q << "\t != " << beamCharge << " : rejected" << endl;
	if(sevt->track[goodTracks[0]->trackID].q!=beamCharge) return ke2_failCut(sbur, sevt, 2); //!(track charge = beam charge)

	if(DEBUG_3) cout << "~~~~ Cut 3 ~~~~" << endl;
	if(DEBUG_3) cout << "Track quality :\t" << sevt->track[goodTracks[0]->trackID].quality << "\t < 0.7 : rejected" << endl;
	htrackQuality->Fill(sevt->track[goodTracks[0]->trackID].quality);
	if(sevt->track[goodTracks[0]->trackID].quality<0.7) return ke2_failCut(sbur, sevt, 3); //!(track quality >0.7)

	if(DEBUG_4) cout << "~~~~ Cut 4 ~~~~" << endl;
	if(DEBUG_4) cout << "Track momentum :\t" << goodTracks[0]->pMag << "\t <13 || > 65 : rejected" << endl;
	htrackMomentum->Fill(goodTracks[0]->pMag);
	if(goodTracks[0]->pMag<13 || goodTracks[0]->pMag>65) return ke2_failCut(sbur, sevt, 4); //!(13Gev<p<65GeV)

	if(DEBUG_5) cout << "~~~~ Cut 5 ~~~~" << endl;
	reassociateCluster(sevt);
	badClusters = ke2_vetoLKr(sevt);
	if(DEBUG_5) cout << "Bad clusters :\t" << badClusters << "\t > 0 : rejected" << endl;
	if(badClusters>0) return ke2_failCut(sbur, sevt, 5); //!no bad cluster LKr Veto

	outTree->Fill();

	// z cuts---
	// check zmax
	if(DEBUG_6) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(DEBUG_6) cout << "Vertex Z :\t " << goodTracks[0]->vertex.Z() << "\t > 9000 : rejected" << endl;
	hvertexZ->Fill(goodTracks[0]->vertex.Z());
	if(goodTracks[0]->vertex.Z() > 9000) return ke2_failCut(sbur, sevt, 6); // !z<z_max

	if(DEBUG_7) cout << "~~~~ Cut 7 ~~~~" << endl;
	minZ = ke2_getZMinBin(ke2_getMomentumBin(goodTracks[0]->pMag));
	if(DEBUG_7) cout << "Vertex Z :\t" << goodTracks[0]->vertex.Z() << "\t < " << minZ << " : rejected" << endl;
	if(goodTracks[0]->vertex.Z() < minZ) return ke2_failCut(sbur, sevt, 7); // !z>z_min(p)

	if(DEBUG_8) cout << "~~~~ Cut 8 ~~~~" << endl;
	if(DEBUG_8) cout << "CDA :\t" << goodTracks[0]->cda << "\t > 3.5 : rejected" << endl;
	hcda->Fill(goodTracks[0]->cda);
	if(goodTracks[0]->cda>3.5) return ke2_failCut(sbur, sevt, 8); // !vertex CDA < 3.5

	if(DEBUG_9) cout << "~~~~ Cut 9 ~~~~" << endl;
	lkrTrackPos = propagateAfter(Geom->Lkr.z, goodTracks[0]);
	hlkrAcceptance->Fill(lkrTrackPos.X(), lkrTrackPos.Y());
	lkrAcceptance = LKr_acc(sbur->nrun, lkrTrackPos.X(), lkrTrackPos.Y(), 8);
	if(DEBUG_9) cout << "Mauro condition :\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=0) return ke2_failCut(sbur, sevt, 9); // !satisfy Mauro's condition

	if(DEBUG_10) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(ke2_isCloseToHotCell()) return ke2_failCut(sbur, sevt, 10); // track or cluster <10cm of the track is close to a hot LKr cell

	if(DEBUG_12) cout << "~~~~ Cut 12 ~~~~" << endl;
	excludedZone = ke2_isTrackInExcludedZone();
	if(DEBUG_12) cout << "Beam Halo :\t" << excludedZone << "\t == " << true << ": rejected" << endl;
	//Only for K-
	//if(goodTracks[0]->t.q==-1 && excludedZone) return failCut(sbur, sevt, 12); //beam halo suppression

	if(DEBUG_13) cout << "~~~~ Cut 13 ~~~~" << endl;
	dchAcceptance = ke2_isTrackInDCHAcceptance();
	if(DEBUG_13) cout << "DCH acceptance :\t" << dchAcceptance << "\t != " << true << ": rejected" << endl;
	if(!dchAcceptance) return ke2_failCut(sbur, sevt, 13); //trak not in DCH acceptance

	if(DEBUG_15) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(DEBUG_15) cout << "Dead Cell :\t" << sevt->track[goodTracks[0]->trackID].dDeadCell << "\t < 2 : rejected" << endl;
	hdeadCellDist->Fill(sevt->track[goodTracks[0]->trackID].dDeadCell);
	if(sevt->track[goodTracks[0]->trackID].dDeadCell<=2) return ke2_failCut(sbur, sevt, 15); //!track distance to deadCell > 2cm

	if(DEBUG_16) cout << "~~~~ Cut 16 ~~~~" << endl;
	if(DEBUG_16) cout << "Track time :\t" << fabs(sevt->track[goodTracks[0]->trackID].time-sbur->tOffst.Dch) << "\t > 20 : rejected" << endl;
	htimeDchOffset->Fill(fabs(sevt->track[goodTracks[0]->trackID].time-sbur->tOffst.Dch));
	if(fabs(sevt->track[goodTracks[0]->trackID].time-sbur->tOffst.Dch) > 20) return ke2_failCut(sbur, sevt, 16); //!track time wrt DCH offset

	// 3) Ke2 specific seletion
	if(DEBUG_KE2){
		cout << "Track energy " << goodTracks[0]->energy << endl;
		cout << "Track momentum " << goodTracks[0]->pMag << endl;
		cout << "Track momentum " << printVector3(goodTracks[0]->momentum) << endl;
		cout << "Track uncorrected momentum " << goodTracks[0]->unPMag << endl;
		cout << "Track uncorrected momentum " << printVector3(goodTracks[0]->unMomentum) << endl;
		cout << "Kaon momentum magnitude  " << kaonP << endl;
		cout << "Kaon momentum " << printVector3(kaonMomentum) << endl;

	}

	MK = 0.493677;
	if(ke2_ke2Selection(sbur, sevt, missMass2(MK, Me, kaonMomentum*kaonP, goodTracks[0]->momentum*goodTracks[0]->pMag))==1) return 0;
	//if(ke2Selection(sbur, sevt, missMass2Man(MK, Me, kaonP, abcog_params.pkdxdzp, abcog_params.pkdydzp, goodTracks[0]->energy, goodTracks[0]->corrdxdz, goodTracks[0]->corrdydz))==1) return 0;
	//makePlots(sbur, sevt, t);
	ke2_passSelection(sbur,sevt);
	/*----------- End of user C code -----------*/
	return 0;
}
