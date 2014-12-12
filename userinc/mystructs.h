/*
 * mystructs.h
 *
 *  Created on: 4 Feb 2014
 *      Author: ncl
 */
#define __CINT_NICO__ 1
#ifndef __CINT_NICO__
#include "reader.h"
#else
class trak;
class cluster;
#endif

#include "TVector3.h"
#include <vector>
#include <iostream>
using namespace std;
#ifndef MYSTRUCTS_H_
#define MYSTRUCTS_H_


class CorrectedTrack: public TObject {
public:
	CorrectedTrack(){
		energy = 0;
		pMag = 0;
		cda = 0;
		trackID = -1;
		clusterEnergy = 0;
		unPMag = 0;
		unCda = 0;
		corrdxdz = 0.;
		corrdydz = 0.;
	};

	/*CorrectedTrack(const CorrectedTrack &obj): TObject(obj){
		energy = obj.energy;
		momentum = obj.momentum;
		pMag = obj.pMag;
		vertex = obj.vertex;
		cda = obj.cda;
		middlePos = obj.middlePos;
		corrdxdz = obj.corrdxdz;
		corrdydz = obj.corrdydz;
		detBPos = obj.detBPos;
		detPos = obj.detPos;
		//t = obj.t;
		trackID = obj.trackID;
		clusterEnergy = obj.clusterEnergy;
		unMomentum = obj.unMomentum;
		unBMomentum = obj.unBMomentum;
		unPMag = obj.unPMag;
		V0 = obj.V0;
		unCda = obj.unCda;
	}*/
	//All momenta are unit vectors
	//Corrected values
	double energy;
	TVector3 momentum;
	double pMag;
	TVector3 vertex;
	double cda;
	TVector3 middlePos;

	double corrdxdz;
	double corrdydz;

	//Original values (no correction needed)
	TVector3 detBPos;
	TVector3 detPos;
	//trak t; //!
	int trackID;

	//Original values (correction needed)
	double clusterEnergy;
	TVector3 unMomentum;
	TVector3 unBMomentum;
	double unPMag;

	//Vertex
	TVector3 V0;
	double unCda;

	ClassDefNV(CorrectedTrack, 1);
};


class CorrectedCluster: public TObject{
public:
	CorrectedCluster(){
		clusterID =-1;
		energy = 0;
	};
	TVector3 position;

	TVector3 unPosition;

	double energy;

	int clusterID;
	//cluster c; //!

	ClassDefNV(CorrectedCluster, 1);
};


typedef struct eventID_t{
	eventID_t(int r, int b, int t){rnum=r;bnum=b;timestamp=t;};
	int rnum;
	int bnum;
	int timestamp;
} eventID;

class pi0dEvent: public TObject{
public:
	pi0dEvent():
	gamma(NULL),
	ep(NULL),
	em(NULL),
	pip(NULL),
	mee(0),
	mPi0(0),
	mK(0){};
	CorrectedCluster *gamma;
	CorrectedTrack *ep, *em, *pip;
	TVector3 pGamma;
	TVector3 pTotal;
	double mee;
	double mPi0;
	double mK;
	double weight;
	int runNumber;
	int burstNumber;
	int timestamp;
	double zVtx;
	TVector3 posLKr;
	TVector3 posDCH1;
	double x;
	double xTrue;
	double meeTrue;
	TVector3 epPTrue, emPTrue;
	double epETrue, emETrue;

	ClassDefNV(pi0dEvent, 1)
};

typedef struct cutsValues_t{
	int triggerMask;
	int numVertex3;
	int minZVertex;
	int maxZVertex;
	int maxChi2Vertex;
	int maxExtraTracks;
	int maxTrackTime;
	bool boolBadTrack;
	int numBadTrackCombi;
	int numPiCandidates;
	bool boolBadECandidates;
	int minTrackMomentum;
	int maxTrackMomentum;
	int numAddGoodCluster;
	int lkrAcceptance;
	int minGammaEnergy;
	int minDeadCellDist;
	int minGammaDCHRadius;
	int minTotalMomentum;
	int maxTotalMomentum;
	double maxPt;
	double maxPi0MassDiff;
	double minKaonMassDiff;
	double maxKaonMassDiff;
} cutsValues;

#endif /* MYSTRUCTS_H_ */
