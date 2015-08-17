/*
 * ScanCuts.h
 *
 *  Created on: Feb 17, 2015
 *      Author: ncl
 */

#ifndef SCANCUTS_H_
#define SCANCUTS_H_

#include <TObject.h>
class Cuts: public TObject{
public:
	Cuts();

	void print();
public:
	int triggerMask;
	int numVertex3;
	int minZVertex;
	int maxZVertex;
	int maxChi2Vertex;
	int maxExtraTracks;
	int maxTrackTime;
	bool boolBadTrack;
	int numBadTrackCombi;
	int numXCandidates;
	bool boolBadECandidates;
	int minTrackMomentum;
	int maxTrackMomentum;
	int numAddGoodCluster;
	int lkrAcceptance;
	int minGammaEnergy;
	int minDeadCellDist;
	int minGammaDCHRadius;
	struct k2pi_t{
		int minTotalMomentum;
		int maxTotalMomentum;
		double maxPt;
		double maxPi0MassDiff;
		double minKaonMassDiff;
		double maxKaonMassDiff;
	} k2pi;
	struct kmu3_t{
		int maxTotalMomentum;
		double minPt;
		double maxPt;
		double maxPi0MassDiff;
		double maxMissMassSq;
	} kmu3;

	ClassDefNV(Cuts, 1)
};

class ScanCuts : public Cuts{
public:
	ScanCuts();
	~ScanCuts();

	bool addParseCutsFile(std::string fileName);
	void generateLists(int number);
	void print();
	void print(int index);

	int loadNextList();
	int loadList(int i);
	void loadDefault();

	int getNLists() {
		return cutsLists.size();
	}
	int getCurrentList() const {
		return currentList;
	}

	int getDefaultIndex() const {
		return defaultIndex;
	}

	void setDefaultIndex(int defaultIndex) {
		this->defaultIndex = defaultIndex;
	}

private:
	bool parseCuts(std::string fileName);
	int defaultIndex;
	int currentList;				//! current directory
	std::vector<Cuts> cutsLists; 	//! current directory


	ClassDefNV(ScanCuts, 1)
};

#endif /* SCANCUTS_H_ */
