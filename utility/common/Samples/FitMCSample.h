/*
 * FitSample.h
 *
 *  Created on: Dec 23, 2015
 *      Author: nlurkin
 */

#ifndef COMMON_FITMCSAMPLE_H_
#define COMMON_FITMCSAMPLE_H_

#include "../Interface/SubSample.h"
#include "../Interface/MCSample.h"

class FitMCSample: public SubSample, public MCSample {
public:
	typedef struct bContent_t{
		double d1,d2,d3,dNew,dAlpha,dBeta,dGamma;

		struct bContent_t &operator +=(struct bContent_t other){
			d1 += other.d1;
			d2 += other.d2;
			d3 += other.d3;
			dNew += other.dNew;
			dAlpha += other.dAlpha;
			dBeta += other.dBeta;
			dGamma += other.dGamma;

			return *this;
		}
	} bContent;

	FitMCSample();
	virtual ~FitMCSample();

	virtual void processEvent(ROOTPhysicsEvent *eventBrch, ROOTBurst *burstBrch,
			ROOTRawEvent *rawBrch, ROOTCorrectedEvent *corrBrch,
			ROOTFileHeader *headerBrch, ROOTMCEvent *mcEvent, NGeom *geomBrch,
			std::vector<bool> *cutsPass, const ConfigFile *cfg, const RunWeights *weights);
	virtual void doGet(TDirectory* inputFD, TFile* tempFD);
	virtual void doWrite();
	virtual void doSetName();
	virtual void initHisto(int nbins, double* bins, const ConfigFile *cfg);
	virtual void scaleToData(bContent totalMC, double nData);
	virtual void setPlotStyle(std::vector<int> color);
	virtual TH1D* getMainHisto() { return dAlpha; }
	virtual double getFFIntegral(double a);
	virtual void renameHisto() {};

	bContent getIntegrals();

	void scale();

	bContent getBinContent(int bin);
	bContent getBinError(int bin);
	virtual SubSample* Add(const SubSample* other);

	const TH1D* getD1() const {
		return d1;
	}

	const TH1D* getD2() const {
		return d2;
	}

	const TH1D* getD3() const {
		return d3;
	}

	const TH1D* getAlpha() const {
		return dAlpha;
	}

	const TH1D* getBeta() const {
		return dBeta;
	}

	const TH1D* getGamma() const {
		return dGamma;
	}

	const TH1D* getNew() const {
		return dNew;
	}

private:
	TH1D *d1, *d2, *d3;
	TH1D *dNew;
	TH1D *dAlpha, *dBeta, *dGamma;
};

#endif /* COMMON_FITMCSAMPLE_H_ */
