/*
 * Drawer.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#include "Drawer.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TNamed.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../Samples/CombineDataSample.h"
#include "../Samples/CombineMCSample.h"
#include "../Samples/CombineSample.h"
#include "../Samples/FitDataSample.h"
#include "../Samples/FitMCSample.h"
#include "StackDrawer.h"
#include "StackRatioDrawer.h"

using namespace std;

Drawer::Drawer() {
	// TODO Auto-generated constructor stub

}

Drawer::~Drawer() {
	// TODO Auto-generated destructor stub
}

void Drawer::drawFitResult(MinuitFitter* fit, std::vector<Sample*> mcSamples,
		std::vector<Sample*> dataSamples, Sample* finalMCSample,
		Sample* finalDataSample) {
	StackRatioDrawer d;

	d.setName(fit->getName());
	d.setTitle(fit->getName());

	for(Sample* sample : mcSamples){
		FitMCSample *ssample = static_cast<FitMCSample*>(sample->getSubSample(0));
		TH1D *dAlpha_c = (TH1D*) ssample->getAlpha()->Clone(TString::Format("dAlpha_c%s%i", fit->getName().c_str(), sample->getIndex()));
		dAlpha_c->Scale(fit->getNorm());
		TH1D *dBeta_c = (TH1D*) ssample->getBeta()->Clone(TString::Format("dBeta_c%s%i", fit->getName().c_str(), sample->getIndex()));
		dBeta_c->Scale(fit->getNorm() * 2. * fit->getFormFactor());
		TH1D *dGamma_c = (TH1D*) ssample->getGamma()->Clone(TString::Format("dGamma_c%s%i", fit->getName().c_str(), sample->getIndex()));
		dGamma_c->Scale(fit->getNorm() * pow(fit->getFormFactor(), 2));
		d.AddHisto1(dAlpha_c, TString::Format("%s #alpha", sample->getLegend().c_str()).Data());
		d.AddHisto1(dBeta_c, TString::Format("%s #beta", sample->getLegend().c_str()).Data());
		d.AddHisto1(dGamma_c, TString::Format("%s #gamma", sample->getLegend().c_str()).Data());
	}

	for(Sample* sample : dataSamples){
			FitDataSample *ssample = static_cast<FitDataSample*>(sample->getSubSample(0));
			TH1D *dSig_c = (TH1D*) ssample->getSig()->Clone(TString::Format("dSig_c%s%i", fit->getName().c_str(), sample->getIndex()));
			cout << dSig_c->GetBinContent(1) << endl;
			d.AddHisto2(dSig_c, sample->getLegend().c_str());
		}

	TH1D* ratio = StackRatioDrawer::buildRatio(finalMCSample->getSubSample(0)->getMainHisto(), finalDataSample->getSubSample(0)->getMainHisto());
	d.SetSecondary(ratio);

	d.draw();
}

void Drawer::drawFitPreparation(std::vector<Sample*> mcSamples,
		std::vector<Sample*> dataSamples, std::string title) {

	StackDrawer d1,d2,d3,dNew,dAlpha,dBeta,dGamma, dSig;

	d1.setName(title);
	d1.setTitle(title);
	d2.setName(title);
	d2.setTitle(title);
	d3.setName(title);
	d3.setTitle(title);
	dNew.setName(title);
	dNew.setTitle(title);
	dAlpha.setName(title);
	dAlpha.setTitle(title);
	dBeta.setName(title);
	dBeta.setTitle(title);
	dGamma.setName(title);
	dGamma.setTitle(title);

	for(auto sample : mcSamples) {
		FitMCSample *ssample = static_cast<FitMCSample*>(sample->getSubSample(0));
		d1.AddHisto((TH1D*) ssample->getD1()->Clone(), sample->getLegend());
		d2.AddHisto((TH1D*) ssample->getD2()->Clone(), sample->getLegend());
		d3.AddHisto((TH1D*) ssample->getD3()->Clone(), sample->getLegend());
		dNew.AddHisto((TH1D*) ssample->getNew()->Clone(), sample->getLegend());
		dAlpha.AddHisto((TH1D*) ssample->getAlpha()->Clone(), sample->getLegend());
		dBeta.AddHisto((TH1D*) ssample->getBeta()->Clone(), sample->getLegend());
		dGamma.AddHisto((TH1D*) ssample->getGamma()->Clone(), sample->getLegend());
	}

	for(auto sample : dataSamples){
		FitDataSample *ssample = static_cast<FitDataSample*>(sample->getSubSample(0));
		dSig.AddHisto((TH1D*) ssample->getSig()->Clone(), sample->getLegend());
	}

	TCanvas *c1 = new TCanvas("cInput", "Inputs", 1600, 800);
	c1->Divide(2, 2);
	d1.generate(static_cast<TPad*>(c1->GetPad(1)));
	d2.generate(static_cast<TPad*>(c1->GetPad(2)));
	d3.generate(static_cast<TPad*>(c1->GetPad(3)));
	dSig.generate(static_cast<TPad*>(c1->GetPad(4)));

	TCanvas *c2 = new TCanvas("cInputsNew", "Inputs New", 1600, 800);
	dNew.generate(c2);

	TCanvas *c3 = new TCanvas("cABG", "Alpha, Beta, gamma", 1600, 800);
	c3->Divide(2, 2);
	dAlpha.generate(static_cast<TPad*>(c3->GetPad(1)));
	dBeta.generate(static_cast<TPad*>(c3->GetPad(1)));
	dGamma.generate(static_cast<TPad*>(c3->GetPad(1)));
}

void Drawer::drawCombineStack(std::vector<Sample*> mcSamples,
		std::vector<Sample*> dataSamples, Sample* finalMCSample, Sample* finalDataSample) {
	vector<StackRatioDrawer*> d;

	for(auto sample : mcSamples){
		CombineMCSample *ssample = static_cast<CombineMCSample*>(sample->getSubSample(0));
		vector<TH1D*>d1 = ssample->getD1();
		for (unsigned int i = 0; i < d1.size(); ++i) {
			if(d.size()<=i) d.push_back(new StackRatioDrawer());
			d[i]->AddHisto1(d1[i], sample->getLegend(), "f");
		}
	}

	for(auto sample : dataSamples){
		CombineDataSample *ssample = static_cast<CombineDataSample*>(sample->getSubSample(0));
		vector<TH1D*>d1 = ssample->getD1();
		for (unsigned int i = 0; i < d1.size(); ++i) {
			if(d.size()<=i) d.push_back(new StackRatioDrawer());
			d[i]->AddHisto1(d1[i], sample->getLegend(), "lep");
		}
	}
	vector<TH1D*> data= static_cast<CombineDataSample*>(finalDataSample->getSubSample(0))->getD1();
	vector<TH1D*> mc= static_cast<CombineMCSample*>(finalMCSample->getSubSample(0))->getD1();
	for(unsigned int iCanvas=0; iCanvas<data.size(); ++iCanvas){
		//TCanvas *c = new TCanvas(TString::Format("c%i", iCanvas), data[iCanvas]->GetTitle());
		TH1D* ratio = StackRatioDrawer::buildRatio(mc[iCanvas], data[iCanvas]);
		d[iCanvas]->SetSecondary(ratio);
		d[iCanvas]->draw();
	}
}




