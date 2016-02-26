/*
 * Drawer.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#include "Drawer.h"

#include <Rtypes.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TAttText.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TString.h>
#include <cmath>
#include <string>
#include <vector>

#include "../Samples/CombineDataSample.h"
#include "../Samples/CombineMCSample.h"
#include "../Samples/CombineSample.h"
#include "../Samples/FitDataSample.h"
#include "../Samples/FitMCSample.h"
#include "ScanDrawer.h"
#include "StackDrawer.h"
#include "StackRatioDrawer.h"

using namespace std;

Drawer::Drawer() {
	// TODO Auto-generated constructor stub

}

Drawer::~Drawer() {
	// TODO Auto-generated destructor stub
}

void Drawer::drawFitResult(std::vector<MinuitFitter*> vFit,
		std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples,
		Sample* finalMCSample, Sample* finalDataSample, int index) {
	StackRatioDrawer d;

	MinuitFitter* fit = vFit[index];

	d.setName(fit->getName());
	d.setTitle(fit->getName());

	for (Sample* sample : mcSamples) {
		FitMCSample *ssample = static_cast<FitMCSample*>(sample->getSubSample(
				index));
		TH1D *dAlpha_c = (TH1D*) ssample->getAlpha()->Clone(
				TString::Format("dAlpha_c%s%i", fit->getName().c_str(),
						sample->getIndex()));
		dAlpha_c->Scale(fit->getNorm());
		TH1D *dBeta_c = (TH1D*) ssample->getBeta()->Clone(
				TString::Format("dBeta_c%s%i", fit->getName().c_str(),
						sample->getIndex()));
		dBeta_c->Scale(fit->getNorm() * 2. * fit->getFormFactor());
		TH1D *dGamma_c = (TH1D*) ssample->getGamma()->Clone(
				TString::Format("dGamma_c%s%i", fit->getName().c_str(),
						sample->getIndex()));
		dGamma_c->Scale(fit->getNorm() * pow(fit->getFormFactor(), 2));
		d.AddHisto1(dAlpha_c,
				TString::Format("%s #alpha", sample->getLegend().c_str()).Data());
		d.AddHisto1(dBeta_c,
				TString::Format("%s #beta", sample->getLegend().c_str()).Data());
		d.AddHisto1(dGamma_c,
				TString::Format("%s #gamma", sample->getLegend().c_str()).Data());
	}

	for (Sample* sample : dataSamples) {
		FitDataSample *ssample =
				static_cast<FitDataSample*>(sample->getSubSample(index));
		TH1D *dSig_c = (TH1D*) ssample->getSig()->Clone(
				TString::Format("dSig_c%s%i", fit->getName().c_str(),
						sample->getIndex()));
		d.AddHisto2(dSig_c, sample->getLegend().c_str());
	}

	TH1D* ratio = StackRatioDrawer::buildRatio(
			finalMCSample->getSubSample(index)->getMainHisto(),
			finalDataSample->getSubSample(index)->getMainHisto());
	d.SetSecondary(ratio);

	TCanvas* c1 = new TCanvas(fit->getName().c_str(), fit->getName().c_str());
	d.generate(c1);
	TPaveText *fitR = new TPaveText(0.47, 0.78, 0.77, 0.94, "NDC BR");
	fitR->AddText("Fit result");
	fitR->AddLine(0., 0.7, 1., 0.7);
	fitR->AddText(Form("G = %f #pm %f", fit->getNorm(), fit->getNormErr()));
	fitR->AddText(
			Form("FF = %f #pm %f", fit->getFormFactor(),
					fit->getFormFactorErr()));
	fitR->SetTextAlign(12);
	TF1 *f = new TF1("f", "[0]*(1+[1]*2.0*x+[1]*[1]*x*x)", 0., 1.);
	f->SetParameter(0, fit->getNorm());
	f->SetParameter(1, fit->getFormFactor());
	f->SetLineColor(kRed);
	c1->cd(1);
	c1->cd(1)->SetLogx(true);
	c1->cd(1)->SetLogy(true);
	fitR->Draw();
	c1->cd(2);
	f->Draw("LSAME");
	ratio->SetMarkerColor(kRed);
	c1->cd(2)->SetLogx(true);
}

void Drawer::drawFitPreparation(std::vector<Sample*> mcSamples,
		std::vector<Sample*> dataSamples, std::string title) {

	StackDrawer d1, d2, d3, dNew, dAlpha, dBeta, dGamma, dSig;

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

	for (auto sample : mcSamples) {
		FitMCSample *ssample =
				static_cast<FitMCSample*>(sample->getSubSample(0));
		d1.AddHisto((TH1D*) ssample->getD1()->Clone(), sample->getLegend());
		d2.AddHisto((TH1D*) ssample->getD2()->Clone(), sample->getLegend());
		d3.AddHisto((TH1D*) ssample->getD3()->Clone(), sample->getLegend());
		dNew.AddHisto((TH1D*) ssample->getNew()->Clone(), sample->getLegend());
		dAlpha.AddHisto((TH1D*) ssample->getAlpha()->Clone(),
				sample->getLegend());
		dBeta.AddHisto((TH1D*) ssample->getBeta()->Clone(),
				sample->getLegend());
		dGamma.AddHisto((TH1D*) ssample->getGamma()->Clone(),
				sample->getLegend());
	}

	for (auto sample : dataSamples) {
		FitDataSample *ssample =
				static_cast<FitDataSample*>(sample->getSubSample(0));
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

void Drawer::drawFitScan(std::vector<MinuitFitter*> fit,
		std::vector<Sample*> mcSamples, std::vector<Sample*> dataSamples,
		Sample* finalMCSample, Sample* finalDataSample, vector<int> use) {

	ScanDrawer d;
	d.setDefaultCutValue(std::find(use.begin(), use.end(), mcSamples[0]->getMainSubSample())-use.begin());
	const ScanCuts *defaultCuts = mcSamples[0]->getSubSample(
			mcSamples[0]->getMainSubSample())->getCutDef();
	double diff;

	for (auto i : use) {
		if (i != mcSamples[0]->getMainSubSample())
			diff = mcSamples[0]->getSubSample(i)->getCutDef()->getDiff(
					defaultCuts), fit[i]->getFormFactor();
		else
			diff = defaultCuts->getDiff(
					mcSamples[0]->getSubSample((i==0 ? i+1 : i-1))->getCutDef());

		d.addScanValue(diff, fit[i]->getFormFactor(),
				fit[i]->getFormFactorErr(),
				finalDataSample->getSubSample(i)->getSelSize(), mcSamples[0]->getSubSample(i)->getSelSize(), mcSamples.size()>0 ? mcSamples[1]->getSubSample(i)->getSelSize() : 0);
	}
	d.computeUncorrError();
	d.draw();
	d.print();
}

void Drawer::drawCombineStack(std::vector<Sample*> mcSamples,
		std::vector<Sample*> dataSamples, Sample* finalMCSample,
		Sample* finalDataSample) {
	vector<StackRatioDrawer*> d;

	int mainSubSample = mcSamples[0]->getMainSubSample();

	for (auto sample : mcSamples) {
		CombineMCSample *ssample =
				static_cast<CombineMCSample*>(sample->getSubSample(mainSubSample));
		vector<TH1D*> d1 = ssample->getD1();
		for (unsigned int i = 0; i < d1.size(); ++i) {
			if (d.size() <= i){
				d.push_back(new StackRatioDrawer());
				d[i]->setName(d1[i]->GetName());
				d[i]->setTitle(d1[i]->GetTitle());
			}
			d[i]->AddHisto1(d1[i], sample->getLegend(), "f");
		}
	}

	for (auto sample : dataSamples) {
		CombineDataSample *ssample =
				static_cast<CombineDataSample*>(sample->getSubSample(mainSubSample));
		vector<TH1D*> d1 = ssample->getD1();
		for (unsigned int i = 0; i < d1.size(); ++i) {
			if (d.size() <= i)
				d.push_back(new StackRatioDrawer());
			d[i]->AddHisto2(d1[i], sample->getLegend(), "lep");
		}
	}
	vector<TH1D*> data =
			static_cast<CombineDataSample*>(finalDataSample->getSubSample(mainSubSample))->getD1();
	vector<TH1D*> mc =
			static_cast<CombineMCSample*>(finalMCSample->getSubSample(mainSubSample))->getD1();
	for (unsigned int iCanvas = 0; iCanvas < data.size(); ++iCanvas) {
		TH1D* ratio = StackRatioDrawer::buildRatio(mc[iCanvas], data[iCanvas]);
		d[iCanvas]->SetSecondary(ratio);
		d[iCanvas]->draw();
		d[iCanvas]->save();
	}
}

