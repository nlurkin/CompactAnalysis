/*
 * CombineDrawer.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: nlurkin
 */

#include "CombineDrawer.h"
#include <TCanvas.h>
#include <iostream>

using namespace std;
CombineDrawer::CombineDrawer() {
	fLegend = new TLegend(.75, 0.6, 0.98, 0.82);
}

CombineDrawer::~CombineDrawer() {
	// TODO Auto-generated destructor stub
}

void CombineDrawer::draw() {
	for (unsigned int iCanvas = 0; iCanvas < fStack.size(); iCanvas++) {
		cout << fStack[iCanvas]->GetTitle() << endl;
		TCanvas *c1 = new TCanvas(TString::Format("c%i", iCanvas),
				fStack[iCanvas]->GetTitle());
		TH1D* ratio = buildRatio(-1, nullptr, fSum[iCanvas], fDataSum[iCanvas]);

		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(3);
		pad1->SetGrid();
		pad1->Draw();             // Draw the upper pad: pad1
		pad1->cd();               // pad1 becomes the current pad
		fStack[iCanvas]->Draw("HIST");               // Draw h1
		fDataStack[iCanvas]->Draw("SAME E P");         // Draw h2 on top of h1
		//mc->GetYaxis()->SetLabelSize(0.05);
		fLegend->Draw();

		c1->cd();
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0);
		pad2->SetBottomMargin(0.2);
		pad2->SetGrid(); // vertical grid
		pad2->Draw();
		pad2->cd();
		ratio->SetStats(0);      // No statistics on lower plot
		ratio->Draw("ep");
		ratio->SetMarkerColor(kRed);
		ratio->GetYaxis()->SetRangeUser(0.89, 1.11);

		// Y axis mc plot settings
		fStack[iCanvas]->GetYaxis()->SetTitleSize(20);
		fStack[iCanvas]->GetYaxis()->SetTitleFont(43);
		fStack[iCanvas]->GetYaxis()->SetTitleOffset(1.55);

		// Ratio plot (ratio) settings
		ratio->SetTitle(""); // Remove the ratio title

		// Y axis ratio plot settings
		ratio->GetYaxis()->SetTitle("Data/MC");
		ratio->GetYaxis()->SetNdivisions(505);
		ratio->GetYaxis()->SetTitleSize(20);
		ratio->GetYaxis()->SetTitleFont(43);
		ratio->GetYaxis()->SetTitleOffset(1);
		ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		ratio->GetYaxis()->SetLabelSize(15);

		// X axis ratio plot settings
		ratio->GetXaxis()->SetTitleSize(20);
		ratio->GetXaxis()->SetTitleFont(43);
		ratio->GetXaxis()->SetTitleOffset(4.);
		ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		ratio->GetXaxis()->SetLabelSize(15);

		c1->SaveAs(TString(fStack[iCanvas]->GetName()) + ".png");

		c1->Close();
		delete c1;
		++iCanvas;
	}
}

void CombineDrawer::addHistoMC(unsigned int index, TH1* histo) {
	THStack *hStack;
	TH1D* sumHisto;
	if (fStack.size() <= index) {
		hStack = new THStack(histo->GetName(), histo->GetTitle());
		sumHisto = new TH1D(Form("%s_mcsum", histo->GetName()),
				histo->GetTitle(), histo->GetNbinsX(),
				histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
		fStack.push_back(hStack);
		fSum.push_back(sumHisto);
	} else {
		hStack = fStack[index];
		sumHisto = fSum[index];
	}

	hStack->Add(histo);
	sumHisto->Add(histo, 1.);
}

void CombineDrawer::addLegendMC(TH1* histo, string legend) {
	fLegend->AddEntry(histo, legend.c_str(), "f");
}

void CombineDrawer::addHistoData(unsigned int index, TH1* histo) {
	THStack *hStack;
	TH1D* sumHisto;
	if (fDataStack.size() <= index) {
		hStack = new THStack(histo->GetName(), histo->GetTitle());
		sumHisto = new TH1D(Form("%s_datasum", histo->GetName()),
				histo->GetTitle(), histo->GetNbinsX(),
				histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
		fDataStack.push_back(hStack);
		fDataSum.push_back(sumHisto);
	} else {
		hStack = fDataStack[index];
		sumHisto = fDataSum[index];
	}

	hStack->Add(histo);
	sumHisto->Add(histo, 1.);
}

void CombineDrawer::addLegendData(TH1* histo, string legend) {
	fLegend->AddEntry(histo, legend.c_str(), "lep");
}
