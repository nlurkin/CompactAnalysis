/*
 * DistribRatioDrawer.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#include "StackRatioDrawer.h"

#include <Rtypes.h>
#include <TAttFill.h>
#include <TAttMarker.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TPad.h>
#include <TString.h>
#include <cstdlib>
#include <string>

using namespace std;

StackRatioDrawer::StackRatioDrawer() :
fSecondary(nullptr){
	fStack2 = new THStack("Signal", "Signal Stack");
}

StackRatioDrawer::~StackRatioDrawer() {
}

void StackRatioDrawer::draw(){
	fCanvas = new TCanvas(fName.c_str(), fTitle.c_str());
	generate(fCanvas);
}
void StackRatioDrawer::generate(TPad* pad) {
	pad->Divide(1, 2);
	pad->cd(1);
	TPad* pad1 = static_cast<TPad*>(pad->GetPad(1));
	pad1->SetPad(0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(3);
	pad1->SetGrid();
	pad1->cd();               // pad1 becomes the current pad
	fStack->Draw("HIST");    // Draw h1
	fStack2->Draw("SAME E P");    // Draw h2 on top of h1
	fLegend->Draw();

	TPad* pad2 = static_cast<TPad*>(pad->GetPad(2));
	pad2->SetPad(0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->cd();
	fSecondary->SetStats(0);      // No statistics on lower plot
	fSecondary->Draw("ep");
	fSecondary->SetMarkerColor(kRed);
	fSecondary->GetYaxis()->SetRangeUser(0.89, 1.11);

	// Y axis mc plot settings
	fStack->GetYaxis()->SetTitleSize(20);
	fStack->GetYaxis()->SetTitleFont(43);
	fStack->GetYaxis()->SetTitleOffset(1.55);

	// Ratio plot (ratio) settings
	fSecondary->SetTitle(""); // Remove the ratio title

	// Y axis ratio plot settings
	fSecondary->GetYaxis()->SetTitle("Data/MC");
	fSecondary->GetYaxis()->SetNdivisions(505);
	fSecondary->GetYaxis()->SetTitleSize(20);
	fSecondary->GetYaxis()->SetTitleFont(43);
	fSecondary->GetYaxis()->SetTitleOffset(1);
	fSecondary->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	fSecondary->GetYaxis()->SetLabelSize(15);

	// X axis ratio plot settings
	fSecondary->GetXaxis()->SetTitleSize(20);
	fSecondary->GetXaxis()->SetTitleFont(43);
	fSecondary->GetXaxis()->SetTitleOffset(4.);
	fSecondary->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	fSecondary->GetXaxis()->SetLabelSize(15);


	TFile *fd = TFile::Open(Form("%s.root", fName.c_str()), "UPDATE");
	fd->cd();
	fStack->Write("mcStack");
	fStack2->Write("sigStack");
	fLegend->Write("legend");
	fd->Close();
}

void StackRatioDrawer::AddHisto1(TH1* h, std::string legend, std::string option) {
	AddHisto(h, legend, option.c_str());
}

void StackRatioDrawer::AddHisto2(TH1* h, std::string legend, std::string option) {
	fStack2->Add(h);
	fLegend->AddEntry(h, legend.c_str(), option.c_str());
}


TH1D* StackRatioDrawer::buildRatio(TH1D* mc, TH1D* data) {
	TH1D* r;
	int rnd = rand() % 99999;
	if(data->GetXaxis()->GetXbins()->fArray)
		r = new TH1D(Form("ratio%i", rnd), "ratio", data->GetXaxis()->GetNbins(), data->GetXaxis()->GetXbins()->fArray);
	else
		r = new TH1D(Form("ratio%i", rnd), "ratio", data->GetXaxis()->GetNbins(), data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());

	mc->SetFillColor(8);

	mc->Sumw2();
	r->Sumw2();
	r->Divide(data, mc, 1, 1, "B");
	return r;
}
