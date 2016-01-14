/*
 * DistribRatioDrawer.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#include "StackRatioDrawer.h"

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
}

void StackRatioDrawer::AddHisto1(TH1* h, std::string legend) {
	AddHisto(h, legend);
}

void StackRatioDrawer::AddHisto2(TH1* h, std::string legend) {
	fStack2->Add(h);
	fLegend->AddEntry(h, legend.c_str());
}
