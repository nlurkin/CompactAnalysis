/*
 * FitResultDrawer.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: nlurkin
 */

#include "FitResultDrawer.h"

#include <TCanvas.h>
#include <TPaveText.h>
#include <TF1.h>

using namespace std;

FitResultDrawer::FitResultDrawer() {
	fFit = new THStack("MC", "Stack MC");
	fFitSig = new THStack("Signal", "Stack Signal");
	fLegFit = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegFitSig = new TLegend(.75, 0.6, 0.98, 0.82);
}

FitResultDrawer::~FitResultDrawer() {
	// TODO Auto-generated destructor stub
}

void FitResultDrawer::draw() {
	TCanvas *c2 = new TCanvas(("fit" + fTitle).c_str(), ("Procedure " + fTitle).c_str());
	TPaveText *fitR = new TPaveText(0.47, 0.78, 0.77, 0.94, "NDC BR");
	fitR->AddText("Fit result");
	fitR->AddLine(0., 0.7, 1., 0.7);
	fitR->AddText(Form("G = %f #pm %f", fNorm, fA));
	fitR->AddText(Form("FF = %f #pm %f", fA, fAErr));
	fitR->SetTextAlign(12);

	TH1D* ratio = buildRatio(fNBins, fBinning, fMC, fSig);
	TF1 *f = new TF1("f", "[0]*(1+[1]*2.0*x+[1]*[1]*x*x)", 0., 1.);
	f->SetParameter(0, fNorm);
	f->SetParameter(1, fA);
	f->SetLineColor(kRed);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(3);
	pad1->SetGrid();
	pad1->Draw();
	pad1->cd();
	fFitSig->Draw("S E P");
	fFit->Draw("HIST SAME");
	fFitSig->Draw("SAMES E P");
	fLegFit->Draw();
	fitR->Draw();
	pad1->SetLogx(true);
	pad1->SetLogy(true);
	c2->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();
	//ratio->Fit("f", "R");
	//ratio->SetStats(0);      // No statistics on lower plot
	ratio->Draw("ep");
	f->Draw("LSAME");
	ratio->SetMarkerColor(kRed);
	pad2->SetLogx(true);

	// Y axis mc plot settings
	//stack->GetYaxis()->SetTitleSize(20);
	//stack->GetYaxis()->SetTitleFont(43);
	//stack->GetYaxis()->SetTitleOffset(1.55);

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

}
