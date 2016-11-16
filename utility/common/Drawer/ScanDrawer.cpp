/*
 * ScanDrawer.cpp
 *
 *  Created on: Jan 15, 2016
 *      Author: nlurkin
 */

#include "ScanDrawer.h"

#include <TAttFill.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TLegend.h>
#include <TPad.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <TColor.h>

using namespace std;

ScanDrawer::ScanDrawer() {
	// TODO Auto-generated constructor stub

}

ScanDrawer::~ScanDrawer() {
	// TODO Auto-generated destructor stub
}

void ScanDrawer::generateResult(TPad* pad) {

	//Creation
	TGraphErrors *scanWErr = new TGraphErrors();
	TGraphErrors *scanwUncorr = new TGraphErrors();
	TGraphErrors *scanDefault = new TGraphErrors();
	scanWErr->SetName("ScanWErr");
	scanwUncorr->SetName("scanwUncorr");
	scanDefault->SetName("scanDefault");

	//Points
	for (unsigned int i = 0; i < fScanValues.size(); ++i) {
		scanWErr->SetPoint(i, fScanValues[i], fResultValues[i]);
		scanwUncorr->SetPoint(i, fScanValues[i], fResultValues[i]);
		scanWErr->SetPointError(i, 0, fResultErrors[i]);
		scanwUncorr->SetPointError(i, 0, fUncorrErrors[i]);
	}

	if (*min_element(fScanValues.begin(), fScanValues.end()) == 0) {
		int size = fScanValues.size();
		scanWErr->SetPoint(size, -0.5, 0);
	}

	scanDefault->SetPoint(fDefaultCutValue, fScanValues[fDefaultCutValue],
			0.04);
	scanDefault->SetPointError(fDefaultCutValue, 0, 10);


	//Style
	scanWErr->SetMarkerStyle(23);
	//scanWErr->SetMarkerStyle(0);
	scanWErr->SetMarkerColor(4);
	//scanWErr->SetMarkerColor(TColor::GetColor(137,116,232));
	scanWErr->SetMarkerSize(1);
	scanWErr->SetLineColor(4);
	//scanWErr->SetLineColor(TColor::GetColor(137,116,232));
	scanWErr->SetLineWidth(3);
	scanWErr->SetFillStyle(0);

	scanwUncorr->SetMarkerStyle(0);
	scanwUncorr->SetMarkerColor(4);
	scanwUncorr->SetLineColor(kRed);
	//scanwUncorr->SetLineColor(TColor::GetColor(232,126,116));
	scanwUncorr->SetLineWidth(4);
	scanwUncorr->SetFillStyle(0);

	scanDefault->SetMarkerStyle(0);
	scanDefault->SetLineStyle(7);
	scanDefault->SetLineWidth(3);
	scanDefault->SetLineColor(8);
	scanDefault->SetFillStyle(0);

	THStack* mcStack = getStackFromFile("mcStack");
	THStack* dataStack = getStackFromFile("sigStack");
	TLegend* legend = getLegendFromFile();

	//Plotting
	if (mcStack || dataStack) {
		pad->Divide(1, 2);
		pad->cd(1);
		TPad* pad1 = static_cast<TPad*>(pad->GetPad(1));
		pad1->SetPad(0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(3);
		pad1->SetGrid();
		pad1->cd();               // pad1 becomes the current pad

		if (mcStack)
			mcStack->Draw("HIST");
		if (dataStack)
			dataStack->Draw("SAME E P");
		if (legend)
			legend->Draw();
		if (mcStack)
			mcStack->GetXaxis()->SetRangeUser(scanWErr->GetXaxis()->GetXmin(),
					scanWErr->GetXaxis()->GetXmax());

		TPad* pad2 = static_cast<TPad*>(pad->GetPad(2));
		pad2->SetPad(0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0.2);
		pad2->SetBottomMargin(0.2);
		pad2->SetGrid(); // vertical grid
		pad2->cd();
	} else {
		pad->SetTopMargin(0.08);
		pad->SetBottomMargin(0.2);
		pad->SetRightMargin(0.05);
		pad->SetGrid();
		pad->cd();
	}
	scanWErr->Draw("AP");
	scanDefault->Draw("PSAME");
	scanWErr->Draw("PSAME");
	scanwUncorr->Draw("PSAME");

	scanWErr->SetTitle("FF Slope fit result");
	scanWErr->SetTitle("");

	scanWErr->GetYaxis()->SetTitle("FF Slope a");
	scanWErr->GetYaxis()->SetNdivisions(505);
	scanWErr->GetYaxis()->SetTitleSize(25);
	scanWErr->GetYaxis()->SetTitleFont(43);
	scanWErr->GetYaxis()->SetTitleOffset(0.5);
	scanWErr->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	scanWErr->GetYaxis()->SetLabelSize(15);
	scanWErr->GetYaxis()->SetRangeUser(0.0249, 0.051);

	scanWErr->GetXaxis()->SetTitle("Cut value");
	scanWErr->GetXaxis()->SetTitleSize(25);
	scanWErr->GetXaxis()->SetTitleFont(43);
	scanWErr->GetXaxis()->SetTitleColor(kBlack);
	scanWErr->GetXaxis()->SetLabelFont(43);
	scanWErr->GetXaxis()->SetLabelSize(15);
	scanWErr->GetXaxis()->SetTitleOffset(0.7);

	pad->SaveAs("fitResult.pdf");
	pad->SaveAs("fitResult.png");

}


void ScanDrawer::generateNSelected(TPad* pad) {
	//Creation
	TGraph *scanSelected = new TGraph();
	TGraphErrors *scanDefault = new TGraphErrors();
	scanSelected->SetName("NSelected");
	scanDefault->SetName("scanDefault");

	//Points
	for (unsigned int i = 0; i < fScanValues.size(); ++i) {
		scanSelected->SetPoint(i, fScanValues[i], fNSelected[i]);
		scanDefault->SetPoint(i, fScanValues[i], 0);
	}

	scanDefault->SetPoint(fDefaultCutValue, fScanValues[fDefaultCutValue],
			fNSelected[fDefaultCutValue]);
	scanDefault->SetPointError(fDefaultCutValue, 0, 99999999);

	scanSelected->Draw("AP");
	scanDefault->Draw("PSAME");
	scanSelected->Draw("PLSAME");

	//Style
	scanSelected->SetMarkerStyle(20);
	scanSelected->SetMarkerColor(4);
	scanSelected->SetLineColor(4);
	scanSelected->SetLineWidth(3);
	scanSelected->SetFillStyle(0);

	scanDefault->SetMarkerStyle(0);
	scanDefault->SetLineStyle(7);
	scanDefault->SetLineWidth(3);
	scanDefault->SetLineColor(8);
	scanDefault->SetFillStyle(0);

	//Plotting
	pad->SetTopMargin(0.08);
	pad->SetBottomMargin(0.2);
	pad->SetRightMargin(0.05);
	pad->SetGrid(); // vertical grid
	pad->cd();
	//scanSelected->SetTitle("Number of selected events");
	scanSelected->GetYaxis()->SetTitle("Selected events");
	scanSelected->GetYaxis()->SetTitleSize(25);
	scanSelected->GetYaxis()->SetTitleFont(43);
	scanSelected->GetYaxis()->SetTitleOffset(0.7);
	scanSelected->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	scanSelected->GetYaxis()->SetLabelSize(15);

	scanSelected->GetXaxis()->SetTitle("Cut value");
	scanSelected->GetXaxis()->SetTitleSize(25);
	scanSelected->GetXaxis()->SetTitleFont(43);
	scanSelected->GetXaxis()->SetTitleColor(kBlack);
	scanSelected->GetXaxis()->SetLabelFont(43);
	scanSelected->GetXaxis()->SetLabelSize(15);
	scanSelected->GetXaxis()->SetTitleOffset(0.7);

	pad->SaveAs("nSelected.pdf");
	pad->SaveAs("nSelected.png");

}

void ScanDrawer::draw() {
	TCanvas *c1 = new TCanvas("plot", "", 600, 300);
	generateResult(c1);
	TCanvas *c2 = new TCanvas("Selected", "", 600, 300);
	generateNSelected(c2);
}

void ScanDrawer::addScanValue(double val, double result, double err, int ndata,
		int npi, int nmu) {
	fScanValues.push_back(val);
	fResultValues.push_back(result);
	fResultErrors.push_back(err);
	fNSelected.push_back(ndata);
	fPiSelected.push_back(npi);
	fMuSelected.push_back(nmu);
}

void ScanDrawer::computeUncorrError() {
	double deflt = fResultErrors[fDefaultCutValue];

	for (auto err : fResultErrors)
		fUncorrErrors.push_back(sqrt(fabs(pow(err, 2) - pow(deflt, 2))));
}

void ScanDrawer::print() {
	for (unsigned int i = 0; i < fScanValues.size(); ++i) {
		cout << fixed << setprecision(4) << setw(10) << fScanValues[i]
				<< setw(10) << fResultValues[i] * 100 << setw(10)
				<< fResultErrors[i] * 100 << setw(10) << fUncorrErrors[i] * 100
				<< setw(10) << fNSelected[i] << setw(10) << fPiSelected[i]
				<< setw(10) << fMuSelected[i] << endl;
	}
}

THStack* ScanDrawer::getStackFromFile(std::string name) {
	TFile *fd = TFile::Open("distribRef.root", "READ");
	if (!fd)
		return nullptr;
	THStack *r = static_cast<THStack*>(fd->Get(name.c_str()));
	fd->Close();
	return r;
}

TLegend* ScanDrawer::getLegendFromFile() {
	TFile *fd = TFile::Open("distribRef.root", "READ");
	if (!fd)
		return nullptr;
	TLegend *r = static_cast<TLegend*>(fd->Get("legend"));
	fd->Close();
	return r;
}
