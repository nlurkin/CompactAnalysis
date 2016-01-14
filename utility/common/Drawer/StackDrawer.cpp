/*
 * StackDrawer.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: nlurkin
 */

#include "StackDrawer.h"

StackDrawer::StackDrawer() :
fCanvas(nullptr){
	fStack = new THStack("stack", "stack");
	fLegend = new TLegend(.75, 0.6, 0.98, 0.82);
}

StackDrawer::~StackDrawer() {
}


void StackDrawer::save() {
	if (fCanvas)
		fCanvas->SaveAs(TString(fStack->GetName()) + ".png");
}

void StackDrawer::generate(TPad *pad) {
	pad->cd();
	fStack->Draw("HIST");
	fLegend->Draw();
}

void StackDrawer::draw() {
}

void StackDrawer::free() {
	if (fCanvas) {
		fCanvas->Close();
		delete fCanvas;
	}
}

void StackDrawer::AddHisto(TH1* h, std::string legend, std::string option) {
	fStack->Add(h);
	fLegend->AddEntry(h, legend.c_str(), option.c_str());
}
