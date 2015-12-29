/*
 * InputFitDrawer.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: nlurkin
 */

#include "InputFitDrawer.h"

#include <TCanvas.h>


InputFitDrawer::InputFitDrawer() {
	fStd1 = new THStack("std1", "Stack FF1");
	fStdx = new THStack("stdx", "Stack FFx");
	fStdxx = new THStack("stdxx", "Stack FFxx");
	fStdNew = new THStack("stdNew", "Stack FFNew");
	fStdAlpha = new THStack("stdAlpha", "Stack Alpha");
	fStdBeta = new THStack("stdBeta", "Stack Beta");
	fStdGamma = new THStack("stdGamma", "Stack Gamma");
	fSig = new THStack("signal", "Stack Signal");
	fLeg1 = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegx = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegxx = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegNew = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegAlpha = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegBeta = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegGamma = new TLegend(.75, 0.6, 0.98, 0.82);
	fLegSig = new TLegend(.75, 0.6, 0.98, 0.82);
}

InputFitDrawer::~InputFitDrawer() {
	// TODO Auto-generated destructor stub
}

void InputFitDrawer::draw() {
	//Plot input histo
	TCanvas *c1 = new TCanvas("cInput", "Inputs", 1600, 800);
	c1->Divide(2, 2);
	c1->cd(4);
	fSig->Draw();
	c1->cd(1);
	fStd1->Draw("HIST");
	fLeg1->Draw();
	c1->cd(2);
	fStdx->Draw("HIST");
	fLegx->Draw();
	c1->cd(3);
	fStdxx->Draw("HIST");
	fLegxx->Draw();

	new TCanvas("cInputsNew", "Inputs New", 1600, 800);
	fSig->DrawClone();
	fStdNew->Draw("HIST");
	fLegNew->Draw();

	TCanvas *c3 = new TCanvas("cABG", "Alpha, Beta, gamma", 1600, 800);
	c3->Divide(2, 2);
	c3->cd(1);
	fStdAlpha->Draw("HIST");
	fLegAlpha->Draw();
	c3->cd(2);
	fStdBeta->Draw("HIST");
	fLegBeta->Draw();
	c3->cd(3);
	fStdGamma->Draw("HIST");
	fLegGamma->Draw();
}
