/*
 * SubSample.cpp
 *
 *  Created on: Jan 11, 2016
 *      Author: nlurkin
 */

#include "SubSample.h"
#include <iostream>
#include <TTree.h>

using namespace std;

void initFitStruct(fitStruct &s) {
	s.n1 = 0;
	s.nx = 0;
	s.nxx = 0;
	s.selEvents = 0;
	s.totEvents = 0;
}

void sumTreeFitStruct(fitStruct &in, TTree *t, fitStruct &out, double factor) {
	for (int i = 0; i < t->GetEntries(); i++) {
		t->GetEntry(i);
		out.n1 += factor * in.n1;
		out.nx += factor * in.nx;
		out.nxx += factor * in.nxx;
		out.selEvents += factor * in.selEvents;
		out.totEvents += factor * in.totEvents;
	}

	std::cout << "Total events: \t" << out.totEvents << std::endl;
	std::cout << "Sel.  events: \t" << out.selEvents << std::endl;
	std::cout << "n1    events: \t" << out.n1 << std::endl;
	std::cout << "nx    events: \t" << out.nx << std::endl;
	std::cout << "nxx   events: \t" << out.nxx << std::endl;
}

SubSample::SubSample() :
		fFitTree(nullptr) {
	fFitBrch.n1 = 0;
	fFitBrch.nx = 0;
	fFitBrch.nxx = 0;
	fFitBrch.selEvents = 0;
	fFitBrch.totEvents = 0;
}

SubSample::~SubSample() {
	// TODO Auto-generated destructor stub
}

void SubSample::initNewFile(int totalChanEvents, int selEvents) {
	fFitBrch.totEvents += totalChanEvents;
	fFitBrch.selEvents += selEvents;
}

void SubSample::scale(TH1 *histo, double scaleFactor) {
	cout << "Scaling " << histo->GetName() << " " << histo->GetEntries() << " "
			<< histo->Integral() << " " << fBr << " " << fFitBrch.totEvents
			<< " " << scaleFactor << endl;
	histo->Scale(fBr / (fFitBrch.totEvents * scaleFactor));
}

SubSample *SubSample::Add(const SubSample *other) {
	fFitBrch += other->fFitBrch;

	return this;
}

void SubSample::initOutput() {
	fFitTree = new TTree("fitStruct", "fitStruct tree");
	fFitTree->Branch("fitStruct", &fFitBrch, "totEvents/I:selEvents:n1:nx:nxx");
	fFitBrch.n1 = 0;
	fFitBrch.nx = 0;
	fFitBrch.nxx = 0;
	fFitBrch.selEvents = 0;
	fFitBrch.totEvents = 0;
}

void SubSample::writeTree() {
	fFitTree->Fill();
	fFitTree->Write();
}
