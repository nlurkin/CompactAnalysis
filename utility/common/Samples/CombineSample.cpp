/*
 * CombineMCSample.cpp
 *
 *  Created on: Jan 2, 2016
 *      Author: nlurkin
 */

#include "CombineSample.h"

#include <iostream>
#include "../userinc/funLib.h"
#include <TFile.h>
#include <iomanip>
#include "Functions.h"
#include <TStyle.h>
#include "../Drawer/CombineDrawer.h"

using namespace std;

CombineSample::CombineSample() {
}

CombineSample::CombineSample(int index, ConfigFile* cfg) :
		Sample(index, cfg) {
}

CombineSample::~CombineSample() {
	// TODO Auto-generated destructor stub
}

void CombineSample::doFill(TFile* inputFD, TFile* tempFD) {
	//Get the TTree
	//Input
	ROOTPhysicsEvent *eventBrch = new ROOTPhysicsEvent();
	ROOTBurst *burstBrch = new ROOTBurst();
	ROOTRawEvent *rawBrch = new ROOTRawEvent();
	ROOTCorrectedEvent *corrBrch = new ROOTCorrectedEvent();
	ROOTFileHeader *headerBrch = new ROOTFileHeader();
	ROOTMCEvent *mcEvent = 0;
	TTree *t = (TTree*) inputFD->Get("event");
	TTree *th = (TTree*) inputFD->Get("header");
	if (t->GetListOfBranches()->Contains("mc"))
		mcEvent = new ROOTMCEvent();
	NGeom *geomBrch = new NGeom();

	t->SetBranchAddress("pi0dEvent", &eventBrch);
	t->SetBranchAddress("rawBurst", &burstBrch);
	t->SetBranchAddress("rawEvent", &rawBrch);
	t->SetBranchAddress("corrEvent", &corrBrch);
	th->SetBranchAddress("header", &headerBrch);
	if (mcEvent)
		t->SetBranchAddress("mc", &mcEvent);
	th->SetBranchAddress("geom", &geomBrch);

	tempFD->cd();
	//Set event nb
	int nevt = t->GetEntries();
	int totalChanEvents = 0;
	for (int i = 0; i < th->GetEntries(); i++) {
		th->GetEntry(i);
		totalChanEvents += headerBrch->NProcessedEvents;
	}

	th->GetEntry(0);

	fFitBrch.totEvents += totalChanEvents;
	fFitBrch.selEvents += nevt;

	//Read events and fill histo
	cout << "Filling " << nevt << endl;
	for (int i = 0; i < nevt; ++i) {
		if (i % 10000 == 0)
			cout << setprecision(2) << i * 100. / (double) nevt << "% " << i
					<< "/" << nevt << "\r";
		cout.flush();
		t->GetEntry(i);
		if (!fCfg->testUseRun(burstBrch->nrun, burstBrch->period))
			continue;

		if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch,
				fFitBrch))
			continue;

		fillHisto(eventBrch, rawBrch, corrBrch, mcEvent, geomBrch, burstBrch);
	}
}

void CombineSample::fillHisto(ROOTPhysicsEvent *evt, ROOTRawEvent *rawEvt,
		ROOTCorrectedEvent *corrEvent, ROOTMCEvent *mcEvent, NGeom *rootGeom,
		ROOTBurst *rootBurst, double weight) {
	int i = -1;
	int iMap = -1;
	TVector3 propPos, propPos2, propPos3;

	d1.at(++i)->Fill(evt->kaon.P.M(), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);

	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);

	d1.at(++i)->Fill(evt->pic.vertex.Z(), weight);
	d1.at(++i)->Fill(evt->pic.vertex.Z(), weight);
	d1.at(++i)->Fill(evt->pic.vertex.Z(), weight);
	d1.at(++i)->Fill(rawEvt->vtx[evt->pic.parentVertex].charge, weight);
	d1.at(++i)->Fill(rawEvt->vtx[evt->pic.parentVertex].cda, weight);
	d1.at(++i)->Fill(evt->kaon.P.Perp2(corrEvent->kaonMomentum), weight);
	d1.at(++i)->Fill(evt->kaon.P.Vect().Mag(), weight);

	d1.at(++i)->Fill(evt->pi0.P.M(), weight);

	//Photon
	d1.at(++i)->Fill(evt->gamma.P.E(), weight);
	d1.at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.X(),
			weight);
	d1.at(++i)->Fill(corrEvent->pCluster[evt->gamma.parentCluster].position.Y(),
			weight);
	d1.at(++i)->Fill(
			distance2D(corrEvent->pCluster[evt->gamma.parentCluster].position,
					TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(evt->gamma.P.Vect().Mag(), weight);

	//e+/e-
	d1.at(++i)->Fill(evt->ep.P.Vect().Mag(), weight);
	d1.at(++i)->Fill(evt->ep.P.Vect().Unit().X(), weight);
	d1.at(++i)->Fill(evt->ep.P.Vect().Unit().Y(), weight);
	d1.at(++i)->Fill(evt->ep.P.Vect().Unit().Z(), weight);
	d1.at(++i)->Fill(evt->ep.P.E(), weight);
	d1.at(++i)->Fill(
			corrEvent->pTrack[evt->ep.parentTrack].E
					/ corrEvent->pTrack[evt->ep.parentTrack].p, weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);

	d1.at(++i)->Fill(evt->em.P.Vect().Mag(), weight);
	d1.at(++i)->Fill(evt->em.P.Vect().Unit().X(), weight);
	d1.at(++i)->Fill(evt->em.P.Vect().Unit().Y(), weight);
	d1.at(++i)->Fill(evt->em.P.Vect().Unit().Z(), weight);
	d1.at(++i)->Fill(evt->em.P.E(), weight);
	d1.at(++i)->Fill(
			corrEvent->pTrack[evt->em.parentTrack].E
					/ corrEvent->pTrack[evt->em.parentTrack].p, weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);

	d1.at(++i)->Fill(evt->mee, weight);

	//pi+
	d1.at(++i)->Fill(evt->pic.P.Vect().Mag(), weight);
	d1.at(++i)->Fill(evt->pic.P.Vect().Unit().X(), weight);
	d1.at(++i)->Fill(evt->pic.P.Vect().Unit().Y(), weight);
	d1.at(++i)->Fill(evt->pic.P.Vect().Unit().Z(), weight);
	d1.at(++i)->Fill(evt->pic.P.E(), weight);
	d1.at(++i)->Fill(
			corrEvent->pTrack[evt->pic.parentTrack].E
					/ corrEvent->pTrack[evt->pic.parentTrack].p, weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, propPos2), weight);
	d1.at(++i)->Fill(distance2D(propPos, propPos3), weight);
	d1.at(++i)->Fill(distance2D(propPos2, propPos3), weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	propPos3 = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(distance2D(propPos, propPos2), weight);
	d1.at(++i)->Fill(distance2D(propPos, propPos3), weight);
	d1.at(++i)->Fill(distance2D(propPos2, propPos3), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(
			distance2D(propPos,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos2,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos3,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	propPos3 = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(
			distance2D(propPos,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos2,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos3,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	propPos3 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	d1.at(++i)->Fill(
			distance2D(propPos,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos2,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);
	d1.at(++i)->Fill(
			distance2D(propPos3,
					corrEvent->pCluster[evt->gamma.parentCluster].position),
			weight);

	double ELKr_ep = 0;
	double ELKr_em = 0;
	bool goodPBWall;
	NPhysicsTrack t_ep = corrEvent->pTrack[evt->ep.parentTrack];
	NPhysicsTrack t_em = corrEvent->pTrack[evt->em.parentTrack];

	propPos = propagateAfter(rootGeom->Lkr.z, t_ep, *rawEvt);
	goodPBWall = true;
	if (rootBurst->pbWall && (propPos.Y() > -33.575 && propPos.Y() < -11.850))
		goodPBWall = false;
	if (t_ep.lkr_acc == 0 && goodPBWall
			&& rawEvt->track[t_ep.trackID].dDeadCell > 2.)
		ELKr_ep = t_ep.p;

	propPos = propagateAfter(rootGeom->Lkr.z, t_em, *rawEvt);
	goodPBWall = true;
	if (rootBurst->pbWall && (propPos.Y() > -33.575 && propPos.Y() < -11.850))
		goodPBWall = false;
	if (t_em.lkr_acc == 0 && goodPBWall
			&& rawEvt->track[t_em.trackID].dDeadCell > 2.)
		ELKr_em = t_em.p;

	d1.at(++i)->Fill(ELKr_ep, weight);
	d1.at(++i)->Fill(ELKr_em, weight);
	d1.at(++i)->Fill(evt->gamma.P.E(), weight);
	d1.at(++i)->Fill(ELKr_ep + ELKr_em + evt->gamma.P.E(), weight);

	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	propPos2 = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	i += int((propPos.Y() + 150.) / 30.) * 4;
	d1.at(++i)->Fill(distance2D(propPos, TVector3(0, 0, 0)), weight);
	d1.at(++i)->Fill(propPos.X(), weight);
	d1.at(++i)->Fill(propPos.Y(), weight);
	d1.at(++i)->Fill(distance2D(propPos, propPos2), weight);

	if (mcEvent)
		dMap.at(++iMap)->Fill(mcEvent->xTrue, evt->x, weight);
	else
		++iMap;
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Lkr.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	dMap.at(++iMap)->Fill(
			corrEvent->pCluster[evt->gamma.parentCluster].position.X(),
			corrEvent->pCluster[evt->gamma.parentCluster].position.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[0].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateBefore(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[1].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[2].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->ep.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->em.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagateAfter(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pTrack[evt->pic.parentTrack], *rawEvt);
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
	propPos = propagate(rootGeom->Dch[3].PosChamber.z,
			corrEvent->pCluster[evt->gamma.parentCluster].position,
			evt->gamma.P.Vect());
	dMap.at(++iMap)->Fill(propPos.X(), propPos.Y(), weight);
}

void CombineSample::getHisto(TFile* inputFD, TFile* tempFD, TString name) {
	TH1D* xxx = (TH1D*) inputFD->Get(name);
	xxx->SetName(TString::Format("%s%i", name.Data(), fIndex));
	tempFD->cd();
	d1.push_back((TH1D*) xxx->Clone());
}
void CombineSample::getHisto2(TFile* inputFD, TFile* tempFD, TString name) {
	TH2D* xxx = (TH2D*) inputFD->Get(name);
	xxx->SetName(TString::Format("%s%i", name.Data(), fIndex));
	tempFD->cd();
	dMap.push_back((TH2D*) xxx->Clone());
}

void CombineSample::doGetHisto(TFile* inputFD, TFile* tempFD) {
	getHisto(inputFD, tempFD, "mK");

	getHisto(inputFD, tempFD, "R_DCH1_ep");
	getHisto(inputFD, tempFD, "X_DCH1_ep");
	getHisto(inputFD, tempFD, "Y_DCH1_ep");
	getHisto(inputFD, tempFD, "R_DCH1_em");
	getHisto(inputFD, tempFD, "X_DCH1_em");
	getHisto(inputFD, tempFD, "Y_DCH1_em");
	getHisto(inputFD, tempFD, "R_DCH1_pip");
	getHisto(inputFD, tempFD, "X_DCH1_pip");
	getHisto(inputFD, tempFD, "Y_DCH1_pip");
	getHisto(inputFD, tempFD, "R_DCH1_gamma");
	getHisto(inputFD, tempFD, "X_DCH1_gamma");
	getHisto(inputFD, tempFD, "Y_DCH1_gamma");

	getHisto(inputFD, tempFD, "R_DCH2_ep");
	getHisto(inputFD, tempFD, "X_DCH2_ep");
	getHisto(inputFD, tempFD, "Y_DCH2_ep");
	getHisto(inputFD, tempFD, "R_DCH2_em");
	getHisto(inputFD, tempFD, "X_DCH2_em");
	getHisto(inputFD, tempFD, "Y_DCH2_em");
	getHisto(inputFD, tempFD, "R_DCH2_pip");
	getHisto(inputFD, tempFD, "X_DCH2_pip");
	getHisto(inputFD, tempFD, "Y_DCH2_pip");
	getHisto(inputFD, tempFD, "R_DCH2_gamma");
	getHisto(inputFD, tempFD, "X_DCH2_gamma");
	getHisto(inputFD, tempFD, "Y_DCH2_gamma");

	getHisto(inputFD, tempFD, "R_DCH3_ep");
	getHisto(inputFD, tempFD, "X_DCH3_ep");
	getHisto(inputFD, tempFD, "Y_DCH3_ep");
	getHisto(inputFD, tempFD, "R_DCH3_em");
	getHisto(inputFD, tempFD, "X_DCH3_em");
	getHisto(inputFD, tempFD, "Y_DCH3_em");
	getHisto(inputFD, tempFD, "R_DCH3_pip");
	getHisto(inputFD, tempFD, "X_DCH3_pip");
	getHisto(inputFD, tempFD, "Y_DCH3_pip");
	getHisto(inputFD, tempFD, "R_DCH3_gamma");
	getHisto(inputFD, tempFD, "X_DCH3_gamma");
	getHisto(inputFD, tempFD, "Y_DCH3_gamma");

	getHisto(inputFD, tempFD, "R_DCH4_ep");
	getHisto(inputFD, tempFD, "X_DCH4_ep");
	getHisto(inputFD, tempFD, "Y_DCH4_ep");
	getHisto(inputFD, tempFD, "R_DCH4_em");
	getHisto(inputFD, tempFD, "X_DCH4_em");
	getHisto(inputFD, tempFD, "Y_DCH4_em");
	getHisto(inputFD, tempFD, "R_DCH4_pip");
	getHisto(inputFD, tempFD, "X_DCH4_pip");
	getHisto(inputFD, tempFD, "Y_DCH4_pip");
	getHisto(inputFD, tempFD, "R_DCH4_gamma");
	getHisto(inputFD, tempFD, "X_DCH4_gamma");
	getHisto(inputFD, tempFD, "Y_DCH4_gamma");

	getHisto(inputFD, tempFD, "Zvtx");
	//	getHisto(inputFD, tempFD, "Zvtx_low");
	//	getHisto(inputFD, tempFD, "Zvtx_high");
	getHisto(inputFD, tempFD, "Qvtx");
	getHisto(inputFD, tempFD, "CDAvtx");
	getHisto(inputFD, tempFD, "Pt2");
	getHisto(inputFD, tempFD, "P");

	getHisto(inputFD, tempFD, "Mpi0");

	//Photon
	getHisto(inputFD, tempFD, "gEnergy");
	getHisto(inputFD, tempFD, "gPositionX");
	getHisto(inputFD, tempFD, "gPositionY");
	getHisto(inputFD, tempFD, "gRadius");
	getHisto(inputFD, tempFD, "gP");

	//e+/e-
	getHisto(inputFD, tempFD, "epPMag");
	getHisto(inputFD, tempFD, "epPx");
	getHisto(inputFD, tempFD, "epPy");
	getHisto(inputFD, tempFD, "epPz");
	getHisto(inputFD, tempFD, "epEnergy");
	getHisto(inputFD, tempFD, "epeop");
	getHisto(inputFD, tempFD, "epLKrX");
	getHisto(inputFD, tempFD, "epLKrY");
	getHisto(inputFD, tempFD, "epLKrR");

	getHisto(inputFD, tempFD, "emPMag");
	getHisto(inputFD, tempFD, "emPx");
	getHisto(inputFD, tempFD, "emPy");
	getHisto(inputFD, tempFD, "emPz");
	getHisto(inputFD, tempFD, "emEnergy");
	getHisto(inputFD, tempFD, "emeop");
	getHisto(inputFD, tempFD, "emLKrX");
	getHisto(inputFD, tempFD, "emLKrY");
	getHisto(inputFD, tempFD, "emLKrR");

	getHisto(inputFD, tempFD, "mee");

	//pi+
	getHisto(inputFD, tempFD, "pipPMag");
	getHisto(inputFD, tempFD, "pipPx");
	getHisto(inputFD, tempFD, "pipPy");
	getHisto(inputFD, tempFD, "pipPz");
	getHisto(inputFD, tempFD, "pipEnergy");
	getHisto(inputFD, tempFD, "pieop");

	getHisto(inputFD, tempFD, "t_epem_DCH");
	getHisto(inputFD, tempFD, "t_eppip_DCH");
	getHisto(inputFD, tempFD, "t_empip_DCH");
	getHisto(inputFD, tempFD, "t_epem_LKr");
	getHisto(inputFD, tempFD, "t_eppip_LKr");
	getHisto(inputFD, tempFD, "t_empip_LKr");

	getHisto(inputFD, tempFD, "t_gep_DCH");
	getHisto(inputFD, tempFD, "t_gem_DCH");
	getHisto(inputFD, tempFD, "t_gpip_DCH");
	getHisto(inputFD, tempFD, "t_gep_LKr");
	getHisto(inputFD, tempFD, "t_gem_LKr");
	getHisto(inputFD, tempFD, "t_gpip_LKr");

	getHisto(inputFD, tempFD, "undeft_gep_LKr");
	getHisto(inputFD, tempFD, "undeft_gem_LKr");
	getHisto(inputFD, tempFD, "undeft_gpip_LKr");

	getHisto(inputFD, tempFD, "L3_E_LKr_ep");
	getHisto(inputFD, tempFD, "L3_E_LKr_em");
	getHisto(inputFD, tempFD, "L3_E_LKr_gamma");
	getHisto(inputFD, tempFD, "L3_E_LKr");

	getHisto(inputFD, tempFD, "R_DCH1_ep_0");
	getHisto(inputFD, tempFD, "X_DCH1_ep_0");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_0");
	getHisto(inputFD, tempFD, "t_epem_DCH1_0");
	getHisto(inputFD, tempFD, "R_DCH1_ep_1");
	getHisto(inputFD, tempFD, "X_DCH1_ep_1");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_1");
	getHisto(inputFD, tempFD, "t_epem_DCH1_1");
	getHisto(inputFD, tempFD, "R_DCH1_ep_2");
	getHisto(inputFD, tempFD, "X_DCH1_ep_2");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_2");
	getHisto(inputFD, tempFD, "t_epem_DCH1_2");
	getHisto(inputFD, tempFD, "R_DCH1_ep_3");
	getHisto(inputFD, tempFD, "X_DCH1_ep_3");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_3");
	getHisto(inputFD, tempFD, "t_epem_DCH1_3");
	getHisto(inputFD, tempFD, "R_DCH1_ep_4");
	getHisto(inputFD, tempFD, "X_DCH1_ep_4");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_4");
	getHisto(inputFD, tempFD, "t_epem_DCH1_4");
	getHisto(inputFD, tempFD, "R_DCH1_ep_5");
	getHisto(inputFD, tempFD, "X_DCH1_ep_5");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_5");
	getHisto(inputFD, tempFD, "t_epem_DCH1_5");
	getHisto(inputFD, tempFD, "R_DCH1_ep_6");
	getHisto(inputFD, tempFD, "X_DCH1_ep_6");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_6");
	getHisto(inputFD, tempFD, "t_epem_DCH1_6");
	getHisto(inputFD, tempFD, "R_DCH1_ep_7");
	getHisto(inputFD, tempFD, "X_DCH1_ep_7");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_7");
	getHisto(inputFD, tempFD, "t_epem_DCH1_7");
	getHisto(inputFD, tempFD, "R_DCH1_ep_8");
	getHisto(inputFD, tempFD, "X_DCH1_ep_8");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_8");
	getHisto(inputFD, tempFD, "t_epem_DCH1_8");
	getHisto(inputFD, tempFD, "R_DCH1_ep_9");
	getHisto(inputFD, tempFD, "X_DCH1_ep_9");
	getHisto(inputFD, tempFD, "Y_DCH1_ep_9");
	getHisto(inputFD, tempFD, "t_epem_DCH1_9");
	getHisto2(inputFD, tempFD, "xMap");
	getHisto2(inputFD, tempFD, "LKr_XY_ep");
	getHisto2(inputFD, tempFD, "LKr_XY_em");
	getHisto2(inputFD, tempFD, "LKr_XY_pip");
	getHisto2(inputFD, tempFD, "LKr_XY_gamma");
	getHisto2(inputFD, tempFD, "DCH1_XY_ep");
	getHisto2(inputFD, tempFD, "DCH1_XY_em");
	getHisto2(inputFD, tempFD, "DCH1_XY_pip");
	getHisto2(inputFD, tempFD, "DCH1_XY_gamma");
	getHisto2(inputFD, tempFD, "DCH2_XY_ep");
	getHisto2(inputFD, tempFD, "DCH2_XY_em");
	getHisto2(inputFD, tempFD, "DCH2_XY_pip");
	getHisto2(inputFD, tempFD, "DCH2_XY_gamma");
	getHisto2(inputFD, tempFD, "DCH3_XY_ep");
	getHisto2(inputFD, tempFD, "DCH3_XY_em");
	getHisto2(inputFD, tempFD, "DCH3_XY_pip");
	getHisto2(inputFD, tempFD, "DCH3_XY_gamma");
	getHisto2(inputFD, tempFD, "DCH4_XY_ep");
	getHisto2(inputFD, tempFD, "DCH4_XY_em");
	getHisto2(inputFD, tempFD, "DCH4_XY_pip");
	getHisto2(inputFD, tempFD, "DCH4_XY_gamma");
}

void CombineSample::doWrite() {
	for (auto histo : d1)
		histo->Write();

	for (auto histo : dMap)
		histo->Write();
}

void CombineSample::doSetName() {
//	for (auto histo : d1)
//		histo->SetName(TString::Format("%s%i", histo->GetName(), fIndex));
//
//	for (auto histo : dMap)
//		histo->SetName(TString::Format("%s%i", histo->GetName(), fIndex));
}

void CombineSample::addHisto(TString name, int bins, double min, double max) {
	TH1D* xxx1 = new TH1D(name, "sample 1", bins, min, max);
	xxx1->Sumw2();
	d1.push_back(xxx1);
}

void CombineSample::addHisto(TString name, int binsx, double minx, double maxx,
		int binsy, double miny, double maxy) {
	TH2D* xxx1 = new TH2D(name, "sample 1", binsx, minx, maxx, binsy, miny,
			maxy);
	xxx1->Sumw2();
	dMap.push_back(xxx1);
}

void CombineSample::initHisto(int, double*) {
	//1
	addHisto("mK", 100, 0.47, 0.52);

	//13
	addHisto("R_DCH1_ep", 150, 0, 150);
	addHisto("X_DCH1_ep", 300, -150, 150);
	addHisto("Y_DCH1_ep", 300, -150, 150);
	addHisto("R_DCH1_em", 150, 0, 150);
	addHisto("X_DCH1_em", 300, -150, 150);
	addHisto("Y_DCH1_em", 300, -150, 150);
	addHisto("R_DCH1_pip", 150, 0, 150);
	addHisto("X_DCH1_pip", 300, -150, 150);
	addHisto("Y_DCH1_pip", 300, -150, 150);
	addHisto("R_DCH1_gamma", 150, 0, 150);
	addHisto("X_DCH1_gamma", 300, -150, 150);
	addHisto("Y_DCH1_gamma", 300, -150, 150);

	//25
	addHisto("R_DCH2_ep", 150, 0, 150);
	addHisto("X_DCH2_ep", 300, -150, 150);
	addHisto("Y_DCH2_ep", 300, -150, 150);
	addHisto("R_DCH2_em", 150, 0, 150);
	addHisto("X_DCH2_em", 300, -150, 150);
	addHisto("Y_DCH2_em", 300, -150, 150);
	addHisto("R_DCH2_pip", 150, 0, 150);
	addHisto("X_DCH2_pip", 300, -150, 150);
	addHisto("Y_DCH2_pip", 300, -150, 150);
	addHisto("R_DCH2_gamma", 150, 0, 150);
	addHisto("X_DCH2_gamma", 300, -150, 150);
	addHisto("Y_DCH2_gamma", 300, -150, 150);

	//37
	addHisto("R_DCH3_ep", 150, 0, 150);
	addHisto("X_DCH3_ep", 300, -150, 150);
	addHisto("Y_DCH3_ep", 300, -150, 150);
	addHisto("R_DCH3_em", 150, 0, 150);
	addHisto("X_DCH3_em", 300, -150, 150);
	addHisto("Y_DCH3_em", 300, -150, 150);
	addHisto("R_DCH3_pip", 150, 0, 150);
	addHisto("X_DCH3_pip", 300, -150, 150);
	addHisto("Y_DCH3_pip", 300, -150, 150);
	addHisto("R_DCH3_gamma", 150, 0, 150);
	addHisto("X_DCH3_gamma", 300, -150, 150);
	addHisto("Y_DCH3_gamma", 300, -150, 150);

	//49
	addHisto("R_DCH4_ep", 150, 0, 150);
	addHisto("X_DCH4_ep", 300, -150, 150);
	addHisto("Y_DCH4_ep", 300, -150, 150);
	addHisto("R_DCH4_em", 150, 0, 150);
	addHisto("X_DCH4_em", 300, -150, 150);
	addHisto("Y_DCH4_em", 300, -150, 150);
	addHisto("R_DCH4_pip", 150, 0, 150);
	addHisto("X_DCH4_pip", 300, -150, 150);
	addHisto("Y_DCH4_pip", 300, -150, 150);
	addHisto("R_DCH4_gamma", 150, 0, 150);
	addHisto("X_DCH4_gamma", 300, -150, 150);
	addHisto("Y_DCH4_gamma", 300, -150, 150);

	//54
	addHisto("Zvtx", 100, -2000, 9000);
	addHisto("Zvtx_low", 100, -2100, 0);
	addHisto("Zvtx_high", 100, 5000, 9500);
	addHisto("Qvtx", 10, -5, 5);
	addHisto("CDAvtx", 100, 0, 10);
	addHisto("Pt2", 100, 0, 0.001);
	addHisto("P", 100, 68, 80);

	//55
	addHisto("Mpi0", 60, 0.120, 0.150);

	//Photon
	//60
	addHisto("gEnergy", 80, 0, 80);
	addHisto("gPositionX", 300, -150, 150);
	addHisto("gPositionY", 300, -150, 150);
	addHisto("gRadius", 150, 0, 150);
	addHisto("gP", 60, 0, 60);

	//e+/e-
	//69
	addHisto("epPMag", 60, 0, 60);
	addHisto("epPx", 60, -0.015, 0.015);
	addHisto("epPy", 60, -0.015, 0.015);
	addHisto("epPz", 110, 0, 1.1);
	addHisto("epEnergy", 60, 0, 60);
	addHisto("epeop", 100, 0, 1.5);
	addHisto("epLKrX", 300, -150, 150);
	addHisto("epLKrY", 300, -150, 150);
	addHisto("epLKrR", 300, -150, 150);

	//78
	addHisto("emPMag", 60, 0, 60);
	addHisto("emPx", 60, -0.015, 0.015);
	addHisto("emPy", 60, -0.015, 0.015);
	addHisto("emPz", 110, 0, 1.1);
	addHisto("emEnergy", 60, 0, 60);
	addHisto("emeop", 100, 0, 1.5);
	addHisto("emLKrX", 300, -150, 150);
	addHisto("emLKrY", 300, -150, 150);
	addHisto("emLKrR", 300, -150, 150);

	//79
	addHisto("mee", 140, 0, 0.14);

	//pi+
	//85
	addHisto("pipPMag", 60, 0, 60);
	addHisto("pipPx", 60, -0.015, 0.015);
	addHisto("pipPy", 60, -0.015, 0.015);
	addHisto("pipPz", 110, 0, 1.1);
	addHisto("pipEnergy", 60, 0, 60);
	addHisto("pieop", 100, 0, 1.5);

	addHisto("t_epem_DCH", 150, 0, 150);
	addHisto("t_eppip_DCH", 150, 0, 150);
	addHisto("t_empip_DCH", 150, 0, 150);
	addHisto("t_epem_LKr", 400, 0, 400);
	addHisto("t_eppip_LKr", 400, 0, 400);
	addHisto("t_empip_LKr", 400, 0, 400);

	addHisto("t_gep_DCH", 150, 0, 150);
	addHisto("t_gem_DCH", 150, 0, 150);
	addHisto("t_gpip_DCH", 150, 0, 150);
	addHisto("t_gep_LKr", 400, 0, 400);
	addHisto("t_gem_LKr", 400, 0, 400);
	addHisto("t_gpip_LKr", 400, 0, 400);

	addHisto("undeft_gep_LKr", 400, 0, 400);
	addHisto("undeft_gem_LKr", 400, 0, 400);
	addHisto("undeft_gpip_LKr", 400, 0, 400);

	addHisto("L3_E_LKr_ep", 160, 0, 80);
	addHisto("L3_E_LKr_em", 160, 0, 80);
	addHisto("L3_E_LKr_gamma", 160, 0, 80);
	addHisto("L3_E_LKr", 160, 0, 80);

	addHisto("R_DCH1_ep_0", 150, 0, 150);
	addHisto("X_DCH1_ep_0", 300, -150, 150);
	addHisto("Y_DCH1_ep_0", 300, -150, 150);
	addHisto("t_epem_DCH1_0", 150, 0, 150);
	addHisto("R_DCH1_ep_1", 150, 0, 150);
	addHisto("X_DCH1_ep_1", 300, -150, 150);
	addHisto("Y_DCH1_ep_1", 300, -150, 150);
	addHisto("t_epem_DCH1_1", 150, 0, 150);
	addHisto("R_DCH1_ep_2", 150, 0, 150);
	addHisto("X_DCH1_ep_2", 300, -150, 150);
	addHisto("Y_DCH1_ep_2", 300, -150, 150);
	addHisto("t_epem_DCH1_2", 150, 0, 150);
	addHisto("R_DCH1_ep_3", 150, 0, 150);
	addHisto("X_DCH1_ep_3", 300, -150, 150);
	addHisto("Y_DCH1_ep_3", 300, -150, 150);
	addHisto("t_epem_DCH1_3", 150, 0, 150);
	addHisto("R_DCH1_ep_4", 150, 0, 150);
	addHisto("X_DCH1_ep_4", 300, -150, 150);
	addHisto("Y_DCH1_ep_4", 300, -150, 150);
	addHisto("t_epem_DCH1_4", 150, 0, 150);
	addHisto("R_DCH1_ep_5", 150, 0, 150);
	addHisto("X_DCH1_ep_5", 300, -150, 150);
	addHisto("Y_DCH1_ep_5", 300, -150, 150);
	addHisto("t_epem_DCH1_5", 150, 0, 150);
	addHisto("R_DCH1_ep_6", 150, 0, 150);
	addHisto("X_DCH1_ep_6", 300, -150, 150);
	addHisto("Y_DCH1_ep_6", 300, -150, 150);
	addHisto("t_epem_DCH1_6", 150, 0, 150);
	addHisto("R_DCH1_ep_7", 150, 0, 150);
	addHisto("X_DCH1_ep_7", 300, -150, 150);
	addHisto("Y_DCH1_ep_7", 300, -150, 150);
	addHisto("t_epem_DCH1_7", 150, 0, 150);
	addHisto("R_DCH1_ep_8", 150, 0, 150);
	addHisto("X_DCH1_ep_8", 300, -150, 150);
	addHisto("Y_DCH1_ep_8", 300, -150, 150);
	addHisto("t_epem_DCH1_8", 150, 0, 150);
	addHisto("R_DCH1_ep_9", 150, 0, 150);
	addHisto("X_DCH1_ep_9", 300, -150, 150);
	addHisto("Y_DCH1_ep_9", 300, -150, 150);
	addHisto("t_epem_DCH1_9", 150, 0, 150);

	addHisto("xMap", 1000, 0, 1, 1000, 0, 1);
	addHisto("LKr_XY_ep", 600, -120, 120, 600, -120, 120);
	addHisto("LKr_XY_em", 600, -120, 120, 600, -120, 120);
	addHisto("LKr_XY_pip", 600, -120, 120, 600, -120, 120);
	addHisto("LKr_XY_gamma", 600, -120, 120, 600, -120, 120);
	addHisto("DCH1_XY_ep", 500, -100, 100, 500, -100, 100);
	addHisto("DCH1_XY_em", 500, -100, 100, 500, -100, 100);
	addHisto("DCH1_XY_pip", 500, -100, 100, 500, -100, 100);
	addHisto("DCH1_XY_gamma", 500, -100, 100, 500, -100, 100);
	addHisto("DCH2_XY_ep", 500, -100, 100, 500, -100, 100);
	addHisto("DCH2_XY_em", 500, -100, 100, 500, -100, 100);
	addHisto("DCH2_XY_pip", 500, -100, 100, 500, -100, 100);
	addHisto("DCH2_XY_gamma", 500, -100, 100, 500, -100, 100);
	addHisto("DCH3_XY_ep", 500, -100, 100, 500, -100, 100);
	addHisto("DCH3_XY_em", 500, -100, 100, 500, -100, 100);
	addHisto("DCH3_XY_pip", 500, -100, 100, 500, -100, 100);
	addHisto("DCH3_XY_gamma", 500, -100, 100, 500, -100, 100);
	addHisto("DCH4_XY_ep", 500, -100, 100, 500, -100, 100);
	addHisto("DCH4_XY_em", 500, -100, 100, 500, -100, 100);
	addHisto("DCH4_XY_pip", 500, -100, 100, 500, -100, 100);
	addHisto("DCH4_XY_gamma", 500, -100, 100, 500, -100, 100);
}

void CombineSample::scale() {
	cout << "Rescaling" << endl;
	cout << fBr << " " << fFitBrch.selEvents << " " << fFitBrch.totEvents
			<< endl;
	// Rescale histo
	for (unsigned int i = 0; i < d1.size(); ++i) {
		Sample::scale(d1.at(i), 1.);
	}
	for (unsigned int i = 0; i < dMap.size(); ++i) {
		Sample::scale(dMap.at(i), 1.);
	}
}

double CombineSample::getFFIntegral(double) {
	return 0;
}

void CombineSample::setPlotStyle(std::vector<int> color) {
	for(auto plot : d1)
		plot->SetFillColor(gStyle->GetColorPalette(color[0]));

	for(auto plot : dMap)
		plot->SetFillColor(gStyle->GetColorPalette(color[0]));
}

void CombineSample::populateStack(HistoDrawer *drawer) {
	CombineDrawer *myDrawer = static_cast<CombineDrawer*>(drawer);

	for (unsigned int i = 0; i < d1.size(); ++i) {
		myDrawer->addHistoMC(i, (TH1D*) (d1[i]->Clone()), fLegend);
	}
}

void CombineSample::populateFit(HistoDrawer *, double, double) {
	//FitResultDrawer *myDrawer = static_cast<FitResultDrawer*>(drawer);
}

TH1D* CombineSample::getMainHisto() {
	return nullptr;
}

