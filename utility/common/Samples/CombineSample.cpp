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

using namespace std;

CombineSample::CombineSample() {
}

CombineSample::~CombineSample() {
}

void CombineSample::processEvent(ROOTPhysicsEvent *eventBrch,
		ROOTBurst *burstBrch, ROOTRawEvent *rawBrch,
		ROOTCorrectedEvent *corrBrch, ROOTFileHeader *, ROOTMCEvent *mcEvent,
		NGeom *geomBrch, std::vector<bool> *cutsPass, const ConfigFile *cfg,
		const RunWeights * weights) {

	if (cutsPass) {
		if (!cutsPass->at(fScanID)) {
			fFitBrch.selEvents--;
			return;
		}
	}
	if (!cfg->testUseRun(burstBrch->nrun, burstBrch->period))
		return;

	if (!testAdditionalCondition(eventBrch, corrBrch, geomBrch, rawBrch, burstBrch,
			fFitBrch))
		return;

	fillHisto(eventBrch, rawBrch, corrBrch, mcEvent, geomBrch, burstBrch, weights);
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

void CombineSample::getHisto(TDirectory* inputFD, TFile* tempFD, TString name) {
	TH1D* xxx = (TH1D*) inputFD->Get(name);
	xxx->SetName(TString::Format("%s%i", name.Data(), fIndex));
	tempFD->cd();
	d1.push_back((TH1D*) xxx->Clone());
}
void CombineSample::getHisto2(TDirectory* inputFD, TFile* tempFD, TString name) {
	TH2D* xxx = (TH2D*) inputFD->Get(name);
	xxx->SetName(TString::Format("%s%i", name.Data(), fIndex));
	tempFD->cd();
	dMap.push_back((TH2D*) xxx->Clone());
}

void CombineSample::doGetHisto(TDirectory* inputFD, TFile* tempFD) {
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
	getHisto(inputFD, tempFD, "Zvtx_low");
	getHisto(inputFD, tempFD, "Zvtx_high");
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

//	getHisto2(inputFD, tempFD, "xMap");
//	getHisto2(inputFD, tempFD, "LKr_XY_ep");
//	getHisto2(inputFD, tempFD, "LKr_XY_em");
//	getHisto2(inputFD, tempFD, "LKr_XY_pip");
//	getHisto2(inputFD, tempFD, "LKr_XY_gamma");
//	getHisto2(inputFD, tempFD, "DCH1_XY_ep");
//	getHisto2(inputFD, tempFD, "DCH1_XY_em");
//	getHisto2(inputFD, tempFD, "DCH1_XY_pip");
//	getHisto2(inputFD, tempFD, "DCH1_XY_gamma");
//	getHisto2(inputFD, tempFD, "DCH2_XY_ep");
//	getHisto2(inputFD, tempFD, "DCH2_XY_em");
//	getHisto2(inputFD, tempFD, "DCH2_XY_pip");
//	getHisto2(inputFD, tempFD, "DCH2_XY_gamma");
//	getHisto2(inputFD, tempFD, "DCH3_XY_ep");
//	getHisto2(inputFD, tempFD, "DCH3_XY_em");
//	getHisto2(inputFD, tempFD, "DCH3_XY_pip");
//	getHisto2(inputFD, tempFD, "DCH3_XY_gamma");
//	getHisto2(inputFD, tempFD, "DCH4_XY_ep");
//	getHisto2(inputFD, tempFD, "DCH4_XY_em");
//	getHisto2(inputFD, tempFD, "DCH4_XY_pip");
//	getHisto2(inputFD, tempFD, "DCH4_XY_gamma");
}

void CombineSample::doWrite() {
	for (auto histo : d1)
		histo->Write();

	for (auto histo : dMap)
		histo->Write();
}

void CombineSample::doSetName() {
//	renameHisto();
}

void CombineSample::renameHisto() {
	for (auto histo : d1)
		histo->SetName(TString::Format("%s%i", histo->GetName(), fIndex));

	for (auto histo : dMap)
		histo->SetName(TString::Format("%s%i", histo->GetName(), fIndex));
}

void CombineSample::addHisto(TString name, int bins, double min, double max) {
	TH1D* xxx1 = new TH1D(name, name, bins, min, max);
	xxx1->Sumw2();
	d1.push_back(xxx1);
}

void CombineSample::addHisto(TString name, int binsx, double minx, double maxx,
		int binsy, double miny, double maxy) {
	TH2D* xxx1 = new TH2D(name, name, binsx, minx, maxx, binsy, miny, maxy);
	xxx1->Sumw2();
	dMap.push_back(xxx1);
}

void CombineSample::initHisto(int, double*, const ConfigFile *) {
	//1
	addHisto("mK", 200, 0.45, 0.53);

	//13
	addHisto("R_DCH1_ep", 	75, 0, 150);
	addHisto("X_DCH1_ep", 	150, -150, 150);
	addHisto("Y_DCH1_ep",	150, -150, 150);
	addHisto("R_DCH1_em", 	75, 0, 150);
	addHisto("X_DCH1_em", 	150, -150, 150);
	addHisto("Y_DCH1_em", 	150, -150, 150);
	addHisto("R_DCH1_pip", 	75, 0, 150);
	addHisto("X_DCH1_pip", 	150, -150, 150);
	addHisto("Y_DCH1_pip", 	150, -150, 150);
	addHisto("R_DCH1_gamma",75, 0, 150);
	addHisto("X_DCH1_gamma",150, -150, 150);
	addHisto("Y_DCH1_gamma",150, -150, 150);

	//25
	addHisto("R_DCH2_ep", 	75,  0, 150);
	addHisto("X_DCH2_ep", 	150, -150, 150);
	addHisto("Y_DCH2_ep", 	150, -150, 150);
	addHisto("R_DCH2_em", 	75,  0, 150);
	addHisto("X_DCH2_em", 	150, -150, 150);
	addHisto("Y_DCH2_em", 	150, -150, 150);
	addHisto("R_DCH2_pip", 	75,  0, 150);
	addHisto("X_DCH2_pip", 	150, -150, 150);
	addHisto("Y_DCH2_pip", 	150, -150, 150);
	addHisto("R_DCH2_gamma",75,  0, 150);
	addHisto("X_DCH2_gamma",150, -150, 150);
	addHisto("Y_DCH2_gamma",150, -150, 150);

	//37
	addHisto("R_DCH3_ep", 	75,  0, 150);
	addHisto("X_DCH3_ep", 	150, -150, 150);
	addHisto("Y_DCH3_ep", 	150, -150, 150);
	addHisto("R_DCH3_em", 	75,  0, 150);
	addHisto("X_DCH3_em", 	150, -150, 150);
	addHisto("Y_DCH3_em", 	150, -150, 150);
	addHisto("R_DCH3_pip", 	75,  0, 150);
	addHisto("X_DCH3_pip", 	150, -150, 150);
	addHisto("Y_DCH3_pip", 	150, -150, 150);
	addHisto("R_DCH3_gamma",75,  0, 150);
	addHisto("X_DCH3_gamma",150, -150, 150);
	addHisto("Y_DCH3_gamma",150, -150, 150);

	//49
	addHisto("R_DCH4_ep", 	75,  0, 150);
	addHisto("X_DCH4_ep", 	150, -150, 150);
	addHisto("Y_DCH4_ep", 	150, -150, 150);
	addHisto("R_DCH4_em", 	75,  0, 150);
	addHisto("X_DCH4_em", 	150, -150, 150);
	addHisto("Y_DCH4_em", 	150, -150, 150);
	addHisto("R_DCH4_pip", 	75,  0, 150);
	addHisto("X_DCH4_pip", 	150, -150, 150);
	addHisto("Y_DCH4_pip", 	150, -150, 150);
	addHisto("R_DCH4_gamma",75,  0, 150);
	addHisto("X_DCH4_gamma",150, -150, 150);
	addHisto("Y_DCH4_gamma",150, -150, 150);

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
	addHisto("gPositionX", 150, -150, 150);
	addHisto("gPositionY", 150, -150, 150);
	addHisto("gRadius", 75, 0, 150);
	addHisto("gP", 60, 0, 60);

	//e+/e-
	//69
	addHisto("epPMag", 60, 0, 60);
	addHisto("epPx", 60, -0.015, 0.015);
	addHisto("epPy", 60, -0.015, 0.015);
	addHisto("epPz", 110, 0, 1.1);
	addHisto("epEnergy", 60, 0, 60);
	addHisto("epeop", 100, 0, 1.5);
	addHisto("epLKrX", 150, -150, 150);
	addHisto("epLKrY", 150, -150, 150);
	addHisto("epLKrR", 75, 0, 150);

	//78
	addHisto("emPMag", 60, 0, 60);
	addHisto("emPx", 60, -0.015, 0.015);
	addHisto("emPy", 60, -0.015, 0.015);
	addHisto("emPz", 110, 0, 1.1);
	addHisto("emEnergy", 60, 0, 60);
	addHisto("emeop", 100, 0, 1.5);
	addHisto("emLKrX", 150, -150, 150);
	addHisto("emLKrY", 150, -150, 150);
	addHisto("emLKrR", 75, 00, 150);

	//79
	addHisto("mee", 140, 0, 0.14);

	//pi+
	//85
	addHisto("pipPMag", 70, 0, 70);
	addHisto("pipPx", 60, -0.015, 0.015);
	addHisto("pipPy", 60, -0.015, 0.015);
	addHisto("pipPz", 110, 0, 1.1);
	addHisto("pipEnergy", 70, 0, 70);
	addHisto("pieop", 100, 0, 1.5);

	addHisto("t_epem_DCH", 75, 0, 150);
	addHisto("t_eppip_DCH", 100, 0, 200);
	addHisto("t_empip_DCH", 100, 0, 200);
	addHisto("t_epem_LKr", 125, 0, 250);
	addHisto("t_eppip_LKr", 125, 0, 250);
	addHisto("t_empip_LKr", 125, 0, 250);

	addHisto("t_gep_DCH", 100, 0, 200);
	addHisto("t_gem_DCH", 100, 0, 200);
	addHisto("t_gpip_DCH", 100, 0, 200);
	addHisto("t_gep_LKr", 125, 0, 250);
	addHisto("t_gem_LKr", 125, 0, 250);
	addHisto("t_gpip_LKr", 125, 0, 250);

	addHisto("undeft_gep_LKr", 100, 0, 200);
	addHisto("undeft_gem_LKr", 100, 0, 200);
	addHisto("undeft_gpip_LKr", 125, 0, 250);

	addHisto("L3_E_LKr_ep", 120, 0, 80);
	addHisto("L3_E_LKr_em", 120, 0, 80);
	addHisto("L3_E_LKr_gamma", 120, 0, 80);
	addHisto("L3_E_LKr", 120, 0, 80);

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
	// Rescale histo
	for (unsigned int i = 0; i < d1.size(); ++i) {
		SubSample::scale(d1.at(i), 1.);
	}
//	for (unsigned int i = 0; i < dMap.size(); ++i) {
//		SubSample::scale(dMap.at(i), 1.);
//	}
}

double CombineSample::getFFIntegral(double) {
	return 0;
}

void CombineSample::setPlotStyle(std::vector<int> color) {
	for (auto plot : d1)
		plot->SetFillColor(gStyle->GetColorPalette(color[0]));

//	for (auto plot : dMap)
//		plot->SetFillColor(gStyle->GetColorPalette(color[0]));
}

TH1D* CombineSample::getMainHisto() {
	return nullptr;
}

SubSample* CombineSample::Add(const SubSample* other) {
	SubSample::Add(other);
	const CombineSample *myOther = static_cast<const CombineSample*>(other);
	for (unsigned int i = 0; i < d1.size(); i++) {
		d1[i]->Add(myOther->d1[i], 1.);
	}
//	for (unsigned int i = 0; i < dMap.size(); i++) {
//		dMap[i]->Add(myOther->dMap[i], 1.);
//	}
	return this;
}
