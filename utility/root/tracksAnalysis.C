#include <iostream>
#include <sstream>

double Mpic = 0.13957018;
double Me = 0.000510998928;
double Mpi0 = 0.1349766;
double MK = 0.493677;

using namespace std;

int pid_opposite_sign(TLorentzVector em, TLorentzVector ep, TLorentzVector gamma, TLorentzVector pip, int vtxCharge){
	TLorentzVector tem;
	TLorentzVector t1ep, t1x;
	TLorentzVector ee1;
	TLorentzVector k1;
	double Mx;
	int nCandidates = 0;
	Mx = Mpic;

	int goodTrack1, goodTrack2;

	int nNegative = 0;

	if(vtxCharge==1){
		t1x.SetVectM(em.Vect(), Mx);
		tem.SetVectM(pip.Vect(), Me);
		t1ep.SetVectM(ep.Vect(), Me);
	}
	else{
		t1x.SetVectM(ep.Vect(), Mx);
		t1ep.SetVectM(pip.Vect(), Me);
		tem.SetVectM(em.Vect(), Me);
	}

	ee1 = tem+t1ep+gamma;
	k1 = ee1+t1x;

	double Mk;
	double diffpi01;
	double diffk1;

	diffpi01 = fabs(ee1.M()-Mpi0);

	Mk = MK;

	diffk1 = fabs(k1.M()-Mk);

	double x1 = pow((tem+t1ep).M()/Mpi0,2.);

	double y1 = 2*ee1*(t1ep-tem)/(pow(Mpi0,2)*(1-x1));

	//cout << "Track1 pi0mass: " << ee1.M() << " kmass: " << k1.M() << endl;
	//cout << " M_pi0:" << ee1.M() << " <" << 0.115 << " || " << ee1.M() << " >" << 0.145 << " && " << endl;
	//cout << " M_K:" << k1.M() << " <" << 0.465 << " || " << k1.M() << " >" << 0.510 << " : rejected" << endl;
	//cout << " x:" << x1 << " <0 || " << x1 << " >1 : rejected" << endl;
	//cout << " y:" << y1 << " <0 || " << y1 << " >1 : rejected" << endl;
	if( ee1.M() > 0.115 && ee1.M() < 0.145
			&& k1.M() > 0.465 && k1.M() < 0.510
			&& x1<1 && x1>0
			&& fabs(y1)<1 && fabs(y1)>0){
		nCandidates++;
	}

	return nCandidates;
}

void tracksAnalysis(){

	gSystem->Load("../../obj/libexportClasses.so");
	gSystem->Load("../../obj/libfunLib.so");
	ROOTMCEvent *mc = new ROOTMCEvent();
	ROOTPhysicsEvent *corr = new ROOTPhysicsEvent();
	TFile *fd = TFile::Open("extra_pi.root");
	TTree *evt = (TTree*)fd->Get("event");

	evt->SetBranchAddress("mc", &mc);
	evt->SetBranchAddress("pi0dEvent", &corr);

	int nValidCandidates = 0;
	int nValid;
	int nGoodAssoc = 0;
	int nGoodSign=0;
	int nWrongSign=0;
	TLorentzVector realPip, realem, realep;

	double assocLimit = 3;
	for(int i=0; i<evt->GetEntries(); i++){
		evt->GetEntry(i);

		int charge = mc->k.pdgID<0 ? -1 : 1;
		int goodAssoc = 0;
		int goodem=0;
		int goodpi=0;
		int goodep=0;
		nValid = pid_opposite_sign(mc->em.P, mc->ep.P, mc->gamma.P, mc->pip.P, charge);

		if(nValid==0){
			if(charge==1){
				realPip = corr->em.P;
				realem = corr->pic.P;
				realep = corr->ep.P;
			}
			else{
				realPip = corr->ep.P;
				realem = corr->em.P;
				realep = corr->pic.P;
			}
			TVector3 diff = mc->em.P.Vect() - realem.Vect();
			if(diff.Mag()<assocLimit) goodem++;
			TVector3 diff = mc->ep.P.Vect() - realep.Vect();
			if(diff.Mag()<assocLimit) goodep++;
			TVector3 diff = mc->pip.P.Vect() - realPip.Vect();
			if(mc->pip.P.Vect().Mag()>0.001){
				if(diff.Mag()<assocLimit) goodpi++;
			}
			else goodpi++;

			goodAssoc = goodep + goodem + goodpi;

			if(goodAssoc!=3){
				if(goodAssoc==1){
				if(charge==1){
					if(goodem==0) nWrongSign++;
					else nGoodSign++;
				}
				else{
					if(goodep==0) nWrongSign++;
					else nGoodSign++;
				}
				}

				if(goodAssoc==2 && mc->pip.P.Vect().Mag()>0.01){
				TVector3 diff = mc->em.P.Vect() - realem.Vect();
				std::cout << (goodem==1 ? "" : "->") << "em (" << mc->em.P.X() << "," << mc->em.P.Y() << "," << mc->em.P.Z() << ") (" << realem.X() << "," << realem.Y() << "," << realem.Z() << ") (" << diff.X() << "," << diff.Y() << "," << diff.Z() << ") " << diff.Mag() << std::endl;
				TVector3 diff = mc->ep.P.Vect() - realep.Vect();
				std::cout << (goodep==1 ? "" : "->") << "ep (" << mc->ep.P.X() << "," << mc->ep.P.Y() << "," << mc->ep.P.Z() << ") (" << realep.X() << "," << realep.Y() << "," << realep.Z() << ") (" << diff.X() << "," << diff.Y() << "," << diff.Z() << ") " << diff.Mag()  << std::endl;
				TVector3 diff = mc->pip.P.Vect() - realPip.Vect();
				std::cout << (goodpi==1 ? "" : "->") << "pip(" << mc->pip.P.X() << "," << mc->pip.P.Y() << "," << mc->pip.P.Z() << ") (" << realPip.X() << "," << realPip.Y() << "," << realPip.Z() << ") (" << diff.X() << "," << diff.Y() << "," << diff.Z() << ") " << diff.Mag()  << std::endl;
				std::cout << mc->pip.P.Vect().Mag() << endl;
				TVector3 diff = mc->k.P.Vect() - corr->kaon.P.Vect();
				std::cout << "k " << mc->k.pdgID << " (" << mc->k.P.X() << "," << mc->k.P.Y() << "," << mc->k.P.Z() << ") (" << corr->kaon.P.X() << "," << corr->kaon.P.Y() << "," << corr->kaon.P.Z() << ") (" << diff.X() << "," << diff.Y() << "," << diff.Z() << ")" << std::endl;
				cout << endl << endl;
				}

			}
		}
		nValidCandidates += nValid;
		if(goodAssoc==3) nGoodAssoc++;
	}

	cout << "Valid candidates: " << nValidCandidates << endl;
	cout << "Good association: " << nGoodAssoc << endl;
	cout << "Good sign: " << nGoodSign << endl;
	cout << "Wrong sign: " << nWrongSign << endl;

}
