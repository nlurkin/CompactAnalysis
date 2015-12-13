/// Compact includes
#include "funLib.h"

#include "SelectionFunctions.h"
#include "pid_res.h"

CompactIO io;
OptionsParser options;

///### TTree objects
ROOTRawEvent &rawEvent = io.getRawEvent();
ROOTCorrectedEvent &corrEvent = io.getCorrEvent();
ROOTBurst &rootBurst = io.getRootBurst();
ROOTFileHeader &rootFileHeader = io.getInputFileHeader();
NGeom &rootGeom = io.getRootGeom();
ROOTMCEvent &rootMC = io.getRootMc();
ROOTPhysicsEvent &rootPhysics = io.getRootPhysics();
ROOTFileHeader &outputFileHeader = io.getOutputFileHeader();

// ### Database objects
NAbcog_params abcog_params;

// ### K2pi kmu3 selectors
//double Mx;
//int xPDGId;
//NRecoParticle *xRootPhysics;

const int stdNbins = 1000;
// ### Histograms for Michal pre_pid
TH1D meeTrue = TH1D("meeTrue", "m_{e^{+}e^{-}} (Good combination);m_{e^{+}e^{-}}", stdNbins, 0, 6);
TH1D m2pi0True = TH1D("m2pi0True", "m^{2}_{#pi^{0}} (Good combination);m^{2}_{#pi^{0}}", stdNbins, -1, 1);
TH2D mee2pi0True = TH2D("mee2pi0True", "m^{2}_{#pi^{0}} vs. m_{e^{+}e^{-}} (Good combination);m_{e^{+}e^{-}};m^{2}_{#pi^{0}}", stdNbins, 0, 6, stdNbins, -1, 1);

TH1D meeFalse = TH1D("meeFalse", "m_{e^{+}e^{-}} (Wrong combination);m_{e^{+}e^{-}}", stdNbins, 0, 6);
TH1D m2pi0False = TH1D("m2pi0False", "m_{#pi^{0}} (Wrong combination);m_{#pi^{0}}", stdNbins, -1, 1);
TH2D mee2pi0False = TH2D("mee2pi0False", "m^{2}_{#pi^{0}} vs. m_{e^{+}e^{-}} (Wrong combination);m_{e^{+}e^{-}};m^{2}_{#pi^{0}}", stdNbins, 0, 6, stdNbins, -1, 1);

TH1D meeTotal = TH1D("meeTotal", "m_{e^{+}e^{-}} (Both combinations);m_{e^{+}e^{-}}", stdNbins, 0, 6);
TH1D m2pi0Total = TH1D("m2pi0Total", "m_{#pi^{0}} (Both combinations);m_{#pi^{0}}", stdNbins, -1, 1);
TH2D mee2pi0Total = TH2D("mee2pi0Total", "m^{2}_{#pi^{0}} vs. m_{e^{+}e^{-}} (Both combinations);m_{e^{+}e^{-}};m^{2}_{#pi^{0}}", stdNbins, 0, 6, stdNbins, -1, 1);

TH1D m2pi0DiffTrue = TH1D("m2pi0DiffTrue", "m^{2}_{#pi^{0}}^{reco}-0.02 (Good combination);m^{2}_{#pi^{0}}^{reco}-0.02", stdNbins, -1, 1);

TH1D m2pi0DiffFalse = TH1D("m2pi0DiffFalse", "m^{2}_{#pi^{0}}^{reco}-0.02 (Wrong combination);m^{2}_{#pi^{0}}^{reco}-0.02", stdNbins, -1, 1);

// ### Histograms for pid
TH1D meegTrue = TH1D("meegTrue", "m_{#pi^{0}} (Good combination);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTrue = TH1D("mkTrue", "m_{K} (Good combination);m_{K}", stdNbins, 0.15, 0.8);
TH2D meegexTrue = TH2D("meegexTrue", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Good combination);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meegkTrue = TH2D("meegkTrue", "m_{K} vs. m_{#pi^{0}} (Good combination);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meegFalse = TH1D("meegFalse", "m_{#pi^{0}} (Wrong combination);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkFalse = TH1D("mkFalse", "m_{K} (Wrong combination);m_{K}", stdNbins, 0.15, 0.8);
TH2D meegexFalse = TH2D("meegexFalse", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Wrong combination);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meegkFalse = TH2D("meegkFalse", "m_{K} vs. m_{#pi^{0}} (Wrong combination);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meegTotal = TH1D("meegTotal", "m_{#pi^{0}} (Both combinations);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTotal = TH1D("mkTotal", "m_{K} (Both combinations);m_{K}", stdNbins, 0.15, 0.8);
TH2D meegexTotal = TH2D("meegexTotal", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Both combinations);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meegkTotal = TH2D("meegkTotal", "m_{K} vs. m_{#pi^{0}} (Both combinations);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meegDiffTrue = TH1D("meegDiffTrue", "m_{#pi^{0}}^{reco}-M_{#pi^{0}} (Good combination);m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
TH1D mkDiffTrue = TH1D("mkDiffTrue", "m_{K}^{reco}-M_{K} (Good combination);m_{K}^{reco}-M_{K}", stdNbins, -0.4, 0.4);

TH1D meegDiffFalse = TH1D("meegDiffFalse", "m_{#pi^{0}}^{reco}-M_{#pi^{0}} (Wrong combination);m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
TH1D mkDiffFalse = TH1D("mkDiffFalse", "m_{K}^{reco}-M_{K} (Wrong combination);m_{K}^{reco}-M_{K}", stdNbins, -0.4, 0.4);

TH1D nxCandidates = TH1D("nxCandidates", "nxCandidates", 20, 0, 20);
TH1D nxCandidatesNew = TH1D("nxCandidatesNew", "nxCandidatesNew", 20, 0, 20);

TH1D nBeamSign = TH1D("nBeamSign", "nBeamSign", 20, 0, 20);

TH2D xTruexFalse = TH2D("xTruexFalse", "xTruexFalse", 2*stdNbins, 0, 2, 2*stdNbins, 0, 2);
TH2D xTruexFalseMany = TH2D("xTruexFalseMany", "xTruexFalseMany", 2*stdNbins, 0, 2, 2*stdNbins, 0, 2);
TH2D xTruexFalseNo = TH2D("xTruexFalseNo", "xTruexFalseNo", 2*stdNbins, 0, 2, 2*stdNbins, 0, 2);
TH1D xMCNoID = TH1D("xMCNoID", "xMCNoID", 2*stdNbins, 0, 2);
TH1D xMCManyID = TH1D("xMCManyID", "xMCManyID", 2*stdNbins, 0, 2);

TH2D xTruexMCMany = TH2D("xTruexMCMany", "xTruexMCMany", 2*stdNbins, 0, 2, 2*stdNbins, 0, 2);
TH2D xTruexMCNo = TH2D("xTruexMCNo", "xTruexMCNo", 2*stdNbins, 0, 2, 2*stdNbins, 0, 2);

TH2D combi2Deeg = TH2D("combi2Deeg", "combi2Deeg", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D combi2Dk = TH2D("combi2Dk", "combi2Dk", stdNbins, 0.15, 0.8, stdNbins, 0.15, 0.8);
TH2D combi2Dk_exclu = TH2D("combi2Dk_exclu", "combi2Dk_exclu", stdNbins, 0.15, 0.8, stdNbins, 0.15, 0.8);

TH1D eopx = TH1D("eopx", "eopx", 1500, 0, 1.5);
TH1D eope = TH1D("eope", "eope", 1500, 0, 1.5);

TH1D eope_goodx = TH1D("eope_goodx", "eope_goodx", 1500, 0, 1.5);
TH1D eope_badx = TH1D("eope_badx", "eope_badx", 1500, 0, 1.5);
TH1D eopx_bade = TH1D("eopx_bade", "eopx_bade", 1500, 0, 1.5);
TH1D eopx_goode = TH1D("eopx_goode", "eopx_goode", 1500, 0, 1.5);

TH1D eoplowest = TH1D("eoplowest", "eoplowest", 1500, 0, 1.5);
TH1D eopsecond = TH1D("eopsecond", "eopsecond", 1500, 0, 1.5);
TH1D eophighest = TH1D("eophighest", "eophighest", 1500, 0, 1.5);

int em, ep, xPart;
bool flBad;
double xTrue, xFalse;
//bool isK2PiType = true;
//bool isKMu3Type = true;

int michal_prepid(int &xCandidate, OptionsParser::ESelectionType t){
	TLorentzVector tem;
	TLorentzVector t1ep, t1pi;
	TLorentzVector t2ep, t2pi;
	TLorentzVector pi01, pi02;
	double x1, x2;

	int goodTrack1, goodTrack2;

	int nCandidates = 0;

	int nNegative = 0;
	int vtxCharge = rawEvent.vtx[corrEvent.goodVertexID].charge;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].q==-1*vtxCharge){
		goodTrack1 = 1;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 1;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p, Me);
		nNegative++;
	}

	nBeamSign.Fill(nNegative);
	if(nNegative!=1) return 0;

	TLorentzVector myKaon;
	myKaon.SetVectM(corrEvent.kaonMomentum*corrEvent.kaonP, MK);

	//Try e = goodTrack1, pi = goodTrack2
	t1ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Me);
	t1pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Mpic);
	x1 = (tem+t1ep).M2()/pow(Mpi0,2);
	pi01 = myKaon - t1pi;

	//Try e = goodTrack2, pi = goodTrack1
	t2pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Mpic);
	t2ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Me);
	x2 = (tem+t2ep).M2()/pow(Mpi0,2);
	pi02 = myKaon - t2pi;

	double diffpi01, diffpi02;

	if(t==OptionsParser::K2PI){
		diffpi01 = fabs(pi01.M2() - 0.02);
		diffpi02 = fabs(pi02.M2() - 0.02);
	}
	else if(t==OptionsParser::KMU3){
		//Bypass the missing mass cut
		diffpi01 = 0;
		diffpi02 = 0;
	}

	if( (x1<1) && (diffpi01<0.01) ){
		nCandidates++;
		xCandidate = goodTrack2;
	}

	if( (x2<1) && (diffpi02<0.01) ){
		nCandidates++;
		xCandidate = goodTrack1;
	}

	meeTotal.Fill(x1);
	meeTotal.Fill(x2);
	m2pi0Total.Fill(pi01.M2());
	m2pi0Total.Fill(pi02.M2());
	mee2pi0Total.Fill(x1, pi01.M2());
	mee2pi0Total.Fill(x2, pi02.M2());

	if(!flBad){
		//Good MC association, fill the plots
		if(xPart==goodTrack2){
			meeTrue.Fill(x1);
			m2pi0True.Fill(pi01.M2());
			mee2pi0True.Fill(x1, pi01.M2());

			m2pi0DiffTrue.Fill(pi01.M2()-0.02);

			meeFalse.Fill(x2);
			m2pi0False.Fill(pi02.M2());
			mee2pi0False.Fill(x2, pi02.M2());

			m2pi0DiffFalse.Fill(pi02.M2()-0.02);
		}
		else if(xPart==goodTrack1){
			meeTrue.Fill(x2);
			m2pi0True.Fill(pi02.M2());
			mee2pi0True.Fill(x2, pi02.M2());

			m2pi0DiffTrue.Fill(pi02.M2()-0.02);

			meeFalse.Fill(x1);
			m2pi0False.Fill(pi01.M2());
			mee2pi0False.Fill(x1, pi01.M2());

			m2pi0DiffFalse.Fill(pi01.M2()-0.02);
		}
	}

	return nCandidates;
}

int pid(int &xCandidate, TLorentzVector &gamma, OptionsParser::ESelectionType t){
	TLorentzVector tem;
	TLorentzVector t1ep, t1x;
	TLorentzVector t2ep, t2x;
	TLorentzVector ee1, ee2;
	TLorentzVector k1, k2;
	double Mx;
	int nCandidates = 0;
	if(t==OptionsParser::K2PI) Mx = Mpic;
	else if(t==OptionsParser::KMU3) Mx = Mmu;

	int goodTrack1, goodTrack2;

	int nNegative = 0;
	int vtxCharge = rawEvent.vtx[corrEvent.goodVertexID].charge;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].q==-1*vtxCharge){
		goodTrack1 = 1;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 2;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p, Me);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 1;
		tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p, Me);
		nNegative++;
	}

	nBeamSign.Fill(nNegative);
	if(nNegative!=1) return 0;

	//Try e = goodTrack1, x = goodTrack2
	t1ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Me);
	t1x.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Mx);
	ee1 = tem+t1ep+gamma;
	k1 = ee1+t1x;

	//Try e = goodTrack2, x = goodTrack1
	t2x.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Mx);
	t2ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Me);
	ee2 = tem+t2ep+gamma;
	k2 = ee2+t2x;

	double Mk;
	double diffpi01, diffpi02;
	double diffk1, diffk2;

	diffpi01 = fabs(ee1.M()-Mpi0);
	diffpi02 = fabs(ee2.M()-Mpi0);

	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1) Mk = MK;
	//else Mk = abcog_params.mkn;
	else Mk = MK;

	diffk1 = fabs(k1.M()-Mk);
	diffk2 = fabs(k2.M()-Mk);

	//diffk1 = 0;
	//diffk2 = 0;

	if(t==OptionsParser::K2PI){
		if(options.isOptDebug()){
			cout << "Track1 pi0mass: " << ee1.M() << " kmass: " << k1.M() << endl;
			cout << " diffpi0:" << diffpi01 << " >" << io.cutsDefinition.k2pi.maxPi0MassDiff << " && " << endl;
			cout << " diffk:" << diffk1 << " >" << io.cutsDefinition.k2pi.maxKaonMassDiff << " : rejected" << endl;
		}
		if( (diffpi01<io.cutsDefinition.k2pi.maxPi0MassDiff) && diffk1<io.cutsDefinition.k2pi.maxKaonMassDiff){
			nCandidates++;
			xCandidate = goodTrack2;
		}
		if(options.isOptDebug()){
			cout << "Track2 pi0mass: " << ee2.M() << " kmass: " << k2.M() << endl;
			cout << " diffpi0:" << diffpi02 << " >" << io.cutsDefinition.k2pi.maxPi0MassDiff << " && " << endl;
			cout << " diffk:" << diffk2 << " >" << io.cutsDefinition.k2pi.maxKaonMassDiff << " : rejected" << endl;
		}
		if( (diffpi02<io.cutsDefinition.k2pi.maxPi0MassDiff) && diffk2<io.cutsDefinition.k2pi.maxKaonMassDiff){
			nCandidates++;
			xCandidate = goodTrack1;
		}
	}
	else if(t==OptionsParser::KMU3){
		if( (diffpi01<io.cutsDefinition.kmu3.maxPi0MassDiff)){
			nCandidates++;
			xCandidate= goodTrack2;
		}
		if( (diffpi02<io.cutsDefinition.kmu3.maxPi0MassDiff)){
			nCandidates++;
			xCandidate= goodTrack1;
		}
	}

	meegTotal.Fill(ee1.M());
	meegTotal.Fill(ee2.M());
	mkTotal.Fill(k1.M());
	mkTotal.Fill(k2.M());
	meegexTotal.Fill(ee1.M(), (tem+t1x).M());
	meegexTotal.Fill(ee2.M(), (tem+t2x).M());
	meegkTotal.Fill(ee1.M(), k1.M());
	meegkTotal.Fill(ee2.M(), k2.M());

	combi2Deeg.Fill(ee1.M(), ee2.M());
	combi2Dk.Fill(k1.M(), k2.M());

	if(t==OptionsParser::K2PI){
		if(diffpi01<io.cutsDefinition.k2pi.maxPi0MassDiff && diffpi02<io.cutsDefinition.k2pi.maxPi0MassDiff){
			combi2Dk_exclu.Fill(k1.M(), k2.M());
		}
	}
	else if(t==OptionsParser::KMU3){
		if(diffpi01<io.cutsDefinition.kmu3.maxPi0MassDiff){
			combi2Dk_exclu.Fill(k1.M(), k2.M());
		}
	}

	if(!flBad){
		//Good MC association, fill the plots
		if(xPart==goodTrack2){
			meegTrue.Fill(ee1.M());
			mkTrue.Fill(k1.M());
			meegexTrue.Fill(ee1.M(), (tem+t1x).M());
			meegkTrue.Fill(ee1.M(), k1.M());

			meegDiffTrue.Fill(ee1.M()-Mpi0);
			mkDiffTrue.Fill(k1.M()-Mk);

			meegFalse.Fill(ee2.M());
			mkFalse.Fill(k2.M());
			meegexFalse.Fill(ee2.M(), (tem+t2x).M());
			meegkFalse.Fill(ee2.M(), k2.M());

			meegDiffFalse.Fill(ee2.M()-Mpi0);
			mkDiffFalse.Fill(k2.M()-Mk);

			xTrue = pow((tem+t1ep).M()/Mpi0, 2.);
			xFalse = pow((tem+t2ep).M()/Mpi0, 2.);
		}
		else if(xPart==goodTrack1){
			meegTrue.Fill(ee2.M());
			mkTrue.Fill(k2.M());
			meegexTrue.Fill(ee2.M(), (tem+t2x).M());
			meegkTrue.Fill(ee2.M(), k2.M());

			meegDiffTrue.Fill(ee2.M()-Mpi0);
			mkDiffTrue.Fill(k2.M()-Mk);

			meegFalse.Fill(ee1.M());
			mkFalse.Fill(k1.M());
			meegexFalse.Fill(ee1.M(), (tem+t1x).M());
			meegkFalse.Fill(ee1.M(), k1.M());

			meegDiffFalse.Fill(ee1.M()-Mpi0);
			mkDiffFalse.Fill(k1.M()-Mk);

			xTrue = pow((tem+t2ep).M()/Mpi0,2.);
			xFalse = pow((tem+t1ep).M()/Mpi0, 2.);
		}
	}

	return nCandidates;
}

int pid_opposite_sign(int &xCandidate, TLorentzVector &gamma, OptionsParser::ESelectionType t){
	TLorentzVector tem;
	TLorentzVector t1ep, t1x;
	TLorentzVector ee1;
	TLorentzVector k1;
	double Mx;
	int nCandidates = 0;
	if(t==OptionsParser::K2PI) Mx = Mpic;
	else if(t==OptionsParser::KMU3) Mx = Mmu;

	int goodTrack1, goodTrack2;

	int nNegative = 0;
	int vtxCharge = rawEvent.vtx[corrEvent.goodVertexID].charge;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].q==-1*vtxCharge){
		goodTrack1 = 1;
		goodTrack2 = 2;
		xCandidate = 0;
		t1x.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p, Mpic);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 2;
		xCandidate = 1;
		t1x.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p, Mpic);
		nNegative++;
	}
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].q==-1*vtxCharge){
		goodTrack1 = 0;
		goodTrack2 = 1;
		xCandidate = 2;
		t1x.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p, Mpic);
		nNegative++;
	}

	if(nNegative!=1) return 0;

	//Try e = goodTrack1, x = goodTrack2
	t1ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Me);
	tem.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Me);
	ee1 = tem+t1ep+gamma;
	k1 = ee1+t1x;

	double Mk;
	double diffpi01;
	double diffk1;

	diffpi01 = fabs(ee1.M()-Mpi0);

	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1) Mk = MK;
	else Mk = MK;

	diffk1 = fabs(k1.M()-Mk);


	if(t==OptionsParser::K2PI){
		if(options.isOptDebug()){
			cout << "Track1 pi0mass: " << ee1.M() << " kmass: " << k1.M() << endl;
			cout << " diffpi0:" << diffpi01 << " >" << io.cutsDefinition.k2pi.maxPi0MassDiff << " && " << endl;
			cout << " diffk:" << diffk1 << " >" << io.cutsDefinition.k2pi.maxKaonMassDiff << " : rejected" << endl;
		}
		if( (diffpi01<io.cutsDefinition.k2pi.maxPi0MassDiff) && diffk1<io.cutsDefinition.k2pi.maxKaonMassDiff){
			nCandidates++;
		}
	}
	else if(t==OptionsParser::KMU3){
		if( (diffpi01<io.cutsDefinition.kmu3.maxPi0MassDiff)){
			nCandidates++;
		}
	}

	meegTotal.Fill(ee1.M());
	mkTotal.Fill(k1.M());
	meegexTotal.Fill(ee1.M(), (tem+t1x).M());
	meegkTotal.Fill(ee1.M(), k1.M());

	if(!flBad){
		//Good MC association, fill the plots
		if(xPart==goodTrack2){
			meegTrue.Fill(ee1.M());
			mkTrue.Fill(k1.M());
			meegexTrue.Fill(ee1.M(), (tem+t1x).M());
			meegkTrue.Fill(ee1.M(), k1.M());

			meegDiffTrue.Fill(ee1.M()-Mpi0);
			mkDiffTrue.Fill(k1.M()-Mk);

			xTrue = pow((tem+t1ep).M()/Mpi0, 2.);
		}
		else if(xPart==goodTrack1){
			meegFalse.Fill(ee1.M());
			mkFalse.Fill(k1.M());
			meegexFalse.Fill(ee1.M(), (tem+t1x).M());
			meegkFalse.Fill(ee1.M(), k1.M());

			meegDiffFalse.Fill(ee1.M()-Mpi0);
			mkDiffFalse.Fill(k1.M()-Mk);

			xFalse = pow((tem+t1ep).M()/Mpi0, 2.);
		}
	}
	return nCandidates;
}

bool pi0d_L3Trigger(NPhysicsTrack &t){
	TVector3 propPos = propagateAfter(rootGeom.Lkr.z, t);
	bool lkrAcceptance = t.lkr_acc;
	bool goodAcceptance = false;
	if(options.isOptDebug()) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t == 0 && " << t.p << " >=6 && " <<  t.E/t.p << " >0.8: ok" << endl;
	if(lkrAcceptance==0 && t.p>=6 && t.E/t.p>0.8 && rawEvent.track[t.trackID].dDeadCell>2) goodAcceptance=true;

	// Track position on LKr with Pb Wall
	bool goodPBWall = true;
	if(rootBurst.pbWall){
		if(options.isOptDebug()) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
		if(propPos.Y()>-33.575 && propPos.Y() < -11.850) goodPBWall = false;
	}

	if(goodAcceptance && goodPBWall) return true;
	return false;
}

//Selection
int pi0d_tracksAcceptance(){
	//bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(rootGeom.Dch[0].PosChamber.x,rootGeom.Dch[0].PosChamber.y,rootGeom.Dch[0].PosChamber.z);
	TVector3 dch2(rootGeom.Dch[1].PosChamber.x,rootGeom.Dch[1].PosChamber.y,rootGeom.Dch[1].PosChamber.z);
	TVector3 dch4(rootGeom.Dch[3].PosChamber.x,rootGeom.Dch[3].PosChamber.y,rootGeom.Dch[3].PosChamber.z);

	int ntrackLkr = 0;
	for(unsigned int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		//propPos = propagateAfter(rootGeom.Lkr.z, t);
		//lkrAcceptance = t.lkr_acc;
		//bool goodAcceptance = false;
		//if(options.isOptDebug()) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t == 0 && " << t.p << " >=6 && " <<  t.E/t.p << " >0.8: ok" << endl;
		//to remove
		//if(lkrAcceptance!=0) badTrack = true;
		//if(lkrAcceptance==0 && t.p>=6 && t.E/t.p>0.8 && rawEvent.track[t.trackID].dDeadCell>2) goodAcceptance=true;
		//goodAcceptance = true;

		// Track position on LKr with Pb Wall
		//bool goodPBWall = true;
		//if(rootBurst.pbWall){
		//	if(options.isOptDebug()) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			//to remove
			//if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		//	if(propPos.Y()>-33.575 && propPos.Y() < -11.850) goodPBWall = false;
		//}

		//if(goodAcceptance && goodPBWall) ntrackLkr++;

		propPos = propagateBefore(rootGeom.Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		if(options.isOptDebug()) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateBefore(rootGeom.Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		if(options.isOptDebug()) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(rootGeom.Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		if(options.isOptDebug()) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;
	}

	//At least 1 e+/e- track in lkr acceptance
	//if(ntrackLkr==0) badTrack = true;
	return badTrack;
}

int pi0d_trackCombinationVeto_loose(){
	int ntracks = 3;

	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			int trackID1 = corrEvent.goodTracks[i];
			int trackID2 = corrEvent.goodTracks[j];
			NPhysicsTrack t1 = corrEvent.pTrack[trackID1];
			NPhysicsTrack t2 = corrEvent.pTrack[trackID2];
			bad = false;
			if(options.isOptDebug()) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(options.isOptDebug()) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(rootBurst.isData){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				if(options.isOptDebug()) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			/*propPos1 = propagateAfter(rootGeom.Lkr.z, t1);
			propPos2 = propagateAfter(rootGeom.Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
			//to remove
			if(RLKr<=20) bad = true;*/

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_trackCombinationVeto_tight(NRecoParticle &xParticle){
	int ntracks = 3;
	TVector3 propPos1, propPos2;
	double RDCH1, RLKr;
	double tDiff;

	bool bad = false;

	int badCombis = 0;

	for(int i=0; i<ntracks-1; i++){
		for(int j=i+1; j<ntracks; j++){
			int trackID1 = corrEvent.goodTracks[i];
			int trackID2 = corrEvent.goodTracks[j];
			NPhysicsTrack t1 = corrEvent.pTrack[trackID1];
			NPhysicsTrack t2 = corrEvent.pTrack[trackID2];
			bad = false;
			if(options.isOptDebug()) cout << "\tTrying combination :\t" << i << " " << j << endl;

			// Track-to-Track distance in DCH1 plane >1cm
			propPos1 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t1);
			propPos2 = propagateBefore(rootGeom.Dch[0].PosChamber.z, t2);

			RDCH1 = distance2D(propPos1, propPos2);
			if(options.isOptDebug()) cout << "\t\tR_DCH1 :\t" << RDCH1 << "\t <2: rejected" << endl;
			if(RDCH1<=2) bad = true;

			// Track DCH times
			// |t1-t2|<15
			if(rootBurst.isData){
				tDiff = fabs(rawEvent.track[t1.trackID].time - rawEvent.track[t2.trackID].time);
				if(options.isOptDebug()) cout << "\t\t|t_i-t_j| :\t" << tDiff << "\t >15: rejected" << endl;
				if(tDiff>=15) bad = true;
			}

			// Track separation in LKr plane >15
			/*propPos1 = propagateAfter(rootGeom.Lkr.z, t1);
			propPos2 = propagateAfter(rootGeom.Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(trackID1==xParticle.parentTrack || trackID2==xParticle.parentTrack){
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
				if(RLKr<=20) bad = true;
			}*/

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_goodClusters_loose(){
	TVector3 propPos;
	double distance;
	double trackR;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) 	conditions = 7;
	else 					conditions = 6;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		//Ignore clusters behind Pb Wall
		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		//Ignore clusters with low energy
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;
		if(c.E < 3.) continue;

		NPhysicsTrack t1 = corrEvent.pTrack[corrEvent.goodTracks[0]];
		NPhysicsTrack t2 = corrEvent.pTrack[corrEvent.goodTracks[1]];
		NPhysicsTrack t3 = corrEvent.pTrack[corrEvent.goodTracks[2]];

		// separation from x impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, t1);
		distance = distance2D(propPos, c.position);
		trackR = distance2D(propPos, TVector3(0,0,0));
		if(options.isOptDebug()) cout << "\t\td_LKr_1 :\t\t" << distance << "\t > 20 || R_LKr_1 :\t" <<  trackR << "<10 : ++" << endl;
		if(distance>20 || trackR<10) cond++;

		// separation from undeflected x trajectories >20cm
		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[t1.trackID].bDetPos, rawEvent.track[t1.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_1 :\t" << distance << "\t > " << io.cutsDefinition.unDeflectedElDist << " : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, t2);
		distance = distance2D(propPos, c.position);
		trackR = distance2D(propPos, TVector3(0,0,0));
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 || R_LKr_1 :\t" <<  trackR << "<10 : ++" << endl;
		if(distance>20 || trackR<10) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[t2.trackID].bDetPos, rawEvent.track[t2.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_2 :\t" << distance << "\t > " << io.cutsDefinition.unDeflectedElDist << " : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, t3);
		distance = distance2D(propPos, c.position);
		trackR = distance2D(propPos, TVector3(0,0,0));
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 || R_LKr_1 :\t" <<  trackR << "<10 : ++" << endl;
		if(distance>20 || trackR<10) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[t3.trackID].bDetPos, rawEvent.track[t3.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_3 :\t" << distance << "\t > " << io.cutsDefinition.unDeflectedElDist << " : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(options.isOptDebug()) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(options.isOptDebug()) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			rootPhysics.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

int pi0d_goodClusters_tight(NRecoParticle &xParticle, ROOTPhysicsEvent &event){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) conditions = 7;
	else conditions=6;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		//Ignore clusters behind Pb Wall
		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		//Ignore clusters with low energy
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;
		if(c.E < 3.) continue;

		NPhysicsTrack x = corrEvent.pTrack[xParticle.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[event.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[event.em.parentTrack];

		// separation from x impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, x);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_x :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>20 || distance2D(propPos, TVector3(0,0,0))<10) cond++;
		//cond++;

		// separation from undeflected x impact point >20cm
		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[x.trackID].bDetPos, rawEvent.track[x.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_x :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;
		//cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, ep);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20 || distance2D(propPos, TVector3(0,0,0))<10) cond++;
		//cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, em);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20 || distance2D(propPos, TVector3(0,0,0))<10) cond++;
		//cond++;

		// separation from undeflected e+ e- trajectories >20cm
		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;
		//cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>io.cutsDefinition.unDeflectedElDist) cond++;
		//cond++;

		// |t_g - t_vtx|<10ns
		if(rootBurst.isData){
			tDiff = fabs(rawEvent.cluster[c.clusterID].time - rawEvent.vtx[corrEvent.goodVertexID].time);
			if(options.isOptDebug()) cout << "\t\t|t_g - t_vtx|:\t\t" << tDiff << "\t < 10 : ++" << endl;
			if(tDiff<10) cond++;
		}

		if(options.isOptDebug()) cout << "\tConditions :\t\t" << cond << "\t == " << conditions << " : Good cluster" << endl;
		if(cond==conditions){
			goodClusters++;
			event.gamma.parentCluster = i;
		}
	}
	return goodClusters;
}

bool pi0d_failCut(int i){
	if(!options.isDoScan())corrEvent.failedCond = i;
	if(options.isOptDebug()) cout << "Event is not passing selection " << i << endl;
	if(io.isDoOutput()) io.output.f1() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << " " << i << endl;
	return false;
}
void pi0d_passSelection(){
	if(options.isOptDebug()) cout << "Event is passing selection" << endl;
	if(io.isDoOutput()) io.output.f2() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << endl;
}

void savePlots(){
	//Michal ones
	meeTrue.Write();
	m2pi0True.Write();
	mee2pi0True.Write();

	meeTotal.Write();
	m2pi0Total.Write();
	mee2pi0Total.Write();

	m2pi0DiffTrue.Write();

	meeFalse.Write();
	m2pi0False.Write();
	mee2pi0False.Write();

	m2pi0DiffFalse.Write();

	//PID ones
	meegTrue.Write();
	mkTrue.Write();
	meegexTrue.Write();
	meegkTrue.Write();

	meegTotal.Write();
	mkTotal.Write();
	meegexTotal.Write();
	meegkTotal.Write();

	meegDiffTrue.Write();
	mkDiffTrue.Write();

	meegFalse.Write();
	mkFalse.Write();
	meegexFalse.Write();
	meegkFalse.Write();

	meegDiffFalse.Write();
	mkDiffFalse.Write();

	nxCandidates.Write();
	nxCandidatesNew.Write();
	nBeamSign.Write();
	xTruexFalse.Write();
	xTruexFalseMany.Write();
	xTruexFalseNo.Write();
	xMCNoID.Write();
	xMCManyID.Write();

	xTruexMCMany.Write();
	xTruexMCNo.Write();

	combi2Deeg.Write();
	combi2Dk.Write();
	combi2Dk_exclu.Write();

	eopx.Write();
	eope.Write();

	eope_goodx.Write();
	eope_badx.Write();
	eopx_bade.Write();
	eopx_goode.Write();

	eoplowest.Write();
	eopsecond.Write();
	eophighest.Write();
}

bool pi0d_failCutInc(int i, bool assoc, bool good, bool bad, struct alt_pid_res *pid_res){
	if(pid_res!=NULL) pid_res->incNonPass(assoc, good, bad);
	if(!options.isDoScan())corrEvent.failedCond = i;
	if(options.isOptDebug()) cout << "Event is not passing selection " << i << endl;
	if(io.isDoOutput()) io.output.f1() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << " " << i << endl;
	return false;
}

vector<float> sortEOP(){
	double eop;
	vector<float> ret;

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		int goodTrackID = corrEvent.goodTracks[i];
		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;
		ret.push_back(eop);
	}

	sort(ret.begin(), ret.end());

	return ret;
}

int pi0d_identifyPiEOP(int &piCandidate, bool &badElectron){
	double eop;

	int piCandidatesNb = 0;
	badElectron = false;

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		int goodTrackID = corrEvent.goodTracks[i];
		if(options.isOptDebug()) cout << "\tTrying track :\t\t" << i <<  endl;

		eop = corrEvent.pTrack[goodTrackID].E/corrEvent.pTrack[goodTrackID].p;

		if(options.isOptDebug()) cout << "\tE over P :\t\t" << eop << "\t\t <0.85 : pi+ candidate; >1.15 : badElectron" <<  endl;
		if(eop<0.85){
			piCandidatesNb++;
			piCandidate = i;
		}
		else if(eop>1.15) badElectron = true;

	}

	return piCandidatesNb;
}


bool associateMCTracks(struct alt_pid_res &pid_res, struct alt_pid_res *pid_res_mu){
	double limit = 2e-4;
	bool flBad = false;

	int ep1=-1, ep2=-1, ep3=-1;
	int em1=-1, em2=-1, em3=-1;
	//Find out the real tracks
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()<limit) ep1=0;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()<limit) ep2=1;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()<limit) ep3=2;

	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()<limit) em1=0;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()<limit) em2=1;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()<limit) em3=2;

	if((ep1==-1 && ep2==-1 && ep3!=-1) || (ep1==-1 && ep2!=-1 && ep3==-1) || (ep1!=-1 && ep2==-1 && ep3==-1)){
		//OK
		if(ep1!=-1) ep=ep1;
		if(ep2!=-1) ep=ep2;
		if(ep3!=-1) ep=ep3;
		pid_res.good[pid_res.currentID]++;
		if(pid_res_mu) pid_res_mu->good[pid_res_mu->currentID]++;
	}
	else{
		flBad = true;
		pid_res.bad[pid_res.currentID]++;
		if(pid_res_mu) pid_res_mu->bad[pid_res_mu->currentID]++;
	}

	if((em1==-1 && em2==-1 && em3!=-1) || (em1==-1 && em2!=-1 && em3==-1) || (em1!=-1 && em2==-1 && em3==-1)){
		//OK
		if(em1!=-1) em=em1;
		if(em2!=-1) em=em2;
		if(em3!=-1) em=em3;
		pid_res.good[pid_res.currentID]++;
		if(pid_res_mu) pid_res_mu->good[pid_res_mu->currentID]++;
	}
	else{
		pid_res.bad[pid_res.currentID]++;
		if(pid_res_mu) pid_res_mu->bad[pid_res_mu->currentID]++;
		flBad = true;
	}
	//if(em==ep){ // Both are identified to the same tracks
		//pid_res.good -= 2;
		//pid_res.bad += 2;
		//flBad = true;
	//}

	if(!flBad){
		if((ep==0 && em==1) || (ep==1 && em==0)) xPart=2;
		else if((ep==0 && em==2) || (ep==2 && em==0)) xPart=1;
		else if((ep==1 && em==2) || (ep==2 && em==1)) xPart=0;
		pid_res.incAssociated();
		if(pid_res_mu) pid_res_mu->incAssociated();
	}
	else{
		ep = -1;
		em = -1;
		xPart = -1;
	}

	return flBad;
}
