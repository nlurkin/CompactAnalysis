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

// ### Histograms
const int stdNbins = 1000;
TH1D meeTrue = TH1D("meeTrue", "m_{#pi^{0}} (Good combination);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTrue = TH1D("mkTrue", "m_{K} (Good combination);m_{K}", stdNbins, 0.15, 0.8);
TH2D meeexTrue = TH2D("meeexTrue", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Good combination);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekTrue = TH2D("meekTrue", "m_{K} vs. m_{#pi^{0}} (Good combination);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeFalse = TH1D("meeFalse", "m_{#pi^{0}} (Wrong combination);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkFalse = TH1D("mkFalse", "m_{K} (Wrong combination);m_{K}", stdNbins, 0.15, 0.8);
TH2D meeexFalse = TH2D("meeexFalse", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Wrong combination);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekFalse = TH2D("meekFalse", "m_{K} vs. m_{#pi^{0}} (Wrong combination);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeTotal = TH1D("meeTotal", "m_{#pi^{0}} (Both combinations);m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTotal = TH1D("mkTotal", "m_{K} (Both combinations);m_{K}", stdNbins, 0.15, 0.8);
TH2D meeexTotal = TH2D("meeexTotal", "m_{e^{#pm}#x^{#mp}} vs. m_{#pi^{0}} (Both combinations);m_{#pi^{0}};m_{e^{#pm}#x^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekTotal = TH2D("meekTotal", "m_{K} vs. m_{#pi^{0}} (Both combinations);m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeDiffTrue = TH1D("meeDiffTrue", "m_{#pi^{0}}^{reco}-M_{#pi^{0}} (Good combination);m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
TH1D mkDiffTrue = TH1D("mkDiffTrue", "m_{K}^{reco}-M_{K} (Good combination);m_{K}^{reco}-M_{K}", stdNbins, -0.4, 0.4);

TH1D meeDiffFalse = TH1D("meeDiffFalse", "m_{#pi^{0}}^{reco}-M_{#pi^{0}} (Wrong combination);m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
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

TH2D combi2Dpi = TH2D("combi2Dpi", "combi2Dpi", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
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

	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1) Mk = abcog_params.mkp;
	else Mk = abcog_params.mkn;

	diffk1 = fabs(k1.M()-Mk);
	diffk2 = fabs(k2.M()-Mk);

	diffk1 = 0;
	diffk2 = 0;

	if(t==OptionsParser::K2PI){
		if( (diffpi01<io.cutsDefinition.k2pi.maxPi0MassDiff) && diffk1<io.cutsDefinition.k2pi.maxKaonMassDiff){
			nCandidates++;
			xCandidate = goodTrack2;
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

	meeTotal.Fill(ee1.M());
	meeTotal.Fill(ee2.M());
	mkTotal.Fill(k1.M());
	mkTotal.Fill(k2.M());
	meeexTotal.Fill(ee1.M(), (tem+t1x).M());
	meeexTotal.Fill(ee2.M(), (tem+t2x).M());
	meekTotal.Fill(ee1.M(), k1.M());
	meekTotal.Fill(ee2.M(), k2.M());

	combi2Dpi.Fill(ee1.M(), ee2.M());
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
			meeTrue.Fill(ee1.M());
			mkTrue.Fill(k1.M());
			meeexTrue.Fill(ee1.M(), (tem+t1x).M());
			meekTrue.Fill(ee1.M(), k1.M());

			meeDiffTrue.Fill(ee1.M()-Mpi0);
			mkDiffTrue.Fill(k1.M()-Mk);

			meeFalse.Fill(ee2.M());
			mkFalse.Fill(k2.M());
			meeexFalse.Fill(ee2.M(), (tem+t2x).M());
			meekFalse.Fill(ee2.M(), k2.M());

			meeDiffFalse.Fill(ee2.M()-Mpi0);
			mkDiffFalse.Fill(k2.M()-Mk);

			xTrue = pow((tem+t1ep).M()/Mpi0, 2.);
			xFalse = pow((tem+t2ep).M()/Mpi0, 2.);
		}
		else if(xPart==goodTrack1){
			meeTrue.Fill(ee2.M());
			mkTrue.Fill(k2.M());
			meeexTrue.Fill(ee2.M(), (tem+t2x).M());
			meekTrue.Fill(ee2.M(), k2.M());

			meeDiffTrue.Fill(ee2.M()-Mpi0);
			mkDiffTrue.Fill(k2.M()-Mk);

			meeFalse.Fill(ee1.M());
			mkFalse.Fill(k1.M());
			meeexFalse.Fill(ee1.M(), (tem+t1x).M());
			meekFalse.Fill(ee1.M(), k1.M());

			meeDiffFalse.Fill(ee1.M()-Mpi0);
			mkDiffFalse.Fill(k1.M()-Mk);

			xTrue = pow((tem+t2ep).M()/Mpi0,2.);
			xFalse = pow((tem+t1ep).M()/Mpi0, 2.);
		}
	}

	return nCandidates;
}

//Selection
int pi0d_tracksAcceptance(){
	bool lkrAcceptance;
	TVector3 propPos;
	bool badTrack = false;
	double radius;
	TVector3 dch1(rootGeom.Dch[0].PosChamber.x,rootGeom.Dch[0].PosChamber.y,rootGeom.Dch[0].PosChamber.z);
	TVector3 dch2(rootGeom.Dch[1].PosChamber.x,rootGeom.Dch[1].PosChamber.y,rootGeom.Dch[1].PosChamber.z);
	TVector3 dch4(rootGeom.Dch[3].PosChamber.x,rootGeom.Dch[3].PosChamber.y,rootGeom.Dch[3].PosChamber.z);

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); ++i){
		int iGoodTrack = corrEvent.goodTracks[i];
		NPhysicsTrack t = corrEvent.pTrack[iGoodTrack];

		propPos = propagateAfter(rootGeom.Lkr.z, t);
		lkrAcceptance = t.lkr_acc;
		if(options.isOptDebug()) cout << "LKr acceptance :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
		if(lkrAcceptance!=0) badTrack = true;

		// Track position on LKr with Pb Wall
		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\t\tPbWall y_LKr :\t\t-33.575 < " << propPos.Y() << " < -11.850: rejected" << endl;
			if(propPos.Y()>-33.575 && propPos.Y() < -11.850) badTrack = true;
		}

		propPos = propagateCorrBefore(rootGeom.Dch[0].PosChamber.z, t);
		radius = distance2D(dch1, propPos);
		if(options.isOptDebug()) cout << "DCH1 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateCorrBefore(rootGeom.Dch[1].PosChamber.z, t);
		radius = distance2D(dch2, propPos);
		if(options.isOptDebug()) cout << "DCH2 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;

		propPos = propagateAfter(rootGeom.Dch[3].PosChamber.z, t);
		radius = distance2D(dch4, propPos);
		if(options.isOptDebug()) cout << "DCH4 radius :\t\t" << radius << "\t <12 || > 110 : rejected" << endl;
		if(radius<12 || radius>110) badTrack = true;
	}

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
			propPos1 = propagateAfter(rootGeom.Lkr.z, t1);
			propPos2 = propagateAfter(rootGeom.Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
			if(RLKr<=20) bad = true;

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
			propPos1 = propagateAfter(rootGeom.Lkr.z, t1);
			propPos2 = propagateAfter(rootGeom.Lkr.z, t2);

			RLKr = distance2D(propPos1, propPos2);
			if(trackID1==xParticle.parentTrack || trackID2==xParticle.parentTrack){
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <50: rejected" << endl;
				if(RLKr<=50) bad = true;
			}
			else{
				if(options.isOptDebug()) cout << "\t\tR_LKr :\t\t" << RLKr << "\t <20: rejected" << endl;
				if(RLKr<=20) bad = true;
			}

			if(bad) badCombis++;
		}
	}

	return badCombis;
}

int pi0d_goodClusters_loose(){
	TVector3 propPos;
	double distance;
	double tDiff;

	int cond;
	int goodClusters = 0;

	int conditions;

	if(rootBurst.isData) conditions = 4;
	else conditions=3;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack t1 = corrEvent.pTrack[corrEvent.goodTracks[0]];
		NPhysicsTrack t2 = corrEvent.pTrack[corrEvent.goodTracks[1]];
		NPhysicsTrack t3 = corrEvent.pTrack[corrEvent.goodTracks[2]];

		// separation from x impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, t1);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_1 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, t2);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, t3);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_2 :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

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

	if(rootBurst.isData) conditions = 6;
	else conditions=5;

	if(options.isOptDebug()) cout << "\tNumber of vclusters :\t" << corrEvent.pCluster.size() << endl;
	if(options.isOptDebug()) cout << "\tNumber of clusters :\t" << rawEvent.Ncluster << endl;

	for(unsigned int i=0; i<corrEvent.pCluster.size(); i++){
		NPhysicsCluster c = corrEvent.pCluster[i];
		cond = 0;

		if(options.isOptDebug()) cout << "\tTrying cluster :\t" << i << endl;

		if(rootBurst.pbWall){
			if(options.isOptDebug()) cout << "\tPbWall distance y_cluster :\t-33.575 < " << c.position.Y() << " < -11.850 : reject" << endl;
			if(c.position.Y()>-33.575 && c.position.Y() < -11.850) continue;
		}
		if(options.isOptDebug()) cout << "\tEnergy :\t" << c.E << endl;

		NPhysicsTrack x = corrEvent.pTrack[xParticle.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[event.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[event.em.parentTrack];

		// separation from x impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, x);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_x :\t\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		// separation from e+ e- impact point >10cm
		propPos = propagateAfter(rootGeom.Lkr.z, ep);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e+ :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		propPos = propagateAfter(rootGeom.Lkr.z, em);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_e- :\t\t" << distance << "\t > 20 : ++" << endl;
		if(distance>20) cond++;

		// separation from undeflected e+ e- trajectories >20cm
		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

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
	meeTrue.Write();
	mkTrue.Write();
	meeexTrue.Write();
	meekTrue.Write();

	meeTotal.Write();
	mkTotal.Write();
	meeexTotal.Write();
	meekTotal.Write();

	meeDiffTrue.Write();
	mkDiffTrue.Write();

	meeFalse.Write();
	mkFalse.Write();
	meeexFalse.Write();
	meekFalse.Write();

	meeDiffFalse.Write();
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

	combi2Dpi.Write();
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


bool associateMCTracks(struct alt_pid_res &pid_res){
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
	}
	else{
		flBad = true;
		pid_res.bad[pid_res.currentID]++;
	}

	if((em1==-1 && em2==-1 && em3!=-1) || (em1==-1 && em2!=-1 && em3==-1) || (em1!=-1 && em2==-1 && em3==-1)){
		//OK
		if(em1!=-1) em=em1;
		if(em2!=-1) em=em2;
		if(em3!=-1) em=em3;
		pid_res.good[pid_res.currentID]++;
	}
	else{
		pid_res.bad[pid_res.currentID]++;
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
	}
	else{
		ep = -1;
		em = -1;
		xPart = -1;
	}

	return flBad;
}
