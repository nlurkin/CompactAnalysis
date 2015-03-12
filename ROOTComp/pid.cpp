#define OUTSIDECOMPACT

/// Compact includes
#include "funLib.h"

/// Std include
#include <iomanip>
using namespace std;

// Local includes
#include "CompactIO.h"
#include "OptionsParser.h"

#include <TH2D.h>
#include <TTree.h>

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

const int stdNbins = 1000;

TH1D meeTrue = TH1D("meeTrue", "meegTrue;m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTrue = TH1D("mkTrue", "mkTrue;m_{K}", stdNbins, 0.15, 0.8);
TH2D meeepiTrue = TH2D("meeepiTrue", "meegepiTrue;m_{#pi^{0}};m_{e^{#pm}#pi^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekTrue = TH2D("meekTrue", "meegkTrue;m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeFalse = TH1D("meeFalse", "meegFalse;m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkFalse = TH1D("mkFalse", "mkFalse;m_{K}", stdNbins, 0.15, 0.8);
TH2D meeepiFalse = TH2D("meeepiFalse", "meegepiFalse;m_{#pi^{0}};m_{e^{#pm}#pi^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekFalse = TH2D("meekFalse", "meegkFalse;m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeTotal = TH1D("meeTotal", "meegTotal;m_{#pi^{0}}", stdNbins, 0, 0.6);
TH1D mkTotal = TH1D("mkTotal", "mkTotal;m_{K}", stdNbins, 0.15, 0.8);
TH2D meeepiTotal = TH2D("meeepiTotal", "meegepiTotal;m_{#pi^{0}};m_{e^{#pm}#pi^{#mp}}", stdNbins, 0, 0.6, stdNbins, 0, 0.6);
TH2D meekTotal = TH2D("meekTotal", "meegkTotal;m_{#pi^{0}};m_{K}", stdNbins, 0, 0.6, stdNbins, 0, 1);

TH1D meeDiffTrue = TH1D("meeDiffTrue", "meegDiffTrue;m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
TH1D mkDiffTrue = TH1D("mkDiffTrue", "mkDiffTrue;m_{K}^{reco}-M_{K}", stdNbins, -0.4, 0.4);

TH1D meeDiffFalse = TH1D("meeDiffFalse", "meegDiffFalse;m_{#pi^{0}}^{reco}-M_{#pi^{0}}", stdNbins, -0.2, 0.5);
TH1D mkDiffFalse = TH1D("mkDiffFalse", "mkDiffFalse;m_{K}^{reco}-M_{K}", stdNbins, -0.4, 0.4);


TH1D nPiCandidates = TH1D("nPiCandidates", "nPiCandidates", 20, 0, 20);
TH1D nPiCandidatesNew = TH1D("nPiCandidatesNew", "nPiCandidatesNew", 20, 0, 20);

TH1D nBeamSign = TH1D("nBeamSign", "nBeamSign", 20, 0, 20);

TH2D xTruexFalse = TH2D("xTruexFalse", "xTruexFalse", 2*stdNbins, 0, 2, 2*stdNbins, 0, 20);

struct alt_pid_res{
	void incTotal(){ total.events++;}
	void incPrelim(){ total.prelim.events++;}
	void incNoID(bool assoc){
		if(assoc) total.prelim.associated.noID.events++;
		else total.prelim.noID.events++;
	}
	void incManyID(bool assoc){
		if(assoc) total.prelim.associated.manyID.events++;
		else total.prelim.manyID.events++;
	}
	void incAssociated(){
		total.prelim.associated.events++;
	}
	void incIded(bool assoc){
		if(assoc) total.prelim.associated.ided.events++;
		else total.prelim.ided.events++;
	}
	void incGoodId(){
		total.prelim.associated.ided.good.events++;
	}
	void incBadId(){
		total.prelim.associated.ided.bad.events++;
	}
	void incPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.pass.events++;
			else if(bad) total.prelim.associated.ided.bad.pass.events++;
		}
		else{
			total.prelim.ided.pass.events++;
		}
	}
	void incNonPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.noPass.events++;
			else if(bad) total.prelim.associated.ided.bad.noPass.events++;
		}
		else{
			total.prelim.ided.noPass.events++;
		}
	}

	//MC association
	int good=0, bad=0;

	struct total_t{
		int events = 0;

		struct prelim_t{
			int events = 0;
			struct noID_t{
				int events = 0;
			} noID;
			struct manyID_t{
				int events = 0;
			} manyID;

			struct associated_t{
				int events = 0;
				struct noID_t{
					int events = 0;
				} noID;
				struct manyID_t{
					int events = 0;
				} manyID;

				struct ided_t{
					int events = 0;
					struct good_t{
						int events = 0;
						struct pass_t{
							int events = 0;
						} pass;
						struct noPass_t{
							int events = 0;
						} noPass;
					} good;
					struct bad_t{
						int events = 0;
						struct pass_t{
							int events = 0;
						} pass;
						struct noPass_t{
							int events = 0;
						} noPass;
					} bad;
				} ided;
			} associated;

			struct ided_t{
				int events = 0;
				struct pass_t{
					int events = 0;
				} pass;
				struct noPass_t{
					int events = 0;
				}noPass;
			} ided;
		} prelim;
	} total;
} pid_res;

int em, ep, pic;
bool flBad;
double xTrue, xFalse;

string printDiv(int a, int b, int c=0, int d=0, int e=0, int f=0){
	stringstream s;

	s.precision(3);
	s << std::fixed;
	s << a;
	if(b!=0) s << "\t\t" << a*100./(double)b << "%";
	if(c!=0) s << "\t\t" << a*100./(double)c << "%";
	if(d!=0) s << "\t\t" << a*100./(double)d << "%";
	if(e!=0) s << "\t\t" << a*100./(double)e << "%";
	if(f!=0) s << "\t\t" << a*100./(double)f << "%";
	return s.str();
}

bool associateMCTracks(){
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
		pid_res.good++;
	}
	else{
		flBad = true;
		pid_res.bad++;
	}

	if((em1==-1 && em2==-1 && em3!=-1) || (em1==-1 && em2!=-1 && em3==-1) || (em1!=-1 && em2==-1 && em3==-1)){
		//OK
		if(em1!=-1) em=em1;
		if(em2!=-1) em=em2;
		if(em3!=-1) em=em3;
		pid_res.good++;
	}
	else{
		pid_res.bad++;
		flBad = true;
	}

	if(!flBad){
		if((ep==0 && em==1) || (ep==1 && em==0)) pic=2;
		else if((ep==0 && em==2) || (ep==2 && em==0)) pic=1;
		else if((ep==1 && em==2) || (ep==2 && em==1)) pic=0;
		pid_res.incAssociated();
	}
	else{
		ep = -1;
		em = -1;
		pic = -1;
	}

	return flBad;
}

int pi0d_identifyPi(int &piCandidate, bool &badElectron){
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

int pid(int &piCandidate, TLorentzVector &gamma){
	TLorentzVector tem;
	TLorentzVector t1ep, t1pi;
	TLorentzVector t2ep, t2pi;
	TLorentzVector ee1, ee2;
	TLorentzVector k1, k2;

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

	//Try e = goodTrack1, pi = goodTrack2
	//cout << endl << corrEvent.pTrack.size() << " " << corrEvent.goodTracks.size() << " " << goodTrack1 << " " << goodTrack2 << endl;
	t1ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Me);
	t1pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Mpic);
	ee1 = tem+t1ep+gamma;
	k1 = ee1+t1pi;

	//Try e = goodTrack2, pi = goodTrack1
	t2pi.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack1]].p, Mpic);
	t2ep.SetVectM(corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[goodTrack2]].p, Me);
	ee2 = tem+t2ep+gamma;
	k2 = ee2+t2pi;

	double Mk;
	double diffpi01, diffpi02;
	double diffk1, diffk2;
	double diff1, diff2;

	diffpi01 = fabs(ee1.M()-Mpi0);
	diffpi02 = fabs(ee2.M()-Mpi0);

	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1) Mk = abcog_params.mkp;
	else Mk = abcog_params.mkn;

	diffk1 = fabs(k1.M()-Mk);
	diffk2 = fabs(k2.M()-Mk);

	diffk1 = 0;
	diffk2 = 0;

	diff1 = diffpi01;// + diffk1;
	diff2 = diffpi02;// + diffk2;

	double pi0DiffLimit = io.cutsDefinition.pid_pi0Diff;
	double kDiffLimit = io.cutsDefinition.pid_kDiff;

	double pi0MSupLim = Mpi0+0.03;
	double kMSupLim = Mk+0.04;
	double kMLowLim1 = 0.3 + ((Mk-0.32)/(Mpi0+0.03))*ee1.M();
	double kMLowLim2 = 0.3 + ((Mk-0.32)/(Mpi0+0.03))*ee2.M();
	//take the smallest ee mass as the ep em, the other is pi
	//if(ee1.M()<(Mpi0+pi0DiffLimit) && k1.M()<(Mk+kDiffLimit)){
	if( (ee1.M()<pi0MSupLim) && (k1.M()<kMSupLim) && (k1.M() > kMLowLim1)){
		nCandidates++;
		piCandidate = goodTrack2;
	}
	//if((ee2.M()-Mpi0)<pi0DiffLimit && fabs(k2.M()-Mk)<kDiffLimit){
	if( (ee2.M()<pi0MSupLim) && (k2.M()<kMSupLim) && (k2.M() > kMLowLim2)){
		nCandidates++;
		piCandidate = goodTrack1;
	}

	/*if(diff1<diff2 && diff1<diff3) piCandidate = 2;
	else if(diff2<diff1 && diff2<diff3) piCandidate = 1;
	else if(diff3<diff2 && diff3<diff1) piCandidate = 0;*/

	meeTotal.Fill(ee1.M());
	meeTotal.Fill(ee2.M());
	mkTotal.Fill(k1.M());
	mkTotal.Fill(k2.M());
	meeepiTotal.Fill(ee1.M(), (tem+t1pi).M());
	meeepiTotal.Fill(ee2.M(), (tem+t2pi).M());
	meekTotal.Fill(ee1.M(), k1.M());
	meekTotal.Fill(ee2.M(), k2.M());

	if(!flBad){
		//Good MC association, fill the plots
		if(pic==goodTrack2){
			meeTrue.Fill(ee1.M());
			mkTrue.Fill(k1.M());
			meeepiTrue.Fill(ee1.M(), (tem+t1pi).M());
			meekTrue.Fill(ee1.M(), k1.M());

			meeDiffTrue.Fill(ee1.M()-Mpi0);
			mkDiffTrue.Fill(k1.M()-Mk);

			meeFalse.Fill(ee2.M());
			mkFalse.Fill(k2.M());
			meeepiFalse.Fill(ee2.M(), (tem+t2pi).M());
			meekFalse.Fill(ee2.M(), k2.M());

			meeDiffFalse.Fill(ee2.M()-Mpi0);
			mkDiffFalse.Fill(k2.M()-Mk);

			xTrue = pow((tem+t1ep).M()/Mpi0, 2.);
			xFalse = pow((tem+t2ep).M()/Mpi0, 2.);
		}
		else if(pic==goodTrack1){
			meeTrue.Fill(ee2.M());
			mkTrue.Fill(k2.M());
			meeepiTrue.Fill(ee2.M(), (tem+t2pi).M());
			meekTrue.Fill(ee2.M(), k2.M());

			meeDiffTrue.Fill(ee2.M()-Mpi0);
			mkDiffTrue.Fill(k2.M()-Mk);

			meeFalse.Fill(ee1.M());
			mkFalse.Fill(k1.M());
			meeepiFalse.Fill(ee1.M(), (tem+t1pi).M());
			meekFalse.Fill(ee1.M(), k1.M());

			meeDiffFalse.Fill(ee1.M()-Mpi0);
			mkDiffFalse.Fill(k1.M()-Mk);

			xTrue = pow((tem+t2ep).M()/Mpi0,2.);
			xFalse = pow((tem+t1ep).M()/Mpi0, 2.);
		}
	}

	return nCandidates;
}

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

int pi0d_trackCombinationVeto(){
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
			if(trackID1==rootPhysics.pic.parentTrack || trackID2==rootPhysics.pic.parentTrack){
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

		// separation from pi impact point >30cm
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

		// separation from undeflected e+ e- trajectories >20cm
		/*propPos = propagate(rootGeom.Lkr.z, rawEvent.track[ep.trackID].bDetPos, rawEvent.track[ep.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e+ :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;

		propPos = propagate(rootGeom.Lkr.z, rawEvent.track[em.trackID].bDetPos, rawEvent.track[em.trackID].bMomentum);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tUndeflected R_LKr_e- :\t" << distance << "\t > 50 : ++" << endl;
		if(distance>50) cond++;*/

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

int pi0d_goodClusters_tight(){
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

		NPhysicsTrack pi = corrEvent.pTrack[rootPhysics.pic.parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[rootPhysics.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[rootPhysics.em.parentTrack];

		// separation from pi impact point >30cm
		propPos = propagateAfter(rootGeom.Lkr.z, pi);
		distance = distance2D(propPos, c.position);
		if(options.isOptDebug()) cout << "\t\tR_LKr_pi :\t\t" << distance << "\t > 50 : ++" << endl;
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
			rootPhysics.gamma.parentCluster = i;
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
bool pi0d_failCutInc(int i, bool assoc, bool good, bool bad){
	pid_res.incNonPass(assoc, good, bad);
	if(!options.isDoScan())corrEvent.failedCond = i;
	if(options.isOptDebug()) cout << "Event is not passing selection " << i << endl;
	if(io.isDoOutput()) io.output.f1() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << " " << i << endl;
	return false;
}
void pi0d_passSelection(){
	if(options.isOptDebug()) cout << "Event is passing selection" << endl;
	if(io.isDoOutput()) io.output.f2() << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << endl;
}

bool nico_pi0DalitzSelect(){
	bool badAcceptance;

	int badCombis=0;
	int piTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int piCandNb;
	bool badElectron;

	int goodClusters;
	//int goodClusterID=-1;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	double pt;

	vector<double> vMass;
	vector<TVector3> vP;

	if(corrEvent.failedCond>=0) {
		outputFileHeader.NFailedEvents--;
		pi0d_failCut(corrEvent.failedCond);
		return false;
	}
	if(options.isOptDebug()) cout << endl;

	// 6) Track DCH time
	if(options.isOptDebug()) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(rootBurst.isData){
		if(options.isOptDebug()) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
		if(options.isOptDebug()) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
		if(options.isOptDebug()) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6); return false;}
	}
	else{
		if(options.isOptDebug()) cout << "\tMC: Not applicable" << endl;
	}

	// 7) Track acceptance veto
	if(options.isOptDebug()) cout << "~~~~ Cut 7 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(options.isOptDebug()) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==io.cutsDefinition.boolBadTrack) {pi0d_failCut(7); return false;}

	// 8) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return false;}


	// 12) Exactly 1 good LKr cluster (loose)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_loose();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(12); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}
	TLorentzVector tempGamma;
	tempGamma.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rawEvent.vtx[corrEvent.goodVertexID].position).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);

	// 13) Photon candidate in LKr acceptance
	if(options.isOptDebug()) cout << "~~~~ Cut 13 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options.isOptDebug()) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=io.cutsDefinition.lkrAcceptance) {pi0d_failCut(13); return false;}

	// 14) E_gamma>3GeV
	if(options.isOptDebug()) cout << "~~~~ Cut 14 ~~~~" << endl;
	if(options.isOptDebug()) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << tempGamma.E() << "\t <= 3 : rejected" << endl;
	if(tempGamma.E()<=io.cutsDefinition.minGammaEnergy) {pi0d_failCut(14); return false;}

	// 15) D_deadcell>2cm
	if(options.isOptDebug()) cout << "~~~~ Cut 15 ~~~~" << endl;
	if(options.isOptDebug()) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t1) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15); return false;}

	// 16) Photon DCH1 intercept >13cm
	if(options.isOptDebug()) cout << "~~~~ Cut 16 ~~~~" << endl;
	propPos = propagate(rootGeom.Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, tempGamma.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options.isOptDebug()) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=io.cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16); return false;}


	//Start Kinematic cuts
	TVector3 totalP; // = t1 + t2 + t3 + tGamma
	totalP =  corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p \
			+ tempGamma.Vect();
	//Pt
	pt = totalP.Perp2(corrEvent.kaonMomentum);

	// 11) Tracks momenta
	if(options.isOptDebug()) cout << "~~~~ Cut 11 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_1 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[0]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_2 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[1]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_3 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[2]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[0]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[0]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[1]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[1]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[2]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[2]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11); return false;}

	// 17) Total momentum 70<p<78
	if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t <70 || >78 : rejected" << endl;
	if(totalP.Mag()<io.cutsDefinition.minTotalMomentum || totalP.Mag()>io.cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return false;}

	// 18) Transverse momentum^2 < 5E-4
	if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
	if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.maxPt << " : rejected" << endl;
	if(pt>=io.cutsDefinition.maxPt) {pi0d_failCut(18); return false;}


	// PID
	flBad = false;
	ep=-1;
	em=-1;
	pic=-1;
	xTrue = 999;
	xFalse = 999;

	bool good=false, bad=false, ided=false;

	flBad = associateMCTracks();

	pid_res.incPrelim();


	// 9) Identify candidates
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	piCandNb = pid(piTrack, tempGamma);
	badElectron = false;
	//piCandNb = pi0d_identifyPi(piTrack, badElectron);
	if(piCandNb==0) pid_res.incNoID(!flBad);
	else if(piCandNb==1){
		if(!flBad){
			if(piTrack == pic){
				pid_res.incGoodId();
				good = true;
			}
			else{
				xTruexFalse.Fill(xTrue, xFalse);
				pid_res.incBadId();
				bad = true;
			}
		}
		pid_res.incIded(!flBad);
		ided = true;
	}
	else if(piCandNb>1) pid_res.incManyID(!flBad);
	nPiCandidatesNew.Fill(piCandNb);
	if(options.isOptDebug()) cout << "Number of pi track candidates :\t" << piCandNb << "\t != 1: rejected" << endl;
	if(piCandNb!=io.cutsDefinition.numPiCandidates) {pi0d_failCutInc(9, !flBad, good, bad); return false;}

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if((int)i==piTrack) continue;

		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	//Create physics event from tracks (e+, e-, pi+-)
	rootPhysics.em.parentTrack = corrEvent.goodTracks[emTrack];
	rootPhysics.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.em.parentVertex = corrEvent.goodVertexID;
	rootPhysics.ep.parentTrack = corrEvent.goodTracks[epTrack];
	rootPhysics.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.ep.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pic.parentTrack = corrEvent.goodTracks[piTrack];
	rootPhysics.pic.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pic.parentVertex = corrEvent.goodVertexID;


	rootPhysics.em.P.SetVectM(corrEvent.pTrack[rootPhysics.em.parentTrack].momentum*corrEvent.pTrack[rootPhysics.em.parentTrack].p, Me);
	rootPhysics.ep.P.SetVectM(corrEvent.pTrack[rootPhysics.ep.parentTrack].momentum*corrEvent.pTrack[rootPhysics.ep.parentTrack].p, Me);
	rootPhysics.pic.P.SetVectM(corrEvent.pTrack[rootPhysics.pic.parentTrack].momentum*corrEvent.pTrack[rootPhysics.pic.parentTrack].p, Mpic);
	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		rootPhysics.pic.pdgID = 211;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		rootPhysics.pic.pdgID = -211;
	}

	// 10) Bad electron cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options.isOptDebug()) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==io.cutsDefinition.boolBadECandidates) {pi0d_failCutInc(10, !flBad, good, bad); return false;}

	// 12) Exactly 1 good LKr cluster (tight)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_tight();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCutInc(12, !flBad, good, bad); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}

	//Add cluster information to physics event (gamma), build reco physics (pi0, kaon)
	rootPhysics.gamma.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.gamma.parentVertex = corrEvent.goodVertexID;
	rootPhysics.gamma.P.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rootPhysics.gamma.vertex).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);
	rootPhysics.pi0.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.pi0.parentVertex = corrEvent.goodVertexID;
	rootPhysics.pi0.P = rootPhysics.em.P + rootPhysics.ep.P + rootPhysics.gamma.P;
	rootPhysics.kaon.P = rootPhysics.pic.P + rootPhysics.pi0.P;

	//Start mass cuts
	// 19) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=io.cutsDefinition.maxPi0MassDiff) {pi0d_failCutInc(19, !flBad, good, bad); return false;}

	// 20) 0.475 < M_pieeg < 0.510
	if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
	if(rootPhysics.kaon.P.M()<io.cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>io.cutsDefinition.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad); return false;}

	//pi0dalitz variables
	rootPhysics.mee = (rootPhysics.em.P + rootPhysics.ep.P).M();
	rootPhysics.x = pow(rootPhysics.mee/Mpi0, 2.);
	rootPhysics.y = 2.*(rootPhysics.em.P.E() - rootPhysics.ep.P.E())/(Mpi0*(1-rootPhysics.x));

	pid_res.incPass(!flBad, good, bad);

	pi0d_passSelection();
	return true;
}

bool newEvent(int i, int &nevt){
	rootPhysics.clear();

	if(options.isOptDebug()) cout << "--------------------------------------------" << endl;
	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}
	if(i==0) cout << endl;
	if(options.isOptDebug()) cout << endl << "--------------------------------------------" << endl;

	pid_res.incTotal();
	++nevt;
	abcog_params = rootBurst.abcog_params;
	if(!options.isDoScan()) return nico_pi0DalitzSelect(); // Normal thing
	else{
		bool result;
		bool defaultResult = false;
		bool globalResult = false;
		io.output.resetResult();
		for(int i=io.cutsDefinition.loadList(0); i!=-1; i=io.cutsDefinition.loadNextList()){
			result = nico_pi0DalitzSelect();
			io.output.newResult(result);
			globalResult |= result;
			if(i==io.cutsDefinition.getDefaultIndex()) defaultResult = result;
		}
		corrEvent.failedCond = defaultResult;
		return globalResult;
	}
}

int main(int argc, char **argv){
	options.parse(argc, argv, io);

	if(options.isDoScan()) io.cutsDefinition.generateLists(options.getScan());
	if(!io.cutsDefinition.addParseCutsFile(options.getCutsFile())) return -1;
	io.cutsDefinition.print();
	if(!options.isDoScan()) io.cutsDefinition.loadDefault(); // No scan, load the default values

	options.printSummary(io);
	options.parseFilter();

	io.openAll(options.isDoScan());

	TTree *pid = new TTree("pid", "pid");
	pid->Branch("pid_res.good", &pid_res.good, "good/I");
	pid->Branch("pid_res.bad", &pid_res.bad, "bad/I");
	pid->Branch("pid_res.total", &pid_res.total.events, "total/I");
	pid->Branch("pid_res.total.prelim", &pid_res.total.prelim.events, "prelim/I");
	pid->Branch("pid_res.total.prelim.associated", &pid_res.total.prelim.associated.events, "associated/I");
	pid->Branch("pid_res.total.prelim.associated.noID", &pid_res.total.prelim.associated.noID.events, "noID/I");
	pid->Branch("pid_res.total.prelim.associated.manyID", &pid_res.total.prelim.associated.manyID.events, "manyID/I");
	pid->Branch("pid_res.total.prelim.associated.ided", &pid_res.total.prelim.associated.ided.events, "ided/I");
	pid->Branch("pid_res.total.prelim.associated.ided.good", &pid_res.total.prelim.associated.ided.good.events, "good/I");
	pid->Branch("pid_res.total.prelim.associated.ided.bad", &pid_res.total.prelim.associated.ided.bad.events, "bad/I");
	pid->Branch("pid_res.total.prelim.noID", &pid_res.total.prelim.noID.events, "noID/I");
	pid->Branch("pid_res.total.prelim.manyID", &pid_res.total.prelim.manyID.events, "manyID/I");
	pid->Branch("pid_res.total.prelim.associated.ided.good.pass", &pid_res.total.prelim.associated.ided.good.pass.events, "pass/I");
	pid->Branch("pid_res.total.prelim.associated.ided.good.noPass", &pid_res.total.prelim.associated.ided.good.noPass.events, "noPass/I");
	pid->Branch("pid_res.total.prelim.associated.ided.bad.pass", &pid_res.total.prelim.associated.ided.bad.pass.events, "pass/I");
	pid->Branch("pid_res.total.prelim.associated.ided.bad.noPass", &pid_res.total.prelim.associated.ided.bad.noPass.events, "noPass/I");
	pid->Branch("pid_res.total.prelim.ided", &pid_res.total.prelim.ided.events, "ided/I");
	pid->Branch("pid_res.total.prelim.ided.pass", &pid_res.total.prelim.ided.pass.events, "pass/I");
	pid->Branch("pid_res.total.prelim.ided.noPass", &pid_res.total.prelim.ided.noPass.events, "noPass/I");

	int nevent = 0;

	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io.input.getNEvents() << endl;
	for(int i=io.input.firstEvent(outputFileHeader); !io.input.eof() && (options.getMaxEvents()<0 || nevent<options.getMaxEvents() ); i = io.input.nextEvent(outputFileHeader)){
		if(newEvent(i, nevent)){
			io.output.fillEvent();
			outputFileHeader.NPassedEvents++;
		}
		else{
			outputFileHeader.NFailedEvents++;
			if(options.isExportAllEvents()) io.output.fillEvent();
		}
	}
	pid->Fill();
	cout << endl << endl;

	cout << "Processed events ->\t" << outputFileHeader.NProcessedEvents << endl;
	cout << "Failed events    ->\t" << outputFileHeader.NFailedEvents << endl;
	cout << "Passed events    ->\t" << outputFileHeader.NPassedEvents << endl;

	cout << endl << "MC Association (track)" << endl << "--------------" << endl;
	cout << "Good = " << pid_res.good << " " << pid_res.good*100./(double)(pid_res.good+pid_res.bad) << endl;
	cout << "Bad  = " << pid_res.bad << " " << pid_res.bad*100./(double)(pid_res.good+pid_res.bad) << endl;

	cout << "\t\t\t\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	cout << "Total: \t\t\t" << printDiv(pid_res.total.events, pid_res.total.events) << endl;
	cout << "|-Prelim: \t\t" << printDiv(pid_res.total.prelim.events, pid_res.total.events, pid_res.total.prelim.events) << endl;
	cout << "  |-NoID: \t\t" << printDiv(pid_res.total.prelim.noID.events, pid_res.total.events, pid_res.total.prelim.events) << endl;
	cout << "  |-ManyID: \t\t" << printDiv(pid_res.total.prelim.manyID.events, pid_res.total.events, pid_res.total.prelim.events) << endl;
	cout << "  |-Ided: \t\t" << printDiv(pid_res.total.prelim.ided.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.ided.events) << endl;
	cout << "    |-Pass: \t\t" << printDiv(pid_res.total.prelim.ided.pass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.ided.events) << endl;
	cout << "    |-NoPass: \t\t" << printDiv(pid_res.total.prelim.ided.noPass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.ided.events) << endl;
	cout << "  |-Associated: \t" << printDiv(pid_res.total.prelim.associated.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events) << endl;
	cout << "  | |-NoID: \t\t" << printDiv(pid_res.total.prelim.associated.noID.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events) << endl;
	cout << "  | |-ManyID: \t\t" << printDiv(pid_res.total.prelim.associated.manyID.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events) << endl;
	cout << "  | |-Ided: \t\t" << printDiv(pid_res.total.prelim.associated.ided.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events) << endl;
	cout << "  |   |-Good: \t\t" << printDiv(pid_res.total.prelim.associated.ided.good.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.good.events) << endl;
	cout << "  |   | |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.pass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.good.events) << endl;
	cout << "  |   | |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.noPass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.good.events) << endl;
	cout << "  |   |-Bad: \t\t" << printDiv(pid_res.total.prelim.associated.ided.bad.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.bad.events) << endl;
	cout << "  |     |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.pass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.bad.events) << endl;
	cout << "  |     |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.noPass.events, pid_res.total.events, pid_res.total.prelim.events, pid_res.total.prelim.associated.events, pid_res.total.prelim.associated.ided.events, pid_res.total.prelim.associated.ided.bad.events) << endl;


	meeTrue.Write();
	mkTrue.Write();
	meeepiTrue.Write();
	meekTrue.Write();

	meeTotal.Write();
	mkTotal.Write();
	meeepiTotal.Write();
	meekTotal.Write();

	meeDiffTrue.Write();
	mkDiffTrue.Write();

	meeFalse.Write();
	mkFalse.Write();
	meeepiFalse.Write();
	meekFalse.Write();

	meeDiffFalse.Write();
	mkDiffFalse.Write();

	nPiCandidates.Write();
	nPiCandidatesNew.Write();
	nBeamSign.Write();
	xTruexFalse.Write();


	pid->Write();


	io.closeAll();
	return 0;
}


