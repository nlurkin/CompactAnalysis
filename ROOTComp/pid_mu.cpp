#define OUTSIDECOMPACT

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

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

#define PRINTVAR(v) #v << "= " << v << " "

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
double Mx;
int xPDGId;
NRecoParticle *xRootPhysics;

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

struct alt_pid_res{
	void incTotal(){ total.events[currentID]++;}
	void incPrelim(){ total.prelim.events[currentID]++;}
	void incNoID(bool assoc){
		if(assoc) total.prelim.associated.noID.events[currentID]++;
		else total.prelim.noID.events[currentID]++;
	}
	void incManyID(bool assoc){
		if(assoc) total.prelim.associated.manyID.events[currentID]++;
		else total.prelim.manyID.events[currentID]++;
	}
	void incAssociated(){
		total.prelim.associated.events[currentID]++;
	}
	void incIded(bool assoc){
		if(assoc) total.prelim.associated.ided.events[currentID]++;
		else total.prelim.ided.events[currentID]++;
	}
	void incGoodId(){
		total.prelim.associated.ided.good.events[currentID]++;
	}
	void incBadId(){
		total.prelim.associated.ided.bad.events[currentID]++;
	}
	void incPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.pass.events[currentID]++;
			else if(bad) total.prelim.associated.ided.bad.pass.events[currentID]++;
		}
		else{
			total.prelim.ided.pass.events[currentID]++;
		}
	}
	void incNonPass(bool assoc, bool good, bool bad){
		if(assoc){
			if(good) total.prelim.associated.ided.good.noPass.events[currentID]++;
			else if(bad) total.prelim.associated.ided.bad.noPass.events[currentID]++;
		}
		else{
			total.prelim.ided.noPass.events[currentID]++;
		}
	}

	void NewEntry(){
		currentID++;
	}
	void ResetEntry(){
		currentID=0;
	}

	void Init(int N){
		good.resize(N, 0);
		bad.resize(N, 0);
		total.events.resize(N, 0);
		total.prelim.events.resize(N, 0);
		total.prelim.noID.events.resize(N, 0);
		total.prelim.manyID.events.resize(N, 0);
		total.prelim.associated.events.resize(N, 0);
		total.prelim.associated.noID.events.resize(N, 0);
		total.prelim.associated.manyID.events.resize(N, 0);
		total.prelim.associated.ided.events.resize(N, 0);
		total.prelim.associated.ided.good.events.resize(N, 0);
		total.prelim.associated.ided.good.pass.events.resize(N, 0);
		total.prelim.associated.ided.good.noPass.events.resize(N, 0);
		total.prelim.associated.ided.bad.events.resize(N, 0);
		total.prelim.associated.ided.bad.noPass.events.resize(N, 0);
		total.prelim.associated.ided.bad.pass.events.resize(N, 0);
		total.prelim.ided.events.resize(N, 0);
		total.prelim.ided.pass.events.resize(N, 0);
		total.prelim.ided.noPass.events.resize(N, 0);
	}

	int currentID=0;
	//MC association
	vector<int> good, bad;

	struct total_t{
		vector<int> events;

		struct prelim_t{
			vector<int> events;
			struct noID_t{
				vector<int> events;
			} noID;
			struct manyID_t{
				vector<int> events;
			} manyID;

			struct associated_t{
				vector<int> events;
				struct noID_t{
					vector<int> events;
				} noID;
				struct manyID_t{
					vector<int> events;
				} manyID;

				struct ided_t{
					vector<int> events;
					struct good_t{
						vector<int> events;
						struct pass_t{
							vector<int> events;
						} pass;
						struct noPass_t{
							vector<int> events;
						} noPass;
					} good;
					struct bad_t{
						vector<int> events;
						struct pass_t{
							vector<int> events;
						} pass;
						struct noPass_t{
							vector<int> events;
						} noPass;
					} bad;
				} ided;
			} associated;

			struct ided_t{
				vector<int> events;
				struct pass_t{
					vector<int> events;
				} pass;
				struct noPass_t{
					vector<int> events;
				}noPass;
			} ided;
		} prelim;
	} total;
} pid_res;

int em, ep, xPart;
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
	//double limit = 1e-4;
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

	/*if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[0]].p)<limit) ep1=0;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[1]].p)<limit) ep2=1;
	if((rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[2]].p)<limit) ep3=2;

	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[0]].p)<limit) em1=0;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[1]].p)<limit) em2=1;
	if((rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[2]].p)<limit) em3=2;*/

	//cout << (rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[0]].p) << endl;
	//cout << (rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[1]].p) << endl;
	//cout << (rootMC.ep.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()*fabs(rootMC.ep.P.E()-corrEvent.pTrack[corrEvent.goodTracks[2]].p) << endl;
	//cout << (rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[0]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[0]].p) << endl;
	//cout << (rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[1]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[1]].p) << endl;
	//cout << (rootMC.em.P.Vect().Unit()-corrEvent.pTrack[corrEvent.goodTracks[2]].momentum).Mag()*fabs(rootMC.em.P.E()-corrEvent.pTrack[corrEvent.goodTracks[2]].p) << endl;
	//cout << PRINTVAR(ep1) << PRINTVAR(ep2) << PRINTVAR(ep3) << endl;
	//cout << PRINTVAR(em1) << PRINTVAR(em2) << PRINTVAR(em3) << endl;
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
		//cout << PRINTVAR(ep) << PRINTVAR(em) << endl;
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

int pid(int &xCandidate, TLorentzVector &gamma){
	TLorentzVector tem;
	TLorentzVector t1ep, t1x;
	TLorentzVector t2ep, t2x;
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
	
	//Try e = goodTrack1, x = goodTrack2
	//cout << endl << corrEvent.pTrack.size() << " " << corrEvent.goodTracks.size() << " " << goodTrack1 << " " << goodTrack2 << endl;
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

	//double pi0MSupLim = Mpi0+0.03;
	//double kMSupLim = Mk+0.04;
	//double kMLowLim1 = 0.3 + ((Mk-0.32)/(Mpi0+0.03))*ee1.M();
	//double kMLowLim2 = 0.3 + ((Mk-0.32)/(Mpi0+0.03))*ee2.M();
	//take the smallest ee mass as the ep em, the other is pi
	//if(ee1.M()<(Mpi0+pi0DiffLimit) && k1.M()<(Mk+kDiffLimit)){
	//if( (ee1.M()<pi0MSupLim) && (k1.M()<kMSupLim) && (k1.M() > kMLowLim1)){
	if(options.getSelectionType()==OptionsParser::K2PI){
		if( (diffpi01<io.cutsDefinition.maxPi0MassDiff) && diffk1<io.cutsDefinition.maxKaonMassDiff){
			nCandidates++;
			xCandidate = goodTrack2;
		}
		//if((ee2.M()-Mpi0)<pi0DiffLimit && fabs(k2.M()-Mk)<kDiffLimit){
		//if( (ee2.M()<pi0MSupLim) && (k2.M()<kMSupLim) && (k2.M() > kMLowLim2)){
		if( (diffpi02<io.cutsDefinition.maxPi0MassDiff) && diffk2<io.cutsDefinition.maxKaonMassDiff){
			nCandidates++;
			xCandidate = goodTrack1;
		}
	}
	else{
		if( (diffpi01<io.cutsDefinition.maxPi0MassDiffMu)){
			nCandidates++;
			xCandidate = goodTrack2;
		}
		//if((ee2.M()-Mpi0)<pi0DiffLimit && fabs(k2.M()-Mk)<kDiffLimit){
		//if( (ee2.M()<pi0MSupLim) && (k2.M()<kMSupLim) && (k2.M() > kMLowLim2)){
		if( (diffpi02<io.cutsDefinition.maxPi0MassDiffMu)){
			nCandidates++;
			xCandidate = goodTrack1;
		}
	}

	/*if(diff1<diff2 && diff1<diff3) piCandidate = 2;
	else if(diff2<diff1 && diff2<diff3) piCandidate = 1;
	else if(diff3<diff2 && diff3<diff1) piCandidate = 0;*/

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

	if(diffpi01<io.cutsDefinition.maxPi0MassDiff && diffpi02<io.cutsDefinition.maxPi0MassDiff){
		combi2Dk_exclu.Fill(k1.M(), k2.M());
	}

	//cout << PRINTVAR(flBad) << PRINTVAR(xPart) << PRINTVAR(goodTrack1) << PRINTVAR(goodTrack2) << endl;
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

int pi0d_trackCombinationVeto_tight(){
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
			if(trackID1==xRootPhysics->parentTrack || trackID2==xRootPhysics->parentTrack){
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

		NPhysicsTrack x = corrEvent.pTrack[xRootPhysics->parentTrack];
		NPhysicsTrack ep = corrEvent.pTrack[rootPhysics.ep.parentTrack];
		NPhysicsTrack em = corrEvent.pTrack[rootPhysics.em.parentTrack];

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

bool nico_pi0DalitzSelect(){
	bool badAcceptance;

	int badCombis=0;
	int xTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int xCandNb;
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
//	badCombis = pi0d_trackCombinationVeto(); //Wrong - depends on pi+ id
	badCombis = pi0d_trackCombinationVeto_loose();
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

	if(options.getSelectionType()==OptionsParser::K2PI){
		// 17) Total momentum 70<p<78
		if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
		if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t <70 || >78 : rejected" << endl;
		if(totalP.Mag()<io.cutsDefinition.minTotalMomentum || totalP.Mag()>io.cutsDefinition.maxTotalMomentum) {pi0d_failCut(17); return false;}

		// 18) Transverse momentum^2 < 5E-4
		if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
		if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.maxPt << " : rejected" << endl;
		if(pt>=io.cutsDefinition.maxPt) {pi0d_failCut(18); return false;}
	}
	else{
		// 17) Total momentum p<78
		if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
		if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t >78 : rejected" << endl;
		if(totalP.Mag()>io.cutsDefinition.maxTotalMuMomentum) {pi0d_failCut(17); return false;}
		
		// 18) Transverse momentum^2 < 5E-4
		if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
		if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.maxMuPt << " || <= " << io.cutsDefinition.minMuPt << " : rejected" << endl;
		if(pt <= io.cutsDefinition.minMuPt || pt>=io.cutsDefinition.maxMuPt) {pi0d_failCut(18); return false;}
	}


	// PID
	flBad = false;
	ep=-1;
	em=-1;
	xPart=-1;
	xTrue = 999;
	xFalse = 999;

	bool good=false, bad=false, ided=false;

	flBad = associateMCTracks();
	pid_res.incPrelim();


	vector<float> eopSort = sortEOP();
	eoplowest.Fill(eopSort[0]);
	if(eopSort[0]<0.85){
		eopsecond.Fill(eopSort[1]);
		if(eopSort[1]>0.85 && eopSort[1]<1.15){
			eophighest.Fill(eopSort[2]);
		}
	}

	// 9) Identify candidates
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	xCandNb = pid(xTrack, tempGamma);
	badElectron = false;
	//piCandNb = pi0d_identifyPi(piTrack, badElectron);
	if(xCandNb==0){
		pid_res.incNoID(!flBad);
		xMCNoID.Fill(rootMC.xTrue);
		xTruexFalseNo.Fill(xTrue, xFalse);
		xTruexMCNo.Fill(xTrue, rootMC.xTrue);
	}
	else if(xCandNb==1){
		if(!flBad){
			if(xTrack == xPart){
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
	else if(xCandNb>1){
		pid_res.incManyID(!flBad);
		xMCManyID.Fill(rootMC.xTrue);
		xTruexFalseMany.Fill(xTrue, xFalse);
		xTruexMCMany.Fill(xTrue, rootMC.xTrue);
	}
	nxCandidatesNew.Fill(xCandNb);
	if(options.isOptDebug()) cout << "Number of x track candidates :\t" << xCandNb << "\t != 1: rejected" << endl;
	if(xCandNb!=io.cutsDefinition.numPiCandidates) {pi0d_failCutInc(9, !flBad, good, bad); return false;}

	for(unsigned int i=0; i<corrEvent.goodTracks.size(); i++){
		if((int)i==xTrack) continue;

		if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==+1 && epTrack==-1) epTrack=i;
		else if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[i]].trackID].q==-1 && emTrack==-1) emTrack=i;
		else{
			if(epTrack==-1) epTrack=i;
			if(emTrack==-1) emTrack=i;
		}
	}

	//Create physics event from tracks (e+, e-, x+-)
	rootPhysics.em.parentTrack = corrEvent.goodTracks[emTrack];
	rootPhysics.em.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.em.parentVertex = corrEvent.goodVertexID;
	rootPhysics.ep.parentTrack = corrEvent.goodTracks[epTrack];
	rootPhysics.ep.vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	rootPhysics.ep.parentVertex = corrEvent.goodVertexID;
	xRootPhysics->parentTrack = corrEvent.goodTracks[xTrack];
	xRootPhysics->vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	xRootPhysics->parentVertex = corrEvent.goodVertexID;


	rootPhysics.em.P.SetVectM(corrEvent.pTrack[rootPhysics.em.parentTrack].momentum*corrEvent.pTrack[rootPhysics.em.parentTrack].p, Me);
	rootPhysics.ep.P.SetVectM(corrEvent.pTrack[rootPhysics.ep.parentTrack].momentum*corrEvent.pTrack[rootPhysics.ep.parentTrack].p, Me);
	xRootPhysics->P.SetVectM(corrEvent.pTrack[xRootPhysics->parentTrack].momentum*corrEvent.pTrack[xRootPhysics->parentTrack].p, Mx);
	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		xRootPhysics->pdgID = xPDGId;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		xRootPhysics->pdgID = -xPDGId;
	}

	//plot e/p
	double xeop = corrEvent.pTrack[xRootPhysics->parentTrack].E/xRootPhysics->P.Vect().Mag();
	double epeop = corrEvent.pTrack[rootPhysics.ep.parentTrack].E/rootPhysics.ep.P.Vect().Mag();
	double emeop = corrEvent.pTrack[rootPhysics.em.parentTrack].E/rootPhysics.em.P.Vect().Mag();
	eopx.Fill(xeop);
	eope.Fill(epeop);
	eope.Fill(emeop);

	if(xeop < 0.85){
		eope_goodx.Fill(epeop);
		eope_goodx.Fill(emeop);
	}
	else{
		eope_badx.Fill(epeop);
		eope_badx.Fill(emeop);
	}
	if(epeop > 0.85 && epeop < 1.15 && emeop > 0.85 && emeop < 1.15) eopx_goode.Fill(xeop);
	else eopx_bade.Fill(xeop);

	// 10) Bad electron cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 10 ~~~~" << endl;
	if(options.isOptDebug()) cout << "Bad electron tracks eop :\t" << badElectron << "\t == " << true << ": rejected" << endl;
	if(badElectron==io.cutsDefinition.boolBadECandidates) {pi0d_failCutInc(10, !flBad, good, bad); return false;}

	/*// 8) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto_tight();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8); return false;}*/
	
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
	rootPhysics.kaon.P = xRootPhysics->P + rootPhysics.pi0.P;

	//Start mass cuts
	// 19) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=io.cutsDefinition.maxPi0MassDiff) {pi0d_failCutInc(19, !flBad, good, bad); return false;}


	if(options.getSelectionType()==OptionsParser::K2PI){
		// 20) 0.475 < M_pieeg < 0.510
		if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
		if(options.isOptDebug()) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
		//if(rootPhysics.kaon.P.M()<io.cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>io.cutsDefinition.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad); return false;}
		if(fabs(rootPhysics.kaon.P.M() - abcog_params.mkp) > io.cutsDefinition.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad); return false;}
	}
	else{
		// 20) Missing mass square < 0.01 MeV^2
		double mmasssq = (rootPhysics.kaon.P - rootPhysics.mu.P - rootPhysics.pi0.P).M2();
		if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
		if(options.isOptDebug()) cout << "M_miss^2:\t\t" << mmasssq << "\t >0.01: rejected" << endl;
		if(mmasssq > io.cutsDefinition.maxMissMassSq) {pi0d_failCutInc(20, !flBad, good, bad); return false;}
		// 20) 0.475 < M_mueeg < 0.510
		//if(rootPhysics.kaon.P.M()<io.cutsDefinition.minKaonMassDiff || rootPhysics.kaon.P.M()>io.cutsDefinition.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad); return false;}
		//if(fabs(rootPhysics.kaon.P.M() - abcog_params.mkp) > io.cutsDefinition.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad); return false;}
	}
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

	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}

	if(options.isOptDebug()) cout << "--------------------------------------------" << endl;
	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(i==0) cout << endl;
	if(options.isOptDebug()) cout << endl << "--------------------------------------------" << endl;

	pid_res.incTotal();
	++nevt;
	abcog_params = rootBurst.abcog_params;
	if(!options.isDoScan()) return nico_pi0DalitzSelect(i); // Normal thing
	else{
		bool result;
		bool defaultResult = false;
		bool globalResult = false;
		io.output.resetResult();
		pid_res.ResetEntry();
		for(int i=io.cutsDefinition.loadList(0); i!=-1; i=io.cutsDefinition.loadNextList()){
			result = nico_pi0DalitzSelect(i);
			io.output.newResult(result);
			globalResult |= result;
			if(i==io.cutsDefinition.getDefaultIndex()) defaultResult = result;
			pid_res.NewEntry();
		}
		corrEvent.failedCond = defaultResult;
		return globalResult;
	}
}

int main(int argc, char **argv){
	if(!options.parse(argc, argv, io)) return 0;

	if(options.isDoScan()) io.cutsDefinition.generateLists(options.getScan());
	if(!io.cutsDefinition.addParseCutsFile(options.getCutsFile())) return -1;
	io.cutsDefinition.print();
	if(!options.isDoScan()){
		io.cutsDefinition.loadDefault(); // No scan, load the default values
		pid_res.Init(1);
	}
	else pid_res.Init(options.getScan());

	options.printSummary(io);
	options.parseFilter();

	io.openAll(options.isDoScan());

	if(options.getSelectionType()==OptionsParser::K2PI){
		Mx = Mpic;
		xRootPhysics = &rootPhysics.pic;
	}
	else {
		Mx = Mmu;
		xRootPhysics = &rootPhysics.mu;
	}

	TTree *pid = new TTree("pid", "pid");
	pid->Branch("pid_res.good", &pid_res.good);
	pid->Branch("pid_res.bad", &pid_res.bad);
	pid->Branch("pid_res.total", &pid_res.total.events);
	pid->Branch("pid_res.total.prelim", &pid_res.total.prelim.events);
	pid->Branch("pid_res.total.prelim.associated", &pid_res.total.prelim.associated.events);
	pid->Branch("pid_res.total.prelim.associated.noID", &pid_res.total.prelim.associated.noID.events);
	pid->Branch("pid_res.total.prelim.associated.manyID", &pid_res.total.prelim.associated.manyID.events);
	pid->Branch("pid_res.total.prelim.associated.ided", &pid_res.total.prelim.associated.ided.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good", &pid_res.total.prelim.associated.ided.good.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad", &pid_res.total.prelim.associated.ided.bad.events);
	pid->Branch("pid_res.total.prelim.noID", &pid_res.total.prelim.noID.events);
	pid->Branch("pid_res.total.prelim.manyID", &pid_res.total.prelim.manyID.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good.pass", &pid_res.total.prelim.associated.ided.good.pass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.good.noPass", &pid_res.total.prelim.associated.ided.good.noPass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad.pass", &pid_res.total.prelim.associated.ided.bad.pass.events);
	pid->Branch("pid_res.total.prelim.associated.ided.bad.noPass", &pid_res.total.prelim.associated.ided.bad.noPass.events);
	pid->Branch("pid_res.total.prelim.ided", &pid_res.total.prelim.ided.events);
	pid->Branch("pid_res.total.prelim.ided.pass", &pid_res.total.prelim.ided.pass.events);
	pid->Branch("pid_res.total.prelim.ided.noPass", &pid_res.total.prelim.ided.noPass.events);

	int nevent = 0;

	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io.input.getNEvents() << endl;
	for(int i=io.input.firstEvent(outputFileHeader, options.getStartEvent()); !io.input.eof() && (options.getMaxEvents()<0 || nevent<options.getMaxEvents() ); i = io.input.nextEvent(outputFileHeader)){
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
	cout << "Good = " << pid_res.good[0] << " " << pid_res.good[0]*100./(double)(pid_res.good[0]+pid_res.bad[0]) << endl;
	cout << "Bad  = " << pid_res.bad[0] << " " << pid_res.bad[0]*100./(double)(pid_res.good[0]+pid_res.bad[0]) << endl;

	cout << "\t\t\t\t\t0\t\t1\t\t2\t\t3\t\t4" << endl;
	cout << "Total: \t\t\t" << printDiv(pid_res.total.events[0], pid_res.total.events[0]) << endl;
	cout << "|-Prelim: \t\t" << printDiv(pid_res.total.prelim.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << endl;
	cout << "  |-NoID: \t\t" << printDiv(pid_res.total.prelim.noID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << endl;
	cout << "  |-ManyID: \t\t" << printDiv(pid_res.total.prelim.manyID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0]) << endl;
	cout << "  |-Ided: \t\t" << printDiv(pid_res.total.prelim.ided.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << endl;
	cout << "    |-Pass: \t\t" << printDiv(pid_res.total.prelim.ided.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << endl;
	cout << "    |-NoPass: \t\t" << printDiv(pid_res.total.prelim.ided.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.ided.events[0]) << endl;
	cout << "  |-Associated: \t" << printDiv(pid_res.total.prelim.associated.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << endl;
	cout << "  | |-NoID: \t\t" << printDiv(pid_res.total.prelim.associated.noID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << endl;
	cout << "  | |-ManyID: \t\t" << printDiv(pid_res.total.prelim.associated.manyID.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0]) << endl;
	cout << "  | |-Ided: \t\t" << printDiv(pid_res.total.prelim.associated.ided.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0]) << endl;
	cout << "  |   |-Good: \t\t" << printDiv(pid_res.total.prelim.associated.ided.good.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << endl;
	cout << "  |   | |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << endl;
	cout << "  |   | |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.good.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.good.events[0]) << endl;
	cout << "  |   |-Bad: \t\t" << printDiv(pid_res.total.prelim.associated.ided.bad.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << endl;
	cout << "  |     |-Pass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.pass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << endl;
	cout << "  |     |-NoPass: \t" << printDiv(pid_res.total.prelim.associated.ided.bad.noPass.events[0], pid_res.total.events[0], pid_res.total.prelim.events[0], pid_res.total.prelim.associated.events[0], pid_res.total.prelim.associated.ided.events[0], pid_res.total.prelim.associated.ided.bad.events[0]) << endl;


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

	pid->Write();


	io.closeAll();
	return 0;
}


