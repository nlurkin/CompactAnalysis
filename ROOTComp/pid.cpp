#define OUTSIDECOMPACT

/// Std include
#include <iomanip>
using namespace std;

/// Compact includes
#include "funLib.h"

// Local includes
#include "SelectionFunctions.h"
#include "pid_res.h"

struct alt_pid_res pid_res;

NRecoParticle *xParticle;
double Mx;
int xPDGId;

bool nico_pi0DalitzSelect(){
	bool badAcceptance;
	
	int badCombis=0;
	int xTrack = -1;
	int epTrack = -1;
	int emTrack = -1;

	int xCandNb;
	bool badElectron;

	int goodClusters;

	int lkrAcceptance;

	TVector3 propPos;
	double radius;

	double pt;

	vector<double> vMass;
	vector<TVector3> vP;

	int firstCutIndex = 0;

	if(corrEvent.failedCond>=0) {
		outputFileHeader.NFailedEvents--;
		pi0d_failCut(corrEvent.failedCond);
		return false;
	}
	if(options.isOptDebug()) cout << endl;


	// 1) Track DCH time
	if(options.isOptDebug()) cout << "~~~~ Cut 1 ~~~~" << endl;
	if(rootBurst.isData){
		if(options.isOptDebug()) cout << "|t_1| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6+firstCutIndex); return false;}
		if(options.isOptDebug()) cout << "|t_2| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6+firstCutIndex); return false;}
		if(options.isOptDebug()) cout << "|t_3| :\t\t\t\t" << fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch) << "\t\t > 25: rejected" << endl;
		if(fabs(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].time - rootBurst.tOffst.Dch)>=io.cutsDefinition.maxTrackTime) {pi0d_failCut(6+firstCutIndex); return false;}
	}
	else{
		if(options.isOptDebug()) cout << "\tMC: Not applicable" << endl;
	}

	// 2) Track acceptance veto
	if(options.isOptDebug()) cout << "~~~~ Cut 2 ~~~~" << endl;
	badAcceptance = pi0d_tracksAcceptance();
	if(options.isOptDebug()) cout << "Track acceptance :\t\t" << badAcceptance << "\t == " << true << ": rejected" << endl;
	if(badAcceptance==io.cutsDefinition.boolBadTrack) {pi0d_failCut(7+firstCutIndex); return false;}


	// 3) Track combination veto loose
	if(options.isOptDebug()) cout << "~~~~ Cut 3 ~~~~" << endl;
//	badCombis = pi0d_trackCombinationVeto(); //Wrong - depends on pi+ id
	badCombis = pi0d_trackCombinationVeto_loose();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8+firstCutIndex); return false;}

	// 4) Michal pre identification
	int nPreCandidates;
	int xPreSelected;
	if(options.isOptDebug()) cout << "~~~~ Cut 4 ~~~~" << endl;
	nPreCandidates = michal_prepid(xPreSelected, OptionsParser::KMU3);
	if(options.isOptDebug()) cout << "Michal pre-id:\t\t " << nPreCandidates << "\t == 0 : rejected" << endl;
	if(nPreCandidates==0) {pi0d_failCut(9+firstCutIndex); return false;}

	// 4) Exactly 1 good LKr cluster
	if(options.isOptDebug()) cout << "~~~~ Cut 4 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_loose();
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCut(12+firstCutIndex); return false;}

	if(rootPhysics.gamma.parentCluster==-1){
		return 0;
	}
	TLorentzVector tempGamma;
	tempGamma.SetVectM((corrEvent.pCluster[rootPhysics.gamma.parentCluster].position - rawEvent.vtx[corrEvent.goodVertexID].position).Unit()*corrEvent.pCluster[rootPhysics.gamma.parentCluster].E, 0.0);

	// 5) Photon candidate in LKr acceptance
	if(options.isOptDebug()) cout << "~~~~ Cut 5 ~~~~" << endl;
	lkrAcceptance = rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].lkr_acc;
	if(options.isOptDebug()) cout << "Mauro condition :\t\t" << lkrAcceptance << "\t != 0 : rejected" << endl;
	if(lkrAcceptance!=io.cutsDefinition.lkrAcceptance) {pi0d_failCut(13+firstCutIndex); return false;}

	// 6) E_gamma>3GeV
	if(options.isOptDebug()) cout << "~~~~ Cut 6 ~~~~" << endl;
	if(options.isOptDebug()) cout << "E_g :\t\t\t\t" << fixed << setprecision(7) << tempGamma.E() << "\t <= 3 : rejected" << endl;
	if(tempGamma.E()<=io.cutsDefinition.minGammaEnergy) {pi0d_failCut(14+firstCutIndex); return false;}

	// 7) D_deadcell>2cm
	if(options.isOptDebug()) cout << "~~~~ Cut 7 ~~~~" << endl;
	if(options.isOptDebug()) cout << "d_deadcell :\t\t\t" << rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.cluster[corrEvent.pCluster[rootPhysics.gamma.parentCluster].clusterID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15+firstCutIndex); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t1) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[0]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15+firstCutIndex); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[1]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15+firstCutIndex); return false;}

	if(options.isOptDebug()) cout << "d_deadcell(t2) :\t\t\t" << rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell << "\t <= 2 : rejected" << endl;
	if(rawEvent.track[corrEvent.pTrack[corrEvent.goodTracks[2]].trackID].dDeadCell<=io.cutsDefinition.minDeadCellDist) {pi0d_failCut(15+firstCutIndex); return false;}

	// 8) Photon DCH1 intercept >13cm
	if(options.isOptDebug()) cout << "~~~~ Cut 8 ~~~~" << endl;
	propPos = propagate(rootGeom.Dch[0].PosChamber.z, corrEvent.pCluster[rootPhysics.gamma.parentCluster].position, tempGamma.Vect());
	radius = sqrt(pow(propPos.X(),2) + pow(propPos.Y(),2));
	if(options.isOptDebug()) cout << "R_gDCH1 :\t\t\t" << radius << "\t <= 13 : rejected" << endl;
	if(radius<=io.cutsDefinition.minGammaDCHRadius) {pi0d_failCut(16+firstCutIndex); return false;}


	//Start Kinematic cuts
	TVector3 totalP; // = t1 + t2 + t3 + tGamma
	totalP =  corrEvent.pTrack[corrEvent.goodTracks[0]].momentum*corrEvent.pTrack[corrEvent.goodTracks[0]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[1]].momentum*corrEvent.pTrack[corrEvent.goodTracks[1]].p \
			+ corrEvent.pTrack[corrEvent.goodTracks[2]].momentum*corrEvent.pTrack[corrEvent.goodTracks[2]].p \
			+ tempGamma.Vect();
	//Pt
	pt = totalP.Perp2(corrEvent.kaonMomentum);
	
	// 9) Tracks momenta
	if(options.isOptDebug()) cout << "~~~~ Cut 9 ~~~~" << endl;
	if(options.isOptDebug()) cout << "p_1 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[0]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_2 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[1]].p << "\t <5 || > 60 : rejected" << endl;
	if(options.isOptDebug()) cout << "p_3 :\t\t\t\t" << corrEvent.pTrack[corrEvent.goodTracks[2]].p << "\t <5 || > 60 : rejected" << endl;
	if(corrEvent.pTrack[corrEvent.goodTracks[0]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[0]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11+firstCutIndex); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[1]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[1]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11+firstCutIndex); return false;}
	if(corrEvent.pTrack[corrEvent.goodTracks[2]].p<=io.cutsDefinition.minTrackMomentum || corrEvent.pTrack[corrEvent.goodTracks[2]].p>=io.cutsDefinition.maxTrackMomentum) {pi0d_failCut(11+firstCutIndex); return false;}

	if(options.getSelectionType()==OptionsParser::K2PI){
		// 10) Total momentum 70<p<78
		if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
		if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t <70 || >78 : rejected" << endl;
		if(totalP.Mag()<io.cutsDefinition.k2pi.minTotalMomentum || totalP.Mag()>io.cutsDefinition.k2pi.maxTotalMomentum) {pi0d_failCut(17+firstCutIndex); return false;}

		// 11) Transverse momentum^2 < 5E-4
		if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
		if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.k2pi.maxPt << " : rejected" << endl;
		if(pt>=io.cutsDefinition.k2pi.maxPt) {pi0d_failCut(18+firstCutIndex); return false;}
	}
	else{
		// 17) Total momentum p<78
		if(options.isOptDebug()) cout << "~~~~ Cut 17 ~~~~" << endl;
		if(options.isOptDebug()) cout << "p_tot :\t\t\t" << totalP.Mag() << "\t >78 : rejected" << endl;
		if(totalP.Mag()>io.cutsDefinition.kmu3.maxTotalMomentum) {pi0d_failCut(17+firstCutIndex); return false;}
		
		// 18) Transverse momentum^2 < 5E-4
		if(options.isOptDebug()) cout << "~~~~ Cut 18 ~~~~" << endl;
		if(options.isOptDebug()) cout << "P_t^2 :\t\t" << pt << "\t >= " << io.cutsDefinition.kmu3.maxPt << " || <= " << io.cutsDefinition.kmu3.minPt << " : rejected" << endl;
		if(pt <= io.cutsDefinition.kmu3.minPt || pt>=io.cutsDefinition.kmu3.maxPt) {pi0d_failCut(18+firstCutIndex); return false;}
	}


	// PID
	flBad = false;
	ep=-1;
	em=-1;
	xPart=-1;
	xTrue = 999;
	xFalse = 999;

	bool good=false, bad=false;

	flBad = associateMCTracks(pid_res, NULL);
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
	xCandNb = pid(xTrack, tempGamma, options.getSelectionType());
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
	}
	else if(xCandNb>1){
		pid_res.incManyID(!flBad);
		xMCManyID.Fill(rootMC.xTrue);
		xTruexFalseMany.Fill(xTrue, xFalse);
		xTruexMCMany.Fill(xTrue, rootMC.xTrue);
	}
	nxCandidatesNew.Fill(xCandNb);
	if(options.isOptDebug()) cout << "Number of x track candidates :\t" << xCandNb << "\t != 1: rejected" << endl;
	if(xCandNb!=io.cutsDefinition.numXCandidates) {pi0d_failCutInc(9, !flBad, good, bad, &pid_res); return false;}

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
	xParticle->parentTrack = corrEvent.goodTracks[xTrack];
	xParticle->vertex = rawEvent.vtx[corrEvent.goodVertexID].position;
	xParticle->parentVertex = corrEvent.goodVertexID;


	rootPhysics.em.P.SetVectM(corrEvent.pTrack[rootPhysics.em.parentTrack].momentum*corrEvent.pTrack[rootPhysics.em.parentTrack].p, Me);
	rootPhysics.ep.P.SetVectM(corrEvent.pTrack[rootPhysics.ep.parentTrack].momentum*corrEvent.pTrack[rootPhysics.ep.parentTrack].p, Me);
	xParticle->P.SetVectM(corrEvent.pTrack[xParticle->parentTrack].momentum*corrEvent.pTrack[xParticle->parentTrack].p, Mx);
	// Select event charge
	if(rawEvent.vtx[corrEvent.goodVertexID].charge==1){
		rootPhysics.kaon.pdgID = 321;
		xParticle->pdgID = xPDGId;
	}
	else{
		rootPhysics.kaon.pdgID = -321;
		xParticle->pdgID = -xPDGId;
	}

	//plot e/p
	double xeop = corrEvent.pTrack[xParticle->parentTrack].E/xParticle->P.Vect().Mag();
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
	if(badElectron==io.cutsDefinition.boolBadECandidates) {pi0d_failCutInc(10, !flBad, good, bad, &pid_res); return false;}

	/*// 11) Track combination veto
	if(options.isOptDebug()) cout << "~~~~ Cut 11 ~~~~" << endl;
	badCombis = pi0d_trackCombinationVeto_tight();
	if(options.isOptDebug()) cout << "Bad track combination :\t\t" << badCombis << "\t != 0: rejected" << endl;
	if(badCombis!=io.cutsDefinition.numBadTrackCombi) {pi0d_failCut(8+firstCutIndex); return false;}*/
	
	// 12) Exactly 1 good LKr cluster (tight)
	if(options.isOptDebug()) cout << "~~~~ Cut 12 ~~~~" << endl;
	goodClusters = pi0d_goodClusters_tight(*xParticle, rootPhysics);
	if(options.isOptDebug()) cout << "Good LKr clusters :\t\t" << goodClusters << "\t != 1 : rejected" << endl;
	if(goodClusters!=io.cutsDefinition.numAddGoodCluster) {pi0d_failCutInc(12, !flBad, good, bad, &pid_res); return false;}

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
	rootPhysics.kaon.P = xParticle->P + rootPhysics.pi0.P;

	//Start mass cuts
	// 19) |M_eeg - M_pi0|<8 MeV
	if(options.isOptDebug()) cout << "~~~~ Cut 19 ~~~~" << endl;
	if(options.isOptDebug()) cout << "M_ee :\t\t" << rootPhysics.mee << endl;
	if(options.isOptDebug()) cout << "|M_eeg - M_pi0| :\t\t" << fabs(rootPhysics.pi0.P.M()-Mpi0) << "\t >= 0.008 : rejected" << endl;
	if(fabs(rootPhysics.pi0.P.M()-Mpi0)>=io.cutsDefinition.k2pi.maxPi0MassDiff) {pi0d_failCutInc(19, !flBad, good, bad, &pid_res); return false;}


	if(options.getSelectionType()==OptionsParser::K2PI){
		// 20) 0.475 < M_pieeg < 0.510
		if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
		if(options.isOptDebug()) cout << "M_pieeg :\t\t" << rootPhysics.kaon.P.M() << "\t <0.475 || >0.510: rejected" << endl;
		if(fabs(rootPhysics.kaon.P.M() - abcog_params.mkp) > io.cutsDefinition.k2pi.maxKaonMassDiff) {pi0d_failCutInc(20, !flBad, good, bad, &pid_res); return false;}
	}
	else{
		// 20) Missing mass square < 0.01 MeV^2
		double mmasssq = (rootPhysics.kaon.P - rootPhysics.mu.P - rootPhysics.pi0.P).M2();
		if(options.isOptDebug()) cout << "~~~~ Cut 20 ~~~~" << endl;
		if(options.isOptDebug()) cout << "M_miss^2:\t\t" << mmasssq << "\t >0.01: rejected" << endl;
		if(mmasssq > io.cutsDefinition.kmu3.maxMissMassSq) {pi0d_failCutInc(20, !flBad, good, bad, &pid_res); return false;}
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
	//Clear the event
	rootPhysics.clear();

	// Filter events
	if(options.getBadEventsList().size()>0){
		if(!isFilteredEvent(rootBurst.nrun, rootBurst.time, rawEvent.timeStamp, options.getBadEventsList())) return false;
	}

	// Debugging info
	if(options.isOptDebug()) cout << "--------------------------------------------" << endl;
	if(options.getPeriodKeep()!= 0 && rootBurst.period!=options.getPeriodKeep()) return false;
	if(i==0) cout << "First event: ";
	if(i % options.getOutputModulo() == 0) cout << i << " " << rootBurst.nrun << " " << rootBurst.time << " " << rawEvent.timeStamp << "                 \r" << flush;
	if(i==0) cout << endl;
	if(options.isOptDebug()) cout << endl << "--------------------------------------------" << endl;

	pid_res.incTotal();
	++nevt;
	abcog_params = rootBurst.abcog_params;
	if(!options.isDoScan()) return nico_pi0DalitzSelect(); // Normal thing
	else{
		// Scan cuts
		bool result;
		bool defaultResult = false;
		bool globalResult = false;

		// Reset result
		io.output.resetResult();
		pid_res.ResetEntry();

		// Loop over all cuts definition
		for(int i=io.cutsDefinition.loadList(0); i!=-1; i=io.cutsDefinition.loadNextList()){
			result = nico_pi0DalitzSelect();
			io.output.newResult(result);

			// Set as success if pass at least one cut
			globalResult |= result;
			if(i==io.cutsDefinition.getDefaultIndex()) defaultResult = result;
			pid_res.NewEntry();
		}
		corrEvent.failedCond = defaultResult;
		return globalResult;
	}
}

int main(int argc, char **argv){

	// Get options

	if(!options.parse(argc, argv, io)) return 0;

	// Generate cuts
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

	// Define some values specific to each selection
	if(options.getSelectionType()==OptionsParser::K2PI){
		Mx = Mpic;
		xParticle = &rootPhysics.pic;
	}
	else {
		Mx = Mmu;
		xParticle = &rootPhysics.mu;
	}

	// Create the output pid_res tree
	TTree *pid = new TTree("pid", "pid");
	branchResTree(pid, pid_res);

	// Initialise event loop
	int nevent = 0;
	outputFileHeader.NPassedEvents = 0;
	cout << "Entries in the tree: " << io.input.getNEvents() << endl;

	// Event loop
	for(int i=io.input.firstEvent(outputFileHeader, options.getStartEvent());
			!io.input.eof() && (options.getMaxEvents()<0 || nevent<options.getMaxEvents() );
			i = io.input.nextEvent(outputFileHeader)){
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

	// Print pid result structure
	cout << endl << endl;
	printResStruct(pid_res, NULL);

	pid->Write();
	io.closeAll();
	return 0;
}


