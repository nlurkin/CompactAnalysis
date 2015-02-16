
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <math.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <iterator>
#include <bitset>
#include <iomanip>
using namespace std;


int nico_pi0DalitzAna(superBurst *sbur,superCmpEvent *sevt){
	//TH1D *nTracks = (TH1D*)gFile->Get("NTracks");
	//TH1D *nVertices = (TH1D*)gFile->Get("NVertices");
	//TH1D *DCH1Rad = (TH1D*)gFile->Get("DCH1Rad");

	//TH1D *sel_nTracks = (TH1D*)gFile->Get("sel_NTracks");
	//TH1D *sel_nVertices = (TH1D*)gFile->Get("sel_NVertices");
	//TH1D *sel_DCH1Rad = (TH1D*)gFile->Get("sel_DCH1Rad");

	//TH1D *eemass = (TH1D*)gFile->Get("eeMass");
	//TH1D *peegmass = (TH1D*)gFile->Get("peegMass");
	//TH1D *eegmass = (TH1D*)gFile->Get("eegMass");

	//TH1D *xDistrib = (TH1D*)gFile->Get("xDistrib");

	//sel_nTracks->Add(nTracks);
	//sel_nVertices->Add(nVertices);
	//sel_DCH1Rad->Add(DCH1Rad);

	//double x = pow(fullEvent.mee/Mpi0,2);
	double weight=1.;
	/*if(mcOnly){
		if(ffWeightType!=-1) weight /= (1+0.032*x);
		if(ffWeightType==1) weight *= x;
		else if(ffWeightType==2) weight *= x*x;
	}*/

	//eemass->Fill(fullEvent.mee, weight*fullEvent.weight);
	//eegmass->Fill(fullEvent.mPi0, weight*fullEvent.weight);
	//peegmass->Fill(fullEvent.mK, weight*fullEvent.weight);

	//xDistrib->Fill(x, weight*fullEvent.weight);
	return 0;
}
