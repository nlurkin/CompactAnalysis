/**************************************************************/
/* COmPACT user routine: user_superMcEvent(superMcEvent *evt) */
/*                                                            */
/* User routine called everytime an event `*evt' is           */
/* loaded. A return value of greater than zero denotes        */
/* an error condition has occured.                            */
/*                                   BH 13/2/98   RWM 11/7/97 */
/**************************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"  /* to get procedure calls for fortran routines */
#include "F77_ana.h" /* mapping of fortran common blocks for analysis routines */
#include "funLib.h"

int user_superMcEvent(superBurst *sbur,superMcEvent *evt) {
/* WARNING: do not alter things before this line */
	vector<double> vMass;
	vector<TVector3> vP;

	//Test MC
	if(!mcBranched){
		outTree->Branch("mc" ,"ROOTMCEvent", &rootMC);
		mcBranched = true;
	}
	rootMC.Clear();

	if(!opts["filter"].empty())
			if(!isFilteredEvent(sbur->nrun, sbur->time, evt->scmpevt.timeStamp, badEventsList)) return 0;

	if(optDebug){
		cout << "======SuperMCDecay======" << endl;
		cout << "KType: " << evt->decay.Ktype << "\t DType: " << evt->decay.Dtype << endl;
		cout << "Vertex: (" << evt->decay.dvertex[0] << "," << evt->decay.dvertex[1] << "," << evt->decay.dvertex[2] << ")" << endl;
		cout << "p_k: (" << evt->decay.p[0] << "," << evt->decay.p[1] << "," << evt->decay.p[2] << "," << evt->decay.p[3] << ")" << endl;
		cout << "======------------======" << endl;
		cout << "======  Particle  ======" << endl;
	}

	bool ep = false;
	bool em = false;
	bool k = false;
	bool pip = false;
	bool g = false;

	for(int i=0; i<= evt->Npart; i++){
		if(abs(evt->part[i].type)==512 && !k){
			if(evt->part[i].type==-512) rootMC.k.pdgID = -321;
			rootMC.k.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.k.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2]);
			rootMC.k.decay.SetXYZ(evt->part[i].dvertex[0], evt->part[i].dvertex[1], evt->part[i].dvertex[2]);
			k = true;
		}
		if(evt->part[i].type==64 && !ep){
			rootMC.ep.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.ep.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2]);
			rootMC.ep.decay.SetXYZ(evt->part[i].dvertex[0], evt->part[i].dvertex[1], evt->part[i].dvertex[2]);
			ep = true;
		}
		if(evt->part[i].type==-64 && !em){
			rootMC.em.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.em.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2]);
			rootMC.em.decay.SetXYZ(evt->part[i].dvertex[0], evt->part[i].dvertex[1], evt->part[i].dvertex[2]);
			em = true;
		}
		if(abs(evt->part[i].type)==8 && !pip){
			rootMC.pip.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.pip.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2]);
			rootMC.pip.decay.SetXYZ(evt->part[i].dvertex[0], evt->part[i].dvertex[1], evt->part[i].dvertex[2]);
			pip = true;
		}
		if(evt->part[i].type==16 && !g){
			if(evt->part[i].type==-16) rootMC.k.pdgID = -16;
			rootMC.gamma.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.gamma.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2]);
			rootMC.gamma.decay.SetXYZ(evt->part[i].dvertex[0], evt->part[i].dvertex[1], evt->part[i].dvertex[2]);
			g = true;
		}
		if(optDebug){
			cout << " ###### " << endl;
			cout << "\tType: " << evt->part[i].type << endl;
			cout << "\tp_i:   (" << evt->part[i].p[0] << "," << evt->part[i].p[1] << "," << evt->part[i].p[2] << "," << evt->part[i].p[3] << ")" << endl;
			cout << "\tProd:  (" << evt->part[i].pvertex[0] << "," << evt->part[i].pvertex[1] << "," << evt->part[i].pvertex[2] << ")" << endl;
			cout << "\tDecay: (" << evt->part[i].dvertex[0] << "," << evt->part[i].dvertex[1] << "," << evt->part[i].dvertex[2] << ")" << endl;
			cout << " ##### " << endl;
		}
	}
	if(optDebug)
		cout << "======------------======" << endl;

	rootMC.xTrue = pow((rootMC.em.P+rootMC.ep.P).M()/Mpi0, 2.);
	evt->scmpevt.DETstatus[0].LV3ABTrig = 1;



	rootBurst.isData = false;
	rootBurst.isMC = true;
	user_superCmpEvent(sbur, &evt->scmpevt);

/*---------- Add user C code here ----------*/
  /*static int nuserevt=0;

   if(nuserevt<20)
    {
      printSuperMcEvent(evt,fprt);
      nuserevt++;
    }
*/

/*----------- End of user C code -----------*/
  return 0;
}
