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

	rootMC.Clear();

	bool ep = false;
	bool em = false;
	for(int i=0; i<= evt->Npart; i++){
		if(evt->part[i].type==64 && !ep){
			rootMC.ep.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.ep.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2])
			ep = true;
		}
		if(evt->part[i].type==-64 && !em){
			rootMC.em.P.SetXYZT(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3], evt->part[i].p[0]);
			rootMC.em.vertex.SetXYZ(evt->part[i].pvertex[0], evt->part[i].pvertex[1], evt->part[i].pvertex[2])
			em = true;
		}
	}
	rootMC.xTrue = pow((rootMC.em.P+rootMC.ep.P).M()/Mpi0, 2.);


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
