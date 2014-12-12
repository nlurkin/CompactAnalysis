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

int user_superMcEvent(superBurst *sbur,superMcEvent *evt) {
/* WARNING: do not alter things before this line */
	vector<double> vMass;
	vector<TVector3> vP;

	bool ep = false;
	bool em = false;
	for(int i=0; i<= evt->Npart; i++){
		//cout << "\t" << evt->part[i].type << endl;
		//cout << evt->part[i].pvertex[2] << endl;
		if(evt->part[i].type==64 && !ep){
			fullEvent.epPTrue = TVector3(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3]).Unit();
			fullEvent.epETrue = evt->part[i].p[0];
			ep = true;
		}
		if(evt->part[i].type==-64 && !em){
			fullEvent.emPTrue = TVector3(evt->part[i].p[1], evt->part[i].p[2], evt->part[i].p[3]).Unit();
			fullEvent.emETrue = evt->part[i].p[0];
			em = true;
		}
	}
	vMass.push_back(Me);
	vP.push_back(fullEvent.epPTrue*fullEvent.epETrue);
	vMass.push_back(Me);
	vP.push_back(fullEvent.emPTrue*fullEvent.emETrue);

	fullEvent.meeTrue = sqrt(invMass2(vMass, vP));
	fullEvent.xTrue = pow(fullEvent.meeTrue/Mpi0, 2.);


	user_superCmpEvent(sbur, &evt->scmpevt);

	/*if(fullEvent.ep){
		cout << fullEvent.mee << " " << fullEvent.meeTrue << endl;
		cout << printVector3(fullEvent.ep->momentum) << " " << printVector3(fullEvent.epPTrue) << endl;
		cout << printVector3(fullEvent.em->momentum) << " " << printVector3(fullEvent.emPTrue) << endl;
	}*/
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
