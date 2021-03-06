/***********************************************************/
/* COmPACT user routine: user_superBurst(superBurst *sbur) */
/*                                                         */
/* User routine called everytime a SuperCOmPACT burst      */
/* `*bur' is loaded.                                       */
/* A return value of greater than zero denotes an error    */
/* condition has occured.                                  */
/*                                             BH 2/3/98   */
/***********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "compactLib.h"

int user_superBurst(superBurst *sbur) {
	rootBurst.clear();
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	if(sbur->brtype==2){
		rootBurst.isMC = true;
		rootBurst.isData = false;
	}
	else if(sbur->brtype==1){
		rootBurst.isMC = false;
		rootBurst.isData = true;
	}
	sbur->BadB.Skip = 0; /* see user_superBurst.example.c to learn to use it */
	if(		rootBurst.isData && (
			sbur->BadB.Phys != 0 ||
			sbur->BadB.Dch != 0)
	) sbur->BadB.Skip = 1;
	//Exclude Physics and DCH bad bursts

	loadEOPData(sbur);
	defineBeamCharge(sbur);
	/*----------- End of user C code -----------*/
	return 0;
}
