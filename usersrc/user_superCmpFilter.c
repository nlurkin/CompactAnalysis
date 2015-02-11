/******************************************************************/
/* COmPACT user routine: user_superCmpFilter(superCmpEvent *sevt) */
/*                                                                */
/* This routine acts as a filter for all superCmpEvent data.      */
/* If it returns a value <0 the event will be sent to             */
/* an output file if one has been opened.                         */
/*                                                    BH 10/2/98  */
/******************************************************************/

#include "cmpio.h"
#include "reader.h"
#include "user.h"

int user_superCmpFilter(superBurst *sbur,superCmpEvent *sevt,hyperCmpEvent *hevt) {
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	/* Remember to fill hyperCOmPACT structure if you want hcmp output */

	/*----------- End of user C code -----------*/
	/*ask trigger + 3 tracks and 2 with eop>0.5*/
	return -1;//nico_pi0DalitzFilter(sbur, sevt);
}
