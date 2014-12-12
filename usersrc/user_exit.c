/********************************************************/
/* COmPACT user routine: user_exit()                    */
/*                                                      */
/* User routine called once all data has been read. This*/
/* allows output files to be closed and any other user  */
/* resources to be tidied up.                           */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "TFile.h"
#include <iostream>
#include "TH1D.h"

int user_exit() {
	/* WARNING: do not alter things before this line */
	/*---------- Add user C code here ----------*/
	cout << endl;
	if(outTree) outTree->Write();
	gFile->Write();
	gFile->Purge();
	gFile->Close();

	delete gFile;

	if(!noOutput){
		fclose(fprt);
		fclose(fprt2);
	}
	/*----------- End of user C code -----------*/
	return 0;
}
