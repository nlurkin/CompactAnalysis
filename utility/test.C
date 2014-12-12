
void create(){
	TFile *fd = TFile::Open("test.root", "RECREATE");

	TTree *tree = new TTree("event", "event");
	int *eventBrch = new pi0dEvent();
	tree->Branch("toto", "pi0dEvent", &eventBrch);

	eventBrch->mK = 1;
	tree->Fill();
	eventBrch->mK = 2;
	tree->Fill();
	eventBrch->mK = 3;
	tree->Fill();

	fd->Write();
}

void read(){
	TFile *fd = TFile::Open("test.root");
	TTree *tree = (TTree*)fd->Get("event");
	pi0dEvent *event = 0;

	tree->SetBranchAddress("toto", &event);
	event = new pi0dEvent();
	tree->Show(0);
}

int test(){
	gSystem->Load("../obj/libmystructs.so");

	read();
}
