double getMax(double a, double b){
	if(a>=b) return a;
	else return b;
}

void simpleCombine(){
	gSystem->Load("../obj/libmystructs.so");
	TString MCFile = "/afs/cern.ch/user/n/nlurkin/work/MCPi3.root";
	TString dataFile = "/afs/cern.ch/user/n/nlurkin/work/Datap3.root";
	int max = 100000000000;
	TH1D *mcTotP = new TH1D("mcTotP", "mcTotP", 100, 0, 100);
	TH1D *mcMK = new TH1D("mcMK", "mcMK", 100, 0.4, 0.52);
	TH1D *mcPipE = new TH1D("mcPipE", "mcPipE", 60, 0, 60);

	TH1D *dataTotP = new TH1D("dataTotP", "dataTotP", 100, 0, 100);
	TH1D *dataMK = new TH1D("dataMK", "dataMK", 100, 0.4, 0.52);
	TH1D *dataPipE = new TH1D("dataPipE", "dataPipE", 60, 0, 60);

	pi0dEvent *eventBrch = new pi0dEvent();

	TFile *fd = TFile::Open(MCFile, "r");
	TTree *treeMC = (TTree*)fd->Get("event");
	treeMC->SetBranchAddress("pi0dEvent", &eventBrch);

	for(int i=0; i<treeMC->GetEntries() && i<max; i++){
		if(i % 1000 == 0){
			cout << i << "/" << treeMC->GetEntries() << "\r";
			cout.flush();
		}
		treeMC->GetEntry(i);

		mcTotP->Fill(eventBrch->pTotal.Mag());
		mcMK->Fill(eventBrch->mK);
		mcPipE->Fill(eventBrch->pip->clusterEnergy);
	}

	fd->Close();
	cout << endl;

	TFile *fd = TFile::Open(dataFile, "r");
	TTree *treeData = (TTree*)fd->Get("event");
	treeData->SetBranchAddress("pi0dEvent", &eventBrch);

	for(int i=0; i<treeData->GetEntries() && i<max; i++){
		if(i % 1000 ==0){
			cout << i << "/" << treeData->GetEntries() << "\r";
			cout.flush();
		}
		treeData->GetEntry(i);

		dataTotP->Fill(eventBrch->pTotal.Mag());
		dataMK->Fill(eventBrch->mK);
		dataPipE->Fill(eventBrch->pip->clusterEnergy);
	}

	THStack *hStack = new THStack("pTotal", "pTotal");

	//Scale MC to Data
	double totalMC = 0;
	totalMC += mcTotP->Integral();
	double factor = ((double)(dataTotP->Integral()))/totalMC;

	mcTotP->Scale(factor);
	mcMK->Scale(factor);
	mcPipE->Scale(factor);

	//Style data
	mcTotP->SetLineColor(kRed);

	TH1D* r = new TH1D("ratio_TotP", "ratio_TotP", 100, 0, 100);
	r->Sumw2();
	r->Divide(dataTotP, mcTotP, 1, 1, "B");
	r->SetMarkerColor(kRed);
	r->SetMaximum(r->GetMaximum()*1.1);
	r->SetMinimum(r->GetMinimum()*0.9);

	TCanvas *c1 = new TCanvas("cTotP", "TotP");
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetGrid();
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	mcTotP->Draw("HIST");               // Draw h1
	dataTotP->Draw("SAME E P");         // Draw h2 on top of h1
	mcTotP->GetYaxis()->SetRangeUser(0,1.2*getMax(dataTotP->GetMaximum(), mcTotP->GetMaximum()));

	c1->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();
	r->SetStats(0);      // No statistics on lower plot
	r->Draw("ep");
	r->GetYaxis()->SetLabelSize(0.05);
	r->SetMarkerColor(kRed);

	//Style data
	mcMK->SetLineColor(kRed);

	TH1D* r1 = new TH1D("ratio_MK", "ratio_MK", 100, 0.4, 0.52);
	r1->Sumw2();
	r1->Divide(dataMK, mcMK, 1, 1, "B");
	r1->SetMarkerColor(kRed);
	r1->SetMaximum(r1->GetMaximum()*1.1);
	r1->SetMinimum(r1->GetMinimum()*0.9);

	TCanvas *c1 = new TCanvas("cMK", "MK");
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetGrid();
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	mcMK->Draw("HIST");               // Draw h1
	dataMK->Draw("SAME E P");         // Draw h2 on top of h1
	mcMK->GetYaxis()->SetRangeUser(0,1.2*getMax(dataMK->GetMaximum(), mcMK->GetMaximum()));

	c1->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();
	r1->SetStats(0);      // No statistics on lower plot
	r1->Draw("ep");
	r1->GetYaxis()->SetLabelSize(0.05);
	r1->SetMarkerColor(kRed);

	//Style data
	mcPipE->SetLineColor(kRed);

	TH1D* r2 = new TH1D("ratio_PipE", "ratio_PipE", 60, 0, 60);
	r2->Sumw2();
	r2->Divide(dataPipE, mcPipE, 1, 1, "B");
	r2->SetMarkerColor(kRed);
	r2->SetMaximum(r2->GetMaximum()*1.1);
	r2->SetMinimum(r2->GetMinimum()*0.9);

	TCanvas *c1 = new TCanvas("cPipE", "PipE");
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetGrid();
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	mcPipE->Draw("HIST");               // Draw h1
	dataPipE->Draw("SAME E P");         // Draw h2 on top of h1
	mcPipE->GetYaxis()->SetRangeUser(0,1.2*getMax(dataPipE->GetMaximum(), mcPipE->GetMaximum()));

	c1->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid(); // vertical grid
	pad2->Draw();
	pad2->cd();
	r2->SetStats(0);      // No statistics on lower plot
	r2->Draw("ep");
	r2->GetYaxis()->SetLabelSize(0.05);
	r2->SetMarkerColor(kRed);
}
