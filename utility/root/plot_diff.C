#include <TDirectory.h>
#include <TArrow.h>
#include <TString.h>
#include <TLegend.h>

TFile *fdtmp;

void prepare_dflt() {
	TFile *fdflt = TFile::Open("/afs/cern.ch/user/n/nlurkin/work/MinCompactExtRange/Dflt/data.root", "READ");
	TTree *tdflt = (TTree*) fdflt->Get("event");
	fdtmp->cd();
	tdflt->Draw("pi0dEvent.x>>hdflt (100, 0, 1)", "", "goff");
	fdflt->Close();
	TFile *fdflt = TFile::Open("/afs/cern.ch/user/n/nlurkin/work/MinCompactExtRange/Dflt/pi.root", "READ");
	TTree *tdflt = (TTree*) fdflt->Get("event");
	fdtmp->cd();
	tdflt->Draw("pi0dEvent.x>>hdflt_pi (100, 0, 1)", "", "goff");
	fdflt->Close();
}

void do_it(TString file, TString output, bool alsoPi) {
	TFile *fextra = TFile::Open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "READ");
	TTree *textra = (TTree*) fextra->Get("event");
	fdtmp->cd();
	textra->Draw("pi0dEvent.x>>hextra(100, 0, 1)", "", "goff");
	fextra->Close();

	TH1D *hextra = (TH1D*) gDirectory->Get("hextra");
	TH1D *hdflt = (TH1D*) gDirectory->Get("hdflt");
	TH1D* n = new TH1D("diff", "", 100, 0, 1);
	n->Add(hextra, 1);
	n->Add(hdflt, -1);

	TH1D* n_pi;
	if(alsoPi){
		TFile *fextra = TFile::Open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/pi.root", "READ");
		TTree *textra = (TTree*) fextra->Get("event");
		fdtmp->cd();
		textra->Draw("pi0dEvent.x>>hextra_pi(100, 0, 1)", "", "goff");
		fextra->Close();

		TH1D *hextra_pi = (TH1D*) gDirectory->Get("hextra_pi");
		TH1D *hdflt_pi = (TH1D*) gDirectory->Get("hdflt_pi");
		n_pi = new TH1D("diff_pi", "", 100, 0, 1);
		n_pi->Add(hextra_pi, 1);
		n_pi->Add(hdflt_pi, -1);
	}

	TCanvas *c = new TCanvas("c", "c", 600, 600);
	c->SetRightMargin(0.05);
	c->SetLeftMargin(0.15);
	//n_pi->Draw();
	n->SetMarkerColor(kRed);
	n->SetLineColor(kRed);
	n->SetMarkerStyle(8);
	n->SetMarkerSize(0.8);
	n->SetFillColor(kRed);
	n->SetStats(false);
	if(alsoPi){
		n->Draw("EP");
		n_pi->Scale(n->Integral()/n_pi->Integral());
		n_pi->SetFillColor(TColor::GetColor(137,116,232));
		n_pi->SetLineColor(TColor::GetColor(137,116,232));
		n_pi->Draw("HIST SAME");
		TLegend *l = new TLegend(0.75, 0.75, 0.95, 0.9);
		l->SetEntrySeparation(0.01);
		l->SetTextSize(20);
		l->SetTextFont(43);
		l->AddEntry(n, "Data", "lp");
		l->AddEntry(n_pi, "K^{#pm}#rightarrow#pi^{#pm}#pi^{0}_{D}", "lf");
		l->Draw();
	}
	else
		n->Draw();

	n->GetYaxis()->SetTitle("Events/(0.01)");
	n->GetYaxis()->SetTitleSize(25);
	n->GetYaxis()->SetTitleFont(43);
	n->GetYaxis()->SetTitleOffset(1.70);

	// X axis ratio plot settings
	n->GetXaxis()->SetTitle("x");
	n->GetXaxis()->SetTitleSize(25);
	n->GetXaxis()->SetTitleColor(1);
	n->GetXaxis()->SetTitleFont(43);
	n->GetXaxis()->SetLabelFont(43);
	n->GetXaxis()->SetLabelSize(15);

	c->SaveAs(output + ".png");
	c->SaveAs(output + ".pdf");
}

void do_time_track(TString file, TString name) {
	fdtmp->cd();
	TH1D* htime = new TH1D("htime", "", 160, -80, 80);

	//for(int x =0 ; x<233; ++x){
	//TFile *ftime = TFile::Open(TString::Format("xroot://eosna62.cern.ch//eos/na62/user/n/nlurkin/%s/Data/job_data%i.root", file.Data(), x), "READ");
	TFile *ftime = TFile::Open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "READ");
	TTree *ttime = (TTree*) ftime->Get("event");
	fdtmp->cd();
	ttime->Draw("rawEvent.track.time-rawBurst.tOffst.Dch>>htimetmp(160, -80, 80)", "", "goff");
	ftime->Close();

	TH1D *htimetmp = (TH1D*) gDirectory->Get("htimetmp");
	htime->Add(htimetmp, 1);
	//}

	TCanvas *c = new TCanvas("c", "c", 600, 600);
	c->SetRightMargin(0.05);
	c->SetLeftMargin(0.15);
	c->SetLogy(true);
	htime->SetTitle("");
	htime->SetMarkerColor(kRed);
	htime->SetLineColor(kRed);
	htime->SetMarkerStyle(8);
	htime->SetMarkerSize(0.8);
	htime->SetFillColor(kRed);
	htime->SetStats(false);
	htime->Draw();

	htime->GetYaxis()->SetTitle("Events/(0.01)");
	htime->GetYaxis()->SetTitleSize(25);
	htime->GetYaxis()->SetTitleFont(43);
	htime->GetYaxis()->SetTitleOffset(1.70);
	htime->GetYaxis()->SetRangeUser(10, 5e5);

	// X axis ratio plot settings
	htime->GetXaxis()->SetTitle("t_{track}-TOffset_{DCH}");
	htime->GetXaxis()->SetTitleSize(25);
	htime->GetXaxis()->SetTitleColor(1);
	htime->GetXaxis()->SetTitleFont(43);
	htime->GetXaxis()->SetLabelFont(43);
	htime->GetXaxis()->SetLabelSize(15);

	c->SaveAs(name + ".png");
	c->SaveAs(name + ".pdf");
}

void do_time_vtx(TString file, TString name) {
	TFile *ftime = TFile::Open("xroot://castorpublic.cern.ch//castor/cern.ch/user/n/nlurkin/" + file + "/data.root", "READ");
	TTree *ttime = (TTree*) ftime->Get("event");
	fdtmp->cd();
	ttime->Draw("rawEvent.vtx.time-rawBurst.tOffst.Dch>>htime(160, -80, 80)", "", "goff");
	ftime->Close();

	TH1D *htime = (TH1D*) gDirectory->Get("htime");

	TCanvas *c = new TCanvas("c", "c", 600, 600);
	c->SetRightMargin(0.05);
	c->SetLeftMargin(0.15);
	c->SetLogy(true);
	htime->SetTitle("");
	htime->SetMarkerColor(kRed);
	htime->SetLineColor(kRed);
	htime->SetMarkerStyle(8);
	htime->SetMarkerSize(0.8);
	htime->SetFillColor(kRed);
	htime->SetStats(false);
	htime->Draw();
	TArrow *a1 = new TArrow(25, 15000, 25, 1000, 0.01, "|>");
	a1->Draw();
	TArrow *a2 = new TArrow(-25, 15000, -25, 1000, 0.01, "|>");
	a2->Draw();

	htime->GetYaxis()->SetTitle("Events/(0.01)");
	htime->GetYaxis()->SetTitleSize(25);
	htime->GetYaxis()->SetTitleFont(43);
	htime->GetYaxis()->SetTitleOffset(1.70);
	htime->GetYaxis()->SetRangeUser(0.5, 1e5);

	// X axis ratio plot settings
	htime->GetXaxis()->SetTitle("t_{vtx} - TOffset_{DCH}");
	htime->GetXaxis()->SetTitleSize(25);
	htime->GetXaxis()->SetTitleColor(1);
	htime->GetXaxis()->SetTitleFont(43);
	htime->GetXaxis()->SetLabelFont(43);
	htime->GetXaxis()->SetLabelSize(15);

	c->SaveAs(name + ".png");
	c->SaveAs(name + ".pdf");
}


void plot_diff() {
	fdtmp = TFile::Open("tmp.root", "RECREATE");
	prepare_dflt();

	do_it("Stage2_MinCompactExtRange_correct/extratrack", "x_extra_tracks", false);
	do_it("Stage2_MinCompactExtRange_correct/vtx_charge", "x_vertex_qvtx", false);
	do_it("Stage2_MinCompactExtRange_correct/track_time", "x_out_of_time", false);
	do_it("Stage2_MinCompactExtRange/rndpid", "two_hypo", true);
	do_time_track("Stage2_MinCompactExtRange_correct/track_time", "track_time");
	do_time_vtx("Stage2_MinCompactExtRange_correct/track_time", "vtx_time");
}
