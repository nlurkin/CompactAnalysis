
double alpha=1./137.;
double r=5.73E-5;
double deltaX=0;

double F(double x, double a){
	return (1+a*x);
}

double fun(double x, double a){
	return ((2*alpha)/(2*TMath::Pi()))*(pow(1-x,3)/x)*(1+pow(r,2)/(2*x))*sqrt(1-pow(r,2)/x)*pow(F(x,a), 2)*(1+deltaX);
}

distrib(){
	int bins=1000;
	double step = 1./(double)bins;
	double x;
	double a = 0.3;

	TH1D *d1 = new TH1D("d1", "Decay Rate", bins+1, 0-step/2., 1+step/2.);
	TH1D *d2 = new TH1D("d2", "Decay Rate", bins+1, 0-step/2., 1+step/2.);
	double start = 0.005/step;
	for(int i=start; i<bins; i++){
		x = step*i;
		d1->Fill(x, fun(x, 0.)*100000.);
		d2->Fill(x, fun(x, a)*100000.);
	}
	
	d1->SetLineColor(kRed);
	d2->SetLineColor(kBlue);
	d2->GetXaxis()->SetTitle("x");
	d1->Draw();
	d2->Draw("SAME");
}
