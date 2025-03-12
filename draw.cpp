#include <iostream>
#include <fstream>

using namespace std;

int draw(){

	TH1F *h[1];
    
	h[0] = new TH1F( "h0", "The mass distribution of HNL;m(GeV);Events", 100, 0, 0.1);
	h[0]->Sumw2();
	TTree tree("tree", "temp tree");
	tree.ReadFile("mass.txt", "x/F");
	tree.Draw("x >> h0", "683.768", "goff");

	//h[1] = new TH1F( "h1", "The E distribution of tau-;E(GeV);Events", 35, 1, 8);
        //h[1]->Sumw2();
	//tree.ReadFile("meebw.txt", "x/F");
        //tree.Draw("x >> h1", "", "goff");

	//h[2] = new TH1F( "h2", "The E distribution of tau-;E(GeV);Events", 35, 1, 8);
        //h[2]->Sumw2();
	//tree.ReadFile("hee.txt", "x/F");
        //tree.Draw("x >> h2", "", "goff");

        //归一化
        //h[0]->Scale(1.0 / h[0]->Integral());
        //h[1]->Scale(1.0 / h[1]->Integral());
	//h[2]->Scale(1.0 / h[2]->Integral());

	//h[0]->SetLineColor(kRed);
	//h[1]->SetLineColor(kBlue);
	//h[2]->SetLineColor(kGreen);
	//h[0]->SetLineWidth(2);
	//h[1]->SetLineWidth(2);
        //h[2]->SetLineWidth(2);

        TCanvas *c1 = new TCanvas("c1","Multi-Histogram", 800, 600);
	
	//关闭默认统计框
	//h[0]->SetStats(0);
	//h[1]->SetStats(0);
	//h[2]->SetStats(0);

	//h[0]->Draw("HIST E");
	//h[1]->Draw("HIST E SAME");
	//h[2]->Draw("HIST E SAME");
        h[0]->Draw("HIST");
        //h[1]->Draw("HIST SAME");
        //h[2]->Draw("HIST SAME");

        //TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        //leg->AddEntry(h[0], Form("%s: #mu=%.2f, #sigma=%.2f", "bwcutoff=15", h[0]->GetMean(), h[0]->GetRMS()), "l");
	//leg->AddEntry(h[1], Form("%s: #mu=%.2f, #sigma=%.2f", "bwcutoff=10^7", h[1]->GetMean(), h[1]->GetRMS()), "l");
	//leg->AddEntry(h[2], Form("%s: #mu=%.2f, #sigma=%.2f", "bwcutoff=10^7*", h[2]->GetMean(), h[2]->GetRMS()), "l");

	//leg->Draw();

        c1->Print("mass.png");
	
	return 0;

}
