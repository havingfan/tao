#include <iostream>
#include <fstream>

using namespace std;

int draw(){

	TH1F *h[6];
    
	h[0] = new TH1F( "h0", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
	//h[0]->Sumw2();
	TTree tree("tree", "temp tree");
	tree.ReadFile("mass0.3.txt", "x/F");
	tree.Draw("x >> h0", "", "goff");

	h[1] = new TH1F( "h1", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
        //h[1]->Sumw2();
	TTree tree1("tree", "temp tree");
	tree1.ReadFile("mass0.5.txt", "x/F");
        tree1.Draw("x >> h1", "", "goff");

	h[2] = new TH1F( "h2", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
        //h[2]->Sumw2();
	TTree tree2("tree", "temp tree");
	tree2.ReadFile("mass0.7.txt", "x/F");
        tree2.Draw("x >> h2", "", "goff");

        h[3] = new TH1F( "h3", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
        //h[3]->Sumw2();
	TTree tree3("tree", "temp tree");
        tree3.ReadFile("mass1.0.txt", "x/F");
        tree3.Draw("x >> h3", "", "goff");

	h[4] = new TH1F( "h4", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
        //h[4]->Sumw2();
	TTree tree4("tree", "temp tree");
        tree4.ReadFile("mass1.3.txt", "x/F");
        tree4.Draw("x >> h4", "", "goff");

	h[5] = new TH1F( "h5", "The Mass distribution of u d~;Mass Fraction;Entries per bin", 100, 0, 1);
        //h[5]->Sumw2();
	TTree tree5("tree", "temp tree");
        tree5.ReadFile("mass1.6.txt", "x/F");
        tree5.Draw("x >> h5", "", "goff");

        //归一化
        //h[0]->Scale(1.0 / h[0]->Integral());
        //h[1]->Scale(1.0 / h[1]->Integral());
	//h[2]->Scale(1.0 / h[2]->Integral());

	h[0]->SetLineColor(kCyan);
	h[1]->SetLineColor(kMagenta);
	h[2]->SetLineColor(kGray);
	h[3]->SetLineColor(kBlue);
        h[4]->SetLineColor(kGreen);
        h[5]->SetLineColor(kRed);
	h[0]->SetLineWidth(2);
	h[1]->SetLineWidth(2);
        h[2]->SetLineWidth(2);
        h[3]->SetLineWidth(2);
        h[4]->SetLineWidth(2);
        h[5]->SetLineWidth(2);

        //将直方图数值缩放
	//h[0]->Scale(1e-1);
	//h[1]->Scale(1e-1);
        //h[2]->Scale(1e-1);
	//h[3]->Scale(1e-1);
        //h[4]->Scale(1e-1);
        //h[5]->Scale(1e-1);

	//设置标签格式（显示整数）
	//h[0]->GetYaxis()->SetLabelFormat(%d);
        //h[0]->GetYaxis()->SetLabelFormat(%d);

        TCanvas *c1 = new TCanvas("c1","Multi-Histogram", 800, 600);
	
	//关闭默认统计框
	h[0]->SetStats(0);
	h[1]->SetStats(0);
	h[2]->SetStats(0);
        h[3]->SetStats(0);
        h[4]->SetStats(0);
        h[5]->SetStats(0);

	//h[0]->Draw("HIST E");
	//h[1]->Draw("HIST E SAME");
	//h[2]->Draw("HIST E SAME");
        h[5]->Draw("");
        h[1]->Draw("SAME");
        h[2]->Draw("SAME");
        h[3]->Draw("SAME");
        h[4]->Draw("SAME");
        h[0]->Draw("SAME");

	//标签
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        //leg->AddEntry(h[0], Form("%s: #mu=%.2f, #sigma=%.2f", "bwcutoff=15", h[0]->GetMean(), h[0]->GetRMS()), "l");
	leg->AddEntry(h[5], Form("m_N3=1.6GeV", h[5]->GetRMS()), "l");
	leg->AddEntry(h[4], Form("m_N3=1.3GeV", h[4]->GetRMS()), "l");
	leg->AddEntry(h[3], Form("m_N3=1.0GeV", h[3]->GetRMS()), "l");
        leg->AddEntry(h[2], Form("m_N3=0.7GeV", h[2]->GetRMS()), "l");
        leg->AddEntry(h[1], Form("m_N3=0.5GeV", h[1]->GetRMS()), "l");
        leg->AddEntry(h[0], Form("m_N3=0.3GeV", h[0]->GetRMS()), "l");

	leg->Draw();

        c1->Print("mass.png");
	
	return 0;

}
