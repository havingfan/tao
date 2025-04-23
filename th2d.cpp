#include <fstream>
#include <vector>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

int th2d() {

	int q = 0;

	//创建二维直方图（定义范围和分箱数）
	//TH2D *h = new TH2D("h", "The Decay Length Distribution of n3; m_n3(GeV); vtan3; Decay Length", nx, *xmin, *xmax, ny, *ymin, *ymax);
	TCanvas *c = new TCanvas("c", "DATA", 800, 600);
	TH2D *h = new TH2D("h", "ta- > pi+ pi- pi- vt(Electron Tag);Mass Fraction;Energy Fraction", 50, 0.0, 1.0, 50, 0.0, 1.0);

	//填充数据
        ifstream xfile("smm.txt"), yfile("sme.txt");
        string linex, liney;
        while (getline(xfile, linex) && getline(yfile, liney)) {
                       h->Fill(stod(linex),stod(liney));
	               q += 1;}

        //设置等高线数值
	//int nContours = 3;
	//h2->SetContour(nContours);
	//h2->SetContourLevel(0, 0.1000);
	//h2->SetContourLevel(1, 1.0000);

	//设置调色板和样式
	//gStyle->SetPalette(kRainBow);//使用彩虹色系
        //gStyle->SetOptStat(0);//关闭统计信息框
	h->SetContour(20);//设置颜色分级数

	//绘制
        h->Draw("COLZ");//COLZ表示颜色映射+显示颜色条

	//改变统计框位置
	gStyle->SetOptStat(1111);
        TPaveStats *stats = (TPaveStats*)h->FindObject("stats");
	if (stats) {
		stats->SetX1NDC(0.65);//左下角X坐标
                stats->SetY1NDC(0.15);//左下角Y坐标
		stats->SetX2NDC(0.85);//右上角X坐标
                stats->SetY2NDC(0.35);//右上角Y坐标
		c->Modified();
		c->Update();}

	//显示网格线
	//c->SetGridx();
	//c->SetGridy();
	//c->Update();
	
	//边界
	c->SetLeftMargin(0.1);
	c->SetRightMargin(0.1);
	c->SetTopMargin(0.1);
	c->SetBottomMargin(0.1);

	c->SaveAs("me.png");
	cout << q << endl;
       	return 0;	     
}
