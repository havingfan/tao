#include <iostream>
#include <fstream>

using namespace std;

int tgraph(){
      
      //创建数据点
      const int n = 6;
      double x[n] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      double y[n] = {2.725911e15, 2.089096e14, 3.849165e13, 6.429107e12, 1.726936e12, 6.105890e11};
      //double y_err[n] = {0.0001226, 0.000834, 0.0007503, 0.0002974, 2.377e-05, 1.571e-07};

      //创建TGraph对象
      TGraph *graph = new TGraph(n, x, y);
      //TGraphErrors *graph = new TGraphErrors(n, x, y, nullptr, y_err);

      //设置图形样式
      graph->SetTitle("The Life Distribution of n3;m_n3(GeV);Life(1/GeV)");
      //graph->SetLineColor(kBlue);
      graph->SetLineWidth(2);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1);
      //graph->SetMarkerColor(kBlue);
      
      //绘制误差棒样式
      //graph->SetFillColor(kBlue)
      //graph-.SetFillStyle(3001);

      //创建画布并绘制
      TCanvas *c1 = new TCanvas("c1", "Error Bars", 800, 600);
      c1->SetLogy();

      graph->Draw("APL");

      //添加误差棒显示
      //graph->Draw("E1 same");

      //保存图像
      c1->SaveAs("life.png");

      return 0;

}
