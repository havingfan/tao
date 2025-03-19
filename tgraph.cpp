#include <iostream>
#include <fstream>

using namespace std;

int tgraph(){
      
      //创建数据点
      const int n1 = 6;
      //double x1[n1] = {pow(0.3, 5), pow(0.5, 5), pow(0.7, 5), 1, pow(1.3, 5), pow(1.6, 5)};
      //double y1[n1] = {3.668499e-16, 4.78676e-15, 2.597966e-14, 1.555426e-13, 5.790601e-13, 1.637763e-12};
      double x1[n1] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      double y1[n1] = {0.0001188, 0.000576, 0.0005154, 0.0001855, 2.439e-05, 2.384e-07};
      double y_err[n] = {3.538e-07, 1.366e-06, 1.154e-06, 4.196e-07, 7.637e-08, 7.071e-10};

      const int n2 = 6;
      double x2[n2] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      double y2[n2] = {2.213e+05, 1.089e+04, 1570, 216.3, 53.42, 18.29};
      double y

      //const int n3 = 6;
      //double x3[n3] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      //double y3[n3] = {2.205e+07, 1.085e+06, 1.579e+05, 2.189e+04, 5356, 1830};

      //const int n4 = 6;
      //double x4[n4] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      //double y4[n4] = {2.211e+09, 1.1e+08, 1.568e+07, 2.171e+06, 5.358e+05, 1.827e+05};

      //const int n5 = 6;
      //double x5[n5] = {0.3, 0.5, 0.7, 1.0, 1.3, 1.6};
      //double y5[n5] = {2.207e+11, 1.088e+10, 1.572e+09, 2.183e+08, 5.339e+07, 1.828e+07};

      //创建TGraph对象
      TGraph *graph1 = new TGraph(n1, x1, y1);
      //TGraphErrors *graph = new TGraphErrors(n, x, y, nullptr, y_err);
      //TGraph *graph2 = new TGraph(n2, x2, y2);
      //TGraph *graph3 = new TGraph(n3, x3, y3);
      //TGraph *graph4 = new TGraph(n4, x4, y4);
      //TGraph *graph5 = new TGraph(n5, x5, y5);

      //设置图形样式
      graph1->SetTitle("The Width Distribution of N3;m_N3^5(GeV^5);Width(GeV)");
      //graph1->SetLineColor(kRed);
      //graph1->SetLineWidth(2);
      graph1->SetMarkerStyle(20);
      graph1->SetMarkerSize(1);
      //graph1->SetMarkerColor(kRed);
      //graph2->SetLineColor(kBlue);
      //graph2->SetLineWidth(2);
      //graph2->SetMarkerStyle(20);
      //graph2->SetMarkerSize(1);
      //graph2->SetMarkerColor(kBlue);
      //graph3->SetLineColor(kGreen);
      //graph3->SetLineWidth(2);
      //graph3->SetMarkerStyle(20);
      //graph3->SetMarkerSize(1);
      //graph3->SetMarkerColor(kGreen);
      //graph4->SetLineColor(kOrange);
      //graph4->SetLineWidth(2);
      //graph4->SetMarkerStyle(20);
      //graph4->SetMarkerSize(1);
      //graph4->SetMarkerColor(kOrange);
      //graph5->SetLineColor(kCyan);
      //graph5->SetLineWidth(2);
      //graph5->SetMarkerStyle(20);
      //graph5->SetMarkerSize(1);
      //graph5->SetMarkerColor(kCyan);
      
      //定义线性函数y=kx+b
      TF1 *fitFunc = new TF1("fitFunc", "[0]*x", x1[0], x1[n1-1]);
      fitFunc->SetParNames("k");
      fitFunc->SetLineColor(kRed);
      
      //执行拟合
      graph1->Fit(fitFunc, "R");

      //获取拟合参数
      double k = fitFunc->GetParameter(0);
      //double b = fitFunc->GetParameter(1);
      double k_err = fitFunc->GetParError(0);
      //double b_err = fitFunc->GetParError(1);

      //绘制误差棒样式
      //graph->SetFillColor(kBlue)
      //graph->SetFillStyle(3001);

      //创建TMultiGraph并添加子图
      //TMultiGraph *mg = new TMultiGraph();
      //mg->Add(graph1);
      //mg->Add(graph2);
      //mg->Add(graph3);
      //mg->Add(graph4);
      //mg->Add(graph5);
      //mg->SetTitle("The Decay Length Disstribution of N3;m_N3(GeV);Decay Length(m)");

      //创建画布并绘制
      TCanvas *c1 = new TCanvas("c1", "Error Bars", 800, 600);
      c1->SetLogx();
      c1->SetLogy();

      //graph1->Draw("APL");
      //mg->Draw("APL");
      graph1->Draw("AP");
      fitFunc->Draw("same");

      //添加拟合公式标签
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextColor(kRed);
      latex.DrawLatexNDC(0.2, 0.8, Form("y = (%.5fe-13 #pm %.5fe-17)x", k*1e+13, k_err*1e+17));

      //添加误差棒显示
      //graph->Draw("E1 same");
   
      //添加图例
      //TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      //leg->AddEntry(graph1, "VtaN3=3e-2", "l");
      //leg->AddEntry(graph2, "VtaN3=3e-3", "l");
      //leg->AddEntry(graph3, "VtaN3=3e-4", "l");
      //leg->AddEntry(graph4, "VtaN3=3e-5", "l");
      //leg->AddEntry(graph5, "VtaN3=3e-6", "l");
      //leg->Draw();

      //保存图像
      c1->SaveAs("decay.png");

      return 0;

}
