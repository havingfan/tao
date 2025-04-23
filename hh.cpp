#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

//定义计算函数模板
double compute(double n, double s0, double b) {
                double s = n - b;	
                if (s > s0 ){
			return 0;}
		else{
	                return -2.0*(n*log((s0+b)/n)-s0-b+n);
		}}

//生成单个泊松随机数
int gen_poisson_single(double lambda) {
	static mt19937 engine(random_device{}());
	poisson_distribution<int> dist(lambda);
	return dist(engine);}

int hh() {

	for (int i = 0; i < 50000; i+=1) {
	      double num3 = gen_poisson_single(17);
	      cout << compute(num3, 7, 10) << endl;}
	
	return 0;
}
