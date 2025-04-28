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
#include <ctime>
#include <stdexcept>
using namespace std;

//定义计算函数模板
double compute(double n, double s0, double b, double sm) {
                //if (s > s0 ){
		//	return 0;}
		//else{
	                //return n*log(sm+b)-sm-b-lgamma(n+1);        
		        return -2.0*(n*log((s0+b)/(sm+b))-s0-b+sm+b);
                }//}

//生成单个泊松随机数
int gen_poisson_single(double lambda) {
	thread_local static mt19937 engine(random_device{}());
	poisson_distribution<int> dist(lambda);
	return dist(engine);}

//牛顿迭代法
//计算方程值a/(x+b)-c
double compute_f(const vector<double>& a, const vector<double>& b, const vector<double>& c, double x) {
        double sum = 0.0;
	for (size_t i = 0; i <a.size(); ++i) {
		double denom = x +b[i];
	        if (fabs(denom) < 1e-12) {
		        throw runtime_error("Division by zero in compute_f");
	        }
		sum += a[i]/denom-c[i];
	}
	return sum;
}

//计算方程导数值-a/(x+b)^2
double compute_df(const vector<double>& a, const vector<double>& b, double x) {
	double sum = 0.0;
	for (size_t i = 0; i < a.size(); ++i) {
		double denom = x + b[i];
	        if (fabs(denom) < 1e-12) {
		        throw runtime_error("Division by zero in compute_df");
                }
	        sum += -a[i]/(denom*denom);
	}
        return sum;
}

//牛顿迭代法
double newton_raphson(const vector<double>& a, const vector<double>& b, const vector<double>& c, double x0, double tol, int max_iter) {
	double x = x0;
	for (int i = 0; i < max_iter; ++i) {
			try{
				double fx = compute_f(a, b, c, x);
				if (fabs(fx) < tol) {
					return x;
				}
				double dfx = compute_df(a, b, x);
				if (fabs(dfx) < 1e-12) {
					throw runtime_error("Zero derivative in Newton-Raphson");
				}
				x = x - fx / dfx;
			} catch (const runtime_error& e) {
				cerr << "Error: " << e.what() << " at x = " << x << endl;
				return NAN;
			}
		}
		cerr << "Newton-Raphson did not converge within " << max_iter << "iterations." << endl;
		return NAN;
	}

int hh() {

	double u = 4e-4;
	double y = 0;
	for (int i = 0; i < 5000; i+=1) {
	      double num31 = gen_poisson_single(29051+253328*u/(1-u)*62.807624/128.81029);
	      double num32 = gen_poisson_single(224277+253328*u/(1-u)*66.002666/128.81029);
	      
	      //求y
              vector<double> a = {num31, num32};
	      vector<double> b = {29051/(62.807624/128.81029), 224277/(66.002666/128.81029)};
	      vector<double> c = {62.807624/128.81029, 66.002666/128.81029};
 
	     //初始值
	      double x0 = 150;
	      double tol = 1e-9;
	      int max_iter = 100;
	      double y = newton_raphson(a, b, c, x0, tol, max_iter);
	      if(!isnan(y)) {
		      bool valid = true;
		      for (size_t i = 0; i < b.size(); ++i) {
			      if (fabs(y + b[i]) < 1e-9) {
				      valid = false;
				      break;
			      }
		      }
		      if (valid) {
			     // cout << "Root found: y = " << y << endl;
		      } else {
			      cout << "Root causes a division by zero." << endl;}
	      } else {
		      cout << "Failed to find a root." << endl;
	      }
	      
	      cout << compute(num31, 253328*u/(1-u)*62.807624/128.81029, 29051, y*62.807624/128.81029)+
		     compute(num32, 253328*u/(1-u)*66.002666/128.81029, 224277, y*66.002666/128.81029) << endl;}
	
	return 0;
}
