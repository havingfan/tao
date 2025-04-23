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
template <typename T>
T compute(T n, T s0, T b) {
                double s = n - b;	
                if (s > s0 ){
			return 0;}
		else{
			//cout << s0 << endl;
                        //cout << b << endl;
                        //cout << n << endl;
                        //cout << -2.0*(n*log((s0+b)/n)-s0-b+n) << endl;
	                return -2.0*(n*log((s0+b)/n)-s0-b+n);
		}}

//生成单个泊松随机数
int gen_poisson_single(double lambda) {
	static mt19937 engine(random_device{}());
	poisson_distribution<int> dist(lambda);
	return dist(engine);}

int blfnm() {

	int qn = 0;
	for (int i = 0; i < 50000; i+=1) {

	//打开两个数据文件
        ifstream file1("b01.txt");
	ifstream file2("s01.txt");
        if(!file1.is_open() || !file2.is_open()) {
                cerr << "Wrong:Can't open the file!" << endl;
               return 1;
        }

        //逐行读取并处理数据
        string line1, line2;
        int line_number = 0;
        vector<double> results;
        while (getline(file1, line1) && getline(file2, line2)) {
                line_number++;
                stringstream ss1(line1), ss2(line2);
                double num1, num2, num3, num11, num21, num110, num210;
                if(!(ss1 >> num1) || !(ss2 >> num2)) {
                        cerr << "line" << line_number << ":type wrong!" << endl;
                        return 1;}
                if (ss1.peek() != EOF || ss2.peek() != EOF) {
                        cerr << "line" << line_number << "extra data!" << endl;
                        return 1;}
		if ( num1 != 0) {
			double x = 5.46e-4;
			//double coeff_b = 424.0*0.919*9.31*1782.0*(1-x)/253328.0;
			//double coeff_s = 424.0*0.919*9.31*1782.0*x/(253328.0*0.007682/15.108);
                        double coeff_s = 253328*x/(1-x)/(253328.0*0.007682/15.108);
                        num21 = coeff_s*num2;
			for(int j=0; j< 1000; j++){
			      num3 = gen_poisson_single(num1+num21);
			      if(num3 !=0){
				      break;}
			}
                results.push_back(compute(num3, num21, num1));}}
	
	//输出结果
	double q = 0;
	for (size_t i = 0; i < results.size(); ++i){
		q += results[i];}
	cout << fixed << q << endl;
	//if (q >= 34.4){
	  //      qn+=1;}
	file1.close();
	file2.close();
	}
	cout << qn << endl;
	return 0;
}
