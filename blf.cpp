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
#include <stdexcept>

using namespace std;

//定义计算函数模板
template <typename T>
T compute(T n, T s0, T b, T sm) {
                //double s = n - b;	
                //double q = -2.0*(n*log((s0+b)/(sm+b))-s0-b+sm+b);
	        if (sm > s0){
			return 0;}
		else{
	                //return n*log(sm+b)-sm-b-lgamma(n+1);//似然
			//cout << s0 << endl;
			//cout << b << endl;
			//cout << sm << endl;
			//cout << -2.0*(n*log((s0+b)/(sm+b))-s0-b+sm+b) << endl;
			return -2.0*(n*log((s0+b)/(sm+b))-s0-b+sm+b);//似然比
                }}

//生成单个泊松随机数
int gen_poisson_single(double lambda) {
	static mt19937 engine(random_device{}());
	poisson_distribution<int> dist(lambda);
	return dist(engine);}

int blf() {

	for (int i = 0; i < 5000; i+=1) {
        

		double p = 0;
		double p_all = 0;
		double u = 4e-4;

		//以输入文件的数字为平均值，生成泊松随机数
		//读取文件
		vector<double> input_numbers;
		ifstream in ("b01.txt");
		ifstream in01 ("s01.txt");

		if (!in.is_open()) {
			cerr << "WRONG:Can't open the input file!" << endl;
			return 1;}

                string line, line01;
		while (getline(in, line) && getline(in01, line01)) {
			try {
				stringstream ss0(line);
				stringstream ss01(line01);
			        double num, num01;	
				ss0 >> num;
				ss01 >> num01;
				input_numbers.push_back(num+num01*253328*u/(1-u)/(253328.0*0.007682/15.108));
			        //input_numbers.push_back(num+num01);
			} catch ( const invalid_argument& e) {
				cerr << "WARING:Skip the invalid line '" << line << "'" << endl;
			}
		}
		in.close();

		//生成随机数并写入临时文件
		ofstream out("temp.txt");
		if (!out.is_open()) {
			cerr << "WRONG:Can't open the output file!" << endl;
			return 2;}

		for (int num001 : input_numbers) {
			if (num001 < 0) {
				cerr << "WARING:Skip the negative number:"<< num001  << endl;
				continue;
			}
		        //for(int j=0; j< 1000; j++){
                               p = gen_poisson_single(num001);
                          //     if(p !=0){
                            //                break;}}
			p_all = p_all + p;
			out << p << endl;
		}
		out.close();

		vector<double> results;

		//if (p_all-253328.0<253328.0*u/(1-u)) {
		//if (p_all-253328.0>0) {
		//compute the value of q
	        //打开数据文件
                ifstream file1("b01.txt");
	        ifstream file2("s01.txt");
	        ifstream file3("temp.txt");
                if(!file1.is_open() || !file2.is_open()) {
                       cerr << "Wrong:Can't open the file!" << endl;
                       return 1;
                }

                //逐行读取并处理数据
                string line1, line2, line3;
                int line_number = 0;
                while (getline(file1, line1) && getline(file2, line2) && getline(file3, line3)) {
                       line_number++;
                       stringstream ss1(line1), ss2(line2), ss3(line3);
                       double num1, num2, num3, num11, num21, num110, num211;
                       if(!(ss1 >> num1) || !(ss2 >> num2) || !(ss3 >> num3)) {
                               cerr << "line" << line_number << ":type wrong!" << endl;
                               return 1;}
                       if (ss1.peek() != EOF || ss2.peek() != EOF) {
                               cerr << "line" << line_number << "extra data!" << endl;
                               return 1;}
		       if ( num1 != 0) {
			       //double x = 4e-4;
			       double y = 1;
			       //double coeff_b = 424.0*0.919*9.31*1782.0*(1-x)/253328.0;
			       //double coeff_s = 424.0*0.919*9.31*1782.0*x/(253328.0*0.007682/15.108);
                               double coeff_s = 253328*u/(1-u)/(253328.0*0.007682/15.108);
			       num21 = coeff_s*num2;
			       num211 = (p_all-253328.0)*y*num2/(253328.0*0.007682/15.108); 
			       //num211 = p_all*y*(num2+num1)/(253328.0+253328.0*0.007682/15.108)-num1;  
			      //num211 = p_all*num1/253328.0-num1;  
                               results.push_back(compute(num3, num21, num1, num211));}
	}  //}
        //else {
        //              results.push_back(0);
	//}
	
	                //输出结果
	                double q = 0;
	                for (size_t i = 0; i < results.size(); ++i){
		               q += results[i];
	                       //cout << results[i] << endl;
	                }
                        //if (q <0){
			//	cout << 0 <<endl;}
			//else {
	                        cout << fixed << q << endl;
				//}
	                //if (q >= 34.4){
	                //      qn+=1;}
	                //file1.close();
	                //file2.close();
	                //file3.close();
		}
	//cout << qn << endl;
	
        //清空临时文件，即存储泊松随机数的文件
	ofstream clear("temp.txt", ios::trunc);
	clear.close();

        return 0;
}
