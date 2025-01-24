//##########################################################
// Topic: Rebuild tao+ and tao-
// (e+ e- --- tao+ tao- --- pi+ pi- nutao nutao)
// (The momentum of electrins in a straight line is needed.)
// Author: Yirui Liu
// Data: 2025.1.20
// #########################################################

#include <iostream>
#include <cmath>
#include <Math/Vector4D.h>
#include <map>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

int rebuildtau(){

  //E11, E12, E31, p31x, p31z,  p31z, E32, p32x, p32y,  p32z are needed.
  double E11 = 0;//e+
  double p11x = 0;
  double p11y = 0;
  double p11z = 0;
  double E12 = 0;//e-
  double p12x = 0;
  double p12y = 0;
  double p12z = 0;
  double E21 = 0;//tau+
  double E211 = 0;
  double p21 = 0;
  double m21 = 1.777;
  double p21x1 = 0;
  double p21x2 = 0;
  double p21y1 = 0;
  double p21y2 = 0;
  double p21z1 = 0;
  double p21z2 = 0;
  double p21x = 0;
  double p21y = 0;
  double p21z = 0;
  double p211x = 0;
  double p211y = 0;
  double p211z = 0;
  double E22 = 0;//tao-
  double E221 = 0;
  double p22 = 0;
  double m22= 1.777;
  double p22x1 = 0;
  double p22x2 = 0;
  double p22y1 = 0;
  double p22y2 = 0;
  double p22z1 = 0;
  double p22z2 = 0;
  double p22x = 0;
  double p22y = 0;
  double p22z = 0;
  double p221x = 0;
  double p221y = 0;
  double p221z= 0;
  double E31 = 0;//pi+
  double p31 = 0;
  double p31x = 0;
  double p31y = 0;
  double p31z = 0;
  double E32 = 0;//pi-
  double p32 = 0;
  double p32x = 0;
  double p32y = 0;
  double p32z = 0;
  double p41 = 0;//nu tau(+)
  double m41 = 0;
  double p42 = 0;//nu tao(-)
  double m42 = 0;



  //Obtain them from LHE.
  E11 = 2.3392666187;/*1.8423046603;2.3400000078;*///e+
  p11x = 0.0;
  p11y = 0.0;
  p11z = 2.3392665629;/*1.8423045895;2.3291061821;*/
  E12 = 2.3398366652;/*2.3324773522;2.3399999966;*///e-
  p12x = 0.0;
  p12y = 0.0;
  p12z = -2.3398366094;/*-2.3324772962;-2.3399999408;*/
  E211 = 2.3393746552;/*1.9968124544;2.3333242172;*///tau+
  p211x = -0.40747578281;/*0.48719567922;1.1619564199;*/
  p211y = -0.22104076697;/*0.56134465058;0.81543415519;*/
  p211z = 1.4491551184;/*0.52636907695;0.52115297180;*/
  E221 = 2.3397286287;/*2.1779695581;2.3357820175;*///tao-
  p221x = 0.40747578281;/*-0.48719567922;-1.1619564199;*/
  p221y = 0.22104076697;/*-0.56134465058;-0.81543415519;*/
  p221z = -1.4497251649;/*-1.0165417837;-0.53204673052;*/
  E31 = 0.66683407219;/*1.2671709430;1.3362044263;*///pi+
  p31x = 0.63300752413;/*0.003691537754;1.0521962699;*/
  p31y = -0.028059421976;/*1.1011448981;-0.13549998777;*/
  p31z = 0.1539544814;/*0.61131713856;0.80030318695;*/
  p31 = sqrt(pow(p31x, 2) + pow(p31y, 2) + pow(p31z, 2));
  E32 = 1.003537174;/*1.0559859264;1.6904112985;*///pi-
  p32x = 0.027112355107;/*-0.95681000885;-1.2986792091;*/
  p32y = 0.91950677434;/*0.22027291060;-1.0721477634;*/
  p32z = -0.37600429658;/*-0.36279593145;0.044068160123;*/
  p32 = sqrt(pow(p32x, 2) + pow(p32y, 2) + pow(p32z, 2));
  
  //Coefficients of the equation.
  double a11 = 0;
  double a12 = 0;
  double b11 = 0;
  double b12 = 0;
  double a = 0;
  double b = 0;
  double c = 0;

  //The reference frame is selected in the momentum center frame of e+ and e-.
  //Obtain E21 and the size of p.
  E21 = (pow(E11 + E12, 2) + pow (m21, 2) - pow(m22, 2)) / 2 /(E11 +E12);
  E22 = (pow(E11 + E12, 2) - pow (m21, 2) + pow(m22, 2)) / 2 /(E11 +E12);
  p21 = sqrt((pow((E11 + E12), 4) + pow((pow(m21, 2) - pow(m22, 2)), 2) - 2 * (pow(m21, 2) + pow(m22, 2)) * pow((E11 + E12) ,2)) / 4 / pow(E11 + E12 ,2));
  p22 = p21;
  p41 = sqrt(pow((E21 - E31),2) - pow(m41 ,2));
  p42 = sqrt(pow((E22 - E32),2) - pow(m42 ,2));
 
  //Obtain px, py, pz of tao+.
  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32z + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31z) / (p31y * p32z - p32y * p31z);
  a12 = (-p31x * p32z + p32x * p31z) / (p31y * p32z - p32y * p31z);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32y + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31y) / (p31z * p32y - p32z * p31y);
  b12 = (-p31x * p32y + p32x * p31y) / (p31z * p32y - p32z * p31y);
  a = 1 + pow(a12, 2) + pow(b12, 2);
  b = 2 * a12 * a11 + 2 * b12 * b11;
  c = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);
  p21x1 = (- b + sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
  p21x2 = (- b - sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
  //cout << p21x1 <<endl;
  //cout << p21x2 <<endl;
  
  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32z + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31z) / (p31x * p32z - p32x * p31z);
  a12 = (-p31y * p32z + p32y * p31z) / (p31x * p32z - p32x * p31z);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32x + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31x) / (p31z * p32x - p32z * p31x);
  b12 = (-p31y * p32x + p32y * p31x) / (p31z * p32x - p32z * p31x);
  a = 1 + pow(a12, 2) + pow(b12, 2);
  b = 2 * a12 * a11 + 2 * b12 * b11;
  c = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);
  p21y1 = (- b + sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
  p21y2 = (- b - sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
 
  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32x + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31x) / (p31y * p32x - p32y * p31x);
  a12 = (-p31z * p32x + p32z * p31x) / (p31y * p32x - p32y * p31x);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32y + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31y) / (p31x * p32y - p32x * p31y);
  b12 = (-p31z * p32y + p32z * p31y) / (p31x * p32y - p32x * p31y);
  a = 1 + pow(a12, 2) + pow(b12, 2);
  b = 2 * a12 * a11 + 2 * b12 * b11;
  c = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);
  p21z1 = (- b + sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
  p21z2 = (- b - sqrt( pow(b, 2) - 4 *a *c)) / 2 / a;
  
  //Find out the right px, py, pz.
  double diff1 = fabs(pow(p21x1, 2) + pow(p21y1, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x1,p21y1,p21z1
  double diff2 = fabs(pow(p21x1, 2) + pow(p21y1, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x1,p21y1,p21z2
  double diff3 = fabs(pow(p21x1, 2) + pow(p21y2, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x1,p21y2,p21z1
  double diff4 = fabs(pow(p21x1, 2) + pow(p21y2, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x1,p21y2,p21z2
  double diff5 = fabs(pow(p21x2, 2) + pow(p21y1, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x2,p21y1,p21z1
  double diff6 = fabs(pow(p21x2, 2) + pow(p21y1, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x2,p21y1,p21z2
  double diff7 = fabs(pow(p21x2, 2) + pow(p21y2, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x2,p21y2,p21z1
  double diff8 = fabs(pow(p21x2, 2) + pow(p21y2, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x2,p21y2,p21z2
  double eps = 1e-6;
  vector<double> PX(0);
  vector<double> PY(0);
  vector<double> PZ(0);
  vector<double> PXX(0);
  vector<double> PYY(0);
  vector<double> PZZ(0);
  if(diff1 < eps){
	  if(PX.size() == 0){
	  PX.push_back(p21x1);
	  PY.push_back(p21y1);
	  PZ.push_back(p21z1);}
	  else{PXX.push_back(p21x1);
	       PYY.push_back(p21y1);
	       PZZ.push_back(p21z1);}}
  if(diff2 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x1);
          PY.push_back(p21y1);
          PZ.push_back(p21z2);}
          else{PXX.push_back(p21x1);
               PYY.push_back(p21y1);
               PZZ.push_back(p21z2);}}
  if(diff3 < eps){
	  if(PX.size() == 0){
          PX.push_back(p21x1);
          PY.push_back(p21y2);
          PZ.push_back(p21z1);}
          else{PXX.push_back(p21x1);
               PYY.push_back(p21y2);
               PZZ.push_back(p21z1);}}
  if(diff4 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x1);
          PY.push_back(p21y2);
          PZ.push_back(p21z2);}
          else{PXX.push_back(p21x1);
               PYY.push_back(p21y2);
               PZZ.push_back(p21z2);}}
  if(diff5 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x2);
          PY.push_back(p21y1);
          PZ.push_back(p21z1);}
          else{PXX.push_back(p21x2);
               PYY.push_back(p21y1);
               PZZ.push_back(p21z1);}}
  if(diff6 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x2);
          PY.push_back(p21y1);
          PZ.push_back(p21z2);}
          else{PXX.push_back(p21x2);
               PYY.push_back(p21y1);
               PZZ.push_back(p21z2);}}
  if(diff7 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x2);
          PY.push_back(p21y2);
          PZ.push_back(p21z1);}
          else{PXX.push_back(p21x2);
               PYY.push_back(p21y2);
               PZZ.push_back(p21z1);}}
  if(diff8 < eps){
          if(PX.size() == 0){
          PX.push_back(p21x2);
          PY.push_back(p21y2);
          PZ.push_back(p21z2);}
          else{PXX.push_back(p21x2);
               PYY.push_back(p21y2);
               PZZ.push_back(p21z2);}}
  if(PX.size() == 1 && PXX.size() == 1){
         if(fabs(PX[0] - PY[0]) < fabs(PXX[0] - PYY[0])){
		p21x = PX[0];
		p21y = PY[0];
		p21z = PZ[0];}
	    else{
		    p21x = PXX[0];
                    p21y = PYY[0];
                    p21z = PZZ[0];}

  //Lorentz transformation
  TLorentzVector L0;
  L0.SetPxPyPzE(p11x+p12x, p11y+p12y, p11z+p12z, E11+E12);
  TLorentzVector L1;
  L1.SetPxPyPzE(p211x, p211y, p211z, E211);
  L1.Boost(-L0.BoostVector());
  TLorentzVector L2;
  L2.SetPxPyPzE(p221x, p221y, p221z, E221);
  L2.Boost(-L0.BoostVector());
  
  //Output the data.
  //tau+
  std::cout << "tao+" << std::endl;
  //std::cout << "px(real)(GeV):   " << L1.Px() << std::endl;
  //std::cout << "px(reconstructed)(GeV):   " << p21x << std::endl;
  std::cout << "Error of px(GeV):   " << (p21x -L1.Px()) / L1.Px() *100 << "%" << std::endl;
  //std::cout << "py(real)(GeV):   " << L1.Py() << std::endl;
  //std::cout << "py(reconstructed)(GeV):   " << p21y << std::endl;
  std::cout << "Error of py(GeV):   " << (p21y -L1.Py()) / L1.Py() *100 << "%" <<  std::endl;
  //std::cout << "pz(real)(GeV):   " << L1.Pz() << std::endl;
  //std::cout << "pz(reconstructed)(GeV):   " << p21z << std::endl;
  std::cout << "Error of pz(GeV):   " << (p21z -L1.Pz()) / L1.Pz() *100 << "%" << std::endl;
  //std::cout << "E(real)(GeV):   " << L1.E() << std::endl;
  //std::cout << "E(reconstructed)(GeV):   " << E21 << std::endl;
  std::cout << "Error of E(GeV):   " << (E21 -L1.E()) / L1.E() *100 << "%" << std::endl;
  
  //tau-
  std::cout << "tao-" << std::endl;
  //std::cout << "px(real)(GeV):   " << L2.Px() << std::endl;
  //std::cout << "px(reconstructed)(GeV):   " << -p21x << std::endl;
  std::cout << "Error of px(GeV):   " << (-p21x -L2.Px()) / L2.Px() *100 << "%" << std::endl;
  //std::cout << "py(real)(GeV):   " << L2.Py() << std::endl;
  //std::cout << "py(reconstructed)(GeV):   " << -p21y << std::endl;
  std::cout << "Error of py(GeV):   " << (-p21y -L2.Py()) / L2.Py() *100 << "%" <<  std::endl;
  //std::cout << "pz(real)(GeV):   " << L2.Pz() << std::endl;
  //std::cout << "pz(reconstructed)(GeV):   " << -p21z << std::endl;
  std::cout << "Error of pz(GeV):   " << (-p21z -L2.Pz()) / L2.Pz() *100 << "%" << std::endl;
  //std::cout << "E(real)(GeV):   " << L2.E() << std::endl;
  //std::cout << "E(reconstructed)(GeV):   " << E22 << std::endl;
  std::cout << "Error of E(GeV):   " << (E22 - L2.E()) / L2.E() *100 << "%" << std::endl;
  }
     else{
     cout << "WRONG!" << endl;}
  return 0;
}
