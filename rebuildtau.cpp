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
  int m = 0;
  double E111 = 0;//e+
  double p111x = 0;
  double p111y = 0;
  double p111z = 0;
  double E121 = 0;//e-
  double p121x = 0;
  double p121y = 0;
  double p121z = 0;
  double E2111 = 0;//tau+
  double p2111x = 0;
  double p2111y = 0;
  double p2111z = 0;
  double p211 = 0;
  double E2211 = 0;//tao-
  double p2211x = 0;
  double p2211y = 0;
  double p2211z = 0;
  double p221 = 0;
  double E311 = 0;//pi+
  double p311x = 0;
  double p311y = 0;
  double p311z = 0;
  double p311 = 0;
  double E321 = 0;//pi-
  double p321x = 0;
  double p321y = 0;
  double p321z = 0;
  double p321 = 0;

  //Obtain them from LHE.
  int N;
  vector<double> *energy = 0;
  vector<double> *px = 0;
  vector<double> *py = 0;
  vector<double> *pz = 0;
  vector<int> *status = 0;
  vector<int> *pid = 0;
  vector<int> *mother1 = 0;
  vector<int> *mother2 = 0;
  int numParticles = 0;
  int n0=0;
  int n01=0;
  int n1=0;
  int n11=0;
  int n2=0;
  int n21=0;
  int n3 = 0;
  int n4 = 0;

  string dir = "/home/yrliu/Desktop/";
  string name = "lhe.root";
  string filename = dir + name;

  TBranch *b_energy;
  TBranch *b_pid;
  TBranch *b_px;
  TBranch *b_py;
  TBranch *b_pz;
  TBranch *b_numParticles;
  TBranch *b_status;
  TBranch *b_mother1;
  TBranch *b_mother2;

  TFile *f = new TFile( filename.c_str() );
  TTree *my_tree = (TTree*)f->Get("lhedata");

  my_tree->SetBranchAddress("energy", &energy, &b_energy );
  my_tree->SetBranchAddress("pid", &pid, &b_pid );
  my_tree->SetBranchAddress("px", &px, &b_px );
  my_tree->SetBranchAddress("py", &py, &b_py );
  my_tree->SetBranchAddress("pz", &pz, &b_pz );
  my_tree->SetBranchAddress("numParticles", &numParticles, &b_numParticles );
  my_tree->SetBranchAddress("status", &status, &b_status );
  my_tree->SetBranchAddress("mother1", &mother1, &b_mother1 );
  my_tree->SetBranchAddress("mother2", &mother2, &b_mother2 );

  //Get one event and analysis it.
  N = (int) my_tree->GetEntries();
  for(int i=0; i<N+1; i++){
        my_tree->GetEntry( i );

  //Determine whether or not the event is what we need.
  if(numParticles == 8){
        unordered_set<int> set1 = {pid->at(0), pid->at(1), pid->at(2), pid->at(3), pid->at(4), pid->at(5), pid->at(6), pid->at(7)};
        unordered_set<int> set2 = {11, -11, -15, -16, 211, 15, 16, -211};
        if(set1 ==set2){

              //Match data and particles one by one.
              map< int, int > Particles{{pid->at(0), 0}, {pid->at(1), 1}, {pid->at(2), 2}, {pid->at(3), 3}, {pid->at(4), 4}, {pid->at(5), 5}, {pid->at(6), 6}, {pid->at(7), 7}};
              n0 = Particles.find(-15)->second;
              n1 = Particles.find(-16)->second;
              n2 = Particles.find(211)->second;
              n01 = Particles.find(15)->second;
              n11 = Particles.find(16)->second;
              n21 = Particles.find(-211)->second;
              n3 = Particles.find(-11)->second;
              n4 = Particles.find(11)->second;

	      //Further screening
              //if( status->at(n0) == 2 && status->at(n1) == 1 && status->at(n2) == 1 && status->at(n3) == -1 &&
              //status->at(n01) == 2 && status->at(n11) == 1 && status->at(n21) == 1 && status->at(n4) == -1 &&
 	      //px->at(n3) == 0 && py->at(n3) == 0 && px->at(n4) == 0 && py->at(n4) == 0 &&
	      //mother1->at(n0) == 1 && mother2->at(n0) == 2 && mother1->at(n1) == 3 && mother2->at(n1) == 3 &&
	      //mother1->at(n2) == 3 && mother2->at(n2) == 3 && mother1->at(n01) == 1 && mother2->at(n01) == 2 &&
	      //mother1->at(n11) == 4 && mother2->at(n11) == 4 && mother1->at(n21) == 4 && mother2->at(n21) == 4 &&
	      //mother1->at(n3) == 0 && mother2->at(n3) == 0 && mother1->at(n4) == 0 && mother2->at(n4) ==0){
                    E111 = energy->at(n3);//e+
                    p111x = px->at(n3);
                    p111y = py->at(n3);
                    p111z = pz->at(n3);
                    E121 = energy->at(n4);//e-
                    p121x = px->at(n4);
                    p121y = py->at(n4);
                    p121z = pz->at(n4);
                    E2111 = energy->at(n0);//tau+
                    p2111x = px->at(n0);
                    p2111y = py->at(n0);
                    p2111z = pz->at(n0);
                    E2211 = energy->at(n01);//tao-
                    p2211x = px->at(n01);
                    p2211y = py->at(n01);
                    p2211z = pz->at(n01);
                    E311 = energy->at(n2);//pi+
                    p311x = px->at(n2);
                    p311y = py->at(n2);
                    p311z = pz->at(n2);
                    p311 = sqrt(pow(p31x, 2) + pow(p31y, 2) + pow(p31z, 2));
                    E321 = energy->at(n21);//pi-
                    p321x = px->at(n21);
                    p321y = py->at(n21);
                    p321z = pz->at(n21);
                    p321 = sqrt(pow(p32x, 2) + pow(p32y, 2) + pow(p32z, 2));
                    
		    //Lorentz transformation
                    TLorentzVector L0;
                    L0.SetPxPyPzE(p111x+p121x, p111y+p121y, p111z+p121z, E111+E121);
		    TLorentzVector L11;
                    L11.SetPxPyPzE(p111x, p111y, p111z, E111);
                    L11.Boost(-L0.BoostVector());
		    TLorentzVector L12;
                    L12.SetPxPyPzE(p121x, p121y, p121z, E121);
                    L12.Boost(-L0.BoostVector());
		    TLorentzVector L21;
                    L21.SetPxPyPzE(p2111x, p2111y, p2111z, E2111);
                    L21.Boost(-L0.BoostVector());
                    TLorentzVector L22;
                    L22.SetPxPyPzE(p2211x, p2211y, p2211z, E2211);
                    L22.Boost(-L0.BoostVector());
                    TLorentzVector L31;
                    L31.SetPxPyPzE(p311x, p311y, p311z, E311);
                    L31.Boost(-L0.BoostVector());
                    TLorentzVector L32;
                    L32.SetPxPyPzE(p321x, p321y, p321z, E321);
                    L32.Boost(-L0.BoostVector());

                    E11 = L11.E();//e+
                    p11x = L11.Px();
                    p11y = L11.Py();
                    p11z = L11.Pz();
                    E12 = L12.E();//e-
                    p12x = L12.Px();
                    p12y = L12.Py();
                    p12z = L12.Pz();
                    E211 = L21.E();//tau+
                    p211x = L21.Px();
                    p211y = L21.Py();
                    p211z = L21.Pz();
		    p211 = sqrt(pow(p211x, 2) + pow(p211y, 2) + pow(p211z, 2));
                    E221 = L22.E();//tao-
                    p221x = L22.Px();
                    p221y = L22.Py();
                    p221z = L22.Pz();
		    p221 = sqrt(pow(p221x, 2) + pow(p221y, 2) + pow(p221z, 2));
		    E31 = L31.E();//pi+
                    p31x = L31.Px();
                    p31y = L31.Py();
                    p31z = L31.Pz();
                    p31 = sqrt(pow(p31x, 2) + pow(p31y, 2) + pow(p31z, 2));
                    E32 = L31.E();//pi-
                    p32x = L32.Px();
                    p32y = L32.Py();
                    p32z = L32.Pz();
                    p32 = sqrt(pow(p32x, 2) + pow(p32y, 2) + pow(p32z, 2));

  //Coefficients of the equation.
  double a11 = 0;
  double a12 = 0;
  double b11 = 0;
  double b12 = 0;
  double ax = 0;
  double bx = 0;
  double cx = 0;
  double ay = 0;
  double by = 0;
  double cy = 0;
  double az = 0;
  double bz = 0;
  double cz = 0;
  double p21x1 = 0;
  double p21x2 = 0;
  double p21y1 = 0;
  double p21y2 = 0;
  double p21z1 = 0;
  double p21z2 = 0;

  //The reference frame is selected in the momentum center frame of e+ and e-.
  //Obtain E21 and the size of p.
  E21 = (E11 + E12) / 2;
  E22 = (E11 + E12) / 2;
  p21 = sqrt(pow(E21, 2) - pow(m21, 2));
  p22 = sqrt(pow(E22, 2) - pow(m22, 2));
  p41 = fabs(E21 - E31);
  p42 = fabs(E22 - E32);

  //Obtain px, py, pz of tao+.
  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32z + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31z) / (p31y * p32z - p32y * p31z);
  a12 = (-p31x * p32z + p32x * p31z) / (p31y * p32z - p32y * p31z);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32y + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31y) / (p31z * p32y - p32z * p31y);
  b12 = (-p31x * p32y + p32x * p31y) / (p31z * p32y - p32z * p31y);
  ax = 1 + pow(a12, 2) + pow(b12, 2);
  bx = 2 * a12 * a11 + 2 * b12 * b11;
  cx = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);

  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32z + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31z) / (p31x * p32z - p32x * p31z);
  a12 = (-p31y * p32z + p32y * p31z) / (p31x * p32z - p32x * p31z);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32x + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31x) / (p31z * p32x - p32z * p31x);
  b12 = (-p31y * p32x + p32y * p31x) / (p31z * p32x - p32z * p31x);
  ay = 1 + pow(a12, 2) + pow(b12, 2);
  by = 2 * a12 * a11 + 2 * b12 * b11;
  cy = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);

  a11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32x + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31x) / (p31y * p32x - p32y * p31x);
  a12 = (-p31z * p32x + p32z * p31x) / (p31y * p32x - p32y * p31x);
  b11 = ((pow(p31, 2) + pow(p21, 2) - pow(p41, 2)) /2 * p32y + (pow(p32, 2) + pow(p22, 2) - pow(p42, 2)) /2 * p31y) / (p31x * p32y - p32x * p31y);
  b12 = (-p31z * p32y + p32z * p31y) / (p31x * p32y - p32x * p31y);
  az = 1 + pow(a12, 2) + pow(b12, 2);
  bz = 2 * a12 * a11 + 2 * b12 * b11;
  cz = pow(a11, 2) + pow(b11, 2) - pow(p21, 2);
  if (pow(bz, 2) - 4 *az *cz < 0){
     p21x = - bx / 2 / ax;
     p21y = - by / 2 / ay;
     p21z = - bz / 2 / az;}
     //cout << "p21x1" << p21x << endl;
     //cout << "p21x2" << 0 << endl;
     //cout << "p21y1" << p21y << endl;
     //cout << "p21y2" << 0 << endl;
     //cout << "p21z1" << p21z << endl;
     //cout << "p21z2" << 0 << endl;}
  else{
     p21x1 = (- bx + sqrt( pow(bx, 2) - 4 *ax *cx)) / 2 / ax;
     p21x2 = (- bx - sqrt( pow(bx, 2) - 4 *ax *cx)) / 2 / ax;
     p21y1 = (- by + sqrt( pow(by, 2) - 4 *ay *cy)) / 2 / ay;
     p21y2 = (- by - sqrt( pow(by, 2) - 4 *ay *cy)) / 2 / ay;
     p21z1 = (- bz + sqrt( pow(bz, 2) - 4 *az *cz)) / 2 / az;
     p21z2 = (- bz - sqrt( pow(bz, 2) - 4 *az *cz)) / 2 / az;
     //cout << "p21x1" << p21x1 << endl;
     //cout << "p21x2" << p21x2 << endl;
     //cout << "p21y1" << p21y1 << endl;
     //cout << "p21y2" << p21y2 << endl;
     //cout << "p21z1" << p21z1 << endl;
     //cout << "p21z2" << p21z2 << endl;

  //Find out the right px, py, pz.
  double diff1 = fabs(pow(p21x1, 2) + pow(p21y1, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x1,p21y1,p21z1
  double diff2 = fabs(pow(p21x1, 2) + pow(p21y1, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x1,p21y1,p21z2
  double diff3 = fabs(pow(p21x1, 2) + pow(p21y2, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x1,p21y2,p21z1
  double diff4 = fabs(pow(p21x1, 2) + pow(p21y2, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x1,p21y2,p21z2
  double diff5 = fabs(pow(p21x2, 2) + pow(p21y1, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x2,p21y1,p21z1
  double diff6 = fabs(pow(p21x2, 2) + pow(p21y1, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x2,p21y1,p21z2
  double diff7 = fabs(pow(p21x2, 2) + pow(p21y2, 2) + pow(p21z1, 2) - pow(p21, 2));//p21x2,p21y2,p21z1
  double diff8 = fabs(pow(p21x2, 2) + pow(p21y2, 2) + pow(p21z2, 2) - pow(p21, 2));//p21x2,p21y2,p21z2
  
  //min diff and min |px-py|
  double min1 = 0;
  double min2 = 0;
  int nmin1 = 0;
  int nmin2 = 0;
  double minx1 = 0;
  double minx2 = 0;
  double miny1 = 0;
  double miny2 = 0;
  double minz1 = 0;
  double minz2 = 0;
  double MIN[8] = { diff1, diff2, diff3, diff4, diff5, diff6, diff7, diff8 };
  if ( MIN[0] > MIN[1] ){
	  min1 = MIN[1];
	  min2 = MIN[0];
          nmin1 = 2;
          nmin2 = 1;}
  else {
	  min1 = MIN[0];
	  min2 = MIN[1];
          nmin1 = 1;
          nmin2 = 2;}
  for ( int i=2; i<7; i++ ){
          if ( MIN[i]<min1 ){
		  min2 = min1;
		  min1 = MIN[i];
		  nmin2 = nmin1;
	          nmin1 = i + 1;}
	  else if ( MIN[i]<min2 ){
  		  min2 = MIN[i];
	          nmin2 = i + 1;}}
  if ( nmin1==1 ){ 
	  minx1 = p21x1;
	  miny1 = p21y1;
	  minz1 = p21z1;}
  if ( nmin1==2 ){ 
	  minx1 = p21x1;
	  miny1 = p21y1;
	  minz1 = p21z2;}
  if ( nmin1==3 ){ 
	  minx1 = p21x1;
	  miny1 = p21y2;
	  minz1 = p21z1;}
  if ( nmin1==4 ){
	  minx1 = p21x1;
	  miny1 = p21y2;
	  minz1 = p21z2;}
  if ( nmin1==5 ){ 
	  minx1 = p21x2;
	  miny1 = p21y1;
	  minz1 = p21z1;}
  if ( nmin1==6 ){
	  minx1 = p21x2;
	  miny1 = p21y1;
	  minz1 = p21z2;}
  if ( nmin1==7 ){ 
	  minx1 = p21x2;
	  miny1 = p21y2;
	  minz1 = p21z1;}
  if ( nmin1==8 ){ 
	  minx1 = p21x2;
	  miny1 = p21y2;
	  minz1 = p21z2;}
  if ( nmin2==1 ){ 
	  minx2 = p21x1;
	  miny2 = p21y1;
	  minz2 = p21z1;}
  if ( nmin2==2 ){ 
	  minx2 = p21x1;
	  miny2 = p21y1;
	  minz2 = p21z2;}
  if ( nmin2==3 ){ 
	  minx2 = p21x1;
	  miny2 = p21y2;
	  minz2 = p21z1;}
  if ( nmin2==4 ){ 
	  minx2 = p21x1;
	  miny2 = p21y2;
	  minz2 = p21z2;}
  if ( nmin2==5 ){ 
	  minx2 = p21x2;
	  miny2 = p21y1;
	  minz2 = p21z1;}
  if ( nmin2==6 ){ 
	  minx2 = p21x2;
	  miny2 = p21y1;
	  minz2 = p21z2;}
  if ( nmin2==7 ){ 
	  minx2 = p21x2;
	  miny2 = p21y2;
	  minz2 = p21z1;}
  if ( nmin2==8 ){ 
	  minx2 = p21x2;
	  miny2 = p21y2;
	  minz2 = p21z2;}
  //if(fabs(minx1 - miny1) < fabs(minx2 - minx2)){
    //            p21x = minx1;
      //          p21y = miny1;
        //        p21z = minz1;}
          //      else{
            //        p21x = minx2;
              //      p21y = miny2;
                //    p21z = minz2;} 
  p21x = (minx1+minx2)/2;
  p21y = (miny1+miny2)/2;
  p21z = (minz1+minz2)/2;

  //cout << "minx1" << minx1 << endl;
  //cout << "minx2" << minx2 << endl;
  //cout << "miny1" << miny1 << endl;
  //cout << "miny2" << miny2 << endl;
  //cout << "minz1" << minz1 << endl;
  //cout << "minz2" << minz2 << endl;
 }
 
  m += 1;

  //cout << "Px:" << L1.Px() << endl;
  //cout << "Py:" << L1.Py() << endl;
  //cout << "Pz:" << L1.Pz() << endl;
  //Output the data.
  //tau+
  //cout << "tao+" << endl;
  //cout << "px(real)(GeV):   " << p211x << endl;
  //cout << "px(reconstructed)(GeV):   " << p21x << endl;
  cout << (p21x - p211x) /p211x *100 << "\n";
  //cout << "py(real)(GeV):   " << p211y << endl;
  //cout << "py(reconstructed)(GeV):   " << p21y << endl;
  //cout << (p21y - p211y) / p211y *100 << "\n";
  //cout << "pz(real)(GeV):   " << p211z << endl;
  //cout << "pz(reconstructed)(GeV):   " << p21z << endl;
  //cout << (p21z - p211z) / p211z *100 << "\n";
  //cout << "E(real)(GeV):   " << E211 << endl;
  //cout << "E(reconstructed)(GeV):   " << E21 << endl;
  //cout << (E21 E211) / E211 *100 << "\n";
  //cout << (-pow(p211, 2) + pow(p21x, 2) + pow(p21y, 2) + pow(p21z, 2)) /pow(p21, 2) * 100 << endl;

  //tau-
  //cout << "tao-" << endl;
  //cout << "px(real)(GeV):   " << p221x << endl;
  //cout << "px(reconstructed)(GeV):   " << -p21x << endl;
  //std::cout << (-p21x -p221x) / p221x *100 << "\n";
  //cout << "py(real)(GeV):   " << p221y << endl;
  //cout << "py(reconstructed)(GeV):   " << -p21y << endl;
  //cout << (-p21y -p221y()) / p221y *100 << "\n";
  //cout << "pz(real)(GeV):   " << p221z << endl;
  //cout << "pz(reconstructed)(GeV):   " << -p21z << endl;
  //cout << (-p21z -p221z) / p221z *100 << "\n";
  //cout << "E(real)(GeV):   " << E221 << endl;
  //cout << "E(reconstructed)(GeV):   " << E22 << endl;
  //cout << (E22 - E221) / E221 *100 << "\n";
  //cout << (-pow(p221, 2) + pow(p22x, 2) + pow(p22y, 2) + pow(p22z, 2)) /pow(p22, 2) * 100 << endl;

  }}}

  // Total nummber of the events and the nummber of the events we could analysisby this means
  cout << m << endl;
  return 0;
}
