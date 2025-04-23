//##########################
//Topic:Select HNL Events
//Author: Yirui Liu
// Data: 2025.3.28
//##########################

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Math/Vector4D.h>
#include <algorithm>
#include <map>

using namespace std;

int select(){

//Define variables.
  int N=0;  
  vector<double> *p_e=0;
  vector<double> *p_x=0;
  vector<double> *p_y=0;
  vector<double> *p_z=0;
  vector<int> *p_id=0;
  vector<int> *mother1=0;
  vector<int> *mother2=0;
  vector<double> *mass=0;
  int npar=0;
  int n1=0;
  int n2=0;
  int n3=0;
  int n4=0;
  int n31=0;
  int n32=0;
  int n33=0;
  int n34=0;
  int n41=0;
  int n42=0;
  int n43=0;
  int n44=0;
  int q0e=0;
  int q0m=0;
  int q1e=0;
  int q1m=0;
  int q0ef=0;
  int q0mf=0;
  int q1ef=0;
  int q1mf=0;
  int NUM[10][10]={0};

//Get data form ".root" file.
  string dir = "/home/yrliu/Desktop/";
  string name = "sm.root";
  string filename = dir + name;

  TBranch *b_p_e;
  TBranch *b_p_id;
  TBranch *b_p_x;
  TBranch *b_p_y;
  TBranch *b_p_z;
  TBranch *b_npar;
  TBranch *b_mother1;
  TBranch *b_mother2;
  TBranch *b_mass;

  TFile *f = new TFile( filename.c_str() );
  TTree *my_tree = (TTree*)f->Get("lhedata");
 
  my_tree->SetBranchAddress("energy", &p_e, &b_p_e );
  my_tree->SetBranchAddress("pid", &p_id, &b_p_id );
  my_tree->SetBranchAddress("px", &p_x, &b_p_x );
  my_tree->SetBranchAddress("py", &p_y, &b_p_y );
  my_tree->SetBranchAddress("pz", &p_z, &b_p_z );
  my_tree->SetBranchAddress("numParticles", &npar, &b_npar );
  my_tree->SetBranchAddress("mother1", &mother1, &b_mother1 );
  my_tree->SetBranchAddress("mother2", &mother2, &b_mother2 );
  my_tree->SetBranchAddress("mass", &mass, &b_mass);

//Get one event and analysis it.
  N = (int) my_tree->GetEntries();
  for(int i=0; i<N; i++){
        my_tree->GetEntry( i );
        
        //Determine whether or not the event is what we need.
	//3p---tau+---e tag
	if(npar == 11){
  	        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10)};
                unordered_set<int> set2 = {11, -11, -15, 15, 211, 211, -211, -16, 16, -12, 11};
                if(set1 == set2){
	
	                //Match data and particles one by one.
                        n42 = 0;
                        for(int i=0; i < npar; i+=1){
                                if(p_id->at(i)==-11 && mother1->at(i)==0){
                                        n1 = i;}//e+
                                if(p_id->at(i)==11 && mother1->at(i)==0){
                                        n2 = i;}//e-
                                if(p_id->at(i)==-15 && mother1->at(i)==1){
                                        n3 = i;}//tau+
                                if(p_id->at(i)==15 && mother1->at(i)==1){
                                        n4 = i;}//tau-
                                if(p_id->at(i)==-211 && mother1->at(i)==4){
                                        n41 = i;}//pi-
                                if(p_id->at(i)==211 && mother1->at(i)==4){
                                        if(n42==0){
                                                n42 = i;//pi+
                                        } else {
                                                n43 = i;}}//pi+
                                if(p_id->at(i)==-16 && mother1->at(i)==4){
                                        n44 = i;}//vt~
                                if(p_id->at(i)==16 && mother1->at(i)==3){
                                        n31 = i;}//vt
                                if(p_id->at(i)==-12 && mother1->at(i)==3){
                                        n32 = i;}//ve~
                                if(p_id->at(i)==11 && mother1->at(i)==3){
                                        n33 = i;}}//e-

			//Lorentz transformation
                        TLorentzVector L01;
                        L01.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        TLorentzVector L02;
                        L02.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));
                        
			TLorentzVector total = L01 +L02;
			TVector3 beta_cm = total.BoostVector();

                        TLorentzVector L1;
                        L1.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        L1.Boost(-beta_cm);

                        TLorentzVector L2;
                        L2.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));
                        L2.Boost(-beta_cm);

			TLorentzVector L3;
                        L3.SetPxPyPzE(p_x->at(n3), p_y->at(n3), p_z->at(n3), p_e->at(n3));
                        L3.Boost(-beta_cm);
                        
			TLorentzVector L4;
                        L4.SetPxPyPzE(p_x->at(n4), p_y->at(n4), p_z->at(n4), p_e->at(n4));
                        L4.Boost(-beta_cm);

			TLorentzVector L31;
                        L31.SetPxPyPzE(p_x->at(n31), p_y->at(n31), p_z->at(n31), p_e->at(n31));
                        L31.Boost(-beta_cm);

                        TLorentzVector L32;
                        L32.SetPxPyPzE(p_x->at(n32), p_y->at(n32), p_z->at(n32), p_e->at(n32));
                        L32.Boost(-beta_cm);

			TLorentzVector L33;
                        L33.SetPxPyPzE(p_x->at(n33), p_y->at(n33), p_z->at(n33), p_e->at(n33));
                        L33.Boost(-beta_cm);

                        TLorentzVector L34;
                        L34.SetPxPyPzE(p_x->at(n34), p_y->at(n34), p_z->at(n34), p_e->at(n34));
                        L34.Boost(-beta_cm);

                        TLorentzVector L41;
                        L41.SetPxPyPzE(p_x->at(n41), p_y->at(n41), p_z->at(n41), p_e->at(n41));
                        L41.Boost(-beta_cm);

                        TLorentzVector L42;
                        L42.SetPxPyPzE(p_x->at(n42), p_y->at(n42), p_z->at(n42), p_e->at(n42));
                        L42.Boost(-beta_cm);

                        TLorentzVector L43;
                        L43.SetPxPyPzE(p_x->at(n43), p_y->at(n43), p_z->at(n43), p_e->at(n43));
                        L43.Boost(-beta_cm);

			//selection
		//	if( sqrt(pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2))/L4.P() < 0.5 ){
                  //      q0e+=1;}

		//	double mh2 = pow(L41.E()+L42.E()+L43.E(),2)-pow(L41.Px()+L42.Px()+L43.Px(),2)-pow(L41.Py()+L42.Py()+L43.Py(),2)-pow(L41.Pz()+L42.Pz()+L43.Pz(),2);

                  //      if( sqrt(mh2)/mass->at(n3) < 0.52 || sqrt(mh2)/mass->at(n3) > 0.8){
                    //            if((L41.E()+L42.E()+L43.E())/L3.E() < 0.96){
                                //cout << (L41.E()+L42.E()+L43.E())/L4.E() << endl;
                      //          q0e += 1;}}

                        //selection--thrust
/*			double max = 0;
		        if ( pow(L31.Px() + L32.Px() + L33.Px(), 2) + pow(L31.Py() + L32.Py() + L33.Py(), 2) + pow(L31.Pz() + L32.Pz() + L33.Pz(), 2) > pow(-L31.Px() + L32.Px() + L33.Px(), 2) + pow(-L31.Py() + L32.Py() + L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() + L33.Pz(), 2) ){
				max = pow(L31.Px() + L32.Px() + L33.Px(), 2) + pow(L31.Py() + L32.Py() + L33.Py(), 2) + pow(L31.Pz() + L32.Pz() + L33.Pz(), 2);}
			else {
				max = pow(-L31.Px() + L32.Px() + L33.Px(), 2) + pow(-L31.Py() + L32.Py() + L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() + L33.Pz(), 2);}
			if ( pow(L31.Px() - L32.Px() + L33.Px(), 2) + pow(L31.Py() - L32.Py() + L33.Py(), 2) + pow(L31.Pz() - L32.Pz() + L33.Pz(), 2) > max ){
				max = pow(L31.Px() - L32.Px() + L33.Px(), 2) + pow(L31.Py() - L32.Py() + L33.Py(), 2) + pow(L31.Pz() - L32.Pz() + L33.Pz(), 2);}
			if ( pow(- L31.Px() - L32.Px() + L33.Px(), 2) + pow(-L31.Py() - L32.Py() + L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() + L33.Pz(), 2) > max ){
                                max = pow(- L31.Px() - L32.Px() + L33.Px(), 2) + pow(- L31.Py() - L32.Py() + L33.Py(), 2) + pow(- L31.Pz() - L32.Pz() + L33.Pz(), 2);}
			if ( pow(L31.Px() + L32.Px() - L33.Px(), 2) + pow(L31.Py() + L32.Py() - L33.Py(), 2) + pow(L31.Pz() + L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(L31.Px() + L32.Px() - L33.Px(), 2) + pow(L31.Py() + L32.Py() - L33.Py(), 2) + pow(L31.Pz() + L32.Pz() - L33.Pz(), 2);}
			if ( pow(-L31.Px() + L32.Px() - L33.Px(), 2) + pow(-L31.Py() + L32.Py() - L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(-L31.Px() + L32.Px() - L33.Px(), 2) + pow(-L31.Py() + L32.Py() - L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() - L33.Pz(), 2);}
			if ( pow(L31.Px() - L32.Px() - L33.Px(), 2) + pow(L31.Py() - L32.Py() - L33.Py(), 2) + pow(L31.Pz() - L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(L31.Px() - L32.Px() - L33.Px(), 2) + pow(L31.Py() - L32.Py() - L33.Py(), 2) + pow(L31.Pz() - L32.Pz() - L33.Pz(), 2);}
			if ( pow(-L31.Px() - L32.Px() - L33.Px(), 2) + pow(-L31.Py() - L32.Py() - L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(-L31.Px() - L32.Px() - L33.Px(), 2) + pow(-L31.Py() - L32.Py() - L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() - L33.Pz(), 2);}
                       
			if( sqrt(max)/(sqrt(pow(L31.Px(), 2) + pow(L31.Py(), 2) + pow(L31.Pz(), 2)) + sqrt(pow(L32.Px(), 2) + pow(L32.Py(), 2) + pow(L32.Pz(), 2)) + sqrt(pow(L33.Px(), 2) + pow(L33.Py(), 2) + pow(L33.Pz(), 2))) > 0.85){							
				q0e += 1;}
			
			//selection---pt loss
			//if( fabs(L31.Pt() + L32.Pt() + L33.Pt() - L3.Pt()) > 0.009 * 10.564){
                          //      q0e += 1;}

			//selection---3p mass
                       // if(pow(p_e->at(n31)+p_e->at(n32)+p_e->at(n33),2)-pow(p_x->at(n31)+p_x->at(n32)+p_x->at(n33),2)-pow(p_y->at(n31)+p_y->at(n32)+p_y->at(n33),2)-pow(p_z->at(n31)+p_z->at(n32)+p_z->at(n33),2)<=1.776*1.776){
                       // q0e += 1;}
	
//		double mh2 = pow(p_e->at(n31)+p_e->at(n32)+p_e->at(n33),2)-pow(p_x->at(n31)+p_x->at(n32)+p_x->at(n33),2)-pow(p_y->at(n31)+p_y->at(n32)+p_y->at(n33),2)-pow(p_z->at(n31)+p_z->at(n32)+p_z->at(n33),2);

//                        if(sqrt(mh2)/mass->at(n3) > 0.5614 + 0.1293 || sqrt(mh2)/mass->at(n3) < 0.5614 - 0.1293 ){
 //                       q0e += 1;}

//			cout << sqrt(mh2)/mass->at(n3) << endl;
*/		}}
	
	//3p---tau+---mu tag
	if(npar == 11){
		unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10)};
                unordered_set<int> set2 = {11, -11, -15, 15, 211, 211, -211, -16, 16, -14, 13};
                if(set1 == set2){
  
		        //Match data and particles one by one.
                        n42 = 0;
                        for(int i=0; i < npar; i+=1){
                                if(p_id->at(i)==-11 && mother1->at(i)==0){
                                        n1 = i;}//e+
                                if(p_id->at(i)==11 && mother1->at(i)==0){
                                        n2 = i;}//e-
                                if(p_id->at(i)==-15 && mother1->at(i)==1){
                                        n3 = i;}//tau+
                                if(p_id->at(i)==15 && mother1->at(i)==1){
                                        n4 = i;}//tau-
                                if(p_id->at(i)==-211 && mother1->at(i)==4){
                                        n41 = i;}//pi-
                                if(p_id->at(i)==211 && mother1->at(i)==4){
                                        if(n42==0){
                                                n42 = i;//pi+
                                        } else {
                                                n43 = i;}}//pi+
                                if(p_id->at(i)==-16 && mother1->at(i)==4){
                                        n44 = i;}//vt~
                                if(p_id->at(i)==16 && mother1->at(i)==3){
                                        n31 = i;}//vt
                                if(p_id->at(i)==-14 && mother1->at(i)==3){
                                        n32 = i;}//vm~
                                if(p_id->at(i)==13 && mother1->at(i)==3){
                                        n33 = i;}}//m-

			//Lorentz transformation
                        TLorentzVector L01;
                        L01.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        TLorentzVector L02;
                        L02.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));

                        TLorentzVector total = L01 +L02;
                        TVector3 beta_cm = total.BoostVector();

                        TLorentzVector L1;
                        L1.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        L1.Boost(-beta_cm);

                        TLorentzVector L2;
                        L2.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));
                        L2.Boost(-beta_cm);

                        TLorentzVector L3;
                        L3.SetPxPyPzE(p_x->at(n3), p_y->at(n3), p_z->at(n3), p_e->at(n3));
                        L3.Boost(-beta_cm);

                        TLorentzVector L4;
                        L4.SetPxPyPzE(p_x->at(n4), p_y->at(n4), p_z->at(n4), p_e->at(n4));
                        L4.Boost(-beta_cm);

                        TLorentzVector L31;
                        L31.SetPxPyPzE(p_x->at(n31), p_y->at(n31), p_z->at(n31), p_e->at(n31));
                        L31.Boost(-beta_cm);

			TLorentzVector L32;
                        L32.SetPxPyPzE(p_x->at(n32), p_y->at(n32), p_z->at(n32), p_e->at(n32));
                        L32.Boost(-beta_cm);

                        TLorentzVector L33;
                        L33.SetPxPyPzE(p_x->at(n33), p_y->at(n33), p_z->at(n33), p_e->at(n33));
                        L33.Boost(-beta_cm);

                        TLorentzVector L34;
                        L34.SetPxPyPzE(p_x->at(n34), p_y->at(n34), p_z->at(n34), p_e->at(n34));
                        L34.Boost(-beta_cm);

                        TLorentzVector L41;
                        L41.SetPxPyPzE(p_x->at(n41), p_y->at(n41), p_z->at(n41), p_e->at(n41));
                        L41.Boost(-beta_cm);

                        TLorentzVector L42;
                        L42.SetPxPyPzE(p_x->at(n42), p_y->at(n42), p_z->at(n42), p_e->at(n42));
                        L42.Boost(-beta_cm);

                        TLorentzVector L43;
                        L43.SetPxPyPzE(p_x->at(n43), p_y->at(n43), p_z->at(n43), p_e->at(n43));
                        L43.Boost(-beta_cm);  	

			//selection
		//	if( sqrt(pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2))/L4.P() < 0.5 ){
                  //      q0m+=1;}

	//		double mh2 = pow(L41.E()+L42.E()+L43.E(),2)-pow(L41.Px()+L42.Px()+L43.Px(),2)-pow(L41.Py()+L42.Py()+L43.Py(),2)-pow(L41.Pz()+L42.Pz()+L43.Pz(),2);

          //              if( sqrt(mh2)/mass->at(n3) < 0.52 || sqrt(mh2)/mass->at(n3) >0.8){                                if((L41.E()+L42.E()+L43.E())/L3.E() < 0.96){
                                //cout << (L41.E()+L42.E()+L43.E())/L4.E() << endl;
            //                    q0m += 1;}}

			//selection--thrust
                       /* double max = 0;
                        if ( pow(L31.Px() + L32.Px() + L33.Px(), 2) + pow(L31.Py() + L32.Py() + L33.Py(), 2) + pow(L31.Pz() + L32.Pz() + L33.Pz(), 2) > pow(-L31.Px() + L32.Px() + L33.Px(), 2) + pow(-L31.Py() + L32.Py() + L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() + L33.Pz(), 2) ){
                                max = pow(L31.Px() + L32.Px() + L33.Px(), 2) + pow(L31.Py() + L32.Py() + L33.Py(), 2) + pow(L31.Pz() + L32.Pz() + L33.Pz(), 2);}
                        else {
                                max = pow(-L31.Px() + L32.Px() + L33.Px(), 2) + pow(-L31.Py() + L32.Py() + L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() + L33.Pz(), 2);}
                        if ( pow(L31.Px() - L32.Px() + L33.Px(), 2) + pow(L31.Py() - L32.Py() + L33.Py(), 2) + pow(L31.Pz() - L32.Pz() + L33.Pz(), 2) > max ){
                                max = pow(L31.Px() - L32.Px() + L33.Px(), 2) + pow(L31.Py() - L32.Py() + L33.Py(), 2) + pow(L31.Pz() - L32.Pz() + L33.Pz(), 2);}
                        if ( pow(- L31.Px() - L32.Px() + L33.Px(), 2) + pow(-L31.Py() - L32.Py() + L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() + L33.Pz(), 2) > max ){
                                max = pow(- L31.Px() - L32.Px() + L33.Px(), 2) + pow(- L31.Py() - L32.Py() + L33.Py(), 2) + pow(- L31.Pz() - L32.Pz() + L33.Pz(), 2);}
                        if ( pow(L31.Px() + L32.Px() - L33.Px(), 2) + pow(L31.Py() + L32.Py() - L33.Py(), 2) + pow(L31.Pz() + L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(L31.Px() + L32.Px() - L33.Px(), 2) + pow(L31.Py() + L32.Py() - L33.Py(), 2) + pow(L31.Pz() + L32.Pz() - L33.Pz(), 2);}
                        if ( pow(-L31.Px() + L32.Px() - L33.Px(), 2) + pow(-L31.Py() + L32.Py() - L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(-L31.Px() + L32.Px() - L33.Px(), 2) + pow(-L31.Py() + L32.Py() - L33.Py(), 2) + pow(-L31.Pz() + L32.Pz() - L33.Pz(), 2);}
                        if ( pow(L31.Px() - L32.Px() - L33.Px(), 2) + pow(L31.Py() - L32.Py() - L33.Py(), 2) + pow(L31.Pz() - L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(L31.Px() - L32.Px() - L33.Px(), 2) + pow(L31.Py() - L32.Py() - L33.Py(), 2) + pow(L31.Pz() - L32.Pz() - L33.Pz(), 2);}
			if ( pow(-L31.Px() - L32.Px() - L33.Px(), 2) + pow(-L31.Py() - L32.Py() - L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() - L33.Pz(), 2) > max ){
                                max = pow(-L31.Px() - L32.Px() - L33.Px(), 2) + pow(-L31.Py() - L32.Py() - L33.Py(), 2) + pow(-L31.Pz() - L32.Pz() - L33.Pz(), 2);}

                        if( sqrt(max)/(sqrt(pow(L31.Px(), 2) + pow(L31.Py(), 2) + pow(L31.Pz(), 2)) + sqrt(pow(L32.Px(), 2) + pow(L32.Py(), 2) + pow(L32.Pz(), 2)) + sqrt(pow(L33.Px(), 2) + pow(L33.Py(), 2) + pow(L33.Pz(), 2))) > 0.85){
                                q0m += 1;}*/

			//selection---pt loss
    		      //  if( fabs(L31.Pt() + L32.Pt() + L33.Pt() + L43.Pt() - L1.Pt() - L2.Pt() - L3.Pt()) > 0.009 * 10.564){
                        //        q0m += 1;}

			//selection---3p mass
                       // if(pow(p_e->at(n31)+p_e->at(n32)+p_e->at(n33),2)-pow(p_x->at(n31)+p_x->at(n32)+p_x->at(n33),2)-pow(p_y->at(n31)+p_y->at(n32)+p_y->at(n33),2)-pow(p_z->at(n31)+p_z->at(n32)+p_z->at(n33),2)<=1.776*1.776){
                       // q0m += 1;}

		//	double mh2 = pow(p_e->at(n31)+p_e->at(n32)+p_e->at(n33),2)-pow(p_x->at(n31)+p_x->at(n32)+p_x->at(n33),2)-pow(p_y->at(n31)+p_y->at(n32)+p_y->at(n33),2)-pow(p_z->at(n31)+p_z->at(n32)+p_z->at(n33),2);

  //                      if(sqrt(mh2)/mass->at(n3) > 0.5614 + 0.1293  || sqrt(mh2)/mass->at(n3) <0.5614 - 0.1293 ){
    //                    q0m += 1;}

		//	cout << sqrt(mh2)/mass->at(n3) << endl;
		}}
         
        //3p---tau- ---e tag
	if(npar == 11){
                unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10)};
                unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, -211, 16, -16, 12, -11};
                if(set1 == set2){
               
			//Match data and particles one by one.
                        n42 = 0;
			for(int i=0; i < npar; i+=1){
				if(p_id->at(i)==-11 && mother1->at(i)==0){
					n1 = i;}//e+
				if(p_id->at(i)==11 && mother1->at(i)==0){
                                        n2 = i;}//e-
                                if(p_id->at(i)==-15 && mother1->at(i)==1){
                                        n3 = i;}//tau+
				if(p_id->at(i)==15 && mother1->at(i)==1){
                                        n4 = i;}//tau-
                                if(p_id->at(i)==211 && mother1->at(i)==4){
                                        n41 = i;}//pi+
				if(p_id->at(i)==-211 && mother1->at(i)==4){
					if(n42==0){
						n42 = i;//pi-
					} else {
						n43 = i;}}//pi-
                                if(p_id->at(i)==16 && mother1->at(i)==4){
                                        n44 = i;}//vt
                                if(p_id->at(i)==-16 && mother1->at(i)==3){
                                        n31 = i;}//vt~
                                if(p_id->at(i)==12 && mother1->at(i)==3){
                                        n32 = i;}//ve
				if(p_id->at(i)==-11 && mother1->at(i)==3){
                                        n33 = i;}}//e+

			//Lorentz transformation
                        TLorentzVector L01;
                        L01.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        TLorentzVector L02;
                        L02.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));

                        TLorentzVector total = L01 + L02;
			TVector3 beta_cm = total.BoostVector();

                        TLorentzVector L1;
                        L1.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        L1.Boost(-beta_cm);

                        TLorentzVector L2;
                        L2.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));
                        L2.Boost(-beta_cm);

                        TLorentzVector L3;
                        L3.SetPxPyPzE(p_x->at(n3), p_y->at(n3), p_z->at(n3), p_e->at(n3));
                        L3.Boost(-beta_cm);

                        TLorentzVector L4;
                        L4.SetPxPyPzE(p_x->at(n4), p_y->at(n4), p_z->at(n4), p_e->at(n4));
                        L4.Boost(-beta_cm);

                        TLorentzVector L31;
                        L31.SetPxPyPzE(p_x->at(n31), p_y->at(n31), p_z->at(n31), p_e->at(n31));
			L31.Boost(-beta_cm);

			TLorentzVector L32;
                        L32.SetPxPyPzE(p_x->at(n32), p_y->at(n32), p_z->at(n32), p_e->at(n32));
                        L32.Boost(-beta_cm);

                        TLorentzVector L33;
                        L33.SetPxPyPzE(p_x->at(n33), p_y->at(n33), p_z->at(n33), p_e->at(n33));
			L33.Boost(-beta_cm);

                        TLorentzVector L41;
                        L41.SetPxPyPzE(p_x->at(n41), p_y->at(n41), p_z->at(n41), p_e->at(n41));
			L41.Boost(-beta_cm);

                        TLorentzVector L42;
                        L42.SetPxPyPzE(p_x->at(n42), p_y->at(n42), p_z->at(n42), p_e->at(n42));
			L42.Boost(-beta_cm);

                        TLorentzVector L43;
                        L43.SetPxPyPzE(p_x->at(n43), p_y->at(n43), p_z->at(n43), p_e->at(n43));
			L43.Boost(-beta_cm);

			TLorentzVector L44;
                        L44.SetPxPyPzE(p_x->at(n44), p_y->at(n44), p_z->at(n44), p_e->at(n44));
                        L44.Boost(-beta_cm);
                        
			//selection
			double mh2 = pow(L41.E()+L42.E()+L43.E(),2)-pow(L41.Px()+L42.Px()+L43.Px(),2)-pow(L41.Py()+L42.Py()+L43.Py(),2)-pow(L41.Pz()+L42.Pz()+L43.Pz(),2);

		//cout << L41.Pz() + L42.Pz() + L43.Pz() + L33.Pz() << endl;}
		
		//		double mh3 = (4 * p_e->at(n4) * (mh2 + pow(mass->at(n4),2) - pow(mass->at(n44),2))-sqrt(16*(pow(p_e->at(n4),2)-pow(mass->at(n4),2))*(pow(-mh2 + pow(mass->at(n4),2) - pow(mass->at(n44),2),2) - 4*pow(mass->at(n44),2)*mh2)))/8/pow(mass->at(n4),2) - p_e->at(n4) + sqrt(pow(mass->at(n44),2)+pow(mass->at(n4)/2*(mh2 - pow(mass->at(n4),2) - pow(mass->at(n44),2))/pow(mass->at(n4),2)*sqrt(pow(p_e->at(n4)/mass->at(n4),2)-1)+p_e->at(n4)/2*sqrt((1-pow((sqrt(mh2)+mass->at(n44))/mass->at(n4),2))*(1-pow((sqrt(mh2)-mass->at(n44))/mass->at(n4),2))),2));
//			cout << mh3 << endl;
                        //selection--thrust
                      /*  double max = 0;
                        if ( pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2) > pow(-L41.Px() + L42.Px() + L43.Px(), 2) + pow(-L41.Py() + L42.Py() + L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() + L43.Pz(), 2) ){
				max = pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2);}
                        else {
                                max = pow(-L41.Px() + L42.Px() + L43.Px(), 2) + pow(-L41.Py() + L42.Py() + L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() + L43.Pz(), 2);}
                        if ( pow(L41.Px() - L42.Px() + L43.Px(), 2) + pow(L41.Py() - L42.Py() + L43.Py(), 2) + pow(L41.Pz() - L42.Pz() + L43.Pz(), 2) > max ){
                                max = pow(L41.Px() - L42.Px() + L43.Px(), 2) + pow(L41.Py() - L42.Py() + L43.Py(), 2) + pow(L41.Pz() - L42.Pz() + L43.Pz(), 2);}
                        if ( pow(- L41.Px() - L42.Px() + L43.Px(), 2) + pow(-L41.Py() - L42.Py() + L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() + L43.Pz(), 2) > max ){
                                max = pow(- L41.Px() - L42.Px() + L43.Px(), 2) + pow(- L41.Py() - L42.Py() + L43.Py(), 2) + pow(- L41.Pz() - L42.Pz() + L43.Pz(), 2);}
                        if ( pow(L41.Px() + L42.Px() - L43.Px(), 2) + pow(L41.Py() + L42.Py() - L43.Py(), 2) + pow(L41.Pz() + L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(L41.Px() + L42.Px() - L43.Px(), 2) + pow(L41.Py() + L42.Py() - L43.Py(), 2) + pow(L41.Pz() + L42.Pz() - L43.Pz(), 2);}
                        if ( pow(-L41.Px() + L42.Px() - L43.Px(), 2) + pow(-L41.Py() + L42.Py() - L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(-L41.Px() + L42.Px() - L43.Px(), 2) + pow(-L41.Py() + L42.Py() - L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() - L43.Pz(), 2);}
                        if ( pow(L41.Px() - L42.Px() - L43.Px(), 2) + pow(L41.Py() - L42.Py() - L43.Py(), 2) + pow(L41.Pz() - L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(L41.Px() - L42.Px() - L43.Px(), 2) + pow(L41.Py() - L42.Py() - L43.Py(), 2) + pow(L41.Pz() - L42.Pz() - L43.Pz(), 2);}
                        if ( pow(-L41.Px() - L42.Px() - L43.Px(), 2) + pow(-L41.Py() - L42.Py() - L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() - L43.Pz(), 2) > max ){
				max = pow(-L41.Px() - L42.Px() - L43.Px(), 2) + pow(-L41.Py() - L42.Py() - L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() - L43.Pz(), 2);}

 			cout << sqrt(max)/(sqrt(pow(L41.Px(), 2) + pow(L41.Py(), 2) + pow(L41.Pz(), 2)) + sqrt(pow(L42.Px(), 2) + pow(L42.Py(), 2) + pow(L42.Pz(), 2)) + sqrt(pow(L43.Px(), 2) + pow(L43.Py(), 2) + pow(L43.Pz(), 2))) << endl;
	*/
			//selection---pt lose
			//cout <<  fabs(L41.Pt() + L42.Pt() + L43.Pt() + L33.Pt()) << endl;
		        //cout << (L41.E()+L42.E()+L43.E())/L4.E() << endl;
			//q0ef+=1;
			//cout << sqrt(mh2)/mass->at(n4) << endl;
			double mh24 = sqrt(mh2)/mass->at(n4);
			double eh24 = (L41.E()+L42.E()+L43.E())/L4.E();
			for (float i = 0.1; i < 1.1; i+=0.1) {
				if(mh24 <= i && mh24 >i-0.1 ){
			                for (float j = 0.1; j < 1.1; j+=0.1) {
                                                if(eh24 <= j && eh24 > (j-0.1) ) {
						        int j1 = j/0.1-1; 
							int i1 = i/0.1-1;
							NUM[j1][i1]+=1;}}}}


        }}


        //3p---tau- ---m tag
	if(npar == 11){
        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10)};
        unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, -211, 16, -16, 14, -13};
                if(set1 == set2){

			//Match data and particles one by one.
                        n42 = 0;
                        for(int i=0; i < npar; i+=1){
                                if(p_id->at(i)==-11 && mother1->at(i)==0){
                                        n1 = i;}//e+
                                if(p_id->at(i)==11 && mother1->at(i)==0){
                                        n2 = i;}//e-
                                if(p_id->at(i)==-15 && mother1->at(i)==1){
                                        n3 = i;}//tau+
                                if(p_id->at(i)==15 && mother1->at(i)==1){
                                        n4 = i;}//tau-
                                if(p_id->at(i)==211 && mother1->at(i)==4){
                                        n41 = i;}//pi+
                                if(p_id->at(i)==-211 && mother1->at(i)==4){
                                        if(n42==0){
                                                n42 = i;//pi-
                                        } else {
                                                n43 = i;}}//pi-
                                if(p_id->at(i)==16 && mother1->at(i)==4){
                                        n44 = i;}//vt
                                if(p_id->at(i)==-16 && mother1->at(i)==3){
                                        n31 = i;}//vt~
                                if(p_id->at(i)==14 && mother1->at(i)==3){
                                        n32 = i;}//vm
                                if(p_id->at(i)==-13 && mother1->at(i)==3){
                                        n33 = i;}}//mu+

			//Lorentz transformation
                        TLorentzVector L01;
                        L01.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        TLorentzVector L02;
                        L02.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));

                        TLorentzVector total = L01 + L02;
                        TVector3 beta_cm = total.BoostVector();

                        TLorentzVector L1;
                        L1.SetPxPyPzE(p_x->at(n1), p_y->at(n1), p_z->at(n1), p_e->at(n1));
                        L1.Boost(-beta_cm);

                        TLorentzVector L2;
                        L2.SetPxPyPzE(p_x->at(n2), p_y->at(n2), p_z->at(n2), p_e->at(n2));
                        L2.Boost(-beta_cm);

                        TLorentzVector L3;
                        L3.SetPxPyPzE(p_x->at(n3), p_y->at(n3), p_z->at(n3), p_e->at(n3));
                        L3.Boost(-beta_cm);

                        TLorentzVector L4;
                        L4.SetPxPyPzE(p_x->at(n4), p_y->at(n4), p_z->at(n4), p_e->at(n4));
                        L4.Boost(-beta_cm);

                        TLorentzVector L31;
                        L31.SetPxPyPzE(p_x->at(n31), p_y->at(n31), p_z->at(n31), p_e->at(n31));
                        L31.Boost(-beta_cm);

			TLorentzVector L32;
                        L32.SetPxPyPzE(p_x->at(n32), p_y->at(n32), p_z->at(n32), p_e->at(n32));
                        L32.Boost(-beta_cm);

                        TLorentzVector L33;
                        L33.SetPxPyPzE(p_x->at(n33), p_y->at(n33), p_z->at(n33), p_e->at(n33));
                        L33.Boost(-beta_cm);

                        TLorentzVector L41;
                        L41.SetPxPyPzE(p_x->at(n41), p_y->at(n41), p_z->at(n41), p_e->at(n41));
                        L41.Boost(-beta_cm);

                        TLorentzVector L42;
                        L42.SetPxPyPzE(p_x->at(n42), p_y->at(n42), p_z->at(n42), p_e->at(n42));
                        L42.Boost(-beta_cm);

                        TLorentzVector L43;
                        L43.SetPxPyPzE(p_x->at(n43), p_y->at(n43), p_z->at(n43), p_e->at(n43));
                        L43.Boost(-beta_cm);

                        TLorentzVector L44;
                        L44.SetPxPyPzE(p_x->at(n44), p_y->at(n44), p_z->at(n44), p_e->at(n44));
                        L44.Boost(-beta_cm);

			//selection
		//	if( sqrt(pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2))/L4.P() < 0.5 ){
                  //      q0mf+=1;}

		//	double mh2 = pow(L41.E()+L42.E()+L43.E(),2)-pow(L41.Px()+L42.Px()+L43.Px(),2)-pow(L41.Py()+L42.Py()+L43.Py(),2)-pow(L41.Pz()+L42.Pz()+L43.Pz(),2);

                  //      if( sqrt(mh2)/mass->at(n4) < 0.52 || sqrt(mh2)/mass->at(n4) > 0.8){
                    //            if((L41.E()+L42.E()+L43.E())/L4.E() < 0.96){
                                //cout << (L41.E()+L42.E()+L43.E())/L4.E() << endl;
                      //          q0mf += 1;}}

			//selection--thrust
               /*         double max = 0;
                        if ( pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2) > pow(-L41.Px() + L42.Px() + L43.Px(), 2) + pow(-L41.Py() + L42.Py() + L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() + L43.Pz(), 2) ){
                                max = pow(L41.Px() + L42.Px() + L43.Px(), 2) + pow(L41.Py() + L42.Py() + L43.Py(), 2) + pow(L41.Pz() + L42.Pz() + L43.Pz(), 2);}
                        else {
                                max = pow(-L41.Px() + L42.Px() + L43.Px(), 2) + pow(-L41.Py() + L42.Py() + L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() + L43.Pz(), 2);}
                        if ( pow(L41.Px() - L42.Px() + L43.Px(), 2) + pow(L41.Py() - L42.Py() + L43.Py(), 2) + pow(L41.Pz() - L42.Pz() + L43.Pz(), 2) > max ){
                                max = pow(L41.Px() - L42.Px() + L43.Px(), 2) + pow(L41.Py() - L42.Py() + L43.Py(), 2) + pow(L41.Pz() - L42.Pz() + L43.Pz(), 2);}
                        if ( pow(- L41.Px() - L42.Px() + L43.Px(), 2) + pow(-L41.Py() - L42.Py() + L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() + L43.Pz(), 2) > max ){
                                max = pow(- L41.Px() - L42.Px() + L43.Px(), 2) + pow(- L41.Py() - L42.Py() + L43.Py(), 2) + pow(- L41.Pz() - L42.Pz() + L43.Pz(), 2);}
                        if ( pow(L41.Px() + L42.Px() - L43.Px(), 2) + pow(L41.Py() + L42.Py() - L43.Py(), 2) + pow(L41.Pz() + L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(L41.Px() + L42.Px() - L43.Px(), 2) + pow(L41.Py() + L42.Py() - L43.Py(), 2) + pow(L41.Pz() + L42.Pz() - L43.Pz(), 2);}
                        if ( pow(-L41.Px() + L42.Px() - L43.Px(), 2) + pow(-L41.Py() + L42.Py() - L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(-L41.Px() + L42.Px() - L43.Px(), 2) + pow(-L41.Py() + L42.Py() - L43.Py(), 2) + pow(-L41.Pz() + L42.Pz() - L43.Pz(), 2);}
                        if ( pow(L41.Px() - L42.Px() - L43.Px(), 2) + pow(L41.Py() - L42.Py() - L43.Py(), 2) + pow(L41.Pz() - L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(L41.Px() - L42.Px() - L43.Px(), 2) + pow(L41.Py() - L42.Py() - L43.Py(), 2) + pow(L41.Pz() - L42.Pz() - L43.Pz(), 2);}
			if ( pow(-L41.Px() - L42.Px() - L43.Px(), 2) + pow(-L41.Py() - L42.Py() - L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() - L43.Pz(), 2) > max ){
                                max = pow(-L41.Px() - L42.Px() - L43.Px(), 2) + pow(-L41.Py() - L42.Py() - L43.Py(), 2) + pow(-L41.Pz() - L42.Pz() - L43.Pz(), 2);}

                        if(sqrt(max)/(sqrt(pow(L41.Px(), 2) + pow(L41.Py(), 2) + pow(L41.Pz(), 2)) + sqrt(pow(L42.Px(), 2) + pow(L42.Py(), 2) + pow(L42.Pz(), 2)) + sqrt(pow(L43.Px(), 2) + pow(L43.Py(), 2) + pow(L43.Pz(), 2))) >0.85){
                                q0mf += 1;}
*/
			//selection---pt loss
			//if( fabs(L41.Pt() + L42.Pt() + L43.Pt() - L4.Pt()) > 0.009 * 10.564){
                              //  q0mf += 1;}

			//selection---3p mass
//                        if(pow(p_e->at(n41)+p_e->at(n42)+p_e->at(n43),2)-pow(p_x->at(n41)+p_x->at(n42)+p_x->at(n43),2)-pow(p_y->at(n41)+p_y->at(n42)+p_y->at(n43),2)-pow(p_z->at(n41)+p_z->at(n42)+p_z->at(n43),2)<=1.776*1.776){
  //                      q0mf += 1;}

		//	double mh2 = pow(p_e->at(n41)+p_e->at(n42)+p_e->at(n43),2)-pow(p_x->at(n41)+p_x->at(n42)+p_x->at(n43),2)-pow(p_y->at(n41)+p_y->at(n42)+p_y->at(n43),2)-pow(p_z->at(n41)+p_z->at(n42)+p_z->at(n43),2);

		//	cout << sqrt(mh2)/mass->at(n3) << endl;
  //                      if(sqrt(mh2)/mass->at(n4) > 0.5614 + 0.1293 || sqrt(mh2)/mass->at(n4) < 0.5614 - 0.1293 ){
    //                    q0mf += 1;}

		}}

        //4p---tau+---e tag
/*	if(npar == 12){
        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10), p_id->at(11)};
        unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, 211, 111, -16, 16, -12, 11};
        if(set1 == set2){
                q1e +=1;}}

	//4p---tau+---mu tag
        if(npar == 12){
        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10), p_id->at(11)};
        unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, 211, 111, -16, 16, -14, 13};
        if(set1 == set2){
                q1m +=1;}}

	//4p---tau- ---e tag
        if(npar == 12){
        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10), p_id->at(11)};
        unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, -211, 111, 16, -16, 12, -11};
        if(set1 == set2){
                q1ef +=1;}}

	//4p---tau- ---mu tag
        if(npar == 12){
        unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10), p_id->at(11)};
        unordered_set<int> set2 = {11, -11, -15, 15, 211, -211, -211, 111, -16, 16, 14, -13};
        if(set1 == set2){
                q1mf +=1;}}
*/


  }

//Output the nummber of events of each bin.
  double qn1 = 0;  
  for (int j = 0; j < 10; j+=1) {
          for (int i = 0; i < 10; i+=1) {
		  qn1 += NUM[j][i]/*253328.0*0.007682/15.108/1000000.0*/;
                  cout << NUM[j][i]/*253328.0*0.007682/15.108/1000000.0*/ << endl;}}

//  cout << q0e << endl;
//  cout << q0m << endl;
//  cout << q0ef << endl;
//  cout << q0mf << endl;

  cout << qn1 << endl;
  cout << endl;
  return 0;
}
