//##########################
// Author: Yirui Liu
// Data: 2025.1.18
//##########################

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Math/Vector4D.h>
#include <algorithm>
#include <map>

using namespace std;

int ana(){

//Define variables.
  int N;  
  vector<double> *p_E=0;
  vector<double> *p_Pt=0;
  vector<double> *p_Eta=0;
  vector<double> *p_Phi=0;
  vector<double> *p_status=0;
  vector<int> *p_id=0;
  vector<int> *p_VID=0;
  vector<int> *p_ChildVID=0;
  vector<double> *mass=0;
  int npar=0;
  float eventweight=0;
  int n0=0;
  int n01=0;
  int n1=0;
  int n11=0;
  int n2=0;
  int n21=0;
  int n31=0;
  int n32=0;

  int q=0;
  double q1=0;
  
//Get data form ".root" file.
  string dir = "/home/yrliu/Desktop/";
  string name = "lhe.root";
  string filename = dir + name;

  TBranch *b_p_E;
  TBranch *b_p_id;
  TBranch *b_p_Pt;
  TBranch *b_p_Eta;
  TBranch *b_p_Phi;
  TBranch *b_npar;
  TBranch *b_p_status;
  TBranch *b_p_VID;
  TBranch *b_p_ChildVID;
  TBranch *b_eventweight;
  TBranch *b_mass;

  TFile *f = new TFile( filename.c_str() );
  TTree *my_tree = (TTree*)f->Get("lhedata");
 
  my_tree->SetBranchAddress("energy", &p_E, &b_p_E );
  my_tree->SetBranchAddress("pid", &p_id, &b_p_id );
  my_tree->SetBranchAddress("px", &p_Pt, &b_p_Pt );
  my_tree->SetBranchAddress("py", &p_Eta, &b_p_Eta );
  my_tree->SetBranchAddress("pz", &p_Phi, &b_p_Phi );
  my_tree->SetBranchAddress("numParticles", &npar, &b_npar );
  my_tree->SetBranchAddress("status", &p_status, &b_p_status );
  my_tree->SetBranchAddress("mother1", &p_VID, &b_p_VID );
  my_tree->SetBranchAddress("mother2", &p_ChildVID, &b_p_ChildVID );
  my_tree->SetBranchAddress("eventweight", &eventweight, &b_eventweight);
  my_tree->SetBranchAddress("mass", &mass, &b_mass);

//Get one event and analysis it.
  N = (int) my_tree->GetEntries();
  for(int i=0; i<N; i++){
        my_tree->GetEntry( i );
        
	

        //Determine whether or not the event is what we need.
      //  if(npar == 8){       
  //	unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5), p_id->at(6), p_id->at(7), p_id->at(8), p_id->at(9), p_id->at(10), p_id->at(11)};
      //  unordered_set<int> set2 = {11, -11, -15, -16, 211, 15, 16, -211};
//	   if(set1 ==set2){

	      //Match data and particles one by one.
              map< int, int > Particles{{p_id->at(0), 0}, {p_id->at(1), 1}, {p_id->at(2), 2}, {p_id->at(3), 3}, {p_id->at(4), 4}, {p_id->at(5), 5}, {p_id->at(6), 6}, {p_id->at(7), 7}, {p_id->at(8), 8}, {p_id->at(9), 9}, {p_id->at(10), 10}, {p_id->at(11), 11}};
  //            n0 = Particles.find(-15)->second;
//	      n1 = Particles.find(-1)->second;
//	      n2 = Particles.find(2)->second;
	      n01 = Particles.find(9900016)->second;
//	      n11 = Particles.find(16)->second;
//	      n21 = Particles.find(-211)->second;
//              n31 = Particles.find(-11)->second;
//	      n32 = Particles.find(11)->second;
/*
		      //Further screening
		      if( p_status->at(n0) == 2 && p_status->at(n1) == 1 && p_status->at(n2) == 1 &&
		      p_status->at(n01) == 2 && p_status->at(n11) == 1 && p_status->at(n21) == 1 &&
		      p_ChildVID->at(n0) == p_VID->at(n1) && p_ChildVID->at(n0) == p_VID->at(n2) &&
		      p_ChildVID->at(n01) == p_VID->at(n11) && p_ChildVID->at(n01) == p_VID->at(n21)){           	
			    
			    //Lorentz transformation
			    TLorentzVector L0;
			    L0.SetPtEtaPhiE(p_Pt->at(n0), p_Eta->at(n0), p_Phi->at(n0), p_E->at(n0));
			    TLorentzVector L1;
			    L1.SetPtEtaPhiE(p_Pt->at(n1), p_Eta->at(n1), p_Phi->at(n1), p_E->at(n1));
			    L1.Boost(-L0.BoostVector());
			    TLorentzVector L2;
			    L2.SetPtEtaPhiE(p_Pt->at(n2), p_Eta->at(n2), p_Phi->at(n2), p_E->at(n2));
			    L2.Boost(-L0.BoostVector());
			    
			    TLorentzVector L01;
			    L01.SetPtEtaPhiE(p_Pt->at(n01), p_Eta->at(n01), p_Phi->at(n01), p_E->at(n01));
			    TLorentzVector L11;
			    L11.SetPtEtaPhiE(p_Pt->at(n11), p_Eta->at(n11), p_Phi->at(n11), p_E->at(n11));
			    L11.Boost(-L01.BoostVector());
			    TLorentzVector L21;
			    L21.SetPtEtaPhiE(p_Pt->at(n21), p_Eta->at(n21), p_Phi->at(n21), p_E->at(n21));
			    L21.Boost(-L01.BoostVector());
		*/	   
			    //Fill in data.
			    //cout << sqrt(pow(p_E->at(n1)+p_E->at(n2), 2) - pow(p_Eta->at(n1) + p_Eta->at(n2), 2) - pow(p_Pt->at(n1) + p_Pt->at(n2), 2) -pow(p_Phi->at(n1) + p_Phi->at(n2), 2))/mass->at(n0) << endl;
			    //cout << sqrt(pow(p_E->at(n01),2) / pow(mass->at(n01),2) -1) * 1e7 * 1.97327/1.473993 << endl;
                            q1 = q1 + sqrt(pow(p_E->at(n01),2) / pow(mass->at(n01),2) -1) * 1e0 * 1.97327/3.668499/10000;
			    q += 1;
		          
		} //}  }
  cout << q1 << endl;
  cout << q << endl;
  cout << endl;
  return 0;
}
