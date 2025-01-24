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
  int npar=0;
  int n0=0;
  int n01=0;
  int n1=0;
  int n11=0;
  int n2=0;
  int n21=0;
  
//Get data form ".root" file.
  string dir = "/home/yrliu/Desktop/";
  string name = "events.hepmc.root";
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

  TFile *f = new TFile( filename.c_str() );
  TTree *my_tree = (TTree*)f->Get("truth");
 
  my_tree->SetBranchAddress("p_E", &p_E, &b_p_E );
  my_tree->SetBranchAddress("p_id", &p_id, &b_p_id );
  my_tree->SetBranchAddress("p_Pt", &p_Pt, &b_p_Pt );
  my_tree->SetBranchAddress("p_Eta", &p_Eta, &b_p_Eta );
  my_tree->SetBranchAddress("p_Phi", &p_Phi, &b_p_Phi );
  my_tree->SetBranchAddress("npar", &npar, &b_npar );
  my_tree->SetBranchAddress("p_status", &p_status, &b_p_status );
  my_tree->SetBranchAddress("p_VID", &p_VID, &b_p_VID );
  my_tree->SetBranchAddress("p_ChildVID", &p_ChildVID, &b_p_ChildVID );

//Prepare for drawing.
 /* TCanvas *c11 = new TCanvas();
  TH1I *h11 = new TH1I( "h11", "The Px distribution of nu tau1", 40, -1, 1);
  TH1I *h12 = new TH1I( "h12", "The Px distribution of pi+", 40, -1, 1);
  TH1I *h21 = new TH1I( "h21", "The Py distribution of nu tau1", 40, -1, 1);
  TH1I *h22 = new TH1I( "h22", "The Py distribution of pi+", 40, -1, 1);
  TH1I *h31 = new TH1I( "h31", "The Pz distribution of nu tau1", 40, -1, 1);
  TH1I *h32 = new TH1I( "h32", "The Pz distribution of pi+", 40, -1, 1);
  TH1I *h41 = new TH1I( "h41", "The E distribution of nu tau1", 30, 0.86, 0.9);
  TH1I *h42 = new TH1I( "h42", "The E distribution of pi+", 30, 0.88, 0.92);
  TH1I *h111 = new TH1I( "h111", "The Px distribution of nu tau2", 40, -1, 1);
  TH1I *h121 = new TH1I( "h121", "The Px distribution of pi-", 40, -1, 1);
  TH1I *h211 = new TH1I( "h211", "The Py distribution of nu tau2", 40, -1, 1);
  TH1I *h221 = new TH1I( "h221", "The Py distribution of pi-", 40, -1, 1);
  TH1I *h311 = new TH1I( "h311", "The Pz distribution of nu tau2", 40, -1, 1);
  TH1I *h321 = new TH1I( "h321", "The Pz distribution of pi-", 40, -1, 1);
  TH1I *h411 = new TH1I( "h411", "The E distribution of nu tau2", 30, 0.86, 0.9);
  TH1I *h421 = new TH1I( "h421", "The E distribution of pi-", 30, 0.88, 0.92);
 */
  
//Get one event and analysis it.
  N = (int) my_tree->GetEntries();
  for(int i=0; i<N+1; i++){
        my_tree->GetEntry( i );

        //Determine whether or not the event is what we need.
	if(npar == 6){       
  	unordered_set<int> set1 = {p_id->at(0), p_id->at(1), p_id->at(2), p_id->at(3), p_id->at(4), p_id->at(5)};
        unordered_set<int> set2 = {-15, -16, 211, 15, 16, -211};
	   if(set1 ==set2){

	      //Match data and particles one by one.
              map< int, int > Particles{{p_id->at(0), 0}, {p_id->at(1), 1}, {p_id->at(2), 2}, {p_id->at(3), 3}, {p_id->at(4), 4}, {p_id->at(5), 5}};
              n0 = Particles.find(-15)->second;
	      n1 = Particles.find(-16)->second;
	      n2 = Particles.find(211)->second;
	      n01 = Particles.find(15)->second;
	      n11 = Particles.find(16)->second;
	      n21 = Particles.find(-211)->second;

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
		   
		    //Fill in data.   
                    cout << L1.Px() << "\n";
		    /*h11->Fill(L1.Px());
		    h12->Fill(L2.Px());
		    h21->Fill(L1.Py());
		    h22->Fill(L2.Py());
		    h31->Fill(L1.Pz());
		    h32->Fill(L2.Pz());
		    h41->Fill(L1.E());
		    h42->Fill(L2.E());
                    h111->Fill(L11.Px());
                    h121->Fill(L21.Px());
                    h211->Fill(L11.Py());
                    h221->Fill(L21.Py());
                    h311->Fill(L11.Pz());
                    h321->Fill(L21.Pz());
                    h411->Fill(L11.E());
                    h421->Fill(L21.E());  */      
    }}}}
 
    //Draw pictures.
  /*  h11->GetXaxis()->SetTitle("Px(GeV)");
    h11->GetYaxis()->SetTitle("Events");
    h11->Draw();
    c11->Print("c11.png");
    cout << endl;
    h12->GetXaxis()->SetTitle("Px(GeV)");
    h12->GetYaxis()->SetTitle("Events");
    h12->Draw();
    c11->Print("c12.png");
    cout << endl;
    h21->GetXaxis()->SetTitle("Py(GeV)");
    h21->GetYaxis()->SetTitle("Events");
    h21->Draw();
    c11->Print("c21.png");
    cout << endl;
    h22->GetXaxis()->SetTitle("Py(GeV)");
    h22->GetYaxis()->SetTitle("Events");
    h22->Draw();
    c11->Print("c22.png");
    cout << endl;
    h31->GetXaxis()->SetTitle("Pz(GeV)");
    h31->GetYaxis()->SetTitle("Events");
    h31->Draw();
    c11->Print("c31.png");
    cout << endl;
    h32->GetXaxis()->SetTitle("Pz(GeV)");
    h32->GetYaxis()->SetTitle("Events");
    h32->Draw();
    c11->Print("c32.png");
    cout << endl;
    h41->GetXaxis()->SetTitle("E(GeV)");
    h41->GetYaxis()->SetTitle("Events");
    h41->Draw();
    c11->Print("c41.png");
    cout << endl;
    h42->GetXaxis()->SetTitle("E(GeV)");
    h42->GetYaxis()->SetTitle("Events");
    h42->Draw();
    c11->Print("c42.png");
    cout << endl;
    h111->GetXaxis()->SetTitle("Px(GeV)");
    h111->GetYaxis()->SetTitle("Events");
    h111->Draw();
    c11->Print("c111.png");
    cout << endl;
    h121->GetXaxis()->SetTitle("Px(GeV)");
    h121->GetYaxis()->SetTitle("Events");
    h121->Draw();
    c11->Print("c121.png");
    cout << endl;
    h211->GetXaxis()->SetTitle("Py(GeV)");
    h211->GetYaxis()->SetTitle("Events");
    h211->Draw();
    c11->Print("c211.png");
    cout << endl;
    h221->GetXaxis()->SetTitle("Py(GeV)");
    h221->GetYaxis()->SetTitle("Events");
    h221->Draw();
    c11->Print("c221.png");
    cout << endl;
    h311->GetXaxis()->SetTitle("Pz(GeV)");
    h311->GetYaxis()->SetTitle("Events");
    h311->Draw();
    c11->Print("c311.png");
    cout << endl;
    h321->GetXaxis()->SetTitle("Pz(GeV)");
    h321->GetYaxis()->SetTitle("Events");
    h321->Draw();
    c11->Print("c321.png");
    cout << endl;
    h411->GetXaxis()->SetTitle("E(GeV)");
    h411->GetYaxis()->SetTitle("Events");
    h411->Draw();
    c11->Print("c411.png");
    cout << endl;
    h421->GetXaxis()->SetTitle("E(GeV)");
    h421->GetYaxis()->SetTitle("Events");
    h421->Draw();
    c11->Print("c421.png");
    cout << endl;
*/
  cout << endl;
  return 0;
}
