#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Math/Vector4D.h>
int ana(){
 
  string dir = "/home/yrliu/Desktop/";
  string name = "events.hepmc.root";
  string filename = dir + name;
  
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
  //PE0--tao+,PE01--tao-,PE1--nutau+,PE11--nutao-;PE2--pi+,PE21--pi-
  vector<double> PE0, PE1, PE2, PE01, PE11, PE21, PPT0, PPT1, PPT2, PPT01, PPT11, PPT21, PETA0, PETA1, PETA2, PETA01, PETA11, PETA21, PPHI0, PPHI1, PPHI2, PPHI01, PPHI11, PPHI21, PVID0, PVID01, PCHILDVID1, PCHILDVID11, PCHILDVID2, PCHILDVID21;
  double Pt0 = 0;
  double Pt1 = 0;
  double Pt2 = 0;
  double Pt11 = 0;
  double Pt22 = 0;
  double Eta0 = 0;
  double Eta1 = 0;
  double Eta2 = 0;
  double Eta11 = 0;
  double Eta22 = 0;
  double Phi0 = 0;
  double Phi1 = 0;
  double Phi2 = 0;
  double Phi11 = 0;
  double Phi22 = 0;
  double E0 = 0;
  double E1 = 0;
  double E2 = 0;
  double E11 = 0;
  double E22 = 0;
  double diffmass = 0;
  int VID0 = 0;
  int ChildVID1 = 0;
  int ChildVID2 = 0;  
  
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

  TCanvas *c11 = new TCanvas();
  TH1I *h11 = new TH1I( "h11", "The Pt distribution of nu tau1",30, 0, 3);
  TH1I *h12 = new TH1I( "h12", "The Pt distribution of pi+", 30, 0, 3);
  TH1I *h21 = new TH1I( "h21", "The Eta distribution of nu tau1", 40, -4, 4);
  TH1I *h22 = new TH1I( "h22", "The Eta distribution of pi+", 40, -4, 4);
  TH1I *h31 = new TH1I( "h31", "The Phi distribution of nu tau1", 40, -4, 4);
  TH1I *h32 = new TH1I( "h32", "The Phi distribution of pi+", 40, -4, 4);
  TH1I *h41 = new TH1I( "h41", "The E distribution of nu tau1", 35, 0, 3.5);
  TH1I *h42 = new TH1I( "h42", "The E distribution of pi+", 35, 0, 3.5);
  TH1I *h5 = new TH1I( "h5", "The difference of M1 distribution",45, 0, 1);
  TH1I *h111 = new TH1I( "h111", "The Pt distribution of nu tau2",30, 0, 3);
  TH1I *h121 = new TH1I( "h121", "The Pt distribution of pi-", 30, 0, 3);
  TH1I *h211 = new TH1I( "h211", "The Eta distribution of nu tau2", 40, -4, 4);
  TH1I *h221 = new TH1I( "h221", "The Eta distribution of pi-", 40, -4, 4);
  TH1I *h311 = new TH1I( "h311", "The Phi distribution of nu tau2", 40, -4, 4);
  TH1I *h321 = new TH1I( "h321", "The Phi distribution of pi-", 40, -4, 4);
  TH1I *h411 = new TH1I( "h411", "The E distribution of nu tau2", 35, 0, 3.5);
  TH1I *h421 = new TH1I( "h421", "The E distribution of pi-", 35, 0, 3.5);
  TH1I *h51 = new TH1I( "h51", "The difference of M2 distribution",35, 0, 2);


//Pick out three kinds of particles and put them in six vectors.
  N = (int) my_tree->GetEntries();
  int sum1 = 0;
  int sum2 = 0;
  for(int i=0; i<N+1; i++){
	my_tree->GetEntry( i );
	if(npar == 6 ){
	for(int j = 0; j < npar; j++){
	switch (p_id->at(j) )
	{ 
	   case -15 :
	      if(p_status->at(j) == 2 ){
              PE0.push_back(p_E->at(j));
              PPT0.push_back(p_Pt->at(j));
              PETA0.push_back(p_Eta->at(j));
              PPHI0.push_back(p_Phi->at(j));
	      PVID0.push_back(p_VID->at(j));}
	      break;
	   case 16 :
	      if(p_status->at(j) == 1 ){
	      PE1.push_back(p_E->at(j));
              PPT1.push_back(p_Pt->at(j));
              PETA1.push_back(p_Eta->at(j));
              PPHI1.push_back(p_Phi->at(j));
	      PCHILDVID1.push_back(p_ChildVID->at(j));}
	      break;
           case 211 :
	      if(p_status->at(j) == 1 ){
	      PE2.push_back(p_E->at(j));
              PPT2.push_back(p_Pt->at(j));
              PETA2.push_back(p_Eta->at(j));
              PPHI2.push_back(p_Phi->at(j));
	      PCHILDVID2.push_back(p_ChildVID->at(j));}
	      break;
	   case 15 :
	      if(p_status->at(j) == 2 ){
              PE01.push_back(p_E->at(j));
              PPT01.push_back(p_Pt->at(j));
              PETA01.push_back(p_Eta->at(j));
              PPHI01.push_back(p_Phi->at(j));
	      PVID01.push_back(p_ChildVID->at(j));}
	      break;
	   case -16 :
              if(p_status->at(j) == 1 ){
              PE11.push_back(p_E->at(j));
              PPT11.push_back(p_Pt->at(j));
              PETA11.push_back(p_Eta->at(j));
              PPHI11.push_back(p_Phi->at(j));
              PCHILDVID11.push_back(p_ChildVID->at(j));}
              break;
	   case -211 :
	      if(p_status->at(j) == 1 ){
              PE21.push_back(p_E->at(j));
              PPT21.push_back(p_Pt->at(j));
              PETA21.push_back(p_Eta->at(j));
              PPHI21.push_back(p_Phi->at(j));
	      PCHILDVID21.push_back(p_ChildVID->at(j));}
	      break;
           default:
              goto endLoop;
	}    	
      }
	
        for(int m=0; m<PPT0.size(); m++){
	          Pt0 = PPT0[i];
		  Eta0 = PETA0[i];
		  Phi0 = PPHI0[i];
		  E0 = PE0[i];
		  VID0 = PVID0[i];
	          for(int i=0; i<PPT1.size(); i++){
                       Pt1 = PPT1[i];
	               Eta1 = PETA1[i];
                       Phi1 = PPHI1[i];
		       E1 = PE1[i];
		       ChildVID1 = PCHILDVID1[i];
	               for(int i=0; i<PPT2.size(); i++){
                            Pt2 = PPT2[i]; 
                            Eta2 = PETA2[i];
		            Phi2 = PPHI2[i];
			    E2 = PE2[i];
			    ChildVID2 = PCHILDVID2[i];
			    if (ChildVID1 == VID0 && ChildVID2 == VID0 ){
                 		 TLorentzVector PtEtaPhiEVector(0,0,0,0);
				 PtEtaPhiEVector.SetPtEtaPhiE(Pt0,Eta0,Phi0,E0);
				 TLorentzVector L1(0,0,0,0);
				 L1.SetPtEtaPhiE(Pt1,Eta1,Phi1,E1);
				 Pt11 = L1.Pt();
                                 Eta11 = L1.Eta();
                                 Phi11 = L1.Phi();
                                 E11 = L1.E();
				 TLorentzVector L2(0,0,0,0);
				 L2.SetPtEtaPhiE(Pt2,Eta2,Phi2,E2);
                                 Pt22 = L2.Pt();
                                 Eta22 = L2.Eta();
                                 Phi22 = L2.Phi();
                                 E22 = L2.E();
				 TLorentzVector diff = L1 + L2;
				 //diffmass = diff.Mag();
				 diffmass = L2.Mag();
				 h11->Fill(Pt11);
		                 h12->Fill(Pt22);
		                 h21->Fill(Eta11);
		                 h22->Fill(Eta22);
		                 h31->Fill(Phi11);
		                 h32->Fill(Phi22);
		                 h41->Fill(E11);
		                 h42->Fill(E22);
		                 h5->Fill(diffmass);
			         sum1+=1;}}}}
	for(int m=0; m<PPT01.size(); m++){
                  Pt0 = PPT01[i];
                  Eta0 = PETA01[i];
                  Phi0 = PPHI01[i];
                  E0 = PE01[i];
                  VID0 = PVID01[i];
                  for(int i=0; i<PPT11.size(); i++){
                       Pt1 = PPT11[i];
                       Eta1 = PETA11[i];
                       Phi1 = PPHI11[i];
                       E1 = PE11[i];
                       ChildVID1 = PCHILDVID11[i];
                       for(int i=0; i<PPT21.size(); i++){
                            Pt2 = PPT21[i];
                            Eta2 = PETA21[i];
                            Phi2 = PPHI21[i];
                            E2 = PE21[i];
                            ChildVID2 = PCHILDVID21[i];
                            if (ChildVID1 == VID0 && ChildVID2 == VID0 ){
                                 TLorentzVector PtEtaPhiEVector(0,0,0,0);
				 PtEtaPhiEVector.SetPtEtaPhiE(Pt0,Eta0,Phi0,E0);
                                 TLorentzVector L1(0,0,0,0);
                                 L1.SetPtEtaPhiE(Pt1,Eta1,Phi1,E1);
                                 Pt11 = L1.Pt();
                                 Eta11 = L1.Eta();
                                 Phi11 = L1.Phi();
                                 E11 = L1.E();
                                 TLorentzVector L2(0,0,0,0);
                                 L2.SetPtEtaPhiE(Pt2,Eta2,Phi2,E2);
                                 Pt22 = L2.Pt();
                                 Eta22 = L2.Eta();
                                 Phi22 = L2.Phi();
                                 E22 = L2.E();
                                 TLorentzVector diff = L1 + L2;
                                 //diffmass = diff.Mag();
				 diffmass  = L2.Mag();
                                 h111->Fill(Pt11);
                                 h121->Fill(Pt22);
                                 h211->Fill(Eta11);
                                 h221->Fill(Eta22);
                                 h311->Fill(Phi11);
                                 h321->Fill(Phi22);
                                 h411->Fill(E11);
                                 h421->Fill(E22);
                                 h51->Fill(diffmass);
                                 sum2+=1;}}}}
    endLoop:
              PE0.clear(); PE1.clear(); PE2.clear();
              PPT0.clear(); PPT1.clear(); PPT2.clear();
              PETA0.clear(); PETA1.clear(); PETA2.clear();
              PPHI0.clear(); PPHI1.clear(); PPHI2.clear();
              PE01.clear(); PE11.clear(); PE21.clear();
              PPT01.clear(); PPT11.clear(); PPT21.clear();
              PETA01.clear(); PETA11.clear(); PETA21.clear();
              PPHI01.clear(); PPHI11.clear(); PPHI21.clear();
              PVID0.clear(); PVID01.clear(); PCHILDVID1.clear();
              PCHILDVID11.clear(); PCHILDVID2.clear(); PCHILDVID21.clear();
    }}
    h11->GetXaxis()->SetTitle("Pt");
    h11->GetYaxis()->SetTitle("Events");
    h11->Draw();
    c11->Print("c11.png");
    cout << endl;
    h12->GetXaxis()->SetTitle("Pt");
    h12->GetYaxis()->SetTitle("Events");
    h12->Draw();
    c11->Print("c12.png");
    cout << endl;
    h21->GetXaxis()->SetTitle("Eta");
    h21->GetYaxis()->SetTitle("Events");
    h21->Draw();
    c11->Print("c21.png");
    cout << endl;
    h22->GetXaxis()->SetTitle("Eta");
    h22->GetYaxis()->SetTitle("Events");
    h22->Draw();
    c11->Print("c22.png");
    cout << endl;
    h31->GetXaxis()->SetTitle("Phi");
    h31->GetYaxis()->SetTitle("Events");
    h31->Draw();
    c11->Print("c31.png");
    cout << endl;
    h32->GetXaxis()->SetTitle("Phi");
    h32->GetYaxis()->SetTitle("Events");
    h32->Draw();
    c11->Print("c32.png");
    cout << endl;
    h41->GetXaxis()->SetTitle("E");
    h41->GetYaxis()->SetTitle("Events");
    h41->Draw();
    c11->Print("c41.png");
    cout << endl;
    h42->GetXaxis()->SetTitle("E");
    h42->GetYaxis()->SetTitle("Events");
    h42->Draw();
    c11->Print("c42.png");
    cout << endl;
    h5->GetXaxis()->SetTitle("Difference of M1");
    h5->GetYaxis()->SetTitle("Events");
    h5->Draw();
    c11->Print("c5.png");
    cout << endl;
    h111->GetXaxis()->SetTitle("Pt");
    h111->GetYaxis()->SetTitle("Events");
    h111->Draw();
    c11->Print("c111.png");
    cout << endl;
    h121->GetXaxis()->SetTitle("Pt");
    h121->GetYaxis()->SetTitle("Events");
    h121->Draw();
    c11->Print("c121.png");
    cout << endl;
    h211->GetXaxis()->SetTitle("Eta");
    h211->GetYaxis()->SetTitle("Events");
    h211->Draw();
    c11->Print("c211.png");
    cout << endl;
    h221->GetXaxis()->SetTitle("Eta");
    h221->GetYaxis()->SetTitle("Events");
    h221->Draw();
    c11->Print("c221.png");
    cout << endl;
    h311->GetXaxis()->SetTitle("Phi");
    h311->GetYaxis()->SetTitle("Events");
    h311->Draw();
    c11->Print("c311.png");
    cout << endl;
    h321->GetXaxis()->SetTitle("Phi");
    h321->GetYaxis()->SetTitle("Events");
    h321->Draw();
    c11->Print("c321.png");
    cout << endl;
    h411->GetXaxis()->SetTitle("E");
    h411->GetYaxis()->SetTitle("Events");
    h411->Draw();
    c11->Print("c411.png");
    cout << endl;
    h421->GetXaxis()->SetTitle("E");
    h421->GetYaxis()->SetTitle("Events");
    h421->Draw();
    c11->Print("c421.png");
    cout << endl;
    h51->GetXaxis()->SetTitle("Difference of M1");
    h51->GetYaxis()->SetTitle("Events");
    h51->Draw();
    c11->Print("c51.png");
    cout << endl;
    cout << "tau+pi+nutao:" << sum1 << endl;
    cout << "tau-pi-nutau:" << sum2 << endl;

  return 0;
}
