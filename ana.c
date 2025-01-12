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
  int npar=0;
  TLorentzVector PtEtaPhiEVector;
  vector<double> PE0, PE1, PE2, PPT0, PPT1, PPT2, PETA0, PETA1, PETA2, PPHI0, PPHI1, PPHI2;
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
  
  TBranch *b_p_E;
  TBranch *b_p_id;
  TBranch *b_p_Pt;
  TBranch *b_p_Eta;
  TBranch *b_p_Phi;
  TBranch *b_npar;
  TBranch *b_p_status;

  TFile *f = new TFile( filename.c_str() );
  TTree *my_tree = (TTree*)f->Get("truth");
 
  my_tree->SetBranchAddress("p_E", &p_E, &b_p_E );
  my_tree->SetBranchAddress("p_id", &p_id, &b_p_id );
  my_tree->SetBranchAddress("p_Pt", &p_Pt, &b_p_Pt );
  my_tree->SetBranchAddress("p_Eta", &p_Eta, &b_p_Eta );
  my_tree->SetBranchAddress("p_Phi", &p_Phi, &b_p_Phi );
  my_tree->SetBranchAddress("npar", &npar, &b_npar );
  my_tree->SetBranchAddress("p_status", &p_status, &b_p_status );
  
  TCanvas *c11 = new TCanvas();
  TH1I *h11 = new TH1I( "h11", "The Pt distribution of nu tau",30, 0, 3);
  TH1I *h12 = new TH1I( "h12", "The Pt distribution of pi", 30, 0, 3);
  TH1I *h21 = new TH1I( "h21", "The Eta distribution of nu tau", 40, -4, 4);
  TH1I *h22 = new TH1I( "h22", "The Eta distribution of pi", 40, -4, 4);
  TH1I *h31 = new TH1I( "h31", "The Phi distribution of nu tau", 40, -4, 4);
  TH1I *h32 = new TH1I( "h32", "The Phi distribution of pi", 40, -4, 4);
  TH1I *h41 = new TH1I( "h41", "The E distribution of nu tau", 35, 0, 3.5);
  TH1I *h42 = new TH1I( "h42", "The E distribution of pi", 35, 0, 3.5);

//Pick out three kinds of particles and put them in three vectors.
  N = (int) my_tree->GetEntries();
  for(int i=0; i<N+1; i++){
	my_tree->GetEntry( i );
        //cout << p_Pt->at(0) << "\n";
	// h1->Fill(p_E->at(0));
	for(int j = 0; j < npar; j++){
	switch (p_id->at(j) )
	{ 
	   case -15 :
	      if(p_status->at(j) == 2 ){
              PE0.push_back(p_E->at(j));
              PPT0.push_back(p_Pt->at(j));
              PETA0.push_back(p_Eta->at(j));
              PPHI0.push_back(p_Phi->at(j));}
	   case 16 :
	      if(p_status->at(j) == 1 ){
	      PE1.push_back(p_E->at(j));
              PPT1.push_back(p_Pt->at(j));
              PETA1.push_back(p_Eta->at(j));
              PPHI1.push_back(p_Phi->at(j));}
           case 211 :
	      if(p_status->at(j) == 1 ){
	      PE2.push_back(p_E->at(j));
              PPT2.push_back(p_Pt->at(j));
              PETA2.push_back(p_Eta->at(j));
              PPHI2.push_back(p_Phi->at(j));}
           default:
              continue;
	}    	
      }
//Pick out the partilles we need by using the conservation of momentum.
        int N0 = PPT0.size();
	int N1 = PPT1.size();
	int N2 = PPT2.size();
    	int MIN = N0 > N1 ? (N1 > N2 ? N2 : N1) : (N0 > N2 ? N2 : N0);
        double sum = 5000000;
        double sum1 = 0;
	if (N0 == MIN){
            
             for(int m=0; m<PPT0.size(); m++){
	          Pt0 = PPT0[i];
		  Eta0 = PETA0[i];
		  Phi0 = PPHI0[i];
		  E0 = PE0[i];
	          for(int i=0; i<PPT1.size(); i++){
                       Pt1 = PPT1[i];
	               Eta1 = PETA1[i];
                       Phi1 = PPHI1[i];
		       E1 = PE1[i];
	               for(int i=0; i<PPT2.size(); i++){
                            Pt2 = PPT2[i]; 
                            Eta2 = PETA2[i];
		            Phi2 = PPHI2[i];
			    E2 = PE2[i];
			    sum1 = abs(Pt1 + Pt2 - Pt0) + abs(Eta1 + Eta2 - Eta0) +abs(Phi1 + Phi2 -Phi0);
		            if (sum1 < sum){
				 sum = sum1;
				 PtEtaPhiEVector.SetPtEtaPhiE(Pt1,Eta1,Phi1,E1);
				 Pt11 = PtEtaPhiEVector.Pt();
                                 Eta11 = PtEtaPhiEVector.Eta();
                                 Phi11 = PtEtaPhiEVector.Phi();
                                 E11 = PtEtaPhiEVector.E();
				 PtEtaPhiEVector.SetPtEtaPhiE(Pt2,Eta2,Phi2,E2);
                                 Pt22 = PtEtaPhiEVector.Pt();
                                 Eta22 = PtEtaPhiEVector.Eta();
                                 Phi22 = PtEtaPhiEVector.Phi();
                                 E22 = PtEtaPhiEVector.E();}
                        }
	           }
	           h11->Fill(Pt11);
		   h12->Fill(Pt22);
		   h21->Fill(Eta11);
		   h22->Fill(Eta22);
		   h31->Fill(Phi11);
		   h32->Fill(Phi22);
		   h41->Fill(E11);
		   h42->Fill(E22);
                   int j = 0;
		   for(int i = 0; i<PPT1.size(); i++){
	                if(PPT1[i] != Pt11){
		             PPT1[j] = PPT1[i];
			     PETA1[j] = PETA1[i];
			     PPHI1[j] = PPHI1[i];
			     PE1[j] = PE1[i];
			     j++;} }
		   j = 0;
		   for(int i = 0; i<PPT2.size(); i++){
                        if(PPT2[i] != Pt22){
                             PPT2[j] = PPT2[i];
                             PETA2[j] = PETA2[i];
                             PPHI2[j] = PPHI2[i];
                             PE2[j] = PE2[i];
                             j++;} } 
	     }
	}
	sum = 5000000;
        if (N1 == MIN && N0 != MIN){

             for(int m=0; m<PPT1.size(); m++){
                  Pt1 = PPT1[i];
                  Eta1 = PETA1[i];
                  Phi1 = PPHI1[i];
                  E1 = PE1[i];
                  for(int i=0; i<PPT0.size(); i++){
                       Pt0 = PPT0[i];
                       Eta0 = PETA0[i];
                       Phi0 = PPHI0[i];
                       E0 = PE0[i];
                       for(int i=0; i<PPT2.size(); i++){
                            Pt2 = PPT2[i];
                            Eta2 = PETA2[i];
                            Phi2 = PPHI2[i];
                            E2 = PE2[i];
                            sum1 = abs(Pt1 + Pt2 - Pt0) + abs(Eta1 + Eta2 - Eta0) +abs(Phi1 + Phi2 -Phi0);
                            if (sum1 < sum){
                                 sum = sum1;
				 sum = sum1;
                                 PtEtaPhiEVector.SetPtEtaPhiE(Pt1,Eta1,Phi1,E1);
                                 Pt11 = PtEtaPhiEVector.Pt();
                                 Eta11 = PtEtaPhiEVector.Eta();
                                 Phi11 = PtEtaPhiEVector.Phi();
                                 E11 = PtEtaPhiEVector.E();
                                 PtEtaPhiEVector.SetPtEtaPhiE(Pt2,Eta2,Phi2,E2);
                                 Pt22 = PtEtaPhiEVector.Pt();
                                 Eta22 = PtEtaPhiEVector.Eta();
                                 Phi22 = PtEtaPhiEVector.Phi();
                                 E22 = PtEtaPhiEVector.E();}
                        }
                   }
                   h11->Fill(Pt11);
                   h12->Fill(Pt22);
                   h21->Fill(Eta11);
                   h22->Fill(Eta22);
                   h31->Fill(Phi11);
                   h32->Fill(Phi22);
                   h41->Fill(E11);
                   h42->Fill(E22);
                   int j = 0;
                   for(int i = 0; i<PPT1.size(); i++){
                        if(PPT1[i] != Pt11){
                             PPT1[j] = PPT1[i];
                             PETA1[j] = PETA1[i];
		             PPHI1[j] = PPHI1[i];
                             PE1[j] = PE1[i];
			     j++;} }
                   j = 0;
                   for(int i = 0; i<PPT2.size(); i++){
                        if(PPT2[i] != Pt22){
                             PPT2[j] = PPT2[i];
                             PETA2[j] = PETA2[i];
                             PPHI2[j] = PPHI2[i];
                             PE2[j] = PE2[i];
                             j++;} }
             }
        }
        if (N1 != MIN && N0 != MIN){

             for(int m=0; m<PPT2.size(); m++){
                  Pt2 = PPT2[i];
                  Eta2 = PETA2[i];
                  Phi2 = PPHI2[i];
                  E2 = PE2[i];
                  for(int i=0; i<PPT0.size(); i++){
                       Pt0 = PPT0[i];
                       Eta0 = PETA0[i];
                       Phi0 = PPHI0[i];
                       E0 = PE0[i];
                       for(int i=0; i<PPT1.size(); i++){
                            Pt1 = PPT1[i];
                            Eta1 = PETA1[i];
                            Phi1 = PPHI1[i];
                            E1 = PE1[i];
			    sum1 = abs(Pt1 + Pt2 - Pt0) + abs(Eta1 + Eta2 - Eta0) +abs(Phi1 + Phi2 -Phi0);
                            if (sum1 < sum){
                                 sum = sum1;
				 sum = sum1;
                                 PtEtaPhiEVector.SetPtEtaPhiE(Pt1,Eta1,Phi1,E1);
                                 Pt11 = PtEtaPhiEVector.Pt();
                                 Eta11 = PtEtaPhiEVector.Eta();
                                 Phi11 = PtEtaPhiEVector.Phi();
                                 E11 = PtEtaPhiEVector.E();
                                 PtEtaPhiEVector.SetPtEtaPhiE(Pt2,Eta2,Phi2,E2);
                                 Pt22 = PtEtaPhiEVector.Pt();
                                 Eta22 = PtEtaPhiEVector.Eta();
                                 Phi22 = PtEtaPhiEVector.Phi();
                                 E22 = PtEtaPhiEVector.E();}
                        }
                   }
                   h11->Fill(Pt11);
                   h12->Fill(Pt22);
                   h21->Fill(Eta11);
                   h22->Fill(Eta22);
                   h31->Fill(Phi11);
                   h32->Fill(Phi22);
                   h41->Fill(E11);
                   h42->Fill(E22);
                   int j = 0;
                   for(int i = 0; i<PPT1.size(); i++){
                        if(PPT1[i] != Pt11){
                             PPT1[j] = PPT1[i];
                             PETA1[j] = PETA1[i];
                             PPHI1[j] = PPHI1[i];
                             PE1[j] = PE1[i];
                             j++;} }
                   j = 0;
                   for(int i = 0; i<PPT2.size(); i++){
                        if(PPT2[i] != Pt22){
                             PPT2[j] = PPT2[i];
                             PETA2[j] = PETA2[i];
                             PPHI2[j] = PPHI2[i];
                             PE2[j] = PE2[i];
                             j++;} }
             }
        }
        
    PE0.clear(); PE1.clear(); PE2.clear();
    PPT0.clear(); PPT1.clear(); PPT2.clear();
    PETA0.clear(); PETA1.clear(); PETA2.clear();
    PPHI0.clear(); PPHI1.clear(); PPHI2.clear();
    }
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
   
  return 0;
}
