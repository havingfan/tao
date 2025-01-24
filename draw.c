#include <iostream>
#include <fstream>

using namespace std;

int draw(){

        ifstream read("nutao1px.txt");

	int N = 5089;
	double a = 1;

	TCanvas *c1 = new TCanvas();
        TH1I *h1 = new TH1I( "h11", "The Px distribution of nu tau1", 40, -1, 1);
        
	for( int i = 0; i < N; i++ ){
		read >> a;
		h1->Fill(a);
	}
        
	h1->GetXaxis()->SetTitle("Px(GeV)");
        h1->GetYaxis()->SetTitle("Events");
	h1->Draw();
        c1->Print("nutao1px.png");

	return 0;

}

