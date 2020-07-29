#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"


void Prettify(TString test1){

    TString filename1 = test1 + ".root";
    TFile *f1 = TFile::Open(filename1);
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050);
    TString current_histo1; 
    TH2F *h2;


    h2 = (TH2F*) f1 -> Get("LRT_h_NPixvTDLength");
    h2 -> Draw();
    Float_t ymax = h2 -> GetMaximum();
    cout << ymax << endl;
    
    TLine *line1 = new TLine(33.25,0,33.25,20);
    line1->SetLineColor(kRed);
    line1->Draw("same");
    TLine *line2 = new TLine(50.5,0,50.5,20);
    line2->SetLineColor(kRed);
    line2->Draw("same");
    TLine *line3 = new TLine(88.5,0,88.5,20);
    line3->SetLineColor(kRed);
    line3->Draw("same");
    TLine *line4 = new TLine(122.5,0,122.5,20);
    line4->SetLineColor(kRed);
    line4->Draw("same");

    c1 -> Print("plots/testing.png");


}
    