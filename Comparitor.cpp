#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"


void Comparitor(TString test1, TString test2){
    TString filename1 = test1 + ".root";
    TString filename2 = test2 + ".root";
    cout << filename1 << " " << filename2 << endl;

    TFile *f1 = TFile::Open(filename1);
    TFile *f2 = TFile::Open(filename2);
    
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050);
    TLegend *legend;
    TString efftitle1;
    TString efftitle2;
    TString current_histo1; 
    TString current_histo2; 
    TH1F *h1;
    TH2F *h2;
    TH1F *h1_effn;
    TH1F *h1_effd;
    TH1F *h2_effn;
    TH1F *h2_effd;
    bool eff1consistent;
    bool eff2consistent;

    TEfficiency* pEff1;
    TEfficiency* pEff2;


    //Comparing standard stuff 
    //1d
    std::vector<TString> h_name_h1;

    //2d
    std::vector<TString> h_name_h2;


    //List of Efficiency Histos to compare
    std::vector<TString> h_eff_names;
    h_eff_names.push_back("FTF_h_d0eff");
    h_eff_names.push_back("LRT_h_d0eff");
    h_eff_names.push_back("trig_h_d0eff");
    h_eff_names.push_back("trig_h_etaeff");
    h_eff_names.push_back("trig_h_pTeff");
    h_eff_names.push_back("offl_h_d0eff");
    h_eff_names.push_back("offl_h_etaeff");
    h_eff_names.push_back("offl_h_pTeff");
    h_eff_names.push_back("trig_h_TDLeff");
    h_eff_names.push_back("tgof_h_d0eff");
    h_eff_names.push_back("trig_h_z0eff");
    h_eff_names.push_back("trig_h_pTefflog");
    h_eff_names.push_back("offl_h_z0eff");
    h_eff_names.push_back("offl_h_TDLeff");
    h_eff_names.push_back("tgof_h_etaeff");
    h_eff_names.push_back("tgof_h_pTeff");
    h_eff_names.push_back("tgof_h_z0eff");
    h_eff_names.push_back("tgof_h_TDLeff");
    

    //List of Lables
    std::vector<TString> h_eff_xlabel;
    h_eff_xlabel.push_back("d0 (mm) [FTF Only]");
    h_eff_xlabel.push_back("d0 (mm) [LRT Only]");
    h_eff_xlabel.push_back("d0 (mm)");
    h_eff_xlabel.push_back("eta");
    h_eff_xlabel.push_back("pT (GeV)");
    h_eff_xlabel.push_back("d0 (mm)");
    h_eff_xlabel.push_back("eta");
    h_eff_xlabel.push_back("pT (GeV)");
    h_eff_xlabel.push_back("Transverse Decay Length (mm)");
    h_eff_xlabel.push_back("d0 (mm) [respect to offline]");
    h_eff_xlabel.push_back("z0 (mm)");
    h_eff_xlabel.push_back("pT (GeV) [Log]");
    h_eff_xlabel.push_back("z0 (mm)");
    h_eff_xlabel.push_back("Transverse Decay Length (mm)");
    h_eff_xlabel.push_back("eta");
    h_eff_xlabel.push_back("pT (GeV)");
    h_eff_xlabel.push_back("z0 (mm)");
    h_eff_xlabel.push_back("Transverse Decay Length (mm)");
    
    //Comparing for Efficiencies
    for(Int_t i = 0; i < h_eff_names.size(); i++){
        h1_effn = (TH1F*) f1 -> Get(h_eff_names[i] + "_n");
        h1_effd = (TH1F*) f1 -> Get(h_eff_names[i] + "_d");
        h2_effn = (TH1F*) f2 -> Get(h_eff_names[i] + "_n");
        h2_effd = (TH1F*) f2 -> Get(h_eff_names[i] + "_d");

        //cout << h_eff_names[i] + "_n" << endl;
        if(TEfficiency::CheckConsistency(*h1_effn,*h1_effd)){
            eff1consistent = true;
            pEff1 = new TEfficiency(*h1_effn, *h1_effd);
            efftitle1 = "Efficiency of " + h_eff_xlabel[i] + "; " + h_eff_xlabel[i] + "; Efficiency";
        } else{
            eff1consistent = false;
            cout << h_eff_names[i] << " histogram is not consistent" << endl;
        }

        if(TEfficiency::CheckConsistency(*h1_effn,*h1_effd)){
            eff2consistent = true;
            pEff2 = new TEfficiency(*h1_effn, *h1_effd);
            efftitle2 = "Efficiency of " + h_eff_xlabel[i] + "; " + h_eff_xlabel[i] + "; Efficiency";
            
        } else{
            eff2consistent = false;
            cout << h_eff_names[i] << " histogram is not consistent" << endl;
        }

        if (eff1consistent && eff2consistent){
            ///pEff1 properties
            pEff1->SetTitle(efftitle1); 
            pEff1->SetLineColor(2);
            pEff1->SetMarkerStyle(40);
            pEff1->SetMarkerColor(2);
            pEff1->SetMarkerSize(1.5);
            

            ///pEff2 properties
            pEff2->SetTitle(efftitle2); 
            pEff2->SetLineColor(4);
            pEff2->SetMarkerStyle(44);
            pEff2->SetMarkerColor(4);
            pEff2->SetMarkerSize(1.5);
   
            //Legends
            legend = new TLegend(0.7,0.8,0.9,0.9);
            legend->SetHeader(test1 + " v " + test2 + " Efficiency Markers","C"); // option "C" allows to center the header
            legend->AddEntry(pEff1, test1);
            legend->AddEntry(pEff2,test2);
            //legend->AddEntry("gr","Graph with error bars","lep");
   
            //Draw
            pEff1 -> Draw("APE0");
            pEff2 -> Draw("same");
            legend-> Draw("same");
            c1 -> Print("plots/" + test1+ "_v_"+test2 +"_"+ h_eff_names[i] + ".png");
        }
    }
}
