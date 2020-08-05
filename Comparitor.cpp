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
    TString filename1 = "roots/" + test1 + ".root";
    TString filename2 = "roots/" + test2 + ".root";

    //TString filename1 = "/scratch/baines/signal_tau1nsAOD/AOD.pool.root";
    //TString filename2 = "/scratch/baines/signal_tau1nsAOD/test13/AOD.pool.root";
    cout << filename1 << " " << filename2 << endl;

    TFile *f_1 = TFile::Open(filename1);
    TFile *f_2 = TFile::Open(filename2);
    
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050);
    gStyle -> SetOptStat(0);
    TLegend *legend;
    TString efftitle1;
    TString efftitle2;
    TString current_histo1; 
    TString current_histo2; 
    TH1F *h1_1;
    TH1F *h1_2;
    TH2F *h2_1;
    TH2F *h2_2;
    TH1F *heffn_1;
    TH1F *heffd_1;
    TH1F *heffn_2;
    TH1F *heffd_2;
    bool eff1consistent;
    bool eff2consistent;
    Double_t histmax;

    TEfficiency* pEff1;
    TEfficiency* pEff2;


    //Comparing standard stuff 
    //1d
    std::vector<TString> h_name_h1;
    h_name_h1.push_back("LRT_h_NPix");
    h_name_h1.push_back("LRT_h_NSct");
    h_name_h1.push_back("LRT_h_NCluster");
    h_name_h1.push_back("LRT_h_Nblayer");
    h_name_h1.push_back("LRT_h_NcontribPix");
    
    h_name_h1.push_back("LRT_h_d0fakes");
    h_name_h1.push_back("LRT_h_d0max_event");

    //TH1 Loop
    for(Int_t k = 0; k < h_name_h1.size(); k++){
        current_histo1 = h_name_h1[k];
        h1_1 = (TH1F*) f_1 -> Get(current_histo1);
        h1_2 = (TH1F*) f_2 -> Get(current_histo1);

        histmax = h1_1 -> GetMaximum();
        if(histmax < h1_2 -> GetMaximum()) histmax = h1_2 -> GetMaximum();

        h1_1 -> SetAxisRange(0, histmax + 10, "Y");
        h1_1->SetLineColor(2);
        
        h1_1->Draw();
        h1_2->Draw("same");

        legend = new TLegend(0.75,0.8,0.9,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h1_1, test1);
        legend->AddEntry(h1_2, test2);
        legend-> Draw("same"); 

        c1 -> Print("plots/" + test1+ "_v_"+test2 +"_"+ current_histo1 + ".png");
    }

    //2d
    std::vector<TString> h_name_h2;
    h_name_h2.push_back("LRT_h_NPixvd0");
    h_name_h2.push_back("LRT_h_NPixvTDLength");
    h_name_h2.push_back("LRT_h_NSctvd0");
    h_name_h2.push_back("LRT_h_NSctvTDLength");
    h_name_h2.push_back("LRT_h_NClustervd0");
    
    h_name_h2.push_back("LRT_h_NClustervTDLength");

    //TH1 Loop
    for(Int_t j = 0; j < h_name_h2.size(); j++){
        current_histo2 = h_name_h2[j];
        h2_1 = (TH2F*) f_1 -> Get(current_histo2);
        h2_2 = (TH2F*) f_2 -> Get(current_histo2);

        //histmax = h2_1 -> GetMaximum();
        //if(histmax < h2_2 -> GetMaximum()) histmax = h1_2 -> GetMaximum();

        //h1_1 -> SetAxisRange(0, histmax + 10, "Y");
        h2_1->SetLineColor(2);
        h2_1->Draw();
        h2_2->Draw("samebox");

        legend = new TLegend(0.75,0.8,0.9,0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h2_1, test1);
        legend->AddEntry(h2_2, test2);
        legend-> Draw("same"); 



        c1 -> Print("plots/" + test1+ "_v_"+test2 +"_"+ current_histo2 + ".png");
    }


    //List of Efficiency Histos to compare
    std::vector<TString> h_eff_names;
    h_eff_names.push_back("trig_h_d0eff");
    h_eff_names.push_back("trig_h_TDLeff");
    
    h_eff_names.push_back("FTF_h_d0eff");
    h_eff_names.push_back("LRT_h_d0eff");
    
    h_eff_names.push_back("trig_h_etaeff");
    h_eff_names.push_back("trig_h_pTeff");
    h_eff_names.push_back("offl_h_d0eff");
    h_eff_names.push_back("offl_h_etaeff");
    h_eff_names.push_back("offl_h_pTeff");
    
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
    h_eff_xlabel.push_back("d0 (mm)");
    h_eff_xlabel.push_back("Transverse Decay Length (mm)");
    
    h_eff_xlabel.push_back("d0 (mm) [FTF Only]");
    h_eff_xlabel.push_back("d0 (mm) [LRT Only]");
    
    h_eff_xlabel.push_back("eta");
    h_eff_xlabel.push_back("pT (GeV)");
    h_eff_xlabel.push_back("d0 (mm)");
    h_eff_xlabel.push_back("eta");
    h_eff_xlabel.push_back("pT (GeV)");
    
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
        heffn_1 = (TH1F*) f_1 -> Get(h_eff_names[i] + "_n");
        heffd_1 = (TH1F*) f_1 -> Get(h_eff_names[i] + "_d");
        heffn_2 = (TH1F*) f_2 -> Get(h_eff_names[i] + "_n");
        heffd_2 = (TH1F*) f_2 -> Get(h_eff_names[i] + "_d");
        //auto gr1 = new TGraphAsymmErrors(heffn_1, heffd_1);
        //auto gr2 = new TGraphAsymmErrors(heffn_2, heffd_2);
        //gr2 ->SetLineColor(2);

        //cout << h_eff_names[i] + "_n" << endl;
        if(TEfficiency::CheckConsistency(*heffn_1,*heffd_1)){
            eff1consistent = true;
            pEff1 = new TEfficiency(*heffn_1, *heffd_1);
            efftitle1 = "Efficiency of " + h_eff_xlabel[i] + "; " + h_eff_xlabel[i] + "; Efficiency";
            
        } else{
            eff1consistent = false;
            cout << h_eff_names[i] << " histogram is not consistent" << endl;
        }

        if(TEfficiency::CheckConsistency(*heffn_2,*heffd_2)){
            eff2consistent = true;
            pEff2 = new TEfficiency(*heffn_2, *heffd_2);
            efftitle2 = "Efficiency of " + h_eff_xlabel[i] + "; " + h_eff_xlabel[i] + "; Efficiency";
            
        } else{
            eff2consistent = false;
            cout << h_eff_names[i] << " histogram is not consistent" << endl;
        }

        if (eff1consistent && eff2consistent){
            ///pEff1 properties
            pEff1->SetTitle(efftitle1); 
            pEff1->SetLineColor(2);
            pEff1->SetMarkerStyle(23);
            pEff1->SetMarkerColor(2);
            pEff1->SetMarkerSize(1);
            

            ///pEff2 properties
            pEff2->SetTitle(efftitle2); 
            pEff2->SetLineColor(4);
            pEff2->SetMarkerStyle(26);
            pEff2->SetMarkerColor(4);
            pEff2->SetMarkerSize(1);
   
            //Legends
            legend = new TLegend(0.75,0.8,0.9,0.9);
            legend->SetHeader("Efficiency Markers","C"); // option "C" allows to center the header
            legend->AddEntry(pEff1, test1);
            legend->AddEntry(pEff2,test2);
            //legend->AddEntry("gr","Graph with error bars","lep");
   
            //Draw          
            pEff1 -> Draw();


            gPad->Update(); 
            auto graph = pEff1->GetPaintedGraph(); 
            graph->SetMinimum(0);
            graph->SetMaximum(1); 
            gPad->Update(); 

            pEff2 -> Draw("same");
            legend-> Draw("same"); 
            

            //gr1 -> Draw("AP0");
            //gr2 -> Draw("same");
            c1 -> Print("plots/" + test1+ "_v_"+test2 +"_"+ h_eff_names[i] + ".png");
        }
    }
}
