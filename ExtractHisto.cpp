#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"

//test

void ExtractHisto(){
    TString filename = "MyxAODAnalysis.outputs.root";
    TFile *f = TFile::Open(filename);
    TCanvas *c1 = new TCanvas("c1"," Efficiency ",50,50,1680,1050);
    TString current_histo; 
    TH1F *h1;
    TH2F *h2;

    TH1F *h1_effn;
    TH1F *h1_effd;



    std::vector<TString> h_name_h1;

    h_name_h1.push_back("h_childDecayLength");


    h_name_h1.push_back("h_offpass");
    h_name_h1.push_back("h_ftfpass");
    h_name_h1.push_back("h_lrtpass");

    
    h_name_h1.push_back("truth_h_dr");
    h_name_h1.push_back("truth_h_d0");
    h_name_h1.push_back("truth_h_eta");
    h_name_h1.push_back("truth_h_pT");
    h_name_h1.push_back("truth_h_phi");
    h_name_h1.push_back("truth_h_dphi");
    h_name_h1.push_back("truth_h_deta");


    h_name_h1.push_back("FTF_h_dr");
    h_name_h1.push_back("FTF_h_dr_0-10");
    h_name_h1.push_back("FTF_h_dr_20-30");
    h_name_h1.push_back("FTF_h_d0");
    h_name_h1.push_back("FTF_h_z0");
    h_name_h1.push_back("FTF_h_dz");
    h_name_h1.push_back("FTF_h_dz_0-10");
    h_name_h1.push_back("FTF_h_dz_20-30");
    h_name_h1.push_back("FTF_h_eta");
    h_name_h1.push_back("FTF_h_pT");
    h_name_h1.push_back("FTF_h_phi");
    h_name_h1.push_back("FTF_h_dphi");
    h_name_h1.push_back("FTF_h_dphi_0-10");
    h_name_h1.push_back("FTF_h_dphi_20-30");
    h_name_h1.push_back("FTF_h_deta");

    h_name_h1.push_back("FTF_h_NPix");
    h_name_h1.push_back("FTF_h_NSct");
    h_name_h1.push_back("FTF_h_Nblayer");
    h_name_h1.push_back("FTF_h_NcontribPix");


    h_name_h1.push_back("LRT_h_dr");
    h_name_h1.push_back("LRT_h_dr_0-10");
    h_name_h1.push_back("LRT_h_dr_20-30");
    h_name_h1.push_back("LRT_h_d0");
    h_name_h1.push_back("LRT_h_z0");
    h_name_h1.push_back("LRT_h_dz");
    h_name_h1.push_back("LRT_h_dz_0-10");
    h_name_h1.push_back("LRT_h_dz_20-30");
    h_name_h1.push_back("LRT_h_eta");
    h_name_h1.push_back("LRT_h_pT");
    h_name_h1.push_back("LRT_h_phi");
    h_name_h1.push_back("LRT_h_dphi");
    h_name_h1.push_back("LRT_h_dphi_0-10");
    h_name_h1.push_back("LRT_h_dphi_20-30");
    h_name_h1.push_back("LRT_h_deta");

    h_name_h1.push_back("LRT_h_NPix");
    h_name_h1.push_back("LRT_h_NSct");
    h_name_h1.push_back("LRT_h_Nblayer");
    h_name_h1.push_back("LRT_h_NcontribPix");


    h_name_h1.push_back("offline_h_dr");
    h_name_h1.push_back("offline_h_d0");
    h_name_h1.push_back("offline_h_dz");
    h_name_h1.push_back("offline_h_eta");
    h_name_h1.push_back("offline_h_pT");
    h_name_h1.push_back("offline_h_phi");
    h_name_h1.push_back("offline_h_dphi");
    h_name_h1.push_back("offline_h_deta");

    h_name_h1.push_back("FTF_h_d0fakes_n");
    h_name_h1.push_back("FTF_h_d0fakes_d");
    h_name_h1.push_back("FTF_h_d0fakes");

    h_name_h1.push_back("LRT_h_d0fakes_n");
    h_name_h1.push_back("LRT_h_d0fakes_d");
    h_name_h1.push_back("LRT_h_d0fakes");

    h_name_h1.push_back("FTF_h_etafakes");
    h_name_h1.push_back("LRT_h_etafakes");



    /*
    h_name_h1.push_back("trig_h_d0eff");
    h_name_h1.push_back("trig_h_etaeff");
    h_name_h1.push_back("trig_h_pTeff");

    h_name_h1.push_back("offl_h_d0eff");
    h_name_h1.push_back("offl_h_etaeff");
    h_name_h1.push_back("offl_h_pTeff");
    */
    
    std::vector<TString> h_name_h1_log;

    
    
    h_name_h1_log.push_back("FTF_h_dr");
    h_name_h1_log.push_back("FTF_h_dr_0-10");
    h_name_h1_log.push_back("FTF_h_dr_20-30");
    h_name_h1_log.push_back("FTF_h_dz");
    h_name_h1_log.push_back("FTF_h_dz_0-10");
    h_name_h1_log.push_back("FTF_h_dz_20-30");
    h_name_h1_log.push_back("FTF_h_dphi");
    h_name_h1_log.push_back("FTF_h_dphi_0-10");
    h_name_h1_log.push_back("FTF_h_dphi_20-30");
    h_name_h1_log.push_back("FTF_h_deta");



    h_name_h1_log.push_back("LRT_h_dr");
    h_name_h1_log.push_back("LRT_h_dr_0-10");
    h_name_h1_log.push_back("LRT_h_dr_20-30");
    h_name_h1_log.push_back("LRT_h_dz");
    h_name_h1_log.push_back("LRT_h_dz_0-10");
    h_name_h1_log.push_back("LRT_h_dz_20-30");
    h_name_h1_log.push_back("LRT_h_dphi");
    h_name_h1_log.push_back("LRT_h_dphi_0-10");
    h_name_h1_log.push_back("LRT_h_dphi_20-30");
    h_name_h1_log.push_back("LRT_h_deta");


    std::vector<TString> h_name_h2;
    h_name_h2.push_back("compare/h_d0truthvtrack_offline");
    h_name_h2.push_back("compare/h_d0truthvtrack_FTF");
    h_name_h2.push_back("compare/h_d0truthvtrack_LRT");
    h_name_h2.push_back("compare/h_d0truthvtrack_All");



    h_name_h2.push_back("compare/h_d0diff_offline");
    h_name_h2.push_back("compare/h_d0diff_FTF");
    h_name_h2.push_back("compare/h_d0diff_LRT");
    
    h_name_h2.push_back("offline_h_phivTDLength");
    h_name_h2.push_back("offline_h_dphivTDLength");
    h_name_h2.push_back("truth_h_phivTDLength");
    h_name_h2.push_back("truth_h_dphivTDLength");
    h_name_h2.push_back("FTF_h_phivTDLength");
    h_name_h2.push_back("FTF_h_dphivTDLength");
    h_name_h2.push_back("LRT_h_phivTDLength");
    h_name_h2.push_back("LRT_h_dphivTDLength");
    
    h_name_h2.push_back("FTF_h_NPixvd0");
    h_name_h2.push_back("FTF_h_NPixvTDLength");
    h_name_h2.push_back("LRT_h_NPixvd0");
    h_name_h2.push_back("LRT_h_NPixvTDLength");

    h_name_h2.push_back("FTF_h_dzvz");
    h_name_h2.push_back("LRT_h_dzvz");




   
    
    //TH1 Loop
    for(Int_t k = 0; k < h_name_h1.size(); k++){
        current_histo = h_name_h1[k];
        h1 = (TH1F*) f -> Get(current_histo);

        h1->Draw();
        c1 -> Print("plots/" + current_histo + ".png");
    }
    
    //TH1 Log Loop
    gPad->SetLogy(1);

    for(Int_t k1 = 0; k1 < h_name_h1_log.size(); k1++){
        current_histo = h_name_h1_log[k1];
        h1 = (TH1F*) f -> Get(current_histo);

        h1->Draw();
        c1 -> Print("plots/" + current_histo + "_log.png");
    }

    gPad->SetLogy(0);
    //TH2Loop
    for(Int_t j = 0; j < h_name_h2.size(); j++){
        current_histo = h_name_h2[j];
        h2 = (TH2F*) f -> Get(current_histo);

        h2->Draw();
        c1 -> Print("plots/" + current_histo + ".png");
    }
    
    

    ///TEFFICIENCY LETS GO 

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

    TString efftitle;

    

    for(Int_t i = 0; i < h_eff_names.size(); i++){
        h1_effn = (TH1F*) f -> Get(h_eff_names[i] + "_n");
        h1_effd = (TH1F*) f -> Get(h_eff_names[i] + "_d");

        //cout << h_eff_names[i] + "_n" << endl;
        if(TEfficiency::CheckConsistency(*h1_effn,*h1_effd)){
            ///cout << "Histograms Consistent " << endl;
            TEfficiency* pEff = new TEfficiency(*h1_effn, *h1_effd);
            efftitle = "Efficiency of " + h_eff_xlabel[i] + "; " + h_eff_xlabel[i] + "; Efficiency";
            pEff->SetTitle(efftitle); 
            pEff->SetLineColor(4);
            pEff->SetMarkerStyle(20);
            pEff->SetMarkerColor(1);
            pEff -> Draw("APE0");
            if (h_eff_names[i].EndsWith("log")) {
                gPad->SetLogx(1);
            } else{
                gPad->SetLogx(0);
            }
            //pEff->GetPaintedGraph()->GetXaxis()->SetTitle("x3");
            c1 -> Print("plots/effs/" + h_eff_names[i] + ".png");
        } else{
            cout << h_eff_names[i] << " histogram is not consistent" << endl;
        }
    }





    

    



}       