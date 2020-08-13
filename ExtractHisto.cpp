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
    h_name_h1.push_back("LRT_h_NCluster");
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

    h_name_h1.push_back("trigger/h_LRT_d0max_event");
    h_name_h1.push_back("trigger/h_FTF_d0max_event");
    h_name_h1.push_back("trigger/h_muonsig_d0max_event");
    h_name_h1.push_back("trigger/h_LRT_trigd_d0max_event");
    h_name_h1.push_back("trigger/h_muonprt_d0max_event");

    h_name_h1.push_back("muon/t_muoncounthistogram");
    h_name_h1.push_back("muon/t_multid0");
    h_name_h1.push_back("muon/t_multideta");
    h_name_h1.push_back("muon/t_multidphi");
    h_name_h1.push_back("muon/t_doubledeta");
    h_name_h1.push_back("muon/t_doubledphi");

    






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
    h_name_h2.push_back("LRT_h_drvTDLength");
    h_name_h2.push_back("LRT_h_d0vTDLength");
    h_name_h2.push_back("LRT_hNoCut_drvTDLength");
    
    
    h_name_h2.push_back("FTF_h_NPixvd0");
    h_name_h2.push_back("FTF_h_NPixvTDLength");
    h_name_h2.push_back("LRT_h_NPixvd0");
    h_name_h2.push_back("LRT_h_NPixvTDLength");
    h_name_h2.push_back("LRT_h_NSctvd0");
    h_name_h2.push_back("LRT_h_NSctvTDLength");
    h_name_h2.push_back("LRT_h_NClustervd0");
    h_name_h2.push_back("LRT_h_NClustervTDLength");

    h_name_h2.push_back("FTF_h_dzvz");
    h_name_h2.push_back("LRT_h_dzvz");

    h_name_h2.push_back("muon/t_multidetavdphi");
    h_name_h2.push_back("muon/t_doubledetavdphi");




   
    
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
        if(h_name_h2[j].EndsWith("drvTDLength")){
            gStyle->SetOptStat(111111);
        } else {
            gStyle->SetOptStat(1111);
        }

        
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

    h_eff_names.push_back("muon/h_d0eff_standalone");
    h_eff_names.push_back("muon/h_d0eff_combined");

    


        
    
    


    
    

    
    
    
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

    h_eff_xlabel.push_back("d0 (mm)");
    h_eff_xlabel.push_back("d0 (mm)");
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


    

    //Lmao time to superimpose offline with trigger
    
    bool eff1consistent;
    bool eff2consistent;
    TH1F *heffn_1;
    TH1F *heffd_1;
    TH1F *heffn_2;
    TH1F *heffd_2;
    TEfficiency* pEff1;
    TEfficiency* pEff2;
    TLegend *legend;
    TString efftitle1;
    TString efftitle2;

    std::vector<TString> t_eff_impose;
    t_eff_impose.push_back("trig_h_d0eff");
    t_eff_impose.push_back("trig_h_etaeff");
    t_eff_impose.push_back("trig_h_pTeff");
    t_eff_impose.push_back("trig_h_TDLeff");

    std::vector<TString> o_eff_impose;
    o_eff_impose.push_back("offl_h_d0eff");
    o_eff_impose.push_back("offl_h_etaeff");
    o_eff_impose.push_back("offl_h_pTeff");
    o_eff_impose.push_back("offl_h_TDLeff");
    
    std::vector<TString> to_eff_xlabel;
    to_eff_xlabel.push_back("d0 (mm)");
    to_eff_xlabel.push_back("eta");
    to_eff_xlabel.push_back("pT (GeV)");
    to_eff_xlabel.push_back("Transverse Decay Length (mm)");



    for(Int_t i = 0; i < t_eff_impose.size(); i++){
        heffn_1 = (TH1F*) f -> Get(t_eff_impose[i] + "_n");
        heffd_1 = (TH1F*) f -> Get(t_eff_impose[i] + "_d");
        heffn_2 = (TH1F*) f -> Get(o_eff_impose[i] + "_n");
        heffd_2 = (TH1F*) f -> Get(o_eff_impose[i] + "_d");

        if(TEfficiency::CheckConsistency(*heffn_1,*heffd_1)){
            eff1consistent = true;
            pEff1 = new TEfficiency(*heffn_1, *heffd_1);
            efftitle1 = "Efficiency of " + to_eff_xlabel[i] + "; " + to_eff_xlabel[i] + "; Efficiency";
        
        } else{
            eff1consistent = false;
            cout << t_eff_impose[i] << " histogram is not consistent" << endl;
        }

        if(TEfficiency::CheckConsistency(*heffn_2,*heffd_2)){
            eff2consistent = true;
            pEff2 = new TEfficiency(*heffn_2, *heffd_2);
            efftitle2 = "Efficiency of " + to_eff_xlabel[i] + "; " + to_eff_xlabel[i] + "; Efficiency";
            
        } else{
            eff2consistent = false;
            cout << o_eff_impose[i] << " histogram is not consistent" << endl;
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
            legend->AddEntry(pEff1, "Trigger");
            legend->AddEntry(pEff2,"Offline");
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
            c1 -> Print("plots/compare/" + t_eff_impose[i]+ "_v_"+o_eff_impose[i]+ ".png");
        }
    }

    


    



}       