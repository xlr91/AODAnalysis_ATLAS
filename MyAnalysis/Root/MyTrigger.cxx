#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyTrigger.h>
#include <xAODEventInfo/EventInfo.h>
#include <cmath>
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TMath.h"



MyTrigger :: MyTrigger (const std::string& name,
                                  ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm (name, pSvcLocator)
{
    declareProperty( "d0min" , m_d0min  = 50, "Minimum value of d0");
    declareProperty( "dphimin" , m_dphimin  = 1, "Minimum value of d0");
    declareProperty( "dpt" , m_dpt  = 1000, "Minimum value of d0");
}

bool MyTrigger::triggercond1(const xAOD::TrackParticle* LRT_track, Float_t t_d0min){
    bool pass = false;
    if(LRT_track -> d0() > t_d0min){
        pass = true;
    }
    return pass;
}

bool MyTrigger::triggercond2(const xAOD::TrackParticle* LRT_track1, const xAOD::TrackParticle* LRT_track2, Float_t t_d0min, Float_t t_dphimin, Float_t t_dpt){
    bool pass = false;
    if(LRT_track1 -> d0() > t_d0min && LRT_track2 -> d0() > t_d0min){
        if(abs(LRT_track1 -> phi() - LRT_track2 -> phi()) > t_dphimin){
            pass = true;
        }
    }
    return pass;
}


StatusCode MyTrigger :: initialize ()
{
    //store all the histograms under MyTrigger/
    ANA_MSG_INFO("Initialize trigger");
    ANA_CHECK(book(TH1F("MyTrigger/h_truthpass", "Event has or not rhadron to 2muon", 2, 0, 2)));
    ANA_CHECK(book(TH1F("MyTrigger/h_trackpass", "Event passed the trigger cuts", 2, 0, 2)));
    ANA_CHECK(book(TH1F("MyTrigger/h_passclassification", "0: no. total events, 1: no. truth pass, 2: no. track pass, 3: both truth and track pass", 4, 0, 4)));

    ANA_CHECK(book(TH1F("MyTrigger/t_eff_n", "Efficiency_function_of_d0_n", 1, 0, 1)));
    ANA_CHECK(book(TH1F("MyTrigger/t_eff_d", "Efficiency_function_of_d0_n", 1, 0, 1)));
    /*
    ANA_CHECK(book(TH1F("MyTrigger/t_eff_d0_n", "Efficiency_function_of_d0_n", 25, -300, 300)));
    ANA_CHECK(book(TH1F("MyTrigger/t_eff_d0_d", "Efficiency_function_of_d0_n", 25, -300, 300)));
    ANA_CHECK(book(TH1F("MyTrigger/t_eff_d0", "Efficiency_function_of_d0_n", 25, -300, 300)));
    ANA_CHECK(book(TH1F("MyTrigger/t_eff_n", "Efficiency_function_of_d0_n", 1, 0, 1)));
    */
    return StatusCode::SUCCESS;
}



StatusCode MyTrigger :: execute ()
{
    ANA_MSG_INFO("Executing Trigger");
    
    int muoncount;
    bool truthpass = false;
    bool trackpass = false;
    bool trigger_cuts;
    bool trigger_cuts2;


    //here we have the truth loop
    //here we need to check if they decay into the two muons thing
    //its similar to the thing one that you have on the analysis one

    const xAOD::TruthParticleContainer* truthparticles = nullptr;
    const xAOD::TruthParticle* child = nullptr; 
    const xAOD::TruthParticle* gchild = nullptr; 
    const xAOD::TrackParticleContainer* trig_LRTparticles = nullptr;

    ANA_CHECK (evtStore() -> retrieve (truthparticles, "TruthParticles"));
    ANA_CHECK (evtStore() -> retrieve (trig_LRTparticles, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_FullScanlrt_FTF"));

    //truth loop
    for (const xAOD::TruthParticle* truth : *truthparticles) {
        if (truth -> absPdgId() == 1000006 && truth -> nChildren() > 1){
            for (size_t ichild=0; ichild< truth->nChildren() ; ichild++) {
                child =  truth ->child(ichild);
                muoncount = 0;
                for (size_t igchild=0; igchild< child->nChildren() ; igchild++) {
                    gchild = child -> child(igchild);
                    if (gchild == nullptr) continue;
                    
                    if(gchild -> absPdgId() == 13){
                        muoncount = muoncount + 1; 
                    }

                }
                ANA_MSG_INFO("muoncounttriggr " << muoncount);
                if (muoncount == 2) truthpass = true;
            }
        }
    }
    
    for (const xAOD::TrackParticle* LRT_T : *trig_LRTparticles){
        trigger_cuts = triggercond1(LRT_T, m_d0min);
        if (trigger_cuts){
            for (const xAOD::TrackParticle* LRT_T_2 : *trig_LRTparticles){
                trigger_cuts2 = triggercond2(LRT_T, LRT_T_2, m_d0min, m_dphimin, m_dpt);
                if(trigger_cuts2){
                    trackpass = true;         
                }
            }
        }
    }
    
    hist("MyTrigger/h_passclassification") -> Fill(0);

    

    if(truthpass){
        hist("MyTrigger/h_passclassification") -> Fill(1);
        hist("MyTrigger/t_eff_d") -> Fill(0);
    }
    if(trackpass){
        hist("MyTrigger/h_passclassification") -> Fill(2);
    }

    if(truthpass && trackpass){
        hist("MyTrigger/h_passclassification") -> Fill(3);
        hist("MyTrigger/t_eff_n") -> Fill(0);
    }
    
    //if bool has been triggered, 'record the event'
    
    ANA_MSG_INFO("truth event: " << truthpass);
    ANA_MSG_INFO("track event: " << trackpass);
    hist("MyTrigger/h_truthpass") -> Fill(truthpass);
    hist("MyTrigger/h_trackpass") -> Fill(trackpass);
    return StatusCode::SUCCESS;
}



StatusCode MyTrigger :: finalize ()
{
    ANA_MSG_INFO("ayy it works lmao");
    return StatusCode::SUCCESS;
}