

#Job Properties
test = 14 #3
J_CUTNUMBER = 2 #2
DRCUTnumber = 0.01
DZCUTnumber = 150
ETACUTnumber = 0.1
Mflag = False

#original AOD file to test on (Offline only)
if test == 0:
    testFile = "/afs/cern.ch/user/e/emuhamma/public/runtest/AOD.pool.root"
    OFlag = True
    Tflag = True

#original AOD file to test on (Offline only)
if test == 1:
    testFile = "/eos/atlas/atlasdatadisk/rucio/mc16_13TeV/05/e6/AOD.19016897._000001.pool.root.1"
    OFlag = True
    Tflag = False
    
#original MET file to test on (trigger only)
if test == 2:
    testFile = "/afs/cern.ch/user/b/baines/work/public/LRT/METtest/AOD.pool.root"
    OFlag = False
    Tflag = True
    
#50 mm d0
if test == 3:
    testFile = "/scratch/baines/signal_tau1nsAOD/AOD.pool.root" 
    OFlag = True
    Tflag = True

#300 mm d0 (can chagne)
if test == 4:
    testFile = "/scratch/baines/signal_tau1nsAOD/test13/AOD.pool.root"
    OFlag = False
    Tflag = True

#muons with PU 40 with lrttest=4 (50mm d0Max)
if test == 5:
    testFile = "/scratch/baines/muonaod/AOD.pool.root"
    OFlag = False
    Tflag = True

#Test14
if test == 6:
    testFile = "/scratch/emuhamma/testAOD/test14/AOD.pool.root"
    OFlag = True
    Tflag = True

#Test15
if test == 7:
    testFile = "/scratch/emuhamma/testAOD/test15/AOD.pool.root"
    OFlag = True
    Tflag = True

#Test16
if test == 8:
    testFile = "/scratch/emuhamma/testAOD/test16/AOD.pool.root"
    OFlag = True
    Tflag = True

#Test17
if test == 9:
    testFile = "/scratch/emuhamma/testAOD/test17/AOD.pool.root"
    OFlag = True
    Tflag = True

#Test18-21
if test == 10:
    testFile = "/scratch/emuhamma/testAOD/test21/AOD.pool.root"
    OFlag = True
    Tflag = True

#Test22-26
if test == 11:
    testFile = "/scratch/emuhamma/testAOD/test22/AOD.pool.root"
    OFlag = True
    Tflag = True

#John's new thing
if test == 12:
    testFile = "/scratch/baines/mrtest/AOD.pool.root"
    OFlag = False
    Tflag = True

#multiple signal samples
if test == 13:
    testFile = ["/scratch/emuhamma/test21/signalstat/signalstat1/AOD.pool.root", 
                "/scratch/emuhamma/test21/signalstat/signalstat2/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat3/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat4/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat5/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat6/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat7/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat8/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat9/AOD.pool.root",
                "/scratch/emuhamma/test21/signalstat/signalstat10/AOD.pool.root"]
    OFlag = True
    Tflag = True

#signal making trigger
if test == 14:  
    OFlag = True
    Tflag = True
    Mflag = True
    testFile = "/scratch/emuhamma/test21/signalmuon/AOD.pool.root"

 #prompt making trigger
if test == 15:  
    OFlag = True
    Tflag = True
    testFile = "/scratch/emuhamma/test21/prompt/AOD.pool.root"

#signal with muons  trigger
if test == 16:  
    OFlag = True
    Tflag = False
    Mflag = True
    testFile = "/eos/atlas/atlasdatadisk/rucio/mc16_13TeV/05/e6/AOD.19016897._000001.pool.root.1"

#See: https://twiki.cern.ch/twiki/bin/viewauth/AtlasComputing/SoftwareTutorialxAODAnalysisInCMake for more details about anything here


#override next line on command line with: --filesInput=XXX
jps.AthenaCommonFlags.FilesInput = [testFile] 

if type(testFile) is list:
    jps.AthenaCommonFlags.FilesInput = testFile


#Specify AccessMode (read mode) ... ClassAccess is good default for xAOD
jps.AthenaCommonFlags.AccessMode = "ClassAccess" 

jps.AthenaCommonFlags.HistOutputs = ["ANALYSIS:MyxAODAnalysis.outputs.root"]
svcMgr.THistSvc.MaxFileSize=-1 #speeds up jobs that output lots of histograms


# Create the algorithm's configuration.
from AnaAlgorithm.DualUseConfig import createAlgorithm

alg = createAlgorithm ( 'MyxAODAnalysis', 'AnalysisAlg' )
alg.OfflineRead = OFlag
alg.TriggerRead = Tflag
alg.MuonRead = Mflag
alg.cutnumber = J_CUTNUMBER
alg.dzcut = DZCUTnumber
alg.drcut = DRCUTnumber
alg.etacut = ETACUTnumber
athAlgSeq += alg

#add in the new trigger algorithm
algtrig = createAlgorithm ( 'MyTrigger', 'AnalysisTrig' )
athAlgSeq += algtrig

# limit the number of events (for testing purposes)
#theApp.EvtMax = 500
theApp.EvtMax = 500000000 # 50000000

#Msg limits
MessageSvc.defaultLimit = 200000  # all messages 

# optional include for reducing printout from athena
include("AthAnalysisBaseComps/SuppressLogging.py")
