
#See: https://twiki.cern.ch/twiki/bin/viewauth/AtlasComputing/SoftwareTutorialxAODAnalysisInCMake for more details about anything here

#testFile = os.getenv("ALRB_TutorialData") + '/r9315/mc16_13TeV.410501.PowhegPythia8EvtGen_A14_ttbar_hdamp258p75_nonallhad.merge.AOD.e5458_s3126_r9364_r9315/AOD.11182705._000001.pool.root.1'
testFile = "/eos/atlas/atlasdatadisk/rucio/mc16_13TeV/05/e6/AOD.19016897._000001.pool.root.1"

#override next line on command line with: --filesInput=XXX
jps.AthenaCommonFlags.FilesInput = [testFile] 

#Specify AccessMode (read mode) ... ClassAccess is good default for xAOD
jps.AthenaCommonFlags.AccessMode = "ClassAccess" 

jps.AthenaCommonFlags.HistOutputs = ["ANALYSIS:MyxAODAnalysis.outputs.root"]
svcMgr.THistSvc.MaxFileSize=-1 #speeds up jobs that output lots of histograms


# Create the algorithm's configuration.
from AnaAlgorithm.DualUseConfig import createAlgorithm
alg = createAlgorithm ( 'MyxAODAnalysis', 'AnalysisAlg' )

# later on we'll add some configuration options for our algorithm that go here
#alg.ElectronPtCut = 30000.0
#alg.SampleName = 'Zee'
# Add our algorithm to the main alg sequence
athAlgSeq += alg

# limit the number of events (for testing purposes)
#theApp.EvtMax = 500
theApp.EvtMax = 1000

#Msg limits
MessageSvc.defaultLimit = 1000  # all messages 

# optional include for reducing printout from athena
include("AthAnalysisBaseComps/SuppressLogging.py")