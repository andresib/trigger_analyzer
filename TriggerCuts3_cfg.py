import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      # Multiple file should be comma separated
      # This is the format for using a remote file
      #'root://cmseos.fnal.gov//store/user/cmsdas/2022/pre_exercises/Set4/Input/DoubleMuon/slimMiniAOD_data_MuEle_1.root',
     #'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/786F4462-B19E-E811-8D2E-FA163EE048DE.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/183D2890-AA9E-E811-9E64-FA163E1E3B31.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/52D67EA2-7C9F-E811-99CC-02163E01A047.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/6AF02A7D-A19E-E811-B222-FA163E761FF5.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/72792556-36A0-E811-B8A1-02163E019ED6.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/786F4462-B19E-E811-8D2E-FA163EE048DE.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/90802BC5-AC9E-E811-B6C8-FA163E8CE6BD.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/96F72773-9E9E-E811-BF77-FA163EF29411.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/A6B4FF76-A59E-E811-B387-FA163EE653E7.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/A8611447-A89E-E811-8CBD-FA163ED5D1CD.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/F82666A8-A39E-E811-AFA0-FA163E75AD21.root',

      # The format for using a local file can be found in the commented line below
      # 'file:slimMiniAOD_data_MuEle_1.root'
  )
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.analyzeBasicPat = cms.EDAnalyzer("TriggerCuts3",
  triggerR = cms.untracked.InputTag( 'TriggerResults', '', 'HLT' ),
  triggerO = cms.untracked.InputTag("slimmedPatTrigger"),
  triggerNames = cms.untracked.vstring("HLT_DiJet110_35_Mjj650_PFMET110_v", "HLT_DiJet110_35_Mjj650_PFMET120_v", "HLT_DiJet110_35_Mjj650_PFMET130_v"),
 #  triggerNames = cms.untracked.vstring("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v", "HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_v", "HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_v"),
  muonSrc = cms.untracked.InputTag("slimmedMuons"),                  
  elecSrc = cms.untracked.InputTag("slimmedElectrons"),
  tauSrc = cms.untracked.InputTag("slimmedTaus"),
  jetSrc = cms.untracked.InputTag("slimmedJets"),
  metSrc = cms.untracked.InputTag("slimmedMETs"),
)

#process.dump=cms.EDAnalyzer('EventContentAnalyzer' )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('finalcuts3.root')
                                   )

process.p = cms.Path(process.analyzeBasicPat)

#process.Tracer = cms.Service("Tracer")
#process.p = cms.Path(process.analyzeBasicPat*process.dump)
process.p = cms.Path(process.analyzeBasicPat)