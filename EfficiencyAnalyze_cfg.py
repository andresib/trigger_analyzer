import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      # Multiple file should be comma separated
      # This is the format for using a remote file
      #'root://cmseos.fnal.gov//store/user/cmsdas/2022/pre_exercises/Set4/Input/DoubleMuon/slimMiniAOD_data_MuEle_1.root',
      'root://eoscms.cern.ch//eos/cms/store/data/Run2018D/MET/MINIAOD/PromptReco-v2/000/321/122/00000/786F4462-B19E-E811-8D2E-FA163EE048DE.root',
      # The format for using a local file can be found in the commented line below
      # 'file:slimMiniAOD_data_MuEle_1.root'
  )
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.analyzeBasicPat = cms.EDAnalyzer("EfficiencyAnalyzer",
  muonSrc = cms.untracked.InputTag("slimmedMuons"),                  
  elecSrc = cms.untracked.InputTag("slimmedElectrons"),
  tauSrc = cms.untracked.InputTag("slimmedTaus"),
  jetSrc = cms.untracked.InputTag("slimmedJets"),
  metSrc = cms.untracked.InputTag("slimmedMETs"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('fullevents_cuts2denom.root')
                                   )

process.p = cms.Path(process.analyzeBasicPat)

