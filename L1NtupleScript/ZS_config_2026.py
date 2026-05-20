# L1Ntuple for ZDC studies with 2023 data
# Hannah Bossi, <hannah.bossi@cern.ch>
# 9/17/2023
import FWCore.ParameterSet.Config as cms
import os
import sys

#from Configuration.Eras.Era_Run3_pp_on_PbPb_cff import Run3_pp_on_PbPb
from Configuration.Eras.Era_Run3_pp_on_PbPb_2025_cff import Run3_pp_on_PbPb_2025
process = cms.Process('RAW2DIGI', Run3_pp_on_PbPb_2025)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')  #
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi') #
process.load('FWCore.MessageService.MessageLogger_cfi') #
process.load('Configuration.EventContent.EventContent_cff') #
process.load('Configuration.StandardSequences.GeometryRecoDB_cff') #
process.load('Configuration.StandardSequences.MagneticField_cff') # 
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff') #
process.load('Configuration.StandardSequences.EndOfProcess_cff') # 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') # 


from Configuration.AlCa.GlobalTag import GlobalTag
# use the prompt queue for now, will need to update when the full GT becomes available
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_Prompt_Queue', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v4', '')

# needed to supress error from cmssw 14
# process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))

# To change the number of events, change this part
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(600),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 2026 pp run for ZS
        '/store/data/Run2026C/ZeroBias/RAW/v1/000/403/392/00000/792e4bd1-1c12-460e-9bbe-5e32613908b5.root'
        #below 2024 special run
        # '/store/hidata/HIRun2024B/HIZeroBias0/RAW/v1/000/388/702/00000/b1aa4bcf-753b-40ab-a896-b30781b28d0f.root'
        # '/store/hidata/HIRun2024A/HIZeroBias0/RAW/v1/000/387/892/00000/e690293d-75ab-4db2-8650-6ba318af0512.root'
        # '/store/hidata/OORun2025/IonPhysics0/RAW/v1/000/394/175/00000/a3ae4f1d-9c53-482c-ab2c-df24a82f1bc1.root'
        # '/store/hidata/OORun2025/IonPhysics0/RAW/v1/000/394/175/00000/93400011-286f-48db-b79d-a06d241f30c3.root'
        #pp below
        # '/store/data/Run2025F/ZeroBias/RAW/v1/000/397/527/00000/024142b8-3c1f-43a9-9cad-0eb1f7f5bfdf.root'
    )
)

# Output definition
from Configuration.Applications.ConfigBuilder import MassReplaceInputTag

# Additional output definition

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi) #
process.endjob_step = cms.EndPath(process.endOfProcess) #

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step, process.endjob_step) #


from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask #
associatePatAlgosToolsTask(process) #


# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul


# =========Change below for TP change
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 
# from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAWsimEcalTP

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)
# process = L1TReEmulFromRAWsimEcalTP(process)
# ========= 


# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMU #

#call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process) #

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
# from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2025_ZS_tests
# process = L1TSettingsToCaloParams_2025_ZS_tests(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
# from L1Trigger.Configuration.customiseSettings import Tight_ZS_iEta_26_27_28_masked #
# process = Tight_ZS_iEta_26_27_28_masked(process) #


from L1Trigger.Configuration.customiseSettings import Nominal_ZS #
process = Nominal_ZS(process) #


# Automatic addition of the customisation function from L1Trigger.Configuration.customiseUtils
from L1Trigger.Configuration.customiseUtils import L1TGlobalMenuXML  #

#call to customisation function L1TGlobalMenuXML imported from L1Trigger.Configuration.customiseUtils
process = L1TGlobalMenuXML(process) #

# End of customisation functions


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete #
process = customiseEarlyDelete(process) #
# End adding early deletion


# #########################
# # ZDC RecHit Producer && Analyzer
# #########################
# # to prevent crash related to HcalSeverityLevelComputerRcd record
# process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi") #
# process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersPbPb_cff') #
# process.zdcreco_step = cms.Path(process.zdcrecoRun3) #
# process.zdcanalyzer_step = cms.Path(process.zdcanalyzer) #
# process.schedule.append(process.zdcreco_step) #
# process.schedule.append(process.zdcanalyzer_step) #
# #=======================================================================


# # ======================================================================
# # ======================== add in the emulator =========================
# # unpacked etsums
# process.l1UpgradeTree.sumZDCTag = cms.untracked.InputTag("gtStage2Digis","EtSumZDC") #
# process.l1UpgradeTree.sumZDCToken = cms.untracked.InputTag("gtStage2Digis","EtSumZDC") #

# # l1 emulator sums
# process.l1UpgradeEmuTree.sumZDCTag = cms.untracked.InputTag("etSumZdcProducer") #
# process.l1UpgradeEmuTree.sumZDCToken = cms.untracked.InputTag("etSumZdcProducer") #

# # do not change these settings
# process.etSumZdcProducer = cms.EDProducer('L1TZDCProducer', #
#                                           hcalTPDigis = cms.InputTag("simHcalTriggerPrimitiveDigis"), #
# #                                          hcalTPDigis = cms.InputTag("hcalDigis"), 
#                                           bxFirst = cms.int32(-2), #
#                                           bxLast = cms.int32(3) #
# )

# #Via Hannah, for the simHcal collection
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis:ZDC") #
# process.simHcalTriggerPrimitiveDigis.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis:ZDC") #


# process.etSumZdc = cms.Path(process.etSumZdcProducer) #
# process.schedule.append(process.etSumZdc) #
# #======================================================================


#UNCOMMENT HERE TO WORK WITH THE LATEST GREATEST
# MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")
