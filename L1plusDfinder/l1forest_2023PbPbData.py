### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
import glob as glob
process = cms.Process('HiForest',Run3_2023)

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('/store/hidata/HIRun2023A/HIForward0/MINIAOD/PromptReco-v2/000/375/055/00000/4eacc5c8-5e8a-4fc1-80ee-0237cd27e364.root'),
    secondaryFileNames = cms.untracked.vstring(
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/040646b7-0187-47b3-b31b-7d9adb974224.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/051566d8-ea25-4f29-aa02-053ab3b974a1.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/0e30af1b-c718-41e9-a721-83e8f7b346ef.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/2b0f75f3-e0f1-4692-905f-484ac0f4f96a.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/2e8373e5-4087-4c51-ab4f-e09a5b2bcb52.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/363851ec-8263-4819-8d02-821533c0e494.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/3716b848-004f-4567-a205-984955d695b1.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/4707a4a6-2dbe-4228-b96c-8240f519ff9d.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/581e1c6b-371c-4e55-8f2a-e6d64d3128ac.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/592f6602-4858-4f3d-ac2c-41dc1bf82651.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/594f2666-a048-44c1-a8a9-307a902ab534.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/5ac205dd-4bdb-4839-8932-0dd4207bba2d.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/5e843e83-1f81-4148-bd99-8428080900b6.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/655f405a-1543-4330-8b57-035c54671078.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/65e29a49-759d-4087-b02a-090fdfcb574b.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/679e3265-8b9b-46a9-8bbd-5a50493755e5.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/6c1a1813-6f0d-4e40-a452-30f5a9225891.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/79c06a0c-8a81-4f5c-a382-414add9d8206.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/86a33962-ab65-4284-b48c-09629bceae1f.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/8aed86e5-0cd7-4a2a-9574-03e28c147f31.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/8f39f3b6-47c9-432b-adf4-ed552aa57b1a.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/92ca43a5-87af-41ae-a0c7-6c632ca36434.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/959fe1c0-bd71-4434-bdbd-52f3e72226fa.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/ac23bde5-31ab-40d4-bcb1-0ea324fb9218.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/ba3e318d-090c-4bd5-89d6-2183fb85be95.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/c4429ccf-360d-44e9-86a2-564ee4aec539.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/d212fee4-9b0e-48fd-b975-df254d76186f.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/d377433f-eb1e-4532-aab2-1e9b501ccdb3.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/d6311f69-9842-4ce7-af17-5463fa5f2ebf.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/dad79f54-8d45-41b0-88dd-d36fa97e2796.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/dc38a722-290d-4bc5-a667-e707c3871f07.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/e65dbe5f-7cc3-420b-9937-6c10e11f80af.root',
        '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/055/00000/f69d6c20-0575-4527-89bb-10d4362a2156.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_gtStage2Digis_*_RECO',
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

#####################################################################################
# Load Global Tag, geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

#####################################################################################
# Define tree output
#####################################################################################

# process.TFileService = cms.Service("TFileService",
#     fileName = cms.string("HiForestMiniAOD.root"))

############################
# Event analysis
############################

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doHFfilters = cms.bool(False)

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
process.hltobject.triggerNames = trigger_list_data_2023_skimmed

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.particleFlowAnalyser.ptMin = cms.double(0.0)

#########################
# Photons, electrons, and muons
#########################

process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.ggHiNtuplizer.useValMapIso = cms.bool(False)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")

#########################
# Jets
#########################

#process.load("HeavyIonsAnalysis.JetAnalysis.ak2PFJetSequence_ppref_data_cff")
#process.load("HeavyIonsAnalysis.JetAnalysis.ak3PFJetSequence_ppref_data_cff")
process.load("HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff")
#process.load('HeavyIonsAnalysis.JetAnalysis.ak4CaloJetSequence_pp_data_cff')

# The following series of analyzers is to hack in a calorimeter jet correction
process.hltAK4CaloRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
                                                     algorithm = cms.string('AK4Calo'),
                                                     level = cms.string('L2Relative')
)
process.hltAK4CaloAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
                                                     algorithm = cms.string('AK4Calo'),
                                                     level = cms.string('L3Absolute')
                                                 )
process.hltAK4CaloCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
                                             correctors = cms.VInputTag("hltAK4CaloRelativeCorrector", "hltAK4CaloAbsoluteCorrector")
)
process.hltAK4CaloJetsCorrected = cms.EDProducer("CorrectedCaloJetProducer",
                                                 correctors = cms.VInputTag("hltAK4CaloCorrector"),
                                                 src = cms.InputTag("slimmedCaloJets")
)

#process.ak4CaloJetAnalyzer.jetTag = cms.InputTag("hltAK4CaloJetsCorrected")
# End calorimeter jet correction hack

#########################
# Tracks
#########################

process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")

#########################
# ZDC
#########################

process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcdigi.SOI = cms.untracked.int32(2)
process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False
process.zdcanalyzer.nZdcTs = cms.int32(6)

from CondCore.CondDB.CondDB_cfi import *
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("HcalElectronicsMapRcd"),
            tag = cms.string("HcalElectronicsMap_2021_v2.0_data")
        )
    ),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        authenticationMethod = cms.untracked.uint32(1)
)

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource(
    'HcalTextCalibrations',
    input = cms.VPSet(
        cms.PSet(
            object = cms.string('ElectronicsMap'),
            file = cms.FileInPath("L1Trigger/L1TZDC/data/emap_2023_newZDC_v3.txt")
        )
    )
)

#########################
# Main analysis list
#########################

process.forest = cms.Path(
    process.HiForestInfo +
    process.hiEvtAnalyzer +
    process.hltanalysis +
    # process.hltobject +
    process.l1object +
    process.trackSequencePP +
#    process.hltAK4CaloRelativeCorrector + 
#    process.hltAK4CaloAbsoluteCorrector +
#    process.hltAK4CaloCorrector +
#    process.hltAK4CaloJetsCorrected +
#    process.ak4CaloJetAnalyzer +
#    process.akPu4CaloJetAnalyzer +
    process.particleFlowAnalyser +
    process.ggHiNtuplizer +
    # process.zdcdigi +
    # process.QWzdcreco +
    process.zdcanalyzer +
    process.muonSequencePP
    # process.unpackedMuons +
    # process.muonAnalyzer
)

#########################
# Event selection
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)

process.pAna = cms.EndPath(process.skimanalysis)

#########################
# Jet customization
#########################

# Select the types of jets filled
addR2Jets = False
addR3Jets = False
addR3FlowJets = False
addR4Jets = True
addR4FlowJets = False
addUnsubtractedR4Jets = False

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = True        # Add jet phi and eta for WTA axis

# This is only for non-reclustered jets
addCandidateTagging = False

if addR4Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupPprefJets

    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections                                                                                    
        print("ADD R4 JETS")
        process.jetsR4 = cms.Sequence()
        setupPprefJets('ak04PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
        process.ak04PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.ak04PFpatJetCorrFactors.primaryVertices = "offlineSlimmedPrimaryVertices"
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.ak4PFJetAnalyzer.jetTag = 'ak04PFpatJets'
        process.ak4PFJetAnalyzer.jetName = 'ak04PF'
        process.ak4PFJetAnalyzer.doSubEvent = False # Need to disable this, since there is some issue with the gen jet constituents. More debugging needed is want to use constituents. 
        process.forest += process.extraJetsData * process.jetsR4 * process.ak4PFJetAnalyzer

#########################
# D finder
#########################

AddCaloMuon = False
runOnMC = False
HIFormat = False
UseGenPlusSim = False
VtxLabel = "offlineSlimmedPrimaryVertices"
TrkLabel = "packedPFCandidates"
GenLabel = "prunedGenParticles"
TrkChi2Label = "packedPFCandidateTrackChi2"
useL1Stage2 = True
HLTProName = "HLT"

from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X

finderMaker_75X(process, AddCaloMuon, runOnMC, HIFormat, UseGenPlusSim, VtxLabel, TrkLabel, TrkChi2Label, GenLabel, useL1Stage2, HLTProName)

process.Dfinder.MVAMapLabel = cms.InputTag(TrkLabel, "MVAValues")
process.Dfinder.makeDntuple = cms.bool(True)
process.Dfinder.tkPtCut = cms.double(1.0) # before fit
process.Dfinder.tkEtaCut = cms.double(2.4) # before fit
process.Dfinder.dPtCut = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0) # before fit
process.Dfinder.VtxChiProbCut = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.0, 0.0, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.05, 0.05)
process.Dfinder.dCutSeparating_PtVal = cms.vdouble(5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5.)
process.Dfinder.tktkRes_svpvDistanceCut_lowptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.tktkRes_svpvDistanceCut_highptD = cms.vdouble(0., 0., 0., 0., 0., 0., 0., 0., 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0.)
process.Dfinder.svpvDistanceCut_lowptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.svpvDistanceCut_highptD = cms.vdouble(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 0., 0., 0., 0., 0., 0., 2.5, 2.5)
process.Dfinder.Dchannel = cms.vint32(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
process.Dfinder.dropUnusedTracks = cms.bool(True)
process.Dfinder.detailMode = cms.bool(False)
process.Dfinder.printInfo = cms.bool(False)

process.dfinder = cms.Path(process.DfinderSequence)

#####################################################################################
# L1 emulation
#####################################################################################

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.forest,process.dfinder,process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW

# Call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMU

# Call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParamsHI_2024_v0_0_0

# Call to customisation function L1TSettingsToCaloParamsHI_2024_v0_4_2 imported from L1Trigger.Configuration.customiseSettings
process = L1TSettingsToCaloParamsHI_2024_v0_0_0(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseUtils
from L1Trigger.Configuration.customiseUtils import L1TGlobalMenuXML

# Call to customisation function L1TGlobalMenuXML imported from L1Trigger.Configuration.customiseUtils
process = L1TGlobalMenuXML(process)



# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

process.l1UpgradeTree.sumZDCTag = cms.untracked.InputTag("etSumZdcProducer")

#process.HFAdcana = cms.EDAnalyzer("HFAdcToGeV",
#    digiLabel = cms.untracked.InputTag("hcalDigis"),
#    minimized = cms.untracked.bool(True),
#    fillhf = cms.bool(False)
#)

#process.HFAdc = cms.Path(process.HFAdcana)
#process.schedule.append(process.HFAdc)

process.hcalDigis.saveQIE10DataNSamples = cms.untracked.vint32(6) 
process.hcalDigis.saveQIE10DataTags = cms.untracked.vstring( "MYDATA" )

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel

process.hltfilter = hltHighLevel.clone(
   HLTPaths = [
        # Double muons
        'HLT_HIUPC_DoubleMuCosmic*_MaxPixelCluster1000_v*',
        'HLT_HIUPC_DoubleMuOpen*_NotMBHF*AND_v*',

        # Not MB
        'HLT_HIUPC_NotMBHF*_v*',

        # Jet triggers
        'HLT_HIUPC_SingleJet*_ZDC1n*XOR_*MaxPixelCluster*',
        'HLT_HIUPC_SingleJet*_NotMBHF2AND_*MaxPixelCluster*',

        # Single muon
        'HLT_HIUPC_SingleMu*_NotMBHF*_MaxPixelCluster*',

        # ZDC 1n or, low pixel clusters
        'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
        'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*',

        # ZB, single pixel track
        'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*'
   ]
)
#process.hltfilter.andOr = cms.bool(True)  # True = OR, False = AND between the HLT paths
#process.hltfilter.throw = cms.bool(False) # throw exception on unknown path names

process.filterSequence = cms.Sequence(
    process.hltfilter * 
    process.primaryVertexFilter
)

process.superFilterPath = cms.Path(process.filterSequence)
process.skimanalysis.superFilters = cms.vstring("superFilterPath")

for path in process.paths:
   getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")
