import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
from RecoTauTag.HLTProducers.l2TauNNProducer_cfi import *
from RecoTauTag.HLTProducers.l2TauTagFilter_cfi import *

def update(process):
    thWp = {
            'WP_Tight': 0.180858813224404,
            'WP_Medium': 0.12267940863785043,
            'WP_Loose': 0.08411243185219064,
    }

    working_point = "Tight"
    rateWP=("WP_{}").format(working_point)
    graphPath = 'RecoTauTag/TrainingFiles/L2TauNNTag/L2TauTag_Run3v1.pb'

    normalizationDict = 'RecoTauTag/TrainingFiles/L2TauNNTag/NormalizationDict.json'


    process = customizeHLTforPatatrackTriplets(process)
    process.hltL2TauTagNNProducer = l2TauNNProducer.clone(
        debugLevel = 0,
        L1Taus= cms.VPSet(
            cms.PSet(
                L1CollectionName = cms.string('DoubleTau'),
                L1TauTrigger =cms.InputTag('hltL1sDoubleTauBigOR'),
            ),
        ),
        hbheInput = cms.InputTag("hltHbhereco"),
        hoInput = cms.InputTag("hltHoreco"),
        ebInput =cms.InputTag("hltEcalRecHit:EcalRecHitsEB"),
        eeInput =cms.InputTag("hltEcalRecHit:EcalRecHitsEE"),
        pataVertices = cms.InputTag("hltPixelVerticesSoA"),
        pataTracks = cms.InputTag("hltPixelTracksSoA"),
        BeamSpot = cms.InputTag("hltOnlineBeamSpot"),
        graphPath = cms.string(graphPath),
        normalizationDict = cms.string(normalizationDict)
    )
    process.hltL2DoubleTauTagNNFilter = l2TauTagFilter.clone(
        nExpected = 2,
        L1TauSrc = cms.InputTag('hltL1sDoubleTauBigOR'),
        L2Outcomes = ('hltL2TauTagNNProducer', 'DoubleTau'),
        DiscrWP = cms.double(thWp[rateWP])
    ) 
    # L2 updated Sequence
    process.hltL2TauTagNNSequence = cms.Sequence(process.HLTDoCaloSequence + process.hltL1sDoubleTauBigOR + process.hltL2TauTagNNProducer)


    # Regional -> Global customization
    process.hltHpsPFTauTrackPt1DiscriminatorReg.PFTauProducer = cms.InputTag("hltHpsPFTauProducer")
    process.hltHpsDoublePFTau35Reg.inputTag = cms.InputTag( "hltHpsPFTauProducer")
    process.hltHpsSelectedPFTausTrackPt1Reg.src = cms.InputTag( "hltHpsPFTauProducer")
    process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminatorReg.particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminatorReg.particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    process.hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg.src = cms.InputTag( "hltHpsPFTauProducer" )

    # re-define path with l2 updated sequence
    #process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4 = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR +
    #process.hltPreDoubleMediumChargedIsoPFTauHPS35Trk1eta2p1Reg +  process.hltL2TauTagNNSequence + process.hltL2DoubleTauTagNNFilter + process.HLTGlobalPFTauHPSSequence +
    #process.HLTHPSDoublePFTauPt35Eta2p1Trk1Reg + process.HLTHPSMediumChargedIsoPFTauSequenceReg +
    #process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationReg +
    #process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg +
    #process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg +
    #process.HLTEndSequence, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.HLTL2TauJetsL1TauSeededSequence)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltDoubleL2Tau26eta2p2)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.HLTL2p5IsoTauL1TauSeededSequence)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltDoubleL2IsoTau26eta2p2 )
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.HLTRegionalPFTauHPSSequence )

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(3, process.hltL2TauTagNNSequence)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(4, process.hltL2DoubleTauTagNNFilter)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(5, process.HLTGlobalPFTauHPSSequence)


    #process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    return process
