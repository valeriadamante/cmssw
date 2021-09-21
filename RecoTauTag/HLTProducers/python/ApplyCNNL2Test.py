import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
from TauMLTools.Production.l2TauNNProducer_cfi import *
from TauMLTools.Production.l2TauTagFilter_cfi import *

def update(process):
    thWp = {
        'do_0':{
            'opt_threshold_3': 0.180858813224404,
            'opt_threshold_4': 0.12267940863785043,
            'opt_threshold_5': 0.08411243185219064,
        }


    }

    rateValue = 3
    rateWP=("opt_threshold_{}").format(rateValue)
    dropoutValue = '0'
    dropoutWP = ("do_{}").format(dropoutValue)
    graphPath = ('TauMLTools/Analysis/config/graph_model/Saved_model_{}Dropout.pb').format(dropoutValue)

    normalizationDict = 'TauMLTools/Analysis/config/NormalizationDict.json'


    process = customizeHLTforPatatrackTriplets(process)
    process.l2TauNNProducer = l2TauNNProducer.clone(
        debugLevel = 5,
        processName = cms.string('l2TauNNProducer'),

        #L1TauTrigger=cms.InputTag('hltL1sDoubleTauBigOR'),
        L1Taus= cms.VPSet(
            cms.PSet(
                L1CollectionName = cms.string('L1BigOR'),
                L1TauTrigger =cms.InputTag('hltL1sDoubleTauBigOR'),
            ),
        ),
        ecalInputs =cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
        hbheInput = cms.InputTag("hltHbhereco"),
        hoInput = cms.InputTag("hltHoreco"),
        pataVertices = cms.InputTag("hltPixelVerticesSoA"),
        old_pataVertices = cms.InputTag("hltTrimmedPixelVertices"),
        pataTracks = cms.InputTag("hltPixelTracksSoA"),
        BeamSpot = cms.InputTag("hltOnlineBeamSpot"),
        graphPath = cms.string(graphPath),
        normalizationDict = cms.string(normalizationDict)
    )
    process.l2TauTagFilter = l2TauTagFilter.clone(
        debugLevel = 1,
        processName = cms.string('l2TauTagFilter'),
        L2outcomes = ('l2TauNNProducer', 'L1BigOR'),
        discr_threshold = cms.double(thWp[dropoutWP][rateWP])
    )
    # L2 updated Sequence
    process.L2Sequence = cms.Sequence(process.HLTDoCaloSequence +  process.l2TauNNProducer + process.l2TauTagFilter)

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
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4 = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR +
    process.hltPreDoubleMediumChargedIsoPFTauHPS35Trk1eta2p1Reg +  process.L2Sequence + process.HLTGlobalPFTauHPSSequence +
    process.HLTHPSDoublePFTauPt35Eta2p1Trk1Reg + process.HLTHPSMediumChargedIsoPFTauSequenceReg +
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationReg +
    process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg +
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg +
    process.HLTEndSequence, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)


    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4,
    process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    return process
