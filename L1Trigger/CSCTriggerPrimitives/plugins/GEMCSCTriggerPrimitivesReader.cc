#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//user include files below
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "Validation/MuonCSCDigis/interface/CSCStubMatcher.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCluster.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "DataFormats/LCTDebug/interface/LCTDebug.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

using namespace edm;

class GEMCSCTriggerPrimitivesReader : public edm::one::EDAnalyzer<> {
public:
  explicit GEMCSCTriggerPrimitivesReader(const edm::ParameterSet&);
    // geomToken_(esConsumes());
  ~GEMCSCTriggerPrimitivesReader(){};
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  TTree* bookTTree_ALCT();
  TTree* bookTTree_CLCT();
  TTree* bookTTree_LCT();
  void SaveALCTs(const CSCALCTDigiCollection* alcts, bool is_data, bool is_emul);
  void SaveCLCTs(const CSCCLCTDigiCollection* clcts, bool is_data, bool is_emul);
  void SaveLCTs(const CSCCorrelatedLCTDigiCollection* lcts, bool is_data, bool is_emul, const std::vector<LCTDebugobject>* lctdebugvector_);
  GlobalPoint getGlobalPosition(CSCDetId& cscId, const CSCCorrelatedLCTDigi& lct);

  int eventsAnalyzed;

  // TTree variables
  edm::Service<TFileService> fs;
  TTree* tree;
  TTree* ALCT_tree;
  TTree* CLCT_tree;
  TTree* LCT_tree;
  TTree* t;

  int t_RUN;
  long unsigned int t_Event;
  int t_eventsAnalyzed;
  int t_luminosityblock;

  std::vector<bool> t_is_data;
  std::vector<bool> t_is_emul;
  std::vector<int> t_endcap;
  std::vector<int> t_station;
  std::vector<int> t_ring;
  std::vector<int> t_chamber;
  std::vector<bool> t_isValid;
  std::vector<int> t_quality;
  std::vector<int> t_keyWG;
  std::vector<int> t_strip;
  std::vector<int> t_quadStrip;
  std::vector<int> t_eightStrip;
  std::vector<int> t_stripType;
  std::vector<int> t_bend;
  std::vector<int> t_slope;
  std::vector<int> t_bx;
  std::vector<int> t_pattern;
  std::vector<int> t_run3pattern;
  std::vector<int> t_fractional_slope;
  std::vector<int> t_quartStripBit;
  std::vector<int> t_getEighthStripBit;
  std::vector<int> t_cfeb;
  std::vector<int> t_trknmb;
  std::vector<int> t_getKeyStrip;
  std::vector<int> t_fractional_strip;
  std::vector<int> t_full_bx;
  std::vector<int> t_comp_code;
  std::vector<int> t_collisionB;
  std::vector<int> t_accel;
  std::vector<float> t_eta;
  std::vector<float> t_phi;

  std::vector<int> bendinganglevectorlayer1;
  std::vector<int> bendinganglenoalignmentcorrectionvectorlayer1;
  std::vector<bool> layer2boolvector;
  std::vector<int> cluster_keystripvectorlayer1;
  std::vector<int> residualvectorlayer1;
  std::vector<int> residualnoalignmentcorrectionvectorlayer1;
  std::vector<int> clusterrollvectorlayer1;
  std::vector<int> clusterbxvectorlayer1;
  std::vector<bool> isme1avector;
  std::vector<bool> layer1matchvector;

  std::vector<int> bendinganglevectorlayer2;
  std::vector<int> bendinganglenoalignmentcorrectionvectorlayer2;
  std::vector<int> cluster_keystripvectorlayer2;
  std::vector<int> residualvectorlayer2;
  std::vector<int> residualnoalignmentcorrectionvectorlayer2;
  std::vector<int> clusterrollvectorlayer2;
  std::vector<int> clusterbxvectorlayer2;
  std::vector<bool> layer2matchvector;

  // Run number, Event number
  int RUN_;
  long unsigned int Event_;
  int luminosityblock;



  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_d_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_d_token_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lcts_tmb_d_token_;

  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_e_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_e_token_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lcts_tmb_e_token_;

  edm::EDGetTokenT<std::vector<LCTDebugobject>> lctdebugtoken_data;
  edm::EDGetTokenT<std::vector<LCTDebugobject>> lctdebugtoken_emul;

  // // Producer's labels
  std::string lctProducerData_;
  std::string lctProducerEmul_;
  edm::ParameterSet config;

  bool debug;
  const CSCGeometry* cscGeometry_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  // edm::ESHandle<CSCGeometry> CSCGeometry_;
  // const edm::ConsumesCollector iC;
};


GEMCSCTriggerPrimitivesReader::GEMCSCTriggerPrimitivesReader(const edm::ParameterSet& iConfig) : eventsAnalyzed(0) {
    lctProducerData_ = iConfig.getUntrackedParameter<std::string>("CSCLCTProducerData", "muonCSCDigis");
    lctProducerEmul_ = iConfig.getUntrackedParameter<std::string>("CSCLCTProducerEmul", "cscTriggerPrimitiveDigis");//simCscTriggerPrimitivesDigis for simulation
    debug = iConfig.getParameter<bool>("debug");
    if(debug) std::cout<<"lctProducerData: "<<lctProducerData_<<std::endl;
    std::cout<<"lctProducerEmul: "<<lctProducerEmul_<<std::endl;
    config = iConfig;

    //consumes Data
    alcts_d_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCALCTDigi"));
    clcts_d_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCLCTDigi"));
    lcts_tmb_d_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCorrelatedLCTDigi"));
    lctdebugtoken_data = consumes<std::vector<LCTDebugobject>>(edm::InputTag(lctProducerEmul_, "LCTDebugvector"));


    //consumes Emul
    alcts_e_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    clcts_e_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    lcts_tmb_e_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    lctdebugtoken_emul = consumes<std::vector<LCTDebugobject>>(edm::InputTag(lctProducerEmul_, "LCTDebugvector"));
    //book TTree
    ALCT_tree = bookTTree_ALCT();
    CLCT_tree = bookTTree_CLCT();
    LCT_tree = bookTTree_LCT();


    cscGeomToken_ = esConsumes<CSCGeometry, MuonGeometryRecord>();

  }


void
GEMCSCTriggerPrimitivesReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    ++eventsAnalyzed;
    cscGeometry_ = &iSetup.getData(cscGeomToken_);
    RUN_ = iEvent.id().run();
    Event_ = iEvent.id().event();
    luminosityblock = iEvent.id().luminosityBlock();
    if(debug)std::cout<<"RUN: "<<RUN_<<" Event: "<<Event_<<std::endl;

    //get Data
    edm::Handle<CSCALCTDigiCollection> alcts_data;
    edm::Handle<CSCCLCTDigiCollection> clcts_data;
    edm::Handle<CSCCorrelatedLCTDigiCollection> lcts_tmb_data;
    edm::Handle<std::vector<LCTDebugobject>> lctdebugvectorhandle_data;

    iEvent.getByToken(alcts_d_token_, alcts_data);
    iEvent.getByToken(clcts_d_token_, clcts_data);
    iEvent.getByToken(lcts_tmb_d_token_, lcts_tmb_data);
    iEvent.getByToken(lctdebugtoken_data, lctdebugvectorhandle_data);


    //get Emul
    edm::Handle<CSCALCTDigiCollection> alcts_emul;
    edm::Handle<CSCCLCTDigiCollection> clcts_emul;
    edm::Handle<CSCCorrelatedLCTDigiCollection> lcts_tmb_emul;
    edm::Handle<std::vector<LCTDebugobject>> lctdebugvectorhandle_emul;
    iEvent.getByToken(alcts_e_token_, alcts_emul);
    iEvent.getByToken(clcts_e_token_, clcts_emul);
    iEvent.getByToken(lcts_tmb_e_token_, lcts_tmb_emul);
    iEvent.getByToken(lctdebugtoken_emul, lctdebugvectorhandle_emul);

        // Save ALCTs
    if (alcts_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveALCTs(alcts_data.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Data ALCTs not found"<<std::endl;

    if (alcts_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveALCTs(alcts_emul.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Emul ALCTs not found"<<std::endl;
    // Save CLCTs
    if (clcts_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveCLCTs(clcts_data.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Data CLCTs not found"<<std::endl;
    if (clcts_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveCLCTs(clcts_emul.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Emul CLCTs not found"<<std::endl;
    // Save LCTs
    if (lcts_tmb_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveLCTs(lcts_tmb_data.product(), is_data, is_emul, lctdebugvectorhandle_data.product());
    }
    else if(debug) std::cout<<"Data LCTs not found"<<std::endl;
    if (lcts_tmb_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveLCTs(lcts_tmb_emul.product(), is_data, is_emul, lctdebugvectorhandle_emul.product());
    }
    else if(debug) std::cout<<"Emul LCTs not found"<<std::endl;

}


GlobalPoint GEMCSCTriggerPrimitivesReader::getGlobalPosition(CSCDetId& cscId, const CSCCorrelatedLCTDigi& lct){
  CSCDetId keyId(cscId.endcap(), cscId.station(), cscId.ring(), cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);
  float fractional_strip = lct.getFractionalStrip();
  // case ME1/1
  if (cscId.station() == 1 and (cscId.ring() == 4 || cscId.ring() == 1)) {
    int ring = 1;  // Default to ME1/b
    if (lct.getStrip() > CSCConstants::MAX_HALF_STRIP_ME1B) {
      ring = 4;  // Change to ME1/a if the HalfStrip Number exceeds the range of ME1/b
      fractional_strip -= CSCConstants::NUM_STRIPS_ME1B;
    }
    CSCDetId cscId_(cscId.endcap(), cscId.station(), ring, cscId.chamber(), cscId.layer());
    cscId = cscId_;
  }
  // regular cases
  const auto& chamber = cscGeometry_->chamber(cscId);
  const auto& layer_geo = chamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();
  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  const LocalPoint& csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  const GlobalPoint& csc_gp = cscGeometry_->idToDet(keyId)->surface().toGlobal(csc_intersect);
  return csc_gp;
}



void GEMCSCTriggerPrimitivesReader::SaveALCTs(const CSCALCTDigiCollection* alcts, bool is_data, bool is_emul)
{
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
  t_luminosityblock = luminosityblock;
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCALCTDigi> alctV;
          const auto& range = alcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              if(debug) std::cout << "ALCT "  <<(*digiIt) << " data:emul "<< is_data <<":"<< is_emul << std::endl;
              alctV.push_back(*digiIt);
              //Fill TTree
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_accel.push_back((*digiIt).getAccelerator());
              t_collisionB.push_back((*digiIt).getCollisionB());
              t_keyWG.push_back((*digiIt).getKeyWG());
              t_bx.push_back((*digiIt).getBX());
              t_trknmb.push_back((*digiIt).getTrknmb());
            }
          }
        }
      }
    }
  }
  ALCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_accel.clear();
  t_collisionB.clear();
  t_keyWG.clear();
  t_bx.clear();
  t_trknmb.clear();

}

void GEMCSCTriggerPrimitivesReader::SaveCLCTs(const CSCCLCTDigiCollection* clcts, bool is_data, bool is_emul)
{
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
  t_luminosityblock = luminosityblock;
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCCLCTDigi> clctV;
          const auto& range = clcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              if (debug) std::cout << "CLCT " << (*digiIt) <<" data:emul "<< is_data <<":"<< is_emul << std::endl;
              clctV.push_back(*digiIt);
              //Fill TTree
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_pattern.push_back((*digiIt).getPattern());
              t_run3pattern.push_back((*digiIt).getRun3Pattern());
              t_slope.push_back((*digiIt).getSlope());
              t_fractional_slope.push_back((*digiIt).getFractionalSlope());
              t_stripType.push_back((*digiIt).getStripType());
              t_bend.push_back((*digiIt).getBend());
              t_quartStripBit.push_back((*digiIt).getQuartStripBit());
              t_getEighthStripBit.push_back((*digiIt).getEighthStripBit());
              t_cfeb.push_back((*digiIt).getCFEB());
              t_trknmb.push_back((*digiIt).getTrknmb());
              t_getKeyStrip.push_back((*digiIt).getKeyStrip());
              t_fractional_strip.push_back((*digiIt).getFractionalStrip());
              t_full_bx.push_back((*digiIt).getFullBX());
              t_comp_code.push_back((*digiIt).getCompCode());
              t_bx.push_back((*digiIt).getBX());
              t_strip.push_back((*digiIt).getStrip());
            }
          }
        }
      }
    }
  }
  CLCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_pattern.clear();
  t_run3pattern.clear();
  t_slope.clear();
  t_fractional_slope.clear();
  t_stripType.clear();
  t_bend.clear();
  t_quartStripBit.clear();
  t_getEighthStripBit.clear();
  t_cfeb.clear();
  t_trknmb.clear();
  t_getKeyStrip.clear();
  t_fractional_strip.clear();
  t_full_bx.clear();
  t_comp_code.clear();
  t_bx.clear();
  t_strip.clear();
}

void GEMCSCTriggerPrimitivesReader::SaveLCTs(const CSCCorrelatedLCTDigiCollection* lcts, bool is_data, bool is_emul, const std::vector<LCTDebugobject>* lctdebugvector_)
{
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
  t_luminosityblock = luminosityblock;
  // edm::ParameterSet srLUTset;
  // srLUTset.addUntrackedParameter<bool>("ReadLUTs",true);
  // for (unsigned i=0; i<testvector_->size(); i++){
  //   int value = testvector_->at(i);
  //   testvector.push_back(value);
  // }
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCCorrelatedLCTDigi> lctV;
          const auto& range = lcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              GlobalPoint globalposition = getGlobalPosition(detid, *digiIt);
              t_eta.push_back(globalposition.eta());
              t_phi.push_back(globalposition.phi());
              // std::cout<<test.eta()<<std::endl;
              // std::cout<<"phi"<<test.phi()<<std::endl;
              // std::cout << "LCT " << (*digiIt) <<" data:emul "<< is_data <<":"<< is_emul << std::endl;//if (debug)
              lctV.push_back(*digiIt);
              //Fill TTree
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_keyWG.push_back((*digiIt).getKeyWG());
              t_strip.push_back((*digiIt).getStrip());
              t_quadStrip.push_back((*digiIt).getStrip(4));
              t_eightStrip.push_back((*digiIt).getStrip(8));
              t_stripType.push_back((*digiIt).getStripType());
              t_bend.push_back((*digiIt).getBend());
              t_slope.push_back((*digiIt).getSlope());
              t_bx.push_back((*digiIt).getBX());
              t_pattern.push_back((*digiIt).getPattern());
              t_run3pattern.push_back((*digiIt).getRun3Pattern());
              int KeyWG, bx, Bend, KeyStrip, slope;
              KeyWG = (*digiIt).getKeyWG();
              bx = (*digiIt).getBX();
              Bend = (*digiIt).getBend();
              KeyStrip = (*digiIt).getStrip();
              slope = (*digiIt).getSlope();
              bool isme1a = false;
              if ((stat==1) && ((*digiIt).getStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)){
                isme1a = true;
              }
              isme1avector.push_back(isme1a);
              bool debugfound = false;
              for (unsigned debugnum=0; debugnum<lctdebugvector_->size(); debugnum++){

                LCTDebugobject debug = lctdebugvector_->at(debugnum);
                std::vector<int> identifiers = debug.Getidentifiers();
                // std::cout<<"lctdebug identifiers: "<<identifiers[0]<<" "<<identifiers[1]<<" "<<identifiers[2]<<" "<<identifiers[3]<<" "<<identifiers[4]<<std::endl;
                if ((KeyWG==identifiers[0]) && (bx==identifiers[1]) && (Bend==identifiers[2]) && (KeyStrip==identifiers[3]) && (slope==identifiers[4]) && (!debugfound) && (is_emul)){
                  bool layer1match = debug.GetLayer1Match();
                  bool layer2match = debug.GetLayer2Match();
                  layer2boolvector.push_back(debug.Getlayer2bool());
                  layer1matchvector.push_back(layer1match);
                  layer2matchvector.push_back(layer2match);
                  if (layer1match){
                    bendinganglevectorlayer1.push_back(debug.GetbendingangleLayer1());
                    bendinganglenoalignmentcorrectionvectorlayer1.push_back(debug.GetbendinganglenoalignmentcorrectionLayer1());
                    int cluster_keystriplayer1 = debug.GetGEMClusterKeyStripLayer1();
                    cluster_keystripvectorlayer1.push_back(cluster_keystriplayer1);
                    int residualwithalignmentlayer1 = debug.GetresidualLayer1();
                    int residualwithoutalignmentlayer1 = debug.GetresidualnoalignmentcorrectionLayer1();
                    residualvectorlayer1.push_back(residualwithalignmentlayer1);
                    residualnoalignmentcorrectionvectorlayer1.push_back(residualwithoutalignmentlayer1);
                    int clusterbxlayer1 = debug.GetClusterBxLayer1();

                    clusterbxvectorlayer1.push_back(clusterbxlayer1);
                  } else {
                    bendinganglevectorlayer1.push_back(-999);
                    bendinganglenoalignmentcorrectionvectorlayer1.push_back(-999);
                    cluster_keystripvectorlayer1.push_back(-999);
                    residualvectorlayer1.push_back(-999);
                    residualnoalignmentcorrectionvectorlayer1.push_back(-999);
                    clusterbxvectorlayer1.push_back(-999);
                  }
                  if (layer2match){
                    bendinganglevectorlayer2.push_back(debug.GetbendingangleLayer2());
                    bendinganglenoalignmentcorrectionvectorlayer2.push_back(debug.GetbendinganglenoalignmentcorrectionLayer2());
                    int cluster_keystriplayer2 = debug.GetGEMClusterKeyStripLayer2();
                    cluster_keystripvectorlayer2.push_back(cluster_keystriplayer2);
                    int residualwithalignmentlayer2 = debug.GetresidualLayer2();
                    int residualwithoutalignmentlayer2 = debug.GetresidualnoalignmentcorrectionLayer2();
                    residualvectorlayer2.push_back(residualwithalignmentlayer2);
                    residualnoalignmentcorrectionvectorlayer2.push_back(residualwithoutalignmentlayer2);
                    int clusterbxlayer2 = debug.GetClusterBxLayer2();
                    clusterbxvectorlayer2.push_back(clusterbxlayer2);
                  } else{
                    bendinganglevectorlayer2.push_back(-999);
                    bendinganglenoalignmentcorrectionvectorlayer2.push_back(-999);
                    cluster_keystripvectorlayer2.push_back(-999);
                    residualvectorlayer2.push_back(-999);
                    residualnoalignmentcorrectionvectorlayer2.push_back(-999);
                    clusterbxvectorlayer2.push_back(-999);
                  }

                  int clusterrolllayer1 = debug.GetClusterRoll1();
                  int clusterrolllayer2 = debug.GetClusterRoll2();
                  clusterrollvectorlayer1.push_back(clusterrolllayer1);
                  clusterrollvectorlayer2.push_back(clusterrolllayer2);

                  debugfound = true;
                }

              }
              if (!debugfound){
                bendinganglevectorlayer1.push_back(-999);
                bendinganglenoalignmentcorrectionvectorlayer1.push_back(-999);
                layer2boolvector.push_back(false);
                layer1matchvector.push_back(false);
                cluster_keystripvectorlayer1.push_back(-999);
                residualvectorlayer1.push_back(-999);
                residualnoalignmentcorrectionvectorlayer1.push_back(-999);
                clusterrollvectorlayer1.push_back(-999);
                clusterbxvectorlayer1.push_back(-999);

                bendinganglevectorlayer2.push_back(-999);
                bendinganglenoalignmentcorrectionvectorlayer2.push_back(-999);
                layer2matchvector.push_back(false);
                cluster_keystripvectorlayer2.push_back(-999);
                residualvectorlayer2.push_back(-999);
                residualnoalignmentcorrectionvectorlayer2.push_back(-999);
                clusterrollvectorlayer2.push_back(-999);
                clusterbxvectorlayer2.push_back(-999);
              }
            }
          }
        }
      }
    }
  }
  LCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_keyWG.clear();
  t_strip.clear();
  t_quadStrip.clear();
  t_eightStrip.clear();
  t_stripType.clear();
  t_bend.clear();
  t_slope.clear();
  t_bx.clear();
  t_pattern.clear();
  t_run3pattern.clear();
  t_eta.clear();
  t_phi.clear();
  bendinganglenoalignmentcorrectionvectorlayer1.clear();
  bendinganglevectorlayer1.clear();
  layer2boolvector.clear();
  cluster_keystripvectorlayer1.clear();
  residualvectorlayer1.clear();
  residualnoalignmentcorrectionvectorlayer1.clear();
  clusterrollvectorlayer1.clear();
  isme1avector.clear();
  clusterbxvectorlayer1.clear();
  layer1matchvector.clear();

  bendinganglenoalignmentcorrectionvectorlayer2.clear();
  bendinganglevectorlayer2.clear();
  cluster_keystripvectorlayer2.clear();
  residualvectorlayer2.clear();
  residualnoalignmentcorrectionvectorlayer2.clear();
  clusterrollvectorlayer2.clear();
  clusterbxvectorlayer2.clear();
  layer2matchvector.clear();
}

TTree* GEMCSCTriggerPrimitivesReader::bookTTree_ALCT() {
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("ALCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/g");
  t->Branch("LuminosityBlock", &t_luminosityblock, "LuminosityBlock/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("keyWG", "std::vector<int>", &t_keyWG);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("trknmb", "std::vector<int>", &t_trknmb);
  t->Branch("accel", "std::vector<int>", &t_accel);
  t->Branch("collisionB", "std::vector<int>", &t_collisionB);
  return t;
}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_CLCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("CLCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/g");
  t->Branch("LuminosityBlock", &t_luminosityblock, "LuminosityBlock/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("pattern", "std::vector<int>", &t_pattern);
  t->Branch("run3pattern", "std::vector<int>", &t_run3pattern);
  t->Branch("slope", "std::vector<int>", &t_slope);
  t->Branch("fractional_slope", "std::vector<int>", &t_fractional_slope);
  t->Branch("stripType", "std::vector<int>", &t_stripType);
  t->Branch("bend", "std::vector<int>", &t_bend);
  t->Branch("quartStripBit", "std::vector<int>", &t_quartStripBit);
  t->Branch("getEighthStripBit", "std::vector<int>", &t_getEighthStripBit);
  t->Branch("cfeb", "std::vector<int>", &t_cfeb);
  t->Branch("trknmb", "std::vector<int>", &t_trknmb);
  t->Branch("getKeyStrip", "std::vector<int>", &t_getKeyStrip);
  t->Branch("fractional_strip", "std::vector<int>", &t_fractional_strip);
  t->Branch("full_bx", "std::vector<int>", &t_full_bx);
  t->Branch("comp_code", "std::vector<int>", &t_comp_code);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("strip", "std::vector<int>", &t_strip);

  return t;

}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_LCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("LCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/g");
  t->Branch("LuminosityBlock", &t_luminosityblock, "LuminosityBlock/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("keyWG", "std::vector<int>", &t_keyWG);
  t->Branch("strip", "std::vector<int>", &t_strip);
  t->Branch("quadStrip", "std::vector<int>", &t_quadStrip);
  t->Branch("eightStrip", "std::vector<int>", &t_eightStrip);
  t->Branch("stripType", "std::vector<int>", &t_stripType);
  t->Branch("bend", "std::vector<int>", &t_bend);
  t->Branch("slope", "std::vector<int>", &t_slope);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("pattern", "std::vector<int>", &t_pattern);
  t->Branch("run3pattern", "std::vector<int>", &t_run3pattern);
  t->Branch("bendinganglelayer1", "std::vector<int>", &bendinganglevectorlayer1);
  t->Branch("bendinganglenoalignmentcorrectionlayer1", "std::vector<int>", &bendinganglenoalignmentcorrectionvectorlayer1);
  t->Branch("layer2bool", "std::vector<bool>", &layer2boolvector);
  t->Branch("cluster_keystriplayer1","std::vector<int>", &cluster_keystripvectorlayer1);
  t->Branch("residuallayer1","std::vector<int>",&residualvectorlayer1);
  t->Branch("residualnoalignmentcorrectionlayer1","std::vector<int>",&residualnoalignmentcorrectionvectorlayer1);
  t->Branch("clusterrolllayer1","std::vector<int>",&clusterrollvectorlayer1);
  t->Branch("isme1a","std::vector<bool>",&isme1avector);
  t->Branch("clusterbxlayer1","std::vector<int>",&clusterbxvectorlayer1);
  t->Branch("layer1match","std::vector<bool>",&layer1matchvector);

  t->Branch("bendinganglelayer2", "std::vector<int>", &bendinganglevectorlayer2);
  t->Branch("bendinganglenoalignmentcorrectionlayer2", "std::vector<int>", &bendinganglenoalignmentcorrectionvectorlayer2);
  t->Branch("cluster_keystriplayer2","std::vector<int>", &cluster_keystripvectorlayer2);
  t->Branch("residuallayer2","std::vector<int>",&residualvectorlayer2);
  t->Branch("residualnoalignmentcorrectionlayer2","std::vector<int>",&residualnoalignmentcorrectionvectorlayer2);
  t->Branch("clusterrolllayer2","std::vector<int>",&clusterrollvectorlayer2);
  t->Branch("clusterbxlayer2","std::vector<int>",&clusterbxvectorlayer2);
  t->Branch("layer2match","std::vector<bool>",&layer2matchvector);
  t->Branch("eta","std::vector<float>",&t_eta);
  t->Branch("phi","std::vector<float>",&t_phi);
  return t;
}

DEFINE_FWK_MODULE(GEMCSCTriggerPrimitivesReader);