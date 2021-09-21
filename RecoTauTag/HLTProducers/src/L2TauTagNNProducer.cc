/*
Outcomes producedr for L2 hadronic tau selection
*/
// system include files
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <math.h>
#include "Compression.h"
#include "TMath.h"
// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// utilities
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
// Geometry
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
// caloRecHit
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalDetIdCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
// L1 Tau Trigger
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
// Tracks
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "RecoPixelVertexing/PixelTrackFitting/interface/FitUtils.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCurvilinear.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CUDADataFormats/Track/interface/PixelTrackHeterogeneous.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CUDADataFormats/SiPixelCluster/interface/gpuClusteringConstants.h"
//vertices
#include "CUDADataFormats/Vertex/interface/ZVertexSoA.h"
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"


namespace tau_hlt{

  enum class NNInputs {
      nVertices = 0,
      l1Tau_pt,
      l1Tau_eta,
      l1Tau_hwIso,
      EcalEnergySum,
      EcalSize,
      EcalEnergyStdDev,
      EcalDeltaEta,
      EcalDeltaPhi,
      EcalChi2,
      EcalEnergySumForPositiveChi2,
      EcalSizeForPositiveChi2,
      HcalEnergySum,
      HcalSize,
      HcalEnergyStdDev,
      HcalDeltaEta,
      HcalDeltaPhi,
      HcalChi2,
      HcalEnergySumForPositiveChi2,
      HcalSizeForPositiveChi2,
      PatatrackPtSum,
      PatatrackSize,
      PatatrackSizeWithVertex,
      PatatrackPtSumWithVertex,
      PatatrackChargeSum,
      PatatrackDeltaEta,
      PatatrackDeltaPhi,
      PatatrackChi2OverNdof,
      PatatrackNdof,
      PatatrackDxy,
      PatatrackDz
    };

  const std::map<NNInputs, std::string> varNameMap = {
    {NNInputs::nVertices,"nVertices"},
    {NNInputs::l1Tau_pt,"l1Tau_pt"},
    {NNInputs::l1Tau_eta,"l1Tau_eta"},
    {NNInputs::l1Tau_hwIso,"l1Tau_hwIso"},
    {NNInputs::EcalEnergySum,"EcalEnergySum"},
    {NNInputs::EcalSize,"EcalSize"},
    {NNInputs::EcalEnergyStdDev,"EcalEnergyStdDev"},
    {NNInputs::EcalDeltaEta,"EcalDeltaEta"},
    {NNInputs::EcalDeltaPhi,"EcalDeltaPhi"},
    {NNInputs::EcalChi2,"EcalChi2"},
    {NNInputs::EcalEnergySumForPositiveChi2,"EcalEnergySumForPositiveChi2"},
    {NNInputs::EcalSizeForPositiveChi2,"EcalSizeForPositiveChi2"},
    {NNInputs::HcalEnergySum,"HcalEnergySum"},
    {NNInputs::HcalSize,"HcalSize"},
    {NNInputs::HcalEnergyStdDev,"HcalEnergyStdDev"},
    {NNInputs::HcalDeltaEta,"HcalDeltaEta"},
    {NNInputs::HcalDeltaPhi,"HcalDeltaPhi"},
    {NNInputs::HcalChi2,"HcalChi2"},
    {NNInputs::HcalEnergySumForPositiveChi2,"HcalEnergySumForPositiveChi2"},
    {NNInputs::HcalSizeForPositiveChi2,"HcalSizeForPositiveChi2"},
    {NNInputs::PatatrackPtSum,"PatatrackPtSum"},
    {NNInputs::PatatrackSize,"PatatrackSize"},
    {NNInputs::PatatrackSizeWithVertex,"PatatrackSizeWithVertex"},
    {NNInputs::PatatrackPtSumWithVertex,"PatatrackPtSumWithVertex"},
    {NNInputs::PatatrackChargeSum,"PatatrackChargeSum"},
    {NNInputs::PatatrackDeltaEta,"PatatrackDeltaEta"},
    {NNInputs::PatatrackDeltaPhi,"PatatrackDeltaPhi"},
    {NNInputs::PatatrackChi2OverNdof,"PatatrackChi2OverNdof"},
    {NNInputs::PatatrackNdof,"PatatrackNdof"},
    {NNInputs::PatatrackDxy,"PatatrackDxy"},
    {NNInputs::PatatrackDz,"PatatrackDz"}
  };

  struct normDictElement{
    float mean;
    float std;
    float min;
    float max;
  };

  struct L2TauNNProducerCacheData {
    L2TauNNProducerCacheData() : graphDef(nullptr),session(nullptr) {}
    tensorflow::GraphDef* graphDef;
    tensorflow::Session* session;
    std::vector<normDictElement> normVec;
  };

  class L2TauNNProducer : public edm::stream::EDProducer<edm::GlobalCache<L2TauNNProducerCacheData>> {

  public:

    using OutputCollection = std::map<int, std::vector<float>>;

    struct caloRecHitCollections {
      const HBHERecHitCollection  *hbhe;
      const HORecHitCollection *ho;
      const EcalRecHitCollection *eb;
      const EcalRecHitCollection *ee;
      const CaloGeometry *Geometry;
    };

    struct InputDescTau {
      std::string CollectionName;
      edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> input_token;
    };

    static constexpr int nCellEta = 5;
    static constexpr int nCellPhi = 5;
    static constexpr int nVars = 31;
    static constexpr float dR_max = 0.5;
    static constexpr float dR2_max = 0.25;
    static constexpr float dEta_width = 2*dR_max/static_cast<float>(nCellEta);
    static constexpr float dPhi_width = 2*dR_max/static_cast<float>(nCellPhi);

    explicit L2TauNNProducer(const edm::ParameterSet&, const L2TauNNProducerCacheData*);
    ~L2TauNNProducer() override {}
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    static std::unique_ptr<L2TauNNProducerCacheData> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(L2TauNNProducerCacheData*);

  private:
    void checknan(tensorflow::Tensor& tensor, bool printoutTensor, int debugLevel);
    void standardizeTensor(tensorflow::Tensor& tensor);
    std::vector<float> GetTauScore(const tensorflow::Tensor& cellGridMatrix);
    void produce(edm::Event& event, const edm::EventSetup& eventsetup) ;
    void FillL1TauVars(tensorflow::Tensor& cellGridMatrix, const std::vector<l1t::TauRef> allTaus);
    void FillCaloRecHits(tensorflow::Tensor& cellGridMatrix, const std::vector<l1t::TauRef> allTaus, const caloRecHitCollections& caloRecHits);
    void FillPatatracks(tensorflow::Tensor& cellGridMatrix, const std::vector<l1t::TauRef> allTaus, const PixelTrackHeterogeneous& patatracks, const ZVertexHeterogeneous& patavertices,const reco::BeamSpot &beamspot, const MagneticField *magfi);
    std::vector<int> GetVtxGood(const ZVertexHeterogeneous& patavertices, const PixelTrackHeterogeneous& patatracks,const std::vector<int> TrackGood);
    std::vector<float> impactParameter(int it, const PixelTrackHeterogeneous& patatracks, float patatrackPhi,const reco::BeamSpot &beamspot, const MagneticField *magfi );


  private:
    int debugLevel_;
    std::string processName;
    const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> tauTrigger_token;
    std::vector<InputDescTau> L1TauDesc;
    edm::EDGetTokenT<HBHERecHitCollection> hbhe_token;
    edm::EDGetTokenT<HORecHitCollection> ho_token;
    std::vector<edm::InputTag> ecalLabels;
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> Geometry_token;
    edm::EDGetTokenT<ZVertexHeterogeneous> pataVertices_token;
    edm::EDGetTokenT<PixelTrackHeterogeneous> pataTracks_token;
    edm::EDGetTokenT<reco::BeamSpot> BeamSpot_token;
    std::string normalizationDict_;
    std::vector<edm::EDGetTokenT<EcalRecHitCollection>> ecal_tokens;
    std::string inputTensorName_;
    std::string outputTensorName_;
    OutputCollection outputs_;
    const L2TauNNProducerCacheData* L2cacheData;

  };



  std::unique_ptr<L2TauNNProducerCacheData> L2TauNNProducer::initializeGlobalCache(const edm::ParameterSet& cfg) {
    std::unique_ptr<L2TauNNProducerCacheData> cacheData= std::make_unique<L2TauNNProducerCacheData>();

    std::string graphPath = cfg.getParameter<std::string>("graphPath");
    graphPath = edm::FileInPath(graphPath).fullPath();

    cacheData->graphDef = tensorflow::loadGraphDef(graphPath);
    cacheData->session = tensorflow::createSession(cacheData->graphDef);

    tensorflow::setLogging("2");

    boost::property_tree::ptree loadPtreeRoot;
    std::string normalizationDict = cfg.getParameter<std::string>("normalizationDict");
    normalizationDict = edm::FileInPath(normalizationDict).fullPath();
    boost::property_tree::read_json(cfg.getParameter<std::string>("normalizationDict"), loadPtreeRoot);
    for (auto& [key, val] : varNameMap){
      boost::property_tree::ptree var = loadPtreeRoot.get_child(val);
      normDictElement current_element;
      current_element.mean = var.get_child("mean").get_value<float>();
      current_element.std =var.get_child("std").get_value<float>();
      current_element.min =var.get_child("min").get_value<float>();
      current_element.max =var.get_child("max").get_value<float>();
      cacheData->normVec.push_back(current_element);
    }
    return cacheData;
  }
  void L2TauNNProducer::globalEndJob(L2TauNNProducerCacheData* cacheData) {
    if (cacheData->graphDef != nullptr) {
      delete cacheData->graphDef;
    }
    tensorflow::closeSession(cacheData->session);
  }
  void L2TauNNProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<int>("debugLevel", 0)->setComment("set debug level for printing out info");
    desc.add<std::string>("processName","")->setComment("Name of the process");

    edm::ParameterSetDescription l1TausPset;
    l1TausPset.add<std::string>("L1CollectionName", "")->setComment("Name of collections");
    l1TausPset.add<edm::InputTag>("L1TauTrigger", edm::InputTag(""))->setComment("Which trigger should the L1 Taus collection pass");
    desc.addVPSet("L1Taus", l1TausPset);
    //desc.add<edm::InputTag>("L1TauTrigger"); // keep it for debugging
    desc.add<edm::InputTag>("hbheInput", edm::InputTag(""))->setComment("HBHE recHit collection");
    desc.add<edm::InputTag>("hoInput", edm::InputTag(""))->setComment("HO recHit Collection");
    desc.add<std::vector<edm::InputTag>>("ecalInputs")->setComment("EB and EE recHit Collections");
    desc.add<edm::InputTag>("pataVertices", edm::InputTag("hltPixelVerticesSoA"))->setComment("patatrack vertices collection");
    desc.add<edm::InputTag>("pataTracks", edm::InputTag("hltPixelTracksSoA"))->setComment("patatrack collection");
    desc.add<edm::InputTag>("BeamSpot");
    desc.add<std::string>("graphPath","TauMLTools/Analysis/config/graph_model/Saved_model_0Dropout.pb")->setComment("path to the saved CNN");
    desc.add<std::string>("normalizationDict","TauMLTools/Analysis/config/NormalizationDict.json")->setComment("path to the dictionary for variable standardization");
    descriptions.addWithDefaultLabel(desc);
  }


  L2TauNNProducer::L2TauNNProducer(const edm::ParameterSet& cfg, const L2TauNNProducerCacheData* cacheData):
        debugLevel_(cfg.getParameter<int>("debugLevel")),
        processName(cfg.getParameter<std::string>("processName")),
        //tauTrigger_token(consumes<trigger::TriggerFilterObjectWithRefs>(cfg.getParameter<edm::InputTag>("L1TauTrigger"))),  // keep it for debugging
        hbhe_token(consumes<HBHERecHitCollection>(cfg.getParameter<edm::InputTag>("hbheInput"))),
        ho_token(consumes<HORecHitCollection>(cfg.getParameter<edm::InputTag>("hoInput"))),
        ecalLabels(cfg.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
        Geometry_token(esConsumes<CaloGeometry,CaloGeometryRecord>()),
        pataVertices_token(consumes<ZVertexHeterogeneous>(cfg.getParameter<edm::InputTag>("pataVertices"))),
        pataTracks_token(consumes<PixelTrackHeterogeneous>(cfg.getParameter<edm::InputTag>("pataTracks"))),
        BeamSpot_token(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("BeamSpot"))),
        normalizationDict_(cfg.getParameter<std::string>("normalizationDict")),
        inputTensorName_((cacheData->graphDef)->node(0).name()),
        outputTensorName_((cacheData->graphDef)->node((cacheData->graphDef)->node_size()-1).name()),
        L2cacheData(cacheData)
        {
           std::vector<edm::ParameterSet> L1TauCollections_token = cfg.getParameter<std::vector<edm::ParameterSet> >("L1Taus");

          const unsigned nLabels = ecalLabels.size();
          for (unsigned i = 0; i != nLabels; i++)
            ecal_tokens.push_back(consumes<EcalRecHitCollection>(ecalLabels[i]));
          int k = 0 ;
          for (std::vector<edm::ParameterSet>::const_iterator iL1Tau = L1TauCollections_token.begin(); iL1Tau != L1TauCollections_token.end(); ++iL1Tau) {
            InputDescTau toInsert;
            toInsert.CollectionName = iL1Tau->getParameter<std::string>("L1CollectionName");
            toInsert.input_token = consumes<trigger::TriggerFilterObjectWithRefs>(iL1Tau->getParameter<edm::InputTag>("L1TauTrigger"));
            L1TauDesc.push_back(toInsert);
            std::vector<float> outcomes_;
            outputs_.insert(std::pair<int, std::vector<float>>(k, outcomes_));
            k++;

          }
          for (auto& [inp_idx, out_vec] : outputs_ ) {
            produces<std::vector<float>>(L1TauDesc[inp_idx].CollectionName);
          }
          //produces<std::vector<float>>("L1BigOR");  // keep it for debugging

        }

   void L2TauNNProducer::checknan(tensorflow::Tensor& tensor, bool printoutTensor, int debugLevel){
     std::vector<int> tensor_shape;
     for(int d=0; d<tensor.shape().dims(); d++) {
         tensor_shape.push_back(tensor.shape().dim_size(d));
     }
     auto getCell = [&](int tau_idx, int phi_idx, int eta_idx, int var_idx)-> const float& {
       return tensor.tensor<float, 4>()(tau_idx,phi_idx, eta_idx, var_idx);
     };
     for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
       for(int phi_idx =0; phi_idx < tensor_shape.at(1); phi_idx++){
         for(int eta_idx =0; eta_idx < tensor_shape.at(2); eta_idx++){
           for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
             auto nonstd_var = getCell(tau_idx,phi_idx, eta_idx, var_idx);
             if(std::isnan(nonstd_var)){
               std::cout << "var is nan \nvar name= " << varNameMap.at(static_cast<NNInputs>(var_idx))
                << "\t var_idx = " << var_idx << "\t eta_idx = " << eta_idx
                  << "\t phi_idx = " << phi_idx << "\t tau_idx = " << tau_idx << std::endl;
               if(debugLevel>2){
                 std::cout << "other vars in same cell \n";
                 if(var_idx+1 < tensor_shape.at(3)) std::cout << varNameMap.at(static_cast<NNInputs>(var_idx+1))
                   << "\t = " <<getCell(tau_idx,phi_idx, eta_idx, var_idx+1) << std::endl;
                 if(var_idx+2 < tensor_shape.at(3)) std::cout << varNameMap.at(static_cast<NNInputs>(var_idx+2))
                   << "\t = " <<getCell(tau_idx,phi_idx, eta_idx, var_idx+2) << std::endl;
                 if(var_idx+3 < tensor_shape.at(3)) std::cout << varNameMap.at(static_cast<NNInputs>(var_idx+3))
                   << "\t = " <<getCell(tau_idx,phi_idx, eta_idx, var_idx+3) << std::endl;
                 if(var_idx+4 < tensor_shape.at(3)) std::cout << varNameMap.at(static_cast<NNInputs>(var_idx+4))
                   << "\t = " <<getCell(tau_idx,phi_idx, eta_idx, var_idx+4) << std::endl;
               } // end if debugLevel
             } // end if is nan
           } // end loop on var_idx
         } // end loop on eta_idx
       } // end loop on phi_idx
     } // end loop on tau_idx
     if(printoutTensor){
       for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
         for(int phi_idx =0; phi_idx < tensor_shape.at(1); phi_idx++){
           for(int eta_idx =0; eta_idx < tensor_shape.at(2); eta_idx++){
             for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
               // search for a specific tau - needed for debugging
               if(getCell(tau_idx, phi_idx, eta_idx, static_cast<int>(NNInputs::l1Tau_pt)) *256.0 == 165.5){
                 if(debugLevel<5){
                   break;
                 }
                 if(debugLevel>4 && debugLevel<10){
                   std::cout << getCell(tau_idx, phi_idx, eta_idx, var_idx) <<",\n";
                 }
                 else{
                   std::cout << "\nvar name= " << varNameMap.at(static_cast<NNInputs>(var_idx)) <<"\t tau_idx = " << tau_idx<<"\t phi_idx = " << phi_idx << "\t eta_idx = " << eta_idx << " \tvalue =" << getCell(tau_idx,phi_idx, eta_idx, var_idx);
                 }
               } // end if debugLevel
             } // end loop on var_idx
           } // end loop on eta_idx
         } // end loop on phi_idx
       } // end loop on tau_idx
     } // end if printout
   } // end function


   void L2TauNNProducer::standardizeTensor(tensorflow::Tensor& tensor){
     std::vector<int> tensor_shape;
     for(int d=0; d<tensor.shape().dims(); d++) {
         tensor_shape.push_back(tensor.shape().dim_size(d));
     }
     auto getCell = [&](int tau_idx,int phi_idx, int eta_idx, int var_idx)-> float& {
       return tensor.tensor<float, 4>()(tau_idx,phi_idx, eta_idx, var_idx);
     };
     for(int  tau_idx=0; tau_idx < tensor_shape.at(0); tau_idx++){
       for(int  phi_idx=0; phi_idx < tensor_shape.at(1); phi_idx++){
         for(int eta_idx =0; eta_idx < tensor_shape.at(2); eta_idx++){
           for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
             float mean = L2cacheData->normVec.at(var_idx).mean;
             float std = L2cacheData->normVec.at(var_idx).std;
             float min = L2cacheData->normVec.at(var_idx).min;
             float max = L2cacheData->normVec.at(var_idx).max;
             float nonstd_var = getCell(tau_idx, phi_idx, eta_idx, var_idx);
             float std_var = static_cast<float>((nonstd_var-mean)/std) ;
             if(std_var > max ){
               std_var = static_cast<float>(max);
             }
             else if (std_var < min){
               std_var = static_cast<float>(min);
             }
             getCell(tau_idx, phi_idx, eta_idx, var_idx) = std_var;
           } // end loop on var_idx
         } // end loop on eta_idx
       } // end loop on phi_idx
     } // end loop on tau_idx
   } // end function


   void L2TauNNProducer::FillL1TauVars(tensorflow::Tensor &cellGridMatrix, const std::vector<l1t::TauRef> allTaus){
     auto getCell = [&](int tau_idx, int phi_idx, int eta_idx, NNInputs NNInput_idx)-> float& {
       return cellGridMatrix.tensor<float, 4>()(tau_idx,phi_idx, eta_idx, static_cast<int>(NNInput_idx));
     };
     int nTaus = static_cast<int>(allTaus.size());
     for (int tau_idx = 0; tau_idx< nTaus; tau_idx++){
       for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
         for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
           getCell(tau_idx,phi_idx, eta_idx, NNInputs::l1Tau_pt)= allTaus[tau_idx]->polarP4().pt();
           getCell(tau_idx,phi_idx, eta_idx, NNInputs::l1Tau_eta)= allTaus[tau_idx]->polarP4().eta();
           getCell(tau_idx,phi_idx, eta_idx, NNInputs::l1Tau_hwIso)= allTaus[tau_idx]->hwIso();
         } //end loop phi
       } // end loop eta
     } // end loop tau
   } // end function

   void L2TauNNProducer::FillCaloRecHits(tensorflow::Tensor &cellGridMatrix, const std::vector<l1t::TauRef> allTaus, const caloRecHitCollections& caloRecHits){

     auto getCell = [&](int tau_idx, int phi_idx, int eta_idx, NNInputs NNInput_idx)-> float& {
       return cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, static_cast<int>(NNInput_idx));
     };

     int nTaus = static_cast<int>(allTaus.size());
     for (int tau_idx = 0; tau_idx< nTaus; tau_idx++){
      float tauEta = allTaus[tau_idx]->polarP4().eta();

      // caorechit_EE
      for (auto & caloRecHit_ee : *caloRecHits.ee){
        if(caloRecHit_ee.energy()<=0) continue;
        auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ee.id())->getPosition();
        float eeCalEta = position.eta();

        float eeCalEn = caloRecHit_ee.energy();
        float eeCalChi2 = caloRecHit_ee.chi2();
        if(reco::deltaR2(position,allTaus[tau_idx]->polarP4())<dR2_max){
          float deta = eeCalEta-tauEta;
          int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
          float dphi = reco::deltaPhi(position,allTaus[tau_idx]->polarP4());
          int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum)+=(eeCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSize)+=1.;
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergyStdDev)+=(eeCalEn * eeCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaEta)+=(deta*eeCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaPhi)+=(dphi*eeCalEn);
          if(eeCalChi2>=0){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalChi2)+=(eeCalChi2*eeCalEn);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySumForPositiveChi2)+=(eeCalEn);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSizeForPositiveChi2)+=1.;
          }
        }
      } // end of loop over calorechit_ee

      // caorechit_EB
      for (auto & caloRecHit_eb : *caloRecHits.eb){
        if(caloRecHit_eb.energy()<=0) continue;
        auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_eb.id())->getPosition();
        float ebCalEta = position.eta();
        float ebCalEn = caloRecHit_eb.energy();
        float ebCalChi2 = caloRecHit_eb.chi2();
        if(reco::deltaR2(position,allTaus[tau_idx]->polarP4())<dR2_max){
          float deta = ebCalEta-tauEta;
          int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
          float dphi = reco::deltaPhi(position,allTaus[tau_idx]->polarP4());
          int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum)+=(ebCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSize)+=1.;
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergyStdDev)+=(ebCalEn * ebCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaEta)+=(deta*ebCalEn);
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaPhi)+=(dphi*ebCalEn);
          if(ebCalChi2>=0){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalChi2)+=(ebCalChi2*ebCalEn);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySumForPositiveChi2)+=(ebCalEn);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSizeForPositiveChi2)+=1.;
          }
        }
      } // end of loop over calorechit_eb

      // normalize to sum and define stdDev
      for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
        for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
          if( getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum)>0. ){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaEta) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaEta)/getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaPhi) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalDeltaPhi)/getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum);
          }
          if( getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySumForPositiveChi2)>0.){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalChi2) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalChi2)/ getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySumForPositiveChi2);
          }
          if(getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSize)>1.){
            // (stdDev - (enSum*enSum)/size) / (size-1)
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergyStdDev) = (getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergyStdDev)  -  ( getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum) * getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergySum) ) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSize)   ) / (getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalSize) -1 );
          }
          else{
             getCell(tau_idx, phi_idx, eta_idx, NNInputs::EcalEnergyStdDev) = 0.;
          }
        } // end loop on phi
      } // end loop on eta


      // caorechit_HBHE
      for (auto & caloRecHit_hbhe : *caloRecHits.hbhe){
        if(caloRecHit_hbhe.energy()<=0) continue;
        auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_hbhe.id())->getPosition();
        float hbheCalEta = position.eta();
        float hbheCalEn = caloRecHit_hbhe.energy();
        float hbheCalChi2 = caloRecHit_hbhe.chi2();
        if(reco::deltaR2(position,allTaus[tau_idx]->polarP4())<dR2_max){
          float deta = hbheCalEta-tauEta;
          int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
          float dphi = reco::deltaPhi(position,allTaus[tau_idx]->polarP4());
          int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum)+=(hbheCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergyStdDev)+=(hbheCalEn * hbheCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSize)+=1.; //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaEta)+=(deta*hbheCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaPhi)+=(dphi*hbheCalEn); //
          if(hbheCalChi2>=0){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalChi2)+=(hbheCalChi2*hbheCalEn); //
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySumForPositiveChi2)+=(hbheCalEn); //
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSizeForPositiveChi2)+=1.; //
          }
        }
      } // end of loop over calorechit_hbhe

      // caorechit_HO
      for (auto & caloRecHit_ho : *caloRecHits.ho){
        if(caloRecHit_ho.energy()<=0) continue;
        auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ho.id())->getPosition();
        float hoCalEta = position.eta();
        float hoCalEn = caloRecHit_ho.energy();
        if(reco::deltaR2(position,allTaus[tau_idx]->polarP4())<dR2_max){
          float deta = hoCalEta-tauEta;
          int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
          float dphi = reco::deltaPhi(position,allTaus[tau_idx]->polarP4());
          int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum)+=(hoCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergyStdDev)+=(hoCalEn * hoCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSize)+=1.; //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaEta)+=(deta*hoCalEn); //
          getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaPhi)+=(dphi*hoCalEn); //
        }
      } // end of loop over calorechit_ho

      // normalize to sum and define stdDev
      for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
        for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
          if( getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum)>0. ){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaEta) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaEta)/getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaPhi) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalDeltaPhi)/getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum);
          }
          if( getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySumForPositiveChi2)>0.){
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalChi2) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalChi2)/ getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySumForPositiveChi2);
          }
          if(getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSize)>1.){
            // (stdDev - (enSum*enSum)/size) / (size-1)
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergyStdDev) = (getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergyStdDev)  -  ( getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum) * getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergySum) ) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSize)   ) / (getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalSize) -1 );
          }
          else{
             getCell(tau_idx, phi_idx, eta_idx, NNInputs::HcalEnergyStdDev) = 0.;
          }
        } // end loop on phi
      } // end loop on eta

    } // end loop over tau

  } // end function

   std::vector<int> L2TauNNProducer::GetVtxGood(const ZVertexHeterogeneous& patavertices, const PixelTrackHeterogeneous& patatracks, const std::vector<int> TrackGood){
     auto const &patatrack_tsoa = *patatracks.get();
     auto maxTracks = patatrack_tsoa.stride();
     auto const &patavtx_soa = *patavertices.get();
     int nv = patavtx_soa.nvFinal;
     std::vector<int> VtxGood;
     VtxGood.reserve(nv);

     float fractionSumPt2_ = 0.3;
     float minSumPt2_=0.;
     double track_pT_min_ = 1.;
     double track_pT_max_ = 20. ;
     //double track_prob_min_ = -1.; // keep it for future changes
     double track_chi2_max_ = 20.;
     unsigned int maxVtx_ = 100;
     std::vector<double> maxChi2_;
     std::vector<double> pTSquaredSum(nv);

     // loop on vertices and on associated tracks to vertices
     for (int j = nv - 1; j >= 0; --j){
       std::vector<int> trk_ass_to_vtx ;
       auto vtx_idx = patavtx_soa.sortInd[j];
       assert(vtx_idx < nv);
       for(int trk_idx = 0 ; trk_idx < maxTracks; trk_idx++){
         int vtx_ass_to_track = patavtx_soa.idv[trk_idx];
         if (vtx_ass_to_track == int16_t(vtx_idx))
          trk_ass_to_vtx.push_back(trk_idx);
       }
       auto nt = trk_ass_to_vtx.size();
       if (nt == 0) {
        std::cout << "vertex " << vtx_idx << " with no tracks..." << std::endl;
        continue;
       }
       // remove outliers
       if (nt < 2) {
         trk_ass_to_vtx.clear();
         continue;
       }
      for (auto& trk_idx : trk_ass_to_vtx){
        int vtx_ass_to_track = patavtx_soa.idv[trk_idx];
        if(vtx_ass_to_track!=vtx_idx) continue;
           double patatrackPt = patatrack_tsoa.pt[trk_idx];
           if (patatrackPt < track_pT_min_)
             continue;
           if (patatrack_tsoa.chi2(trk_idx) > track_chi2_max_)
             continue;
           if (patatrackPt > track_pT_max_){
             patatrackPt = track_pT_max_;
           }
           pTSquaredSum.at(vtx_idx) += patatrackPt * patatrackPt;
         }
       }

     auto minFOM_fromFrac = pTSquaredSum.at(patavtx_soa.sortInd[nv-1]) * fractionSumPt2_;
     if(minFOM_fromFrac==0){
       for (int j = nv - 1; j >= 0; --j){
         minFOM_fromFrac = pTSquaredSum.at(patavtx_soa.sortInd[j]) * fractionSumPt2_;
         if (minFOM_fromFrac!=0){
           break;
         }
       }
     }

     for (int j = nv - 1; j >= 0; --j){
       auto idx = patavtx_soa.sortInd[j];

       if (VtxGood.size() >= maxVtx_) {
         break;
       }
       if (pTSquaredSum[idx] >= minFOM_fromFrac && pTSquaredSum[idx] > minSumPt2_){
         VtxGood.push_back(idx);
       }
     }
     return VtxGood;
   }


   std::vector<float> L2TauNNProducer::impactParameter(int it, const PixelTrackHeterogeneous& patatracks, float patatrackPhi,const reco::BeamSpot &beamspot, const MagneticField *magfi ){
     auto const &patatrack_tsoa = *patatracks.get();
     auto const &fit = patatrack_tsoa.stateAtBS;
     std::vector<float> impactParameters;
     /* dxy e dz */
     //Rfit::Vector5d ipar, opar;
     //Rfit::Matrix5d icov, ocov;
     riemannFit::Vector5d ipar, opar;
     riemannFit::Matrix5d icov, ocov;
     fit.copyToDense(ipar, icov, it);
     //Rfit::transformToPerigeePlane(ipar, icov, opar, ocov);
     riemannFit::transformToPerigeePlane(ipar, icov, opar, ocov);
     LocalTrajectoryParameters lpar(opar(0), opar(1), opar(2), opar(3), opar(4), 1.);
     float sp = std::sin(patatrackPhi);
     float cp = std::cos(patatrackPhi);
     Surface::RotationType Rotation(sp, -cp, 0, 0, 0, -1.f, cp, sp, 0);
     GlobalPoint BeamSpotPoint(beamspot.x0(), beamspot.y0(), beamspot.z0());
     Plane impPointPlane(BeamSpotPoint, Rotation);
     GlobalTrajectoryParameters gp(impPointPlane.toGlobal(lpar.position()), impPointPlane.toGlobal(lpar.momentum()), lpar.charge(), magfi);
     GlobalPoint vv = gp.position();
     math::XYZPoint pos(vv.x(), vv.y(), vv.z());
     GlobalVector pp = gp.momentum();
     math::XYZVector mom(pp.x(), pp.y(), pp.z());
     auto lambda = M_PI_2 - pp.theta();
     auto phi = pp.phi();
     float patatrackDxy = -vv.x()*std::sin(phi)+ vv.y()*std::cos(phi);
     float patatrackDz = (vv.z()*std::cos(lambda)- (vv.x()*std::cos(phi)+vv.y()*std::sin(phi))*std::sin(lambda))/std::cos(lambda);
     impactParameters.push_back(patatrackDxy);
     impactParameters.push_back(patatrackDz);
     return impactParameters;
   }


   void L2TauNNProducer::FillPatatracks(tensorflow::Tensor &cellGridMatrix, const std::vector<l1t::TauRef> allTaus, const PixelTrackHeterogeneous& patatracks, const ZVertexHeterogeneous& patavertices,const reco::BeamSpot &beamspot, const MagneticField *magfi){

     auto getCell = [&](int tau_idx, int phi_idx, int eta_idx, NNInputs NNInput_idx)-> float& {
       return cellGridMatrix.tensor<float, 4>()(tau_idx,phi_idx, eta_idx, static_cast<int>(NNInput_idx));
     };
     int nTaus = static_cast<int>(allTaus.size());
     for (int tau_idx = 0; tau_idx< nTaus; tau_idx++){
       float tauEta = allTaus[tau_idx]->polarP4().eta();
       float tauPhi = allTaus[tau_idx]->polarP4().phi();

       auto const &patatrack_tsoa = *patatracks.get();
       auto maxTracks = patatrack_tsoa.stride();
       auto const &patavtx_soa = *patavertices.get();
       auto const *quality = patatrack_tsoa.qualityData();

       std::vector<int> TrackGood;
       for (int32_t it = 0; it < maxTracks; ++it) {
          //if(patatrack_tsoa.m_quality[it]<trackQuality::loose)
          //  continue;
          auto q = quality[it];
          if (q != pixelTrack::Quality::loose)
            continue;  // FIXME
          auto nHits = patatrack_tsoa.nHits(it);
          if (nHits == 0)
            break;
          if (nHits < 0)
            continue;
          TrackGood.push_back(it);
       }

       std::vector<int> VtxGood = GetVtxGood(patavertices,patatracks, TrackGood);



       for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
         for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
           getCell(tau_idx, phi_idx, eta_idx, NNInputs::nVertices)= VtxGood.size();
         }
       }

       for (auto& it : TrackGood) {
          //if(patatrack_tsoa.m_quality[it]<trackQuality::loose  ) continue;
          float patatrackPt = patatrack_tsoa.pt[it];
          if(patatrackPt<=0) continue;
          float patatrackPhi = patatrack_tsoa.phi(it);
          float patatrackEta = patatrack_tsoa.eta(it);
          float patatrackCharge = patatrack_tsoa.charge(it);
          float patatrackChi2OverNdof = patatrack_tsoa.chi2(it);
          auto nHits = patatrack_tsoa.nHits(it);
          if (nHits <= 0) continue;
          int patatrackNdof = 2 * nHits - 5;

          int vtx_idx_assTrk = patavtx_soa.idv[it];
          if(reco::deltaR2(patatrackEta,patatrackPhi,tauEta,tauPhi)<dR2_max){
            float deta = patatrackEta-tauEta;
            int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
            float dphi = reco::deltaPhi(patatrackPhi,tauPhi);
            int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum)+=(patatrackPt);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackSize)+=1.;
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChargeSum)+= patatrackCharge;
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaEta)+=(deta*patatrackPt);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaPhi)+=(dphi*patatrackPt);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChi2OverNdof)+=(patatrackChi2OverNdof*patatrackPt);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackNdof)+= (patatrackNdof*patatrackPt);
            std::vector<float> impactParameters = impactParameter(it, patatracks, patatrackPhi,beamspot, magfi );
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDxy)+= (impactParameters.at(0)*patatrackPt);
            getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDz)+= (impactParameters.at(1)*patatrackPt);
            //if(vtx_idx_assTrk>0 ){
            if((std::find(VtxGood.begin(), VtxGood.end(), vtx_idx_assTrk) != VtxGood.end())){
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSumWithVertex)+=(patatrackPt);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackSizeWithVertex)+=1.;
            }
          } // end if deltaR
        } // end of loop over patatracks

        // normalize to sum and define stdDev
        for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
          for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
            if(getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum)>0.){
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaEta) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaEta) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaPhi) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaPhi) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChi2OverNdof) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChi2OverNdof) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackNdof) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackNdof) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDxy) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDxy) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
              getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDz) = getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDz) / getCell(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum);
            } // end if
          } // end loop on phi
        } // end loop on eta
      } // end loop on taus
   } // end function


   std::vector<float> L2TauNNProducer::GetTauScore(const tensorflow::Tensor &cellGridMatrix){
     std::vector<tensorflow::Tensor> pred_tensor;
     std::vector<float> pred_vector;
     tensorflow::run(L2cacheData->session, {{inputTensorName_, cellGridMatrix}}, {outputTensorName_}, &pred_tensor);
     std::vector<int> tensor_shape;
     for(int d=0; d<cellGridMatrix.shape().dims(); d++) {
         tensor_shape.push_back(cellGridMatrix.shape().dim_size(d));
     }
     for(int tau_idx=0; tau_idx < tensor_shape.at(0); tau_idx++){
       pred_vector.push_back(pred_tensor[0].matrix<float>()(tau_idx, 0));
     }

     return pred_vector;
   }


   void L2TauNNProducer::produce(edm::Event& event, const edm::EventSetup& eventsetup){

     std::map<size_t, l1t::TauVectorRef> TauCollectionMap;   // map of l1Taucollections - input_index
     std::vector<std::string> CollectionNamesMap;
     CollectionNamesMap.reserve(L1TauDesc.size());
     l1t::TauVectorRef allTaus;

     for (size_t inp_idx = 0; inp_idx<L1TauDesc.size(); inp_idx++){
       CollectionNamesMap.push_back(L1TauDesc[inp_idx].CollectionName) ;

       edm::Handle<trigger::TriggerFilterObjectWithRefs> l1TriggeredTaus;
       event.getByToken(L1TauDesc[inp_idx].input_token, l1TriggeredTaus);

       l1t::TauVectorRef l1Taus;
       l1TriggeredTaus->getObjects(trigger::TriggerL1Tau, l1Taus);
       TauCollectionMap.insert(std::pair<size_t, l1t::TauVectorRef>(inp_idx, l1Taus));

       for(auto& i : l1Taus){

         if(std::find(allTaus.begin(), allTaus.end(), i)!= allTaus.end()) continue;
         allTaus.push_back(i);
       }

     }

     edm::Handle<EcalRecHitCollection> ebHandle;
     edm::Handle<EcalRecHitCollection> eeHandle;
     for (std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i = ecal_tokens.begin(); i != ecal_tokens.end(); i++) {
       edm::Handle<EcalRecHitCollection> ec_tmp;
       event.getByToken(*i, ec_tmp);
       if (ec_tmp->empty())
         continue;
       if ((ec_tmp->begin()->detid()).subdetId() == EcalBarrel) {
         ebHandle = ec_tmp;
       }
       else if ((ec_tmp->begin()->detid()).subdetId() == EcalEndcap) {
         eeHandle = ec_tmp;
       }
     }

     std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i;
     for (i = ecal_tokens.begin(); i != ecal_tokens.end(); i++) {
       edm::Handle<EcalRecHitCollection> ec;
       event.getByToken(*i, ec);
     }
     edm::Handle<HBHERecHitCollection> hbhe;
     event.getByToken(hbhe_token, hbhe);
     edm::Handle<HORecHitCollection> ho;
     event.getByToken(ho_token, ho);
     edm::Handle<PixelTrackHeterogeneous> pataTracks;
     event.getByToken(pataTracks_token, pataTracks);
     edm::Handle<ZVertexHeterogeneous> pataVertices;
     event.getByToken(pataVertices_token, pataVertices);
     edm::Handle<reco::BeamSpot> bsHandle;
     event.getByToken(BeamSpot_token, bsHandle);
     edm::ESHandle<MagneticField> fieldESH;
     eventsetup.get<IdealMagneticFieldRecord>().get(fieldESH);
     edm::ESHandle<CaloGeometry> Geometry = eventsetup.getHandle(Geometry_token);
     caloRecHitCollections caloRecHits;
     caloRecHits.hbhe= &*hbhe;
     caloRecHits.ho= &*ho;
     caloRecHits.eb= &*ebHandle;
     caloRecHits.ee= &*eeHandle;
     caloRecHits.Geometry = &*Geometry;
     int evt_id = event.id().event();

     int nTaus = allTaus.size();
     tensorflow::Tensor cellGridMatrix(tensorflow::DT_FLOAT, { nTaus, nCellEta, nCellPhi, nVars });
     for (int tau_idx = 0; tau_idx<nTaus ; tau_idx++){
       for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
         for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
           for (int var_idx = 0; var_idx<nVars ; var_idx++){
             cellGridMatrix.tensor<float, 4>()(tau_idx,phi_idx, eta_idx,var_idx) = 0.;
           }
         }
       }
     }
     FillL1TauVars(cellGridMatrix, allTaus);

     FillCaloRecHits(cellGridMatrix, allTaus, caloRecHits);

     FillPatatracks(cellGridMatrix, allTaus, *pataTracks,*pataVertices,*bsHandle, fieldESH.product());

     bool printoutevt=false;
     if(debugLevel_>4){
       if(evt_id==135505){
         printoutevt=true;
       }
     }
     standardizeTensor(cellGridMatrix); // do not remove this when eliminating debugging lines!!
     checknan(cellGridMatrix,printoutevt, debugLevel_);

     std::vector<float> tau_score = GetTauScore(cellGridMatrix);

     for(size_t inp_idx = 0; inp_idx<L1TauDesc.size(); inp_idx++){
       std::vector<float> outcomes_;
       for(int tau_idx = 0; tau_idx < nTaus; tau_idx++){
         std::cout <<  evt_id << " \t " << (allTaus[tau_idx])->polarP4().pt() << " \t "<< tau_score.at(tau_idx) << std::endl;
         //std::cout << "evt " << evt_id << " tau Pt " << (allTaus[tau_idx])->polarP4().pt() << " score "<< tau_score.at(tau_idx) << std::endl;
         if(std::find( TauCollectionMap.at(inp_idx).begin(), TauCollectionMap.at(inp_idx).end() ,allTaus[tau_idx]) != TauCollectionMap.at(inp_idx).end()){
           outcomes_.push_back(tau_score.at(tau_idx));
         }
       }

       auto output = std::make_unique<std::vector<float>>(outcomes_);
       event.put(std::move(output),CollectionNamesMap.at(inp_idx));
     }
   }
  }
using L2TauNNProducer = tau_hlt::L2TauNNProducer;
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L2TauNNProducer);
