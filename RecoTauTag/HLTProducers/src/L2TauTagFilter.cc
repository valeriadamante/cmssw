/*
Outcomes producedr for L2 hadronic tau selection
*/
// system include files
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <math.h>
#include "Compression.h"
#include "TMath.h"
// user include files
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"


namespace tau_hlt{

  class L2TauTagFilter : public edm::stream::EDFilter<> {

  public:
    explicit L2TauTagFilter(const edm::ParameterSet&);
    ~L2TauTagFilter() override {}
    static void fillDescriptions(edm::ConfigurationDescriptions&);

    bool filter(edm::Event& event, const edm::EventSetup& eventsetup) ;

  private:


  private:
    int debugLevel;
    std::string processName;
    //std::vector<float> l2Outcomes;
    edm::EDGetTokenT<std::vector<float>> l2Outcomes_token;
    double discr_threshold;

  };

  void L2TauTagFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<int>("debugLevel", 0)->setComment("set debug level for printing out info");
    desc.add<std::string>("processName","")->setComment("Name of the process");
    desc.add<edm::InputTag>("L2outcomes", edm::InputTag(""))->setComment("L2 CNN outcomes");
    desc.add<double>("discr_threshold", 0.12267940863785043)->setComment("value of discriminator threshold");
    descriptions.addWithDefaultLabel(desc);
  }

  L2TauTagFilter::L2TauTagFilter(const edm::ParameterSet& cfg):
        debugLevel(cfg.getParameter<int>("debugLevel")),
        processName(cfg.getParameter<std::string>("processName")),
        l2Outcomes_token(consumes<std::vector<float>>(cfg.getParameter<edm::InputTag>("L2outcomes"))),
        discr_threshold(cfg.getParameter<double>("discr_threshold"))
        {

        }


   bool L2TauTagFilter::filter(edm::Event& event, const edm::EventSetup& eventsetup){
     bool result = false;
     int nTauPassed = 0;
     edm::Handle<std::vector<float>> l2Outcomes;
     event.getByToken(l2Outcomes_token, l2Outcomes);


     //edm::Handle<trigger::TriggerFilterObjectWithRefs> l1TriggeredTaus;
     //event.getByToken(L1TauDesc[inp_idx].input_token, l1TriggeredTaus);
     int evt_id = event.id().event();

     for(auto& outcome : *l2Outcomes){
        if(outcome > discr_threshold){
          nTauPassed++;
        }
        if(nTauPassed == 2){
          //std::cout << "evt " << evt_id << " has at least two taus" << std::endl;
          return true;
        }
     }
     if(nTauPassed >= 2 ){
       std::cout << "evt " << evt_id << " has at least two taus" << std::endl;
       result = true;
     }

     return result;
   }
  }
using L2TauTagFilter = tau_hlt::L2TauTagFilter;
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L2TauTagFilter);
