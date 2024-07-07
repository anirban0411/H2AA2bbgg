//
// Original Author:  Lorenzo Viliani
//         Created:  Fri, 30 Nov 2018 09:20:09 GMT
//

// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

//
// class declaration
//

class LHEWeightsDumper : public edm::EDAnalyzer {
   public:
      explicit LHEWeightsDumper(const edm::ParameterSet&);
      ~LHEWeightsDumper();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endRun(edm::Run const&, edm::EventSetup const&) override;

      edm::InputTag lheProduct_;
      edm::EDGetTokenT<LHEEventProduct> lheToken_;
};


LHEWeightsDumper::LHEWeightsDumper(const edm::ParameterSet& iConfig):
lheProduct_( iConfig.getParameter<edm::InputTag>("lheProduct") )
{
   lheToken_ = consumes <LHEEventProduct,edm::InEvent> (edm::InputTag(lheProduct_));
   return;
}


LHEWeightsDumper::~LHEWeightsDumper()
{
}


void
LHEWeightsDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle<LHEEventProduct> event;
  
   iEvent.getByToken(lheToken_,event);

   for (unsigned int i=0; i<event->weights().size(); i++) {
     std::cout << "Weight " << i << " ID: " << event->weights()[i].id << " " <<  event->weights()[i].wgt << std::endl;
   }

}


void
LHEWeightsDumper::endRun(edm::Run const& iRun, edm::EventSetup const&)
{
}

void
LHEWeightsDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(LHEWeightsDumper);
