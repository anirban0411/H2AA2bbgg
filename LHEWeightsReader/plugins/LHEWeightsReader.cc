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

class LHEWeightsReader : public edm::EDAnalyzer {
   public:
      explicit LHEWeightsReader(const edm::ParameterSet&);
      ~LHEWeightsReader();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endRun(edm::Run const&, edm::EventSetup const&) override;

      edm::InputTag lheProduct_;
      edm::EDGetTokenT<LHERunInfoProduct> lheToken_;
};


LHEWeightsReader::LHEWeightsReader(const edm::ParameterSet& iConfig):
lheProduct_( iConfig.getParameter<edm::InputTag>("lheProduct") )
{
   lheToken_ = consumes <LHERunInfoProduct,edm::InRun> (edm::InputTag(lheProduct_));
   return;
}


LHEWeightsReader::~LHEWeightsReader()
{
}


void
LHEWeightsReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}


void
LHEWeightsReader::endRun(edm::Run const& iRun, edm::EventSetup const&)
{

   using namespace edm;

   edm::Handle<LHERunInfoProduct> run;
   typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

   iRun.getByToken(lheToken_, run);
   LHERunInfoProduct myLHERunInfoProduct = *(run.product());
   for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
     std::cout << iter->tag() << std::endl;
     std::vector<std::string> lines = iter->lines();
     for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
     }
   }

}

void
LHEWeightsReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(LHEWeightsReader);
