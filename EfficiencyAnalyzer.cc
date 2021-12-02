#include <map>
#include <string>
#include <iostream>


#include "TH1.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

using namespace std;

class EfficiencyAnalyzer : public edm::EDAnalyzer {

public:
   explicit EfficiencyAnalyzer(const edm::ParameterSet&);
   ~EfficiencyAnalyzer();
  
private:

   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
  
   // simple map to contain all histograms; 
   // histograms are booked in the beginJob() 
   // method
   std::map<std::string,TH1F*> histContainer_; 
   // ----------member data ---------------------------     
   edm::EDGetTokenT<pat::MuonCollection> muonCollToken;
   edm::EDGetTokenT<pat::ElectronCollection> elecCollToken;
   edm::EDGetTokenT<pat::TauCollection> tauCollToken;
   edm::EDGetTokenT<pat::JetCollection> jetCollToken;
   edm::EDGetTokenT<pat::METCollection> metCollToken;


   // input tags 
   edm::InputTag muonSrc_;
   edm::InputTag elecSrc_;
   edm::InputTag tauSrc_;
   edm::InputTag jetSrc_;
   edm::InputTag metSrc_;
};


EfficiencyAnalyzer::EfficiencyAnalyzer(const edm::ParameterSet& iConfig):

   histContainer_(),
   muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
   elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elecSrc")),
   tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc")),
   jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc")),
   metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc")){

   muonCollToken = consumes<pat::MuonCollection>(muonSrc_);
   elecCollToken = consumes<pat::ElectronCollection>(elecSrc_);
   tauCollToken = consumes<pat::TauCollection>(tauSrc_);
   jetCollToken = consumes<pat::JetCollection>(jetSrc_);
   metCollToken = consumes<pat::METCollection>(metSrc_);



}

EfficiencyAnalyzer::~EfficiencyAnalyzer(){
}

void
EfficiencyAnalyzer::analyze(const edm::Event& iEvent, 
                   const edm::EventSetup& iSetup){

   // get pat muon collection 
   edm::Handle< std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonCollToken, muons);


   // get pat electron collection 
   edm::Handle< std::vector<pat::Electron>> electrons;
   iEvent.getByToken(elecCollToken, electrons);


   // get pat tau collection 
   edm::Handle< std::vector<pat::Tau>> taus;
   iEvent.getByToken(tauCollToken, taus);


   // get pat jet collection 
   edm::Handle< std::vector<pat::Jet>> jets;
   iEvent.getByToken(jetCollToken, jets);


    // get pat met collection 
   edm::Handle< std::vector<pat::MET>> met;
   iEvent.getByToken(metCollToken, met);

   int cortetau=0;
   int cortetau2=0;
   int cortetauden=0;
   int cortejet=0;
   int cortejet2=0;
   int cortejetden=0;
   histContainer_["eventos"]->Fill(1); 


   // loop over taus
   for (auto itt = taus->cbegin(); itt != taus->cend(); ++itt) {
      histContainer_["tauPt"] ->Fill(itt->pt());
      histContainer_["tauEta"]->Fill(itt->eta());
      histContainer_["tauPhi"]->Fill(itt->phi());    
     
      //cout<<"Tau loop itt: "<<"Pt: "<<itt->pt()<<" Eta: "<<itt->eta()<<" Phi: "<<itt->phi()<<" size: "<<taus->size()<<endl;

      if( itt->pt()>20 && fabs(itt->eta())<2.2 ){
         cortetau++;
         histContainer_["eventos"]->Fill(-1); 

         if( itt->pt()>30 && fabs(itt->eta())<2.1 ){
         cortetau2++;
         histContainer_["eventos"]->Fill(-2);
         }
      }

      if( itt->pt()>30 && fabs(itt->eta())<2.1 ){
         cortetauden++;
         histContainer_["eventos"]->Fill(-5);
         }
   
   }


   // loop over jets
   for (auto itj = jets->cbegin(); itj != jets->cend(); ++itj) {
      histContainer_["jetPt"] ->Fill(itj->pt());
      histContainer_["jetEta"]->Fill(itj->eta());
      histContainer_["jetPhi"]->Fill(itj->phi());    

      //cout<<"Jet loop itj: "<<"Pt: "<<itj->pt()<<" Eta: "<<itj->eta()<<" Phi: "<<itj->phi()<<" size: "<<jets->size()<<endl;
   
   }

//cout<<" Leadingjet0  "<<"Pt: "<<(*jets)[0].pt()<<" Eta: "<<(*jets)[0].eta()<<" Phi: "<<(*jets)[0].phi()<<endl;
//cout<<" Leadingjet1  "<<"Pt: "<<(*jets)[1].pt()<<" Eta: "<<(*jets)[1].eta()<<" Phi: "<<(*jets)[1].phi()<<endl;

   if(jets->size()>1){ 
      histContainer_["dijetMass"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
      histContainer_["jet0Pt"]->Fill((*jets)[0].pt());
      histContainer_["jet1Pt"]->Fill((*jets)[1].pt()); 
      histContainer_["jet0Eta"]->Fill((*jets)[0].eta());
      histContainer_["jet1Eta"]->Fill((*jets)[1].eta()); 
      histContainer_["jet0Phi"]->Fill((*jets)[0].phi());
      histContainer_["jet1Phi"]->Fill((*jets)[1].phi()); 




      if((*jets)[0].eta()*(*jets)[1].eta()<0   &&  fabs((*jets)[0].eta() - (*jets)[1].eta()) > 4.2  ) {

         if(  ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 750 && (*jets)[0].pt()>40 && (*jets)[1].pt()>40 ) { 
            //  histContainer_["dijetMass50_850"]->Fill( ((*jets)[0].p4()+(*jets)[1].p4()).mass()  );    
            //  histContainer_["jet0Pt_50_850"]->Fill((*jets)[0].pt());
            //  histContainer_["jet1Pt_50_850"]->Fill((*jets)[1].pt());
            cortejet++;
            histContainer_["eventos"]->Fill(-3); 
            if(  ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850 && (*jets)[0].pt()>50 && (*jets)[1].pt()>50 ) { 
               cortejet2++;
               histContainer_["eventos"]->Fill(-4); 
            }
         }
      }
   

      if(  ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850 && (*jets)[0].pt()>50 && (*jets)[1].pt()>50 ) { 
         cortejetden++;
         histContainer_["eventos"]->Fill(-6); 
      }



   }

   // loop over met
   for (auto itm = met->cbegin(); itm != met->cend(); ++itm) {
      histContainer_["metPt"] ->Fill(itm->pt());
      histContainer_["metEta"]->Fill(itm->eta());
      histContainer_["metPhi"]->Fill(itm->phi());    
    // cout<<"Met loop itm: "<<"Pt: "<<itm->pt()<<" Eta: "<<itm->eta()<<" Phi: "<<itm->phi()<<" size: "<<met->size()<<endl;

}



if (cortetau>0 && cortejet> 0) {

   histContainer_["eventos"]->Fill(2); 

   histContainer_["tauPt_cut"] ->Fill((*taus)[0].pt());
   histContainer_["tauEta_cut"]->Fill((*taus)[0].eta());
   histContainer_["tauPhi_cut"]->Fill((*taus)[0].phi()); 
   histContainer_["dijetMass_cut"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_cut"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_cut"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_cut"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_cut"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_cut"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_cut"]->Fill((*jets)[1].phi()); 
 
   histContainer_["metPt_cut"] ->Fill((*met)[0].pt());
   histContainer_["metEta_cut"]->Fill((*met)[0].eta());
   histContainer_["metPhi_cut"]->Fill((*met)[0].phi());
}


if (cortetau2>0 && cortejet2> 0) {

   histContainer_["eventos"]->Fill(3); 

   histContainer_["tauPt_cut2"] ->Fill((*taus)[0].pt());
   histContainer_["tauEta_cut2"]->Fill((*taus)[0].eta());
   histContainer_["tauPhi_cut2"]->Fill((*taus)[0].phi()); 
   histContainer_["dijetMass_cut2"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_cut2"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_cut2"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_cut2"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_cut2"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_cut2"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_cut2"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_cut2"] ->Fill((*met)[0].pt());
   histContainer_["metEta_cut2"]->Fill((*met)[0].eta());
   histContainer_["metPhi_cut2"]->Fill((*met)[0].phi());

}

if (cortetauden>0 && cortejetden> 0) {

   histContainer_["eventos"]->Fill(4); 

   histContainer_["tauPt_den"] ->Fill((*taus)[0].pt());
   histContainer_["tauEta_den"]->Fill((*taus)[0].eta());
   histContainer_["tauPhi_den"]->Fill((*taus)[0].phi()); 
   histContainer_["dijetMass_den"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_den"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_den"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_den"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_den"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_den"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_den"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_den"] ->Fill((*met)[0].pt());
   histContainer_["metEta_den"]->Fill((*met)[0].eta());
   histContainer_["metPhi_den"]->Fill((*met)[0].phi());

}

   // Multiplicity
 
   histContainer_["tauMult"]->Fill(taus->size() );
   histContainer_["jetMult"]->Fill(jets->size() );
   histContainer_["metMult"]->Fill(met->size() );

}



void 
EfficiencyAnalyzer::beginJob()
{
   // register to the TFileService
   edm::Service<TFileService> fs;

//  number of events

   histContainer_["eventos"]=fs->make<TH1F>("eventos",   "numero eventos", 40, -10,  10);


   // book histograms for Tau


   histContainer_["tauPt"]=fs->make<TH1F>("tauPt",   "tau Pt", 300, 0, 600);
   histContainer_["tauPt_cut"]=fs->make<TH1F>("tauPt_cut",   "tau Pt", 300, 0, 600);
   histContainer_["tauPt_cut2"]=fs->make<TH1F>("tauPt_cut2",   "tau Pt", 300, 0, 600);
   histContainer_["tauPt_den"]=fs->make<TH1F>("tauPt_den",   "tau Pt", 300, 0, 600);
   histContainer_["tauEta"]=fs->make<TH1F>("tauEta",   "tau Eta",100, -5,  5);
   histContainer_["tauEta_cut"]=fs->make<TH1F>("tauEta_cut",   "tau Eta",100, -5,  5);
   histContainer_["tauEta_cut2"]=fs->make<TH1F>("tauEta_cut2",   "tau Eta",100, -5,  5);
   histContainer_["tauEta_den"]=fs->make<TH1F>("tauEta_den",   "tau Eta",100, -5,  5);
   histContainer_["tauPhi"]=fs->make<TH1F>("tauPhi",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauPhi_cut"]=fs->make<TH1F>("tauPhi_cut",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauPhi_cut2"]=fs->make<TH1F>("tauPhi_cut2",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauPhi_den"]=fs->make<TH1F>("tauPhi_den",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauMult"]=fs->make<TH1F>("tauMult",   "tau multiplicity",     100, 0,  50);

   // book histograms for MET

   histContainer_["metPt"]=fs->make<TH1F>("metPt",   "met Pt", 300, 0, 600);
   histContainer_["metPt_cut"]=fs->make<TH1F>("metPt_cut",   "met Pt", 300, 0, 600);
   histContainer_["metPt_cut2"]=fs->make<TH1F>("metPt_cut2",   "met Pt", 300, 0, 600);
   histContainer_["metPt_den"]=fs->make<TH1F>("metPt_den",   "met Pt", 300, 0, 600);
   histContainer_["metEta"]=fs->make<TH1F>("metEta",   "met Eta",100, -5,  5);
   histContainer_["metEta_cut"]=fs->make<TH1F>("metEta_cut",   "met Eta",100, -5,  5);
   histContainer_["metEta_cut2"]=fs->make<TH1F>("metEta_cut2",   "met Eta",100, -5,  5);
   histContainer_["metEta_den"]=fs->make<TH1F>("metEta_den",   "met Eta",100, -5,  5);
   histContainer_["metPhi"]=fs->make<TH1F>("metPhi",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metPhi_cut"]=fs->make<TH1F>("metPhi_cut",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metPhi_cut2"]=fs->make<TH1F>("metPhi_cut2",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metPhi_den"]=fs->make<TH1F>("metPhi_den",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metMult"]=fs->make<TH1F>("metMult",   "met multiplicity",     100, 0,  50);


   // book histograms for Jets

   histContainer_["jetPt"]=fs->make<TH1F>("jetPt",   "jet Pt", 300, 0, 600);
   histContainer_["jet0Pt"]=fs->make<TH1F>("jet0pt",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet0Pt_cut"]=fs->make<TH1F>("jet0Pt_cut",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet0Pt_cut2"]=fs->make<TH1F>("jet0Pt_cut2",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet0Pt_den"]=fs->make<TH1F>("jet0Pt_den",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet1Pt"]=fs->make<TH1F>("jet1pt",   "jet1 Pt", 300, 0, 600);
   histContainer_["jet1Pt_cut"]=fs->make<TH1F>("jet1Pt_cut",   "jet1 Pt", 300, 0, 600);
   histContainer_["jet1Pt_cut2"]=fs->make<TH1F>("jet1Pt_cut2",   "jet1 Pt", 300, 0, 600);
   histContainer_["jet1Pt_den"]=fs->make<TH1F>("jet1Pt_den",   "jet1 Pt", 300, 0, 600);

   histContainer_["jetEta"]=fs->make<TH1F>("jetEta",   "jet Eta",  100, -5,  5);
   histContainer_["jet0Eta"]=fs->make<TH1F>("jet0Eta",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet0Eta_cut"]=fs->make<TH1F>("jet0Eta_cut",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet0Eta_cut2"]=fs->make<TH1F>("jet0Eta_cut2",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet0Eta_den"]=fs->make<TH1F>("jet0Eta_den",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet1Eta"]=fs->make<TH1F>("jet1Eta",   "jet1 Eta",  100, -5,  5);
   histContainer_["jet1Eta_cut"]=fs->make<TH1F>("jet1Eta_cut",   "jet1 Eta",  100, -5,  5);
   histContainer_["jet1Eta_cut2"]=fs->make<TH1F>("jet1Eta_cut2",   "jet1 Eta",  100, -5,  5);
   histContainer_["jet1Eta_den"]=fs->make<TH1F>("jet1Eta_den",   "jet1 Eta",  100, -5,  5);

   histContainer_["jetPhi"]=fs->make<TH1F>("jetPhi",   "jet Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi"]=fs->make<TH1F>("jet0Phi",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi_cut"]=fs->make<TH1F>("jet0Phi_cut",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi_cut2"]=fs->make<TH1F>("jet0Phi_cut2",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi_den"]=fs->make<TH1F>("jet0Phi_den",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi"]=fs->make<TH1F>("jet1Phi",   "jet1 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi_cut"]=fs->make<TH1F>("jet1Phi_cut",   "jet1 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi_cut2"]=fs->make<TH1F>("jet1Phi_cut2",   "jet1 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi_den"]=fs->make<TH1F>("jet1Phi_den",   "jet1 Phi",     100, -3.5, 3.5);

   histContainer_["jetMult"]=fs->make<TH1F>("jetMult",   "jet multiplicity",     100, 0,  50);
   //histContainer_["dijetMass50_850"]=fs->make<TH1F>("dijetMass850_50", "mass",    300,   0., 1500.);
   //histContainer_["jet0Pt_50_850"]=fs->make<TH1F>("jet0Pt_50_850",   "jet0 Pt", 300, 0, 600);
   //histContainer_["jet1Pt_50_850"]=fs->make<TH1F>("jet1Pt_50_850",   "jet1 Pt", 300, 0, 600);
   histContainer_["dijetMass"]=fs->make<TH1F>("dijetMass", "mass",    300,   0., 1500.);
   histContainer_["dijetMass_cut"]=fs->make<TH1F>("dijetMass_cut", "mass",    300,   0., 1500.);
   histContainer_["dijetMass_cut2"]=fs->make<TH1F>("dijetMass_cut2", "mass",    300,   0., 1500.);
   histContainer_["dijetMass_den"]=fs->make<TH1F>("dijetMass_den", "mass",    300,   0., 1500.);

          


    
}

void 
EfficiencyAnalyzer::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EfficiencyAnalyzer);

