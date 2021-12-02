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

class Dijetmass : public edm::EDAnalyzer {

public:
   explicit Dijetmass(const edm::ParameterSet&);
   ~Dijetmass();
  
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


Dijetmass::Dijetmass(const edm::ParameterSet& iConfig):

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

Dijetmass::~Dijetmass(){
}

void
Dijetmass::analyze(const edm::Event& iEvent, 
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

   // fill pat muon histograms
   for (auto it = muons->cbegin(); it != muons->cend(); ++it) {
      histContainer_["muonPt"] ->Fill(it->pt());
      histContainer_["muonEta"]->Fill(it->eta());
      histContainer_["muonPhi"]->Fill(it->phi());
      cout<<"Muon loop it: "<<"Pt: "<<it->pt()<<" Eta: "<<it->eta()<<" Phi: "<<it->phi()<<" size: "<<muons->size()<<endl;

      if( it->pt()>20 && fabs(it->eta())<2.1 ){
         for (auto it2 = muons->cbegin(); it2 != muons->cend(); ++it2){
            if (it2 > it){
               // check only muon pairs of unequal charge 
               if( it->charge()*it2->charge()<0){
                  if( it2->pt()>20 && fabs(it2->eta())<2.1 ){
                     histContainer_["mumuMass"]->Fill((it->p4()+it2->p4()).mass());

                  }
               }
            }
         }
      }
   }

   // loop over electrons
   for (auto ite = electrons->cbegin(); ite != electrons->cend(); ++ite) {
      histContainer_["elePt"] ->Fill(ite->pt());
      histContainer_["eleEta"]->Fill(ite->eta());
      histContainer_["elePhi"]->Fill(ite->phi());    
      cout<<"Electron loop ite: "<<"Pt: "<<ite->pt()<<" Eta: "<<ite->eta()<<" Phi: "<<ite->phi()<<" size: "<<electrons->size()<<endl;
    }



   // loop over taus
   for (auto itt = taus->cbegin(); itt != taus->cend(); ++itt) {
      histContainer_["tauPt"] ->Fill(itt->pt());
      histContainer_["tauEta"]->Fill(itt->eta());
      histContainer_["tauPhi"]->Fill(itt->phi());    
      histContainer_["taunum"]->Fill(1); 
      cout<<"Tau loop itt: "<<"Pt: "<<itt->pt()<<" Eta: "<<itt->eta()<<" Phi: "<<itt->phi()<<" size: "<<taus->size()<<endl;



      if( itt->pt()>30 && fabs(itt->eta())<2.1 ){
         histContainer_["tauPt30eta2p1"] ->Fill(itt->pt());
         histContainer_["tauEta30eta2p1"]->Fill(itt->eta());
         histContainer_["tauPhi30eta2p1"]->Fill(itt->phi());
         histContainer_["taunum30eta2p1"]->Fill(1);
  
      }
   
   }


   // loop over jets
   for (auto itj = jets->cbegin(); itj != jets->cend(); ++itj) {
      histContainer_["jetPt"] ->Fill(itj->pt());
      histContainer_["jetEta"]->Fill(itj->eta());
      histContainer_["jetPhi"]->Fill(itj->phi());    

      cout<<"Jet loop itj: "<<"Pt: "<<itj->pt()<<" Eta: "<<itj->eta()<<" Phi: "<<itj->phi()<<" size: "<<jets->size()<<endl;
    /*
      if( itj->pt()>50 ){
         histContainer_["jetPt50"] ->Fill(itj->pt());
         histContainer_["jetEta50"]->Fill(itj->eta());
         histContainer_["jetPhi50"]->Fill(itj->phi()); 
         if(jets->size()>1){ 
            histContainer_["dijetmass50"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());    
         }


      }
   */
   }

cout<<" Leadingjet0  "<<"Pt: "<<(*jets)[0].pt()<<" Eta: "<<(*jets)[0].eta()<<" Phi: "<<(*jets)[0].phi()<<endl;
cout<<" Leadingjet1  "<<"Pt: "<<(*jets)[1].pt()<<" Eta: "<<(*jets)[1].eta()<<" Phi: "<<(*jets)[1].phi()<<endl;

if(jets->size()>1){ 
histContainer_["dijetmass"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
histContainer_["jet0"]->Fill((*jets)[0].pt());
histContainer_["jet1"]->Fill((*jets)[1].pt()); 
}


if(  jets->size()>1 && ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850 ) { 
histContainer_["dijetmass850"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());    
histContainer_["jet0_850"]->Fill((*jets)[0].pt());
histContainer_["jet1_850"]->Fill((*jets)[1].pt());
}


if(  jets->size()>1 && ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850 && (*jets)[0].pt()>50 && (*jets)[1].pt()>50 ) { 
histContainer_["dijetmass850_50"]->Fill( ((*jets)[0].p4()+(*jets)[1].p4()).mass()  );    
histContainer_["jet0_50_850"]->Fill((*jets)[0].pt());
histContainer_["jet1_50_850"]->Fill((*jets)[1].pt());
}


if(  jets->size()>1  && (*jets)[0].pt()>50 && (*jets)[1].pt()>50 )     { 
histContainer_["dijetmass50"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
histContainer_["jet0_50"]->Fill((*jets)[0].pt());
histContainer_["jet1_50"]->Fill((*jets)[1].pt());
            }
  

   // loop over met
   for (auto itm = met->cbegin(); itm != met->cend(); ++itm) {
      histContainer_["metPt"] ->Fill(itm->pt());
      histContainer_["metEta"]->Fill(itm->eta());
      histContainer_["metPhi"]->Fill(itm->phi());    
     cout<<"Met loop itm: "<<"Pt: "<<itm->pt()<<" Eta: "<<itm->eta()<<" Phi: "<<itm->phi()<<" size: "<<met->size()<<endl;

}



   // Multiplicity
   histContainer_["eleMult" ]->Fill(electrons->size());
   histContainer_["muonMult"]->Fill(muons->size() );
   histContainer_["tauMult"]->Fill(taus->size() );
   histContainer_["jetMult"]->Fill(jets->size() );
   histContainer_["metMult"]->Fill(met->size() );




   /*
cout<<" MET  "<<"Pt: "<<(*met)[0].pt()<<" Eta: "<<(*met)[0].eta()<<" Phi: "<<(*met)[0].phi()<<endl;
cout<<" Leadingjet0  "<<"Pt: "<<(*jets)[0].pt()<<" Eta: "<<(*jets)[0].eta()<<" Phi: "<<(*jets)[0].phi()<<endl;
cout<<" Leadingjet1  "<<"Pt: "<<(*jets)[1].pt()<<" Eta: "<<(*jets)[1].eta()<<" Phi: "<<(*jets)[1].phi()<<endl;

cout<<" METdentro  "<<"Pt: "<<(*met)[0].pt()<<" Eta: "<<(*met)[0].eta()<<" Phi: "<<(*met)[0].phi()<<endl;
cout<<" Leadingjet0 dentro  "<<"Pt: "<<(*jets)[0].pt()<<" Eta: "<<(*jets)[0].eta()<<" Phi: "<<(*jets)[0].phi()<<endl;
cout<<" Leadingjet1  dentro "<<"Pt: "<<(*jets)[1].pt()<<" Eta: "<<(*jets)[1].eta()<<" Phi: "<<(*jets)[1].phi()<<endl;


cout<< "corte Tau: "<<cortetau<<" corte Jet: "<<cortejet<<endl;

   cout<<"eureka"<<endl;




cout<<"Sin loop... jets->size(): "<<jets->size()<<" Dijetmass: " <<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
   cout<<"Jet loop itj: "<<"jetPt: "<<itj->pt()<<" jetEta: "<<itj->eta()<<" jetPhi: "<<itj->phi()<<endl;
         cout<<"Dentro if taupt> 30 Leadingjet0  "<<"jetPt: "<<(*jets)[0].pt()<<" jetEta: "<<(*jets)[0].eta()<<" jetPhi: "<<(*jets)[0].phi()<<endl;
         cout<<"Dentro if taupt> 30 Leadingjet1  "<<"jetPt: "<<(*jets)[1].pt()<<" jetEta: "<<(*jets)[1].eta()<<" jetPhi: "<<(*jets)[1].phi()<<endl;
  cout<<"Dentro if taupt> 30 ... jets->size(): "<<jets->size()<<" Dijetmass: "<<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
  cout<<"Dentro if taupt> 30 size >1... jets->size(): "<<jets->size()<<" Dijetmass: "<<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
cout<<"Met loop itm: "<<"metPt: "<<itm->pt()<<" metEta: "<<itm->eta()<<" metPhi: "<<itm->phi()<<endl;
 cout<<"final 1 jets->size(): "<<jets->size()<<endl;
 cout<<"((*jets)[0].p4()+(*jets)[1].p4()).mass(): "<<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
cout<<" Leadingjet0  "<<"jetPt: "<<(*jets)[0].pt()<<" jetEta: "<<(*jets)[0].eta()<<" jetPhi: "<<(*jets)[0].phi()<<endl;
cout<<" Leadingjet1  "<<"jetPt: "<<(*jets)[1].pt()<<" jetEta: "<<(*jets)[1].eta()<<" jetPhi: "<<(*jets)[1].phi()<<endl;
 cout<<"final 2 jets->size(): "<<jets->size()<<endl;
  cout<<"((*jets)[0].p4()+(*jets)[1].p4()).mass(): "<<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
 cout<<"final 3 jets->size(): "<<jets->size()<<endl;
  cout<<"((*jets)[0].p4()+(*jets)[1].p4()).mass(): "<<((*jets)[0].p4()+(*jets)[1].p4()).mass()<<endl;
*/
}



void 
Dijetmass::beginJob()
{
   // register to the TFileService
   edm::Service<TFileService> fs;


   histContainer_["mumuMass"]=fs->make<TH1F>("mumuMass", "mass",    90,   30., 120.);
   histContainer_["dijetmass"]=fs->make<TH1F>("dijetmass", "mass",    200,   0., 1000.);
      histContainer_["dijetmass30"]=fs->make<TH1F>("dijetmass30", "mass",    200,   0., 1000.);
         histContainer_["dijetmass50"]=fs->make<TH1F>("dijetmass50", "mass",    200,   0., 1000.);
                  histContainer_["dijetmass850"]=fs->make<TH1F>("dijetmass850", "mass",    200,   0., 1000.);
                                    histContainer_["dijetmass850_50"]=fs->make<TH1F>("dijetmass850_50", "mass",    200,   0., 1000.);
         histContainer_["dijetmass30di850"]=fs->make<TH1F>("dijetmass30di850", "mass",    200,   0., 1000.);
  



   // book histograms for Multiplicity:

   histContainer_["eleMult"]=fs->make<TH1F>("eleMult",   "electron multiplicity", 100, 0,  50);
   histContainer_["muonMult"]=fs->make<TH1F>("muonMult",   "muon multiplicity",     100, 0,  50);
   histContainer_["tauMult"]=fs->make<TH1F>("tauMult",   "tau multiplicity",     100, 0,  50);
   histContainer_["jetMult"]=fs->make<TH1F>("jetMult",   "jet multiplicity",     100, 0,  50);
   histContainer_["metMult"]=fs->make<TH1F>("metMult",   "met multiplicity",     100, 0,  50);
   histContainer_["taunum"]=fs->make<TH1F>("taunum1",   "numero taus antes de cortes", 100, 0,  50);
      histContainer_["taunum30eta2p1"]=fs->make<TH1F>("taunum30eta2p1",   "numero taus despues de cortes pt 30 GeV eta <2.1",     100, 0,  50);

   // book histograms for Pt:

   histContainer_["elePt"]=fs->make<TH1F>("elePt",   "electron Pt", 100, 0,  200);
   histContainer_["muonPt"]=fs->make<TH1F>("muonPt",   "muon Pt", 100, 0, 200);
   histContainer_["tauPt"]=fs->make<TH1F>("tauPt",   "tau Pt", 100, 0, 200);
      histContainer_["tauPt30eta2p1"]=fs->make<TH1F>("tauPt30eta2p1",   "tau Pt", 100, 0, 200);


   histContainer_["jet0_50_850"]=fs->make<TH1F>("jet0_50_850",   "jet0 Pt", 100, 0, 200);
   histContainer_["jet1_50_850"]=fs->make<TH1F>("jet1_50_850",   "jet1 Pt", 100, 0, 200);
   histContainer_["jet0_50"]=fs->make<TH1F>("jet0_50",   "jet0 Pt", 100, 0, 200);
   histContainer_["jet1_50"]=fs->make<TH1F>("jet1_50",   "jet1 Pt", 100, 0, 200);
   histContainer_["jet0_850"]=fs->make<TH1F>("jet0_850",   "jet0 Pt", 100, 0, 200);
   histContainer_["jet1_850"]=fs->make<TH1F>("jet1_850",   "jet1 Pt", 100, 0, 200);
   histContainer_["jet0"]=fs->make<TH1F>("jet0",   "jet0 Pt", 100, 0, 200);
   histContainer_["jet1"]=fs->make<TH1F>("jet1",   "jet1 Pt", 100, 0, 200);


   histContainer_["jetPt"]=fs->make<TH1F>("jetPt",   "jet Pt", 100, 0, 200);
      histContainer_["jetPt30"]=fs->make<TH1F>("jetPt30",   "jet Pt", 100, 0, 200);
           histContainer_["jetPt50"]=fs->make<TH1F>("jetPt50",   "jet Pt", 100, 0, 200);
         histContainer_["jetPt30di850"]=fs->make<TH1F>("jetPt30di850",   "jet Pt", 100, 0, 200);
   histContainer_["metPt"]=fs->make<TH1F>("metPt",   "met Pt", 100, 0, 200);
      histContainer_["metPt30di850"]=fs->make<TH1F>("metPt30di850",   "met Pt", 100, 0, 200);

  
  

   // book histograms for Eta: 
   histContainer_["eleEta"]=fs->make<TH1F>("eleEta",   "electron Eta",100, -5,  5);
   histContainer_["muonEta"]=fs->make<TH1F>("muonEta",   "muon Eta",  100, -5,  5);
   histContainer_["tauEta"]=fs->make<TH1F>("tauEta",   "tau Eta",100, -5,  5);
      histContainer_["tauEta30eta2p1"]=fs->make<TH1F>("tauEta30eta2p1",   "tau Eta",  100, -5,  5);
   histContainer_["jetEta"]=fs->make<TH1F>("jetEta",   "jet Eta",  100, -5,  5);
      histContainer_["jetEta30"]=fs->make<TH1F>("jetEta30",   "jet Eta",  100, -5,  5);
         histContainer_["jetEta50"]=fs->make<TH1F>("jetEta50",   "jet Eta",  100, -5,  5);
         histContainer_["jetEta30di850"]=fs->make<TH1F>("jetEta30di850",   "jet Eta",  100, -5,  5);
   histContainer_["metEta"]=fs->make<TH1F>("metEta",   "met Eta",100, -5,  5);
      histContainer_["metEta30di850"]=fs->make<TH1F>("metEta30di850",   "met Eta",100, -5,  5);

   

   // book histograms for Phi: 
   histContainer_["elePhi"]=fs->make<TH1F>("elePhi",   "electron Phi", 100, -3.5, 3.5);
   histContainer_["muonPhi"]=fs->make<TH1F>("muonPhi",   "muon Phi",     100, -3.5, 3.5);
   histContainer_["tauPhi"]=fs->make<TH1F>("tauPhi",   "tau Phi", 100, -3.5, 3.5);
      histContainer_["tauPhi30eta2p1"]=fs->make<TH1F>("tauPhi30eta2p1",   "tau Phi",     100, -3.5, 3.5);
   histContainer_["jetPhi"]=fs->make<TH1F>("jetPhi",   "jet Phi",     100, -3.5, 3.5);
      histContainer_["jetPhi30"]=fs->make<TH1F>("jetPhi30",   "jet Phi",     100, -3.5, 3.5);
         histContainer_["jetPhi50"]=fs->make<TH1F>("jetPhi50",   "jet Phi",     100, -3.5, 3.5);
         histContainer_["jetPhi30di850"]=fs->make<TH1F>("jetPhi30di850",   "jet Phi",     100, -3.5, 3.5);
   histContainer_["metPhi"]=fs->make<TH1F>("metPhi",   "met Phi", 100, -3.5, 3.5);
      histContainer_["metPhi30di850"]=fs->make<TH1F>("metPhi30di850",   "met Phi", 100, -3.5, 3.5);


    
}

void 
Dijetmass::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Dijetmass);

