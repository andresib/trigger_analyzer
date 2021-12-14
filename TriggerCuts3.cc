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


// trigger
//#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"


using namespace std;

class TriggerCuts3 : public edm::EDAnalyzer {

public:
   explicit TriggerCuts3(const edm::ParameterSet&);
   ~TriggerCuts3();
  
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
  /// input for patTriggerEvent
 // edm::EDGetTokenT< pat::TriggerEvent > triggerEventToken_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_;

   // input tags 
   edm::InputTag muonSrc_;
   edm::InputTag elecSrc_;
   edm::InputTag tauSrc_;
   edm::InputTag jetSrc_;
   edm::InputTag metSrc_;
     /// input for patTrigger
   edm::InputTag triggerO_;
   edm::InputTag triggerR_;
   std::vector<std::string> vstrig ;

};


TriggerCuts3::TriggerCuts3(const edm::ParameterSet& iConfig):

   histContainer_(),
   muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
   elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elecSrc")),
   tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc")),
   jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc")),
   metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc")),
   triggerO_(iConfig.getUntrackedParameter<edm::InputTag>("triggerO")),
   triggerR_(iConfig.getUntrackedParameter<edm::InputTag>("triggerR")),
   vstrig(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerNames" )){

   muonCollToken = consumes<pat::MuonCollection>(muonSrc_);
   elecCollToken = consumes<pat::ElectronCollection>(elecSrc_);
   tauCollToken = consumes<pat::TauCollection>(tauSrc_);
   jetCollToken = consumes<pat::JetCollection>(jetSrc_);
   metCollToken = consumes<pat::METCollection>(metSrc_);
   triggerObjects_ = consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerO_);
   triggerResults_ =  consumes<edm::TriggerResults>(triggerR_);


/////Mapa de eventos
// 1 -> Eventos totales
// -1 -> Tau pt >20 AND eta <2.2
// -2 -> a) Tau pt >20 AND eta <2.2 AND b)  Tau pt >30 AND eta <2.1
// -3 -> Tau pt >30 AND eta <2.1
 
// -4 a) jet1_eta*jet2_eta < 0 AND delta eta >4.2  AND b) dijetmass> 750 GeV AND jet1_pt>40 jet2_pt>40
// -5 a) jet1_eta*jet2_eta < 0 AND delta eta >4.2  AND b) dijetmass> 750 GeV AND jet1_pt>40 jet2_pt>40
//AND c) dijetmass> 850 GeV AND jet1_pt>50 jet2_pt>50
// -6  jet1_eta*jet2_eta < 0 AND delta eta >4.2 AND jet1_pt>50 jet2_pt>50           dijetMass50
// -7  jet1_eta*jet2_eta < 0 AND delta eta >4.2 AND jet1_pt>50 jet2_pt>50  AND dijetmass> 850 GeV AND jet1_pt>50 jet2_pt>50      dijetMass850

// -8 Two jets pt>40 with greater Mjj opposite hemispheres Delta eta >4.2
// -9 Two jets pt>50 with greater Mjj opposite hemispheres Delta eta >4.2
// -10 Two jets pt>40 with greater Mjj>750 opposite hemispheres Delta eta >4.2      dijetMassVBFone40
// -11 Two jets pt>50  with greater Mjj>850 opposite hemispheres Delta eta >4.2     dijetMassVBFone


// -12  passed trigger 0      tauPt_trigger0_alone
// -13  passed trigger 1      tauPt_trigger1_alone
// -14  passed trigger 2      tauPt_trigger2_alone

// 2   -1   AND   -4       tauPt_cut    
// Tau pt >20 AND eta <2.2   AND   a) jet1_eta*jet2_eta < 0 AND delta eta >4.2  AND b) dijetmass> 750 GeV AND jet1_pt>40 jet2_pt>40

// 3   -2    AND   -5      tauPt_cut2
// a) Tau pt >20 AND eta <2.2 AND b)  Tau pt >30 AND eta <2.1   AND 
//    a) jet1_eta*jet2_eta < 0 AND delta eta >4.2  AND b) dijetmass> 750 GeV AND jet1_pt>40 jet2_pt>40
//AND c) dijetmass> 850 GeV AND jet1_pt>50 jet2_pt>50

// 4   -3   AND  -7        tauPt_den_corr
// Tau pt >30 AND eta <2.1    AND 
//  jet1_eta*jet2_eta < 0 AND delta eta >4.2
//  dijetmass> 850 GeV AND jet1_pt>50 jet2_pt>50

//  5  cortetau>0 &&  mjj40 > 750      tauPt_case1
//  6  cortetauden>0 &&  mjj > 850     tauPt_case2


// 7  -3    AND  -7  AND -12     cortetauden>0 && cortejetden_corrected> 0 &&  passtrig0>0         tauPt_trigger0
// 8  -3    AND  -7  AND -13     cortetauden>0 && cortejetden_corrected> 0 &&  passtrig1>0         tauPt_trigger1
// 9  -3    AND  -7  AND -14     cortetauden>0 && cortejetden_corrected> 0 &&  passtrig2>0         tauPt_trigger2

// 10    cortetau>0 &&  mjj40 > 750  &&  passtrig0>0        tauPt_case1_trigger0
// 11    cortetauden>0 &&  mjj > 850  &&  passtrig0>0       tauPt_case2_trigger0
// 12    cortetau>0 &&  mjj40 > 750  &&  passtrig1>0        tauPt_case1_trigger1
// 13    cortetauden>0 &&  mjj > 850  &&  passtrig1>0       tauPt_case2_trigger1
// 14    cortetau>0 &&  mjj40 > 750  &&  passtrig2>0        tauPt_case1_trigger2
// 15    cortetauden>0 &&  mjj > 850  &&  passtrig2>0       tauPt_case2_trigger2



}

TriggerCuts3::~TriggerCuts3(){
}

void
TriggerCuts3::analyze(const edm::Event& iEvent, 
                   const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

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

      // get trigger collection 
   edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   edm::Handle<edm::TriggerResults> triggerResultsHandle_;
   iEvent.getByToken(triggerResults_, triggerResultsHandle_);

 const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle_);

   int cortetau=0;
   int cortetau2=0;
   int cortetauden=0;
   int cortejet=0;
   int cortejet2=0;
   int cortejetden_corrected=0;
   int passtrig0=0;
   int passtrig1=0;
   int passtrig2=0;
   histContainer_["eventos"]->Fill(1); 
   float mjj=-1;
   float tmp_mjj = -1;
   float i_pt = -10;
   float j_pt = -10;
   float j_eta = -10;
   float i_eta = -10;
   float j_phi = -10;
   float i_phi = -10;

   float mjj40=-1;
   float tmp_mjj40 = -1;
   float i_pt40 = -10;
   float j_pt40 = -10;
   float j_eta40 = -10;
   float i_eta40 = -10;
   float j_phi40 = -10;
   float i_phi40 = -10;


   float tau_pt = -10;
   float tau_eta = -10;
   float tau_phi = -10;
   float tmp_taupt = -10;

   float tau_pt20 = -10;
   float tau_eta20 = -10;
   float tau_phi20 = -10;
   float tmp_taupt20 = -10;
   


   //loop over trigger results
  if(triggerResultsHandle_.isValid() == true ){
      for(unsigned int u = 0; u < names.size(); u++){
         for(unsigned int v = 0; v < vstrig.size(); v++){
            if ( names.triggerName(u).find(vstrig.at(v) ) != std::string::npos) {
				   if(triggerResultsHandle_ -> wasrun(u) ){
					   if(triggerResultsHandle_ -> accept(u)){ 
                    if (v==0) {
                        passtrig0++;
                        histContainer_["eventos"]->Fill(-12);   //-8
                          
                     }
                        else if (v==1){
                           passtrig1++;
                           histContainer_["eventos"]->Fill(-13);     //-9
                  
                        }
                        else {
                           passtrig2++;
                           histContainer_["eventos"]->Fill(-14);     //-10
  
                        }
					   }
				   }
			   }
		   }
	   }
  }

   // loop over taus
   for (auto itt = taus->cbegin(); itt != taus->cend(); ++itt) {
      histContainer_["tauPt"] ->Fill(itt->pt());
      histContainer_["tauEta"]->Fill(itt->eta());
      histContainer_["tauPhi"]->Fill(itt->phi());    



      if( itt->pt()>20 && fabs(itt->eta())<2.1 ){
         cortetau++;
         histContainer_["eventos"]->Fill(-1); /////case 1
                  tmp_taupt20=itt->pt();
                     if ( tmp_taupt20 > tau_pt20 ){ 
                        tau_pt20 = tmp_taupt20;
                        histContainer_["eventos"]->Fill(-1.4); //
                        tau_eta20=itt->eta();
                        tau_phi20=itt->phi();

                     }

  
         if( itt->pt()>30 && fabs(itt->eta())<2.1 ){
         cortetau2++;
         histContainer_["eventos"]->Fill(-2);
         }
      }

      if( itt->pt()>30 && fabs(itt->eta())<2.1 ){
         cortetauden++;
         histContainer_["eventos"]->Fill(-3);  /////case 2
         tmp_taupt=itt->pt();
                     if ( tmp_taupt > tau_pt ){ 
                        tau_pt = tmp_taupt;
                        histContainer_["eventos"]->Fill(-3.4); //
                        tau_eta=itt->eta();
                        tau_phi=itt->phi();

                     }



       
         }
   
   }

   // loop over jets
   for (auto itj = jets->cbegin(); itj != jets->cend(); ++itj) {
      histContainer_["jetPt"] ->Fill(itj->pt());
      histContainer_["jetEta"]->Fill(itj->eta());
      histContainer_["jetPhi"]->Fill(itj->phi());       
   }

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
            cortejet++;
            histContainer_["eventos"]->Fill(-4);    /////case 1
            if(  ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850 && (*jets)[0].pt()>50 && (*jets)[1].pt()>50 ) { 
               cortejet2++;
               histContainer_["eventos"]->Fill(-5); 
            }
         }
      }

     if((*jets)[0].eta()*(*jets)[1].eta()<0   &&  fabs((*jets)[0].eta() - (*jets)[1].eta()) > 4.2  ) {
        if ((*jets)[0].pt()>50 && (*jets)[1].pt()>50)
        {
          
      histContainer_["dijetMass50"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
      histContainer_["jet0Pt50"]->Fill((*jets)[0].pt());
      histContainer_["jet1Pt50"]->Fill((*jets)[1].pt()); 
      histContainer_["jet0Eta50"]->Fill((*jets)[0].eta());
      histContainer_["jet1Eta50"]->Fill((*jets)[1].eta()); 
      histContainer_["jet0Phi50"]->Fill((*jets)[0].phi());
      histContainer_["jet1Phi50"]->Fill((*jets)[1].phi()); 
      histContainer_["eventos"]->Fill(-6); /////temporal  

            if(  ((*jets)[0].p4()+(*jets)[1].p4()).mass() > 850   ) { 
            cortejetden_corrected++;
            histContainer_["eventos"]->Fill(-7);     /////case 2
            histContainer_["dijetMass850"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
            histContainer_["jet0Pt850"]->Fill((*jets)[0].pt());
            histContainer_["jet1Pt850"]->Fill((*jets)[1].pt()); 
            histContainer_["jet0Eta850"]->Fill((*jets)[0].eta());
            histContainer_["jet1Eta850"]->Fill((*jets)[1].eta()); 
            histContainer_["jet0Phi850"]->Fill((*jets)[0].phi());
            histContainer_["jet1Phi850"]->Fill((*jets)[1].phi()); 

            }
        }

      }

   }


  for (auto ite_i = jets->cbegin(); ite_i != jets->cend(); ++ite_i) {
      for (auto ite_j = ite_i; ite_j != jets->cend(); ++ite_j) {
         if (ite_j>ite_i){
            if(   ite_i->pt()>50 && ite_j->pt()>50 ) { 
               if(ite_i->eta()*ite_j->eta()<0 ) {
                  if (fabs(ite_i->eta() - ite_j->eta()) > 4.2){
                     tmp_mjj = (ite_i->p4()+ ite_j->p4()).mass();
                     if ( tmp_mjj > mjj ){ 
                        mjj = tmp_mjj;
                        histContainer_["eventos"]->Fill(-9); /////temporal  12
                        i_pt=ite_i->pt();
                        j_pt=ite_j->pt();
                        i_eta=ite_i->eta();
                        j_eta=ite_j->eta();
                        i_phi=ite_i->phi();
                        j_phi=ite_j->phi();
                     }
                  }
               }
            }
            if(   ite_i->pt()>50 && ite_j->pt()>50 ) { 
               if(ite_i->eta()*ite_j->eta()<0 ) {
                  if (fabs(ite_i->eta() - ite_j->eta()) > 4.2){
                     tmp_mjj40 = (ite_i->p4()+ ite_j->p4()).mass();
                     if ( tmp_mjj40 > mjj40 ){ 
                        mjj40 = tmp_mjj40;
                        histContainer_["eventos"]->Fill(-8); /////temporal   11
                        i_pt40=ite_i->pt();
                        j_pt40=ite_j->pt();
                        i_eta40=ite_i->eta();
                        j_eta40=ite_j->eta();
                        i_phi40=ite_i->phi();
                        j_phi40=ite_j->phi();
                     }
                  }
               }
            }
         }
      }   
   }
                     if ( mjj > 850 ){ 
                        histContainer_["dijetMassVBFone"]->Fill(mjj ); 
                        histContainer_["jet0PtVBFone"]->Fill(i_pt);
                        histContainer_["jet1PtVBFone"]->Fill(j_pt); 
                        histContainer_["jet0EtaVBFone"]->Fill(i_eta);
                        histContainer_["jet1EtaVBFone"]->Fill(j_eta); 
                        histContainer_["jet0PhiVBFone"]->Fill(i_phi);
                        histContainer_["jet1PhiVBFone"]->Fill(j_phi);
                        histContainer_["eventos"]->Fill(-11); /////temporal  ///case 2   14
                     }

                     if ( mjj40 > 850 ){ 
                        histContainer_["dijetMassVBFone40"]->Fill(mjj40 ); 
                        histContainer_["jet0PtVBFone40"]->Fill(i_pt40);
                        histContainer_["jet1PtVBFone40"]->Fill(j_pt40); 
                        histContainer_["jet0EtaVBFone40"]->Fill(i_eta40);
                        histContainer_["jet1EtaVBFone40"]->Fill(j_eta40); 
                        histContainer_["jet0PhiVBFone40"]->Fill(i_phi40);
                        histContainer_["jet1PhiVBFone40"]->Fill(j_phi40);
                        histContainer_["eventos"]->Fill(-10); /////temporal  ///case 1     13
                     }

   // loop over met
   for (auto itm = met->cbegin(); itm != met->cend(); ++itm) {
      histContainer_["metPt"] ->Fill(itm->pt());
      histContainer_["metEta"]->Fill(itm->eta());
      histContainer_["metPhi"]->Fill(itm->phi());    
   }

if (cortetau>0 && cortejet> 0) {

   histContainer_["eventos"]->Fill(2);

   histContainer_["tauPt_cut"] ->Fill(tau_pt20);
   histContainer_["tauEta_cut"]->Fill(tau_eta20);
   histContainer_["tauPhi_cut"]->Fill(tau_phi20); 

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



if (cortetauden>0 && cortejetden_corrected> 0) {

   histContainer_["eventos"]->Fill(4); 

   histContainer_["tauPt_den_corr"] ->Fill(tau_pt);
   histContainer_["tauEta_den_corr"]->Fill(tau_eta);
   histContainer_["tauPhi_den_corr"]->Fill(tau_phi); 
   histContainer_["dijetMass_den_corr"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_den_corr"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_den_corr"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_den_corr"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_den_corr"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_den_corr"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_den_corr"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_den_corr"] ->Fill((*met)[0].pt());
   histContainer_["metEta_den_corr"]->Fill((*met)[0].eta());
   histContainer_["metPhi_den_corr"]->Fill((*met)[0].phi());

}

if (cortetauden>0 && cortejetden_corrected> 0 &&  passtrig0>0 ) {

   histContainer_["eventos"]->Fill(7);    //7

   histContainer_["tauPt_trigger0"] ->Fill(tau_pt);
   histContainer_["tauEta_trigger0"]->Fill(tau_eta);
   histContainer_["tauPhi_trigger0"]->Fill(tau_phi); 
   histContainer_["dijetMass_trigger0"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_trigger0"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger0"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger0"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger0"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger0"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger0"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_trigger0"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger0"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger0"]->Fill((*met)[0].phi());

}





if (  passtrig0>0 ) {

   if ( taus->size()>0)
   {
     histContainer_["tauPt_trigger0_alone"] ->Fill((*taus)[0].pt());
     histContainer_["tauEta_trigger0_alone"]->Fill((*taus)[0].eta());
     histContainer_["tauPhi_trigger0_alone"]->Fill((*taus)[0].phi()); 
        histContainer_["tauMult_trigger0_alone"]->Fill(taus->size() );
   }
   

   if (jets->size()>1)
   {
   histContainer_["dijetMass_trigger0_alone"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());      
   histContainer_["jet0Pt_trigger0_alone"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger0_alone"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger0_alone"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger0_alone"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger0_alone"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger0_alone"]->Fill((*jets)[1].phi()); 
      histContainer_["jetMult_trigger0_alone"]->Fill(jets->size() );
   }


     if (met->size()>0)
   {
   histContainer_["metPt_trigger0_alone"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger0_alone"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger0_alone"]->Fill((*met)[0].phi()); 
      histContainer_["metMult_trigger0_alone"]->Fill(met->size() );
   }
   


}




if (cortetauden>0 && cortejetden_corrected> 0 &&  passtrig1>0 ) {

   histContainer_["eventos"]->Fill(8); 

   histContainer_["tauPt_trigger1"] ->Fill(tau_pt);
   histContainer_["tauEta_trigger1"]->Fill(tau_eta);
   histContainer_["tauPhi_trigger1"]->Fill(tau_phi);
   histContainer_["dijetMass_trigger1"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_trigger1"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger1"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger1"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger1"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger1"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger1"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_trigger1"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger1"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger1"]->Fill((*met)[0].phi());

}

if ( passtrig1>0 ) {

   if ( taus->size()>0){
   histContainer_["tauPt_trigger1_alone"] ->Fill((*taus)[0].pt());
   histContainer_["tauEta_trigger1_alone"]->Fill((*taus)[0].eta());
   histContainer_["tauPhi_trigger1_alone"]->Fill((*taus)[0].phi()); 
   histContainer_["tauMult_trigger1_alone"]->Fill(taus->size() );
   }
   if (jets->size()>1) {
   histContainer_["dijetMass_trigger1_alone"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());   
   histContainer_["jet0Pt_trigger1_alone"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger1_alone"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger1_alone"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger1_alone"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger1_alone"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger1_alone"]->Fill((*jets)[1].phi()); 
   histContainer_["jetMult_trigger1_alone"]->Fill(jets->size() );
   }
   if (met->size()>0) {
   histContainer_["metPt_trigger1_alone"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger1_alone"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger1_alone"]->Fill((*met)[0].phi());
   histContainer_["metMult_trigger1_alone"]->Fill(met->size() );
   }

}

if (cortetauden>0 && cortejetden_corrected> 0 &&  passtrig2>0 ) {

   histContainer_["eventos"]->Fill(9); 

   histContainer_["tauPt_trigger2"] ->Fill(tau_pt);
   histContainer_["tauEta_trigger2"]->Fill(tau_eta);
   histContainer_["tauPhi_trigger2"]->Fill(tau_phi);
   histContainer_["dijetMass_trigger2"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
    
   histContainer_["jet0Pt_trigger2"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger2"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger2"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger2"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger2"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger2"]->Fill((*jets)[1].phi()); 

   histContainer_["metPt_trigger2"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger2"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger2"]->Fill((*met)[0].phi());
}
if ( passtrig2>0 ) {

   if ( taus->size()>0){
   histContainer_["tauPt_trigger2_alone"] ->Fill((*taus)[0].pt());
   histContainer_["tauEta_trigger2_alone"]->Fill((*taus)[0].eta());
   histContainer_["tauPhi_trigger2_alone"]->Fill((*taus)[0].phi()); 
   histContainer_["tauMult_trigger2_alone"]->Fill(taus->size() );
   }

   if (jets->size()>1) {
   histContainer_["dijetMass_trigger2_alone"]->Fill(((*jets)[0].p4()+(*jets)[1].p4()).mass());  
   histContainer_["jet0Pt_trigger2_alone"]->Fill((*jets)[0].pt());
   histContainer_["jet1Pt_trigger2_alone"]->Fill((*jets)[1].pt());  
   histContainer_["jet0Eta_trigger2_alone"]->Fill((*jets)[0].eta());
   histContainer_["jet1Eta_trigger2_alone"]->Fill((*jets)[1].eta()); 
   histContainer_["jet0Phi_trigger2_alone"]->Fill((*jets)[0].phi());
   histContainer_["jet1Phi_trigger2_alone"]->Fill((*jets)[1].phi()); 
   histContainer_["jetMult_trigger2_alone"]->Fill(jets->size() );
   }

  if (met->size()>0){
   histContainer_["metPt_trigger2_alone"] ->Fill((*met)[0].pt());
   histContainer_["metEta_trigger2_alone"]->Fill((*met)[0].eta());
   histContainer_["metPhi_trigger2_alone"]->Fill((*met)[0].phi());
   histContainer_["metMult_trigger2_alone"]->Fill(met->size() );
  }
}



if (cortetau>0 &&  mjj40 > 850  &&  passtrig0>0 ) {

   histContainer_["eventos"]->Fill(10);  //-11

   histContainer_["tauPt_case1_trigger0"]->Fill(tau_pt20);
   histContainer_["tauEta_case1_trigger0"]->Fill(tau_eta20);
   histContainer_["tauPhi_case1_trigger0"]->Fill(tau_phi20); 

                        histContainer_["dijetMass_case1_trigger0"]->Fill(mjj40 ); 
                        histContainer_["jet0Pt_case1_trigger0"]->Fill(i_pt40);
                        histContainer_["jet1Pt_case1_trigger0"]->Fill(j_pt40); 
                        histContainer_["jet0Eta_case1_trigger0"]->Fill(i_eta40);
                        histContainer_["jet1Eta_case1_trigger0"]->Fill(j_eta40); 
                        histContainer_["jet0Phi_case1_trigger0"]->Fill(i_phi40);
                        histContainer_["jet1Phi_case1_trigger0"]->Fill(j_phi40);

   histContainer_["metPt_case1_trigger0"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case1_trigger0"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case1_trigger0"]->Fill((*met)[0].phi());

}
if (cortetauden>0 &&  mjj > 850  &&  passtrig0>0 ) {

   histContainer_["eventos"]->Fill(11);   //-12

   histContainer_["tauPt_case2_trigger0"] ->Fill(tau_pt);
   histContainer_["tauEta_case2_trigger0"]->Fill(tau_eta);
   histContainer_["tauPhi_case2_trigger0"]->Fill(tau_phi);

                        histContainer_["dijetMass_case2_trigger0"]->Fill(mjj ); 
                        histContainer_["jet0Pt_case2_trigger0"]->Fill(i_pt);
                        histContainer_["jet1Pt_case2_trigger0"]->Fill(j_pt); 
                        histContainer_["jet0Eta_case2_trigger0"]->Fill(i_eta);
                        histContainer_["jet1Eta_case2_trigger0"]->Fill(j_eta); 
                        histContainer_["jet0Phi_case2_trigger0"]->Fill(i_phi);
                        histContainer_["jet1Phi_case2_trigger0"]->Fill(j_phi);


   histContainer_["metPt_case2_trigger0"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case2_trigger0"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case2_trigger0"]->Fill((*met)[0].phi());

}
if (cortetau>0 &&  mjj40 > 850  &&  passtrig1>0 ) {

   histContainer_["eventos"]->Fill(12);  //-13


   histContainer_["tauPt_case1_trigger1"]->Fill(tau_pt20);
   histContainer_["tauEta_case1_trigger1"]->Fill(tau_eta20);
   histContainer_["tauPhi_case1_trigger1"]->Fill(tau_phi20); 

                        histContainer_["dijetMass_case1_trigger1"]->Fill(mjj40 ); 
                        histContainer_["jet0Pt_case1_trigger1"]->Fill(i_pt40);
                        histContainer_["jet1Pt_case1_trigger1"]->Fill(j_pt40); 
                        histContainer_["jet0Eta_case1_trigger1"]->Fill(i_eta40);
                        histContainer_["jet1Eta_case1_trigger1"]->Fill(j_eta40); 
                        histContainer_["jet0Phi_case1_trigger1"]->Fill(i_phi40);
                        histContainer_["jet1Phi_case1_trigger1"]->Fill(j_phi40);

   histContainer_["metPt_case1_trigger1"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case1_trigger1"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case1_trigger1"]->Fill((*met)[0].phi());

}
if (cortetauden>0 &&  mjj > 850  &&  passtrig1>0 ) {

   histContainer_["eventos"]->Fill(13);  //-14

   histContainer_["tauPt_case2_trigger1"] ->Fill(tau_pt);
   histContainer_["tauEta_case2_trigger1"]->Fill(tau_eta);
   histContainer_["tauPhi_case2_trigger1"]->Fill(tau_phi);

   

                        histContainer_["dijetMass_case2_trigger1"]->Fill(mjj ); 
                        histContainer_["jet0Pt_case2_trigger1"]->Fill(i_pt);
                        histContainer_["jet1Pt_case2_trigger1"]->Fill(j_pt); 
                        histContainer_["jet0Eta_case2_trigger1"]->Fill(i_eta);
                        histContainer_["jet1Eta_case2_trigger1"]->Fill(j_eta); 
                        histContainer_["jet0Phi_case2_trigger1"]->Fill(i_phi);
                        histContainer_["jet1Phi_case2_trigger1"]->Fill(j_phi);

   histContainer_["metPt_case2_trigger1"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case2_trigger1"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case2_trigger1"]->Fill((*met)[0].phi());

}

if (cortetau>0 &&  mjj40 > 850  &&  passtrig2>0 ) {

   histContainer_["eventos"]->Fill(14);  //15

   histContainer_["tauPt_case1_trigger2"] ->Fill(tau_pt20);
   histContainer_["tauEta_case1_trigger2"]->Fill(tau_eta20);
   histContainer_["tauPhi_case1_trigger2"]->Fill(tau_phi20);  

                        histContainer_["dijetMass_case1_trigger2"]->Fill(mjj40 ); 
                        histContainer_["jet0Pt_case1_trigger2"]->Fill(i_pt40);
                        histContainer_["jet1Pt_case1_trigger2"]->Fill(j_pt40); 
                        histContainer_["jet0Eta_case1_trigger2"]->Fill(i_eta40);
                        histContainer_["jet1Eta_case1_trigger2"]->Fill(j_eta40); 
                        histContainer_["jet0Phi_case1_trigger2"]->Fill(i_phi40);
                        histContainer_["jet1Phi_case1_trigger2"]->Fill(j_phi40);

   histContainer_["metPt_case1_trigger2"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case1_trigger2"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case1_trigger2"]->Fill((*met)[0].phi());

}
if (cortetauden>0 &&  mjj > 850  &&  passtrig2>0 ) {

   histContainer_["eventos"]->Fill(15);  //-16

   histContainer_["tauPt_case2_trigger2"] ->Fill(tau_pt);
   histContainer_["tauEta_case2_trigger2"]->Fill(tau_eta);
   histContainer_["tauPhi_case2_trigger2"]->Fill(tau_phi);

                        histContainer_["dijetMass_case2_trigger2"]->Fill(mjj ); 
                        histContainer_["jet0Pt_case2_trigger2"]->Fill(i_pt);
                        histContainer_["jet1Pt_case2_trigger2"]->Fill(j_pt); 
                        histContainer_["jet0Eta_case2_trigger2"]->Fill(i_eta);
                        histContainer_["jet1Eta_case2_trigger2"]->Fill(j_eta); 
                        histContainer_["jet0Phi_case2_trigger2"]->Fill(i_phi);
                        histContainer_["jet1Phi_case2_trigger2"]->Fill(j_phi);

   histContainer_["metPt_case2_trigger2"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case2_trigger2"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case2_trigger2"]->Fill((*met)[0].phi());

}

if (cortetau>0 &&  mjj40 > 850  ) {

   histContainer_["eventos"]->Fill(5);   //-17

   histContainer_["tauPt_case1"] ->Fill(tau_pt20);
   histContainer_["tauEta_case1"]->Fill(tau_eta20);
   histContainer_["tauPhi_case1"]->Fill(tau_phi20);  

                        histContainer_["dijetMass_case1"]->Fill(mjj40 ); 
                        histContainer_["jet0Pt_case1"]->Fill(i_pt40);
                        histContainer_["jet1Pt_case1"]->Fill(j_pt40); 
                        histContainer_["jet0Eta_case1"]->Fill(i_eta40);
                        histContainer_["jet1Eta_case1"]->Fill(j_eta40); 
                        histContainer_["jet0Phi_case1"]->Fill(i_phi40);
                        histContainer_["jet1Phi_case1"]->Fill(j_phi40);

   histContainer_["metPt_case1"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case1"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case1"]->Fill((*met)[0].phi());

}
if (cortetauden>0 &&  mjj > 850   ) {

   histContainer_["eventos"]->Fill(6);  //-18

   histContainer_["tauPt_case2"] ->Fill(tau_pt);
   histContainer_["tauEta_case2"]->Fill(tau_eta);
   histContainer_["tauPhi_case2"]->Fill(tau_phi);

                        histContainer_["dijetMass_case2"]->Fill(mjj ); 
                        histContainer_["jet0Pt_case2"]->Fill(i_pt);
                        histContainer_["jet1Pt_case2"]->Fill(j_pt); 
                        histContainer_["jet0Eta_case2"]->Fill(i_eta);
                        histContainer_["jet1Eta_case2"]->Fill(j_eta); 
                        histContainer_["jet0Phi_case2"]->Fill(i_phi);
                        histContainer_["jet1Phi_case2"]->Fill(j_phi);

   histContainer_["metPt_case2"] ->Fill((*met)[0].pt());
   histContainer_["metEta_case2"]->Fill((*met)[0].eta());
   histContainer_["metPhi_case2"]->Fill((*met)[0].phi());

}





   // Multiplicity
 
   histContainer_["tauMult"]->Fill(taus->size() );
   histContainer_["jetMult"]->Fill(jets->size() );
   histContainer_["metMult"]->Fill(met->size() );

}



void 
TriggerCuts3::beginJob()
{
   // register to the TFileService
   edm::Service<TFileService> fs;

//  number of events

   histContainer_["eventos"]=fs->make<TH1F>("eventos",   "numero eventos", 80, -20,  20);



   // book histograms for Tau


   histContainer_["tauPt"]=fs->make<TH1F>("tauPt",   "tau Pt", 300, 0, 600);
   histContainer_["tauPt_cut"]=fs->make<TH1F>("tauPt_cut",   "tau Pt", 300, 0, 600);
   histContainer_["tauPt_cut2"]=fs->make<TH1F>("tauPt_cut2",   "tau Pt", 300, 0, 600);
 
      histContainer_["tauPt_den_corr"]=fs->make<TH1F>("tauPt_den_corr",   "tau Pt", 300, 0, 600);
   histContainer_["tauEta"]=fs->make<TH1F>("tauEta",   "tau Eta",100, -5,  5);
   histContainer_["tauEta_cut"]=fs->make<TH1F>("tauEta_cut",   "tau Eta",100, -5,  5);
   histContainer_["tauEta_cut2"]=fs->make<TH1F>("tauEta_cut2",   "tau Eta",100, -5,  5);

      histContainer_["tauEta_den_corr"]=fs->make<TH1F>("tauEta_den_corr",   "tau Eta",100, -5,  5);
   histContainer_["tauPhi"]=fs->make<TH1F>("tauPhi",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauPhi_cut"]=fs->make<TH1F>("tauPhi_cut",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauPhi_cut2"]=fs->make<TH1F>("tauPhi_cut2",   "tau Phi", 100, -3.5, 3.5);

     histContainer_["tauPhi_den_corr"]=fs->make<TH1F>("tauPhi_den_corr",   "tau Phi", 100, -3.5, 3.5);
   histContainer_["tauMult"]=fs->make<TH1F>("tauMult",   "tau multiplicity",     100, 0,  50);



   // book histograms for MET

   histContainer_["metPt"]=fs->make<TH1F>("metPt",   "met Pt", 300, 0, 600);
   histContainer_["metPt_cut"]=fs->make<TH1F>("metPt_cut",   "met Pt", 300, 0, 600);
   histContainer_["metPt_cut2"]=fs->make<TH1F>("metPt_cut2",   "met Pt", 300, 0, 600);
 
      histContainer_["metPt_den_corr"]=fs->make<TH1F>("metPt_den_corr",   "met Pt", 300, 0, 600);
   histContainer_["metEta"]=fs->make<TH1F>("metEta",   "met Eta",100, -5,  5);
   histContainer_["metEta_cut"]=fs->make<TH1F>("metEta_cut",   "met Eta",100, -5,  5);
   histContainer_["metEta_cut2"]=fs->make<TH1F>("metEta_cut2",   "met Eta",100, -5,  5);
 
     histContainer_["metEta_den_corr"]=fs->make<TH1F>("metEta_den_corr",   "met Eta",100, -5,  5);
   histContainer_["metPhi"]=fs->make<TH1F>("metPhi",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metPhi_cut"]=fs->make<TH1F>("metPhi_cut",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metPhi_cut2"]=fs->make<TH1F>("metPhi_cut2",   "met Phi", 100, -3.5, 3.5);
 
      histContainer_["metPhi_den_corr"]=fs->make<TH1F>("metPhi_den_corr",   "met Phi", 100, -3.5, 3.5);
   histContainer_["metMult"]=fs->make<TH1F>("metMult",   "met multiplicity",     100, 0,  50);



   // book histograms for Jets
  

   histContainer_["jetPt"]=fs->make<TH1F>("jetPt",   "jet Pt", 300, 0, 600);
   histContainer_["jet0Pt"]=fs->make<TH1F>("jet0pt",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet0Pt_cut"]=fs->make<TH1F>("jet0Pt_cut",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet0Pt_cut2"]=fs->make<TH1F>("jet0Pt_cut2",   "jet0 Pt", 300, 0, 600);

      histContainer_["jet0Pt_den_corr"]=fs->make<TH1F>("jet0Pt_den_corr",   "jet0 Pt", 300, 0, 600);
   histContainer_["jet1Pt"]=fs->make<TH1F>("jet1pt",   "jet1 Pt", 300, 0, 600);
   histContainer_["jet1Pt_cut"]=fs->make<TH1F>("jet1Pt_cut",   "jet1 Pt", 300, 0, 600);
   histContainer_["jet1Pt_cut2"]=fs->make<TH1F>("jet1Pt_cut2",   "jet1 Pt", 300, 0, 600);

      histContainer_["jet1Pt_den_corr"]=fs->make<TH1F>("jet1Pt_den_corr",   "jet1 Pt", 300, 0, 600);

   histContainer_["jetEta"]=fs->make<TH1F>("jetEta",   "jet Eta",  100, -5,  5);
   histContainer_["jet0Eta"]=fs->make<TH1F>("jet0Eta",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet0Eta_cut"]=fs->make<TH1F>("jet0Eta_cut",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet0Eta_cut2"]=fs->make<TH1F>("jet0Eta_cut2",   "jet0 Eta",  100, -5,  5);

      histContainer_["jet0Eta_den_corr"]=fs->make<TH1F>("jet0Eta_den_corr",   "jet0 Eta",  100, -5,  5);
   histContainer_["jet1Eta"]=fs->make<TH1F>("jet1Eta",   "jet1 Eta",  100, -5,  5);
   histContainer_["jet1Eta_cut"]=fs->make<TH1F>("jet1Eta_cut",   "jet1 Eta",  100, -5,  5);
   histContainer_["jet1Eta_cut2"]=fs->make<TH1F>("jet1Eta_cut2",   "jet1 Eta",  100, -5,  5);

      histContainer_["jet1Eta_den_corr"]=fs->make<TH1F>("jet1Eta_den_corr",   "jet1 Eta",  100, -5,  5);

   histContainer_["jetPhi"]=fs->make<TH1F>("jetPhi",   "jet Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi"]=fs->make<TH1F>("jet0Phi",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi_cut"]=fs->make<TH1F>("jet0Phi_cut",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet0Phi_cut2"]=fs->make<TH1F>("jet0Phi_cut2",   "jet0 Phi",     100, -3.5, 3.5);
 
      histContainer_["jet0Phi_den_corr"]=fs->make<TH1F>("jet0Phi_den_corr",   "jet0 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi"]=fs->make<TH1F>("jet1Phi",   "jet1 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi_cut"]=fs->make<TH1F>("jet1Phi_cut",   "jet1 Phi",     100, -3.5, 3.5);
   histContainer_["jet1Phi_cut2"]=fs->make<TH1F>("jet1Phi_cut2",   "jet1 Phi",     100, -3.5, 3.5);

      histContainer_["jet1Phi_den_corr"]=fs->make<TH1F>("jet1Phi_den_corr",   "jet1 Phi",     100, -3.5, 3.5);

   histContainer_["jetMult"]=fs->make<TH1F>("jetMult",   "jet multiplicity",     100, 0,  50);

   histContainer_["dijetMass"]=fs->make<TH1F>("dijetMass", "mass",    500,   0., 2500.);
   histContainer_["dijetMass_cut"]=fs->make<TH1F>("dijetMass_cut", "mass",    500,   0., 2500.);
   histContainer_["dijetMass_cut2"]=fs->make<TH1F>("dijetMass_cut2", "mass",    500,   0., 2500.);

     histContainer_["dijetMass_den_corr"]=fs->make<TH1F>("dijetMass_den_corr", "mass",    500,   0., 2500.);



         histContainer_["tauPt_trigger0"]=fs->make<TH1F>("tauPt_trigger0",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger0"]=fs->make<TH1F>("tauEta_trigger0",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger0"]=fs->make<TH1F>("tauPhi_trigger0",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger0"]=fs->make<TH1F>("metPt_trigger0",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger0"]=fs->make<TH1F>("metEta_trigger0",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger0"]=fs->make<TH1F>("metPhi_trigger0",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger0"]=fs->make<TH1F>("jet0Pt_trigger0",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger0"]=fs->make<TH1F>("jet1Pt_trigger0",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger0"]=fs->make<TH1F>("jet0Eta_trigger0",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger0"]=fs->make<TH1F>("jet1Eta_trigger0",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger0"]=fs->make<TH1F>("jet0Phi_trigger0",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger0"]=fs->make<TH1F>("jet1Phi_trigger0",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger0"]=fs->make<TH1F>("dijetMass_trigger0", "mass",    500,   0., 2500.);
         
         histContainer_["tauPt_trigger1"]=fs->make<TH1F>("tauPt_trigger1",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger1"]=fs->make<TH1F>("tauEta_trigger1",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger1"]=fs->make<TH1F>("tauPhi_trigger1",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger1"]=fs->make<TH1F>("metPt_trigger1",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger1"]=fs->make<TH1F>("metEta_trigger1",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger1"]=fs->make<TH1F>("metPhi_trigger1",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger1"]=fs->make<TH1F>("jet0Pt_trigger1",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger1"]=fs->make<TH1F>("jet1Pt_trigger1",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger1"]=fs->make<TH1F>("jet0Eta_trigger1",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger1"]=fs->make<TH1F>("jet1Eta_trigger1",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger1"]=fs->make<TH1F>("jet0Phi_trigger1",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger1"]=fs->make<TH1F>("jet1Phi_trigger1",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger1"]=fs->make<TH1F>("dijetMass_trigger1", "mass",    500,   0., 2500.);

         histContainer_["tauPt_trigger2"]=fs->make<TH1F>("tauPt_trigger2",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger2"]=fs->make<TH1F>("tauEta_trigger2",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger2"]=fs->make<TH1F>("tauPhi_trigger2",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger2"]=fs->make<TH1F>("metPt_trigger2",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger2"]=fs->make<TH1F>("metEta_trigger2",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger2"]=fs->make<TH1F>("metPhi_trigger2",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger2"]=fs->make<TH1F>("jet0Pt_trigger2",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger2"]=fs->make<TH1F>("jet1Pt_trigger2",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger2"]=fs->make<TH1F>("jet0Eta_trigger2",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger2"]=fs->make<TH1F>("jet1Eta_trigger2",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger2"]=fs->make<TH1F>("jet0Phi_trigger2",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger2"]=fs->make<TH1F>("jet1Phi_trigger2",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger2"]=fs->make<TH1F>("dijetMass_trigger2", "mass",    500,   0., 2500.);


         histContainer_["tauPt_trigger0_alone"]=fs->make<TH1F>("tauPt_trigger0_alone",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger0_alone"]=fs->make<TH1F>("tauEta_trigger0_alone",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger0_alone"]=fs->make<TH1F>("tauPhi_trigger0_alone",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger0_alone"]=fs->make<TH1F>("metPt_trigger0_alone",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger0_alone"]=fs->make<TH1F>("metEta_trigger0_alone",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger0_alone"]=fs->make<TH1F>("metPhi_trigger0_alone",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger0_alone"]=fs->make<TH1F>("jet0Pt_trigger0_alone",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger0_alone"]=fs->make<TH1F>("jet1Pt_trigger0_alone",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger0_alone"]=fs->make<TH1F>("jet0Eta_trigger0_alone",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger0_alone"]=fs->make<TH1F>("jet1Eta_trigger0_alone",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger0_alone"]=fs->make<TH1F>("jet0Phi_trigger0_alone",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger0_alone"]=fs->make<TH1F>("jet1Phi_trigger0_alone",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger0_alone"]=fs->make<TH1F>("dijetMass_trigger0_alone", "mass",    500,   0., 2500.);


         histContainer_["tauMult_trigger0_alone"]=fs->make<TH1F>("tauMult_trigger0_alone",   "tau multiplicity",     100, 0,  50);
         histContainer_["metMult_trigger0_alone"]=fs->make<TH1F>("metMult_trigger0_alone",   "met multiplicity",     100, 0,  50);
         histContainer_["jetMult_trigger0_alone"]=fs->make<TH1F>("jetMult_trigger0_alone",   "jet multiplicity",     100, 0,  50);



         histContainer_["tauPt_trigger1_alone"]=fs->make<TH1F>("tauPt_trigger1_alone",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger1_alone"]=fs->make<TH1F>("tauEta_trigger1_alone",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger1_alone"]=fs->make<TH1F>("tauPhi_trigger1_alone",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger1_alone"]=fs->make<TH1F>("metPt_trigger1_alone",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger1_alone"]=fs->make<TH1F>("metEta_trigger1_alone",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger1_alone"]=fs->make<TH1F>("metPhi_trigger1_alone",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger1_alone"]=fs->make<TH1F>("jet0Pt_trigger1_alone",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger1_alone"]=fs->make<TH1F>("jet1Pt_trigger1_alone",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger1_alone"]=fs->make<TH1F>("jet0Eta_trigger1_alone",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger1_alone"]=fs->make<TH1F>("jet1Eta_trigger1_alone",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger1_alone"]=fs->make<TH1F>("jet0Phi_trigger1_alone",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger1_alone"]=fs->make<TH1F>("jet1Phi_trigger1_alone",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger1_alone"]=fs->make<TH1F>("dijetMass_trigger1_alone", "mass",    500,   0., 2500.);

         histContainer_["tauMult_trigger1_alone"]=fs->make<TH1F>("tauMult_trigger1_alone",   "tau multiplicity",     100, 0,  50);
         histContainer_["metMult_trigger1_alone"]=fs->make<TH1F>("metMult_trigger1_alone",   "met multiplicity",     100, 0,  50);
         histContainer_["jetMult_trigger1_alone"]=fs->make<TH1F>("jetMult_trigger1_alone",   "jet multiplicity",     100, 0,  50);

         histContainer_["tauPt_trigger2_alone"]=fs->make<TH1F>("tauPt_trigger2_alone",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_trigger2_alone"]=fs->make<TH1F>("tauEta_trigger2_alone",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_trigger2_alone"]=fs->make<TH1F>("tauPhi_trigger2_alone",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_trigger2_alone"]=fs->make<TH1F>("metPt_trigger2_alone",   "met Pt", 300, 0, 600);
         histContainer_["metEta_trigger2_alone"]=fs->make<TH1F>("metEta_trigger2_alone",   "met Eta",100, -5,  5);
         histContainer_["metPhi_trigger2_alone"]=fs->make<TH1F>("metPhi_trigger2_alone",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_trigger2_alone"]=fs->make<TH1F>("jet0Pt_trigger2_alone",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_trigger2_alone"]=fs->make<TH1F>("jet1Pt_trigger2_alone",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_trigger2_alone"]=fs->make<TH1F>("jet0Eta_trigger2_alone",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_trigger2_alone"]=fs->make<TH1F>("jet1Eta_trigger2_alone",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_trigger2_alone"]=fs->make<TH1F>("jet0Phi_trigger2_alone",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_trigger2_alone"]=fs->make<TH1F>("jet1Phi_trigger2_alone",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_trigger2_alone"]=fs->make<TH1F>("dijetMass_trigger2_alone", "mass",    500,   0., 2500.);

         histContainer_["tauMult_trigger2_alone"]=fs->make<TH1F>("tauMult_trigger2_alone",   "tau multiplicity",     100, 0,  50);
         histContainer_["metMult_trigger2_alone"]=fs->make<TH1F>("metMult_trigger2_alone",   "met multiplicity",     100, 0,  50);
         histContainer_["jetMult_trigger2_alone"]=fs->make<TH1F>("jetMult_trigger2_alone",   "jet multiplicity",     100, 0,  50);

      

       //  histContainer_["dijetMassVBF"]=fs->make<TH1F>("dijetMassVBF", "mass",    500,   0., 2500.);
       //  histContainer_["jet0PtVBF"]=fs->make<TH1F>("jet0PtVBF",   "jet0 Pt", 300, 0, 600);
       //  histContainer_["jet1PtVBF"]=fs->make<TH1F>("jet1PtVBF",   "jet1 Pt", 300, 0, 600);
       //  histContainer_["jet0EtaVBF"]=fs->make<TH1F>("jet0EtaVBF",   "jet0 Eta",  100, -5,  5);
       //  histContainer_["jet1EtaVBF"]=fs->make<TH1F>("jet1EtaVBF",   "jet1 Eta",  100, -5,  5);
       //  histContainer_["jet0PhiVBF"]=fs->make<TH1F>("jet0PhiVBF",   "jet0 Phi",     100, -3.5, 3.5);
       //  histContainer_["jet1PhiVBF"]=fs->make<TH1F>("jet1PhiVBF",   "jet1 Phi",     100, -3.5, 3.5);

         histContainer_["dijetMassVBFone"]=fs->make<TH1F>("dijetMassVBFone", "mass",    500,   0., 2500.);
         histContainer_["jet0PtVBFone"]=fs->make<TH1F>("jet0PtVBFone",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1PtVBFone"]=fs->make<TH1F>("jet1PtVBFone",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0EtaVBFone"]=fs->make<TH1F>("jet0EtaVBFone",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1EtaVBFone"]=fs->make<TH1F>("jet1EtaVBFone",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0PhiVBFone"]=fs->make<TH1F>("jet0PhiVBFone",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1PhiVBFone"]=fs->make<TH1F>("jet1PhiVBFone",   "jet1 Phi",     100, -3.5, 3.5);

         histContainer_["dijetMassVBFone40"]=fs->make<TH1F>("dijetMassVBFone40", "mass",    500,   0., 2500.);
         histContainer_["jet0PtVBFone40"]=fs->make<TH1F>("jet0PtVBFone40",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1PtVBFone40"]=fs->make<TH1F>("jet1PtVBFone40",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0EtaVBFone40"]=fs->make<TH1F>("jet0EtaVBFone40",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1EtaVBFone40"]=fs->make<TH1F>("jet1EtaVBFone40",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0PhiVBFone40"]=fs->make<TH1F>("jet0PhiVBFone40",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1PhiVBFone40"]=fs->make<TH1F>("jet1PhiVBFone40",   "jet1 Phi",     100, -3.5, 3.5);



         histContainer_["dijetMass50"]=fs->make<TH1F>("dijetMass50", "mass",    500,   0., 2500.);
         histContainer_["jet0Pt50"]=fs->make<TH1F>("jet0Pt50",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt50"]=fs->make<TH1F>("jet1Pt50",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta50"]=fs->make<TH1F>("jet0Eta50",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta50"]=fs->make<TH1F>("jet1Eta50",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi50"]=fs->make<TH1F>("jet0Phi50",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi50"]=fs->make<TH1F>("jet1Phi50",   "jet1 Phi",     100, -3.5, 3.5);


         histContainer_["dijetMass850"]=fs->make<TH1F>("dijetMass850", "mass",    500,   0., 2500.);
         histContainer_["jet0Pt850"]=fs->make<TH1F>("jet0Pt850",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt850"]=fs->make<TH1F>("jet1Pt850",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta850"]=fs->make<TH1F>("jet0Eta850",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta850"]=fs->make<TH1F>("jet1Eta850",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi850"]=fs->make<TH1F>("jet0Phi850",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi850"]=fs->make<TH1F>("jet1Phi850",   "jet1 Phi",     100, -3.5, 3.5);


         histContainer_["tauPt_case1_trigger0"]=fs->make<TH1F>("tauPt_case1_trigger0",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case1_trigger0"]=fs->make<TH1F>("tauEta_case1_trigger0",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case1_trigger0"]=fs->make<TH1F>("tauPhi_case1_trigger0",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case1_trigger0"]=fs->make<TH1F>("metPt_case1_trigger0",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case1_trigger0"]=fs->make<TH1F>("metEta_case1_trigger0",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case1_trigger0"]=fs->make<TH1F>("metPhi_case1_trigger0",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case1_trigger0"]=fs->make<TH1F>("jet0Pt_case1_trigger0",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case1_trigger0"]=fs->make<TH1F>("jet1Pt_case1_trigger0",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case1_trigger0"]=fs->make<TH1F>("jet0Eta_case1_trigger0",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case1_trigger0"]=fs->make<TH1F>("jet1Eta_case1_trigger0",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case1_trigger0"]=fs->make<TH1F>("jet0Phi_case1_trigger0",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case1_trigger0"]=fs->make<TH1F>("jet1Phi_case1_trigger0",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case1_trigger0"]=fs->make<TH1F>("dijetMass_case1_trigger0", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case2_trigger0"]=fs->make<TH1F>("tauPt_case2_trigger0",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case2_trigger0"]=fs->make<TH1F>("tauEta_case2_trigger0",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case2_trigger0"]=fs->make<TH1F>("tauPhi_case2_trigger0",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case2_trigger0"]=fs->make<TH1F>("metPt_case2_trigger0",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case2_trigger0"]=fs->make<TH1F>("metEta_case2_trigger0",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case2_trigger0"]=fs->make<TH1F>("metPhi_case2_trigger0",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case2_trigger0"]=fs->make<TH1F>("jet0Pt_case2_trigger0",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case2_trigger0"]=fs->make<TH1F>("jet1Pt_case2_trigger0",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case2_trigger0"]=fs->make<TH1F>("jet0Eta_case2_trigger0",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case2_trigger0"]=fs->make<TH1F>("jet1Eta_case2_trigger0",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case2_trigger0"]=fs->make<TH1F>("jet0Phi_case2_trigger0",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case2_trigger0"]=fs->make<TH1F>("jet1Phi_case2_trigger0",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case2_trigger0"]=fs->make<TH1F>("dijetMass_case2_trigger0", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case1_trigger1"]=fs->make<TH1F>("tauPt_case1_trigger1",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case1_trigger1"]=fs->make<TH1F>("tauEta_case1_trigger1",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case1_trigger1"]=fs->make<TH1F>("tauPhi_case1_trigger1",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case1_trigger1"]=fs->make<TH1F>("metPt_case1_trigger1",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case1_trigger1"]=fs->make<TH1F>("metEta_case1_trigger1",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case1_trigger1"]=fs->make<TH1F>("metPhi_case1_trigger1",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case1_trigger1"]=fs->make<TH1F>("jet0Pt_case1_trigger1",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case1_trigger1"]=fs->make<TH1F>("jet1Pt_case1_trigger1",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case1_trigger1"]=fs->make<TH1F>("jet0Eta_case1_trigger1",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case1_trigger1"]=fs->make<TH1F>("jet1Eta_case1_trigger1",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case1_trigger1"]=fs->make<TH1F>("jet0Phi_case1_trigger1",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case1_trigger1"]=fs->make<TH1F>("jet1Phi_case1_trigger1",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case1_trigger1"]=fs->make<TH1F>("dijetMass_case1_trigger1", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case2_trigger1"]=fs->make<TH1F>("tauPt_case2_trigger1",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case2_trigger1"]=fs->make<TH1F>("tauEta_case2_trigger1",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case2_trigger1"]=fs->make<TH1F>("tauPhi_case2_trigger1",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case2_trigger1"]=fs->make<TH1F>("metPt_case2_trigger1",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case2_trigger1"]=fs->make<TH1F>("metEta_case2_trigger1",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case2_trigger1"]=fs->make<TH1F>("metPhi_case2_trigger1",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case2_trigger1"]=fs->make<TH1F>("jet0Pt_case2_trigger1",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case2_trigger1"]=fs->make<TH1F>("jet1Pt_case2_trigger1",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case2_trigger1"]=fs->make<TH1F>("jet0Eta_case2_trigger1",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case2_trigger1"]=fs->make<TH1F>("jet1Eta_case2_trigger1",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case2_trigger1"]=fs->make<TH1F>("jet0Phi_case2_trigger1",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case2_trigger1"]=fs->make<TH1F>("jet1Phi_case2_trigger1",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case2_trigger1"]=fs->make<TH1F>("dijetMass_case2_trigger1", "mass",    500,   0., 2500.);



         histContainer_["tauPt_case1_trigger2"]=fs->make<TH1F>("tauPt_case1_trigger2",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case1_trigger2"]=fs->make<TH1F>("tauEta_case1_trigger2",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case1_trigger2"]=fs->make<TH1F>("tauPhi_case1_trigger2",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case1_trigger2"]=fs->make<TH1F>("metPt_case1_trigger2",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case1_trigger2"]=fs->make<TH1F>("metEta_case1_trigger2",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case1_trigger2"]=fs->make<TH1F>("metPhi_case1_trigger2",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case1_trigger2"]=fs->make<TH1F>("jet0Pt_case1_trigger2",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case1_trigger2"]=fs->make<TH1F>("jet1Pt_case1_trigger2",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case1_trigger2"]=fs->make<TH1F>("jet0Eta_case1_trigger2",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case1_trigger2"]=fs->make<TH1F>("jet1Eta_case1_trigger2",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case1_trigger2"]=fs->make<TH1F>("jet0Phi_case1_trigger2",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case1_trigger2"]=fs->make<TH1F>("jet1Phi_case1_trigger2",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case1_trigger2"]=fs->make<TH1F>("dijetMass_case1_trigger2", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case2_trigger2"]=fs->make<TH1F>("tauPt_case2_trigger2",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case2_trigger2"]=fs->make<TH1F>("tauEta_case2_trigger2",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case2_trigger2"]=fs->make<TH1F>("tauPhi_case2_trigger2",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case2_trigger2"]=fs->make<TH1F>("metPt_case2_trigger2",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case2_trigger2"]=fs->make<TH1F>("metEta_case2_trigger2",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case2_trigger2"]=fs->make<TH1F>("metPhi_case2_trigger2",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case2_trigger2"]=fs->make<TH1F>("jet0Pt_case2_trigger2",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case2_trigger2"]=fs->make<TH1F>("jet1Pt_case2_trigger2",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case2_trigger2"]=fs->make<TH1F>("jet0Eta_case2_trigger2",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case2_trigger2"]=fs->make<TH1F>("jet1Eta_case2_trigger2",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case2_trigger2"]=fs->make<TH1F>("jet0Phi_case2_trigger2",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case2_trigger2"]=fs->make<TH1F>("jet1Phi_case2_trigger2",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case2_trigger2"]=fs->make<TH1F>("dijetMass_case2_trigger2", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case1"]=fs->make<TH1F>("tauPt_case1",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case1"]=fs->make<TH1F>("tauEta_case1",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case1"]=fs->make<TH1F>("tauPhi_case1",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case1"]=fs->make<TH1F>("metPt_case1",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case1"]=fs->make<TH1F>("metEta_case1",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case1"]=fs->make<TH1F>("metPhi_case1",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case1"]=fs->make<TH1F>("jet0Pt_case1",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case1"]=fs->make<TH1F>("jet1Pt_case1",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case1"]=fs->make<TH1F>("jet0Eta_case1",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case1"]=fs->make<TH1F>("jet1Eta_case1",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case1"]=fs->make<TH1F>("jet0Phi_case1",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case1"]=fs->make<TH1F>("jet1Phi_case1",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case1"]=fs->make<TH1F>("dijetMass_case1", "mass",    500,   0., 2500.);

         histContainer_["tauPt_case2"]=fs->make<TH1F>("tauPt_case2",   "tau Pt", 300, 0, 600);
         histContainer_["tauEta_case2"]=fs->make<TH1F>("tauEta_case2",   "tau Eta",100, -5,  5);
         histContainer_["tauPhi_case2"]=fs->make<TH1F>("tauPhi_case2",   "tau Phi", 100, -3.5, 3.5);
         histContainer_["metPt_case2"]=fs->make<TH1F>("metPt_case2",   "met Pt", 300, 0, 600);
         histContainer_["metEta_case2"]=fs->make<TH1F>("metEta_case2",   "met Eta",100, -5,  5);
         histContainer_["metPhi_case2"]=fs->make<TH1F>("metPhi_case2",   "met Phi", 100, -3.5, 3.5);
         histContainer_["jet0Pt_case2"]=fs->make<TH1F>("jet0Pt_case2",   "jet0 Pt", 300, 0, 600);
         histContainer_["jet1Pt_case2"]=fs->make<TH1F>("jet1Pt_case2",   "jet1 Pt", 300, 0, 600);
         histContainer_["jet0Eta_case2"]=fs->make<TH1F>("jet0Eta_case2",   "jet0 Eta",  100, -5,  5);
         histContainer_["jet1Eta_case2"]=fs->make<TH1F>("jet1Eta_case2",   "jet1 Eta",  100, -5,  5);
         histContainer_["jet0Phi_case2"]=fs->make<TH1F>("jet0Phi_case2",   "jet0 Phi",     100, -3.5, 3.5);
         histContainer_["jet1Phi_case2"]=fs->make<TH1F>("jet1Phi_case2",   "jet1 Phi",     100, -3.5, 3.5);
         histContainer_["dijetMass_case2"]=fs->make<TH1F>("dijetMass_case2", "mass",    500,   0., 2500.);
     



}

void 
TriggerCuts3::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerCuts3);

