#include "HLTJetMETNtupleProducer.h"
#include <iostream>

bool applyMETFilters_;

edm::EDGetTokenT<bool> badMuonFilterToken_;
edm::EDGetTokenT<bool> badChargedCandidateFilterToken_; 

//edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
edm::EDGetTokenT<edm::TriggerResults> triggerBitsPAT_;

HLTJetMETNtupleProducer::HLTJetMETNtupleProducer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  runJets_ = iConfig.getParameter<bool>("runJets");
  isData_ = iConfig.getParameter<bool>("isData");
  runMets_ = iConfig.getUntrackedParameter<bool>("runMets", false);
  PVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag"));
  MetCollectionToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCollectionTag"));
  if(!isData_){
    GenMetCollectionToken_ = consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("GenMetCollectionTag"));
    GenJetCollectionToken_= consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("GenJetCollectionTag"));
    CaloJetCollectionToken_= consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("CaloJetCollectionTag"));
  }
  applyMETFilters_      = iConfig.getParameter<bool>("applyMETFilters");
  badMuonFilterToken_   = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"));
  badChargedCandidateFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));
  MuonCollectionToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonCollectionTag"));
  ElectronCollectionToken_ = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("ElectronCollectionTag"));
  if(runJets_ || (!isData_)){
    PFJetCollectionToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("PFJetCollectionTag"));
  }
  if(runJets_){
    HLTPFJetCollectionToken_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("HLTPFJetCollectionTag"));
    HLTCaloJetCollectionToken_ = consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("HLTCaloJetCollectionTag"));
  }
  if(runMets_){
    HLTPFMetCollectionToken_ = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("HLTPFMetTag"));
    HLTPFMetType1CollectionToken_ = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("HLTPFMetType1Tag"));
    HLTCaloMetCollectionToken_ = consumes<reco::CaloMETCollection>(iConfig.getParameter<edm::InputTag>("HLTCaloMetTag"));
  }
  eleMvaSpring16WPMediumMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaSpring16WPMediumMap"));
  eleMvaSpring16WPTightMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMvaSpring16WPTightMap"));
  mvaSpring16ValuesMapToken_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaSpring16ValuesMap"));
  mvaSpring16CategoriesMapToken_ = consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaSpring16CategoriesMap"));
  eleSummer16VetoIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleSummer16VetoIdMap"));
  eleSummer16LooseIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleSummer16LooseIdMap"));
  eleSummer16MediumIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleSummer16MediumIdMap"));
  eleSummer16TightIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleSummer16TightIdMap"));
  hlt_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltprocess"));
  triggerPaths_ = iConfig.getUntrackedParameter<std::vector<std::string> >("triggerPaths");
  //triggerBits_          = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")));
  triggerBitsPAT_       = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("PAT")));
  triggerObjects_       = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>( "triggerObjects" ));
}


HLTJetMETNtupleProducer::~HLTJetMETNtupleProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HLTJetMETNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ///event info
  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.getLuminosityBlock().luminosityBlock();
  
  math::XYZPoint pv_position = math::XYZPoint(0.,0.,0.);
  edm::Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(PVToken_, Vertex);
  nPV_ = 0;
  if(Vertex.isValid()) {
    for(unsigned i = 0 ; i < Vertex->size(); i++) {
      nPV_++;
      if(i == 0) {
	PVx_ = (*Vertex)[i].x();
	PVy_ = (*Vertex)[i].y();
	PVz_ = (*Vertex)[i].z();
	pv_position = (*Vertex)[i].position();
      }
    }
  }

  edm::Handle<edm::TriggerResults> triggerBitsPAT;
  iEvent.getByToken(triggerBitsPAT_,triggerBitsPAT);
  edm::TriggerNames namesPAT;
  if( triggerBitsPAT.isValid() )
    {
      namesPAT = iEvent.triggerNames(*triggerBitsPAT);
    }

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(hlt_,hltresults);
  const edm::TriggerNames& TrigNames_ = iEvent.triggerNames(*hltresults);
  const int ntrigs = hltresults->size();
  //store trigger paths that are passed in the event
  triggerResults_.clear();
  for(int itrig=0; itrig<ntrigs; itrig++)
    {
      if(!hltresults->wasrun(itrig) )continue;
      std::string trigName = TrigNames_.triggerName(itrig);
      if(triggerPaths_.size() > 0){
        for(size_t ip = 0; ip < triggerPaths_.size(); ip++){
          if(trigName.find(triggerPaths_[ip]) != std::string::npos){
            if(hltresults->accept(itrig)){
              triggerResults_.push_back(trigName);
            }
          }
        }
      }
    }

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  passedL1MET_ = 0;
  passedHLTCaloMET_ = 0;
  passedHLTCaloMETClean_ = 0;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(TrigNames_);
    //obj.unpackFilterLabels(iEvent, *hltresults);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
      //if (obj.hasFilterLabel("hltL1sETM50IorETM60IorETM70IorETM80IorETM90IorETM100") || obj.hasFilterLabel("hltL1sETM50ToETM120")) passedL1MET_ = 1;
      if (obj.hasFilterLabel("hltL1sAllETMHFSeeds") || obj.hasFilterLabel("hltL1sAllETMHadSeeds") || obj.hasFilterLabel("hltL1sETM50ToETM120")) passedL1MET_ = 1; //2018, 2017, or 2016 filter name
      if (obj.hasFilterLabel("hltMET90")) passedHLTCaloMET_ = 1;
      if (obj.hasFilterLabel("hltMETClean80")) passedHLTCaloMETClean_ = 1;
    }             
  }

  if(runJets_ || (!isData_)){
    //Get offline Reco/PAT Jets
    edm::Handle<pat::JetCollection> offlinepfjets;
    iEvent.getByToken(PFJetCollectionToken_,offlinepfjets);
    pfjetPt_.clear(); pfjetEta_.clear(); pfjetPhi_.clear();
    for(pat::JetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) {  
      if (ipf->pt() > 20){
	pfjetPt_.push_back(ipf->pt());
	pfjetEta_.push_back(ipf->eta()); 
	pfjetPhi_.push_back(ipf->phi()); 
      }
    }
  }//if runJets_ || (!isData_)
  if(!isData_){
    edm::Handle<reco::CaloJetCollection> offlinecalojets;
    iEvent.getByToken(CaloJetCollectionToken_,offlinecalojets);
    calojetPt_.clear(); calojetEta_.clear(); calojetPhi_.clear();
    for(reco::CaloJetCollection::const_iterator icalo = offlinecalojets->begin();icalo != offlinecalojets->end(); ++icalo) {
      if (icalo->pt() > 20){
	calojetPt_ .push_back(icalo->pt());
	calojetEta_.push_back(icalo->eta());
	calojetPhi_.push_back(icalo->phi());
      }
    }
  }
  //Get GenJets
  if(!isData_){
    edm::Handle<reco::GenJetCollection> genjets;
    iEvent.getByToken(GenJetCollectionToken_,genjets);
    genjetPt_.clear(); genjetEta_.clear(); genjetPhi_.clear();
    for(reco::GenJetCollection::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
      if (igen->pt() > 20){
        genjetPt_ .push_back(igen->pt());
        genjetEta_.push_back(igen->eta());
        genjetPhi_.push_back(igen->phi());
      }
    }
  }
  
  //MET filters
  edm::Handle<bool> badMuonFilter;
  //edm::Handle<bool> badChargedCandidateFilter;
  
  iEvent.getByToken(badMuonFilterToken_,badMuonFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_,badChargedCandidateFilter);
  
  bool passMETFilters = 1;
  bool pass_HBHENoiseFilter = 1;
  bool pass_HBHENoiseIsoFilter = 1;
  bool pass_EcalDeadCellTriggerPrimitiveFilter = 1;
  bool pass_goodVertices = 1;
  bool pass_eeBadScFilter = 1;
  bool pass_globalTightHalo2016Filter = 1;
  bool pass_badMuonFilter = *badMuonFilter;
  //bool pass_badChargedCandidateFilter = *badChargedCandidateFilter;
  
  
  if( triggerBitsPAT.isValid() )
    {
      for (unsigned int i = 0, n = triggerBitsPAT->size(); i < n; ++i)
	{
	  std::string triggerName = namesPAT.triggerName(i);
	  
	  std::cout<<"pat filter names "<<triggerName<<std::endl;
	  bool isFired = (triggerBitsPAT->accept(i) ? true : false);
	  
	  if( strcmp(triggerName.c_str(),"Flag_HBHENoiseFilter") == 0 )
	    {
	      if( !isFired ) pass_HBHENoiseFilter = 0;
	    }
	  else if( strcmp(triggerName.c_str(),"Flag_HBHENoiseIsoFilter") == 0 )
	    {
	      if( !isFired ) pass_HBHENoiseIsoFilter = 0;
	       }
	  else if( strcmp(triggerName.c_str(),"Flag_EcalDeadCellTriggerPrimitiveFilter") == 0 )
	    {
	      if( !isFired ) pass_EcalDeadCellTriggerPrimitiveFilter = 0;
	    }
	  else if( strcmp(triggerName.c_str(),"Flag_goodVertices") == 0 )
	    {
	      if( !isFired ) pass_goodVertices = 0;
	    }
	  else if( strcmp(triggerName.c_str(),"Flag_eeBadScFilter") == 0 )
	    {
	      if( !isFired ) pass_eeBadScFilter = 0;
	    }
	  else if( strcmp(triggerName.c_str(),"Flag_globalTightHalo2016Filter") == 0 )
	    {
	      if( !isFired ) pass_globalTightHalo2016Filter = 0;
	    }
	}
    }
   
   //std::cout<<"bad muon filter "<<pass_badMuonFilter<<std::endl;
   //std::cout<<"bad ch. cand filter "<<pass_badChargedCandidateFilter<<std::endl;
  passMETFilters = (pass_HBHENoiseFilter &&
		    pass_HBHENoiseIsoFilter &&
		    pass_EcalDeadCellTriggerPrimitiveFilter &&
		    pass_goodVertices &&
		    pass_eeBadScFilter &&
		    pass_globalTightHalo2016Filter &&
		    pass_badMuonFilter);
                    //&& pass_badChargedCandidateFilter);
  passMETFilter_ = passMETFilters;
  
  if(runJets_){
    //Get hlt jets
    edm::Handle<reco::PFJetCollection> hltpfjets;
    iEvent.getByToken(HLTPFJetCollectionToken_,hltpfjets);
    hltpfjetPt_.clear(); hltpfjetEta_.clear(); hltpfjetPhi_.clear();
    for(reco::PFJetCollection::const_iterator ijet = hltpfjets->begin();ijet != hltpfjets->end(); ++ijet) {  
      if (ijet->pt() > 20){
	hltpfjetPt_.push_back(ijet->pt());
	hltpfjetEta_.push_back(ijet->eta()); 
	hltpfjetPhi_.push_back(ijet->phi());
      }
    }
    edm::Handle<reco::CaloJetCollection> hltcalojets;
    iEvent.getByToken(HLTCaloJetCollectionToken_,hltcalojets);
    hltcalojetPt_.clear(); hltcalojetEta_.clear(); hltcalojetPhi_.clear();
    for(reco::CaloJetCollection::const_iterator ijet = hltcalojets->begin();ijet != hltcalojets->end(); ++ijet) {
      if (ijet->pt() > 20){
	hltcalojetPt_.push_back(ijet->pt());
	hltcalojetEta_.push_back(ijet->eta());
	hltcalojetPhi_.push_back(ijet->phi());
      }
    }
  }//if runJets_

  if(runMets_){
    //hlt met
    edm::Handle<reco::PFMETCollection> hltMETs;
    iEvent.getByToken(HLTPFMetCollectionToken_, hltMETs);
    hltPFMetPt_ = 0.; hltPFMetPhi_ = -999.;
    if(hltMETs.isValid() && hltMETs->size() > 0){
      hltPFMetPt_ = (*hltMETs)[0].pt();
      hltPFMetPhi_ = (*hltMETs)[0].phi();
    }
    edm::Handle<reco::PFMETCollection> hltMETsType1;
    iEvent.getByToken(HLTPFMetType1CollectionToken_, hltMETsType1);
    hltPFMetType1Pt_ = 0.; hltPFMetType1Phi_ = -999.;
    if(hltMETsType1.isValid() && hltMETsType1->size() > 0){
      hltPFMetType1Pt_ = (*hltMETsType1)[0].pt();
      hltPFMetType1Phi_ = (*hltMETsType1)[0].phi();
    }
    edm::Handle<reco::CaloMETCollection> hltCaloMETs;
    iEvent.getByToken(HLTCaloMetCollectionToken_, hltCaloMETs);
    hltCaloMetPt_ = 0.; hltCaloMetPhi_ = -999.;
    if(hltCaloMETs.isValid() && hltCaloMETs->size() > 0){
      hltCaloMetPt_ = (*hltCaloMETs)[0].pt();
      hltCaloMetPhi_ = (*hltCaloMETs)[0].phi();
    }
  }// if runMets_

  edm::Handle<pat::METCollection> METs;
  iEvent.getByToken(MetCollectionToken_, METs);
  
  metPt_ = 0.; metPhi_ = -999.;
  metPx_ = 0.; metPy_ = 0.;
  caloMetPt_ = 0.; caloMetPhi_ = -999.;
  genMetPt_ = 0.; genMetPhi_ = -999.;
  if(METs.isValid() && METs->size() > 0){
    metPt_ = (*METs)[0].pt();
    metPhi_ = (*METs)[0].phi();
    metPx_ = (*METs)[0].px();
    metPy_ = (*METs)[0].py();
    caloMetPt_ = (*METs)[0].caloMETPt();
    caloMetPhi_ = (*METs)[0].caloMETPhi();
    if(!isData_){
      genMetPt_ = (*METs)[0].genMET()->pt();
      genMetPhi_ = (*METs)[0].genMET()->phi();
    }
  }
  /*if(!isData_){
    edm::Handle<reco::GenMETCollection> genMETs;
    iEvent.getByToken(GenMetCollectionToken_, genMETs);
    if(genMETs.isValid() && genMETs->size() > 0){
      genMetPt_ = (*genMETs)[0].pt();
      genMetPhi_ = (*genMETs)[0].phi();
    }
  }
  */

  muonPx_.clear(); muonPy_.clear(); muonPz_.clear(); 
  muonPt_.clear(); muonEta_.clear(); muonPhi_.clear(); muonCharge_.clear();
  muonR04SumChargedHadronPt_.clear(); muonR04SumChargedParticlePt_.clear();
  muonR04SumNeutralHadronEt_.clear(); muonR04SumPhotonEt_.clear();
  muonR04SumPUPt_.clear(); muonIsPF_.clear();
  muonIsGlobal_.clear(); muonIsTracker_.clear(); muonIsICHEPMedium_.clear();
  muonDz_.clear(); muonDxy_.clear(); //muonNormChi2_.clear();
  
  edm::Handle<pat::MuonCollection> Muons;
  iEvent.getByToken(MuonCollectionToken_, Muons);
  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){
	if ((*Muons)[i].pt() < 10.) continue; //require a minimim cut
	muonPx_.push_back((*Muons)[i].px());
	muonPy_.push_back((*Muons)[i].py());
	muonPz_.push_back((*Muons)[i].pz());
	muonPt_.push_back((*Muons)[i].pt());
	muonEta_.push_back((*Muons)[i].eta());
	muonPhi_.push_back((*Muons)[i].phi());
	muonCharge_.push_back((*Muons)[i].charge());
	muonR04SumChargedHadronPt_.push_back((*Muons)[i].pfIsolationR04().sumChargedHadronPt);
	muonR04SumChargedParticlePt_.push_back((*Muons)[i].pfIsolationR04().sumChargedParticlePt);
	muonR04SumNeutralHadronEt_.push_back((*Muons)[i].pfIsolationR04().sumNeutralHadronEt);
	muonR04SumPhotonEt_.push_back((*Muons)[i].pfIsolationR04().sumPhotonEt);
	muonR04SumPUPt_.push_back((*Muons)[i].pfIsolationR04().sumPUPt);
	muonIsPF_.push_back((*Muons)[i].isPFMuon());
	muonIsGlobal_.push_back((*Muons)[i].isGlobalMuon());
	muonIsTracker_.push_back((*Muons)[i].isTrackerMuon());
	muonIsICHEPMedium_.push_back(isICHEPMuon((*Muons)[i]));
	reco::TrackRef bestTrack  = (*Muons)[i].muonBestTrack();
	if (bestTrack.isNonnull()) {
	  muonDxy_.push_back(bestTrack->dxy(pv_position));
	  muonDz_.push_back(bestTrack->dz(pv_position));
	}
	else{
	  muonDxy_.push_back(-9999.);
	  muonDz_.push_back(-9999.);
	}
      }
    }

  edm::Handle<edm::View<pat::Electron> > Electrons;
  iEvent.getByToken(ElectronCollectionToken_, Electrons);

  // cut based                                                                                                                                                                           
  edm::Handle<edm::ValueMap<bool> > summer16_veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > summer16_loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > summer16_medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > summer16_tight_id_decisions;
  iEvent.getByToken(eleSummer16VetoIdMapToken_,summer16_veto_id_decisions);
  iEvent.getByToken(eleSummer16LooseIdMapToken_,summer16_loose_id_decisions);
  iEvent.getByToken(eleSummer16MediumIdMapToken_,summer16_medium_id_decisions);
  iEvent.getByToken(eleSummer16TightIdMapToken_,summer16_tight_id_decisions);

  // mva
  edm::Handle<edm::ValueMap<bool> > spring16_medium_decisions;
  edm::Handle<edm::ValueMap<bool> > spring16_tight_decisions;
  iEvent.getByToken(eleMvaSpring16WPMediumMapToken_,spring16_medium_decisions);
  iEvent.getByToken(eleMvaSpring16WPTightMapToken_,spring16_tight_decisions);

  edm::Handle<edm::ValueMap<float> > mvaSpring16Values;
  edm::Handle<edm::ValueMap<int> > mvaSpring16Categories;
  iEvent.getByToken(mvaSpring16ValuesMapToken_,mvaSpring16Values);
  iEvent.getByToken(mvaSpring16CategoriesMapToken_,mvaSpring16Categories);

  elecPx_.clear(); elecPy_.clear(); elecPz_.clear(); elecPt_.clear(); elecEta_.clear();
  elecPhi_.clear(); elecCharge_.clear(); elecR03SumChargedHadronPt_.clear(); elecR03SumChargedParticlePt_.clear();
  elecR03SumNeutralHadronEt_.clear(); elecR03SumPhotonEt_.clear(); elecR03SumPUPt_.clear();
  elec_mva_value_Spring16_v1_.clear(); elec_mva_category_Spring16_v1_.clear();
  elec_mva_medium_Spring16_v1_.clear(); elec_mva_tight_Spring16_v1_.clear();
  elec_cutId_veto_Summer16_.clear(); elec_cutId_loose_Summer16_.clear();
  elec_cutId_medium_Summer16_.clear(); elec_cutId_tight_Summer16_.clear();
  elec_pass_conversion_.clear(); elec_nmissinginnerhits_.clear();
  elecDxy_.clear(); elecDz_.clear();
  
  if(Electrons.isValid())
    {
      for(unsigned i = 0 ; i < Electrons->size() ; i++){
	const auto el = Electrons->ptrAt(i);

	if (el->pt()<10) continue;

	elecPx_.push_back(el->px());
	elecPy_.push_back(el->py());
	elecPz_.push_back(el->pz());
	elecPt_.push_back(el->pt());
	elecEta_.push_back(el->eta());
	elecPhi_.push_back(el->phi()); 
	elecCharge_.push_back(el->charge());
	elecR03SumChargedHadronPt_.push_back(el->pfIsolationVariables().sumChargedHadronPt);
	elecR03SumChargedParticlePt_.push_back(el->pfIsolationVariables().sumChargedParticlePt);
	elecR03SumNeutralHadronEt_.push_back(el->pfIsolationVariables().sumNeutralHadronEt);
	elecR03SumPhotonEt_.push_back(el->pfIsolationVariables().sumPhotonEt);
	elecR03SumPUPt_.push_back(el->pfIsolationVariables().sumPUPt);
	reco::GsfTrackRef gsfTr_e = el->gsfTrack();
	//elec_nmissinginnerhits_.push_back(gsfTr_e->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
	elec_pass_conversion_.push_back((*Electrons)[i].passConversionVeto());
	elecDxy_.push_back(gsfTr_e->dxy(pv_position));
	elecDz_.push_back(gsfTr_e->dz(pv_position));
	elec_mva_medium_Spring16_v1_.push_back((*spring16_medium_decisions)[el]);
        elec_mva_tight_Spring16_v1_.push_back((*spring16_tight_decisions)[el]);
	elec_mva_value_Spring16_v1_.push_back((*mvaSpring16Values)[el]);
        elec_mva_category_Spring16_v1_.push_back((*mvaSpring16Categories)[el]);
	elec_cutId_veto_Summer16_.push_back((*summer16_veto_id_decisions)[el]);
        elec_cutId_loose_Summer16_.push_back((*summer16_loose_id_decisions)[el]);
        elec_cutId_medium_Summer16_.push_back((*summer16_medium_id_decisions)[el]);
        elec_cutId_tight_Summer16_.push_back((*summer16_tight_id_decisions)[el]);
      }
    }
  /*
     if( (applyMETFilters_ && passMETFilters) || !applyMETFilters_ ) // MET filters    
      {  
	//Fill tree
	tree_->Fill();
      }
  */

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HLTJetMETNtupleProducer::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");

  tree_->Branch("run", &run_, "run/l");
  tree_->Branch("event", &event_, "event/l");
  tree_->Branch("lumi", &lumi_, "lumi/l");

  tree_->Branch("nPV", &nPV_, "nPV/i"); 
  tree_->Branch("PVx", &PVx_, "PVx/F");
  tree_->Branch("PVy", &PVy_, "PVy/F");
  tree_->Branch("PVz", &PVz_, "PVz/F");

  if(runJets_ || (!isData_)){
    tree_->Branch("pfjetPt", "std::vector<float>", &pfjetPt_);
    tree_->Branch("pfjetEta", "std::vector<float>", &pfjetEta_);
    tree_->Branch("pfjetPhi", "std::vector<float>", &pfjetPhi_);
  }
  if(runJets_){
    tree_->Branch("hltpfjetPt", "std::vector<float>", &hltpfjetPt_);
    tree_->Branch("hltpfjetEta", "std::vector<float>", &hltpfjetEta_);
    tree_->Branch("hltpfjetPhi", "std::vector<float>", &hltpfjetPhi_);
    tree_->Branch("hltcalojetPt", "std::vector<float>", &hltcalojetPt_);
    tree_->Branch("hltcalojetEta", "std::vector<float>", &hltcalojetEta_);
    tree_->Branch("hltcalojetPhi", "std::vector<float>", &hltcalojetPhi_);
  }
  if(runMets_){
    tree_->Branch("hltPFMetPt", &hltPFMetPt_, "hltPFMetPt/F");
    tree_->Branch("hltPFMetPhi", &hltPFMetPhi_, "hltPFMetPhi/F");
    tree_->Branch("hltPFMetType1Pt", &hltPFMetType1Pt_, "hltPFMetType1Pt/F");
    tree_->Branch("hltPFMetType1Phi", &hltPFMetType1Phi_, "hltPFMetType1Phi/F");
    tree_->Branch("hltCaloMetPt", &hltCaloMetPt_, "hltCaloMetPt/F");
    tree_->Branch("hltCaloMetPhi", &hltCaloMetPhi_, "hltCaloMetPhi/F");
  }
  if(!isData_){
    tree_->Branch("calojetPt", "std::vector<float>", &calojetPt_);
    tree_->Branch("calojetEta", "std::vector<float>", &calojetEta_);
    tree_->Branch("calojetPhi", "std::vector<float>", &calojetPhi_);
    tree_->Branch("genjetPt", "std::vector<float>", &genjetPt_);
    tree_->Branch("genjetEta", "std::vector<float>", &genjetEta_);
    tree_->Branch("genjetPhi", "std::vector<float>", &genjetPhi_);
    tree_->Branch("genMetPt", &genMetPt_, "genMetPt/F");
    tree_->Branch("genMetPhi", &genMetPhi_, "genMetPhi/F");
  }

  tree_->Branch("metPx", &metPx_, "metPx/F");
  tree_->Branch("metPy", &metPy_, "metPy/F");
  tree_->Branch("metPt", &metPt_, "metPt/F");
  tree_->Branch("metPhi", &metPhi_, "metPhi/F");
  tree_->Branch("caloMetPt", &caloMetPt_, "caloMetPt/F");
  tree_->Branch("caloMetPhi", &caloMetPhi_, "caloMetPhi/F");
  tree_->Branch("passMETFilter", &passMETFilter_, "passMETFilter/O");
  
  tree_->Branch("muonPx", "std::vector<float>", &muonPx_);
  tree_->Branch("muonPy", "std::vector<float>",&muonPy_);
  tree_->Branch("muonPz", "std::vector<float>",&muonPz_);
  tree_->Branch("muonPt", "std::vector<float>",&muonPt_);
  tree_->Branch("muonEta", "std::vector<float>",&muonEta_);
  tree_->Branch("muonPhi", "std::vector<float>",&muonPhi_);
  tree_->Branch("muonCharge", "std::vector<float>",&muonCharge_);
  tree_->Branch("muonR04SumChargedHadronPt", "std::vector<float>",&muonR04SumChargedHadronPt_);
  tree_->Branch("muonR04SumChargedParticlePt", "std::vector<float>",&muonR04SumChargedParticlePt_);
  tree_->Branch("muonR04SumNeutralHadronEt", "std::vector<float>",&muonR04SumNeutralHadronEt_);
  tree_->Branch("muonR04SumPhotonEt", "std::vector<float>",&muonR04SumPhotonEt_);
  tree_->Branch("muonR04SumPUPt", "std::vector<float>",&muonR04SumPUPt_);
  tree_->Branch("muonIsPF", "std::vector<bool>", &muonIsPF_);
  tree_->Branch("muonIsGlobal", "std::vector<bool>", &muonIsGlobal_);
  tree_->Branch("muonIsTracker", "std::vector<bool>", &muonIsTracker_);
  tree_->Branch("muonIsICHEPMedium", "std::vector<bool>", &muonIsICHEPMedium_);
  tree_->Branch("muonDz", "std::vector<float>", &muonDz_);
  tree_->Branch("muonDxy", "std::vector<float>",&muonDxy_);
  //tree_->Branch("muonNormChi2", "std::vector<float>",&muonNormChi2_);

  tree_->Branch("elecPx", "std::vector<float>", &elecPx_);
  tree_->Branch("elecPy", "std::vector<float>",&elecPy_);
  tree_->Branch("elecPz", "std::vector<float>",&elecPz_);
  tree_->Branch("elecPt", "std::vector<float>",&elecPt_);
  tree_->Branch("elecEta", "std::vector<float>",&elecEta_);
  tree_->Branch("elecPhi", "std::vector<float>",&elecPhi_);
  tree_->Branch("elecDxy", "std::vector<float>", &elecDxy_);
  tree_->Branch("elecDz", "std::vector<float>",&elecDz_);
  tree_->Branch("elecCharge", "std::vector<float>",&elecCharge_);
  tree_->Branch("elecR03SumChargedHadronPt", "std::vector<float>",&elecR03SumChargedHadronPt_);
  tree_->Branch("elecR03SumChargedParticlePt", "std::vector<float>",&elecR03SumChargedParticlePt_);
  tree_->Branch("elecR03SumNeutralHadronEt", "std::vector<float>",&elecR03SumNeutralHadronEt_);
  tree_->Branch("elecR03SumPhotonEt", "std::vector<float>",&elecR03SumPhotonEt_);
  tree_->Branch("elecR03SumPUPt", "std::vector<float>",&elecR03SumPUPt_);
  tree_->Branch("elec_mva_medium_Spring16_v1", "std::vector<bool>", &elec_mva_medium_Spring16_v1_);
  tree_->Branch("elec_mva_tight_Spring16_v1", "std::vector<bool>", &elec_mva_tight_Spring16_v1_);
  tree_->Branch("elec_mva_value_Spring16_v1", "std::vector<float>", &elec_mva_value_Spring16_v1_);
  tree_->Branch("elec_mva_category_Spring16_v1", "std::vector<int>", &elec_mva_category_Spring16_v1_);
  tree_->Branch("elec_cutId_veto_Summer16", "std::vector<bool>", &elec_cutId_veto_Summer16_);
  tree_->Branch("elec_cutId_loose_Summer16", "std::vector<bool>", &elec_cutId_loose_Summer16_);
  tree_->Branch("elec_cutId_medium_Summer16", "std::vector<bool>", &elec_cutId_medium_Summer16_);
  tree_->Branch("elec_cutId_tight_Summer16", "std::vector<bool>", &elec_cutId_tight_Summer16_);
  tree_->Branch("elec_nmissinginnerhits", "std::vector<unsigned int>", &elec_nmissinginnerhits_);
  tree_->Branch("elec_pass_conversion", "std::vector<bool>", &elec_pass_conversion_);
  tree_->Branch("triggerResults", "std::vector<std::string>", &triggerResults_);

  tree_->Branch("passedL1MET", &passedL1MET_,"passedL1MET/i"); 
  tree_->Branch("passedHLTCaloMET", &passedHLTCaloMET_, "passedHLTCaloMET/i");
  tree_->Branch("passedHLTCaloMETClean", &passedHLTCaloMETClean_, "passedHLTCaloMETClean/i");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTJetMETNtupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HLTJetMETNtupleProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HLTJetMETNtupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HLTJetMETNtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HLTJetMETNtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTJetMETNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool HLTJetMETNtupleProducer::isICHEPMuon(const reco::Muon& recoMu){
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.49 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}
 
//define this as a plug-in
DEFINE_FWK_MODULE(HLTJetMETNtupleProducer);
