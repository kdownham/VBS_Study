#include "SampleAnalyzer/User/Analyzer/zz_polarization.h"
#include <vector>
#include <algorithm>
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool zz_polarization::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  PHYSICS->mcConfig().Reset();
  Double_t Lower[6] = {140.0, 180.0, 220.0, 260.0, 300.0, 800.0};
  // Initialize a new histogram called "myHisto", with label "ZTZ0" with 20 bins ranging from -1 to 1 on the x axis
  myHisto = new TH1F("myHisto", "Z0ZT", 40, -1., 1.);
  mZZ = new TH1F("mZZ","Z0ZT",5,Lower);
  //mZZ = new TH1F("mZZ","ZZ",50,0.,1000.);

  cout << "END Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void zz_polarization::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  // Write the output of Madanalysis to a fbile with the name "output_zTz0.root"
  TFile* output = new TFile("output_z0zT.root", "RECREATE");
  // Set the title of the X and Y axes, respectively
  myHisto->GetXaxis()->SetTitle("cos#theta^{*}");
  myHisto->GetYaxis()->SetTitle("d#sigma/dcos#theta^{*}");
  mZZ->GetXaxis()->SetTitle("m_{ZZ}");
  mZZ->GetYaxis()->SetTitle("d#sigma/dm_{ZZ}");
  //myHisto->SetMaximum( 1.0 );
  myHisto->SetMinimum( 0 );
  mZZ->SetMinimum( 0 );
  //Scale the histogram to the fraction of events passing selection
  myHisto->Scale(1.0/20000.0);
  mZZ->Scale(1.0/10000.0);
  myHisto->Write();
  mZZ->Write();
  output->Close();
  cout << "END Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool zz_polarization::Execute(SampleFormat& sample, const EventFormat& event)
{
  
  // Initialization of pointers for positive and negative leptons and the z boson

  const MCParticleFormat* z = 0;
  const MCParticleFormat* lepton_p = 0;
  const MCParticleFormat* lepton_n = 0;
  const MCParticleFormat* lepton_p2 = 0;
  const MCParticleFormat* lepton_n2 = 0;
  const MCParticleFormat* z2 = 0;
  //const MCParticleFormat* zz = 0;
  //int lep_candidates = 0;
  int positive = 0;
  int negative = 0;

  std::vector<const MCParticleFormat*> leppt (4);
  //std::vector<double>::size_type sz = leppt.size();
  //std::vector<double>::iterator result;

  // Identification of the particles of interest
  
  if (event.mc()!=0){

  //sumWeights += event.mc()->weight();
  std::cout << "Event weight = " << event.mc()->weight() << endl;

  //std::cout << "Number of particles in this event = " << event.mc()->particles().size() << endl;

  for (unsigned int i=0;i<event.mc()->particles().size();i++){


    const MCParticleFormat* prt=&(event.mc()->particles()[i]);

    //std::cout << "Particle ID = " << prt->pdgid() << endl;

    // Making sure the lepton is a final state particle
    if (!PHYSICS->Id->IsFinalState(prt)) continue;

    // Lepton selection based on PDG-id
    if ( std::abs(prt->pdgid())!=11 && std::abs(prt->pdgid())!=13 ) continue;

    // Getting the mother of the lepton and checking if it is a Z boson
    const MCParticleFormat* mother = prt->mothers()[0];

    if (mother==0) continue;
    if (std::abs(mother->pdgid())!=23) continue;

    // Looking for a positive lepton in the event

      if (prt->pdgid() > 0 && positive < 1){
	  lepton_p = prt;
	  z = mother;
	  positive++;
	  //std::cout << "First positive lepton found" << endl;
     }
     
      if (prt->pdgid() > 0 && positive > 0){
	  lepton_p2 = prt;
	  z2 = mother;
	  //std::cout << "Second positive lepton found" << endl;	  
     }

      if (prt->pdgid() < 0 && negative < 1){
	  lepton_n = prt;
	  negative++;
	  //std::cout << "First negative lepton found" << endl;
     }

      if (prt->pdgid() < 0 && negative > 0){
	  lepton_n2 = prt;
	  //std::cout << "Second negative lepton found" << endl;
     }

  }
  
  }

    //std::cout << "Number of candidate leptons = " << lep_candidates << endl;

    if (lepton_p==0){
        WARNING << "(z > l+ l-) decay chain not found!" << std::endl;
        return false;
       }

    leppt[0] = lepton_p;
    leppt[1] = lepton_n;
    leppt[2] = lepton_p2;
    leppt[3] = lepton_n2;

    SORTER->sort(leppt, PTordering);

    // Need a method here to connect the pT ordered leptons to a specific pdg id
    //

    bool leading_lep_e = std::abs(leppt[0]->pdgid()) == 11 && leppt[0]->momentum().Pt() > 20.0 && std::abs(leppt[0]->eta()) < 2.5;
    bool leading_lep_m = std::abs(leppt[0]->pdgid()) == 13 && leppt[0]->momentum().Pt() > 20.0 && std::abs(leppt[0]->eta()) < 2.4;
    bool second_lep_e = std::abs(leppt[1]->pdgid()) == 11 && leppt[1]->momentum().Pt() > 12.0 && std::abs(leppt[1]->eta()) < 2.5;
    bool second_lep_m = std::abs(leppt[1]->pdgid()) == 13 && leppt[1]->momentum().Pt() > 10.0 && std::abs(leppt[1]->eta()) < 2.4;
    bool third_lep_e = std::abs(leppt[2]->pdgid()) == 11 && leppt[2]->momentum().Pt() > 7.0 && std::abs(leppt[2]->eta()) < 2.5;
    bool third_lep_m = std::abs(leppt[2]->pdgid()) == 13 && leppt[2]->momentum().Pt() > 5.0 && std::abs(leppt[2]->eta()) < 2.4;
    bool fourth_lep_e = std::abs(leppt[3]->pdgid()) == 11 && leppt[3]->momentum().Pt() > 7.0 && std::abs(leppt[3]->eta()) < 2.5;
    bool fourth_lep_m = std::abs(leppt[3]->pdgid()) == 13 && leppt[3]->momentum().Pt() > 5.0 && std::abs(leppt[3]->eta()) < 2.4;

    double z1_mass = 0;
    double z2_mass = 0;

    double mdiff_z1 = std::abs(z->momentum().M() - 91.2);
    double mdiff_z2 = std::abs(z2->momentum().M() - 91.2);

    if (mdiff_z1 < mdiff_z2){
	z1_mass = z->momentum().M();	
        z2_mass = z2->momentum().M();
   }else{    
        z1_mass = z2->momentum().M();
        z2_mass = z->momentum().M();
    }

    //std::cout << "Z1 mass = " << z1_mass << endl;
    //std::cout << "Z2 mass = " << z2_mass << endl;

    bool z1_range = z1_mass < 120.0 && z1_mass > 40.0;
    bool z2_range = z2_mass < 120.0 && z2_mass > 4.0;

    if ( (leading_lep_e || leading_lep_m) ){ 		
	if (second_lep_e || second_lep_m){
	    if (third_lep_e || third_lep_m){
		if (fourth_lep_e || fourth_lep_m){

			// Computing the observable; filling the histogram

    			// Boosting the lepton four-momentum to the Z rest frame
    			MCParticleFormat lepton_new = *lepton_p;
    			lepton_new.ToRestFrame(z);

    			MCParticleFormat lepton2_new = *lepton_p2;
    			lepton2_new.ToRestFrame(z2);

			MALorentzVector pZZ = z->momentum() + z2->momentum();
			MALorentzVector p4l = leppt[0]->momentum() + leppt[1]->momentum() + leppt[2]->momentum() + leppt[3]->momentum();

			std::cout << "pZZ.M() = " << pZZ.M() << endl;
		        std::cout << "m4l = " << p4l.M() << endl;
			//std::cout << "mZZ = " << (z->momentum()+z2->momentum()).M() << endl;

			//zz = z + z2;
			mZZ->Fill( (z->momentum()+z2->momentum()).M(), event.mc()->weight() );

    			// Filling histogram
        
    			myHisto->Fill( cos(lepton_new.angle(z)), event.mc()->weight() );
    			myHisto->Fill( cos(lepton2_new.angle(z2)), event.mc()->weight() );
			//myHisto->Fill( cos(lepton_new.angle(z)) );
			//myHisto->Fill( cos(lepton2_new.angle(z2)) );
		}
	   }
	}
     }

  // ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ***************************************************************************
  
  /*if (event.mc()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Event weight
    double myWeight=1.;
    if (!Configuration().IsNoEventWeight()) myWeight=event.mc()->weight();

    // Initial state
    for (MAuint32 i=0;i<event.mc()->particles().size();i++)
    {
      const MCParticleFormat& part = event.mc()->particles()[i];

      cout << "----------------------------------" << endl;
      cout << "MC particle" << endl;
      cout << "----------------------------------" << endl;

      // display index particle
      cout << "index=" << i+1;

      // display the status code
      cout << "Status Code=" << part.statuscode() << endl;
      if (PHYSICS->Id->IsInitialState(part)) cout << " (Initial state) ";
      else if (PHYSICS->Id->IsFinalState(part)) cout << " (Final state) ";
      else cout << " (Intermediate state) ";
      cout << endl;

      // pdgid
      cout << "pdg id=" << part.pdgid() << endl;
      if (PHYSICS->Id->IsInvisible(part)) cout << " (invisible particle) ";
      else cout << " (visible particle) ";
      cout << endl;

      // display kinematics information
      cout << "px=" << part.px()
                << " py=" << part.py()
                << " pz=" << part.pz()
                << " e="  << part.e()
                << " m="  << part.m() << endl;
      cout << "pt=" << part.pt() 
                << " eta=" << part.eta() 
                << " phi=" << part.phi() << endl;

      // display particle mother id
      if (part.mothers().empty()) 
      {
        cout << "particle with no mother." << endl;
      }
      else
      {
        std::cout << "particle coming from the decay of";
        for(MAuint32 j=0;j<part.mothers().size();j++)
        {
          const MCParticleFormat* mother = part.mothers()[j];
          cout << " " << mother->pdgid();
        }
        std::cout << "." << endl;
      }
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.mc()->MET().pt()
         << " phi=" << event.mc()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.mc()->MHT().pt()
              << " phi=" << event.mc()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.mc()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.mc()->THT() << endl; 
   cout << endl;

  return true;
  }
  */


  // ***************************************************************************
  // Example of analysis with reconstructed objects
  // Concerned samples : 
  //   - LHCO samples
  //   - LHE/STDHEP/HEPMC samples after applying jet-clustering algorithm
  // ***************************************************************************
  /*
  // Event weight
  double myWeight=1.;
  if (!Configuration().IsNoEventWeight() && event.mc()!=0) myWeight=event.mc()->weight();

  Manager()->InitializeForNewEvent(myWeight);

  if (event.rec()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Looking through the reconstructed electron collection
    for (MAuint32 i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      cout << "----------------------------------" << endl;
      cout << "Electron" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << elec.charge() << endl;
      cout << "px=" << elec.px()
                << " py=" << elec.py()
                << " pz=" << elec.pz()
                << " e="  << elec.e()
                << " m="  << elec.m() << endl;
      cout << "pt=" << elec.pt() 
                << " eta=" << elec.eta() 
                << " phi=" << elec.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
                << elec.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed muon collection
    for (MAuint32 i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      cout << "----------------------------------" << endl;
      cout << "Muon" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << mu.charge() << endl;
      cout << "px=" << mu.px()
                << " py=" << mu.py()
                << " pz=" << mu.pz()
                << " e="  << mu.e()
                << " m="  << mu.m() << endl;
      cout << "pt=" << mu.pt() 
                << " eta=" << mu.eta() 
                << " phi=" << mu.phi() << endl;
      cout << "ET/PT isolation criterion =" << mu.ET_PT_isol() << endl;
      cout << "pointer address to the matching MC particle: " 
           << mu.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed hadronic tau collection
    for (MAuint32 i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      cout << "----------------------------------" << endl;
      cout << "Tau" << endl;
      cout << "----------------------------------" << endl;
      cout << "tau: index=" << i+1 
                << " charge=" << tau.charge() << endl;
      cout << "px=" << tau.px()
                << " py=" << tau.py()
                << " pz=" << tau.pz()
                << " e="  << tau.e()
                << " m="  << tau.m() << endl;
      cout << "pt=" << tau.pt() 
                << " eta=" << tau.eta() 
                << " phi=" << tau.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
           << tau.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed jet collection
    for (MAuint32 i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      cout << "----------------------------------" << endl;
      cout << "Jet" << endl;
      cout << "----------------------------------" << endl;
      cout << "jet: index=" << i+1 
           << " charge=" << jet.charge() << endl;
      cout << "px=" << jet.px()
           << " py=" << jet.py()
           << " pz=" << jet.pz()
           << " e="  << jet.e()
           << " m="  << jet.m() << endl;
      cout << "pt=" << jet.pt() 
           << " eta=" << jet.eta() 
           << " phi=" << jet.phi() << endl;
      cout << "b-tag=" << jet.btag()
           << " true b-tag (before eventual efficiency)=" 
           << jet.true_btag() << endl;
      cout << "EE/HE=" << jet.EEoverHE()
           << " ntracks=" << jet.ntracks() << endl;
      cout << endl;
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.rec()->MET().pt()
         << " phi=" << event.rec()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.rec()->MHT().pt()
              << " phi=" << event.rec()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.rec()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.rec()->THT() << endl;
    cout << endl;
  }
  */
  return true;
}

