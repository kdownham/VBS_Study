// -*- C++ -*- 

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"

#include <iostream>
#include <algorithm>
#include <string>

static const float mZ = 91.1876;

namespace Rivet {

class VBS_selection: public Analysis {
public:
    enum CutFlow {
        inSample,
        passingLepAcceptance,
        passingLepVeto,
        passingM2lconstraint,
	passing2jselection,
	passingAll
    };
    float totalEvents = 0;
    float sumSelectedWeights = 0;
    float common_weight = 0;
    float w_1 = 0;
    float w_2 = 0;
    float w_3 = 0;
    float w_4 = 0;
    float w_5 = 0;
    float w_6 = 0;
    float w_7 = 0;
    float w_8 = 0;
    float sumWeightsLeps = 0; 
    float sumWeights2j = 0;
    float sumWeightsCR = 0;
    float sumWeightsVeto = 0;
    float selectedEvents = 0;
    float veto1 = 0;
    float veto2 = 0;
    float veto3 = 0;
    float veto4 = 0;
    float veto5 = 0; 
    float veto6 = 0;
    float sameflavor = 0;
    float dflavor = 0;

    VBS_selection() : Analysis("VBS_selection"){};

    void init() {
	double lepConeSize = 0.4;
        double lepMaxEta = 2.5;

	Cut lepton_cut = (Cuts::abseta < lepMaxEta);

	FinalState fs(-2.5,2.5,0.0*GeV);
        FinalState fsm(-5,5,0.0*GeV);
        addProjection(fs, "FS");
        addProjection(fsm, "FSM");

        ChargedLeptons charged_leptons(fs);
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);

        PromptFinalState prompt_leptons(charged_leptons);
        prompt_leptons.acceptMuonDecays(true);
        prompt_leptons.acceptTauDecays(false);

        PromptFinalState prompt_photons(photons);
        prompt_photons.acceptMuonDecays(true);
        prompt_photons.acceptTauDecays(false);

 	Cut cut_el = Cuts::abseta < 2.5 && Cuts::pT > 7.0*GeV;
	Cut cut_mu = Cuts::abseta < 2.4 && Cuts::pT > 5.0*GeV;


	DressedLeptons dressed_leptons = DressedLeptons(prompt_photons, prompt_leptons, lepConeSize, lepton_cut, true, true);
	addProjection(dressed_leptons, "DressedLeptons");

	addProjection(FastJets(FinalState(-4.7, 4.7), FastJets::ANTIKT, 0.4),"jets");

	addProjection(MissingMomentum(fs), "MET");

	IdentifiedFinalState nu_id;
	nu_id.acceptNeutrinos();
	PromptFinalState neutrinos(nu_id);
	neutrinos.acceptTauDecays(false);
	declare(neutrinos,"neutrinos");

	//ZFinder zeefinder(FinalState(), cut_el, PID::ELECTRON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
	//declare(zeefinder, "ZeeFinder");
	//ZFinder zmmfinder(FinalState(), cut_mu, PID::MUON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
	//declare(zmmfinder, "ZmmFinder");
	
	//WFinder wefinder(FinalState(), cut_el, PID::ELECTRON, 0*GeV, 1000000*GeV, 20*GeV, 0.1);

	// WFinder::ChargedLeptons::PROMPT, WFinder::CLUSTERNODECAY, WFinder::AddPhotons::YES, WFinder::MassWindow::M, 80.4*GeV);
	//declare(wefinder, "WeFinder");

	//WFinder wmfinder(FinalState(), cut_mu, PID::MUON, 0*GeV, 1000000*GeV, 20*GeV, 0.1);
	//declare(wmfinder, "WmFinder");
	
	//VBS variables
	hist_mjj = bookHisto1D("mjj", 6, 500, 3500);
	hist_dEtajj = bookHisto1D("dEtajj", 32, 0, 8);
	hist_mjj_1 = bookHisto1D("mjj_1", 6, 500, 3500);
	hist_dEtajj_1 = bookHisto1D("dEtajj_1", 32, 0, 8);
	hist_mjj_c = bookHisto1D("mjj_c", 7, 0, 3500);
	hist_dEtajj_c = bookHisto1D("dEtajj_c", 32, 0, 8);
	//hist_dRjj = bookHisto1D("dRjj", 60, 0, 15);

	//Cut variables
	hist_pt_l1 = bookHisto1D("leadinglpt", 100, 0, 300);
	hist_pt_l2 = bookHisto1D("secondlpt", 100, 0, 300);
	hist_MET = bookHisto1D("MET", 100, 0, 300);
	hist_ptll = bookHisto1D("ptll", 10, 0, 100);
	hist_mll = bookHisto1D("mll", 50, 0, 300);
	hist_ptj1 = bookHisto1D("jet1_Pt", 100, 0, 1000);
	hist_ptj2 = bookHisto1D("jet2_Pt", 100, 0, 1000);

	hist_ct_Z1 = bookHisto1D("cos_theta_Z1",100,-1.0,1.0);

     }

    void analyze(const Event& event) {
	double weight = event.weight();
	//std::cout << "Event weight = " << weight << endl;
	totalEvents++;

	//const ZFinder& zeefinder = applyProjection<ZFinder>(event, "ZeeFinder");
	//const ZFinder& zmmfinder = applyProjection<ZFinder>(event, "ZmmFinder");
	//const WFinder& wefinder = applyProjection<WFinder>(event, "WeFinder");
	//const WFinder& wmfinder = applyProjection<WFinder>(event, "WmFinder");
	
	Particles leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").particlesByPt(10*GeV);
	if(leptons.size() < 2) {
	   vetoEvent;
	}
	if (leptons.at(0).pT() < 20*GeV) vetoEvent;

	int sumPids = leptons.at(0).pdgId() + leptons.at(1).pdgId();
	//int sumAbsPids = std::abs(leptons.at(0).pdgId()) + std::abs(leptons.at(1).pdgId());
	const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
	FourMomentum leps = leptons.at(0).momentum() + leptons.at(1).momentum();
	float mll = leps.mass();

	
	Jets jets;
        foreach (const Jet& jet, applyProjection<FastJets>(event, "jets").jetsByPt(30*GeV) ) {
            bool isolated = true;
            foreach (const Particle& lepton, leptons) {
                if (deltaR(lepton, jet) < 0.4) {
                    isolated = false;
                    break;
                }
            }
            if (isolated){
                jets.push_back(jet);
		}
        }

        if (jets.size() < 2) {  
            vetoEvent;
        }	

	FourMomentum dijet_system = jets.at(0).momentum() + jets.at(1).momentum();
	float mjj = dijet_system.mass();
	float dEtajj = jets.at(0).momentum().eta() - jets.at(1).momentum().eta();


	//Event Selection
	if (mjj < 500*GeV || std::abs(dEtajj) < 2.5) vetoEvent;

	common_weight += weight;
	
	//Process-specific cuts
	
	//Both true corresponds to ZZ, s_1 = true and s_2 = false corresponds to WZ
	//s_1 = false and s_2 = true corresponds to same-sign WW and both false 
	//corresponds to opposite-sign WW	    

	bool s_1 = false;
	bool s_2 = false;

	int Njets = jets.size();

	//FourMomentum pZ_a, pZ_b, pZ_1, pZ_2;
	//FourMomentum Z_a_l1, Z_a_l2, Z_b_l1, Z_b_l2, pZZ;
	
	FourMomentum Z_l1_l2, Z_l1_l3, Z_l1_l4, Z_l2_l3, Z_l2_l4, Z_l3_l4;
	double mZ_l1_l2 = 0;
	double mZ_l1_l3 = 0;
	double mZ_l1_l4 = 0;
	double mZ_l2_l3 = 0;
	double mZ_l2_l4 = 0;
	double mZ_l3_l4 = 0;
	
	if (s_1 && s_2){
	// ZZ selection
	
	   //std::cout << "Event has reached ZZ selection" << endl;
	   if (leptons.size() < 4) vetoEvent;
	   if (leptons.size() > 4) vetoEvent;

	   double mdiff_12 = 0;
	   double mdiff_13 = 0;
	   double mdiff_14 = 0;
	   double mdiff_23 = 0;
	   double mdiff_24 = 0;
	   double mdiff_34 = 0;
	   double m_a = 0;
	   double m_b = 0;
	   double m_c = 0;
	   double m_d = 0;
	   double m_e = 0;
	   double m_f = 0;
	   double a = 99999999;
	   double b = 99999999;
	   double c = 99999999;
	   double d = 99999999;
	   double e = 99999999;
	   double f = 99999999;

	   double mZ_1 = 0;
	   double mZ_2 = 0;

	   //int z_cand = 0;

	   //std::cout << "Event has passed number of leptons veto" << endl;

	   if (leptons.at(0).pdgId() + leptons.at(1).pdgId() == 0 && leptons.at(2).pdgId() + leptons.at(3).pdgId() == 0){
		Z_l1_l2 = leptons.at(0).momentum() + leptons.at(1).momentum();
		Z_l3_l4 = leptons.at(2).momentum() + leptons.at(3).momentum();
		mZ_l1_l2 = Z_l1_l2.mass();
		mZ_l3_l4 = Z_l3_l4.mass();
		mdiff_12 = std::abs(mZ_l1_l2 - 91.2);
		mdiff_34 = std::abs(mZ_l3_l4 - 91.2);
		if (mdiff_12 < mdiff_34){
		    a = mdiff_12;
		    m_a = mZ_l1_l2;
		    b = mdiff_34;
		    m_b = mZ_l3_l4;
	       }else{
		    a = mdiff_34;
		    m_a = mZ_l3_l4;
		    b = mdiff_12;
		    m_b = mZ_l1_l2;
		}
		//z_cand++;
	   }

	   if (leptons.at(0).pdgId() + leptons.at(2).pdgId() == 0 && leptons.at(1).pdgId() + leptons.at(3).pdgId() == 0){
		Z_l1_l3 = leptons.at(0).momentum() + leptons.at(2).momentum();
		Z_l2_l4 = leptons.at(1).momentum() + leptons.at(3).momentum();
		mZ_l1_l3 = Z_l1_l3.mass();
		mZ_l2_l4 = Z_l2_l4.mass();
		mdiff_13 = std::abs(mZ_l1_l3 - 91.2);
		mdiff_24 = std::abs(mZ_l2_l4 - 91.2);
		//std::cout << "l1 l3 mass difference = " << mdiff_13 << endl;
		//std::cout << "l2 l4 mass difference = " << mdiff_24 << endl;
		if (mdiff_13 < mdiff_24){
		    c = mdiff_13;
		    m_c = mZ_l1_l3;
		    d = mdiff_24;
		    m_d = mZ_l2_l4;
	       }else{
		    c = mdiff_24;
		    m_c = mZ_l2_l4;
		    d = mdiff_13;
		    m_d = mZ_l1_l3;
	        }
		//std::cout << "smallest mass difference from Z mass = " << c << endl;
		//std::cout << "largest mass difference from Z mass = " << d << endl;
	   }


	   if (leptons.at(0).pdgId() + leptons.at(3).pdgId() == 0 && leptons.at(1).pdgId() + leptons.at(2).pdgId() == 0){
	        Z_l1_l4 = leptons.at(0).momentum() + leptons.at(3).momentum();
		Z_l2_l3 = leptons.at(1).momentum() + leptons.at(2).momentum();
		mZ_l1_l4 = Z_l1_l4.mass();
		mZ_l2_l3 = Z_l2_l3.mass();
		mdiff_14 = std::abs(mZ_l1_l4 - 91.2);
		mdiff_23 = std::abs(mZ_l2_l3 - 91.2);
		if (mdiff_14 < mdiff_23){
		    e = mdiff_14;
		    m_e = mZ_l1_l4;
		    f = mdiff_23;
		    m_f = mZ_l2_l3;
	       }else{
		    e = mdiff_23;
		    m_e = mZ_l2_l3;
		    f = mdiff_14;
		    m_f = mZ_l1_l4;
		}
	   }	
	   
	   if (a < c && a < e){
		mZ_1 = m_a;
		mZ_2 = m_b;
		std::cout << "a = " << a << endl;
		std::cout << "b = " << b << endl;
	  }else if (c < a && c < e){
		mZ_1 = m_c;
		mZ_2 = m_d;
		std::cout << "c = " << c << endl;
		std::cout << "d = " << d << endl;
	  }else if (e < a && e < c){
		mZ_1 = m_e;
		mZ_2 = m_f;
		std::cout << "e = " << e << endl;
		std::cout << "f = " << f << endl;
	  }

	  if (mZ_1 < 60.0 || mZ_1 > 120.0) vetoEvent;
	  if (mZ_2 < 60.0 || mZ_2 > 120.0) vetoEvent;
	
	   //if (zeefinder.bosons().size() < 1 && zmmfinder.bosons().size() < 1) vetoEvent;
	
	   /*
	   if (zeefinder.bosons().size() + zmmfinder.bosons().size() < 2) vetoEvent;

	   std::cout << "Event has 2 or more Z bosons" << endl;

	   if (zeefinder.bosons().size() + zmmfinder.bosons().size() > 2) vetoEvent;

	   std::cout << "Event has exactly 2 Z bosons" << endl;

	   if (zmmfinder.bosons().size() > 1 && zeefinder.bosons().size() < 1){
		pZ_a = zmmfinder.bosons()[0];
	 	pZ_b = zmmfinder.bosons()[1];
	 	pZ_1 = pZ_a;
		pZ_2 = pZ_b;
		pZZ = pZ_a + pZ_b;
		Z_a_l1 = zmmfinder.constituents()[0];
		Z_a_l2 = zmmfinder.constituents()[1];
		//Z_a_l1_id = zmmfinder.constituents()[0].pdgId();
		//Z_a_l2_id = zmmfinder.constituents()[1].pdgId();
	   } else if (zmmfinder.bosons().size() < 1 && zeefinder.bosons().size() > 1){
		pZ_a = zeefinder.bosons()[0];
                pZ_b = zeefinder.bosons()[1];
                pZ_1 = pZ_a;
                pZ_2 = pZ_b;
                pZZ = pZ_a + pZ_b;
                Z_a_l1 = zeefinder.constituents()[0];
                Z_a_l2 = zeefinder.constituents()[1];
                //Z_a_l1_id = zeefinder.constituents()[0].pdgId();
                //Z_a_l2_id = zeefinder.constituents()[1].pdgId();
	   } else if (zmmfinder.bosons().size() > 0 && zeefinder.bosons().size() > 0){
		pZ_a = zmmfinder.bosons()[0];
                pZ_b = zeefinder.bosons()[0];
                pZ_1 = pZ_a;
                pZ_2 = pZ_b;
                pZZ = pZ_a + pZ_b;
                Z_a_l1 = zmmfinder.constituents()[0];
                Z_a_l2 = zmmfinder.constituents()[1];
                //Z_a_l1_id = zmmfinder.constituents()[0].pdgId();
                //Z_a_l2_id = zmmfinder.constituents()[1].pdgId();	
		}

	   FourMomentum Z1_l1 = pZ_1 + Z_a_l1;
	   
	   double pZpl = pZ_1.px()*Z_a_l1.px() + pZ_1.py()*Z_a_l1.py() + pZ_1.pz()*Z_a_l1.pz();
	   double pZabs = sqrt(pZ_1.px2() + pZ_1.py2() + pZ_1.pz2());
	   double plabs = sqrt(Z_a_l1.px2() + Z_a_l1.py2() + Z_a_l1.pz2());

	   double ct_Z1_l1 = (pZpl - (pZ_1.E()*Z_a_l1.E())) / (pZabs*plabs);	

	   //const double ct_Z1_l1 = 2 * pZ_1.pT() * Z_a_l1.pT() * sinh(pZ_1.eta() - Z_a_l1.eta()) / Z1_l1.mass() / add_quad(Z1_l1.mass(), Z1_l1.pT());

	   hist_ct_Z1->fill(ct_Z1_l1, weight);
	   */
	   
		    	  
	}else if (s_1 && !s_2){
	// WZ selection
		
	   //std::cout << "Congratulations! You have chosen the WZ selection!" << endl;	

	   FourMomentum Z_12, Z_13, Z_23;

	   double m_12 = 9999999999;
	   double m_13 = 9999999999;
	   double m_23 = 9999999999;
	   double mdiff_1_2 = 0;
	   double mdiff_1_3 = 0;
	   double mdiff_2_3 = 0; 

	   if (leptons.size() < 3) vetoEvent;
	   if (leptons.size() > 3) vetoEvent;

	   // Find the Z candidate
	   // Check if the leading leptons are a Z candidate
	   if (leptons.at(0).pdgId() + leptons.at(1).pdgId() == 0){
		// Form the dilepton momentum (Z candidate momentum)
		Z_12 = leptons.at(0).momentum() + leptons.at(1).momentum();
		// Store the invariant mass of the Z candidate
		m_12 = Z_12.mass();
		// Determine the difference between the Z candidate mass and the nominal Z mass
		mdiff_1_2 = std::abs(m_12 - 91.2);
		}
	   
	   // Check if the leading and third lepton are a Z candidate
	   if (leptons.at(0).pdgId() + leptons.at(2).pdgId() == 0){
		Z_13 = leptons.at(0).momentum() + leptons.at(2).momentum();
		m_13 = Z_13.mass();
		mdiff_1_3 = std::abs(m_13 - 91.2);
		}

	   // Check if the second and third leptons are a Z candidate
	   if (leptons.at(1).pdgId() + leptons.at(2).pdgId() == 0){
		Z_23 = leptons.at(1).momentum() + leptons.at(2).momentum();
		m_23 = Z_23.mass();
		mdiff_2_3 = std::abs(m_23 - 91.2);
		}

	  
	   //Find the Z candidate with mass closest to nominal Z mass
	   if (mdiff_1_2 < mdiff_1_3 && mdiff_1_2 < mdiff_2_3){
		
		// W lepton pT veto
		if (leptons.at(2).pT() < 20) vetoEvent;
		
		// Require Z candidate to be within mass window
		if (m_12 < 60.0 || m_12 > 120.0) vetoEvent;		

	   }else if (mdiff_1_3 < mdiff_1_2 && mdiff_1_3 < mdiff_2_3){

		if (leptons.at(1).pT() < 20) vetoEvent;
		
		if (m_13 < 60.0 || m_13 > 120.0) vetoEvent;

	   }else if (mdiff_2_3 < mdiff_1_2 && mdiff_2_3 < mdiff_1_3){

		if (m_23 < 60.0 || m_23 > 120.0) vetoEvent;

	   }
 

	}else if (!s_1 && s_2){
	// same-sign WW selection
	   if (leptons.size() < 2) vetoEvent;
	   if (leptons.size() > 2) vetoEvent;
	
	   // Opposite sign lepton veto	
	   if (leptons.at(0).pdgId()*leptons.at(1).pdgId() < 0) vetoEvent;

	   if (leptons.at(1).pT() < 20) vetoEvent;

	}else if (!s_1 && !s_2){
	// opposite-sign WW selection
	   if (leptons.size() < 2) vetoEvent;
	   if (leptons.size() > 2) vetoEvent;

	   if (leptons.at(0).pdgId()*leptons.at(1).pdgId() > 0) vetoEvent;

	   if (leptons.at(0).pdgId() + leptons.at(1).pdgId() == 0) vetoEvent;
	   
	   if (leptons.at(1).pT() < 20) vetoEvent;

	   for (int i=0; i <= Njets-1; i++){
		if (jets.at(i).bTagged()){
		    weight = weight*0.15;
			}
		}  
	
	}


		sumSelectedWeights += weight;
		selectedEvents++;
		hist_mjj->fill(mjj/GeV, weight);
		hist_mjj_1->fill(mjj/GeV, weight);
		hist_dEtajj->fill(std::abs(dEtajj), weight);
		hist_dEtajj_1->fill(std::abs(dEtajj), weight);
		hist_pt_l1->fill(leptons.at(0).pT()/GeV, weight);
		hist_pt_l2->fill(leptons.at(1).pT()/GeV, weight);
		hist_ptll->fill(leps.pT()/GeV, weight);
		hist_ptj1->fill(jets.at(0).pT()/GeV, weight);
		hist_ptj2->fill(jets.at(1).pT()/GeV, weight);
		//hist_ct_Z1->fill(ct_Z1_l1, weight);

    }

    void finalize() {
	float efficiency = selectedEvents/totalEvents;
	double xsec = crossSection();

	std::cout << "SumOfWeights() processed = " << sumOfWeights() << endl;
	std::cout << "Selected sumWeights = " << sumSelectedWeights << endl;

	std::cout << "Events with same flavor leptons: " << sameflavor << endl;
	std::cout << "Events with different flavor leptons: " << dflavor << endl;
	std::cout << "Events after veto 1: " << veto1 << endl;
	std::cout << "Events after veto 2: " << veto2 << endl;
	std::cout << "Events after veto 3: " << veto3 << endl;
	std::cout << "Events after veto 4: " << veto4 << endl;
	std::cout << "Events after veto 5: " << veto5 << endl;
	std::cout << "Events after veto 6: " << veto6 << endl;
	
	std::cout << std::endl << "Selected Events = " << selectedEvents << ", Total = " << totalEvents << endl;
	std::cout << "Efficiency = " << efficiency << endl;

	std::cout << "Initial cross section was: " << xsec << std::endl;
	std::cout << "Fiducial cross section is: " << xsec*sumSelectedWeights/sumOfWeights() << std::endl;

	std::cout << "Event Yield: " << 35.9*xsec*sumSelectedWeights/sumOfWeights() << std::endl;

	float sumWeights = sumOfWeights();
	scale(hist_mjj, xsec/sumWeights);
	scale(hist_mjj_1, xsec*35.9/sumWeights);
	scale(hist_mjj_c, xsec/sumWeights);
	scale(hist_dEtajj, xsec/sumWeights);
	scale(hist_dEtajj_1, xsec*35.9/sumWeights);
	scale(hist_dEtajj_c, xsec/sumWeights);
	scale(hist_ptj1, xsec/sumWeights);
	scale(hist_ptj2, xsec/sumWeights);

	}

	//Histograms
	Histo1DPtr hist_mjj, hist_dEtajj, hist_pt_l1, hist_pt_l2, hist_MET, hist_ptll, hist_mll, hist_ptj1, hist_ptj2, hist_mjj_1, hist_mjj_c, hist_dEtajj_1, hist_dEtajj_c, hist_ct_Z1;

     };

  DECLARE_RIVET_PLUGIN(VBS_selection);

}
















