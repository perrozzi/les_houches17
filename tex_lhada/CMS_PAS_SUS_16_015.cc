// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Cutflow.hh"

#include <math.h>
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
namespace Rivet {


  class CMS_PAS_SUS_16_015: public Analysis {
  public:

    ///Cut ids
    enum {kmt2_cut, kdeltaphi_etmiss_jet, kveto, kpreselection,
	  kone_jet_selection, kat_least_two_jet_selection,
	  k1j_loose_1, k1j_loose_2, k1j_medium_1, k1j_medium_2,
	  k1j_medium_3, k1j_medium_4, k1j_medium_5, k1j_medium_6,
	  k1j_medium_7, k1j_medium_8} CutIds;
    /// Constructor
    CMS_PAS_SUS_16_015():      Analysis("CMS_PAS_SUS_16_015"),
      cutflow("CutFlow", {"kmt2_cut", "kdeltaphi_etmiss_jet",
	    "kveto", "kpreselection", "kone_jet_selection",
	    "kat_least_two_jet_selection", "k1j_loose_1",
	    "k1j_loose_2", "k1j_medium_1", "k1j_medium_2",
	    "k1j_medium_3", "k1j_medium_4", "k1j_medium_5",
	    "k1j_medium_6", "k1j_medium_7", "k1j_medium_8"})
    {  }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      VisibleFinalState visfs(fs);

      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "jets_eta47");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "bjets");
      addProjection(MissingMomentum(fs), "met");
    }

    template<typename P>
    double  scalar_pt_sum(const std::vector<P>& momenta){
        double ht  = 0;
        for(const auto& p: momenta){
          ht += p.pt();
        }
        return ht;
    }


    double  btag_eff(const Particles& bjets){
        return 1.;
    }


    template<typename P>
    double  mt2(const P& particle1, const P& particle2, const P& met){
        return P();
    }


    template<typename P1, typename P2>
    double  dphi(const P1& p1, const P2& p2){
        double r = acos(p1.px()*p2.px() + p1.px()+p2.py() )
          / sqrt((p1.px()*p1.px() + p1.py()*p1.py())
		 * (p1.px()*p1.px() + p1.py()*p1.py()));
        return isnan(r) ? 0 : r;
    }


    bool cut_mt2_cut(){
        bool r = true;
        return cutflow.fill(kmt2_cut, r);
    }

    bool cut_deltaphi_etmiss_jet(){
        bool r = true;
        r &= (dphi(met, jets[1 - 1]) > 0.3)  && (dphi(met, jets[2 - 1]) > 0.3)
	  && (dphi(met, jets[3 - 1]) > 0.3)  && (dphi(met, jets[4 - 1]) > 0.3);
        return cutflow.fill(kdeltaphi_etmiss_jet, r);
    }

    bool cut_veto(){
        bool r = true;
        return cutflow.fill(kveto, r);
    }

    bool cut_preselection(){
        bool r = true;
        r &= cut_mt2_cut();
        r &= cut_deltaphi_etmiss_jet();
        r &= cut_veto();
        return cutflow.fill(kpreselection, r);
    }

    bool cut_one_jet_selection(){
        bool r = true;
        r &= cut_preselection();
        r &= jets.size() == 1;
        r &= jets[1 - 1].pt() > 200;
        return cutflow.fill(kone_jet_selection, r);
    }

    bool cut_at_least_two_jet_selection(){
        bool r = true;
        r &= cut_preselection();
        r &= jets.size()> 1;
        r &= cut_mt2_cut()> 200;
        return cutflow.fill(kat_least_two_jet_selection, r);
    }

    bool cut_1j_loose_1(){
        bool r = true;
        r &= cut_one_jet_selection();
        r &= ht> 575;
        return cutflow.fill(k1j_loose_1, r);
    }

    bool cut_1j_loose_2(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= bjets.size() < 3;
        r &= ht> 575;
        r &= ht < 1000;
        r &= cut_mt2_cut()> 200;
        return cutflow.fill(k1j_loose_2, r);
    }

    bool cut_1j_medium_1(){
        bool r = true;
        r &= cut_one_jet_selection();
        r &= bjets.size() < 1;
        r &= ht> 1000;
        return cutflow.fill(k1j_medium_1, r);
    }

    bool cut_1j_medium_2(){
        bool r = true;
        r &= cut_one_jet_selection();
        r &= bjets.size()> 0;
        r &= ht> 575;
        return cutflow.fill(k1j_medium_2, r);
    }

    bool cut_1j_medium_3(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size() < 1;
        r &= ht> 575;
        r &= ht < 1000;
        r &= cut_mt2_cut()> 800;
        return cutflow.fill(k1j_medium_3, r);
    }

    bool cut_1j_medium_4(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size()> 0;
        r &= bjets.size() < 3;
        r &= ht> 575;
        r &= ht < 1000;
        r &= cut_mt2_cut()> 600;
        return cutflow.fill(k1j_medium_4, r);
    }

    bool cut_1j_medium_5(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size() < 2;
        r &= ht> 1000;
        r &= ht < 1500;
        r &= cut_mt2_cut()> 800;
        return cutflow.fill(k1j_medium_5, r);
    }

    bool cut_1j_medium_6(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size() == 2;
        r &= ht> 1000;
        r &= ht < 1500;
        r &= cut_mt2_cut()> 400;
        return cutflow.fill(k1j_medium_6, r);
    }

    bool cut_1j_medium_7(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size() < 2;
        r &= ht> 1500;
        r &= cut_mt2_cut()> 400;
        return cutflow.fill(k1j_medium_7, r);
    }

    bool cut_1j_medium_8(){
        bool r = true;
        r &= cut_at_least_two_jet_selection();
        r &= jets.size() < 4;
        r &= bjets.size() == 2;
        r &= ht> 1500;
        r &= cut_mt2_cut()> 200;
        return cutflow.fill(k1j_medium_8, r);
    }
    /// Perform the per-event analysis
    void analyze(const Event& event) {

       const FastJets& jets_eta47Proj = applyProjection<FastJets>(event, "jets_eta47");
       Jets jets_eta47 = jets_eta47Proj.jetsByPt();
       Jets jets;
       for(const auto& p: jets_eta47){
         if((p.eta() < 2.4)){
           jets.push_back(p);
         }
       }

       const FastJets& bjetsProj = applyProjection<FastJets>(event, "bjets");
       Jets bjets = bjetsProj.jetsByPt();
       const MissingMomentum metProj = applyProjection<MissingMomentum>(event, "met");
       met = metProj.missingMomentum();
       ht = scalar_pt_sum(jets);

       cut_1j_loose_1();
       cut_1j_loose_2();
       cut_1j_medium_1();
       cut_1j_medium_2();
       cut_1j_medium_3();
       cut_1j_medium_4();
       cut_1j_medium_5();
       cut_1j_medium_6();
       cut_1j_medium_7();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       std::cout << "Analsyis cut flow:\n"
                 << "-----------------\n\n"
                 << cutflow << "\n";
    }


  protected:

    //@{
    /** Collections and variables
     */
    //@}
     /** Analysis objects
      * @{
      */
     Jets jets_eta47;

     Jets jets;

     Jets bjets;

     FourMomentum met;

     double  ht;
     /** @}
     */    //@{
    /** Histograms
     */
    //@}

    ///Tracks the event counts after each cut
    Cutflow cutflow;

  };
  DECLARE_RIVET_PLUGIN(CMS_PAS_SUS_16_015);
}



