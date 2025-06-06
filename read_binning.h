#ifndef READ_BINNING_H
#define READ_BINNING_H
#include "TEnv.h"

class read_binning
{
public:
  read_binning(const std::string configfile)
    {
      penv = new TEnv(configfile.c_str());

    }

  Int_t get_nbins(){ return penv->GetValue("nbins", 1); }
  Int_t get_bbins(){ return penv->GetValue("bbins", 1); }
  Int_t get_minentries(){ return penv->GetValue("minentries", 1); }
  Int_t get_measure_bins(){ return penv->GetValue("measure_bins", 1); }
  Double_t get_minimum(){ return penv->GetValue("minimum", 1.0); }
  Double_t get_fixed(){ return penv->GetValue("fixed", 1.0); }
  Int_t get_primer(){ return penv->GetValue("primer", 0); }
  Int_t get_njet_sys(){ return penv->GetValue("NJET", 0); }
  Int_t get_prior_sys(){ return penv->GetValue("PRIOR", 0); }
  Int_t get_vtx_sys(){ return penv->GetValue("VTX", 0); }
  Double_t get_jes_sys(){ return penv->GetValue("JES", 0.0); }
  Double_t get_jer_sys(){ return penv->GetValue("JER", 0.0); }
  Double_t get_first_xj(){ return penv->GetValue("first_xj", 0.0); }
  
  Double_t get_pt_regions(int reg){ return penv->GetValue(Form("pt1_%d", reg), 1.0); }
  Double_t get_pt2_regions(int reg){ return penv->GetValue(Form("pt2_%d", reg), 1.0); }
  
  Double_t get_low_trigger(int reg){ return low_trigger[reg];}
  Double_t get_low_trigger_goal() { return penv->GetValue("low_trigger", 1.0);}
  
  Int_t get_measure_region(int reg){ return measure_region[reg]; }
  Int_t get_subleading_measure_region(int reg){ return subleading_measure_region[reg]; }

  Double_t get_sample_boundary(int ib){ return sample_boundary[ib]; }

  Double_t get_dphicut() { return penv->GetValue("dphi_top", 1.0) * TMath::Pi() / penv->GetValue("dphi_bottom", 2.0); }
  Double_t get_truth_leading_goal(){ return penv->GetValue("truth_leading_goal", 14.0); }
  Double_t get_truth_subleading_goal(){ return penv->GetValue("truth_subleading_goal", 14.0); }

  Int_t get_reco_leading_bin_buffer(){ return penv->GetValue("reco_leading_bin_buffer", 14.0); }
  Int_t get_reco_subleading_bin_buffer(){ return penv->GetValue("reco_subleading_bin_buffer", 14.0); }

  Int_t get_measure_leading_bin_buffer(){ return penv->GetValue("measure_leading_bin_buffer", 14.0); }
  Int_t get_measure_subleading_bin_buffer(){ return penv->GetValue("measure_subleading_bin_buffer", 14.0); }


  int get_truth_leading_bin(){ return truth_leading_bin;}
  int get_truth_subleading_bin(){ return truth_subleading_bin;}

  int get_reco_leading_bin(){ return reco_leading_bin;}
  int get_reco_subleading_bin(){ return reco_subleading_bin;}

  int get_measure_leading_bin(){ return measure_leading_bin;}
  int get_measure_subleading_bin(){ return measure_subleading_bin;}

  
  float get_truth_leading_cut(){ return truth_leading_cut;}
  float get_truth_subleading_cut(){ return truth_subleading_cut;}

  float get_reco_leading_cut(){ return reco_leading_cut;}
  float get_reco_subleading_cut(){ return reco_subleading_cut;}

  float get_measure_leading_cut(){ return measure_leading_cut;}
  float get_measure_subleading_cut(){ return measure_subleading_cut;}

  
  void get_pt_bins(float ipt_bins[])
  {
    Float_t minimum = get_minimum();
    Float_t fixed = get_fixed();
    Int_t beforebin = get_bbins();
    Float_t low_trigger_goal = get_low_trigger_goal();
    
    float alpha = TMath::Power(fixed/minimum, 1/(float)beforebin);

    Float_t truth_leading_goal = get_truth_leading_goal();
    Float_t truth_subleading_goal = get_truth_subleading_goal();
    Int_t nbins = get_nbins();
    int ireg = 0;
    int ireg2 = 0;
    for (int i = 0; i < nbins+1; i++)
      {
	float ipt = minimum*TMath::Power(alpha, (float)i);
	ipt_bins[i] = ipt;
	for (int ib = 0; ib < 4; ib++)
	  {
	    if (sample_boundary[ib] == 0 && ipt >= sample_boundary_goal[ib])
	      {
		sample_boundary[ib] = ipt;
	      }
	  }
	if (low_trigger_goal <= ipt && low_trigger_bin == 0) low_trigger_bin = i;
	if (ireg2 <= get_measure_bins() && subleading_measure_region[ireg2] == 0 && ipt >= get_pt2_regions(ireg2)) subleading_measure_region[ireg2++] = i;
	if (ireg <= get_measure_bins() && measure_region[ireg] == 0 && ipt >= get_pt_regions(ireg)) measure_region[ireg++] = i;
	if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
      }

    low_trigger[0] = ipt_bins[low_trigger_bin];
    low_trigger[1] = ipt_bins[low_trigger_bin+1];
    low_trigger[2] = ipt_bins[low_trigger_bin+2];

    truth_subleading_bin = 0;
    
    reco_leading_bin = truth_leading_bin + get_reco_leading_bin_buffer();
    measure_leading_bin = truth_leading_bin + get_measure_leading_bin_buffer();
    reco_subleading_bin = truth_subleading_bin + get_reco_subleading_bin_buffer();
    measure_subleading_bin = truth_subleading_bin + get_measure_subleading_bin_buffer();

    truth_leading_cut = ipt_bins[truth_leading_bin];
    reco_leading_cut = ipt_bins[reco_leading_bin];
    measure_leading_cut = ipt_bins[measure_leading_bin];
    truth_subleading_cut = ipt_bins[truth_subleading_bin];
    reco_subleading_cut = ipt_bins[reco_subleading_bin];
    measure_subleading_cut = ipt_bins[measure_subleading_bin];


    return;
  }

  void get_xj_bins(float ixj_bins[])
  {
    Int_t nbins = get_nbins();
    Float_t minimum = get_minimum();
    Float_t fixed = get_fixed();
    Int_t beforebin = get_bbins();

    float alpha = TMath::Power(fixed/minimum, 1/(float)beforebin);
    float maximum = minimum*TMath::Power(alpha, nbins );
    float bin_xj10 = 1.0;
    float bin_xj1 = 1.0*(minimum/maximum);

    for (int i = 0; i < nbins+1; i++)
      {
	float ixj = bin_xj1*TMath::Power(alpha, (float)i);
	ixj_bins[i] = ixj;
      }
    return;
  }
 private:
  TEnv *penv{nullptr};

  float sample_boundary_goal[4] = {13.99, 19.99, 29.99, 100};

  float sample_boundary[4] = {0, 0, 0, 100};

  float low_trigger[3] = {0};
  int low_trigger_bin = 0;
  int measure_region[10] = {0};
  int subleading_measure_region[10] = {0};

  float truth_leading_cut = 0;
  float truth_subleading_cut = 0;

  float reco_leading_cut = 0;
  float reco_subleading_cut = 0;

  float measure_leading_cut = 0;
  float measure_subleading_cut = 0;

  
  int truth_leading_bin = 0;
  int truth_subleading_bin = 0;

  int reco_leading_bin = 0;
  int reco_subleading_bin = 0;

  int measure_leading_bin = 0;
  int measure_subleading_bin = 0;

  
};
#endif
