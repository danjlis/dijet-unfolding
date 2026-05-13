#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "dijetfinder.h"

#include "TProfile.h"
#include "TProfile2D.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

static int verbosity = 0;

const int separate_efficiency_correction = 0;

int createResponse_noempty_pp(const std::string configfile = "binning.config", const int full_or_half = 0, const int niterations = 10, const int cone_size = 4,  const int primer = 0, const int useFakes = 1, const int useEfficiencies = 1, const int input_generator = 0, const int unfold_generator = 0, const int get_mapping = 0, const int do_trimming = 1)
{
  gRandom->SetSeed(0);
  gStyle->SetPalette(kRainBow);  
  std::cout << "Using fakes : " << useFakes << std::endl;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  // pt boundaries for samples
  const int nsamples = 7;

  float boundary_r4[nsamples+1];
  boundary_r4[0] = 8;
  boundary_r4[1] = 14;
  boundary_r4[2] = 22;
  boundary_r4[3] = 33;
  boundary_r4[4] = 42;
  boundary_r4[5] = 52;
  boundary_r4[6] = 62;
  boundary_r4[7] = 100;
  float boundary_pileup_r4[nsamples+1];
  boundary_pileup_r4[0] = 8;
  boundary_pileup_r4[1] = 14;
  boundary_pileup_r4[2] = 22;
  boundary_pileup_r4[3] = 33;
  boundary_pileup_r4[4] = 42;
  boundary_pileup_r4[5] = 52;
  boundary_pileup_r4[6] = 100;
  boundary_pileup_r4[7] = 100;

  // read binning from config file
  read_binning rb(configfile.c_str());

  // dijet cuts are found here
  dijetfinder djf(cone_size);
  djf.SetVerbosity(verbosity);
  
  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();

  Double_t vtx_cut = rb.get_vtx_cut();
  Double_t njet_cut = rb.get_njet_cut();

  Int_t zyam_sys = rb.get_zyam_sys();
  Int_t inclusive_sys = rb.get_inclusive_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();

  Double_t pileup_sys = rb.get_pileup_sys();

  Int_t prior_sys = rb.get_prior_sys();
  Int_t full_sys = rb.get_full_sys();
  Int_t emfrac_sys = rb.get_emfrac_sys();
  Int_t use_herwig = rb.get_herwig();
  Int_t trigger_sys = rb.get_trigger_sys();
  Int_t ca_sys = rb.get_crossingangle_sys();
  Int_t philoc_sys = rb.get_philoc_sys();
    
  TF1 *gjer = nullptr;
			 
  TH1D *hjersmear = nullptr;
  std::string sim_string = "PYTHIA";
  if (use_herwig) sim_string = "HERWIG";
  if (pileup_sys > 1) sim_string = "PILEUP";
  else if (pileup_sys > 0) sim_string = "MIX";

  std::cout << "Getting "<< sim_string << std::endl;
  
  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TFile *finjer = new TFile(Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure_%s.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size), sim_string.c_str()), "r");

  std::cout << "done" << std::endl;

  std::string sys_name = "nominal";
  std::string calib_string = "SMEAR";
  double mix = 0;
  if (full_sys)
    {
      sys_name = "FULL";
    }

  if (prior_sys)
    {
      sys_name = "PRIOR";
    }
  
  std::cout << pileup_sys << std::endl;
  bool mixed = 0;
  if (pileup_sys > 1)
    {
      mix = pileup_sys;
      sys_name = "PILEUP";
      std::cout << "Pileup mix: " << mix << std::endl;
    }
  else if (pileup_sys > 0)
    {
      mixed = true;
      mix = pileup_sys;
      sys_name = "PILEUPMIX";
      std::cout << "Pileup mix: " << mix << std::endl;
    }

  if (emfrac_sys)
    {
      sys_name = "EMFRAC";
    }

  if (trigger_sys)
    {
      //std::cout << "trigger sys" << std::endl;
      sys_name = "TRIGGER";
    }

  if (ca_sys)
    {      
      sys_name = "CA" + std::to_string(ca_sys) ;
    }
  
  if (philoc_sys)
    {
      sys_name = "PHILOC" + std::to_string(philoc_sys);;
      std::cout << sys_name << std::endl;
    }

  if (zyam_sys)
    {
      sys_name = "ZYAM";
    }
  if (inclusive_sys)
    {
      sys_name = "INCLUSIVE";
    }
  if (use_herwig)
    {
      sys_name = "HERWIG";
    }

  gjer = (TF1*) finjer->Get("jernom");
  hjersmear = (TH1D*) finjer->Get("h_nominal_smear");

  float nominal_smear = 0.08;

  if (JER_sys != 0)
    {
      if (JER_sys < 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_down_smear");
	  gjer = (TF1*) finjer->Get("jerneg");
	  calib_string = "SMEAR_DOWN";
	  sys_name = "negJER";
	  nominal_smear = 0.00;
	}
      else if (JER_sys > 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_up_smear");
	  gjer = (TF1*) finjer->Get("jerpos");
	  calib_string = "SMEAR_UP";
	  sys_name = "posJER";
	  nominal_smear = 0.15;
	}
    }

  double smearing = gjer->GetParameter(0);
  fsmear->SetParameter(2, smearing);
  if (JES_sys != 0)
    {
      if (JES_sys < 0)
	sys_name = "negJES";
      else if (JES_sys > 0)
      	sys_name = "posJES";
    }
  if (unfold_generator == 2 && input_generator == 1)
    {
      use_herwig = 0;
      sys_name = "PuH";
    }
  if (unfold_generator == 1 && input_generator == 1)
    {
      use_herwig = 0;
      sys_name = "PuP";
    }
  if (unfold_generator == 2 && input_generator == 2)
    {
      use_herwig = 1;
      sys_name = "HuH";
    }
  if (unfold_generator == 1 && input_generator == 2)
    {
      use_herwig = 1;
      sys_name = "HuP";
    }

  std::string sys_name_orig = sys_name;

  if (full_or_half)
    {
      sys_name = "HALF_" + sys_name;
    }

  const int nbins = read_nbins;

  const int nbins_pt = nbins+1;

  float ipt_bins[nbins_pt+1];
  Double_t dpt_bins[nbins_pt+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  ipt_bins[nbins_pt] = 100;

  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      dpt_bins[i] = (Double_t) ipt_bins[i];
      if (verbosity > 1)
	{	  
	  std::cout << ipt_bins[i] << " -- " << " -- " << dpt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}
    }
  dpt_bins[nbins_pt] = 100;

  Int_t max_reco_bin = rb.get_maximum_reco_bin();

  std::string system_string = "pp";

  int sim_version = (use_herwig? 10 : 11);

  int sam = nsamples;
  int first_sample = 8;
  if (pileup_sys > 1)
    {
      sim_version = 12;
      sam = 2*nsamples;
    }

  const int nsamples_new = 2*nsamples;
  std::string j_file[nsamples_new];
  j_file[0] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_8_ana509_MDC2-00000028.root";
  j_file[1] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_12_ana509_MDC2-00000028.root";
  j_file[2] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
  j_file[3] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
  j_file[4] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_40_ana509_MDC2-00000028.root";
  j_file[5] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_50_ana509_MDC2-00000028.root";
  j_file[6] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_60_ana509_MDC2-00000028.root";

  if (mixed && pileup_sys > 0)
    {
      j_file[7] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_8_ana509_MDC2-00000028.root";
      j_file[8] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_12_ana509_MDC2-00000028.root";
      j_file[9] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_20_ana509_MDC2-00000028.root";
      j_file[10] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_30_ana509_MDC2-00000028.root";
      j_file[11] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_40_ana509_MDC2-00000028.root";
      j_file[12] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v12_50_ana509_MDC2-00000028.root";
    }
  
  bool use_sample[nsamples_new];
  for (int i = 0; i < nsamples_new; i++) use_sample[i] = 1;
  if (mixed)
    {      
      use_sample[nsamples_new - 1] = 0;
    }
  else if (pileup_sys > 1)
    {
      use_sample[nsamples_new - 1] = 0;
      for (int i = 0; i < nsamples;i++)
	{
	  use_sample[i] = 0;
	}
      boundary_r4[5] = 100;
    }
  else
    {
      for (int i = 0; i < nsamples;i++)
	{
	  use_sample[nsamples+i] = 0;
	}
    }
  if (use_herwig)
    {
      
      //use_sample[0] = 0;
      use_sample[6] = 0;
      boundary_r4[6] = 100;
      //      boundary_r4[5] = 100;
    }
  
  float n_events[nsamples_new];
  TFile *finsim[nsamples_new];


  TTree *ttree[nsamples_new];
  ULong64_t gl1_scaled[nsamples_new];

  std::vector<float> *truth_jet_pt_ref[nsamples_new] = {0};
  std::vector<float> *truth_jet_pt[nsamples_new] = {0};
  std::vector<float> *truth_jet_eta[nsamples_new] = {0};
  std::vector<float> *truth_jet_phi[nsamples_new] = {0};
  
  std::vector<float> *reco_jet_pt[nsamples_new] = {0};
  std::vector<float> *reco_jet_pt_uncalib[nsamples_new] = {0};
  std::vector<float> *reco_jet_emcal[nsamples_new] = {0};
  std::vector<float> *reco_jet_e[nsamples_new] = {0};
  std::vector<float> *reco_jet_eta[nsamples_new] = {0};
  std::vector<float> *reco_jet_eta_det[nsamples_new] = {0};
  std::vector<float> *reco_jet_phi[nsamples_new] = {0};

  float truth_vertex_z[nsamples_new];
  float mbd_vertex_z[nsamples_new];
  int mbd_hit[nsamples_new];

  for (int j = 0; j < nsamples_new; j++)
    {

      if (!use_sample[j]) continue;

      finsim[j] = new TFile(j_file[j].c_str(),"r");

      
      ttree[j] = (TTree*) finsim[j]->Get("ttree");
      if (cone_size != 4)
	{
	  ttree[j]->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_ref[j]);
	}
      ttree[j]->SetBranchAddress(Form("truth_jet_pt_%d", cone_size), &truth_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_eta_%d", cone_size), &truth_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_phi_%d", cone_size), &truth_jet_phi[j]);
  
      ttree[j]->SetBranchAddress(Form("jet_pt_calib_%d", cone_size), &reco_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("jet_pt_%d", cone_size), &reco_jet_pt_uncalib[j]);
      ttree[j]->SetBranchAddress(Form("jet_emcal_%d", cone_size), &reco_jet_emcal[j]);
      ttree[j]->SetBranchAddress(Form("jet_e_%d", cone_size), &reco_jet_e[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_%d", cone_size), &reco_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_det_%d", cone_size), &reco_jet_eta_det[j]);
      ttree[j]->SetBranchAddress(Form("jet_phi_%d", cone_size), &reco_jet_phi[j]);
      ttree[j]->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z[j]);
      ttree[j]->SetBranchAddress("mbd_hit", &mbd_hit[j]);
      ttree[j]->SetBranchAddress("truth_vertex_z", &truth_vertex_z[j]);
      ttree[j]->SetBranchAddress("gl1_scaled", &gl1_scaled[j]);
      n_events[j] = ttree[j]->GetEntries();
    }

  
  float cs_8 = (1.3878e7);
  float cs_12 = (1.4903e6);
  float cs_20 = (6.2623e4);
  float cs_30 = (2.5298e3);
  float cs_40 = (1.3553e2);
  float cs_50 = (7.3113);
  float cs_60 = (3.3261e-01);

  float scale_factor[nsamples];
  scale_factor[0] = cs_8/cs_60;
  scale_factor[1] = cs_12/cs_60;
  scale_factor[2] = cs_20/cs_60;
  scale_factor[3] = cs_30/cs_60;
  scale_factor[4] = cs_40/cs_60; 
  scale_factor[5] = cs_50/cs_60;
  scale_factor[6] = 1;
  float cs_h[6] = {1.8437e+08, 6.7108e+05, 5.2613e4, 2.0694e3, 1.0510e2, 5.2089};
  float nevents_h[6] = {10001000, 10913000, 10000000., 10000000., 9671056, 11559000}; 
  float scaless[6] = {1./0.592641,1./0.592641, 1, 1, 1, 1};

  if (use_herwig)
    {
      for (int i = 0 ; i < 6; i++)
	{
	  scale_factor[i] = scaless[i]*cs_h[i]*nevents_h[5]/(cs_h[5]*nevents_h[i]);
	}
      //scale_factor[0]*=2;
    }
  

  TH1D *h_flat_data_over_reco_pt1pt2 = nullptr;

  if (get_mapping == 2)
    {
      // Get the pt1pt2 histograms
      TFile *finunfold = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_PRIMER2_%s.root",  rb.get_code_location().c_str(), cone_size, sys_name_orig.c_str()),"r");
      if (!finunfold)
	{
	  std::cout << " no file " << std::endl;
	  return 0;
	}

      h_flat_data_over_reco_pt1pt2 = (TH1D*) finunfold->Get("h_flat_data_over_reco_pt1pt2");
      h_flat_data_over_reco_pt1pt2->SetDirectory(0);
      finunfold->Close();
    }

  int prior_iteration = 2;

  if (!get_mapping && primer == 0 && !unfold_generator && false)
    {
      std::cout << "opening file for prior" << std::endl;
      TFile *fin_iter = new TFile(Form("unfolding_hists/iteration_tune_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()), "r");      

      TH1D *hunc = (TH1D*) fin_iter->Get("h_total_uncertainties");
      int minimum_iter = 0;
      float mini = hunc->GetBinContent(2);;
      for (int ib = 2; ib < 10; ib++)
	{

	  float temp_min = hunc->GetBinContent(ib);

	  if (mini > temp_min)
	    {
	      minimum_iter = ib - 2;
	      mini = temp_min;
	    }
	  if (temp_min > mini)
	    {
	      break;
	    }
	}

      prior_iteration = minimum_iter;

      std::cout << "Using iteration " << minimum_iter << " for prior reweighting" << std::endl;
    }
  
  // Vertex Reweight

  std::vector<std::pair<float, float>> vertex_scales;
  TH2D *h_eta_reweight = nullptr;
  TH2D *h_lead_pt_v_emfrac_reweight = nullptr;
  if (!(get_mapping==1) && primer != 1)
    {
      TFile *fvtx = new TFile(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()),"r");
      TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");
      for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
	{
	  vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
	}
      h_eta_reweight = (TH2D*) fvtx->Get("h_eta_reweight");
      h_lead_pt_v_emfrac_reweight = (TH2D*) fvtx->Get("h_emfrac_reweight");
    }

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  const int mbins = rb.get_measure_bins();
  int measure_bins[10] = {0};
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }


	
  djf.setTruthCuts(truth_leading_cut, truth_subleading_cut);
  djf.setRecoCuts(reco_leading_cut, reco_subleading_cut);
  djf.setPhiLocation(philoc_sys);
  
  if (verbosity > 1)
    {
      std::cout << "Max reco bin: " << max_reco_bin << std::endl;
      std::cout << "Truth1: " << truth_leading_cut << std::endl;
      std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
      std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
      std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
      std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
      std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
    }

  // for trimming from mapping phase of unfold

  TH1D *h_flat_truth_mapping_primer = nullptr;
  TH1D *h_flat_reco_mapping_primer = nullptr;
  TH2D *h_flat_response_mapping_primer[nsamples] = {0};

  int nbins_pt_truth = nbins_pt*nbins_pt;
  std::map<int, int> mapped_pt_bin_truth;
  int nbins_pt_reco = nbins_pt*nbins_pt;
  std::map<int, int> mapped_pt_bin_reco;

  if (!(get_mapping==1))
    {
      TString mappingppath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
      
      mappingppath += "_MAPPING";
    
      mappingppath += "_" + sys_name_orig;
      
      mappingppath += ".root";

      TFile *fr = new TFile(mappingppath.Data(),"r");

      TNtuple *tn_map = (TNtuple*) fr->Get("tn_mapping");

      float ptbin;
      float mapbin_truth;
      float use_truth;
      float mapbin_reco;
      float use_reco;

      tn_map->SetBranchAddress("ptbin", &ptbin);
      tn_map->SetBranchAddress("mapbin_truth", &mapbin_truth);
      tn_map->SetBranchAddress("use_truth", &use_truth);
      tn_map->SetBranchAddress("mapbin_reco", &mapbin_reco);
      tn_map->SetBranchAddress("use_reco", &use_reco);
      nbins_pt_truth = 0;
      nbins_pt_reco = 0;
      for (int i = 0; i < tn_map->GetEntries(); i++)
	{
	  tn_map->GetEntry(i);
	  if (use_truth == 1)
	    {
	      nbins_pt_truth++;
	      mapped_pt_bin_truth[ptbin] = mapbin_truth;
	    }
	  else
	    {
	      mapped_pt_bin_truth[ptbin] = -1;
	    }
	  if (use_reco == 1)
	    {
	      nbins_pt_reco++;
	      mapped_pt_bin_reco[ptbin] = mapbin_reco;
	    }
	  else
	    {
	      mapped_pt_bin_reco[ptbin] = -1;
	    }
	}
      h_flat_truth_mapping_primer = (TH1D*)fr->Get("h_flat_truth_mapping_all_samples");
      h_flat_reco_mapping_primer = (TH1D*)fr->Get("h_flat_reco_mapping_all_samples");
      
      for (int i = 0; i < nsamples; i++)
	{
	  h_flat_response_mapping_primer[i] = (TH2D*) fr->Get(Form("h_flat_response_mapping_%d", i));
	}      
    }

  int nbins_pt_2 = nbins_pt*nbins_pt;
  int nbins_pt_reco_2 = nbins_pt_reco;
  int nbins_pt_truth_2 = nbins_pt_truth;

  TEfficiency *he_dijet_subleading_binned = new TEfficiency("he_dijet_subleading_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Subleading Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_match_binned = new TEfficiency("he_dijet_match_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Match Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_trigger_binned = new TEfficiency("he_dijet_trigger_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Trigger Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_truth_bad_binned = new TEfficiency("he_dijet_truth_bad_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Truth Bad Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_truth_good_binned = new TEfficiency("he_dijet_truth_good_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Truth Good Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  
  TEfficiency *he_dijet_fake_binned = new TEfficiency("he_dijet_fake_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Fake Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_miss_binned = new TEfficiency("he_dijet_miss_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Miss Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_mbd_hit_binned = new TEfficiency("he_dijet_mbd_hit_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", nbins_pt, 0, nbins_pt);
  TEfficiency *he_dijet = new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_mbd_hit = new TEfficiency("he_dijet_mbd_hit",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_reco = new TEfficiency("he_dijet_reco",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_match = new TEfficiency("he_dijet_match",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_match_reco = new TEfficiency("he_dijet_match_reco",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_fake = new TEfficiency("he_dijet_fake",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_trigger = new TEfficiency("he_dijet_trigger",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_vertex = new TEfficiency("he_dijet_vertex",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_vertex_truth = new TEfficiency("he_dijet_vertex_truth",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);

  TH1D *h_truth_lead_sample[nsamples];
  for (int i = 0; i < nsamples; i++)
    {
      h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }
  
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);


  TH2D *h_truth_lead_pt_eta = new TH2D("h_truth_lead_pt_eta", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100, 110, -1.1, 1.1);
  TH2D *h_truth_sublead_pt_eta = new TH2D("h_truth_sublead_pt_eta", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100, 110, -1.1, 1.1);
  TH2D *h_reco_lead_pt_eta = new TH2D("h_reco_lead_pt_eta", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100, 110, -1.1, 1.1);
  TH2D *h_reco_sublead_pt_eta = new TH2D("h_reco_sublead_pt_eta", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100, 110, -1.1, 1.1);

  TH2D *h_truth_lead_pt_phi = new TH2D("h_truth_lead_pt_phi", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100, 128, -1*TMath::Pi(), TMath::Pi());
  TH2D *h_truth_sublead_pt_phi = new TH2D("h_truth_sublead_pt_phi", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100, 128, -1*TMath::Pi(), TMath::Pi());
  TH2D *h_reco_lead_pt_phi = new TH2D("h_reco_lead_pt_phi", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100, 128, -1*TMath::Pi(), TMath::Pi());
  TH2D *h_reco_sublead_pt_phi = new TH2D("h_reco_sublead_pt_phi", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100, 128, -1*TMath::Pi(), TMath::Pi());

  TH2D *h_eta_lead_sublead = new TH2D("h_eta_lead_sublead","", 220, -1.1, 1.1, 220, -1.1, 1.1);
  TH2D *h_eta_phi_lead = new TH2D("h_eta_phi_lead","", 48, -1.1, 1.1, 128, -TMath::Pi(), TMath::Pi());

  TH2D *h_lead_pt_v_emfrac = new TH2D("h_lead_pt_v_emfrac","",nbins_pt, dpt_bins, 12, -0.1, 1.1);
  
  TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_reco_lead = new TH1D("h_match_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_reco_sublead = new TH1D("h_match_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_njets = new TH1D("h_njets", ";n_{jet}; counts", 21, -0.5, 20.5);
  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  
  
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
  TH1D *h_truth_xj_range[3];

  for (int i = 0; i < 3; i++)
    {
      h_truth_xj_range[i] = new TH1D(Form("h_truth_xj_range_%d", i),";A_{J};1/N", nbins, ixj_bins);
    }
  
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

  TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  TH1D *h_flat_truth_pt1pt2_raw = new TH1D("h_truth_flat_pt1pt2_raw",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
  TH1D *h_flat_truth_pt1pt2_scaled_new = new TH1D("h_truth_flat_pt1pt2_scaled_new",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
  TH1D *h_flat_reco_pt1pt2_raw = new TH1D("h_reco_flat_pt1pt2_raw",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_2, 0, nbins_pt_2);

  TH1D *h_flat_reco_all_pt1pt2 = new TH1D("h_reco_flat_all_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH1D *h_flat_reco_no_fakes_pt1pt2 = new TH1D("h_reco_flat_no_fakes_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);

  TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
  TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
  TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt_reco_2, 0, nbins_pt_reco_2, nbins_pt_truth_2, 0, nbins_pt_truth_2);

  h_flat_truth_to_response_pt1pt2->Sumw2();
  h_flat_reco_to_response_pt1pt2->Sumw2();
  h_flat_truth_pt1pt2_raw->Sumw2();
  h_flat_truth_pt1pt2_scaled_new->Sumw2();
  h_flat_reco_pt1pt2_raw->Sumw2();
  h_flat_reco_pt1pt2->Sumw2();
  h_flat_reco_all_pt1pt2->Sumw2();
  h_flat_response_pt1pt2->Sumw2();
  


  TProfile2D *h_pt1pt2_reco_from_pileup = new TProfile2D("h_pt1pt2_reco_from_pileup",";p_{T,1};p_{T,2}; fraction from pileup", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TProfile2D *h_pt1pt2_truth_from_pileup = new TProfile2D("h_pt1pt2_truth_from_pileup",";p_{T,1};p_{T,2}; fraction from pileup", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
						    
  TH1D *h_count_flat_truth_pt1pt2_sample[nsamples];// new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2_sample[nsamples];// new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH2D *h_count_flat_response_pt1pt2 = new TH2D("h_count_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt_reco_2, 0, nbins_pt_reco_2, nbins_pt_truth_2, 0, nbins_pt_truth_2);

  TH1D *h_flat_reco_all_pt1pt2_sample[nsamples];// = new TH1D("h_reco_flat_all_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2_sample[nsamples];// new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2_sample[nsamples];// new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH2D *h_flat_response_pt1pt2_sample[nsamples];// new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);
  TH2D *h_count_flat_response_pt1pt2_sample[nsamples];// new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  for (int i = 0 ; i < nsamples; i++)
    {
      h_flat_reco_all_pt1pt2_sample[i] = new TH1D(Form("h_reco_flat_all_pt1pt2_%d", i),";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);

      h_count_flat_truth_pt1pt2_sample[i] = new TH1D(Form("h_truth_count_flat_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
      h_count_flat_reco_pt1pt2_sample[i] = new TH1D(Form("h_count_reco_flat_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);

      h_flat_truth_to_response_pt1pt2_sample[i] = new TH1D(Form("h_truth_flat_to_response_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
      h_flat_reco_to_response_pt1pt2_sample[i] = new TH1D(Form("h_reco_flat_to_response_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);

      h_flat_response_pt1pt2_sample[i] = new TH2D(Form("h_flat_response_pt1pt2_%d", i), ";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt_reco_2, 0, nbins_pt_reco_2, nbins_pt_truth_2, 0, nbins_pt_truth_2);
      h_count_flat_response_pt1pt2_sample[i] = new TH2D(Form("h_count_flat_response_pt1pt2_%d", i), ";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt_reco_2, 0, nbins_pt_reco_2, nbins_pt_truth_2, 0, nbins_pt_truth_2);

      h_flat_truth_to_response_pt1pt2_sample[i]->Sumw2();
      h_flat_reco_to_response_pt1pt2_sample[i]->Sumw2();
      h_flat_response_pt1pt2_sample[i]->Sumw2();
  
    }

  
  // For the prior sensitivity
  TH1D *h_flatreweight_pt1pt2 = new TH1D("h_flatreweight_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt_truth_2, 0, nbins_pt_truth_2);
  TH2D *h2_flatreweight_pt1pt2 = new TH2D("h2_flatreweight_pt1pt2",";p_{T,1}^{truth} [GeV]; p_{T,2}^{truth} [GeV] ; Reweight Factor", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  
  if (!prior_sys && !primer && !get_mapping)
    {

      TFile *fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()), "r");
      //TH1D *h_unfold_flat = (TH1D*) fun->Get("h_data_flat_pt1pt2");
      TH1D *h_unfold_flat = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
      std::cout << "Got histograms" << std::endl;
      TFile *ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()), "r");

      TH1D *h_truth_flat = (TH1D*) ftr->Get("h_truth_flat_pt1pt2");
      h_truth_flat->SetName("h_truth_flat");
      std::cout << "Got histograms" << std::endl;
      

      // of signal region
      
      double integral_of_signal = 0;
      double integral_of_truth = 0;
      double pt1_unfold_bins[nbins_pt];
      double pt1_truth_bins[nbins_pt];
      for (int i = 0; i < nbins_pt; i++)
	{
	  pt1_unfold_bins[i] = 0;
	  pt1_truth_bins[i] = 0;
	}

      //h_unfold_flat->Scale(1./h_unfold_flat->Integral(),"width");
      //h_truth_flat->Scale(1./h_truth_flat->Integral(),"width");

      for (int ibin = 0; ibin < nbins_pt_2; ibin++)
	{

	  int pt1_bin = ibin/nbins_pt;
	  int pt2_bin = ibin%nbins_pt;

	  int max_bin = std::max(pt1_bin, pt2_bin);
	  int min_bin = std::min(pt1_bin, pt2_bin);
	  //if (max_bin < measure_leading_bin) continue;
	  //if (min_bin < measure_subleading_bin) continue;

	  pt1_unfold_bins[max_bin] += h_unfold_flat->GetBinContent(ibin+1);
	  pt1_truth_bins[max_bin] += h_truth_flat->GetBinContent(ibin+1);

	  integral_of_signal += h_unfold_flat->GetBinContent(ibin+1);
	  integral_of_truth += h_truth_flat->GetBinContent(ibin+1);
	}

      
      for (int ibin = 0; ibin < nbins_pt_2; ibin++)
	{

	  int pt1_bin = ibin/nbins_pt;
	  int pt2_bin = ibin%nbins_pt;
	  int max_bin = std::max(pt1_bin, pt2_bin);
	  int min_bin = std::min(pt1_bin, pt2_bin);
	  if (pt1_unfold_bins[max_bin] == 0 || pt1_truth_bins[max_bin] == 0)
	    {
	      continue;
	    }
	  else
	    {
	      h_unfold_flat->SetBinContent(ibin+1, h_unfold_flat->GetBinContent(ibin+1)/pt1_unfold_bins[max_bin]);
	      h_truth_flat->SetBinContent(ibin+1, h_truth_flat->GetBinContent(ibin+1)/pt1_truth_bins[max_bin]);
	    }
	}

      for (int ibin = 0; ibin < nbins_pt_2; ibin++)
	{
	  float pt1_bin = ibin/nbins_pt;
	  float pt2_bin = ibin%nbins_pt;
	  int max_bin = std::max(pt1_bin, pt2_bin);
	  int gbin = h2_flatreweight_pt1pt2->GetBin(pt1_bin+1, pt2_bin+1);
	  
	  float v = h_unfold_flat->GetBinContent(ibin+1);
	  float b = h_truth_flat->GetBinContent(ibin+1);
	  
	  if (b > 0 && v > 0)
	    {
	      float vb = (v/b);// * lead_factor;
	      if (vb > 5) vb = 5;
	      if (vb < 0.1) vb = 0.1;
	      h_flatreweight_pt1pt2->SetBinContent(ibin+1, vb);
	      h2_flatreweight_pt1pt2->SetBinContent(gbin, vb);
	    }
	  else
	    {
	      continue;//      h_flatreweight_pt1pt2->SetBinContent(ibin+1, 1);
	      continue;//h2_flatreweight_pt1pt2->SetBinContent(gbin, 1);
	    }

	}
      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  double last_weight = 0;
	  for (int id = 0; id <= ib; id++)
	    {
	      int idd = ib - id;
	      int ibb = 1 + ib*nbins_pt + idd;
	      int ibb1 = 1 + idd*nbins_pt + ib;
	      double weight = h_flatreweight_pt1pt2->GetBinContent(ibb);
	      if (weight == 0)
		{
		  h_flatreweight_pt1pt2->SetBinContent(ibb, last_weight);
		  h_flatreweight_pt1pt2->SetBinContent(ibb1, last_weight);
		}
	      else
		{
		  last_weight = weight;
		}
	    }
	}

      TCanvas *cre = new TCanvas("cre","cre", 500, 500);
      cre->SetLeftMargin(0.1);
      cre->SetRightMargin(0.19);
      h2_flatreweight_pt1pt2->Draw("colz");
      dlutility::DrawSPHENIX(0.2, 0.87);
      dlutility::drawText("Prior Reweighting Matrix", 0.2, 0.77);
      cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name_orig.c_str()));
    }

  int nbin_response = nbins_pt_2;
  
  RooUnfoldResponse rooResponse(h_flat_reco_to_response_pt1pt2, h_flat_truth_to_response_pt1pt2);

  TRandom *rng = new TRandom(0);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> mysmearjets;
  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets2;
  std::vector<struct jet> mysmearjets2;
  std::vector<struct jet> myrecojets2;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  std::vector<std::pair<struct jet, struct jet>> matched_dijets2;

  for (int isample = 0; isample < nsamples_new; isample++)
    {
      if (!use_sample[isample]) continue;

      int iisample = isample%nsamples;
      std::cout << "Sample " << isample << std::endl;
      int entries0 = 0;
      int entries2 = ttree[isample]->GetEntries();

      int starting = 1;
      if (verbosity > 5) entries2 = 100000;

      if (iisample != isample)
	{
	  entries2 = n_events[iisample] * mix;
	}
      
      for (int i = entries0; i < entries2; i++)
	{
	  ttree[isample]->GetEntry(i);
	  double event_scale = scale_factor[iisample];
	  if (verbosity > 5)
	    {
	      std::cout << "------------- Event " << i << " --------------" << std::endl;
	    }
	  bool has_vertex = (fabs(mbd_vertex_z[isample]) < vtx_cut);
	  bool has_truth_vertex = (fabs(truth_vertex_z[isample]) < vtx_cut);
	  bool has_mbd_hit = mbd_hit[isample];

	  if ((!has_vertex || !has_mbd_hit) && !full_sys) continue;

	  if (full_sys)
	    {
	      if ((!has_vertex || !has_mbd_hit) && (!has_truth_vertex)) continue;
	    }

	  bool triggered = (trigger_sys ? true : ( ( ( gl1_scaled[isample] >> 22 ) & 0x1) == 0x1) );

	  // Full or half closure 
	  bool fill_response = 1;
	  if (full_or_half && i < (entries2/2))
	    {
	      if (starting)
		{
		  std::cout << "Starting Half" << std::endl;
		  starting = 0;		  
		}
	      fill_response = 0;
	    }

	  // Fill the unfolded matrices?
	  bool fill_unfold = (!full_or_half) || (!fill_response && full_or_half);

	  // make sure sample is good 
	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (cone_size != 4)
	    {
	      maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
	    }

	  if ((isample < nsamples) && (maxpttruth < boundary_r4[iisample] || maxpttruth >= boundary_r4[iisample+1])) continue;
	  if ((isample >= nsamples) && (maxpttruth < boundary_pileup_r4[iisample] || maxpttruth >= boundary_pileup_r4[iisample+1])) continue;

	  // vertex reweighting
	  if (!(get_mapping==1) && primer != 1 && has_vertex)
	    {
	      for (int ib = 0; ib < vertex_scales.size(); ib++)
		{
		  if (mbd_vertex_z[isample] < vertex_scales.at(ib).first)
		    {
		      event_scale *= vertex_scales.at(ib).second;
		      break;
		    }
		}

	    }

	  // collect the truth and reco jets
	  // fully matched dijets
	  matched_dijets.clear();
	  matched_dijets2.clear();
	  // all truth jets
	  mytruthjets.clear();
	  // all reco jets
	  myrecojets.clear();
	  mysmearjets.clear();
	  // top two of each
	  mytruthjets2.clear();
	  myrecojets2.clear();
	  mysmearjets2.clear();

	  int ntruthjets = truth_jet_pt[isample]->size();
	  for (int j = 0; j < ntruthjets;j++)
	    {

	      //if (truth_jet_pt[isample]->at(j) < truth_subleading_cut) continue;
	      //if (fabs(truth_jet_eta[isample]->at(j)) > 1.1) continue;

	      struct jet tempjet;
	      tempjet.istruth = 1;
	      tempjet.pt = truth_jet_pt[isample]->at(j);
	      tempjet.eta = truth_jet_eta[isample]->at(j);
	      tempjet.phi = truth_jet_phi[isample]->at(j);
	      tempjet.id = j;

	      mytruthjets.push_back(tempjet);	  	  

	    }

	  std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  int nrecojets = reco_jet_pt[isample]->size();
	  int nnrecojets = 0;
	  for (int j = 0; j < nrecojets;j++)
	    {
	      // no negative energy jets

	      if (reco_jet_e[isample]->at(j) < 0) continue;

	      if (reco_jet_pt[isample]->at(j) > 5) nnrecojets++;

	      struct jet tempjet;

	      tempjet.istruth = 0;

	      float temppt = reco_jet_pt[isample]->at(j);
	      //float tempptreco = temppt;
	      //float smear1 = gjer->Eval(temppt);//nominal_smear;//hjersmear->GetBinContent(ib );
	      //float jersmear = 0;

	      // if (smear1 > 0)
	      // 	{
	      // 	  fsmear->SetParameter(2, smear1);
	      // 	  jersmear = fsmear->GetRandom();
	      // 	}
	      
	      //tempptreco = temppt + jersmear*temppt;// + JES_sys * temppt;
	
	      tempjet.pt = temppt;
	      tempjet.pt_uncalib = reco_jet_pt[isample]->at(j);
	      tempjet.e = reco_jet_e[isample]->at(j);

	      //if (tempjet.pt < reco_subleading_cut) continue;
	      
	      tempjet.emcal = reco_jet_emcal[isample]->at(j);
	      tempjet.eta = reco_jet_eta[isample]->at(j);
	      tempjet.eta_det = reco_jet_eta_det[isample]->at(j);
	      tempjet.phi = reco_jet_phi[isample]->at(j);
	      tempjet.t = 0;
	      tempjet.id = j;

	      myrecojets.push_back(tempjet);
	    }

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
	  
	  // new smearing technique
	  
	  matched_dijets = djf.match_dijets_smear(myrecojets, mytruthjets);
	  for (auto recojet : myrecojets)
	    {
	      bool mm = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == recojet.id;}) != matched_dijets.end();
	      he_dijet_match->Fill(mm, recojet.pt, event_scale);
	      //	      he_dijet_match_reco->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	    }

	  for (auto match :  matched_dijets)
	    {
	      auto recojet = match.second;
	      auto truthjet = match.first;
	      double smear1 = fsmear->GetRandom();
	      double temppt = recojet.pt;
	      recojet.pt = temppt + truthjet.pt*smear1;  
	      mysmearjets.push_back(recojet);		  
	      //matched_dijets2.push_back(std::make_pair(truthjet, recojet));
	    }
	  //he_dijet_match_reco->Fill(base_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	  //he_dijet_match_reco->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	  
 	  // std::sort(matched_dijets2.begin(), matched_dijets2.end(), [] (auto a, auto b) { return a.first.pt > b.first.pt; });
	  // if (matched_dijets2.size() > 1)
	  //   {
	  //     matched_dijets2 = {matched_dijets2.begin(), matched_dijets2.begin()+2};
	  //   }
	  
	  //matched_dijets = match_dijets_ppg08(myrecojets, mytruthjets);
	  for (auto recojet : myrecojets)
	    {
	      bool mm = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == recojet.id;}) != matched_dijets.end();
	      if (mm) continue;
	      //double smear1 = fsmear->GetRandom();
	      //double temppt = recojet.pt;
	      //recojet.pt = temppt + recojet.pt*smear1;  
	      mysmearjets.push_back(recojet);		  	      
	    }

	  std::sort(mysmearjets.begin(), mysmearjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  matched_dijets2 = djf.match_dijets_response(mysmearjets, mytruthjets);
	  for (auto recojet : mysmearjets)
	    {
	      bool mm = std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.second.id == recojet.id;}) != matched_dijets.end();
	      he_dijet_match_reco->Fill(mm, recojet.pt, event_scale);
	      //	      he_dijet_match_reco->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	    }

	  if (isample == 2 && verbosity >= 5)
	    {
	      int imm = 0;
	      std::cout << "matched: " << std::endl;
	      for (auto mjet : matched_dijets2)
		{
		  std::cout << "Match : " << imm++ << " : " << mjet.first.pt << " / " << mjet.first.eta << " / " << mjet.first.phi << "     " << mjet.second.pt << " / " << mjet.second.eta << " / " << mjet.second.phi << std::endl;
		}
	    }

	  std::sort(matched_dijets2.begin(), matched_dijets2.end(), [] (auto a, auto b) { return a.first.pt > b.first.pt; });

	  if (matched_dijets2.size() > 1)
	    {
	      matched_dijets2 = {matched_dijets2.begin(), matched_dijets2.begin()+2};
	    }

	  
	  bool truth_good = djf.check_dijet_truth(mytruthjets);
	  bool truth_good_base = mytruthjets.size() > 1;


	  if (full_sys)
	    {
	      truth_good &= has_truth_vertex;
	    }

	  if (myrecojets.size() > 1)
	    {
	      int sizet = myrecojets.size();
	      myrecojets2 = {myrecojets.begin(), myrecojets.begin()+sizet};
	      if (verbosity > 5)
		{
		  std::cout << "RECO GOOD" << std::endl;
		}		  
	    }
	  if (mysmearjets.size() > 1)
	    {
	      int sizet = mysmearjets.size();
	      mysmearjets2 = {mysmearjets.begin(), mysmearjets.begin()+2};
	      if (verbosity > 5)
		{
		  std::cout << "RECO GOOD" << std::endl;
		}		  
	    }

	  if (truth_good)
	    {
	      mytruthjets2 = {mytruthjets.begin(), mytruthjets.begin()+2};
	      if (verbosity > 5)
		{
		  std::cout << "TRUTH GOOD" << std::endl;
		}		  
	    }
	  // matched_dijets = djf.match_dijets(myrecojets, mytruthjets);
	  // match the truth dijet with all reco jets for dphi then dR matching
	  	  


	  bool matched = (matched_dijets2.size() > 1);
	  if (matched)
	    {
	      bool found1 = (std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.second.id == mysmearjets.at(0).id ; }) != matched_dijets2.end());
	      bool found2 = (std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.second.id == mysmearjets.at(1).id ; }) != matched_dijets2.end());
	      bool foundt1 = (std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.first.id == mytruthjets.at(0).id ; }) != matched_dijets2.end());
	      bool foundt2 = (std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.first.id == mytruthjets.at(1).id ; }) != matched_dijets2.end());
	      matched &= (foundt1 && foundt2);
	      matched &= (found1 && found2);

	    }

	  float max_truth = 0;
	  float max_reco = 0;
	  float min_truth = 0;
	  float min_reco = 0;

	  float max_phi_truth = 0;
	  float max_phi_reco = 0;
	  float min_phi_truth = 0;
	  float min_phi_reco = 0;
	  float max_eta_truth = 0;
	  float max_eta_reco = 0;
	  float min_eta_truth = 0;
	  float min_eta_reco = 0;

	  /* this variable determines how things are filled */

	  // 0 - truth and reco good and matched
	  // 1 - truth good and no reco
	  // 2 - truth and reco good and no match (fake and miss)
	  // 3 - truth not good and reco good and no match (fake)

	  int fill_fake_miss = -1;

	  bool reco_good = djf.check_dijet_reco(mysmearjets);

	  bool base_good = reco_good;	  
	  // passes njet cut
	  bool njet_pass = (mysmearjets.size() < njet_cut);
	  reco_good &= njet_pass;
	  // passes trigger cut
	  reco_good &= triggered;
	  
	  if (full_sys)
	    {
	      // passes full cut
	      reco_good &= has_vertex;
	      reco_good &= has_mbd_hit;
	    }

	  // passes in range
	  bool in_range = false;
	  if (mysmearjets2.size() > 0)
	    {
	      in_range = (mysmearjets2.at(0).pt < ipt_bins[max_reco_bin]);
	      reco_good &= in_range;
	    }

	  // if truth good in its base, just 2 jets to get reweighting
	  if (truth_good_base)
	    {
	      max_truth = mytruthjets.at(0).pt;
	      min_truth = mytruthjets.at(1).pt;
	      max_phi_truth = mytruthjets.at(0).phi;
	      min_phi_truth = mytruthjets.at(1).phi;
	      max_eta_truth = mytruthjets.at(0).eta;
	      min_eta_truth = mytruthjets.at(1).eta;
	    }
	  
	  if (reco_good)
	    {
	      max_reco = mysmearjets.at(0).pt;
	      min_reco = mysmearjets.at(1).pt;
	      max_phi_reco = mysmearjets.at(0).phi;
	      min_phi_reco = mysmearjets.at(1).phi;
	      max_eta_reco = mysmearjets.at(0).eta;
	      min_eta_reco = mysmearjets.at(1).eta;
	    }

	  if (!(get_mapping==1) && primer != 1 && reco_good && false)
	    {
	      int eta1 = floor((mysmearjets2.at(0).eta + 1.1)/0.1) + 1;
	      int eta2 = floor((mysmearjets2.at(1).eta + 1.1)/0.1) + 1;
	      int gbinn = h_eta_reweight->GetBin(eta1, eta2);
	      event_scale *= h_eta_reweight->GetBinContent(gbinn);
	    }

	  
	  if (!(get_mapping==1) && primer != 1 && reco_good && emfrac_sys)
	    {
	      int ptembin = h_lead_pt_v_emfrac_reweight->FindBin(mysmearjets.at(0).pt, mysmearjets.at(0).emcal);
	      event_scale *= h_lead_pt_v_emfrac_reweight->GetBinContent(ptembin);
	    }

	  
	  if (matched && truth_good && reco_good)
	    {
	      max_truth = matched_dijets2.at(0).first.pt;
	      min_truth = matched_dijets2.at(1).first.pt;

	      max_reco = matched_dijets2.at(0).second.pt;
	      min_reco = matched_dijets2.at(1).second.pt;

	      max_phi_truth = matched_dijets2.at(0).first.phi;
	      min_phi_truth = matched_dijets2.at(1).first.phi;
	      max_eta_truth = matched_dijets2.at(0).first.eta;
	      min_eta_truth = matched_dijets2.at(1).first.eta;

	      max_phi_reco = matched_dijets2.at(0).second.phi;
	      min_phi_reco = matched_dijets2.at(1).second.phi;
	      max_eta_reco = matched_dijets2.at(0).second.eta;
	      min_eta_reco = matched_dijets2.at(1).second.eta;
	    }

	  float pt1_truth_bin = nbins_pt;
	  float pt2_truth_bin = nbins_pt;
	  float pt1_reco_bin = nbins_pt;
	  float pt2_reco_bin = nbins_pt;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;

	  float p1 = max_phi_truth;
	  float p2 = min_phi_truth;
	  float ps1 = max_phi_reco;
	  float ps2 = min_phi_reco;

	  float eta1 = max_eta_truth;
	  float eta2 = min_eta_truth;
	  float etas1 = max_eta_reco;
	  float etas2 = min_eta_reco;

	  float maxi = es1;
	  float mini = es2;
	  if (mini > maxi);
	  {
	    mini = es1;
	    maxi = es2;
	    float tempp = ps1;
	    ps1 = ps2;
	    ps2 = tempp;
	    tempp = etas1;
	    etas1 = etas2;
	    etas2 = tempp;

	  }
	  float maxit = std::max(e1, e2);
	  float minit = std::min(e1, e2);


	  // get bins
	  for (int ib = 0; ib < nbins_pt; ib++)
	    {
	      if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
		{
		  pt1_truth_bin = ib;
		}
	      if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
		{
		  pt2_truth_bin = ib;
		}
	      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
		{
		  pt1_reco_bin = ib;
		}
	      if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
		{
		  pt2_reco_bin = ib;
		}
	    }

	  if (pt1_truth_bin == nbins_pt)
	    {
	      truth_good = false;
	      truth_good_base = false;
	    }

	  if (pt2_truth_bin == nbins_pt)
	    {
	      truth_good_base = false;
	      truth_good = false;
	    }

	  int total_bin_truth_1 = pt1_truth_bin + nbins_pt*pt2_truth_bin;
	  int total_bin_truth_2 = pt2_truth_bin + nbins_pt*pt1_truth_bin;
	  int total_bin_reco_1 = pt1_reco_bin + nbins_pt*pt2_reco_bin;
	  int total_bin_reco_2 = pt2_reco_bin + nbins_pt*pt1_reco_bin;

	  // reweight all truth jet pairs the same since there can't be any monojets

	  double event_scale_fake = event_scale;
	  if (!prior_sys && !get_mapping && !primer)// && reco_good && matched)
	    {
	      int recorrectbin = pt1_truth_bin + nbins_pt*pt2_truth_bin;	      
	      event_scale *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
	      //event_scale_fake *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
	    }
	  double event_scale_new = event_scale;
	  if (!prior_sys && get_mapping==2)// && !primer)// && reco_good && matched)
	    {
	      event_scale_new *= h_flat_data_over_reco_pt1pt2->GetBinContent(total_bin_reco_1);
	    }

	  // trigger efficiency
	  if (base_good && njet_pass)
	    {
	      he_dijet_trigger_binned->Fill(triggered, mysmearjets2.at(0).pt, mysmearjets2.at(1).pt, event_scale);
	      he_dijet_trigger_binned->Fill(triggered, mysmearjets2.at(1).pt, mysmearjets2.at(0).pt, event_scale);
	    }

	  // Fill efficiencies
	  
	  // check if leading truth jet is matched
	  if (mytruthjets2.size() > 1)
	    {
	      if (std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.first.id == mytruthjets2.at(0).id;}) != matched_dijets2.end())
		{
		  
		  bool subleading_matched = std::find_if(matched_dijets2.begin(), matched_dijets2.end(), [=] (auto a) { return a.first.id == mytruthjets2.at(1).id;}) != matched_dijets2.end();
		  
		  he_dijet_match_binned->Fill(subleading_matched, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
		  he_dijet_match_binned->Fill(subleading_matched, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);

		  if (fabs(mytruthjets2.at(0).eta) < 0.7)
		    {
		      bool sub_in = fabs(mytruthjets2.at(1).eta) < 0.7;
		      he_dijet_subleading_binned->Fill(sub_in, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
		      he_dijet_subleading_binned->Fill(sub_in, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);		      
		      if (sub_in)
			{

			  he_dijet_truth_good_binned->Fill(base_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
			  he_dijet_truth_good_binned->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);		      

			}
		      else
			{

			  he_dijet_truth_bad_binned->Fill(base_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
			  he_dijet_truth_bad_binned->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);		      
			  
			}
		    }
		}
	    }
	      
	  if (base_good)
	    {
	      he_dijet_reco->Fill(truth_good, mysmearjets2.at(0).pt, mysmearjets2.at(1).pt, event_scale);
	      he_dijet_reco->Fill(truth_good, mysmearjets2.at(1).pt, mysmearjets2.at(0).pt, event_scale);
	    }
	  if (truth_good)
	    {

	      he_dijet_vertex->Fill(has_vertex, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_vertex->Fill(has_vertex, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	  
	      he_dijet_mbd_hit->Fill(has_mbd_hit, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_mbd_hit->Fill(has_mbd_hit, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);

	      he_dijet_trigger->Fill(triggered, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_trigger->Fill(triggered, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);	  
	  
	      he_dijet->Fill(matched && reco_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	      he_dijet->Fill(matched && reco_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      
	    }

	  //if (!has_vertex || !has_mbd_hit) continue;	  

	  // Get mapping	  
	  if (!(get_mapping==1))
	    {
	      
	      if (mapped_pt_bin_truth[pt1_truth_bin*nbins_pt + pt2_truth_bin] == -1)
		{
		  truth_good = false;
		}
	      else
		{
		  total_bin_truth_1 = mapped_pt_bin_truth[total_bin_truth_1];
		  total_bin_truth_2 = mapped_pt_bin_truth[total_bin_truth_2];
		}
	      
	      if (mapped_pt_bin_reco[pt1_reco_bin*nbins_pt + pt2_reco_bin] == -1)
		{
		  reco_good = false;
		}
	      else
		{
		  total_bin_reco_1 = mapped_pt_bin_reco[total_bin_reco_1];
		  total_bin_reco_2 = mapped_pt_bin_reco[total_bin_reco_2];
		}
	  
	      if (matched && reco_good && truth_good)
		{
		  int binn = h_flat_response_mapping_primer[iisample]->GetBin(1 + pt1_reco_bin*nbins_pt + pt2_reco_bin, 1 + pt1_truth_bin*nbins_pt + pt2_truth_bin);
		  if (h_flat_response_mapping_primer[iisample]->GetBinContent(binn) == 0) continue;
		}	      


	    }

	  if (total_bin_reco_1 <= 0 && reco_good)
	    {
	      std::cout << "BAD" << std::endl;
	      std::cout << "TOTAL: " << total_bin_reco_1 << std::endl;
	      std::cout << "pt1_reco_bin: " << pt1_reco_bin << std::endl;
	      std::cout << "pt2_reco_bin: " << pt2_reco_bin << std::endl;
	      std::cout << "RECO" << std::endl;
	      for (auto recojet : myrecojets)
		{
		  recojet.print();
		}
	      std::cout << "TRUTH" << std::endl;
	      for (auto truthjet : mytruthjets)
		{
		  truthjet.print();
		}
	      std::cout << "SMEAR" << std::endl;
	      for (auto smearjet : mysmearjets)
		{
		  smearjet.print();
		}

	      
	    }

	  if (!reco_good && !truth_good) continue;

	  if (reco_good && truth_good && matched)
	    {
	      fill_fake_miss = 0;

	      he_dijet_miss_binned->Fill(1, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(1, e2, e1, event_scale);

	      he_dijet_fake_binned->Fill(1, es1, es2, event_scale);
	      he_dijet_fake_binned->Fill(1, es2, es1, event_scale);

	    }

	  if (reco_good && truth_good && !matched)
	    {

	      fill_fake_miss = 2;	      
	      
	      he_dijet_fake_binned->Fill(0, es1, es2, event_scale_fake);
	      he_dijet_fake_binned->Fill(0, es2, es1, event_scale_fake);

	      he_dijet_miss_binned->Fill(0, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(0, e2, e1, event_scale);
	    }
	  
	  if (reco_good && !truth_good)
	    {
	      fill_fake_miss = 3;

	      he_dijet_fake_binned->Fill(0, es1, es2, event_scale_fake);
	      he_dijet_fake_binned->Fill(0, es2, es1, event_scale_fake);
	      //if (!useFakes) continue;
	    }
	  if (truth_good && !reco_good)
	    {
	      fill_fake_miss = 1;
	      he_dijet_miss_binned->Fill(0, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(0, e2, e1, event_scale);
	    }
	  
	  if (truth_good)
	    {
	      h_truth_lead->Fill(e1, event_scale);
	      h_truth_sublead->Fill(e2, event_scale);

	      h_truth_lead_pt_eta->Fill(e1,eta1, event_scale);
	      h_truth_sublead_pt_eta->Fill(e2,eta2, event_scale);
	      h_truth_lead_pt_phi->Fill(e1,p1, event_scale);
	      h_truth_sublead_pt_phi->Fill(e2,p2, event_scale);

	      h_truth_lead_sample[iisample]->Fill(e1, event_scale);
	    }

	  if (reco_good)
	    {
	      h_reco_lead->Fill(maxi, event_scale);
	      h_reco_sublead->Fill(mini, event_scale);
	      h_reco_lead_pt_eta->Fill(maxi,etas1, event_scale);
	      h_reco_sublead_pt_eta->Fill(mini,etas2, event_scale);
	      h_reco_lead_pt_phi->Fill(maxi,ps1, event_scale);
	      h_reco_sublead_pt_phi->Fill(mini,ps2, event_scale);
	      
	      
	      h_eta_phi_lead->Fill(mysmearjets.at(0).eta, mysmearjets.at(0).phi, event_scale);
	    }
	  if (verbosity > 5)
	    {
	      std::cout << "FILL FAKE MISS = " << fill_fake_miss << std::endl;
	    }
	  
	  // miss
	  if (fill_fake_miss == 1)
	    {
	      
	      if (fill_response)
		{
		  
		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_2, event_scale);

		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_2);

		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_2, event_scale);

		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_2);
		  h_pt1pt2_truth_from_pileup->Fill(e1, e2, (isample != iisample ? 1 : 0), event_scale);
		  h_pt1pt2_truth_from_pileup->Fill(e2, e1, (isample != iisample ? 1 : 0), event_scale);
		  rooResponse.Miss(total_bin_truth_1, event_scale);
		  rooResponse.Miss(total_bin_truth_2, event_scale);

		}

	      if (fill_unfold)
		{
		  if (minit >= ipt_bins[measure_subleading_bin] && maxit >= measure_leading_cut && maxit < ipt_bins[measure_bins[3]])
		    {
		      h_truth_xj->Fill(minit/maxit, event_scale);
		      for (int irange = 0; irange < 3; irange++)
			{
			  if (minit >= ipt_bins[measure_subleading_bin] &&  maxit >= ipt_bins[measure_bins[irange]] && maxit < ipt_bins[measure_bins[irange+1]])
			    {
			      h_truth_xj_range[irange]->Fill(minit/maxit, event_scale);
			      break;
			    }
			}
		      h_linear_truth_xj->Fill(minit/maxit, event_scale);
		    }

		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_1, event_scale_new);
		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_2, event_scale_new);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_2, event_scale);

		  h_e1e2->Fill(e1, e2, event_scale);
		  h_e1e2->Fill(e2, e1, event_scale);

		}
	      continue;
	    }

	  // Fill
	  if (fill_fake_miss == 0)
	    {
	      
	      if (maxit > truth_leading_cut) h_match_truth_lead->Fill(maxit, event_scale);
	      if (minit >  truth_subleading_cut) h_match_truth_sublead->Fill(minit, event_scale);
	      if (maxi >  reco_leading_cut) h_match_reco_lead->Fill(maxi, event_scale);
	      if (mini >  reco_subleading_cut) h_match_reco_sublead->Fill(mini, event_scale);
	      
	      	      
	      if (fill_response)
		{

		  h_flat_reco_to_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_to_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, event_scale);

		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, event_scale);

		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_2, event_scale);
		  h_flat_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, total_bin_truth_1, event_scale);
		  h_flat_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, total_bin_truth_2, event_scale);

		  h_count_flat_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, total_bin_truth_1);
		  h_count_flat_response_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, total_bin_truth_2);
		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_2);
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_2);	      

		  h_pt1pt2_reco_from_pileup->Fill(maxi, mini, (isample != iisample ? 1 : 0), event_scale);
		  h_pt1pt2_reco_from_pileup->Fill(mini, maxi, (isample != iisample ? 1 : 0), event_scale);

		  h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_2, event_scale);

		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_2, event_scale);

		  h_flat_reco_no_fakes_pt1pt2->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_no_fakes_pt1pt2->Fill(total_bin_reco_2, event_scale);
		  
		  h_pt1pt2_truth_from_pileup->Fill(e1, e2, (isample != iisample ? 1 : 0), event_scale);
		  h_pt1pt2_truth_from_pileup->Fill(e2, e1, (isample != iisample ? 1 : 0), event_scale);

		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_2, event_scale);

		  h_flat_response_pt1pt2->Fill(total_bin_reco_1, total_bin_truth_1, event_scale);
		  h_flat_response_pt1pt2->Fill(total_bin_reco_2, total_bin_truth_2, event_scale);

		  h_count_flat_response_pt1pt2->Fill(total_bin_reco_1, total_bin_truth_1);
		  h_count_flat_response_pt1pt2->Fill(total_bin_reco_2, total_bin_truth_2);
		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_2);
		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_2);	      

		  
		  rooResponse.Fill(total_bin_reco_1,total_bin_truth_1, event_scale);
		  rooResponse.Fill(total_bin_reco_2,total_bin_truth_2, event_scale);

		  
		}
	      if (fill_unfold)
		{
		  if (minit >= ipt_bins[measure_subleading_bin] && maxit >= measure_leading_cut && maxit < ipt_bins[measure_bins[3]])
		    {
		      h_truth_xj->Fill(minit/maxit, event_scale);
		      for (int irange = 0; irange < 3; irange++)
			{
			  if (minit >= ipt_bins[measure_subleading_bin] &&  maxit >= ipt_bins[measure_bins[irange]] && maxit < ipt_bins[measure_bins[irange+1]])
			    {
			      h_truth_xj_range[irange]->Fill(minit/maxit, event_scale);
			      break;
			    }
			}
		      h_linear_truth_xj->Fill(minit/maxit, event_scale);
		    }
		  // if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
		  //   {
		      
		  //     h_truth_xj->Fill(minit/maxit, event_scale);
		  //     for (int irange = 0; irange < 3; irange++)
		  // 	{
		  // 	  if (minit >= measure_subleading_cut &&  maxit >= ipt_bins[measure_bins[irange]] && maxit < ipt_bins[measure_bins[irange+1]])
		  // 	    {
		  // 	      h_truth_xj_range[irange]->Fill(minit/maxit, event_scale);
		  // 	    }
		  // 	}

		  //     h_linear_truth_xj->Fill(minit/maxit, event_scale);
		  //   }

		  h_pt1pt2->Fill(es1, es2, event_scale);
		  h_pt1pt2->Fill(es2, es1, event_scale);
		  h_e1e2->Fill(e1, e2, event_scale);
		  h_e1e2->Fill(e2, e1, event_scale);

		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_1, event_scale_new);
		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_2, event_scale_new);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_2, event_scale);
		  
		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_1, event_scale);
		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_2, event_scale);
		}		
		
		
	      if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
		{
		  h_reco_xj->Fill(mini/maxi, event_scale);
		  h_linear_reco_xj->Fill(mini/maxi, event_scale);
		}

	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(mysmearjets2.at(0).eta, mysmearjets2.at(1).eta, event_scale);
		  h_lead_pt_v_emfrac->Fill(mysmearjets.at(0).pt, mysmearjets.at(0).emcal, event_scale); 	  
		  h_mbd_vertex->Fill(mbd_vertex_z[isample], event_scale);
		  h_njets->Fill(nnrecojets, event_scale);	      
		}
	      continue;
	    } 

	  // both fake and miss
	  if (fill_fake_miss == 2)
	    {
	      
	      if (fill_response)
		{
		  if(useFakes)
		    {
		      rooResponse.Fake(total_bin_reco_1, event_scale_fake);
		      rooResponse.Fake(total_bin_reco_2, event_scale_fake);
		      h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_1, event_scale_fake);
		      h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_2, event_scale_fake);
		    }

		  h_pt1pt2_reco_from_pileup->Fill(maxi, mini, (isample != iisample ? 1 : 0), event_scale_fake);
		  h_pt1pt2_reco_from_pileup->Fill(mini, maxi, (isample != iisample ? 1 : 0), event_scale_fake);
		  
		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_2);
		  
		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_2, event_scale_fake);
		  
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_2);
		  
		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, event_scale_fake);		 		

		  h_pt1pt2_truth_from_pileup->Fill(e1, e2, (isample != iisample ? 1 : 0), event_scale);
		  h_pt1pt2_truth_from_pileup->Fill(e2, e1, (isample != iisample ? 1 : 0), event_scale);


		  rooResponse.Miss(total_bin_truth_1, event_scale);
		  rooResponse.Miss(total_bin_truth_2, event_scale);
		  
		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[iisample]->Fill(total_bin_truth_2, event_scale);

		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2_sample[iisample]->Fill(total_bin_truth_2);

		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(total_bin_truth_2, event_scale);

		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_1);
		  h_count_flat_truth_pt1pt2->Fill(total_bin_truth_2);
		  
		}
	      if (fill_unfold)
		{
		  if (minit >= ipt_bins[measure_subleading_bin] && maxit >= measure_leading_cut && maxit < ipt_bins[measure_bins[3]])
		    {
		      h_truth_xj->Fill(minit/maxit, event_scale);
		      for (int irange = 0; irange < 3; irange++)
			{
			  if (minit >= ipt_bins[measure_subleading_bin] &&  maxit >= ipt_bins[measure_bins[irange]] && maxit < ipt_bins[measure_bins[irange+1]])
			    {
			      h_truth_xj_range[irange]->Fill(minit/maxit, event_scale);
			      break;
			    }
			}
		      h_linear_truth_xj->Fill(minit/maxit, event_scale);
		    }

		  h_pt1pt2->Fill(es1, es2, event_scale);
		  h_pt1pt2->Fill(es2, es1, event_scale);

		  h_e1e2->Fill(e1, e2, event_scale);
		  h_e1e2->Fill(e2, e1, event_scale);

		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_1, event_scale_new);
		  h_flat_truth_pt1pt2_scaled_new->Fill(total_bin_truth_2, event_scale_new);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_1, event_scale);
		  h_flat_truth_pt1pt2_raw->Fill(total_bin_truth_2, event_scale);
		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_2, event_scale_fake);		  
		}
	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(mysmearjets2.at(0).eta, mysmearjets2.at(1).eta, event_scale_fake);
		  h_lead_pt_v_emfrac->Fill(mysmearjets.at(0).pt, mysmearjets.at(0).emcal, event_scale_fake); 
		  h_mbd_vertex->Fill(mbd_vertex_z[isample], event_scale_fake);
		  h_njets->Fill(nnrecojets, event_scale_fake);	      
		}
	      continue;
	    }
	  // only reco fake
	  if (fill_fake_miss == 3)
	    {
	      
	      if (fill_response)
		{
		  h_pt1pt2_reco_from_pileup->Fill(maxi, mini, (isample != iisample ? 1 : 0), event_scale_fake);
		  h_pt1pt2_reco_from_pileup->Fill(mini, maxi, (isample != iisample ? 1 : 0), event_scale_fake);

		  if (useFakes)
		    {
		      rooResponse.Fake(total_bin_reco_1, event_scale_fake);
		      rooResponse.Fake(total_bin_reco_2, event_scale_fake);
		      h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_1, event_scale_fake);
		      h_flat_reco_to_response_pt1pt2->Fill(total_bin_reco_2, event_scale_fake);
		    }
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2_sample[iisample]->Fill(total_bin_reco_2);	      

		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_all_pt1pt2_sample[iisample]->Fill(total_bin_reco_2, event_scale_fake);

		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_1);
		  h_count_flat_reco_pt1pt2->Fill(total_bin_reco_2);	      

		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_all_pt1pt2->Fill(total_bin_reco_2, event_scale_fake);

		}
	      if (fill_unfold)
		{

		  h_pt1pt2->Fill(es1, es2, event_scale_fake);
		  h_pt1pt2->Fill(es2, es1, event_scale_fake);

		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_1, event_scale_fake);
		  h_flat_reco_pt1pt2_raw->Fill(total_bin_reco_2, event_scale_fake);

		}


	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(mysmearjets2.at(0).eta, mysmearjets2.at(1).eta, event_scale_fake);
		  h_lead_pt_v_emfrac->Fill(mysmearjets2.at(0).pt, mysmearjets2.at(0).emcal, event_scale_fake); 
		  h_mbd_vertex->Fill(mbd_vertex_z[isample], event_scale_fake);
		  h_njets->Fill(nnrecojets, event_scale_fake);	      
		}
	      continue;
	    } 

	}
    }


  
  if (get_mapping == 2)
    {
      // now i need to add samples up and make sure that there are 100 dijets total in all samples

      TH1D *h_reweight_denom = (TH1D*) h_flat_truth_pt1pt2_raw->Clone("h_reweight_denom");
      TH1D *h_reweight_numer = (TH1D*) h_flat_truth_pt1pt2_scaled_new->Clone("h_reweight_numer");
      h_reweight_denom->Scale(1./h_reweight_denom->Integral());
      h_reweight_numer->Scale(1./h_reweight_numer->Integral());
      TH1D *h_reweight_truth_pt1pt2 = (TH1D*) h_reweight_numer->Clone("h_reweight_truth_pt1pt2");
      h_reweight_truth_pt1pt2->Divide(h_reweight_denom);
      
      TString mappingpath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
      
      mappingpath += "_REWEIGHT";
    
      mappingpath += "_" + sys_name;
      
      mappingpath += ".root";

      TFile *fr = new TFile(mappingpath.Data(),"recreate");

      h_reweight_truth_pt1pt2->Write();
      h_flat_truth_pt1pt2_raw->Write();
      h_flat_truth_pt1pt2_scaled_new->Write();
      fr->Close();

      return 0;
      
    }
  if (get_mapping == 1)
    {
      // go through and make a mapping of where in my samples we fill to make a response matrix
      // we need a histogram of truth and reco that are trimmed
      TH1D *h_flat_truth_mapping[nsamples];
      TH1D *h_flat_reco_mapping[nsamples];

      TH1D *h_flat_truth_mapping_all_samples = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
      TH1D *h_flat_reco_mapping_all_samples = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
      TH1D *h_flat_truth_count_all_samples = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
      TH1D *h_flat_reco_count_all_samples = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
      
      h_flat_truth_mapping_all_samples->SetName("h_flat_truth_mapping_all_samples");
      h_flat_truth_mapping_all_samples->Reset();
      
      h_flat_truth_count_all_samples->SetName("h_flat_truth_count_all_samples");
      h_flat_truth_count_all_samples->Reset();

      h_flat_reco_mapping_all_samples->SetName("h_flat_reco_mapping_all_samples");
      h_flat_reco_mapping_all_samples->Reset();
      
      h_flat_reco_count_all_samples->SetName("h_flat_reco_count_all_samples");
      h_flat_reco_count_all_samples->Reset();

      TH2D *h_flat_response_mapping[nsamples];


      for (int i = 0; i < nsamples;i++)
	{
	  double count_response_full = 0;
	  double count_response_keep = 0;
	  h_flat_truth_mapping[i] = (TH1D*) h_flat_truth_to_response_pt1pt2_sample[i]->Clone();
	  h_flat_reco_mapping[i] = (TH1D*) h_flat_reco_to_response_pt1pt2_sample[i]->Clone();
	  h_flat_response_mapping[i] = (TH2D*) h_flat_response_pt1pt2_sample[i]->Clone();
	
	  h_flat_truth_mapping[i]->SetName(Form("h_flat_truth_mapping_%d", i));
	  h_flat_reco_mapping[i]->SetName(Form("h_flat_reco_mapping_%d", i));
	  h_flat_response_mapping[i]->SetName(Form("h_flat_response_mapping_%d", i));
	
	  h_flat_truth_mapping[i]->Reset();
	  h_flat_reco_mapping[i]->Reset();
	  h_flat_response_mapping[i]->Reset();


	  for (int ibr = 0; ibr < nbins_pt_2; ibr++)
	    {
	      for (int ibt = 0; ibt < nbins_pt_2; ibt++)
		{
		  
		  float content = h_count_flat_response_pt1pt2_sample[i]->GetBinContent(ibr+1, ibt+1);
		  count_response_full += content;
		  if (content < 5 && do_trimming)
		    {
		      h_count_flat_reco_pt1pt2_sample[i]->SetBinContent(ibr+1, h_count_flat_reco_pt1pt2_sample[i]->GetBinContent(ibr+1) - content);
		      h_count_flat_truth_pt1pt2_sample[i]->SetBinContent(ibt+1, h_count_flat_truth_pt1pt2_sample[i]->GetBinContent(ibt+1) - content);
		      h_flat_response_mapping[i]->SetBinContent(ibr+1, ibt+1, 0);
		    }
		  else
		    {
		      count_response_keep += content;
		      h_flat_response_mapping[i]->SetBinContent(ibr+1, ibt+1, 1);
		    }
		}
	    }

	  std::cout << "Sample: " << i << std::endl;
	  std::cout << "Keep: " << count_response_keep << " / " << count_response_full << " = " << count_response_keep/count_response_full << std::endl;

	  h_flat_truth_count_all_samples->Add(h_count_flat_truth_pt1pt2_sample[i]);
	  h_flat_reco_count_all_samples->Add(h_count_flat_reco_pt1pt2_sample[i]);
	}

      // now i need to add samples up and make sure that there are 100 dijets total in all samples

      TString mappingpath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
      
      mappingpath += "_MAPPING";
    
      mappingpath += "_" + sys_name;
      
      mappingpath += ".root";

      TFile *fr = new TFile(mappingpath.Data(),"recreate");

      TNtuple *tn_mapping = new TNtuple("tn_mapping","mapping for response construction", "ptbin:use_truth:mapbin_truth:use_reco:mapbin_reco");

      int nempty_reco = 0;
      std::vector<int> binnumbers_reco{};
      int nempty_truth = 0;
      std::vector<int> binnumbers_truth{};
      int nrecobins = 0;
      int ntruthbins = 0;
      for (int ib = 0; ib < nbins_pt_2; ib++)
	{
	  float content_truth = h_flat_truth_count_all_samples->GetBinContent(ib+1);
	  float content_reco = h_flat_reco_count_all_samples->GetBinContent(ib+1);
	  //float purity = he_dijet_fake_binned->GetEfficiency(ib+1);

	  bool use_truth = false;//(content_truth < minentries);
	  bool use_reco = false;//(content_reco < minentries);		
	  if (do_trimming)
	    {
	      use_truth = (content_truth < minentries);
	    }
	  float ut = 0;
	  float ur = 0;
	  float bt = 0;
	  float br = 0;
	  if (use_truth)
	    {
	      ut = 0;
	      bt = 0;
	      nempty_truth++;
	      binnumbers_truth.push_back(0);
	      h_flat_truth_mapping_all_samples->SetBinContent(ib+1, 0);
	    }
	  else
	    {
	      ut = 1;
	      bt = ntruthbins;
	      ntruthbins++;
	      binnumbers_truth.push_back(ntruthbins);
	      h_flat_truth_mapping_all_samples->SetBinContent(ib+1, ntruthbins);
	    }
	  if (use_reco)
	    {
	      ur = 1;
	      br = nrecobins;
	      nrecobins++;
	      //nempty_reco++;
	      binnumbers_reco.push_back(nrecobins);
	      h_flat_reco_mapping_all_samples->SetBinContent(ib+1, nrecobins);
	    }
	  else
	    {
	      ur = 1;
	      br = nrecobins;
	      nrecobins++;
	      binnumbers_reco.push_back(nrecobins);
	      h_flat_reco_mapping_all_samples->SetBinContent(ib+1, nrecobins);
	    }
	  tn_mapping->Fill(ib, ut, bt, ur, br);
	}

      
      // now i put the response mappings and the truth mappings into the mapping file
      tn_mapping->Write();
      h_flat_truth_mapping_all_samples->Write();
      h_flat_reco_mapping_all_samples->Write();
      h_count_flat_truth_pt1pt2->Write();
      h_count_flat_reco_pt1pt2->Write();

      for (int i = 0; i < nsamples; i++)
	{
	  h_flat_response_mapping[i]->Write();
	  h_count_flat_truth_pt1pt2_sample[i]->Write();
	  h_count_flat_reco_pt1pt2_sample[i]->Write();
	}
      fr->Write();
      fr->Close();

      return 0;
            
    }

  h_flat_reco_pt1pt2_raw->Scale(.5);
  h_flat_truth_pt1pt2_raw->Scale(.5);

  TH1D *h_dijet_fake_rate = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
  h_dijet_fake_rate->SetName("h_dijet_fake_rate");
  h_dijet_fake_rate->Divide(h_flat_reco_all_pt1pt2);
  
  TH1D *h_flat_reco_pt1pt2_corr = (TH1D*) h_flat_reco_pt1pt2_raw->Clone("h_flat_reco_pt1pt2_corr");
  if (!useFakes)
    {
      h_flat_reco_pt1pt2_corr->Multiply(h_dijet_fake_rate);
    }

  TH1D *h_flat_truth_pt1pt2 = (TH1D*) h_flat_reco_pt1pt2->Clone();
  h_flat_truth_pt1pt2->Reset();
  h_flat_truth_pt1pt2->SetName("h_truth_flat_pt1pt2");

  for (int ib = 0; ib < nbins_pt_2; ib++)
    {
      int bin = mapped_pt_bin_truth[ib];
      if (bin >= 0)
	{
	  h_flat_truth_pt1pt2->SetBinContent(ib+1, h_flat_truth_pt1pt2_raw->GetBinContent(bin+1));
	  h_flat_truth_pt1pt2->SetBinError(ib+1, h_flat_truth_pt1pt2_raw->GetBinError(bin+1));
	}
    }
  for (int ib = 0; ib < nbins_pt_2; ib++)
    {
      int bin = mapped_pt_bin_reco[ib];
      if (bin >= 0)
	{
	  h_flat_reco_pt1pt2->SetBinContent(ib+1, h_flat_reco_pt1pt2_raw->GetBinContent(bin+1));
	  h_flat_reco_pt1pt2->SetBinError(ib+1, h_flat_reco_pt1pt2_raw->GetBinError(bin+1));
	}
    }
  
  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  TH1D* h_flat_refold_pt1pt2[niterations];
  TH1D* h_flat_truth_fold_pt1pt2;
  int niter = 3;


  h_flat_truth_fold_pt1pt2 = (TH1D*) rooResponse.ApplyToTruth(h_flat_truth_pt1pt2_raw, "hTruth_Folded");
  h_flat_truth_fold_pt1pt2->SetName("h_flat_truth_fold_pt1pt2");
  if (useFakes)
    {
      h_flat_truth_fold_pt1pt2->Add(rooResponse.Hfakes());
    }
  else
    {
      // for (int ib = 0; ib < nbins_pt; ib++)
      // 	{
      // 	  for (int ib2 = 0; ib2 < nbins_pt; ib2++)
      // 	    {
      // 	      int ibbf = 1 + ib*nbins_pt + ib2;
      // 	      int ibb = he_dijet_fake_binned->GetGlobalBin(ib+1, ib2+1);
      // 	      double value = he_dijet_fake_binned->GetEfficiency(ibb);
      // 	      h_flat_truth_fold_pt1pt2->SetBinContent(ibbf, h_flat_truth_fold_pt1pt2->GetBinContent(ibbf)/value);
      // 	    }
      // 	}
    }
  for (int iter = 0; iter < niterations; iter++ )
    {
      
      RooUnfoldBayes   unfold (&rooResponse, h_flat_reco_pt1pt2_corr, iter + 1, false, useFakes);    // OR

      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_reco_pt1pt2->Clone();
      h_flat_unfold_pt1pt2[iter]->Reset();
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
      for (int ib = 0; ib < nbins_pt_2; ib++)
	{
	  int bin = mapped_pt_bin_truth[ib];
	  if (bin >= 0)
	    {
	      h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin+1));
	      h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin+1));
	    }
	}

      h_flat_refold_pt1pt2[iter] = (TH1D*) rooResponse.ApplyToTruth(h_flat_unfold_skim[iter], "hRefolded");
      h_flat_refold_pt1pt2[iter]->SetName(Form("h_flat_refold_pt1pt2_iter%d", iter));
      if (useFakes)
	{
	  h_flat_refold_pt1pt2[iter]->Add(rooResponse.Hfakes());
	}
      else
	{
	  // for (int ib = 0; ib < nbins_pt; ib++)
	  //   {
	  //     for (int ib2 = 0; ib2 < nbins_pt; ib2++)
	  // 	{
	  // 	  int ibbf = 1 + ib*nbins_pt + ib2;
	  // 	  int ibb = he_dijet_fake_binned->GetGlobalBin(ib+1, ib2+1);
	  // 	  double value = he_dijet_fake_binned->GetEfficiency(ibb);
	  // 	  h_flat_refold_pt1pt2[iter]->SetBinContent(ibbf, h_flat_refold_pt1pt2[iter]->GetBinContent(ibbf)/value);
	  // 	}
	  //   }
	}
    }
  
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_truth_fold = new TH2D("h_pt1pt2_truth_fold", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];
  TH2D *h_pt1pt2_refold[niterations];

  TH2D *h_a_pt1pt2_reco = new TH2D("h_a_pt1pt2_reco", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_a_pt1pt2_truth_fold = new TH2D("h_a_pt1pt2_truth_fold", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_a_pt1pt2_truth = new TH2D("h_a_pt1pt2_truth", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_a_pt1pt2_unfold[niterations];
  TH2D *h_a_pt1pt2_refold[niterations];
      
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
      h_pt1pt2_refold[iter] = new TH2D("h_pt1pt2_refold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_pt1pt2_refold[iter]->SetName(Form("h_pt1pt2_refold_iter%d", iter));
	  
      h_a_pt1pt2_unfold[iter] = new TH2D("h_a_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_a_pt1pt2_unfold[iter]->SetName(Form("h_a_pt1pt2_unfold_iter%d", iter));
      h_a_pt1pt2_refold[iter] = new TH2D("h_a_pt1pt2_refold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_a_pt1pt2_refold[iter]->SetName(Form("h_a_pt1pt2_refold_iter%d", iter));
	  
    }
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_truth_fold = new TH1D("h_xj_truth_fold", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_truth_direct = new TH1D("h_xj_truth_direct", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  TH1D *h_xj_refold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
      h_xj_refold[iter] = new TH1D(Form("h_xj_refold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
      
  h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");
      
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_truth_fold_pt1pt2, h_pt1pt2_truth_fold, nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins_pt);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins_pt);
      histo_opps::make_sym_pt1pt2(h_flat_refold_pt1pt2[iter], h_pt1pt2_refold[iter], nbins_pt);
    }
      
	  
  histo_opps::make_asym_pt1pt2(h_a_pt1pt2_truth, h_pt1pt2_truth, nbins_pt);
  histo_opps::make_asym_pt1pt2(h_a_pt1pt2_truth_fold, h_pt1pt2_truth_fold, nbins_pt);
  histo_opps::make_asym_pt1pt2(h_a_pt1pt2_reco, h_pt1pt2_reco, nbins_pt);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_asym_pt1pt2(h_a_pt1pt2_unfold[iter], h_pt1pt2_unfold[iter], nbins_pt);
      histo_opps::make_asym_pt1pt2(h_a_pt1pt2_refold[iter], h_pt1pt2_refold[iter], nbins_pt);
    }
      
  gStyle->SetPaintTextFormat("4.0f");

  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 1000, 250);
  cpt1pt2->Divide(4, 1);
  cpt1pt2->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  h_pt1pt2_reco->Draw("colz");
  
  cpt1pt2->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_truth->Draw("colz");
  cpt1pt2->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_unfold[niter]->Draw("colz");

  cpt1pt2->cd(4);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_a_pt1pt2_unfold[niter]->Draw("colz");

  cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2sim_%s_r%02d_%s%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, (primer? "PRIMER" + std::to_string(primer) : "").c_str() , sys_name.c_str()));

  h_e1e2->Scale(0.5);
  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth_fold, h_xj_truth_fold, nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_e1e2, h_xj_truth_direct, nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);

  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_refold[iter], h_xj_refold[iter], nbins_pt, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
    }

  TH1D *h_xj_truth_project[3];
  for (int irange = 0; irange < 3; irange++)
    {
      h_xj_truth_project[irange] = (TH1D*) h_xj_truth_direct->Clone();
      h_xj_truth_project[irange]->SetName(Form("h_xj_truth_project_%d", irange));
      h_xj_truth_project[irange]->Reset();
      histo_opps::project_xj(h_e1e2, h_xj_truth_project[irange], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins);
	      
    }
  TCanvas *cxj_split = new TCanvas("cxj_split","cxj_split", 1800, 600);
  dlutility::systematic_split_canvas(cxj_split, 3);

  for (int irange = 0; irange < 3; irange++)
    {
      dlutility::SetMarkerAtt(h_xj_truth_project[irange], kBlue, 1, 20);
      dlutility::SetLineAtt(h_xj_truth_project[irange], kBlue, 1, 1);
      dlutility::SetMarkerAtt(h_truth_xj_range[irange], kBlack, 1, 20);
      dlutility::SetLineAtt(h_truth_xj_range[irange], kBlack, 1, 1);
      cxj_split->cd(1 + irange*2);
      h_truth_xj_range[irange]->Draw();
      h_xj_truth_project[irange]->Draw("same");
    }
  cxj_split->cd(7);
  dlutility::DrawSPHENIXpp(0.22, 0.84, 0.08);
  dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74, 0, kBlack, 0.08);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59, 0, kBlack, 0.08);
  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49, 0, kBlack, 0.08);
  TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.07);
  leg->SetTextFont(42);
  leg->AddEntry(h_xj_truth_project[0], "Projected");
  leg->AddEntry(h_truth_xj_range[0], "Filled");
  leg->Draw("same");

  cxj_split->Print(Form("%s/unfolding_plots/proj_compare_range_%s_r%02d_%s%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, (primer? "PRIMER" + std::to_string(primer) : "").c_str() , sys_name.c_str()));
	  
  TCanvas *cjetdiv = new TCanvas("cjetdiv","cjetdiv", 700, 500);
  dlutility::createCutCanvas(cjetdiv);
  cjetdiv->cd(1);
  gPad->SetLogy();

  dlutility::SetLineAtt(h_truth_lead, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead, kBlack, 0.8, 24);

  dlutility::SetLineAtt(h_truth_lead_sample[0], kRed, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[0], kRed, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[1], kGreen, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[1], kGreen, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[2], kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[2], kBlue, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[3], kOrange, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[3], kOrange, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[4], kViolet, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[4], kViolet, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[5], kPink, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[5], kPink, 1, 20);
  dlutility::SetLineAtt(h_truth_lead_sample[6], kGreen+1, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead_sample[6], kGreen+1, 1, 20);

  h_truth_lead->SetTitle(";Leading Jet p_{T} [GeV];counts * lumiscale");
  h_truth_lead->SetMinimum(1);
  h_truth_lead->Draw("p");
  h_truth_lead_sample[0]->Draw("same p");
  h_truth_lead_sample[1]->Draw("same p");
  h_truth_lead_sample[2]->Draw("same p");
  h_truth_lead_sample[3]->Draw("same p");
  h_truth_lead_sample[4]->Draw("same p");
  h_truth_lead_sample[5]->Draw("same p");
  h_truth_lead_sample[6]->Draw("same p");
  h_truth_lead->Draw("p same");


  cjetdiv->cd(2);
  if (!use_herwig)
    {
      dlutility::DrawSPHENIXppsize(0.05, 0.9, 0.08, 1, 0, 1, "Pythia8");
    }
  else if (use_herwig == 1)
    {
      dlutility::DrawSPHENIXppsize(0.05, 0.9, 0.08, 1, 0, 1, "HERWIG 7.8");
    }
  dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.05, 0.8, 0, kBlack, 0.08);
  dlutility::drawText("All Dijets", 0.05, 0.75, 0, kBlack, 0.08);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.05, 0.7, 0, kBlack, 0.08);
  dlutility::drawText(Form("|#eta| < %2.1f", 1.1 - 0.1*cone_size), 0.05, 0.65, 0, kBlack, 0.08);

  TLegend *leg3 = new TLegend(0.01, 0.16, 0.7, 0.63);
  leg3->SetLineWidth(0);
  leg3->AddEntry(h_truth_lead, "Combined Sample");
  leg3->AddEntry(h_truth_lead_sample[0], "8 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[1], "12 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[2], "20 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[3], "30 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[4], "40 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[5], "50 GeV Sample");
  leg3->AddEntry(h_truth_lead_sample[6], "60 GeV Sample");
  leg3->SetTextSize(0.08);
  leg3->Draw("same");

  cjetdiv->Print(Form("%s/unfolding_plots/combined_sample_%s_r%02d_%s%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, (primer?Form("PRIMER%d", primer):""), sys_name.c_str()));


  TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
	  
  if (primer > 0)
    {
      responsepath += "_PRIMER" + std::to_string(primer);;
    }

  responsepath += "_" + sys_name;

  responsepath += ".root";
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  rooResponse.Write();
  h_dijet_fake_rate->Write();
  h_count_flat_truth_pt1pt2->Write();
  h_count_flat_reco_pt1pt2->Write();
  if (h2_flatreweight_pt1pt2)
    {
      h2_flatreweight_pt1pt2->Write();
    }
  if (h2_flatreweight_pt1pt2)
    {
      h2_flatreweight_pt1pt2->Write();
    }
  for (int i = 0; i < 6; i++)
    {
      h_count_flat_truth_pt1pt2_sample[i]->Write();
      h_count_flat_reco_pt1pt2_sample[i]->Write();
      h_count_flat_response_pt1pt2_sample[i]->Write();
    }
  
  h_pt1pt2_reco_from_pileup->Write();
  h_pt1pt2_truth_from_pileup->Write();

  he_dijet_subleading_binned->Write();
  he_dijet_match_binned->Write();
  he_dijet_trigger_binned->Write();
  he_dijet_truth_bad_binned->Write();
  he_dijet_truth_good_binned->Write();
  h_lead_pt_v_emfrac->Write();
  h_eta_lead_sublead->Write();
  h_eta_phi_lead->Write();
  h_flat_reco_pt1pt2_raw->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_reco_to_response_pt1pt2->Write();
  h_flat_reco_all_pt1pt2->Write();
  h_flat_reco_no_fakes_pt1pt2->Write();
  h_flat_truth_pt1pt2_raw->Write();
  h_flat_truth_pt1pt2->Write();
  h_flat_truth_to_response_pt1pt2->Write();
  h_flat_truth_fold_pt1pt2->Write();
  h_flat_response_pt1pt2->Write();
  h_mbd_vertex->Write();
  h_njets->Write();
  h_pt1pt2->Write();
  h_e1e2->Write();
  // Full closure drawing
  h_linear_truth_xj->Write();
  h_xj_truth->Write();
  h_xj_truth_fold->Write();
  h_xj_truth_direct->Write();
  h_truth_xj->Write();
  for (int irange = 0; irange < 3; irange++)
    {
      h_truth_xj_range[irange]->Write();
    }
	

  he_dijet_mbd_hit->Write();
  he_dijet_mbd_hit_binned->Write();
  he_dijet_fake_binned->Write();
  he_dijet_miss_binned->Write();
  he_dijet_trigger->Write();
  he_dijet_reco->Write();
  he_dijet_match->Write();
  he_dijet_match_reco->Write();
  he_dijet_vertex->Write();
  he_dijet_vertex_truth->Write();
	  

  for (int iter = 0 ; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter]->Write();
      h_xj_unfold[iter]->Write();
      h_flat_refold_pt1pt2[iter]->Write();
      h_xj_refold[iter]->Write();
    }
  h_truth_lead->Write();
  h_truth_sublead->Write();
  h_reco_lead->Write();
  h_reco_sublead->Write();

  h_truth_lead_pt_phi->Write();
  h_truth_sublead_pt_phi->Write();
  h_reco_lead_pt_phi->Write();
  h_reco_sublead_pt_phi->Write();

  h_truth_lead_pt_eta->Write();
  h_truth_sublead_pt_eta->Write();
  h_reco_lead_pt_eta->Write();
  h_reco_sublead_pt_eta->Write();

  h_xj_reco->Write();
	      
  fr->Close();


  return 0;
}

int main(int argc, char *argv[])
{

  std::string config = "binning.config";
  int full_or_half = 0;
  int niterations = 10;
  int cone_size = 4;
  int primer = 0;
  int set=0;
  int useFakes = 1;
  int useEfficiencies = 1;
  int input_generator = 0;
  int unfold_generator = 0;
  int get_mapping = 0;
  int do_trimming = 1;
  for (int i = 1; i < argc; ++i)
    {
      std::string arg = argv[i];

      if (arg == "-n" && i + 1 < argc)
	{
	  set++;
	  niterations = std::stoi(argv[++i]);  // Convert next argument to int
	}
      else if (arg == "-c" && i + 1 < argc)
	{
	  set++;
	  config = argv[++i];  // Next argument as string
	}
      else if (arg == "-r" && i + 1 < argc)
	{
	  set++;
	  cone_size = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-p" && i + 1 < argc)
	{
	  set++;
	  primer = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-h" && i < argc)
	{
	  set++;
	  full_or_half = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-v" && i+1 < argc)
	{
	  verbosity = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-f" && i+1 < argc)
	{
	  useFakes = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-e" && i+1 < argc)
	{
	  useEfficiencies = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-ig" && i+1 < argc)
	{
	  input_generator = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-ug" && i+1 < argc)
	{
	  unfold_generator = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-m" && i+1 < argc)
	{
	  get_mapping = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-t" && i+1 < argc)
	{
	  do_trimming = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else
	{
	  std::cerr << "Unknown or incomplete argument: " << arg << "\n";
	  return 1;
	}
    }
  if (set < 5)
    {
      std::cout << "Not enough settings: " << std::endl;
      std::cout << "[usage] : " << std::endl;
      std::cout << "    ./createResponse_noempty_pp -c binning.config -r 4 -n 10 -p 1 -h 0 -f 0 -e 0 -t 1 -ig 0 -ug 0" << std::endl;
      return 1;
    }
  
  return createResponse_noempty_pp(config, full_or_half, niterations, cone_size, primer, useFakes, useEfficiencies, input_generator, unfold_generator, get_mapping, do_trimming);
  
}
