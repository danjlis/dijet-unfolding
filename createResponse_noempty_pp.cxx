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
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
static int verbosity = 0;

const int separate_efficiency_correction = 0;

int createResponse_noempty_pp(const std::string configfile = "binning.config", const int full_or_half = 0, const int niterations = 10, const int cone_size = 4,  const int primer = 0, const int useFakes = 1, const int useEfficiencies = 1, const int input_generator = 0, const int unfold_generator = 0)
{
  gRandom->SetSeed(0);
  
  std::cout << "Using fakes : " << useFakes << std::endl;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  float boundary_r4[6];
  boundary_r4[0] = 14;
  boundary_r4[1] = 22;
  boundary_r4[2] = 33;
  boundary_r4[3] = 42;
  boundary_r4[4] = 52;
  boundary_r4[5] = 100;

  float maxreco_r4[6];
  maxreco_r4[0] = 32;
  maxreco_r4[1] = 50;
  maxreco_r4[2] = 60;
  maxreco_r4[3] = 100;
  maxreco_r4[4] = 100;
  // float maxreco_r4[6];
  // maxreco_r4[0] = 39;
  // maxreco_r4[1] = 59;
  // maxreco_r4[2] = 72;
  // maxreco_r4[3] = 100;
  // maxreco_r4[4] = 100;

  read_binning rb(configfile.c_str());

  dijetfinder djf(cone_size);
  djf.SetVerbosity(verbosity);
  
  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();
  //Int_t primer = rb.get_primer();

  //Double_t dphicut = rb.get_dphicut();

  Double_t vtx_cut = rb.get_vtx_cut();
  Double_t njet_cut = rb.get_njet_cut();

  Int_t zyam_sys = rb.get_zyam_sys();
  Int_t inclusive_sys = rb.get_inclusive_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  Int_t prior_sys = rb.get_prior_sys();
  Int_t use_herwig = rb.get_herwig();
  
  TF1 *finjes = nullptr;

  TFile *infilejes = new TFile(Form("%s/r%02d/jes/jes_fits_r%02d_0JES_PYTHIA.root", rb.get_jesr_location().c_str(), cone_size, cone_size), "r");
  finjes = (TF1*) infilejes->Get("fjesprime_2");

  TF1 *gjer = nullptr;
  TH1D *hjersmear = nullptr;

  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TFile *finjer = new TFile(Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure_PYTHIA.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)), "r");
   
  std::string sys_name = "nominal";
  std::string calib_string = "SMEAR";

  if (prior_sys)
    {
      sys_name = "PRIOR";
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
  if (JER_sys != 0)
    {
      if (JER_sys < 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_down_smear");
	  gjer = (TF1*) finjer->Get("jerneg");
	  calib_string = "SMEAR_DOWN";
	  sys_name = "negJER";
	}
      else if (JER_sys > 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_up_smear");
	  gjer = (TF1*) finjer->Get("jerpos");
	  calib_string = "SMEAR_UP";
	  sys_name = "posJER";
	}
    }
  // int nbins_smear = 0;
  // //sWrite
  // if (!hjersmear)
  //   {
  //     std::cout << "no histo for smear in " << Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)) << std::endl;
  //     return 1;
  //   }

  // nbins_smear = hjersmear->GetNbinsX();
  
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
  
  std::string j_file[5];
  j_file[0] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_12_ana509_MDC2-00000028.root";
  j_file[1] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
  j_file[2] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
  j_file[3] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_40_ana509_MDC2-00000028.root";
  j_file[4] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_50_ana509_MDC2-00000028.root";

  
  bool use_sample[5] = {1, 1, 1, 1, 1};
  if (use_herwig)
    {
      use_sample[1] = 0;
      use_sample[2] = 0;
      use_sample[4] = 0;      
    }
  
  float n_events[5];
  TFile *finsim[5];


  TTree *ttree[5];
  ULong64_t gl1_scaled[5];

  std::vector<float> *truth_jet_pt_ref[5] = {0};
  std::vector<float> *truth_jet_pt[5] = {0};
  std::vector<float> *truth_jet_eta[5] = {0};
  std::vector<float> *truth_jet_phi[5] = {0};
  
  std::vector<float> *reco_jet_pt[5] = {0};
  std::vector<float> *reco_jet_pt_uncalib[5] = {0};
  std::vector<float> *reco_jet_emcal[5] = {0};
  std::vector<float> *reco_jet_e[5] = {0};
  std::vector<float> *reco_jet_eta[5] = {0};
  std::vector<float> *reco_jet_eta_det[5] = {0};
  std::vector<float> *reco_jet_phi[5] = {0};

  float truth_vertex_z[5];
  float mbd_vertex_z[5];
  int mbd_hit[5];

  for (int j = 0; j < 5; j++)
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

  
  float cs_12 = (1.49e6);
  float cs_20 = (6.26e4);
  float cs_30 = (2.53e3);
  float cs_40 = (1.355e2);
  float cs_50 = (7.311);

  float csh_10 = 1.5703e-9;
  float csh_30 = 1.473e-12;
  float scale_factor[5];
  scale_factor[0] = cs_12/cs_50;
  scale_factor[1] = cs_20/cs_50;
  scale_factor[2] = cs_30/cs_50;
  scale_factor[3] = cs_40/cs_50; 
  scale_factor[4] = 1;

  if (use_herwig)
    {
      scale_factor[0] = 1.4216108 * ((float) n_events[3]/ ((float)n_events[0])) * csh_10/(csh_30 * 4.);
      scale_factor[3] = 1;
    }
  
  int prior_iteration = 2;

  if (primer == 0 && !unfold_generator && false)
    {
      std::cout << "opening file for prior" << std::endl;
      TFile *fin_iter = new TFile(Form("unfolding_hists/iteration_tune_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");      

      TH1D *hunc = (TH1D*) fin_iter->Get("h_total_uncertainties");
      int minimum_iter = 0;
      float mini = hunc->GetBinContent(2);;
      for (int ib = 3; ib < 10; ib++)
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
  if (primer != 1)
    {
      TFile *fvtx = new TFile(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()),"r");
      TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");
      for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
	{
	  vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
	}
      h_eta_reweight = (TH2D*) fvtx->Get("h_eta_reweight");
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

  if (use_herwig)
    {
      boundary_r4[0] = truth_leading_cut;
      boundary_r4[1] = 30;
      boundary_r4[3] = 30;
      boundary_r4[4] = ipt_bins[nbins + 1];;     
    }

	
  djf.setTruthCuts(truth_leading_cut, truth_subleading_cut);
  djf.setRecoCuts(reco_leading_cut, reco_subleading_cut);

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

  TH2D *h_reco_fakes = new TH2D("h_reco_fakes","", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  TH2D *h_no_truth_fakes = new TH2D("h_no_truth_fakes","", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_match_fakes = new TH2D("h_match_fakes","", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  
  TEfficiency *he_dijet_fake_binned = new TEfficiency("he_dijet_fake_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Fake Rate", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TEfficiency *he_dijet_miss_binned = new TEfficiency("he_dijet_miss_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Miss Rate", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TEfficiency *he_dijet_mbd_hit_binned = new TEfficiency("he_dijet_mbd_hit_binned",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", nbins_pt, 0, nbins_pt);
  TEfficiency *he_dijet = new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_mbd_hit = new TEfficiency("he_dijet_mbd_hit",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_reco = new TEfficiency("he_dijet_reco",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_match = new TEfficiency("he_dijet_match",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_match_reco = new TEfficiency("he_dijet_match_reco",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_fake = new TEfficiency("he_dijet_fake",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_trigger = new TEfficiency("he_dijet_trigger",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_vertex = new TEfficiency("he_dijet_vertex",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);
  TEfficiency *he_dijet_vertex_truth = new TEfficiency("he_dijet_vertex_truth",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100, 50, 0, 100);

  TH1D *h_truth_lead_sample[5];
  for (int i = 0; i < 5; i++)
    {
      h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }
  
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH2D *h_eta_lead_sublead = new TH2D("h_eta_lead_sublead","", 220, -1.1, 1.1, 220, -1.1, 1.1);
  
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

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  TH1D *h_flat_reco_all_pt1pt2 = new TH1D("h_reco_flat_all_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt, nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  h_flat_truth_to_response_pt1pt2->Sumw2();
  h_flat_reco_to_response_pt1pt2->Sumw2();
  h_flat_truth_pt1pt2->Sumw2();
  h_flat_reco_pt1pt2->Sumw2();
  h_flat_response_pt1pt2->Sumw2();
  
  TH2D *h_count_flat_response_pt1pt2 = new TH2D("h_count_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt, nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  TH1D *h_flat_reco_all_pt1pt2_sample[5];// = new TH1D("h_reco_flat_all_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2_sample[5];// new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2_sample[5];// new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2_sample[5];// new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2_sample[5];// new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH2D *h_flat_response_pt1pt2_sample[5];// new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);
  TH2D *h_count_flat_response_pt1pt2_sample[5];// new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  for (int i = 0 ; i < 5; i++)
    {
      h_flat_reco_all_pt1pt2_sample[i] = new TH1D(Form("h_reco_flat_all_pt1pt2_%d", i),";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_count_flat_truth_pt1pt2_sample[i] = new TH1D(Form("h_truth_count_flat_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_count_flat_reco_pt1pt2_sample[i] = new TH1D(Form("h_count_reco_flat_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_flat_truth_to_response_pt1pt2_sample[i] = new TH1D(Form("h_truth_flat_to_response_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_flat_reco_to_response_pt1pt2_sample[i] = new TH1D(Form("h_reco_flat_to_response_pt1pt2_%d", i), ";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_flat_response_pt1pt2_sample[i] = new TH2D(Form("h_flat_response_pt1pt2_%d", i), ";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt, nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      h_count_flat_response_pt1pt2_sample[i] = new TH2D(Form("h_count_flat_response_pt1pt2_%d", i), ";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt, nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
    }


  // for trimming

  TH1D *h_flat_truth_mapping_primer;
  TH1D *h_flat_reco_mapping_primer;
  TH2D *h_count_flat_response_pt1pt2_sample_primer[5];
  TH1D *h_count_flat_reco_pt1pt2_sample_primer[5];
  TH1D *h_count_flat_truth_pt1pt2_sample_primer[5];
  if (primer != 1)
    {
      TString responseppath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);

      responseppath += "_PRIMER1_" + sys_name;
      
      responseppath += ".root";

  
      TFile *fresponse = new TFile(responseppath.Data(),"r");

      h_flat_truth_mapping_primer = (TH1D*) fresponse->Get("h_flat_truth_mapping"); 
      if (!h_flat_truth_mapping_primer)
	{
	  std::cout << "no truth" << std::endl;
	  return 1;
	}
      h_flat_truth_mapping_primer->SetName("h_flat_truth_mapping_primer");

      h_flat_reco_mapping_primer = (TH1D*) fresponse->Get("h_flat_reco_mapping"); 
      if (!h_flat_reco_mapping_primer)
	{
	  std::cout << "no reco" << std::endl;
	  return 1;
	}
      h_flat_reco_mapping_primer->SetName("h_flat_reco_mapping_primer");

      for (int i = 0; i < 5; i++)
	{
	  h_count_flat_truth_pt1pt2_sample_primer[i] = (TH1D*) fresponse->Get(Form("h_truth_count_flat_pt1pt2_%d", i)); 
	  if (!h_count_flat_truth_pt1pt2_sample_primer[i])
	    {
	      std::cout << "no truth" << std::endl;
	      return 1;
	    }
	  h_count_flat_truth_pt1pt2_sample_primer[i]->SetName(Form("h_count_flat_truth_pt1pt2_sample_primer_%d", i));
	  
	  h_count_flat_reco_pt1pt2_sample_primer[i] = (TH1D*) fresponse->Get(Form("h_count_reco_flat_pt1pt2_%d", i)); 
	  if (!h_count_flat_reco_pt1pt2_sample_primer[i])
	    {
	      std::cout << "no reco" << std::endl;
	      return 1;
	    }
	  h_count_flat_reco_pt1pt2_sample_primer[i]->SetName(Form("h_count_flat_reco_pt1pt2_sample_primer_%d", i));
	  
	  h_count_flat_response_pt1pt2_sample_primer[i] = (TH2D*) fresponse->Get(Form("h_count_flat_response_pt1pt2_%d", i)); 
	  if (!h_count_flat_response_pt1pt2_sample_primer[i])
	    {
	      std::cout << "no response" << std::endl;
	      return 1;
	    }
	  h_count_flat_response_pt1pt2_sample_primer[i]->SetName(Form("h_count_flat_response_pt1pt2_sample_primer_%d", i));
	}

    }
  
  // For the prior sensitivity
  TH1D *h_flatreweight_pt1pt2 = new TH1D("h_unfold_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);


  if (!prior_sys && !primer)
    {
      TFile *fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
      //TH1D *h_unfold_flat = (TH1D*) fun->Get("h_data_flat_pt1pt2");
      TH1D *h_unfold_flat = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
      std::cout << "Got histograms" << std::endl;
      TFile *ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");

      TH1D *h_truth_flat = (TH1D*) ftr->Get("h_flat_truth_to_unfold_pt1pt2");
      h_truth_flat->SetName("h_truth_flat");
      std::cout << "Got histograms" << std::endl;
      TH2D *h2_flatreweight_pt1pt2 = new TH2D("h2_flatreweight_pt1pt2",";p_{T,1}^{truth} [GeV]; p_{T,2}^{truth} [GeV] ; Reweight Factor", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

      // of signal region
      
      float integral_of_signal = 0;
      float integral_of_truth = 0;
      for (int ibin = 0; ibin < nbins_pt*nbins_pt; ibin++)
	{
	  int pt1_bin = ibin/nbins_pt;
	  int pt2_bin = ibin%nbins_pt;
	  int max_bin = std::max(pt1_bin, pt2_bin);
	  int min_bin = std::min(pt1_bin, pt2_bin);
	  if (true)//	  if (max_bin >= measure_bins[0] && max_bin < measure_bins[3] && min_bin >= measure_subleading_bin)
	    {
	      integral_of_signal += h_unfold_flat->GetBinContent(ibin+1);
	      integral_of_truth += h_truth_flat->GetBinContent(ibin+1);
	    }	  
	}

      h_unfold_flat->Scale(1./integral_of_signal);
      h_truth_flat->Scale(1./integral_of_truth);

      for (int ibin = 0; ibin < nbins_pt*nbins_pt; ibin++)
	{

	  float pt1_bin = ibin/nbins_pt;
	  float pt2_bin = ibin%nbins_pt;

	  int gbin = h2_flatreweight_pt1pt2->GetBin(pt1_bin+1, pt2_bin+1);

	  float v = h_unfold_flat->GetBinContent(ibin+1);
	  float b = h_truth_flat->GetBinContent(ibin+1);
	  
	  if (b > 0)
	    {
	      float vb = v/b;
	      if (v == 0) vb = 1;
	      h_flatreweight_pt1pt2->SetBinContent(ibin+1, vb);
	      h2_flatreweight_pt1pt2->SetBinContent(gbin, vb);
	    }
	  else
	    {
	      h_flatreweight_pt1pt2->SetBinContent(ibin+1, 1);
	      h2_flatreweight_pt1pt2->SetBinContent(gbin, 1);
	    }

	}

      TCanvas *cre = new TCanvas("cre","cre", 500, 500);
      cre->SetLeftMargin(0.1);
      cre->SetRightMargin(0.19);
      h2_flatreweight_pt1pt2->Draw("colz");
      dlutility::DrawSPHENIX(0.2, 0.87);
      dlutility::drawText("Prior Reweighting Matrix", 0.2, 0.77);
      cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
    }

  int nbin_response = nbins_pt*nbins_pt;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  TRandom *rng = new TRandom(0);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets2;
  std::vector<struct jet> myrecojets2;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  std::vector<std::pair<struct jet, struct jet>> matched_dijets2;

  for (int isample = 0; isample < 5; isample++)
    {
      if (!use_sample[isample]) continue;
      std::cout << "Sample " << isample << std::endl;
      int entries0 = 0;
      int entries2 = ttree[isample]->GetEntries();

      int starting = 1;
      if (verbosity > 5) entries2 = 100000;
      
      for (int i = entries0; i < entries2; i++)
	{
	  ttree[isample]->GetEntry(i);
	  double event_scale = scale_factor[isample];
	  if (verbosity > 5)
	    {
	      std::cout << "------------- Event " << i << " --------------" << std::endl;
	    }
	  
	  bool has_vertex = (fabs(mbd_vertex_z[isample]) < vtx_cut);

	  bool has_mbd_hit = mbd_hit[isample];

	  if (!has_vertex || !has_mbd_hit) continue;	  


	  bool triggered = ( ( ( gl1_scaled[isample] >> 22 ) & 0x1) == 0x1);
	  //if (inrecojets >= njet_cut) continue;
	  // Vertex Rewieghting
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

	  bool fill_unfold = (!full_or_half) || (!fill_response && full_or_half);

	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (cone_size != 4)
	    {
	      maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
	    }

	  float maxptreco = *std::max_element(reco_jet_pt_uncalib[isample]->begin(), reco_jet_pt_uncalib[isample]->end());

	  if (maxpttruth < boundary_r4[isample] || maxpttruth >= boundary_r4[isample+1]) continue;
	  if (maxptreco >= maxreco_r4[isample]) continue;
	  
	  if (primer != 1 && has_vertex)
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
	  // top two of each
	  mytruthjets2.clear();
	  myrecojets2.clear();

	  int ntruthjets = truth_jet_pt[isample]->size();
	  for (int j = 0; j < ntruthjets;j++)
	    {

	      if (truth_jet_pt[isample]->at(j) < truth_subleading_cut) continue;
	      if (fabs(truth_jet_eta[isample]->at(j)) > 0.7) continue;

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

	      if (reco_jet_e[isample]->at(j) < 0) continue;
	      struct jet tempjet;

	      tempjet.istruth = 0;

	      float temppt = reco_jet_pt[isample]->at(j);
	      float tempptreco = temppt;
	      int ib = floor((temppt - 3)/0.1) + 1;
	      float smear1 = 0.0;//hjersmear->GetBinContent(ib );
	      float jersmear = 0;
	      if (smear1 > 0)
		{
		  fsmear->SetParameter(2, smear1);
		  jersmear = fsmear->GetRandom();
		}
	      tempptreco = temppt*(1 + JES_sys) + jersmear*temppt;// + JES_sys * temppt;
	
	      tempjet.pt = tempptreco;
	      tempjet.pt_uncalib = reco_jet_pt[isample]->at(j);
	      if (tempjet.pt > 7) nnrecojets++;
	      if (tempjet.pt < reco_subleading_cut) continue;
	      
	      tempjet.emcal = reco_jet_emcal[isample]->at(j);
	      tempjet.eta = reco_jet_eta[isample]->at(j);
	      tempjet.eta_det = reco_jet_eta_det[isample]->at(j);
	      tempjet.phi = reco_jet_phi[isample]->at(j);
	      tempjet.t = 0;
	      tempjet.id = j;

	      myrecojets.push_back(tempjet);
	    }

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
	  

	  if (verbosity > 5)
	    {
	      std::cout << "----- RECO JETS -----" << std::endl;
	      for (auto rjet : myrecojets)
		{
		  rjet.print();
		}
	      std::cout << "----- TRUTH JETS -----" << std::endl;
	      for (auto tjet : mytruthjets)
		{
		  tjet.print();
		}
	    }


	  bool truth_good = djf.check_dijet_truth(mytruthjets);

	  if (!truth_good && !useFakes)
	    {
	      continue;
	    }

	  if (myrecojets.size() > 1)
	    {
	      int sizet = myrecojets.size();
	      myrecojets2 = {myrecojets.begin(), myrecojets.begin() + sizet};
	      if (verbosity > 5)
		{
		  std::cout << "RECO GOOD" << std::endl;
		}		  
	    }

	  if (mytruthjets.size() > 1)//truth_good)
	    {
	      mytruthjets2 = {mytruthjets.begin(), mytruthjets.begin()+2};
	      if (verbosity > 5)
		{
		  std::cout << "TRUTH GOOD" << std::endl;
		}		  
	    }

	  matched_dijets = djf.match_dijets(myrecojets2, mytruthjets2);

	  std::sort(matched_dijets.begin(), matched_dijets.end(), [] (auto a, auto b) { return a.first.pt > b.first.pt; });

	  bool matched = true;
	  if (matched_dijets.size() > 1)
	    {
	      matched = true;
	      
	      // matched_dijets2 = {matched_dijets.begin(), matched_dijets.begin() + 2};
	      // for (auto mjet : matched_dijets2)
	      // 	{
	      // 	  matched &= mjet.first.id == mytruthjets.at(0).id || mjet.first.id == mytruthjets.at(1).id;
	      // 	  matched &= mjet.first.id == mytruthjets.at(0).id || mjet.first.id == mytruthjets.at(1).id;
	      // 	  matched &= mjet.second.id == myrecojets.at(0).id || mjet.second.id == myrecojets.at(1).id;
	      // 	  matched &= mjet.second.id == myrecojets.at(0).id || mjet.second.id == myrecojets.at(1).id;
	      // 	}
	    }
	  else
	    {
	      matched = false;
	    }
	
	  
	  //	  bool matched = matched_dijets.size() == 2;

	  // smear matched jets
	  if (matched && false && !input_generator)
	    {
	      myrecojets2.clear();
	      
	      for (int ijet = 0 ; ijet < matched_dijets.size(); ijet++)
		{
		  auto mjet = matched_dijets.at(ijet);
		  float temppt = mjet.first.pt;
		  
		  // while (hjersmear->GetBinLowEdge(ib++) < temppt && ib <= nbins_smear)
		  //   {
		  //     continue;
		  //   }

		  float smear1 = 0.08;//gjer->Eval(temppt);//hjersmear->GetBinContent(ib - 1);

		  fsmear->SetParameter(2, smear1);
		  float jersmear = fsmear->GetRandom();
		  float tempptreco = mjet.second.pt + jersmear*temppt + JES_sys * temppt;
		  matched_dijets.at(ijet).second.pt = tempptreco;
		  myrecojets2.push_back(matched_dijets.at(ijet).second);
		}

	      std::sort(myrecojets2.begin(), myrecojets2.end(), [] (auto a, auto b) { return a.pt > b.pt; });
	    }	      

	  float max_truth = 0;
	  float max_reco = 0;
	  float min_truth = 0;
	  float min_reco = 0;

	  /* this variable determines how things are filled */

	  // 0 - truth and reco good and matched
	  // 1 - truth good and no reco
	  // 2 - truth and reco good and no match (fake and miss)
	  // 3 - truth not good and reco good and no match (fake)

	  int fill_fake_miss = 0;

	  bool reco_good = djf.check_dijet_reco(myrecojets2);
	  bool base_good = reco_good;
	  
	  bool njet_pass = (myrecojets.size() < njet_cut);

	  reco_good &= triggered;
	  reco_good &= njet_pass;

	  bool in_range = false;
	  
	  if (myrecojets2.size() > 0)
	    {
	      in_range = (myrecojets2.at(0).pt < ipt_bins[max_reco_bin]);
	      reco_good &= in_range;
	    }
	  
	  
	  //if (!reco_good) continue;
	  
	  if (truth_good)
	    {
	      max_truth = mytruthjets2.at(0).pt;
	      min_truth = mytruthjets2.at(1).pt;
	    }
	  
	  if (reco_good)
	    {
	      max_reco = myrecojets2.at(0).pt;
	      min_reco = myrecojets2.at(1).pt;
	    }
	 
	  if (matched && truth_good && reco_good)
	    {
	      max_truth = matched_dijets.at(0).first.pt;
	      min_truth = matched_dijets.at(1).first.pt;
	      max_reco = matched_dijets.at(0).second.pt;
	      min_reco = matched_dijets.at(1).second.pt;
	    }

	  if (primer != 1 && reco_good)
	    {
	      int eta1 = floor((myrecojets2.at(0).eta + 1.1)/0.1) + 1;
	      int eta2 = floor((myrecojets2.at(1).eta + 1.1)/0.1) + 1;
	      int gbinn = h_eta_reweight->GetBin(eta1, eta2);
	      event_scale *= h_eta_reweight->GetBinContent(gbinn);
	    }

	  double fake_event_scale = event_scale;
	  if (base_good)
	    {
	      he_dijet_reco->Fill(truth_good, myrecojets2.at(0).pt, myrecojets2.at(1).pt, event_scale);
	      he_dijet_reco->Fill(truth_good, myrecojets2.at(1).pt, myrecojets2.at(0).pt, event_scale);
	    }
	  if (truth_good)
	    {
	      if (matched)
		{
		  he_dijet_match_reco->Fill(base_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
		  he_dijet_match_reco->Fill(base_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
		}
	      if(base_good)
		{
	      
		  he_dijet_match->Fill(matched, myrecojets2.at(0).pt, myrecojets2.at(1).pt, event_scale);
		  he_dijet_match->Fill(matched, myrecojets2.at(1).pt, myrecojets2.at(0).pt, event_scale);
	
		}

	      he_dijet_vertex->Fill(has_vertex, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_vertex->Fill(has_vertex, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	  
	      he_dijet_mbd_hit->Fill(has_mbd_hit, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_mbd_hit->Fill(has_mbd_hit, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);

	      he_dijet_trigger->Fill(triggered, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      he_dijet_trigger->Fill(triggered, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);	  
	  
	      he_dijet->Fill(matched && reco_good, mytruthjets2.at(1).pt, mytruthjets2.at(0).pt, event_scale);
	      he_dijet->Fill(matched && reco_good, mytruthjets2.at(0).pt, mytruthjets2.at(1).pt, event_scale);
	      
	    }
	  
	  float pt1_truth_bin = nbins_pt;
	  float pt2_truth_bin = nbins_pt;
	  float pt1_reco_bin = nbins_pt;
	  float pt2_reco_bin = nbins_pt;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;


	  	  
	  float maxi = std::max(es1, es2);
	  float mini = std::min(es1, es2);

	  float maxit = std::max(e1, e2);
	  float minit = std::min(e1, e2);

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
	    }

	  if (pt2_truth_bin == nbins_pt)
	    {
	      truth_good = false;
	    }
	  if (primer != 1)
	    {
	      
	      if (!h_flat_truth_mapping_primer->GetBinContent(1 + pt1_truth_bin*nbins_pt + pt2_truth_bin))
		{
		  truth_good = false;
		}

	      if (!h_flat_reco_mapping_primer->GetBinContent(1 + pt1_reco_bin*nbins_pt + pt2_reco_bin))
		{
		  reco_good = false;
		}

	      if (truth_good && reco_good)
		{
		  int binn = h_count_flat_response_pt1pt2_sample_primer[isample]->GetBin(1 + pt1_reco_bin*nbins_pt + pt2_reco_bin, 1 + pt1_truth_bin*nbins_pt + pt2_truth_bin);
		  if (h_count_flat_response_pt1pt2_sample_primer[isample]->GetBinContent(binn) < 5) continue;
		}
	      if (reco_good && h_count_flat_reco_pt1pt2_sample_primer[isample]->GetBinContent(1 + pt1_reco_bin*nbins_pt + pt2_reco_bin) < 5) continue;
	      if (truth_good && h_count_flat_truth_pt1pt2_sample_primer[isample]->GetBinContent(1 + pt1_truth_bin*nbins_pt + pt2_truth_bin) < 5) continue;;
	      

	    }
	  if (!prior_sys && !primer && truth_good)// && reco_good && matched)
	    {
	      int recorrectbin = pt1_truth_bin*nbins_pt + pt2_truth_bin;	      
	      event_scale *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
	    }

	  if (reco_good && truth_good && matched)
	    {
	      fill_fake_miss = 0;
	      he_dijet_fake_binned->Fill(1, pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
	      he_dijet_fake_binned->Fill(1, pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);
	      
	      he_dijet_miss_binned->Fill(1, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(1, e2, e1, event_scale);
	    }

	  if (reco_good && truth_good && !matched)
	    {
	      h_reco_fakes->Fill(es1, es2);
	      h_reco_fakes->Fill(es2, es1);
	      h_match_fakes->Fill(es1, es2);
	      h_match_fakes->Fill(es2, es1);

	      fill_fake_miss = 2;	      
	      he_dijet_fake_binned->Fill(0, pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
	      he_dijet_fake_binned->Fill(0, pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);
	      he_dijet_miss_binned->Fill(0, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(0, e2, e1, event_scale);

	    }
	  if (reco_good && !truth_good)
	    {
	      fill_fake_miss = 3;
	      h_reco_fakes->Fill(es1, es2);
	      h_reco_fakes->Fill(es2, es1);
	      h_no_truth_fakes->Fill(es1, es2);
	      h_no_truth_fakes->Fill(es2, es1);
	      he_dijet_fake_binned->Fill(0, pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
	      he_dijet_fake_binned->Fill(0, pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);
	    }
	  if (truth_good && !reco_good)
	    {
	      fill_fake_miss = 1;
	      he_dijet_miss_binned->Fill(0, e1, e2, event_scale);
	      he_dijet_miss_binned->Fill(0, e2, e1, event_scale);
	    }

	  h_truth_lead->Fill(e1, event_scale);
	  h_truth_lead_sample[isample]->Fill(maxpttruth, event_scale);
	  
	  if (minit >  truth_subleading_cut) h_truth_sublead->Fill(e2, event_scale);

	  if (maxi >  reco_leading_cut) h_reco_lead->Fill(maxi, event_scale);
	  if (mini >  reco_subleading_cut) h_reco_sublead->Fill(mini, event_scale);

	  
	  // miss
	  if (fill_fake_miss == 1)
	    {
	      
	      if (fill_response)
		{
		  
		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);

		  h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);

		  rooResponse.Miss(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  rooResponse.Miss(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

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

		  h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

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
		  
		  h_flat_reco_to_response_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_to_response_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);
		  h_flat_response_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_response_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_response_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_response_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, pt2_truth_bin + nbins_pt*pt1_truth_bin);
		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);
		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);	      

		  h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  h_flat_reco_all_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_all_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);
		  h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, pt2_truth_bin + nbins_pt*pt1_truth_bin);
		  h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);
		  h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);	      

		  
		  rooResponse.Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin,pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  rooResponse.Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin,pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  
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

		  h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);
		  
		  h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);
		}		
		
		
	      if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
		{
		  h_reco_xj->Fill(mini/maxi, event_scale);
		  h_linear_reco_xj->Fill(mini/maxi, event_scale);
		}

	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(myrecojets2.at(0).eta, myrecojets2.at(1).eta, event_scale);
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
		  h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  rooResponse.Fake(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  rooResponse.Fake(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  rooResponse.Miss(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  rooResponse.Miss(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);
		  
		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);

		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2_sample[isample]->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);

		  h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);

		  h_flat_reco_all_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_all_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		  h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin);
		  h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin);
		  
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

		  // if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_bins[3]])
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

		  h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, event_scale);
		  h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, event_scale);

		  h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		}
	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(myrecojets2.at(0).eta, myrecojets2.at(1).eta, event_scale);
		  h_mbd_vertex->Fill(mbd_vertex_z[isample], event_scale);
		  h_njets->Fill(nnrecojets, event_scale);	      
		}
	      continue;
	    }
	  // only reco fake
	  if (fill_fake_miss == 3)
	    {
	      
	      if (fill_response)
		{

		  h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  rooResponse.Fake(pt1_reco_bin + nbins_pt*pt2_reco_bin, event_scale);
		  rooResponse.Fake(pt2_reco_bin + nbins_pt*pt1_reco_bin, event_scale);

		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);	      

		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_all_pt1pt2_sample[isample]->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		  h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin);
		  h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin);	      

		  h_flat_reco_all_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_all_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		}
	      if (fill_unfold)
		{

		  h_pt1pt2->Fill(es1, es2, event_scale);
		  h_pt1pt2->Fill(es2, es1, event_scale);
		  h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins_pt*pt2_reco_bin, fake_event_scale);
		  h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins_pt*pt1_reco_bin, fake_event_scale);

		}


	      if (reco_good)
		{
		  h_eta_lead_sublead->Fill(myrecojets2.at(0).eta, myrecojets2.at(1).eta, event_scale);
		  h_mbd_vertex->Fill(mbd_vertex_z[isample], event_scale);
		  h_njets->Fill(nnrecojets, event_scale);	      
		}
	      continue;
	    } 

	}
    }


  TH1D *h_dijet_fake_rate = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
  h_dijet_fake_rate->SetName("h_dijet_fake_rate");
  h_dijet_fake_rate->Divide(h_flat_reco_all_pt1pt2);


  TH2D *h_fake_v_entries = new TH2D("h_fake_v_entries",";Fake Rate; entries", 20, 0, 1, 24000, 0, 24000);
  
  TH2D *h_pt1pt2_reco_before = new TH2D("h_pt1pt2_reco_before", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco_before, nbins_pt);

  TH1D *h_flat_reco_before_pt1pt2 = (TH1D*) h_flat_reco_pt1pt2->Clone();
  h_flat_reco_before_pt1pt2->SetName("h_flat_reco_before_pt1pt2");
    
  TH2D *h_pt1pt2_reco_response_before = new TH2D("h_pt1pt2_reco_response_before", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_to_response_pt1pt2, h_pt1pt2_reco_response_before, nbins_pt);

  for (int ib = 0; ib < nbins_pt; ib++)
    {
      for (int ib2 = 0; ib2 < nbins_pt; ib2++)
	{
	  int ibb = 1 + ib*nbins_pt + ib2;
	  double value = he_dijet_fake_binned->GetEfficiency(ibb);

	  h_fake_v_entries->Fill(value, h_count_flat_reco_pt1pt2->GetBinContent(ibb));
	  //h_flat_reco_pt1pt2->SetBinContent(ibb, h_flat_reco_pt1pt2->GetBinContent(ibb)*value);
	}
    }

  h_flat_reco_pt1pt2->Scale(.5);
  h_flat_truth_pt1pt2->Scale(.5);

  // TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  // TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  // TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  // TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  // TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt, nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  // TH1D h_binnumbers_truth_sample[5];
  
  
  // for (int i = 0; i < 5; i++)
  //   {
  //     h_binnumbers_truth_sample[i] = new TH1D(Form("h_binnumbers_truth_sample_%d", i), "", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
      
  //     h_flat_reco_to_response_pt1pt2_sample[isample]->Scale(.5);
  //     h_flat_reco_all_pt1pt2_sample[isample]->Scale(.5);
  //     h_flat_truth_to_response_pt1pt2_sample[isample]->Scale(.5);
  //     h_flat_response_pt1pt2_sample[isample]->Scale(.5);

  //     for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
  // 	{
  // 	  int count_truth = h_count_flat_truth_pt1pt2_sample[i]->GetBinContent(ib+1);
  // 	  if (count_truth > 10)
  // 	    {
	      
  // 	    }
  // 	}
  //   }
  
  int number_of_mins = 1;
  if (minentries == 0) number_of_mins = 20;
  for (int imin = minentries; imin < minentries + number_of_mins; imin++)       
    {
      std::cout << "Minimum Entries == " << imin << std::endl;
      //int nbinsx = h_flat_response_pt1pt2->GetXaxis()->GetNbins();
      //int nbinsy = h_flat_response_pt1pt2->GetYaxis()->GetNbins();

      TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
      TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
      if (minentries)
	{
	  h_flat_truth_mapping->SetName("h_flat_truth_mapping");
	  h_flat_reco_mapping->SetName("h_flat_reco_mapping");
	}
      else
	{
	  h_flat_truth_mapping->SetName(Form("h_flat_truth_mapping_min%d", imin));
	  h_flat_reco_mapping->SetName(Form("h_flat_reco_mapping_min%d", imin));
	}
      
      h_flat_reco_mapping->Reset();
      h_flat_truth_mapping->Reset();
  
      int nempty_reco = 0;
      std::vector<int> binnumbers_reco{};
      int nempty_truth = 0;
      std::vector<int> binnumbers_truth{};
      int nrecobins = 0;
      int ntruthbins = 0;
      for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
	{
	  bool use_reco = true;
	  bool use_truth = true;
	  if (primer ==1)
	    {
	      float content_reco = h_flat_reco_pt1pt2->GetBinContent(ib+1);
	      float error_reco = h_flat_reco_pt1pt2->GetBinError(ib+1);
	      float content_truth = h_flat_truth_pt1pt2->GetBinContent(ib+1);
	      float error_truth = h_flat_truth_pt1pt2->GetBinError(ib+1);
	      float significance_reco = (content_reco > 0 ? content_reco/error_reco : 0);
	      float significance_truth = (content_truth > 0 ? content_truth/error_truth : 0);
	      use_reco = (content_reco < minentries || he_dijet_fake_binned->GetEfficiency(ib + 1) < 0.3);
	      use_truth = (content_truth < minentries);		
	    }
	  else
	    {
	      use_reco = (h_flat_reco_mapping_primer->GetBinContent(ib+1) == 0);
	      use_truth = (h_flat_truth_mapping_primer->GetBinContent(ib+1) == 0);
	    }
	  
	  if (use_reco)
	    {
	      nempty_reco++;
	      binnumbers_reco.push_back(0);
	      h_flat_reco_mapping->SetBinContent(ib+1, 0);
	    }
	  else
	    {
	      nrecobins++;
	      binnumbers_reco.push_back(nrecobins);
	      h_flat_reco_mapping->SetBinContent(ib+1, nrecobins);
	    }
	  if (use_truth)
	    {
	      nempty_truth++;
	      binnumbers_truth.push_back(0);
	      h_flat_truth_mapping->SetBinContent(ib+1, 0);
	    }
	  else
	    {
	      ntruthbins++;
	      binnumbers_truth.push_back(ntruthbins);
	      h_flat_truth_mapping->SetBinContent(ib+1, ntruthbins);
	    }
	}
      // for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
      // 	{
      // 	  if (h_count_flat_reco_pt1pt2->GetBinContent(ib+1) < imin)
      // 	    {
      // 	      nempty_reco++;
      // 	      binnumbers_reco.push_back(0);
      // 	      h_flat_reco_mapping->SetBinContent(ib+1, 0);
      // 	    }
      // 	  else
      // 	    {
      // 	      nrecobins++;
      // 	      binnumbers_reco.push_back(nrecobins);
      // 	      h_flat_reco_mapping->SetBinContent(ib+1, nrecobins);
      // 	    }
      // 	  if (h_count_flat_truth_pt1pt2->GetBinContent(ib+1) < imin)
      // 	    {
      // 	      nempty_truth++;
      // 	      binnumbers_truth.push_back(0);
      // 	      h_flat_truth_mapping->SetBinContent(ib+1, 0);
      // 	    }
      // 	  else
      // 	    {
      // 	      ntruthbins++;
      // 	      binnumbers_truth.push_back(ntruthbins);
      // 	      h_flat_truth_mapping->SetBinContent(ib+1, ntruthbins);
      // 	    }
      // 	}

      TH1D *h_flat_truth_skim = new TH1D(Form("h_flat_truth_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
      TH1D *h_flat_reco_skim = new TH1D(Form("h_flat_reco_skim_min%d", imin),"", nrecobins, 0, nrecobins);
      TH1D *h_flat_truth_to_unfold_skim = new TH1D(Form("h_flat_truth_to_unfold_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
      TH1D *h_flat_reco_to_unfold_skim = new TH1D(Form("h_flat_reco_to_unfold_skim_min%d", imin),"", nrecobins, 0, nrecobins);

      TH2D *h_flat_response_skim = new TH2D(Form("h_flat_response_skim_min%d", imin),"", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

      TH1D *h_flat_truth_to_unfold = (TH1D*) h_flat_truth_pt1pt2->Clone();
      h_flat_truth_to_unfold->SetName("h_flat_truth_to_unfold_pt1pt2");
      h_flat_truth_to_unfold->Reset();
      
      if (minentries)
	{
	  h_flat_response_skim->SetName("h_flat_response_skim");
	  h_flat_truth_skim->SetName("h_flat_truth_skim");
	  h_flat_reco_skim->SetName("h_flat_reco_skim");
	  h_flat_truth_to_unfold_skim->SetName("h_flat_truth_to_unfold_skim");
	  h_flat_reco_to_unfold_skim->SetName("h_flat_reco_to_unfold_skim");
	}

      for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
	{
	  int bintruth = binnumbers_truth.at(ib);      
	  int binreco = binnumbers_reco.at(ib);
	  if (binreco)
	    {
	      h_flat_reco_skim->SetBinContent(binreco, h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1));
	      h_flat_reco_skim->SetBinError(binreco, h_flat_reco_to_response_pt1pt2->GetBinError(ib+1));
	      h_flat_reco_to_unfold_skim->SetBinContent(binreco, h_flat_reco_pt1pt2->GetBinContent(ib+1));
	      h_flat_reco_to_unfold_skim->SetBinError(binreco, h_flat_reco_pt1pt2->GetBinError(ib+1));

	    }

	  if (bintruth)
	    {
	 
	      h_flat_truth_skim->SetBinContent(bintruth, h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1));
	      h_flat_truth_skim->SetBinError(bintruth, h_flat_truth_to_response_pt1pt2->GetBinError(ib+1));
	      h_flat_truth_to_unfold_skim->SetBinContent(bintruth, h_flat_truth_pt1pt2->GetBinContent(ib+1));
	      h_flat_truth_to_unfold_skim->SetBinError(bintruth, h_flat_truth_pt1pt2->GetBinError(ib+1));

	      for (int ibr = 0; ibr < nbins_pt*nbins_pt; ibr++)
		{
		  int binreco2 = binnumbers_reco.at(ibr);
		  if (binreco2)
		    {
		      int rbin = h_flat_response_skim->GetBin(binreco2, bintruth);
		      int tbin = h_flat_response_pt1pt2->GetBin(ibr+1, ib+1);
		      h_flat_response_skim->SetBinContent(rbin, h_flat_response_pt1pt2->GetBinContent(tbin));
		      h_flat_response_skim->SetBinError(rbin, h_flat_response_pt1pt2->GetBinError(tbin));
		    }
		}
	    }
	}

      std::cout <<" Nbins skim reco = "<<h_flat_reco_skim->GetNbinsX()<<std::endl;
      std::cout <<" Nbins skim truth = "<<h_flat_truth_skim->GetNbinsX()<<std::endl;
      std::cout <<" Nbins skim response = "<<h_flat_response_skim->GetXaxis()->GetNbins()<<"  --  " << h_flat_response_skim->GetYaxis()->GetNbins()<<std::endl; 

      RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);

      if (minentries)
	{
	  rooResponsehist.SetName("response_noempty");
	}
      else
	{
	  rooResponsehist.SetName(Form("response_noempty_min%d", imin));
	}

      if (minentries)
	{
	  TH1D* h_flat_unfold_skim[niterations];
	  TH1D* h_flat_unfold_pt1pt2[niterations];
	  TH1D* h_flat_refold_skim[niterations];
	  TH1D* h_flat_refold_pt1pt2[niterations];
	  TH1D* h_flat_truth_fold_skim;
	  TH1D* h_flat_truth_fold_pt1pt2;
	  int niter = 3;

	  for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
	    {
	      int bin = binnumbers_truth.at(ib);
	      if (bin)
		{
		  h_flat_truth_to_unfold->SetBinContent(ib+1, h_flat_truth_to_unfold_skim->GetBinContent(bin));
		  h_flat_truth_to_unfold->SetBinError(ib+1, h_flat_truth_to_unfold_skim->GetBinError(bin));
		}
	    }

	  h_flat_truth_fold_skim = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_truth_to_unfold_skim, "hTruth_Folded");

	  h_flat_truth_fold_pt1pt2 = (TH1D*) h_flat_truth_pt1pt2->Clone();
	  h_flat_truth_fold_pt1pt2->Reset();
	  h_flat_truth_fold_pt1pt2->SetName("h_flat_truth_fold_pt1pt2");

	  for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
	    {
	      int bin = binnumbers_truth.at(ib);
	      if (bin)
		{
		  h_flat_truth_fold_pt1pt2->SetBinContent(ib+1, h_flat_truth_fold_skim->GetBinContent(bin));
		  h_flat_truth_fold_pt1pt2->SetBinError(ib+1, h_flat_truth_fold_skim->GetBinError(bin));
		}
	    }

	  
	  for (int iter = 0; iter < niterations; iter++ )
	    {
      
	      //RooUnfoldBayes   unfold (&rooResponse, h_flat_reco_pt1pt2, iter + 1);    // OR
	      RooUnfoldBayes   unfold (&rooResponsehist, h_flat_reco_to_unfold_skim, iter + 1);    // OR

	      //TH1D *hc = (TH1D*) unfold.Hunfold();
	      //h_flat_unfold_pt1pt2[iter] = (TH1D*) hc->Clone();//h_flat_truth_pt1pt2->Clone();
	      //h_flat_unfold_pt1pt2[iter]->Reset();
	      //h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
	      


	       //std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;

	      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
	      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
	      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
	      h_flat_unfold_pt1pt2[iter]->Reset();
	      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
	      for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
		{
		  int bin = binnumbers_truth.at(ib);
		  if (bin)
		    {
		      h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin));
		      h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin));
		    }
		}

	      h_flat_refold_skim[iter] = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_unfold_skim[iter], "hRefolded");

	      std::cout <<" Nbins skim reco = "<<h_flat_refold_skim[iter]->GetNbinsX()<<std::endl;
	      h_flat_refold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
	      h_flat_refold_pt1pt2[iter]->Reset();
	      h_flat_refold_pt1pt2[iter]->SetName(Form("h_flat_refold_pt1pt2_%d",iter));
	      for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
		{
		  int bin = binnumbers_truth.at(ib);
		  if (bin)
		    {
		      h_flat_refold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_refold_skim[iter]->GetBinContent(bin));
		      h_flat_refold_pt1pt2[iter]->SetBinError(ib+1, h_flat_refold_skim[iter]->GetBinError(bin));
		    }
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

	  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 1000, 300);
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

	  h_truth_lead->SetTitle(";Leading Jet p_{T} [GeV];counts * lumiscale");
	  h_truth_lead->SetMinimum(1);
	  h_truth_lead->Draw("p");
	  h_truth_lead_sample[0]->Draw("same p");
	  h_truth_lead_sample[1]->Draw("same p");
	  h_truth_lead_sample[2]->Draw("same p");
	  h_truth_lead_sample[3]->Draw("same p");
	  h_truth_lead_sample[4]->Draw("same p");
	  h_truth_lead->Draw("p same");


	  TLine  *line1 = new TLine(14, 1, 14, h_truth_lead->GetBinContent(15));
	  dlutility::SetLineAtt(line1, kRed, 1.5, 1);
	  line1->Draw("same");
	  TLine  *lin2 = new TLine(20, 1, 20, h_truth_lead->GetBinContent(21));
	  dlutility::SetLineAtt(lin2, kGreen, 1.5, 1);
	  lin2->Draw("same");
	  TLine  *lin3 = new TLine(30, 1, 30, h_truth_lead->GetBinContent(31));
	  dlutility::SetLineAtt(lin3, kBlue, 1.5, 1);
	  lin3->Draw("same");
	  cjetdiv->cd(2);
	  if (!use_herwig)
	    {
	      dlutility::DrawSPHENIXppsize(0.05, 0.84, 0.08, 1, 0, 1, "Pythia8");
	    }
	  else if (use_herwig == 1)
	    {
	      dlutility::DrawSPHENIXppsize(0.05, 0.84, 0.08, 1, 0, 1, "HERWIG 7.8");
	    }
	  dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.05, 0.74, 0, kBlack, 0.08);
	  dlutility::drawText("All Dijet Pairs", 0.05, 0.69, 0, kBlack, 0.08);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.05, 0.64, 0, kBlack, 0.08);

	  TLegend *leg3 = new TLegend(0.01, 0.13, 0.7, 0.62);
	  leg3->SetLineWidth(0);
	  leg3->AddEntry(h_truth_lead, "Combined Sample");
	  leg3->AddEntry(h_truth_lead_sample[0], "10 GeV Sample");
	  leg3->AddEntry(h_truth_lead_sample[1], "20 GeV Sample");
	  leg3->AddEntry(h_truth_lead_sample[2], "30 GeV Sample");
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
	  rooResponsehist.Write();
	  rooResponse.Write();
	  h_reco_fakes->Write();
	  h_match_fakes->Write();
	  h_no_truth_fakes->Write();

	  h_count_flat_truth_pt1pt2->Write();
	  h_count_flat_reco_pt1pt2->Write();

	  for (int i = 0; i < 5; i++)
	    {
	      h_count_flat_truth_pt1pt2_sample[i]->Write();
	      h_count_flat_reco_pt1pt2_sample[i]->Write();
	      h_count_flat_response_pt1pt2_sample[i]->Write();
	    }
	  h_fake_v_entries->Write();
	  h_eta_lead_sublead->Write();
	  h_flat_reco_pt1pt2->Write();
	  h_flat_reco_before_pt1pt2->Write();
	  h_flat_reco_to_response_pt1pt2->Write();
	  h_flat_reco_all_pt1pt2->Write();
	  h_flat_truth_pt1pt2->Write();
	  h_flat_truth_to_response_pt1pt2->Write();
	  h_flat_truth_fold_pt1pt2->Write();
	  h_flat_response_pt1pt2->Write();
	  h_flat_truth_mapping->Write();
	  h_flat_reco_mapping->Write();
	  h_flat_truth_skim->Write();
	  h_flat_truth_to_unfold_skim->Write();
	  h_flat_truth_to_unfold->Write();
	  h_flat_reco_skim->Write();
	  h_flat_response_skim->Write();
	  h_mbd_vertex->Write();
	  h_njets->Write();
	  h_pt1pt2->Write();
	  h_pt1pt2_reco_before->Write();
	  h_pt1pt2_reco_response_before->Write();
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
	  h_xj_reco->Write();
	      
	  fr->Close();

	}
      else
	{
	  TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_min_" + std::to_string(imin) + ".root";

	  TFile *fr = new TFile(responsepath.Data(),"recreate");
	  rooResponsehist.Write();
	  
	  TH1D *h1 = (TH1D*)h_flat_reco_pt1pt2->Clone();
	  TH1D *h2 = (TH1D*)h_flat_truth_pt1pt2->Clone();
	  TH2D *h3 = (TH2D*)h_flat_response_pt1pt2->Clone();
	  TH1D *h4 = (TH1D*)h_flat_truth_mapping->Clone();
	  TH1D *h5 = (TH1D*)h_flat_reco_mapping->Clone();
	  TH1D *h6 = (TH1D*)h_flat_truth_skim->Clone();
	  TH1D *h7 = (TH1D*)h_flat_reco_skim->Clone();
	  TH2D *h8 = (TH2D*)h_flat_response_skim->Clone();
	  h1->Write();
	  h2->Write();
	  h3->Write();
	  h4->Write();
	  h5->Write();
	  h6->Write();
	  h7->Write();
	  h8->Write();
			  	  
	  fr->Close();
	}
      std::cout << "Reco empty: " << nempty_reco << std::endl;
      std::cout << "Truth empty: " << nempty_truth << std::endl;
    }
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
      std::cout << "    ./createResponse_noempty_pp -c binning.config -r 4 -n 10 -p 1 -h 0 -f 0 -e 0 -ig 0 -ug 0" << std::endl;
      return 1;
    }
  
  return createResponse_noempty_pp(config, full_or_half, niterations, cone_size, primer, useFakes, useEfficiencies, input_generator, unfold_generator);
  
}
