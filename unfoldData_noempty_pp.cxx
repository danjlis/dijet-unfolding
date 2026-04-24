#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "dijetfinder.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

static int verbosity = 0;

int unfoldData_noempty_pp(const std::string configfile = "binning.config", const int niterations = 10, const int cone_size = 4, const int primer = 0, const int unfold_generator = 0, const int input_generator = 0, const int full_or_half = 0)
{
  gStyle->SetPalette(kRainBow);
  gRandom->SetSeed(0);
  if (unfold_generator > 2) return 1;

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

  int version_arr[3] = {8, 11, 10};
  std::string generator_name[3] = {"NA", "PYTHIA", "HERWIG"};

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  
  bool ispp = true;
  std::string system_string = "pp";
    
  read_binning rb(configfile.c_str());
  std::vector<int> runnumbers = rb.getRunnumbers();

  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TFile *finjer = new TFile(Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)), "r");


  TFile *in_trigger = new TFile(Form("%s/efficiencies/trigger_efficiencies.root", rb.get_code_location().c_str()), "r");
  TF2 *ftrigger = (TF2*) in_trigger->Get(Form("f2_jet10_eff_r0%d", cone_size));

  // TFile *in_purity = new TFile(Form("%s/efficiencies/purity_correction_r%02d.root", rb.get_code_location().c_str(), cone_size), "r");
  // TH1D *h_purity = (TH1D*) in_purity->Get(Form("h_purity_r%02d", cone_size));

  dijetfinder djf(cone_size);
  djf.SetVerbosity(verbosity);

  std::string data_file = rb.get_tntuple_location() + "/TREE_DIJET_SKIM_r0" + std::to_string(cone_size) + "_v8_7_ana533_2025p007_v001_gl10-alltime.root";

  int sim_version = version_arr[unfold_generator];
  
  std::string j_file[5];
  
  if (unfold_generator == 0)
    {
      j_file[0] = data_file;
    }
  else
    {
      j_file[0] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_12_ana509_MDC2-00000028.root";
      j_file[1] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
      j_file[2] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
      j_file[3] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_40_ana509_MDC2-00000028.root";
      j_file[4] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_50_ana509_MDC2-00000028.root";
    }

  int sample_skip[5] = {0};
  
  if (unfold_generator == 1)
    {

      for (int i = 0; i < 5; i++)
	{
	  sample_skip[i] = 1;
	}

    }
  else if (unfold_generator == 2)
    {
     
      sample_skip[0] = 1;
      sample_skip[3] = 1;
    }
  else if (unfold_generator == 0)
    {
      sample_skip[0] = 1;
    }
  
  ULong64_t gl1_scaled[5];;
  float mbd_vertex_z[5];
  int mbd_hit[5];  
  int runnumber[5];
  double mbd_avgsigma[5];
  double calib_lead_time;
  double calib_delta_time;
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

  float n_events[5];
  TFile *fin[5];
  TTree *ttree[5];

  
  for (int i = 0; i < 5; i ++)
    {
      if (!sample_skip[i]) continue;
      fin[i] = new TFile(j_file[i].c_str(), "r");
      ttree[i]  = (TTree*) fin[i]->Get("ttree");;

      if (!ttree[i])
	{
	  std::cout << " no data "<< std::endl;
	  return 1;
	}

      if (unfold_generator > 0)
	{
	  if (cone_size != 4)
	    {
	      ttree[i]->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_ref[i]);
	    }
	  ttree[i]->SetBranchAddress(Form("truth_jet_pt_%d", cone_size), &truth_jet_pt[i]);
	  ttree[i]->SetBranchAddress(Form("truth_jet_eta_%d", cone_size), &truth_jet_eta[i]);
	  ttree[i]->SetBranchAddress(Form("truth_jet_phi_%d", cone_size), &truth_jet_phi[i]);
	}
      else
	{
	  ttree[i]->SetBranchAddress("runnumber", &runnumber[i]);
	  ttree[i]->SetBranchAddress("mbd_avgsigma", &mbd_avgsigma[i]);
	}
      ttree[i]->SetBranchAddress(Form("jet_pt_calib_%d", cone_size), &reco_jet_pt[i]);
      ttree[i]->SetBranchAddress(Form("jet_pt_%d", cone_size), &reco_jet_pt_uncalib[i]);
      ttree[i]->SetBranchAddress(Form("jet_emcal_%d", cone_size), &reco_jet_emcal[i]);
      ttree[i]->SetBranchAddress(Form("jet_e_%d", cone_size), &reco_jet_e[i]);
      ttree[i]->SetBranchAddress(Form("jet_eta_%d", cone_size), &reco_jet_eta[i]);
      ttree[i]->SetBranchAddress(Form("jet_eta_det_%d", cone_size), &reco_jet_eta_det[i]);
      ttree[i]->SetBranchAddress(Form("jet_phi_%d", cone_size), &reco_jet_phi[i]);
      ttree[i]->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z[i]);
      ttree[i]->SetBranchAddress("calib_lead_time", &calib_lead_time);
      ttree[i]->SetBranchAddress("calib_delta_time", &calib_delta_time);

      ttree[i]->SetBranchAddress("mbd_hit", &mbd_hit[i]);
      ttree[i]->SetBranchAddress("gl1_scaled", &gl1_scaled[i]);
      n_events[i] = ttree[i]->GetEntries();
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

  if (unfold_generator == 0)
    {
      scale_factor[0] = 1;
    }

  if (unfold_generator == 2)
    {
      scale_factor[0] = 1.4216108 * ((float) n_events[3]/ ((float)n_events[0])) * csh_10/(csh_30 * 4.);
      scale_factor[3] = 1;
    }

  
  Int_t read_nbins = rb.get_nbins();

  //Double_t dphicut = rb.get_dphicut();
  Double_t vtx_cut = rb.get_vtx_cut();
  Double_t njet_cut = rb.get_njet_cut();
  Int_t trigger_sys = rb.get_trigger_sys();
  Int_t philoc_sys = rb.get_philoc_sys();
  Int_t ca_sys = rb.get_crossingangle_sys();
  Int_t full_sys = rb.get_full_sys();
    
  const int nbins = read_nbins;
  const int nbins_pt = read_nbins + 1;
  
  //Int_t njet_sys = rb.get_njet_sys();
  Int_t prior_sys = rb.get_prior_sys();
  Int_t emfrac_sys = rb.get_emfrac_sys();
  Double_t pileup_sys = rb.get_pileup_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  Int_t herwig_sys = rb.get_herwig();
  
  std::string sys_name = "nominal"; 
  if (pileup_sys > 1)
    {
      sys_name = "PILEUP";
    } 
  else  if (pileup_sys > 0)
    {
      sys_name = "PILEUPMIX";
    } 
  if (full_sys)
    sys_name = "FULL";
  if (prior_sys)
    sys_name = "PRIOR";
  if (emfrac_sys)
    sys_name = "EMFRAC";

  if (trigger_sys)
    {
      sys_name = "TRIGGER";
    }
  int runnumber_low = 0;
  int runnumber_high = 100000;
  if (ca_sys)
    {      
      sys_name = "CA" + std::to_string(ca_sys) ;
      if (ca_sys == 1)
	{
	  runnumber_high = 52212;
	}
      else
	{
	  runnumber_low = 52212;
	}
    }

  if (philoc_sys)
    {
      sys_name = "PHILOC" + std::to_string(philoc_sys);;
    }

  if (herwig_sys)
    sys_name = "HERWIG";
  
  if (JER_sys < 0)
    sys_name = "negJER";

  if (JER_sys > 0)
    sys_name = "posJER";

  if (JES_sys < 0)
    sys_name = "negJES";

  if (JES_sys > 0)
    sys_name = "posJES";



  TString responsetestpath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);

 
  
  TFile *fresponsetest;
  
  TH1D *h_flat_data_truth_pt1pt2 = nullptr;

  
  if (unfold_generator == 2 && input_generator == 1)
    {
      sys_name = "PuH";
      responsetestpath += "_PRIMER1_HALF_HERWIG.root";
    }
  if (unfold_generator == 1 && input_generator == 1)
    {

      if (full_or_half)
	{
	  sys_name = "HALF_PuP";
	  responsetestpath += "_PRIMER1_HALF_PuP.root";
	}
      else
	{
	  sys_name = "PuP";
	  responsetestpath += "_PRIMER1_nominal.root";
	}

    }
  if (unfold_generator == 2 && input_generator == 2)
    {
      if (full_or_half)
	{
	  sys_name = "HALF_HuH";
	  responsetestpath += "_PRIMER1_HALF_HuH.root";
	}
      else
	{
	  sys_name = "HuH";
	  responsetestpath += "_PRIMER1_HERWIG.root";
	}
    }
  if (unfold_generator == 1 && input_generator == 2)
    {
      sys_name = "HuP";
      responsetestpath += "_PRIMER1_nominal.root";
    }

  if (unfold_generator && input_generator)
    {
      fresponsetest = new TFile(responsetestpath.Data(),"r");

      h_flat_data_truth_pt1pt2 = (TH1D*) fresponsetest->Get("h_truth_flat_pt1pt2"); 

      if (!h_flat_data_truth_pt1pt2)
	{
	  std::cout << "no truth for data" << std::endl;
	  return 1;
	}
    }
  
  if (h_flat_data_truth_pt1pt2)
    {
      h_flat_data_truth_pt1pt2->SetName("h_flat_data_truth_pt1pt2");
    }

  std::string sys_name_orig = sys_name;
  
  if (primer == 1)
    {
      sys_name = "PRIMER1_" + sys_name;
    }
  else if (primer == 2)
    {
      sys_name = "PRIMER2_" + sys_name;
    }

  TF1 *gjer = (TF1*) finjer->Get("jernom");
  TH1D *hjersmear = (TH1D*) finjer->Get("h_nominal_smear");
  if (JER_sys != 0)
    {
      if (JER_sys < 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_down_smear");
	  gjer = (TF1*) finjer->Get("jerneg");
	}
      else if (JER_sys > 0)
	{
	  hjersmear = (TH1D*) finjer->Get("h_sys_up_smear");
	  gjer = (TF1*) finjer->Get("jerpos");
	}
    }

  
  float ipt_bins[nbins_pt+1];
  double dpt_bins[nbins_pt+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      dpt_bins[i] = (double) ipt_bins[i];
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }
  ipt_bins[nbins_pt] = 100;
  dpt_bins[nbins_pt] = 100;
  
  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
  
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  djf.setRecoCuts(reco_leading_cut, reco_subleading_cut);
  djf.setTruthCuts(truth_leading_cut, truth_subleading_cut);
  djf.setPhiLocation(philoc_sys);
  

  int exbin = 1;
  int measure_bins[10] = {0};
  int mbins = 3;  
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }


  Int_t max_reco_bin = rb.get_maximum_reco_bin();

  std::cout << "Max reco: " << max_reco_bin << " - ( " << ipt_bins[max_reco_bin] << ")" << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  if (unfold_generator == 2)
    {
      boundary_r4[0] = truth_leading_cut;
      boundary_r4[1] = 30;
      boundary_r4[3] = 30;
      boundary_r4[4] = ipt_bins[nbins_pt];;     
    }

  TH1D *h_flat_truth_mapping_primer = nullptr;
  TH1D *h_flat_reco_mapping_primer = nullptr;
  TH2D *h_flat_response_mapping_primer[5] = {0};

  int nbins_pt_truth = nbins_pt*nbins_pt;
  std::map<int, int> mapped_pt_bin_truth;
  int nbins_pt_reco = nbins_pt*nbins_pt;
  std::map<int, int> mapped_pt_bin_reco;

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
      
  for (int i = 0; i < 5; i++)
    {
      h_flat_response_mapping_primer[i] = (TH2D*) fr->Get(Form("h_flat_response_mapping_%d", i));
    }      


  int nbins_pt_2 = nbins_pt*nbins_pt;
  int nbins_pt_reco_2 = nbins_pt_reco;
  int nbins_pt_truth_2 = nbins_pt_truth;
  std::cout <<  "nbins_reco: " << nbins_pt_reco_2 << std::endl;
  std::cout <<  "nbins_truth: " << nbins_pt_truth_2 << std::endl;
  TH1D *h_flat_truth_pt1pt2_primer1 = nullptr;

  if (primer != 1)
    {
      TFile *fresponse1 = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()), "r");
      if (!fresponse1)
	{
	  std::cout << "No file for mbd efficiency" << std::endl;
	  return 1;
	}

      std::cout << Form("response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name_orig.c_str()) << std::endl;
      h_flat_truth_pt1pt2_primer1 = (TH1D*) fresponse1->Get("h_truth_flat_pt1pt2"); 

      h_flat_truth_pt1pt2_primer1->SetName("h_truth_flat_pt1pt2_primer1");

      if (!h_flat_truth_pt1pt2_primer1)
	{
	  std::cout << "no truth" << std::endl;
	  return 1;
	}
    }

  TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);

  responsepath += "_" + sys_name;
  
  responsepath += ".root";

  
  TFile *fresponse = new TFile(responsepath.Data(),"r");
  
  RooUnfoldResponse *rooResponse = (RooUnfoldResponse*) fresponse->Get("response");

  if (!rooResponse)
    {
      std::cout << "no repsonse" << std::endl;
      return 1;
    }

  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fresponse->Get("h_truth_flat_pt1pt2");//h_truth_flat_pt1pt2"); 
  if (!h_flat_truth_pt1pt2)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }

  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fresponse->Get("h_reco_flat_pt1pt2"); 
  if (!h_flat_reco_pt1pt2)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }

  TH1D *h_filled_xj_truth[3];
  for (int i =0 ; i < 3; i++)
    {
      h_filled_xj_truth[i] = (TH1D*) fresponse->Get(Form("h_truth_xj_range_%d", i));
      if (!h_filled_xj_truth[i])
	{
	  std::cout << "no truth xj" << std::endl;
	  return 1;
	}

    }

  TH1D *h_reco_lead = (TH1D*) fresponse->Get("h_reco_lead");
  if (!h_reco_lead)
    {
      std::cout << "no reco lead" << std::endl;
      return 1;
    }
  TH1D *h_reco_sublead = (TH1D*) fresponse->Get("h_reco_sublead");
  if (!h_reco_sublead)
    {
      std::cout << "no reco sublead" << std::endl;
      return 1;
    }
  TH2D *h_reco_eta_phi_lead = (TH2D*) fresponse->Get("h_eta_phi_lead");
  if (!h_reco_eta_phi_lead)
    {
      std::cout << "no reco sublead" << std::endl;
      return 1;
    }
  h_reco_eta_phi_lead->SetName("h_reco_eta_phi_lead");

  TH2D *h_eta_lead_sublead = new TH2D("h_eta_lead_sublead","", 220, -1.1, 1.1, 220, -1.1, 1.1);
  TH2D *h_eta_phi_lead = new TH2D("h_eta_phi_lead","", 48, -1.1, 1.1, 128, -TMath::Pi(), TMath::Pi());

  TH2D *h_lead_pt_v_emfrac = new TH2D("h_lead_pt_v_emfrac","",nbins_pt, dpt_bins, 12, -0.1, 1.1);
  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_njets = new TH1D("h_njets", ";N_{Jet}; counts", 21, -0.5, 20.5);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";x_{J};1/N", nbins, ixj_bins);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T1, data};p_{T2, data}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T1, smear} + p_{T2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH1D *h_flat_data_pt1pt2_eff = new TH1D("h_data_flat_pt1pt2_eff",";p_{T1, smear} + p_{T2, smear}", nbins_pt_reco_2, 0, nbins_pt_reco_2);
  TH1D *h_data_lead = new TH1D("h_data_lead", ";p_{T} [GeV];", 100, 0, 100);
  TH1D *h_data_sublead = new TH1D("h_data_sublead", ";p_{T} [GeV];", 100, 0, 100);
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);

  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets;
  for (int isample = 0; isample < 5; isample++)
    {

      if (!sample_skip[isample]) continue;
      int entries = ttree[isample]->GetEntries();

      int start = 0;
      if (full_or_half) entries /= 2;// entries/2;
      int last_runnumber = 0;
      bool skip_run = false;
      for (int i = start; i < entries; i++)
	{

	  ttree[isample]->GetEntry(i);

	  if (last_runnumber != runnumber[isample])
	    {
	      last_runnumber = runnumber[isample];
	      if (std::find(runnumbers.begin(), runnumbers.end(), runnumber[isample]) == runnumbers.end())
		{
		  std::cout << "Skipping run " << runnumber[isample] << std::endl;
		  skip_run = true;
		}
	      else
		{
		  std::cout << "Analysing run " << runnumber[isample] << std::endl;
		  skip_run = false;
		}
	    }
	  if (skip_run) continue;
	  
	  if (runnumber[isample] >= runnumber_high || runnumber[isample] < runnumber_low) continue;
	  std::cout << "Event : " << i << " \r" << std::flush;
	  myrecojets.clear();
	  mytruthjets.clear();

	  bool has_vertex = (fabs(mbd_vertex_z[isample]) < vtx_cut);
	  
	  bool has_mbd_hit = mbd_hit[isample];
	  if (!has_vertex || !has_mbd_hit)
	    continue;


	  //      if ( !( ((gl1_scaled >> 22) & 0x1 ) == 0x1 ||
	  bool triggered = true;
	  if (unfold_generator)
	    {
	      float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	      if (cone_size != 4)
		{
		  maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
		}

	      if (maxpttruth < boundary_r4[isample] || maxpttruth >= boundary_r4[isample+1]) continue;
	      float maxptreco = *std::max_element(reco_jet_pt_uncalib[isample]->begin(), reco_jet_pt_uncalib[isample]->end());

	      if (maxptreco >= maxreco_r4[isample]) continue;

	      
	      triggered = ((gl1_scaled[isample] >> 34) & 0x1) == 0x1;

	      h_truth_lead->Fill(maxpttruth, scale_factor[isample]);
	    
	    }
	  else
	    {
	      
	      triggered = ((gl1_scaled[isample] >> 18) & 0x1)  == 0x1 ||((gl1_scaled[isample] >> 34) & 0x1)  == 0x1;
	      //if (mbd_avgsigma[isample] >= 0.4) continue;
	    }
      

	  int nrecojets = reco_jet_pt[isample]->size();

	  int nnrecojets = 0;
	  for (int j = 0; j < nrecojets;j++)
	    {

	      if (reco_jet_e[isample]->at(j) < 0) continue;


	      
	      struct jet tempjet;
	      tempjet.istruth = 0;
	      float temppt = reco_jet_pt[isample]->at(j);
	      tempjet.pt = temppt;
	      if (unfold_generator)
		{
		  int ib = floor((temppt - 3)/0.1) + 1;
		  float smear1 = 0.08;//hjersmear->GetBinContent(ib);
		  float jersmear = 0;
		  if (smear1 > 0)
		    {
		      fsmear->SetParameter(2, smear1);
		      jersmear = fsmear->GetRandom();
		    }
		  float tempptreco = temppt + jersmear*temppt;
		  tempjet.pt = tempptreco;
		}
	      if (tempjet.pt > 5) nnrecojets++;
	      if (tempjet.pt < reco_subleading_cut) continue;
	      // if (fabs(reco_jet_eta_det[isample]->at(j)) > 0.7) continue;
	      // if (fabs(reco_jet_eta[isample]->at(j)) > 0.7) continue;
	      tempjet.e = reco_jet_e[isample]->at(j);
	      tempjet.emcal = reco_jet_emcal[isample]->at(j);
	      tempjet.eta = reco_jet_eta[isample]->at(j);
	      tempjet.eta_det = reco_jet_eta_det[isample]->at(j);
	      tempjet.phi = reco_jet_phi[isample]->at(j);
	      tempjet.t = 0;
	      tempjet.id = j;

	      myrecojets.push_back(tempjet);
	    }

	  bool njet_result = (myrecojets.size() < njet_cut);

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  bool reco_good = djf.check_dijet_reco(myrecojets);
	  bool time_good = djf.passes_time_cut(calib_lead_time, calib_delta_time);
	  reco_good &= time_good;
	  reco_good &= triggered;
	  reco_good &= has_mbd_hit;
	  reco_good &= has_vertex;
	  reco_good &= njet_result;
	  
	  if (reco_good)
	    {
	      reco_good &= (myrecojets.at(0).pt < ipt_bins[max_reco_bin]);
	    }
      
	  if (!reco_good) continue;

	  float pt1_reco_bin = nbins_pt;
	  float pt2_reco_bin = nbins_pt;

	  float maxi = myrecojets.at(0).pt;
	  float mini = myrecojets.at(1).pt;
	  float es1 = maxi;
	  float es2 = mini;

	  float trigger_eff = 1./ftrigger->Eval(maxi, mini);
	  
	  //if () continue;

	  for (int ib = 0; ib < nbins_pt; ib++)
	    {
	      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
		{
		  pt1_reco_bin = ib;
		}
	      if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
		{
		  pt2_reco_bin = ib;
		}
	    }


	  if (pt1_reco_bin == nbins_pt) continue;
	  if (pt2_reco_bin == nbins_pt) continue;

	  int total_bin_reco_1 = pt1_reco_bin + nbins_pt*pt2_reco_bin;
	  int total_bin_reco_2 = pt2_reco_bin + nbins_pt*pt1_reco_bin;

	  if (mapped_pt_bin_reco[pt1_reco_bin*nbins_pt + pt2_reco_bin] == -1)
	    {
	      reco_good = false;
	    }
	  else
	    {
	      total_bin_reco_1 = mapped_pt_bin_reco[total_bin_reco_1];
	      total_bin_reco_2 = mapped_pt_bin_reco[total_bin_reco_2];
	    }

	  if (!reco_good) continue;
	  
	  h_pt1pt2->Fill(es1, es2, scale_factor[isample]);
	  h_pt1pt2->Fill(es2, es1, scale_factor[isample]);
	  if (!unfold_generator)
	    {
	      if (trigger_sys)
		{
		  h_flat_data_pt1pt2->Fill(total_bin_reco_1, 1./trigger_eff);
		  h_flat_data_pt1pt2->Fill(total_bin_reco_2, 1./trigger_eff);
		}
	      else
		{
		  h_flat_data_pt1pt2->Fill(total_bin_reco_1);
		  h_flat_data_pt1pt2->Fill(total_bin_reco_2);

		}	      
	      h_flat_data_pt1pt2_eff->Fill(total_bin_reco_1);
	      h_flat_data_pt1pt2_eff->Fill(total_bin_reco_2);
	    }
	  else
	    {
	      h_flat_data_pt1pt2->Fill(total_bin_reco_1, scale_factor[isample]);
	      h_flat_data_pt1pt2->Fill(total_bin_reco_2, scale_factor[isample]);

	    }
	  if (maxi >= ipt_bins[measure_bins[0]] && maxi < ipt_bins[measure_bins[2]] && mini >= measure_subleading_cut)
	    {
	      h_reco_xj->Fill(mini/maxi, scale_factor[isample]);
	    }
	  h_mbd_vertex->Fill(mbd_vertex_z[isample], scale_factor[isample]);
	  h_njets->Fill(nnrecojets, scale_factor[isample]);
	  h_data_lead->Fill(maxi, scale_factor[isample]);
	  h_eta_phi_lead->Fill(myrecojets.at(0).eta, myrecojets.at(0).phi, scale_factor[isample]);
	  h_lead_pt_v_emfrac->Fill(myrecojets.at(0).pt, myrecojets.at(0).emcal,  scale_factor[isample]); 
	  h_data_sublead->Fill(mini, scale_factor[isample]);
	  h_eta_lead_sublead->Fill(myrecojets.at(0).eta, myrecojets.at(1).eta, scale_factor[isample]);
	}
    }

  h_flat_data_pt1pt2->Scale(.5);
  h_flat_data_pt1pt2_eff->Scale(.5);
  //h_flat_data_pt1pt2->Multiply(h_purity);
  TH1D *h_flat_data_pt1pt2_full = new TH1D("h_data_flat_pt1pt2_full",";p_{T1, smear} + p_{T2, smear}", nbins_pt_2, 0, nbins_pt_2);
  for (int ib = 0; ib < nbins_pt_2; ib++)
	{
	  int bin = mapped_pt_bin_reco[ib];
	  if (bin >= 0)
	    {
	      h_flat_data_pt1pt2_full->SetBinContent(ib+1, h_flat_data_pt1pt2->GetBinContent(bin+1));
	      h_flat_data_pt1pt2_full->SetBinError(ib+1, h_flat_data_pt1pt2->GetBinError(bin+1));
	    }
	}

  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 1;
  for (int iter = 0; iter < niterations; iter++ )
    {
     
      RooUnfoldBayes   unfold (rooResponse, h_flat_data_pt1pt2, iter + 1, false, true);    // OR
      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
      std::cout << __LINE__ << std::endl;
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
    }
    
  
  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_truth_primer1 = new TH2D("h_pt1pt2_truth_primer1", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_sim = new TH2D("h_pt1pt2_sim", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_data_all = new TH1D("h_xj_data_all", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_truth_primer1 = new TH1D("h_xj_truth_primer1", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_sim = new TH1D("h_xj_sim", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins_pt);
  if (primer != 1)
    {
      histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2_primer1, h_pt1pt2_truth_primer1, nbins_pt);
    }
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_sim, nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2_full, h_pt1pt2_data, nbins_pt);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins_pt);
    }

  h_pt1pt2_data->SetTitle(";Data p_{T, 1} [GeV]; Data p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  gStyle->SetPaintTextFormat("4.0f");

  TH2D *hct = (TH2D*) h_pt1pt2_truth->Clone();
  hct->SetName("hct");
  hct->Scale(1./hct->Integral());

  TCanvas *ceff = new TCanvas("ceff","ceff", 1000, 300);
  ceff->SetLogy();
  h_flat_data_pt1pt2->Draw();
  h_flat_data_pt1pt2_eff->SetLineColor(kRed);
  h_flat_data_pt1pt2_eff->Draw("same");

  ceff->Print(Form("%s/unfolding_plots/pt1pt2_eff_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, sys_name.c_str()));
  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 2000, 2000);
  cpt1pt2->Divide(2, 2);

  for (int i = 0; i < niterations; i++)
    {
      h_pt1pt2_unfold[i]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");
      
      cpt1pt2->cd(1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.2);
      h_pt1pt2_data->Draw("colz");
      
      cpt1pt2->cd(2);
      gPad->SetRightMargin(0.2);
      gPad->SetLogz();
      h_pt1pt2_truth->SetMinimum(1);
      h_pt1pt2_truth->Draw("colz");
      
      cpt1pt2->cd(3);
      gPad->SetRightMargin(0.2);
      gPad->SetLogz();
      h_pt1pt2_unfold[i]->SetMinimum(0.01);//Draw("colz");
      h_pt1pt2_unfold[i]->Draw("colz");

      TH2D *hc = (TH2D*) h_pt1pt2_unfold[i]->Clone();
      hc->SetName("hc");
      hc->Scale(1./hc->Integral());
      hc->Divide(hct);

      cpt1pt2->cd(4);
      gPad->SetRightMargin(0.2);
      gPad->SetLogz(0);
      hc->Draw("colz");
      
      cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2data_%s_r%02d_iter_%d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, i, sys_name.c_str()));
    }
  // cpt1pt2->cd(1);

  // h_pt1pt2_raw->GetYaxis()->SetTitleSize(0.06);
  // h_pt1pt2_raw->GetXaxis()->SetTitleSize(0.06);

  // h_pt1pt2_raw->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2_raw->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  // TLine *id_l = new TLine(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45, 45);
  // dlutility::SetLineAtt(id_l, kBlack, 2, 2);
  // h_pt1pt2_raw->Draw("colz");
  // id_l->Draw("same");
  // if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  // else dlutility::DrawSPHENIXpp(0.22, 0.85);
  // cpt1pt2->cd(2);
  // h_pt1pt2_ZYAM->GetYaxis()->SetTitleSize(0.06);
  // h_pt1pt2_ZYAM->GetXaxis()->SetTitleSize(0.06);
  // h_pt1pt2_ZYAM->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2_ZYAM->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2_ZYAM->Draw("colz");
  // id_l->Draw("same");
  // if (!ispp) dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.85);

  // cpt1pt2->cd(3);
  // h_pt1pt2->GetYaxis()->SetTitleSize(0.06);
  // h_pt1pt2->GetXaxis()->SetTitleSize(0.06);

  // h_pt1pt2->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2->Draw("colz");
  // id_l->Draw("same");
  // cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2dataZYAM_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, sys_name.c_str()));

  // TCanvas *ctext = new TCanvas("ctext","ctext", 600, 500);
  // ctext->SetRightMargin(0.2);
  
  // TH2D *h_pt1pt2_text = (TH2D*) h_pt1pt2_ZYAM->Clone();
  // h_pt1pt2_text->Divide(h_pt1pt2_raw);

  // h_pt1pt2_text->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  // h_pt1pt2_text->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);

  // gStyle->SetPaintTextFormat("1.2f");
  // h_pt1pt2_text->Draw("colz");
  // h_pt1pt2_text->Draw("text,same");
  // id_l->Draw("same");
  // if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  // else dlutility::DrawSPHENIXpp(0.22, 0.85);
  // //  dlutility::DrawSPHENIX(0.22, 0.85);
  // if (!ispp) dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.75);
  // dlutility::drawText("Background Fraction", 0.22, 0.7);
  // ctext->Print(Form("%s/unfolding_plots/pt1pt2dataTEXT_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, sys_name.c_str()));


  // Define a diverging color palette for TH2D:
  // Blue (low) → Green (central band) → Red (high)

  
  // Color control points (0–1)

  // const unsigned char Number = 10;

  // Double_t Red[Number]   =  {0, 0, 0.0,     0.0, 1,       1, 1, 1, 1, 1};
  // Double_t Green[Number] =  {0, 0, 1,     1, 1,       1, .5, .5, 0, 0};
  // Double_t Blue[Number]  =  {1, 1, 0.0,     0.0, 0.0,       0.0, 0, 0, 0, 0};
  // Double_t  Stops[Number] = {0.00, 0.15, 0.15001, 0.25, 0.2501, 0.4, 0.4001, 0.8, 0.8001, 1.0};

  // TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 20);
  gStyle->SetPalette(kRainBow);
  TCanvas *ceta = new TCanvas("ceta","ceta", 1000, 500);
  ceta->Divide(2,1);
  ceta->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.15);
  h_eta_phi_lead->SetTitle(";#eta;#phi;Counts");
  h_eta_phi_lead->Draw("colz");
  dlutility::DrawSPHENIXpp(0.15, 0.95, 0.03, 0, 1, 1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.15, 0.9, 0, kBlack, 0.03);
  dlutility::drawText("Data", 0.15, 0.85, 0, kBlack, 0.03);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.85, 0.9, 1, kBlack, 0.03);
  dlutility::drawText("|#eta| < 1.1 - #it{R}", 0.85, 0.85, 1, kBlack, 0.03);


  ceta->cd(2);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.15);
  h_reco_eta_phi_lead->SetTitle(";#eta;#phi;Counts");
  h_reco_eta_phi_lead->Draw("colz");
  dlutility::DrawSPHENIXpp(0.15, 0.95, 0.03, 0, 1, 1, 1, "PYTHIA-8");

  dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[reco_leading_bin], ipt_bins[nbins]), 0.15, 0.9, 0, kBlack, 0.03);
  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[reco_subleading_bin]), 0.15, 0.85, 0, kBlack, 0.03);
  dlutility::drawText("Simulation", 0.85, 0.85, 1, kBlack, 0.03);
  ceta->Print(Form("%s/unfolding_plots/etaphi_comp_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
  TCanvas *ccomp = new TCanvas("ccomp","ccomp", 500, 500);

  h_data_lead->Scale(1./h_data_lead->Integral(), "width");
  h_data_sublead->Scale(1./h_data_sublead->Integral(), "width");
  h_reco_lead->Scale(1./h_reco_lead->Integral(), "width");
  h_reco_sublead->Scale(1./h_reco_sublead->Integral(), "width");

  dlutility::SetMarkerAtt(h_data_lead, kBlue, 1, 20);
  dlutility::SetMarkerAtt(h_data_sublead, kBlue, 1, 24);
  dlutility::SetMarkerAtt(h_reco_lead, kRed, 1, 21);
  dlutility::SetMarkerAtt(h_reco_sublead, kRed, 1, 25);
  dlutility::SetLineAtt(h_data_lead, kBlue, 1, 1);
  dlutility::SetLineAtt(h_data_sublead, kBlue, 1, 1);
  dlutility::SetLineAtt(h_reco_lead, kRed, 1, 1);
  dlutility::SetLineAtt(h_reco_sublead, kRed, 1, 1);

  h_reco_lead->SetTitle(";Jet #it{p}_{T}; Arb.");
  h_reco_lead->Draw("p");
  h_reco_sublead->Draw("p same");
  h_data_lead->Draw("p same");
  h_data_sublead->Draw("p same");
  ccomp->SetLogy();
  dlutility::DrawSPHENIXpp(0.56, 0.89, 0.03);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.56, 0.79, 0, kBlack, 0.03);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[reco_leading_bin], ipt_bins[nbins]), 0.56, 0.75, 0, kBlack, 0.03);
  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[reco_subleading_bin]), 0.56, 0.71, 0, kBlack, 0.03);

  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.56, 0.67, 0, kBlack, 0.03);
  dlutility::drawText("|#eta| < 1.1 - #it{R}", 0.56, 0.63, 0, kBlack, 0.03);

  TLegend *ldr = new TLegend(0.6, 0.46, 0.75, 0.61);
  ldr->SetTextSize(0.03);
  ldr->SetLineWidth(0);
  ldr->AddEntry(h_data_lead, "Data - Lead","p");
  ldr->AddEntry(h_data_sublead, "Data - Sublead","p");
  ldr->AddEntry(h_reco_lead, "Reco. Sim - Lead","p");
  ldr->AddEntry(h_reco_sublead, "Reco. Sim - Sublead","p");
  ldr->Draw("same");
  ccomp->Print(Form("%s/unfolding_plots/spectra_comp_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
  TCanvas *cpt2 = new TCanvas("cpt2","cpt2", 1000, 600);
  TH2D *h_asym_data = (TH2D*) h_pt1pt2_data->Clone();
  h_asym_data->SetName("h_asym_data");
  h_asym_data->Reset();
  TH2D *h_asym_sim = (TH2D*) h_pt1pt2_sim->Clone();
  h_asym_sim->SetName("h_asym_sim");
  h_asym_sim->Reset();

  histo_opps::make_asym_pt1pt2(h_asym_data, h_pt1pt2_data, nbins_pt);
  histo_opps::make_asym_pt1pt2(h_asym_sim, h_pt1pt2_sim, nbins_pt);

  cpt2->Divide(3,2);
  TLegend *lcc = new TLegend(0.22, 0.40, 0.75, 0.60);
  lcc->SetTextSize(0.06);
  lcc->SetLineWidth(0);

  h_asym_data->Scale(1./h_asym_data->Integral(),"width");
  h_asym_sim->Scale(1./h_asym_sim->Integral(),"width");
  for (int i = 0; i < 6; i++)
    {
      cpt2->cd(i+1);
      TH1D *hd = nullptr;
      TH1D *hs = nullptr;
      float cut_low = 0;
      float cut_high = 100;

      if (i == 0)
	{
	  hd = h_asym_data->ProjectionY(Form("hd_%d", i), reco_leading_bin, measure_bins[0]);
	  hs = h_asym_sim->ProjectionY(Form("hs_%d", i), reco_leading_bin, measure_bins[0]);
	  hd->GetXaxis()->SetRangeUser(0, ipt_bins[measure_bins[0]]);
	  cut_low = reco_leading_cut;
	  cut_high = ipt_bins[measure_bins[0]];
	}
      else if (i == 4)
	{
	  hd = h_asym_data->ProjectionY(Form("hd_%d", i), measure_bins[3], -1);
	  hs = h_asym_sim->ProjectionY(Form("hs_%d", i), measure_bins[3], -1);
	  cut_low = ipt_bins[measure_bins[3]];
	  cut_high = ipt_bins[nbins_pt];
	}
      else if (i == 5)
	{
	  hd = h_asym_data->ProjectionX(Form("hd_%d", i), 0, -1);
	  hs = h_asym_sim->ProjectionX(Form("hs_%d", i), 0, -1);
	}
      else
	{
	  hd = h_asym_data->ProjectionY(Form("hd_%d", i), measure_bins[i - 1], measure_bins[i]);
	  hs = h_asym_sim->ProjectionY(Form("hs_%d", i), measure_bins[i - 1], measure_bins[i]);
	  hd->GetXaxis()->SetRangeUser(0, ipt_bins[measure_bins[i]]);
	  cut_low = ipt_bins[measure_bins[i-1]];
	  cut_high = ipt_bins[measure_bins[i]];

	}

      hd->Scale(1./hd->Integral(),"width");
      hs->Scale(1./hs->Integral(),"width");

      std::string cut_string = Form("%2.1f #geq #it{p}_{T} < %2.1f GeV", cut_low, cut_high);
      hd->SetTitle("; #it{p}_{T,2} [GeV] ; Self. Norm.");
      
      dlutility::SetMarkerAtt(hd, kBlue, 1, 20);
      dlutility::SetLineAtt(hd, kBlue, 2, 1);

      dlutility::SetMarkerAtt(hs, kRed, 0.9, 21);
      dlutility::SetLineAtt(hs, kRed, 2, 1);

      if (i == 5)      
	{
	  gPad->SetLogy();
	}
      hd->Draw("hist E1");
      hs->Draw("hist E1 same");
      if (i < 6)
	dlutility::drawText(cut_string.c_str(), 0.22, 0.8);
      if (i == 0)
	{
	  lcc->AddEntry(hd, "Data");
	  lcc->AddEntry(hs, "Sim. Reco.");
	}
    }
  // cpt2->cd(6);
  // dlutility::DrawSPHENIXpp(0.22, 0.85, 0.05);
  // dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.74,  0, kBlack, 0.05);
  // dlutility::drawText(Form("Systematic: %s" , sys_name.c_str()), 0.22, 0.69,  0, kBlack, 0.05);
  //lcc->Draw("same");
  cpt2->Print(Form("%s/unfolding_plots/subleading_comp_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));

  
  TCanvas *cunpt1pt2 = new TCanvas("cunpt1pt2","cunpt1pt2", 1200, 1000);

  TH2D *h_ratio_data_reco = (TH2D*) h_pt1pt2_data->Clone();

  double pt1_unfold_bins[nbins_pt];
  double pt1_truth_bins[nbins_pt];
  for (int i = 0; i < nbins_pt; i++)
    {
      pt1_unfold_bins[i] = 0;
      pt1_truth_bins[i] = 0;
    }

  h_ratio_data_reco->SetName("h_ratio_data_reco");

  //h_ratio_data_reco->Scale(1./h_ratio_data_reco->Integral());
  TH2D *h_ratio_sim_reco = (TH2D*) h_pt1pt2_sim->Clone();
  h_ratio_sim_reco->SetName("h_ratio_sim_reco");
  //h_ratio_sim_reco->Scale(1./h_ratio_sim_reco->Integral());

  for (int ibin = 0; ibin < nbins_pt_2; ibin++)
    {
      int pt1_bin = ibin/nbins_pt;
      int pt2_bin = ibin%nbins_pt;
      int max_bin = std::max(pt1_bin, pt2_bin);
      int min_bin = std::min(pt1_bin, pt2_bin);

      pt1_truth_bins[max_bin] += h_ratio_sim_reco->GetBinContent(pt1_bin+1, pt2_bin+1);
      pt1_unfold_bins[max_bin] += h_ratio_data_reco->GetBinContent(pt1_bin+1, pt2_bin+1);
    }
  for (int ibin = 0; ibin < nbins_pt_2; ibin++)
    {
      int pt1_bin = ibin/nbins_pt;
      int pt2_bin = ibin%nbins_pt;
      int max_bin = std::max(pt1_bin, pt2_bin);
      int min_bin = std::min(pt1_bin, pt2_bin);
      int gbin = h_ratio_sim_reco->GetBin(pt1_bin+1, pt2_bin+1);
      h_ratio_sim_reco->SetBinContent(gbin,  h_ratio_sim_reco->GetBinContent(pt1_bin+1, pt2_bin+1)/pt1_truth_bins[max_bin]);
      h_ratio_data_reco->SetBinContent(gbin,  h_ratio_data_reco->GetBinContent(pt1_bin+1, pt2_bin+1)/pt1_unfold_bins[max_bin]);
    }

  h_ratio_data_reco->Divide(h_ratio_sim_reco);
  h_ratio_data_reco->SetTitle(";#it{p}_{T,1};#it{p}_{T,1};Data/Reco");
  
  cunpt1pt2->SetRightMargin(0.2);
  //h_ratio_data_reco->SetMinimum(0);
  //h_ratio_data_reco->SetMaximum(2);
  if (unfold_generator)
    {
      //h_ratio_data_reco->SetMinimum(0.95);
      //h_ratio_data_reco->SetMaximum(1.05);
    }
  h_ratio_data_reco->Draw("colz");

  cunpt1pt2->Print(Form("%s/unfolding_plots/full_unfold_pt1pt2_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));


  histo_opps::normalize_histo(h_reco_xj, nbins);


  TH2D *h_pt1pt2_data_truth = new TH2D("h_pt1pt2_data_truth", ";p_{T1};p_{T2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH1D *h_xj_data_truth = new TH1D("h_xj_data_truth", ";x_{J};",nbins, ixj_bins);
  if (unfold_generator > 0)
    {
      histo_opps::make_sym_pt1pt2(h_flat_data_truth_pt1pt2, h_pt1pt2_data_truth, nbins_pt);
    }
  histo_opps::project_xj(h_pt1pt2_data, h_xj_data_all, nbins_pt, measure_bins[0], measure_bins[2], measure_subleading_bin, nbins - 1);
  TCanvas *cunfold = new TCanvas("cunfold","cunfold", 700, 1000);
  dlutility::ratioPanelCanvas(cunfold);
  
  for (int irange = 0; irange < 3; irange++)
    {
      h_xj_data->Reset();
      h_xj_truth->Reset();
      if (unfold_generator > 0)
	{
	  h_xj_data_truth->Reset();
	}
      if (primer != 1)
	{
	  h_xj_truth_primer1->Reset();
	}
      h_xj_sim->Reset();

      histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);

      if (primer != 1)
	{
	  histo_opps::project_xj(h_pt1pt2_truth_primer1, h_xj_truth_primer1, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
	  histo_opps::normalize_histo(h_xj_truth_primer1, nbins);
	}

      if (unfold_generator > 0)
	{
	  histo_opps::project_xj(h_pt1pt2_data_truth, h_xj_data_truth, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
	  histo_opps::normalize_histo(h_xj_data_truth, nbins);
	}

      histo_opps::project_xj(h_pt1pt2_sim, h_xj_sim, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
      histo_opps::normalize_histo(h_xj_sim, nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold[iter]->Reset();
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
	  histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
	}
      histo_opps::normalize_histo(h_xj_data, nbins);
      histo_opps::normalize_histo(h_xj_truth, nbins);





      dlutility::SetMarkerAtt(h_xj_truth, kBlack, 1, 8);
      dlutility::SetLineAtt(h_xj_truth, kBlack, 1, 1);
      if (primer != 1)
	{

	  dlutility::SetMarkerAtt(h_xj_truth_primer1, kBlack, 1, 24);
	  dlutility::SetLineAtt(h_xj_truth_primer1, kBlack, 1, 1);
	}
      if (unfold_generator > 0)
	{

	  dlutility::SetMarkerAtt(h_xj_data_truth, kGreen+1, 2, 24);
	  dlutility::SetLineAtt(h_xj_data_truth, kGreen+1, 1, 1);
	}

      dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 24);
      dlutility::SetLineAtt(h_xj_data, kBlue, 1, 1);

      dlutility::SetMarkerAtt(h_xj_sim, kRed, 1, 24);
      dlutility::SetLineAtt(h_xj_sim, kRed, 1, 1);

      dlutility::SetLineAtt(h_filled_xj_truth[irange], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_filled_xj_truth[irange], kBlack, 1, 25);

      h_xj_truth->SetMaximum(h_xj_unfold[0]->GetBinContent(h_xj_unfold[0]->GetMaximumBin())*1.3);
      h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
      cunfold->cd(1);
      
      h_xj_truth->Draw();
      //h_filled_xj_truth[irange]->Draw("same");
      
      if (primer != 1)
	{
	  h_xj_truth_primer1->Draw("same");
	}
      if (unfold_generator > 0)
	{
	  h_xj_data_truth->Draw("same");
	}

      h_xj_data->Draw("same");
      h_xj_sim->Draw("same");

      for (int iter = 0; iter < niterations; iter++)
	{
	  if (iter == 0)
	    {
	      dlutility::SetMarkerAtt(h_xj_unfold[iter], kViolet, 1, 1);
	      dlutility::SetLineAtt(h_xj_unfold[iter], kViolet, 2, 1);
	    }
	  else if (iter == niterations - 1)
	    {
	      dlutility::SetMarkerAtt(h_xj_unfold[iter], kGreen, 1, 1);
	      dlutility::SetLineAtt(h_xj_unfold[iter], kGreen, 2, 1);

	    }
	  else
	    {
	      dlutility::SetMarkerAtt(h_xj_unfold[iter], kBlack, 1, 1);
	      dlutility::SetLineAtt(h_xj_unfold[iter], kBlack, 2, 1);
	    }
	  h_xj_unfold[iter]->Draw("same hist");
	}

      dlutility::DrawSPHENIXpp(0.22, 0.85, 0.03);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.75, 0, kBlack, 0.03);
      dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", measure_subleading_cut), 0.22, 0.7, 0, kBlack, 0.03);
      TLegend *lc = new TLegend(0.55, 0.68, 0.75, 0.88);
      lc->SetTextSize(0.03);
      lc->SetLineWidth(0);
      lc->AddEntry(h_xj_data,"Raw Data");
      lc->AddEntry(h_xj_sim,"Sim. Reco.");
      lc->AddEntry(h_xj_unfold[0],"Unfold Iter = 1");
      lc->AddEntry(h_xj_unfold[9],"Unfold Iter = 10");
      lc->AddEntry(h_xj_unfold[1],"Unfold Intermediate Iter");
      if (herwig_sys)
	{
	  lc->AddEntry(h_xj_truth,"Truth Reweighted");
	  if (primer != 1)
	    {
	      lc->AddEntry(h_xj_truth_primer1,"Truth Herwig 7.8");
	    }
	}
      else
	{
	  lc->AddEntry(h_xj_truth,"Truth Reweighted");
	  if (primer != 1)
	    {
	      lc->AddEntry(h_xj_truth_primer1,"Truth PYTHIA8");
	    }

	}
      if (unfold_generator > 0)
	{
	  lc->AddEntry(h_xj_data_truth,Form("Truth %s", generator_name[unfold_generator].c_str()));
	}

      lc->Draw("same");

      cunfold->cd(2);

      TH1D *h_iter_compare[niterations];
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_iter_compare[iter] = (TH1D*) h_xj_unfold[iter]->Clone();

	  if (unfold_generator > 0)
	    {
	      h_iter_compare[iter]->Divide(h_xj_data_truth);
	    }
	  else if (primer != 1)
	    {
	      h_iter_compare[iter]->Divide(h_xj_truth_primer1);
	    }
	  else
	    {
	      h_iter_compare[iter]->Divide(h_xj_truth);
	    }
	  if(iter == 0)
	    {

	      h_iter_compare[iter]->SetTitle(";x_{J}; Unfold/Truth");
	      h_iter_compare[iter]->SetMinimum(0.9);
	      h_iter_compare[iter]->SetMaximum(1.1);
	      h_iter_compare[iter]->Draw();
	    }
	  else
	    {
	      h_iter_compare[iter]->Draw("same");
	    }
	}
      cunfold->Print(Form("%s/unfolding_plots/full_unfold_%s_r%02d_range_%d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, sys_name.c_str()));

      histo_opps::normalize_histo(h_xj_truth, nbins);
      if (primer != 1)
	{
	  histo_opps::normalize_histo(h_xj_truth_primer1, nbins);
	}
      if (unfold_generator > 0)
	{
	  histo_opps::normalize_histo(h_xj_data_truth, nbins);
	}

      histo_opps::normalize_histo(h_xj_sim, nbins);
      histo_opps::normalize_histo(h_xj_data, nbins);
      histo_opps::normalize_histo(h_filled_xj_truth[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
	}

    }


  // TCanvas *cdphi = new TCanvas("cphi","cdphi", 500, 500);

  
  
  // h_dphi_reco->SetMinimum(-20);
  // h_dphi_reco->SetMaximum(h_dphi_reco->GetBinContent(h_dphi_reco->GetMaximumBin())*1.7);
  // h_dphi_reco->Draw("hist");
  // h_dphi_reco_ZYAM->Draw("hist same");
  // h_dphi_reco_SIGNAL->Draw("hist same");
  // h_dphi_reco_subtract->Draw("hist same");
  // //h_dphi_reco_subtract_flow->Draw("hist same");
  // //f_hand_modulation->Draw("same");
  // TLine *lii = new TLine(0, 0, TMath::Pi(), 0);
  // lii->Draw("same");
  // TLine *lis = new TLine(0, average_y, TMath::Pi(), average_y);
  // dlutility::SetLineAtt(lis, kRed, 2, 2);
  // lis->Draw("same");

  // dlutility::DrawSPHENIX(0.22, 0.85);
  // dlutility::drawText(Form("p_{T1} #geq %2.1f GeV", reco_leading_cut), 0.22, 0.75);
  // dlutility::drawText(Form("p_{T2} #geq %2.1f GeV", reco_subleading_cut), 0.22, 0.7);
  // dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.65);
  // TLegend *lc = new TLegend(0.5, 0.63, 0.75, 0.78);
  // lc->SetTextSize(0.03);
  // lc->SetLineWidth(0);
  // lc->AddEntry(h_dphi_reco, "No ZYAM sub.");
  // lc->AddEntry(h_dphi_reco_subtract, "w/ ZYAM sub.");
  // //lc->AddEntry(h_dphi_reco_subtract_flow, "w/ ZYAM sub. + flow");
  // lc->AddEntry(h_dphi_reco_SIGNAL, "Signal region");
  // lc->AddEntry(h_dphi_reco_ZYAM, "ZYAM region");
  // lc->Draw("same");
  // cdphi->Print(Form("%s/unfolding_plots/dphi_ZYAM_%s_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin, cone_size));

  // for (int j = 0; j < 4; j++)
  //   {
  //     dlutility::SetLineAtt(h_dphi_reco_sub[j], kBlack, 2, 1);
  //     dlutility::SetLineAtt(h_dphi_reco_ZYAM_sub[j], kBlack, 2, 4);
  //     h_dphi_reco_ZYAM_sub[j]->SetFillColorAlpha(kYellow, 0.4);
  //     dlutility::SetLineAtt(h_dphi_reco_SIGNAL_sub[j], kBlack, 2, 4);
  //     h_dphi_reco_SIGNAL_sub[j]->SetFillColorAlpha(kSpring+2, 0.4);

  //     TH1D *h_dphi_reco_subtract_sub = (TH1D*) h_dphi_reco_sub[j]->Clone();
  //     h_dphi_reco_subtract_sub->SetName("h_dphi_reco_subtract_sub");
  //     dlutility::SetLineAtt(h_dphi_reco_subtract_sub, kRed, 2, 1);

  //     float nbins_x1 = 0;
  //     float average_y1 = 0;
  //     for (int i = 0; i < h_dphi_reco_ZYAM_sub[j]->GetNbinsX(); i++)
  // 	{
  // 	  if (h_dphi_reco_ZYAM_sub[j]->GetBinCenter(i+1) > dphicut_z_low && h_dphi_reco_ZYAM_sub[j]->GetBinCenter(i+1) < dphicut_z_high)
  // 	    {
  // 	      average_y1 += h_dphi_reco_ZYAM_sub[j]->GetBinContent(i+1);
  // 	      nbins_x1 += 1.0;
  // 	    }
  // 	}

  //     average_y1 /= nbins_x1;
  
  //     for (int i = 0; i < h_dphi_reco_subtract_sub->GetNbinsX(); i++)
  // 	{
  // 	  float v = h_dphi_reco_subtract_sub->GetBinContent(i+1);
  // 	  h_dphi_reco_subtract_sub->SetBinContent(i+1, v - average_y1);
  // 	}

  //     h_dphi_reco_sub[j]->SetMinimum(-5);
  //     h_dphi_reco_sub[j]->SetMaximum(h_dphi_reco_sub[j]->GetBinContent(h_dphi_reco_sub[j]->GetMaximumBin())*1.7);
  //     h_dphi_reco_sub[j]->Draw("hist");
  //     h_dphi_reco_ZYAM_sub[j]->Draw("hist same");
  //     h_dphi_reco_SIGNAL_sub[j]->Draw("hist same");
  //     h_dphi_reco_subtract_sub->Draw("hist same");
  //     //h_dphi_reco_subtract_flow->Draw("hist same");
  //     //f_hand_modulation->Draw("same");
  //     //TLine *lii = new TLine(0, 0, TMath::Pi(), 0);
  //     lii->Draw("same");
  //     TLine *lis2 = new TLine(0, average_y1, TMath::Pi(), average_y1);
  //     dlutility::SetLineAtt(lis2, kRed, 2, 2);
  //     lis2->Draw("same");

  //     dlutility::DrawSPHENIX(0.22, 0.85);
  //     dlutility::drawText(Form("p_{T1} #geq %2.1f GeV", reco_leading_cut), 0.22, 0.75);
  //     if (j == 0)
  // 	dlutility::drawText(Form("%2.1f #leq p_{T2}^{lead} < %2.1f GeV", ipt_bins[1], ipt_bins[2]), 0.22, 0.7);
  //     if (j == 1)
  // 	dlutility::drawText(Form("%2.1f #leq p_{T2}^{lead} < %2.1f GeV", ipt_bins[2], ipt_bins[3]), 0.22, 0.7);
  //     if (j == 2)
  // 	dlutility::drawText(Form("%2.1f #leq p_{T2}^{lead} < %2.1f GeV", ipt_bins[3], ipt_bins[5]), 0.22, 0.7);
  //     if (j == 3)
  // 	dlutility::drawText(Form("%2.1f #leq p_{T2}^{lead} < %2.1f GeV", ipt_bins[5], ipt_bins[7]), 0.22, 0.7);

  //     dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.65);
  //     lc = new TLegend(0.65, 0.7, 0.75, 0.85);
  //     lc->SetTextSize(0.03);
  //     lc->SetLineWidth(0);
  //     lc->AddEntry(h_dphi_reco_sub[j], "No ZYAM sub.");
  //     lc->AddEntry(h_dphi_reco_subtract_sub, "w/ ZYAM sub.");
  //     //lc->AddEntry(h_dphi_reco_subtract_flow, "w/ ZYAM sub. + flow");
  //     lc->AddEntry(h_dphi_reco_SIGNAL_sub[j], "Signal region");
  //     lc->AddEntry(h_dphi_reco_ZYAM_sub[j], "ZYAM region");
  //     lc->Draw("same");
  //     cdphi->Print(Form("%s/unfolding_plots/dphi_ZYAM_AA_sub_%d_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), j, centrality_bin, cone_size));

  //   }

  
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 1, 24);


  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  //h_xj_truth->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_truth->SetMaximum(6);
  h_xj_truth->SetMinimum(0);
  h_xj_truth->Draw("p");
  h_xj_unfold[niter]->Draw("same p");
  h_xj_data->Draw("hist same");
  h_reco_xj->Draw("hist same");
  h_xj_truth->Draw("same hist");
  h_xj_data->Draw("p same");
  h_xj_unfold[niter]->Draw("same hist");
  h_xj_unfold[niter]->Draw("same p");
  if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  else dlutility::DrawSPHENIXpp(0.22, 0.85);
  //dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/3", 0.22, 0.59);
  dlutility::drawText("\\mathscr{L} = 25.7 pb^{-1}", 0.22, 0.54);

  TLegend *leg = new TLegend(0.22, 0.35, 0.4, 0.50);
  leg->SetLineWidth(0);
  leg->AddEntry(h_xj_data, "Data");
  if (herwig_sys)
    {
      leg->AddEntry(h_xj_truth, "Herwig 7.8");
    }
  else
    {
      leg->AddEntry(h_xj_truth, "Pythia8");
    }
  leg->AddEntry(h_xj_unfold[niter], "Unfolded");
  leg->Draw("same");
  dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
    
  cxj->cd(2);

  TH1D *h_data_compare = (TH1D*) h_xj_unfold[niter]->Clone();
  h_data_compare->Divide(h_xj_truth);
  //h_data_compare->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_data_compare->SetTitle(";x_{J}; Unfold / Truth");
  dlutility::SetFont(h_data_compare, 42, 0.06);
  dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);

 
  h_data_compare->Draw("p");
  TLine *line = new TLine(0.3, 1, 1, 1);
  line->SetLineStyle(4);
  line->SetLineColor(kRed + 3);
  line->SetLineWidth(2);
  line->Draw("same");



  dlutility::SetLineAtt(h_xj_data_all, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data_all, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlack, 1, 24);
  
  TCanvas *cproj = new TCanvas("cproj","cproj", 500, 700);
  dlutility::ratioPanelCanvas(cproj);

  histo_opps::normalize_histo(h_xj_data_all, nbins);
  cproj->cd(1);
  h_xj_data_all->SetMaximum(3);
  h_xj_data_all->SetTitle(";x_{J}; #frac{1}{N_{pairs}}#frac{dN_{pair}}{dx_{J}}");
  dlutility::SetFont(h_xj_data_all, 42, 0.04);
  
  h_xj_data_all->Draw();
  h_reco_xj->Draw("same");

  
  if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  else dlutility::DrawSPHENIXpp(0.22, 0.85);
  //  dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
  TLegend *legp = new TLegend(0.63, 0.75, 0.8, 0.85);
  legp->SetLineWidth(0);
  legp->SetTextFont(42);
  legp->SetTextSize(0.04);
  legp->AddEntry(h_reco_xj, "Data Filled");
  legp->AddEntry(h_xj_data_all, "Data Projected");
  legp->Draw("same");

  cproj->cd(2);
  TH1D *h_fillproj_compare = (TH1D*) h_xj_data_all->Clone();

  dlutility::SetFont(h_fillproj_compare, 42, 0.06);
  h_fillproj_compare->Divide(h_reco_xj);

  h_fillproj_compare->SetMaximum(1.2);
  
  h_fillproj_compare->SetMinimum(0.8);
  h_fillproj_compare->SetTitle(";x_{J};Projected / Filled");
  dlutility::SetLineAtt(h_fillproj_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_fillproj_compare, kBlack, 1, 8);

  h_fillproj_compare->Draw();
  TLine *line3 = new TLine(0.1, 1, 1, 1);
  line3->SetLineStyle(4);
  line3->SetLineColor(kRed + 3);
  line3->SetLineWidth(2);
  line3->Draw("same");
  cproj->Print(Form("%s/unfolding_plots/proj_compar_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
  
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_" + system_string + "_r0" + std::to_string(cone_size);

  unfoldpath += "_" + sys_name;
  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");
  h_lead_pt_v_emfrac->Write();
  h_eta_phi_lead->Write();
  h_eta_lead_sublead->Write();
  h_flat_data_pt1pt2->Write();
  h_flat_data_pt1pt2_full->Write();
  h_data_lead->Write();
  h_data_sublead->Write();
  h_truth_lead->Write();
  //h_flat_data_pt1pt2_raw->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_mbd_vertex->Write();
  h_njets->Write();
  
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_unfold_pt1pt2[iter]->Write();
    }
  fout->Close();
  

  return 0;
  
}
int main(int argc, char *argv[])
{

  std::string config = "binning.config";
  int niterations = 10;
  int cone_size = 4;
  int primer = 0;
  int set=0;
  int unfold_generator = 0;
  int input_generator = 0;
  int half = 0;
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
      else if (arg == "-h" && i + 1 < argc)
	{
	  half = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-v" && i+1 < argc)
	{
	  verbosity = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-ig" && i+1 < argc)
	{
	  input_generator = std::stoi(argv[++i]);
	}
      else if (arg == "-ug" && i+1 < argc)
	{
	  unfold_generator = std::stoi(argv[++i]);
	}
      else
	{
	  std::cerr << "Unknown or incomplete argument: " << arg << "\n";
	  return 1;
	}
    }
  if (set < 3)
    {
      std::cout << "Not enough settings: " << std::endl;
      std::cout << "[usage] : " << std::endl;
      std::cout << "    ./unfoldData_noempty_pp -c binning.config -r 4 -n 10 -p 1 -ig 0 -ug 0" << std::endl;
      return 1;
    }
  
  return unfoldData_noempty_pp(config, niterations, cone_size, primer, unfold_generator, input_generator, half);
  
}
