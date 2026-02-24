#include <iostream>
#include <string>

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "dijetfinder.h"

#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"

static int verbosity = 0;

void makeRawDeltaPhi(const int cone_size = 4, const std::string configfile = "binning.config")
{
  if (verbosity > 1) std::cout << "Starting" << std::endl;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  float boundary_r4[6];
  boundary_r4[0] = 14;
  boundary_r4[1] = 17;
  boundary_r4[2] = 22;
  boundary_r4[3] = 35;
  boundary_r4[4] = 52;
  boundary_r4[5] = 100;
  
  
  read_binning rb(configfile.c_str());

  dijetfinder djf(cone_size);
  djf.SetVerbosity(verbosity);
  
  
  std::string j10_file = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v11_10_ana509_MDC2-00000028.root";
  std::string j15_file = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v11_15_ana509_MDC2-00000028.root";
  std::string j20_file = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v11_20_ana509_MDC2-00000028.root";
  std::string j30_file = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v11_30_ana509_MDC2-00000028.root";
  std::string j50_file = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v11_50_ana509_MDC2-00000028.root";

  float n_events[5];
  TFile *finsim[5];
  
  finsim[0] = new TFile(j10_file.c_str(),"r");
  finsim[1] = new TFile(j15_file.c_str(),"r");
  finsim[2] = new TFile(j20_file.c_str(),"r");
  finsim[3] = new TFile(j30_file.c_str(),"r");
  finsim[4] = new TFile(j50_file.c_str(),"r");

  TTree *ttree[5];
  ULong64_t gl1_scaled[5];

  std::vector<float> *truth_jet_pt_ref[5] = {0};
  std::vector<float> *truth_jet_pt[5] = {0};
  std::vector<float> *truth_jet_eta[5] = {0};
  std::vector<float> *truth_jet_phi[5] = {0};
  
  std::vector<float> *reco_jet_pt[5] = {0};
  std::vector<float> *reco_jet_emcal[5] = {0};
  std::vector<float> *reco_jet_e[5] = {0};
  std::vector<float> *reco_jet_eta[5] = {0};
  std::vector<float> *reco_jet_eta_det[5] = {0};
  std::vector<float> *reco_jet_phi[5] = {0};

  float truth_vertex_z[5];
  float mbd_vertex_z[5];

  for (int j = 0; j < 5; j++)
    {
      
      ttree[j] = (TTree*) finsim[j]->Get("ttree");
      ttree[j]->SetName(Form("ttree_%d", j));
      if (cone_size != 4)
	{
	  ttree[j]->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_ref[j]);
	}
      ttree[j]->SetBranchAddress(Form("truth_jet_pt_%d", cone_size), &truth_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_eta_%d", cone_size), &truth_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_phi_%d", cone_size), &truth_jet_phi[j]);
  
      ttree[j]->SetBranchAddress(Form("jet_pt_%d", cone_size), &reco_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("jet_emcal_%d", cone_size), &reco_jet_emcal[j]);
      ttree[j]->SetBranchAddress(Form("jet_e_%d", cone_size), &reco_jet_e[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_%d", cone_size), &reco_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_det_%d", cone_size), &reco_jet_eta_det[j]);
      ttree[j]->SetBranchAddress(Form("jet_phi_%d", cone_size), &reco_jet_phi[j]);
      ttree[j]->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z[j]);
      ttree[j]->SetBranchAddress("truth_vertex_z", &truth_vertex_z[j]);
      ttree[j]->SetBranchAddress("gl1_scaled", &gl1_scaled[j]);
      n_events[j] = ttree[j]->GetEntries();
    }

  
  float cs_10 = (3.997e6);
  float cs_15 = (4.073e5);
  float cs_20 = 6.218e4;
  float cs_30 = (2.502e3);
  float cs_50 = (7.2695);
  
  float scale_factor[5];
  scale_factor[0] = cs_10/cs_50;
  scale_factor[1] = cs_15/cs_50;
  scale_factor[2] = cs_20/cs_50;
  scale_factor[3] = cs_30/cs_50; 
  scale_factor[4] = 1;

  std::string data_file = rb.get_tntuple_location() + "/TREE_DIJET_SKIM_r0" + std::to_string(cone_size) + "_v8_2_ana509_2024p022_v001_gl10-ana468list.root";
  
  ULong64_t gl1_scaled_data;
  float mbd_vertex_z_data;
  
  std::vector<float> *data_jet_pt = 0;
  std::vector<float> *data_jet_emcal = 0;
  std::vector<float> *data_jet_e = 0;
  std::vector<float> *data_jet_eta = 0;
  std::vector<float> *data_jet_eta_det = 0;
  std::vector<float> *data_jet_phi = 0;

  TFile *fin = new TFile(data_file.c_str(), "r");
  TTree *ttreed  = (TTree*) fin->Get("ttree");;

  if (!ttreed)
    {
      std::cout << " no data "<< std::endl;
    }

  ttreed->SetBranchAddress(Form("jet_pt_%d", cone_size), &data_jet_pt);
  ttreed->SetBranchAddress(Form("jet_emcal_%d", cone_size), &data_jet_emcal);
  ttreed->SetBranchAddress(Form("jet_e_%d", cone_size), &data_jet_e);
  ttreed->SetBranchAddress(Form("jet_eta_%d", cone_size), &data_jet_eta);
  ttreed->SetBranchAddress(Form("jet_eta_det_%d", cone_size), &data_jet_eta_det);
  ttreed->SetBranchAddress(Form("jet_phi_%d", cone_size), &data_jet_phi);
  ttreed->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z_data);
  ttreed->SetBranchAddress("gl1_scaled", &gl1_scaled_data);
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;
 

  const int nbinsdphi  = 16;
  float min_dphi = 3*TMath::Pi()/4.;
  float max_dphi = TMath::Pi();
  float stepdphi = (max_dphi - min_dphi)/(float)nbinsdphi;

  float idphi_bins[nbinsdphi+1];
  for (int i = 0; i < nbinsdphi+1; i++)
    {
      idphi_bins[i] = min_dphi + i*stepdphi;
    }
  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  Int_t max_reco_bin = rb.get_maximum_reco_bin();
  
  djf.setRecoCuts(reco_leading_cut, reco_subleading_cut);
  djf.setTruthCuts(truth_leading_cut, truth_subleading_cut);
  
  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int binranges[10] = {0};
  int binrangesmin[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      binranges[ir] = rb.get_measure_region(ir);
      binrangesmin[ir] = rb.get_subleading_measure_region(ir);
      std::cout << binranges[ir] << " -- " <<  binrangesmin[ir] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  TH1D *h_data_dphi = new TH1D("h_data_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_truth_dphi = new TH1D("h_truth_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_dphi = new TH1D("h_reco_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TH1D *h_truth_match_dphi = new TH1D("h_truth_match_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_match_dphi = new TH1D("h_reco_match_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TProfile *h_sim_match_ddphi = new TProfile("h_sim_match_ddphi", ";<#it{p}_{T}> [GeV]; d#Delta#phi", 25, 10, 60, "s");
  TProfile2D *h2_sim_match_ddphi = new TProfile2D("h2_sim_match_ddphi", ";#it{p}_{T,1} [GeV] ;#it{p}_{T,2} [GeV]; d#Delta#phi", 4, 20, 60, 10, 10, 60, "s");

  TH3D *h_truth_pt1pt2dphi = new TH3D("h_truth_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_pt1pt2dphi = new TH3D("h_reco_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_truth_match_pt1pt2dphi = new TH3D("h_truth_match_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi = new TH3D("h_reco_match_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphitruth = new TH3D("h_reco_match_pt1pt2dphitruth",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_data_pt1pt2dphi = new TH3D("h_data_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  // Data

  TH1D *h_data_dphi_counts = new TH1D("h_data_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_truth_dphi_counts = new TH1D("h_truth_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_dphi_counts = new TH1D("h_reco_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TH1D *h_truth_match_dphi_counts = new TH1D("h_truth_match_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_match_dphi_counts = new TH1D("h_reco_match_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TProfile *h_sim_match_ddphi_counts = new TProfile("h_sim_match_ddphi_counts", ";<#it{p}_{T}> [GeV]; d#Delta#phi", 25, 10, 60, "s");
  TProfile2D *h2_sim_match_ddphi_counts = new TProfile2D("h2_sim_match_ddphi_counts", ";#it{p}_{T,1} [GeV] ;#it{p}_{T,2} [GeV]; d#Delta#phi", 4, 20, 60, 10, 10, 60, "s");

  TH3D *h_truth_pt1pt2dphi_counts = new TH3D("h_truth_pt1pt2dphi_counts",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_pt1pt2dphi_counts = new TH3D("h_reco_pt1pt2dphi_counts",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_truth_match_pt1pt2dphi_counts = new TH3D("h_truth_match_pt1pt2dphi_counts",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi_counts = new TH3D("h_reco_match_pt1pt2dphi_counts",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi_countstruth = new TH3D("h_reco_match_pt1pt2dphi_countstruth",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_data_pt1pt2dphi_counts = new TH3D("h_data_pt1pt2dphi_counts",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  // Data

  TH1D *h_correlated_counts_leading_reco[3];
  TH1D *h_correlated_counts_subleading_reco[3];
  TH1D *h_all_counts_leading_reco[3];
  TH1D *h_all_counts_subleading_reco[3];

  TH1D *h_correlated_counts_leading_truth[3];
  TH1D *h_correlated_counts_subleading_truth[3];
  TH1D *h_all_counts_leading_truth[3];
  TH1D *h_all_counts_subleading_truth[3];

  for (int i = 0; i < 3; i++)
    {
      h_correlated_counts_leading_reco[i] = new TH1D(Form("h_correlated_counts_leading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_subleading_reco[i] = new TH1D(Form("h_correlated_counts_subleading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_leading_reco[i] = new TH1D(Form("h_all_counts_leading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_subleading_reco[i] = new TH1D(Form("h_all_counts_subleading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_leading_truth[i] = new TH1D(Form("h_correlated_counts_leading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_subleading_truth[i] = new TH1D(Form("h_correlated_counts_subleading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_leading_truth[i] = new TH1D(Form("h_all_counts_leading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_subleading_truth[i] = new TH1D(Form("h_all_counts_subleading_truth_%d", i), ";;", nbinsdphi, idphi_bins);

    }

  int entriesd = ttreed->GetEntries();

  TFile *infilejes = new TFile(Form("%s/r%02d/jes/jes_fits_r%02d_0JES.root", rb.get_jesr_location().c_str(), cone_size, cone_size), "r");
  TF1 *fJES = (TF1*) infilejes->Get("fjesprime_2");

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);

  TF1 *gjer = nullptr;

  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TFile *finjer = new TFile(Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)), "r");

  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::string sys_name = "nominal";
  gjer = (TF1*) finjer->Get("jernom");
  if (JER_sys != 0)
    {
      if (JER_sys < 0)
	{
	  gjer = (TF1*) finjer->Get("jerneg");
	  sys_name = "negJER";
	}
      else if (JER_sys > 0)
	{
	  gjer = (TF1*) finjer->Get("jerpos");
	  sys_name = "posJER";
	}
    }
  if (JES_sys != 0)
    {
      if (JES_sys < 0)
	{
	  sys_name = "negJES";
	}
      else if (JES_sys > 0)
	{
	  sys_name = "posJES";
	}
    }

  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  float width = 0.1 + JER_sys;
  fgaus->SetParameters(1, 0, width);

  // Vertex Reweight
  Double_t vtx_cut = rb.get_vtx_cut();
  Int_t vtx_sys = rb.get_vtx_sys();
  TFile *fvtx = new TFile(Form("%s/vertex/vertex_reweight_pp_r%02d.root", rb.get_code_location().c_str(), cone_size),"r");
  TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");

  std::vector<std::pair<float, float>> vertex_scales;

  for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
    {
      vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
    }


  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets2;
  std::vector<struct jet> mytruthjets2;
  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  
  for (int i = 0; i < entriesd; i++)
    {
      ttreed->GetEntry(i);

      std::cout << "Event : " << i << " \r" << std::flush;
      myrecojets.clear();
      
      if ( !(fabs(mbd_vertex_z_data) < vtx_cut)) continue;
      if ( !( ((gl1_scaled_data >> 22) & 0x1 ) == 0x1 || ((gl1_scaled_data >> 18) & 0x1)  == 0x1 ) ) continue;

      // if (nrecojets >= njet_cut) continue;//if (trigger != 18 || trigger != 22) continue;
      int nrecojets = data_jet_pt->size();

      for (int j = 0; j < nrecojets;j++)
	{

	  if (data_jet_e->at(j) < 0) continue;


	      
	  struct jet tempjet;
	  tempjet.istruth = 0;
	  float temppt = fJES->Eval(data_jet_pt->at(j));
	  tempjet.pt = temppt;

	  if (tempjet.pt < reco_subleading_cut) continue;

	  tempjet.emcal = data_jet_emcal->at(j);
	  tempjet.eta = data_jet_eta->at(j);
	  tempjet.eta_det = data_jet_eta_det->at(j);
	  tempjet.phi = data_jet_phi->at(j);
	  tempjet.t = 0;
	  tempjet.id = j;

	  myrecojets.push_back(tempjet);
	}

      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      bool reco_good = djf.check_dijet_reco(myrecojets);

      if (reco_good)
	{
	  reco_good &= (myrecojets.at(0).pt < ipt_bins[max_reco_bin]);
	}

      if (!reco_good) continue;

      float maxi = myrecojets.at(0).pt;
      float mini = myrecojets.at(1).pt;
      float es1 = maxi;
      float es2 = mini;      

      float dphi = myrecojets.at(0).phi - myrecojets.at(1).phi;
      float dphi_data = fabs(dphi);
      if (dphi_data > TMath::Pi())
	{
	  dphi_data = 2*TMath::Pi() - dphi_data;
	}

      h_data_pt1pt2dphi->Fill(es1, es2, dphi_data);

      h_data_dphi->Fill(dphi_data);

    }

  // Sim
  for (int isample = 0; isample < 5; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries2 = ttree[isample]->GetEntries();

      for (int i = 0; i < entries2; i++)
	{
	  ttree[isample]->GetEntry(i);
	  double event_scale = scale_factor[isample];
	  if (verbosity > 5)
	    {
	      std::cout << "------------- Event " << i << " --------------" << std::endl;
	    }
	  
	  bool has_vertex = true;
	  if ( !(fabs(mbd_vertex_z[isample]) < vtx_cut)) continue;

	  bool triggered = ( ( ( gl1_scaled[isample] >> 22 ) & 0x1) == 0x1);
	  //if (inrecojets >= njet_cut) continue;
	  // Vertex Rewieghting
	  bool fill_response = 0;

	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (cone_size != 4)
	    {
	      maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
	    }


	  if (maxpttruth < boundary_r4[isample] || maxpttruth >= boundary_r4[isample+1]) continue;

	  for (int ib = 0; ib < vertex_scales.size(); ib++)
	    {
	      if (mbd_vertex_z[isample] < vertex_scales.at(ib).first)
		{
		  event_scale *= vertex_scales.at(ib).second;
		  break;
		}
	    }
	  if (verbosity > 3)
	    {
	      std::cout << __LINE__ << " : in sample" << std::endl;
	    }
	  // collect the truth and reco jets
	  // fully matched dijets
	  matched_dijets.clear();
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
	  for (int j = 0; j < nrecojets;j++)
	    {

	      if (reco_jet_e[isample]->at(j) < 0) continue;


	      
	      struct jet tempjet;
	      tempjet.istruth = 0;
	      float temppt = fJES->Eval(reco_jet_pt[isample]->at(j));

	      tempjet.pt = temppt;
		      

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

	  if (myrecojets.size() > 1)
	    {
	      myrecojets2 = {myrecojets.begin(), myrecojets.begin()+2};
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

	  matched_dijets = djf.match_dijets(myrecojets, mytruthjets2);

	  std::sort(matched_dijets.begin(), matched_dijets.end(), [] (auto a, auto b) { return a.first.pt > b.first.pt; });

	  bool matched = matched_dijets.size() == 2;

	  if (matched)
	    {
	      myrecojets2.clear();
	      for (int ijet = 0 ; ijet < matched_dijets.size(); ijet++)
		{
		  auto mjet = matched_dijets.at(ijet);
		  float temppt = mjet.first.pt;
		  float smear1 = gjer->Eval(temppt);
		  if (JER_sys < 0)
		    {
		      smear1 -= 0.01;
		    }
		  else if (JER_sys > 0)
		    {
		      smear1 += 0.01;
		    }
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
	  int fill_fake_miss = 0;

	  bool reco_good = djf.check_dijet_reco(myrecojets2);
	  reco_good &= triggered;
	  reco_good &= has_vertex;


	  if (reco_good)
	    {
	      reco_good &= (myrecojets2.at(0).pt < ipt_bins[max_reco_bin]);
	    }

	  if (truth_good)
	    {
	      truth_good &= (mytruthjets2.at(0).pt < ipt_bins[nbins]);
	    }

	  if (reco_good && matched)
	    {
	      max_reco = myrecojets2.at(0).pt;
	      min_reco = myrecojets2.at(1).pt;
	    }	      

	  if (truth_good)
	    {
	      max_truth = mytruthjets2.at(0).pt;
	      min_truth = mytruthjets2.at(1).pt;

	    }

	  float dphi_reco = 0;
	  if (reco_good && matched) 
	    {
	      dphi_reco = fabs(myrecojets2.at(0).phi - myrecojets2.at(1).phi);
	      if (dphi_reco > TMath::Pi())
		{
		  dphi_reco = 2*TMath::Pi() - dphi_reco;
		}
	    }
	  float dphi_truth = 0;
	  if (truth_good)
	    {
	      dphi_truth = fabs(mytruthjets2.at(0).phi - mytruthjets2.at(1).phi);
	      if (dphi_truth > TMath::Pi())
		{
		  dphi_truth = 2*TMath::Pi() - dphi_truth;
		}
	    }
	  if (matched)
	    {
	      float ddphi = dphi_truth - dphi_reco;
	      if (ddphi > TMath::Pi()) ddphi -= 2*TMath::Pi();
	      if (ddphi < -1*TMath::Pi()) ddphi += 2*TMath::Pi();
	      h_sim_match_ddphi->Fill((max_truth +  min_truth)/2., ddphi, event_scale);
	    }
	  if (truth_good)
	    {
	      h_truth_dphi->Fill(dphi_truth, event_scale);
	      h_truth_pt1pt2dphi->Fill(max_truth, min_truth, dphi_truth, event_scale);
	    }
	  if (reco_good && matched)
	    {
	      h_reco_dphi->Fill(dphi_reco, event_scale);
	      h_reco_pt1pt2dphi->Fill(max_reco, min_reco, dphi_reco, event_scale);
	    }

	  if (matched && reco_good && truth_good)
	    {
	      h_truth_match_dphi->Fill(dphi_truth, event_scale);
	      h_reco_match_dphi->Fill(dphi_reco, event_scale);
	      h_truth_match_pt1pt2dphi->Fill(max_truth, min_truth, dphi_truth, event_scale);
	      h_reco_match_pt1pt2dphi->Fill(max_reco, min_reco, dphi_reco, event_scale);
	      h_reco_match_pt1pt2dphitruth->Fill(max_reco, min_reco, dphi_truth, event_scale);
	      h_truth_dphi->Fill(dphi_truth, event_scale);
	      h_reco_dphi->Fill(dphi_reco, event_scale);
	      h_truth_match_dphi->Fill(dphi_truth, event_scale);
	      h_reco_match_dphi->Fill(dphi_reco, event_scale);

	      h_truth_match_dphi_counts->Fill(dphi_truth);
	      h_reco_match_dphi_counts->Fill(dphi_reco);
	      h_truth_match_pt1pt2dphi_counts->Fill(max_truth, min_truth, dphi_truth);
	      h_reco_match_pt1pt2dphi_counts->Fill(max_reco, min_reco, dphi_reco);
	      h_reco_match_pt1pt2dphi_countstruth->Fill(max_reco, min_reco, dphi_truth);
	      h_truth_dphi_counts->Fill(dphi_truth);
	      h_reco_dphi_counts->Fill(dphi_reco);
	      h_truth_match_dphi_counts->Fill(dphi_truth);
	      h_reco_match_dphi_counts->Fill(dphi_reco);

	      bool same_dphi = false;
	      for (int ip = 0; ip < nbinsdphi; ip++)
		{
		  
		  if (dphi_truth >= idphi_bins[ip] && dphi_truth < idphi_bins[ip+1] && dphi_reco >= idphi_bins[ip] && dphi_reco < idphi_bins[ip+1])
		    {
		      same_dphi = true;
		      break;
		    }
		}

	      for (int ilead = 0; ilead < 3; ilead++)
		{
		  int same_leading = 0;
		  if (!(max_truth < ipt_bins[binranges[ilead]] || max_truth >= ipt_bins[binranges[ilead+1]]))
		    {
		      same_leading++;
		      h_all_counts_leading_truth[ilead]->Fill(dphi_truth);
		    }
		  if (!(max_reco < ipt_bins[binranges[ilead]] || max_reco >= ipt_bins[binranges[ilead+1]]))
		    {
		      same_leading++;
		      h_all_counts_leading_reco[ilead]->Fill(dphi_reco);
		    }
		  if (same_dphi && same_leading == 2)
		    {
		      h_correlated_counts_leading_reco[ilead]->Fill(dphi_reco);
		      h_correlated_counts_leading_truth[ilead]->Fill(dphi_truth);
		    }

		  for (int isub = 0; isub < 3; isub++)
		    {
		      int same_subleading = 0;
		      if (!(min_truth < ipt_bins[binrangesmin[isub]] || min_truth >= ipt_bins[binrangesmin[isub+1]]))
			{
			  same_subleading++;
			  h_all_counts_subleading_truth[isub]->Fill(dphi_truth);
			}
		      if (!(min_reco < ipt_bins[binrangesmin[isub]] || min_reco >= ipt_bins[binrangesmin[isub+1]]))
			{
			  same_subleading++;
			  h_all_counts_subleading_reco[isub]->Fill(dphi_reco);
			}
		      if (same_dphi && same_subleading == 2)
			{
			  h_correlated_counts_subleading_reco[isub]->Fill(dphi_reco);
			  h_correlated_counts_subleading_truth[isub]->Fill(dphi_truth);
			}
		    }
		}
	      //h_sim_match_ddphi->Fill((max_truth +  min_truth)/2., dphi_truth - dphi_reco, event_scale);
	      h2_sim_match_ddphi->Fill(max_truth,  min_truth, dphi_truth - dphi_reco, event_scale);
	    }
	}
    }

  TString outpath = "dphihists/dphi_hists_pp_r0" + std::to_string(cone_size);
	  
  outpath += "_" + sys_name;

  outpath += ".root";

  TFile *fout = new TFile(outpath.Data(), "recreate");

  h_sim_match_ddphi->Write();
  h2_sim_match_ddphi->Write();
  h_data_pt1pt2dphi->Write();
  h_reco_match_pt1pt2dphi->Write();
  h_reco_match_pt1pt2dphitruth->Write();
  h_truth_match_pt1pt2dphi->Write();

  h_reco_pt1pt2dphi->Write();
  h_truth_pt1pt2dphi->Write();
  h_reco_dphi->Write();
  h_data_dphi->Write();
  h_truth_dphi->Write();

  h_reco_match_dphi->Write();
  h_truth_match_dphi->Write();

  h_sim_match_ddphi_counts->Write();
  h2_sim_match_ddphi_counts->Write();
  h_data_pt1pt2dphi_counts->Write();
  h_reco_match_pt1pt2dphi_counts->Write();
  h_reco_match_pt1pt2dphi_countstruth->Write();
  h_truth_match_pt1pt2dphi_counts->Write();

  h_reco_pt1pt2dphi_counts->Write();
  h_truth_pt1pt2dphi_counts->Write();
  h_reco_dphi_counts->Write();
  h_truth_dphi_counts->Write();

  h_reco_match_dphi_counts->Write();
  h_truth_match_dphi_counts->Write();

  for (int i = 0; i < 3; i++) {
    h_correlated_counts_leading_reco[i]->Write();
    h_correlated_counts_subleading_reco[i]->Write();
    h_correlated_counts_leading_truth[i]->Write();
    h_correlated_counts_subleading_truth[i]->Write();

    h_all_counts_leading_reco[i]->Write();
    h_all_counts_subleading_reco[i]->Write();
    h_all_counts_leading_truth[i]->Write();
    h_all_counts_subleading_truth[i]->Write();

  }
  std::cout << h_reco_dphi->GetMean() <<std::endl;  
  fout->Close();
}

int main(int argc, char *argv[])
{

  std::string config = "binning.config";
  int cone_size = 4;
  int set = 0;
  for (int i = 1; i < argc; ++i)
    {
      std::string arg = argv[i];

      if (arg == "-c" && i + 1 < argc)
	{
	  set++;
	  config = argv[++i];  // Next argument as string
	}
      else if (arg == "-v" && i + 1 < argc)
	{
	  verbosity = std::stoi(argv[++i]);  // Next argument as string
	}
      else if (arg == "-r" && i + 1 < argc)
	{
	  set++;
	  cone_size = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else
	{
	  std::cerr << "Unknown or incomplete argument: " << arg << "\n";
	  return 1;
	}
    }
  if (set < 2)
    {
      std::cout << "Not enough settings: " << std::endl;
      std::cout << "[usage] : " << std::endl;
      std::cout << "    ./makeRawDeltaPhi -c binning.config -r 4 " << std::endl;
      return 1;
    }

  makeRawDeltaPhi(cone_size, config);
  return 0;
  
}
