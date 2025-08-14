#include "read_binning.h"
#include "dlUtility.h"
const bool NUCLEAR = false;//true;
const bool useJES=true;
void makeJER(const int cone_size = 4)
{

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  read_binning rb("binning.config");
  
  TFile *corrFile = new TFile(Form("%s/rootfiles/JES_Calib_Default.root", rb.get_code_location().c_str()), "READ");
  TF1 *fJES = (TF1*)corrFile->Get(Form("f_corr_run21_r0%d_z0_eta0123", cone_size));

  /* std::string data_file = "../unfold/TNTUPLE_DIJET_v6_2_ana462_2024p01015
     _v001_gl10-00047289-00048291.root"; */
  /* std::string j10_file = "../unfold/TREE_MATCH_v5_10_new_ProdA_2024-00000021.root"; */
  /* std::string j20_file = "../unfold/TREE_MATCH_v5_20_new_ProdA_2024-00000021.root"; */
  /* std::string j30_file = "../unfold/TREE_MATCH_v5_30_new_ProdA_2024-00000021.root"; */

  std::string data_file = rb.get_tntuple_location() + "/TNTUPLE_DIJET_r0" + std::to_string(cone_size) + "_v6_6_ana468_2024p012_v001_gl10-all.root";

  std::string j00_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v5_0_new_ProdA_2024-00000021.root";
  std::string j10_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_10_new_ProdA_2024-00000021.root";
  std::string j15_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_15_new_ProdA_2024-00000021.root";
  std::string j20_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_20_new_ProdA_2024-00000021.root";
  std::string j30_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_30_new_ProdA_2024-00000021.root";
  std::string j50_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_50_new_ProdA_2024-00000021.root";

  float maxpttruth[6];
  float pt1_truth[6];
  float pt2_truth[6];
  float dphi_truth[6];
  float pt1_reco[6];
  float pt2_reco[6];
  float mbd_vertex[6];
  float em1_reco[6];
  float em2_reco[6];
  float sim_third_jet_pt[6];
  float truth_third_jet_pt[6];
  float nrecojets[6];
  float dphi_reco[6];
  float match[6];

  float mbd_vertex_data;
  
  float pt1_data;
  float pt2_data;
  float em1_data;
  float em2_data;
  float dphi_data;
  float data_third_jet_pt;
  float ndatajets;
  TFile *find = new TFile(data_file.c_str(), "r");
  TNtuple *tnd  = (TNtuple*) find->Get("tn_dijet");;
  if (!tnd)
    {
      std::cout << " no data "<< std::endl;
    }
  tnd->SetBranchAddress("pt1_reco", &pt1_data);
  tnd->SetBranchAddress("mbd_vertex", &mbd_vertex_data);
  tnd->SetBranchAddress("pt2_reco", &pt2_data);
  tnd->SetBranchAddress("em1_reco", &em1_data);
  tnd->SetBranchAddress("em2_reco", &em2_data);
  tnd->SetBranchAddress("dphi_reco", &dphi_data);
  tnd->SetBranchAddress("njets", &ndatajets);
  tnd->SetBranchAddress("thirdjetpt", &data_third_jet_pt);

  TFile *fin[6];
  fin[0] = new TFile(j00_file.c_str(), "r");
  fin[1] = new TFile(j10_file.c_str(), "r");
  fin[2] = new TFile(j15_file.c_str(), "r");
  fin[3] = new TFile(j20_file.c_str(), "r");
  fin[4] = new TFile(j30_file.c_str(), "r");
  fin[5] = new TFile(j50_file.c_str(), "r");

  TNtuple *tn[6];
  for (int i = 0 ; i < 6; i++)
    {
      tn[i] = (TNtuple*) fin[i]->Get("tn_match");
      tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
      tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
      tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
      tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
      tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
      tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
      tn[i]->SetBranchAddress("em1_reco", &em1_reco[i]);
      tn[i]->SetBranchAddress("em2_reco", &em2_reco[i]);
      tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
      tn[i]->SetBranchAddress("nrecojets", &nrecojets[i]);
      tn[i]->SetBranchAddress("thirdjet_pt", &sim_third_jet_pt[i]);
      tn[i]->SetBranchAddress("thirdjet_pt_truth", &truth_third_jet_pt[i]);
      tn[i]->SetBranchAddress("matched", &match[i]);
      tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
    }

  float cs_00 = 4.197e-2;
  float cs_10 = (3.997e-6);
  float cs_15 = (4.073e-7);
  float cs_20 = 6.218e-8;
  float cs_30 = (2.502e-9);
  float cs_50 = (7.2695e-12);
  float scale_factor[6];
  float nen10 = 10000000.;
  float nen15 = 10000000.;
  float nen20 = 10000000.;
  float nen30 = 9990000.;
  float nen50 = 10000000.;
  scale_factor[0] = (nen50/nen10) * cs_00/(cs_50*10);
  scale_factor[1] = (nen50/nen10) * cs_10/cs_50;
  scale_factor[2] = (nen50/nen15) * cs_15/cs_50;
  scale_factor[3] = (nen50/nen20) * cs_20/cs_50;
  scale_factor[4] = (nen50/nen30) * cs_30/cs_50; 
  scale_factor[5] = 1;

  
  float reco_cut = 10.;
  float third_jet_cut_single = 7;
  float dphicut = 3.*TMath::Pi()/4.;

  float third_jet_cuts[8] = {10, 9, 8, 7, 6, 5, 4, 3}; 
  float sample_boundary[7] = {0,0,0,0,0, 0, 0};

  switch (cone_size)
    {
      case 2:
	
	sample_boundary[0] = 0;
	sample_boundary[1] = 12;
	sample_boundary[2] = 15;
	sample_boundary[3] = 20;
	sample_boundary[4] = 31;
	sample_boundary[5] = 50;
	sample_boundary[6] = 100;
	
	break;
	
      case 3:
	sample_boundary[0] = 0;
	sample_boundary[1] = 13;
	sample_boundary[2] = 16;
	sample_boundary[3] = 21;
	sample_boundary[4] = 33;
	sample_boundary[5] = 51;
	sample_boundary[6] = 100;
	break;

      case 4:
	sample_boundary[0] = 0;
	sample_boundary[1] = 14;
	sample_boundary[2] = 17;
	sample_boundary[3] = 22;
	sample_boundary[4] = 35;
	sample_boundary[5] = 52;
	sample_boundary[6] = 100;
	break;

      case 5:
	sample_boundary[0] = 0;
	sample_boundary[1] = 15;
	sample_boundary[2] = 24;
	sample_boundary[3] = 30;
	sample_boundary[4] = 40;
	sample_boundary[5] = 60;
	sample_boundary[6] = 100;
	break;

      case 6:
	sample_boundary[0] = 0;
	sample_boundary[1] = 17;
	sample_boundary[2] = 26;
	sample_boundary[3] = 35;
	sample_boundary[4] = 45;
	sample_boundary[5] = 63;
	sample_boundary[6] = 100;
	break;
    } 
  TRandom *rng = new TRandom();  

  TH1D *h_lead_sample[6];
  TH1D *h_lead_combined = new TH1D("h_lead_combined","; Leading Truth Jet p_{T} [GeV]; lumi scale * counts",100, 0, 100);
  for (int i = 0; i < 6; i++)
    {
      h_lead_sample[i] = new TH1D(Form("h_lead_sample_%d", i), "; Leading Truth Jet p_{T} [GeV]; lumi scale * counts",100, 0, 100);
    }
  TH1D *h_sigp_sim = new TH1D("h_sigp_sim",";<p_{T}> [GeV];#sigma(p_{T, #psi})", 4, 20, 40);
  TH1D *h_sigp_data = new TH1D("h_sigp_data",";<p_{T}> [GeV];#sigma(p_{T, #psi})", 4, 20, 40);

  TH1D *h_sigv_sim = new TH1D("h_sigv_sim",";<p_{T}> [GeV];#sigma(p_{T, #eta})", 4, 20, 40);
  TH1D *h_sigv_data = new TH1D("h_sigv_data",";<p_{T}> [GeV];#sigma(p_{T, #eta})", 4, 20, 40);

  TH1D *h_jer_sim = new TH1D("h_jer_sim",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
  TH1D *h_jer_data = new TH1D("h_jer_data",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);

  TH1D *h_ptavg_sim = new TH1D("h_ptavg_sim",";<p_{T}> [GeV];Counts", 4, 20, 40);
  TH1D *h_ptavg_data = new TH1D("h_ptavg_data",";<p_{T}> [GeV];Counts", 4, 20, 40);

  TH2D *h_trigger_ptv_pta_sim = new TH2D("h_trigger_ptv_pta_sim",";pta;ptv", 4, 20, 40, 80, -40, 40);
  TH2D *h_trigger_ptv_pta_data = new TH2D("h_trigger_ptv_pta_data",";pta;ptv", 4, 20, 40, 80, -40, 40);

  TH2D *h_trigger_ptp_pta_sim = new TH2D("h_trigger_ptp_pta_sim",";pta;ptp", 4, 20, 40, 80, -40, 40);
  TH2D *h_trigger_ptp_pta_data = new TH2D("h_trigger_ptp_pta_data",";pta;ptp", 4, 20, 40, 80, -40, 40);

  
  TProfile *hp_trigger_ptv_pta_sim = new TProfile("hp_trigger_ptv_pta_sim",";pta;ptv", 4, 20, 40);//,"s");
  TProfile *hp_trigger_ptv_pta_data = new TProfile("hp_trigger_ptv_pta_data",";pta;ptv", 4, 20, 40);//,"s");

  TProfile *hp_trigger_ptp_pta_sim = new TProfile("hp_trigger_ptp_pta_sim",";pta;ptp", 4, 20, 40);//,"s");
  TProfile *hp_trigger_ptp_pta_data = new TProfile("hp_trigger_ptp_pta_data",";pta;ptp", 4, 20, 40);//,"s");

  TProfile *hp_cosdphi_pta_sim = new TProfile("hp_cosdphi_pta_sim",";pta;ptp", 4, 20, 40,"s");
  TProfile *hp_cosdphi_pta_data = new TProfile("hp_cosdphi_pta_data",";pta;ptp", 4, 20, 40,"s");


  // 3jet cut
  TH1D *h_sigp_sim_3jet[8]; // new TH1D("h_sigp_sim",";<p_{T}> [GeV];#sigma(p_{T, #psi})", 10, 20, 40);
  TH1D *h_sigp_data_3jet[8]; // new TH1D("h_sigp_data",";<p_{T}> [GeV];#sigma(p_{T, #psi})", 10, 20, 40);

  TH1D *h_sigv_sim_3jet[8]; // new TH1D("h_sigv_sim",";<p_{T}> [GeV];#sigma(p_{T, #eta})", 10, 20, 40);
  TH1D *h_sigv_data_3jet[8]; // new TH1D("h_sigv_data",";<p_{T}> [GeV];#sigma(p_{T, #eta})", 10, 20, 40);

  TH1D *h_jer_sim_3jet[8]; // new TH1D("h_jer_sim",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_data_3jet[8]; // new TH1D("h_jer_data",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);

  TH1D *h_jer_truth_im_3jet[8]; // new TH1D("h_jer_truth",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_sim_im_3jet[8]; // new TH1D("h_jer_sim",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_data_im_3jet[8]; // new TH1D("h_jer_data",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);

  TH1D *h_jer_truth_im_softcorr_3jet[8]; // new TH1D("h_jer_truth",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_sim_im_softcorr_3jet[8]; // new TH1D("h_jer_sim",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_data_im_softcorr_3jet[8]; // new TH1D("h_jer_data",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);

  TH1D *h_jer_sim_im_softcorr_part_3jet[8]; // new TH1D("h_jer_sim",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);
  TH1D *h_jer_data_im_softcorr_part_3jet[8]; // new TH1D("h_jer_data",";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 10, 20, 40);

  TH1D *h_ptavg_sim_3jet[8]; // new TH1D("h_ptavg_sim",";<p_{T}> [GeV];Counts", 10, 20, 40);
  TH1D *h_ptavg_data_3jet[8]; // new TH1D("h_ptavg_data",";<p_{T}> [GeV];Counts", 10, 20, 40);

  TH2D *h_trigger_ptv_pta_sim_3jet[8]; // new TH2D("h_trigger_ptv_pta_sim",";pta;ptv", 10, 20, 40, 40, -20, 20);
  TH2D *h_trigger_ptv_pta_data_3jet[8]; // new TH2D("h_trigger_ptv_pta_data",";pta;ptv", 10, 20, 40, 40, -20, 20);

  TH2D *h_trigger_ptp_pta_sim_3jet[8]; // new TH2D("h_trigger_ptp_pta_sim",";pta;ptp", 10, 20, 40, 40, -20, 20);
  TH2D *h_trigger_ptp_pta_data_3jet[8]; // new TH2D("h_trigger_ptp_pta_data",";pta;ptp", 10, 20, 40, 40, -20, 20);

  
  TProfile *hp_trigger_ptv_pta_sim_3jet[8]; // new TProfile("hp_trigger_ptv_pta_sim",";pta;ptv", 10, 20, 40);//,"s");
  TProfile *hp_trigger_ptv_pta_data_3jet[8]; // new TProfile("hp_trigger_ptv_pta_data",";pta;ptv", 10, 20, 40);//,"s");

  TProfile *hp_trigger_ptp_pta_sim_3jet[8]; // new TProfile("hp_trigger_ptp_pta_sim",";pta;ptp", 10, 20, 40);//,"s");
  TProfile *hp_trigger_ptp_pta_data_3jet[8]; // new TProfile("hp_trigger_ptp_pta_data",";pta;ptp", 10, 20, 40);//,"s");

  TProfile *hp_cosdphi_pta_sim_3jet[8]; // new TProfile("hp_cosdphi_pta_sim",";pta;ptp", 10, 20, 40);//,"s");
  TProfile *hp_cosdphi_pta_data_3jet[8]; // new TProfile("hp_cosdphi_pta_data",";pta;ptp", 10, 20, 40);//,"s");

  for (int i = 0; i < 8; i++)
  {
    h_sigp_sim_3jet[i] = new TH1D(Form("h_sigp_sim_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T, #psi})", 4, 20, 40);
    h_sigp_data_3jet[i] = new TH1D(Form("h_sigp_data_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T, #psi})", 4, 20, 40);

    h_sigv_sim_3jet[i] = new TH1D(Form("h_sigv_sim_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T, #eta})", 4, 20, 40);
    h_sigv_data_3jet[i] = new TH1D(Form("h_sigv_data_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T, #eta})", 4, 20, 40);
    h_jer_sim_3jet[i] = new TH1D(Form("h_jer_sim_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_data_3jet[i] = new TH1D(Form("h_jer_data_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);

    h_jer_sim_im_3jet[i] = new TH1D(Form("h_jer_sim_im_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_data_im_3jet[i] = new TH1D(Form("h_jer_data_im_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_truth_im_3jet[i] = new TH1D(Form("h_jer_truth_im_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);

    h_jer_sim_im_softcorr_3jet[i] = new TH1D(Form("h_jer_sim_im_softcorr_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_data_im_softcorr_3jet[i] = new TH1D(Form("h_jer_data_im_softcorr_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_truth_im_softcorr_3jet[i] = new TH1D(Form("h_jer_truth_im_softcorr_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);

    h_jer_sim_im_softcorr_part_3jet[i] = new TH1D(Form("h_jer_sim_im_softcorr_part_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);
    h_jer_data_im_softcorr_part_3jet[i] = new TH1D(Form("h_jer_data_im_softcorr_part_3jet%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 4, 20, 40);

    h_ptavg_sim_3jet[i] = new TH1D(Form("h_ptavg_sim_3jet%d", i),";<p_{T}> [GeV];Counts", 4, 20, 40);
    h_ptavg_data_3jet[i] = new TH1D(Form("h_ptavg_data_3jet%d", i),";<p_{T}> [GeV];Counts", 4, 20, 40);

    h_trigger_ptv_pta_sim_3jet[i] = new TH2D(Form("h_trigger_ptv_pta_sim_3jet%d", i),";pta;ptv", 4, 20, 40, 80, -40, 40);
    h_trigger_ptv_pta_data_3jet[i] = new TH2D(Form("h_trigger_ptv_pta_data_3jet%d", i),";pta;ptv", 4, 20, 40, 80, -40, 40);
    h_trigger_ptp_pta_sim_3jet[i] = new TH2D(Form("h_trigger_ptp_pta_sim_3jet%d", i),";pta;ptp", 4, 20, 40, 80, -40, 40);
    h_trigger_ptp_pta_data_3jet[i] = new TH2D(Form("h_trigger_ptp_pta_data_3jet%d", i),";pta;ptp", 4, 20, 40, 80, -40, 40);
    hp_trigger_ptv_pta_sim_3jet[i] = new TProfile(Form("hp_trigger_ptv_pta_sim_3jet%d", i),";pta;ptv", 4, 20, 40);//,"s");
    hp_trigger_ptv_pta_data_3jet[i] = new TProfile(Form("hp_trigger_ptv_pta_data_3jet%d", i),";pta;ptv", 4, 20, 40);//,"s");
    hp_trigger_ptp_pta_sim_3jet[i] = new TProfile(Form("hp_trigger_ptp_pta_sim_3jet%d", i),";pta;ptp", 4, 20, 40);//,"s");
    hp_trigger_ptp_pta_data_3jet[i] = new TProfile(Form("hp_trigger_ptp_pta_data_3jet%d", i),";pta;ptp", 4, 20, 40);//,"s");
    hp_cosdphi_pta_sim_3jet[i] = new TProfile(Form("hp_cosdphi_pta_sim_3jet%d", i),";pta;ptp", 4, 20, 40,"s");
    hp_cosdphi_pta_data_3jet[i] = new TProfile(Form("hp_cosdphi_pta_data_3jet%d", i),";pta;ptp", 4, 20, 40,"s");
  }
  // dijet balance

  TH2D *h2_aj_avgpt_sim_3jet[8];
  TProfile *hp_aj_avgpt_sim_3jet[8];
  TH1D *h_avgpt_sim_3jet[8];

  TH2D *h2_aj_avgpt_truth_3jet[8];
  TProfile *hp_aj_avgpt_truth_3jet[8];
  TH1D *h_avgpt_truth_3jet[8];

  TH2D *h2_aj_avgpt_data_3jet[8];
  TProfile *hp_aj_avgpt_data_3jet[8];
  TH1D *h_avgpt_data_3jet[8];
  
  for (int i = 0; i < 8; i++)
    {
      h2_aj_avgpt_sim_3jet[i] = new TH2D(Form("h2_aj_avgpt_sim_3jet%d", i), ";avg pt; aj",4, 20, 40, 100, -1, 1);
      hp_aj_avgpt_sim_3jet[i] = new TProfile(Form("hp_aj_avgpt_sim_3jet%d", i), ";avg pt; aj",4, 20, 40, "s");
      h_avgpt_sim_3jet[i] = new TH1D(Form("h_avgpt_sim_3jet%d", i), ";avg pt; aj",4, 20, 40); 

      h2_aj_avgpt_truth_3jet[i] = new TH2D(Form("h2_aj_avgpt_truth_3jet%d", i), ";avg pt; aj",4, 20, 40, 100, -1, 1);
      hp_aj_avgpt_truth_3jet[i] = new TProfile(Form("hp_aj_avgpt_truth_3jet%d", i), ";avg pt; aj",4, 20, 40, "s");
      h_avgpt_truth_3jet[i] = new TH1D(Form("h_avgpt_truth_3jet%d", i), ";avg pt; aj",4, 20, 40); 

      h2_aj_avgpt_data_3jet[i] = new TH2D(Form("h2_aj_avgpt_data_3jet%d", i), ";avg pt; aj",4, 20, 40, 100, -1, 1);
      hp_aj_avgpt_data_3jet[i] = new TProfile(Form("hp_aj_avgpt_data_3jet%d", i), ";avg pt; aj",4, 20, 40, "s");
      h_avgpt_data_3jet[i] = new TH1D(Form("h_avgpt_data_3jet%d", i), ";avg pt; aj",4, 20, 40); 
    }
  
  // EM fraction

  const int nbins_em = 2;
  float em_bounds[3] = {-0.2, 0.5, 1.2};
  TH1D *h_em_sigp_sim[5];
  TH1D *h_em_sigp_data[5];

  TH1D *h_em_sigv_sim[5];
  TH1D *h_em_sigv_data[5];

  TH1D *h_em_jer_sim[5];
  TH1D *h_em_jer_data[5];

  TH1D *h_em_ptavg_sim[5];
  TH1D *h_em_ptavg_data[5];

  TH2D *h_em_trigger_ptv_pta_sim[5];
  TH2D *h_em_trigger_ptv_pta_data[5];

  TH2D *h_em_trigger_ptp_pta_sim[5];
  TH2D *h_em_trigger_ptp_pta_data[5];

  
  TProfile *hp_em_trigger_ptv_pta_sim[5];
  TProfile *hp_em_trigger_ptv_pta_data[5];

  TProfile *hp_em_trigger_ptp_pta_sim[5];
  TProfile *hp_em_trigger_ptp_pta_data[5];

  TProfile *hp_em_cosdphi_pta_sim[5];
  TProfile *hp_em_cosdphi_pta_data[5];
  for (int i = 0; i < 5; i++)
    {
      h_em_sigp_sim[i] = new TH1D(Form("h_em_sigp_sim_%d", i),";<p_{T}> [GeV];#sigma(p_{T, #psi})", 5, 20, 40);
      h_em_sigp_data[i] = new TH1D(Form("h_em_sigp_data_%d", i),";<p_{T}> [GeV];#sigma(p_{T, #psi})", 5, 20, 40);

      h_em_sigv_sim[i] = new TH1D(Form("h_em_sigv_sim_%d", i),";<p_{T}> [GeV];#sigma(p_{T, #eta})", 5, 20, 40);
      h_em_sigv_data[i] = new TH1D(Form("h_em_sigv_data_%d", i),";<p_{T}> [GeV];#sigma(p_{T, #eta})", 5, 20, 40);

      h_em_jer_sim[i] = new TH1D(Form("h_em_jer_sim_%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 5, 20, 40);
      h_em_jer_data[i] = new TH1D(Form("h_em_jer_data_%d", i),";<p_{T}> [GeV];#sigma(p_{T})/p_{T}", 5, 20, 40);

      h_em_ptavg_sim[i] = new TH1D(Form("h_em_ptavg_sim_%d", i),";<p_{T}> [GeV];Counts", 5, 20, 40);
      h_em_ptavg_data[i] = new TH1D(Form("h_em_ptavg_data_%d", i),";<p_{T}> [GeV];Counts", 5, 20, 40);

      h_em_trigger_ptv_pta_sim[i] = new TH2D(Form("h_em_trigger_ptv_pta_sim_%d", i),";pta;ptv", 5, 20, 40, 80, -40, 40);
      h_em_trigger_ptv_pta_data[i] = new TH2D(Form("h_em_trigger_ptv_pta_data_%d", i),";pta;ptv", 5, 20, 40, 80, -40, 40);

      h_em_trigger_ptp_pta_sim[i] = new TH2D(Form("h_em_trigger_ptp_pta_sim_%d", i),";pta;ptp", 5, 20, 40, 80, -40, 40);
      h_em_trigger_ptp_pta_data[i] = new TH2D(Form("h_em_trigger_ptp_pta_data_%d", i),";pta;ptp", 5, 20, 40, 80, -40, 40);

  
      hp_em_trigger_ptv_pta_sim[i] = new TProfile(Form("hp_em_trigger_ptv_pta_sim_%d", i),";pta;ptv", 5, 20, 40);//,"s");
      hp_em_trigger_ptv_pta_data[i] = new TProfile(Form("hp_em_trigger_ptv_pta_data_%d", i),";pta;ptv", 5, 20, 40);//,"s");

      hp_em_trigger_ptp_pta_sim[i] = new TProfile(Form("hp_em_trigger_ptp_pta_sim_%d", i),";pta;ptp", 5, 20, 40);//,"s");
      hp_em_trigger_ptp_pta_data[i] = new TProfile(Form("hp_em_trigger_ptp_pta_data_%d", i),";pta;ptp", 5, 20, 40);//,"s");

      hp_em_cosdphi_pta_sim[i] = new TProfile(Form("hp_em_cosdphi_pta_sim_%d", i),";pta;ptp", 5, 20, 40,"s");
      hp_em_cosdphi_pta_data[i] = new TProfile(Form("hp_em_cosdphi_pta_data_%d", i),";pta;ptp", 5, 20, 40,"s");

    }
  
  int entriesd = tnd->GetEntries();

  for (int i = 0; i < entriesd; i++)
    {
      tnd->GetEntry(i);
      //if (fabs(mbd_vertex_data) > 30) continue;
      float trigger_jet = fJES->Eval(pt1_data);
      float probe_jet = fJES->Eval(pt2_data);
      float third_jet = fJES->Eval(data_third_jet_pt);
      if (!useJES)
	{
	  trigger_jet = pt1_data;
	  probe_jet = pt2_data;
	  third_jet = data_third_jet_pt;

	}
      if (rng->Uniform() > 0.5)
	{
	  float tempy = trigger_jet;
	  trigger_jet = probe_jet;
	  probe_jet = tempy;
	}


      bool data_good = (trigger_jet >= reco_cut && probe_jet >= reco_cut && dphi_data >= dphicut && ndatajets == 2);      
      
      float ptv = trigger_jet*cos(dphi_data/2) + probe_jet*cos(dphi_data/2);
      float ptp = trigger_jet*sin(dphi_data/2) - probe_jet*sin(dphi_data/2);

      float pta = (trigger_jet + probe_jet)/2.;
      float aj = (trigger_jet - probe_jet) / (trigger_jet + probe_jet);
      if (rng->Uniform() > 0.5)
	{
	  ptv = -1*ptv;
	}

      if (trigger_jet >= reco_cut && probe_jet >= reco_cut && dphi_data >= dphicut)
	{

	  for (int i3 = 0 ; i3 < 8; i3++)
	    {
	      if (third_jet < third_jet_cuts[i3])
		{
		  h2_aj_avgpt_data_3jet[i3]->Fill(pta, aj);
		  hp_aj_avgpt_data_3jet[i3]->Fill(pta, aj);
		  h_avgpt_data_3jet[i3]->Fill(pta);

		  h_ptavg_data_3jet[i3]->Fill(pta);
		  h_trigger_ptv_pta_data_3jet[i3]->Fill(pta, ptv);
		  hp_trigger_ptv_pta_data_3jet[i3]->Fill(pta, ptv);
		  h_trigger_ptp_pta_data_3jet[i3]->Fill(pta, ptp);
		  hp_trigger_ptp_pta_data_3jet[i3]->Fill(pta, ptp);      
		  hp_cosdphi_pta_data_3jet[i3]->Fill(pta, fabs(cos(dphi_data)));

		}
	    }
	}

      if (!data_good) continue;

      h_ptavg_data->Fill(pta);
      h_trigger_ptv_pta_data->Fill(pta, ptv);
      hp_trigger_ptv_pta_data->Fill(pta, ptv);
      h_trigger_ptp_pta_data->Fill(pta, ptp);
      hp_trigger_ptp_pta_data->Fill(pta, ptp);      
      hp_cosdphi_pta_data->Fill(pta, fabs(cos(dphi_data)));

      
      float ema = (em1_data*pt1_data + em2_data*pt2_data)/(pt1_data+pt2_data);
      for (int iem = 0; iem < nbins_em; iem++)
	{
	  if (ema >= em_bounds[iem] && ema < em_bounds[iem+1])
	  {
	    h_em_ptavg_data[iem]->Fill(pta);
	    h_em_trigger_ptv_pta_data[iem]->Fill(pta, ptv);
	    hp_em_trigger_ptv_pta_data[iem]->Fill(pta, ptv);
	    h_em_trigger_ptp_pta_data[iem]->Fill(pta, ptp);
	    hp_em_trigger_ptp_pta_data[iem]->Fill(pta, ptp);      
	    hp_em_cosdphi_pta_data[iem]->Fill(pta, fabs(cos(dphi_data)));	    
	    break;
	  }
	}
    }

  // Sim
  for (int isample = 1; isample < 6; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries2 = tn[isample]->GetEntries();
      for (int i = 0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  //if (fabs(mbd_vertex[isample]) > 30) continue;
	  h_lead_sample[isample]->Fill(maxpttruth[isample], scale_factor[isample]);
	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;
	  h_lead_combined->Fill(maxpttruth[isample], scale_factor[isample]);
	  float trigger_jet = fJES->Eval(pt1_reco[isample]);
	  float probe_jet = fJES->Eval(pt2_reco[isample]);
	  float third_jet = fJES->Eval(sim_third_jet_pt[isample]);

	  if (!useJES)
	    {
	      trigger_jet = pt1_reco[isample];
	      probe_jet = pt2_reco[isample];
	      third_jet = sim_third_jet_pt[isample];
	    }
	  if (rng->Uniform() > 0.5)
	    {
	      float tempy = trigger_jet;
	      trigger_jet = probe_jet;
	      probe_jet = tempy;
	    }

	  bool reco_good = (trigger_jet >= reco_cut && probe_jet >= reco_cut && dphi_reco[isample] >= dphicut and third_jet < 3);

	  float ptv = trigger_jet*cos(dphi_reco[isample]/2) + probe_jet*cos(dphi_reco[isample]/2);
	  float ptp = trigger_jet*sin(dphi_reco[isample]/2) - probe_jet*sin(dphi_reco[isample]/2);
	  float pta = (trigger_jet + probe_jet)/2.;
	  float aj = (trigger_jet - probe_jet) / (trigger_jet + probe_jet);

	  if (rng->Uniform() > 0.5)
	    {
	      ptv = -1*ptv;
	    }

	  if (trigger_jet >= reco_cut && probe_jet >= reco_cut && dphi_reco[isample] >= dphicut)// &&  match[isample])
	    {
	      for (int i3 = 0 ; i3 < 8; i3++)
		{
		  if (third_jet < third_jet_cuts[i3])
		    {
		      h2_aj_avgpt_sim_3jet[i3]->Fill(pta, aj, scale_factor[isample]);
		      hp_aj_avgpt_sim_3jet[i3]->Fill(pta, aj, scale_factor[isample]);
		      h_avgpt_sim_3jet[i3]->Fill(pta, scale_factor[isample]);

		      h_ptavg_sim_3jet[i3]->Fill(pta, scale_factor[isample]);
		      h_trigger_ptv_pta_sim_3jet[i3]->Fill(pta, ptv, scale_factor[isample]);
		      hp_trigger_ptv_pta_sim_3jet[i3]->Fill(pta, ptv, scale_factor[isample]);
		      h_trigger_ptp_pta_sim_3jet[i3]->Fill(pta, ptp, scale_factor[isample]);
		      hp_trigger_ptp_pta_sim_3jet[i3]->Fill(pta, ptp, scale_factor[isample]);
		      hp_cosdphi_pta_sim_3jet[i3]->Fill(pta, fabs(cos(dphi_reco[isample])), scale_factor[isample]);
		    }
		}
	    }

	  float triggertruth_jet = pt1_truth[isample];
	  float probetruth_jet = pt2_truth[isample];

	  if (rng->Uniform() > 0.5)
	    {
	      triggertruth_jet = pt2_truth[isample];
	      probetruth_jet = pt1_truth[isample];
	    }


	  float ptat = (triggertruth_jet + probetruth_jet)/2.;
	  float ajt = (triggertruth_jet - probetruth_jet) / (triggertruth_jet + probetruth_jet);

	  if (triggertruth_jet >= reco_cut && probetruth_jet >= reco_cut && dphi_truth[isample] >= dphicut &&  match[isample])
	    {

	      for (int i3 = 0 ; i3 < 8; i3++)
		{
		  if (truth_third_jet_pt[isample] < third_jet_cuts[i3])
		    {
		      h2_aj_avgpt_truth_3jet[i3]->Fill(ptat, ajt, scale_factor[isample]);
		      hp_aj_avgpt_truth_3jet[i3]->Fill(ptat, ajt, scale_factor[isample]);
		      h_avgpt_truth_3jet[i3]->Fill(ptat, scale_factor[isample]);
		    }
		}
	    }

	  
	  if (!reco_good) continue;

	  h_ptavg_sim->Fill(pta, scale_factor[isample]);
	  h_trigger_ptv_pta_sim->Fill(pta, ptv, scale_factor[isample]);
	  hp_trigger_ptv_pta_sim->Fill(pta, ptv, scale_factor[isample]);
	  h_trigger_ptp_pta_sim->Fill(pta, ptp, scale_factor[isample]);
	  hp_trigger_ptp_pta_sim->Fill(pta, ptp, scale_factor[isample]);
	  hp_cosdphi_pta_sim->Fill(pta, fabs(cos(dphi_reco[isample])), scale_factor[isample]);

	  float ema = (em1_reco[isample]*pt1_reco[isample] + em2_reco[isample]*pt2_reco[isample])/(pt1_reco[isample]+pt2_reco[isample]);
	  for (int iem = 0; iem < 5; iem++)
	    {
	      if (ema >= em_bounds[iem] && ema < em_bounds[iem+1])
	      {
		h_em_ptavg_sim[iem]->Fill(pta, scale_factor[isample]);
		h_em_trigger_ptv_pta_sim[iem]->Fill(pta, ptv, scale_factor[isample]);
		hp_em_trigger_ptv_pta_sim[iem]->Fill(pta, ptv, scale_factor[isample]);
		h_em_trigger_ptp_pta_sim[iem]->Fill(pta, ptp, scale_factor[isample]);
		hp_em_trigger_ptp_pta_sim[iem]->Fill(pta, ptp, scale_factor[isample]);      
		hp_em_cosdphi_pta_sim[iem]->Fill(pta, fabs(cos(dphi_reco[isample])), scale_factor[isample]);	    
		break;
	      }
	    }

	}
    }


  TH1D *h_sim_3jet_width[4];
  TH1D *h_data_3jet_width[4];
  TH1D *h_truth_3jet_width[4];

  TH1D *h_sim_3jet_softcorr[8];
  TH1D *h_data_3jet_softcorr[8];
  TH1D *h_truth_3jet_softcorr[8];

  for (int i = 0; i < 8; i++)
    {
      h_truth_3jet_softcorr[i] = new TH1D(Form("h_truth_3jet_softcorr_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 4, 20, 40);
      h_sim_3jet_softcorr[i] = new TH1D(Form("h_sim_3jet_softcorr_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 4, 20, 40);
      h_data_3jet_softcorr[i] = new TH1D(Form("h_data_3jet_softcorr_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 4, 20, 40);
    }
  
  TH1D *h_sim_3jet_sigv[4];
  TH1D *h_data_3jet_sigv[4];

  TH1D *h_sim_3jet_sigp[4];
  TH1D *h_data_3jet_sigp[4];

  TH1D *h_sim_3jet_sigma[4];
  TH1D *h_data_3jet_sigma[4];

  TH1D *h_sim_3jet_bi_width[4];
  TH1D *h_data_3jet_bi_width[4];

  TF1 *fg = new TF1("fg","gaus", -20, 20);

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < 4; i++)
    {
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      h_sim_3jet_sigma[i] = new TH1D(Form("h_sim_3jet_sigma_%d", i), ";p_{T,3}; #sigma(p_{T, #psi})", 10, 0.5, 10.5);
      h_data_3jet_sigma[i] = new TH1D(Form("h_data_3jet_sigma_%d", i), ";p_{T,3}; #sigma(p_{T, #psi})", 10, 0.5, 10.5);
      
      h_sim_3jet_sigp[i] = new TH1D(Form("h_sim_3jet_sigp_%d", i), ";p_{T,3}; #sigma(p_{T, #psi})", 10, 0.5, 10.5);
      h_data_3jet_sigp[i] = new TH1D(Form("h_data_3jet_sigp_%d", i), ";p_{T,3}; #sigma(p_{T, #psi})", 10, 0.5, 10.5);

      h_sim_3jet_sigv[i] = new TH1D(Form("h_sim_3jet_sigv_%d", i), ";p_{T,3}; #sigma(p_{T, #eta})", 10, 0.5, 10.5);
      h_data_3jet_sigv[i] = new TH1D(Form("h_data_3jet_sigv_%d", i), ";p_{T,3}; #sigma(p_{T, #eta})", 10, 0.5, 10.5);
      
      h_sim_3jet_bi_width[i] = new TH1D(Form("h_sim_3jet_bi_width_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 10, 0.5, 10.5);
      h_data_3jet_bi_width[i] = new TH1D(Form("h_data_3jet_bi_width_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 10, 0.5, 10.5);

      for (int j = 0; j < 8; j++)
	{

	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  
	  TH1D *h_simv = (TH1D*) h_trigger_ptv_pta_sim_3jet[j]->ProjectionY("hv", i+1, i+1);
	  h_simv->Fit("fg", "R");
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  TH1D *h_simp = (TH1D*) h_trigger_ptp_pta_sim_3jet[j]->ProjectionY("hp", i+1, i+1);
	  h_simp->Fit("fg", "R");
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);
	  
	  h_sim_3jet_sigv[i]->SetBinContent(int(third_jet_cuts[j]), sigv);
	  h_sim_3jet_sigv[i]->SetBinError(int(third_jet_cuts[j]), sig2v);

	  h_sim_3jet_sigp[i]->SetBinContent(int(third_jet_cuts[j]), sigp);
	  h_sim_3jet_sigp[i]->SetBinError(int(third_jet_cuts[j]), sig2p);

	  float pta = hp_trigger_ptv_pta_sim_3jet[i]->GetBinCenter(i+1);
	  
	  float cosdphi = hp_cosdphi_pta_sim_3jet[j]->GetBinContent(i+1);
	  float sigcosdphi = hp_cosdphi_pta_sim_3jet[j]->GetBinError(i+1);
	  float quad_sub = sqrt(sigp*sigp - sigv*sigv);

	  float sig2_quad_sub = TMath::Power(sigp*sig2p/quad_sub, 2) + TMath::Power(sigv*sig2v/quad_sub, 2);
	  float sig_quad_sub = sqrt(sig2_quad_sub);
	  float res = quad_sub /(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig_quad_sub/quad_sub, 2) + TMath::Power(sigcosdphi/cosdphi, 2));
	  
	  if (sigp <= sigv) continue;
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  std::cout << "sim 3jet : " << j <<  " bin : " << i << " : " << "sigp ( " << sigp << " +- "<<sig2p<<") sige ( " << sigv << " +- " << sig2v << ") cosdphi  ( " << cosdphi << " +- " << sigcosdphi << " ) res ( " <<  res << " +- " << sigres << " ) " << std::endl;
	  
	  h_sim_3jet_sigma[i]->SetBinContent(int(third_jet_cuts[j]), quad_sub);
	  h_sim_3jet_sigma[i]->SetBinError(int(third_jet_cuts[j]), sig_quad_sub);
	  
	  h_sim_3jet_bi_width[i]->SetBinContent(int(third_jet_cuts[j]), res);
	  h_sim_3jet_bi_width[i]->SetBinError(int(third_jet_cuts[j]), sigres);

	  h_jer_sim_3jet[j]->SetBinContent(i+1, res);
	  h_jer_sim_3jet[j]->SetBinError(i+1, sigres);
	}
      for (int j = 0; j < 8; j++)
	{
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  TH1D *h_datav = (TH1D*) h_trigger_ptv_pta_data_3jet[j]->ProjectionY("hv", i+1, i+1);
	  h_datav->Fit("fg", "R");
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  TH1D *h_datap = (TH1D*) h_trigger_ptp_pta_data_3jet[j]->ProjectionY("hp", i+1, i+1);
	  h_datap->Fit("fg", "R");
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);

	  /* h_trigger_ptv_pta_data_3jet[j]->FitSlicesY(fg, i+1, i+1); */
	  /* float sigv = fg->GetParameter(2); */
	  /* float sig2v = fg->GetParError(2); */
	  /* h_trigger_ptp_pta_data_3jet[j]->FitSlicesY(fg, i+1, i+1); */
	  /* float sigp = fg->GetParameter(2); */
	  /* float sig2p = fg->GetParError(2); */
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  
	  h_data_3jet_sigv[i]->SetBinContent(int(third_jet_cuts[j]), sigv);//j, third_jet_cuts[j], sigv);
	  h_data_3jet_sigv[i]->SetBinError(int(third_jet_cuts[j]), sig2v);

	  h_data_3jet_sigp[i]->SetBinContent(int(third_jet_cuts[j]), sigp);
	  h_data_3jet_sigp[i]->SetBinError(int(third_jet_cuts[j]), sig2p);
 

	  float pta = hp_trigger_ptv_pta_data_3jet[i]->GetBinCenter(i+1);
	  float sigpta = hp_trigger_ptv_pta_data_3jet[i]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_cosdphi_pta_data_3jet[j]->GetBinContent(i+1);
	  float sigcosdphi = hp_cosdphi_pta_data_3jet[j]->GetBinError(i+1);
	  float quad_sub = sqrt(sigp*sigp - sigv*sigv);
	  float sig2_quad_sub = TMath::Power(sigp*sig2p/quad_sub, 2) + TMath::Power(sigv*sig2v/quad_sub, 2);
	  float sig_quad_sub = sqrt(sig2_quad_sub);
	  float res = quad_sub /(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig_quad_sub/quad_sub, 2) + TMath::Power(sigcosdphi/cosdphi, 2));
	  
	  if (sigp <= sigv) continue;
	  std::cout << "data 3jet : " << j <<  " bin : " << i << " : " << "sigp ( " << sigp << " +- "<<sig2p<<") sige ( " << sigv << " +- " << sig2v << ") cosdphi  ( " << cosdphi << " +- " << sigcosdphi << " ) res ( " <<  res << " +- " << sigres << " ) " << std::endl;
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  h_data_3jet_sigma[i]->SetBinContent(int(third_jet_cuts[j]), quad_sub);
	  h_data_3jet_sigma[i]->SetBinError(int(third_jet_cuts[j]), sig_quad_sub);
	  
	  h_data_3jet_bi_width[i]->SetBinContent(int(third_jet_cuts[j]), res);
	  h_data_3jet_bi_width[i]->SetBinError(int(third_jet_cuts[j]), sigres);

	  h_jer_data_3jet[j]->SetBinContent(i+1, res);
	  h_jer_data_3jet[j]->SetBinError(i+1, sigres);
	  
	}

    }

  TH1D *h_aj_data[4][8];
  TH1D *h_aj_sim[4][8];
  TH1D *h_aj_truth[4][8];

  TF1 *fgg = new TF1("fgg", "gaus", -0.3, 0.3);
  fgg->SetParameters(1, 0, 0.2);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < 4; i++)
    {
      

      h_truth_3jet_width[i] = new TH1D(Form("h_truth_3jet_width_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 10, 0.5, 10.5);
      h_sim_3jet_width[i] = new TH1D(Form("h_sim_3jet_width_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 10, 0.5, 10.5);
      h_data_3jet_width[i] = new TH1D(Form("h_sim_3jet_width_%d", i), ";p_{T,3}; #sigma(p_{T})/p_{T}", 10, 0.5, 10.5);

      
      for (int j = 0; j < 8; j++)
	{
	  h_aj_data[i][j] = (TH1D*) h2_aj_avgpt_data_3jet[j]->ProjectionY(Form("h_aj_data_%d_%d", i, j), i+1, i+1);
	  h_aj_sim[i][j] = (TH1D*) h2_aj_avgpt_sim_3jet[j]->ProjectionY(Form("h_aj_sim_%d_%d", i, j), i+1, i+1);
	  h_aj_truth[i][j] = (TH1D*) h2_aj_avgpt_truth_3jet[j]->ProjectionY(Form("h_aj_truth_%d_%d", i, j), i+1, i+1);

	  h_aj_data[i][j]->Scale(1./h_aj_data[i][j]->Integral(),"width");//;
	  h_aj_sim[i][j]->Scale(1./h_aj_sim[i][j]->Integral(),"width");//;
	  h_aj_truth[i][j]->Scale(1./h_aj_truth[i][j]->Integral(),"width");//;

	  h_aj_sim[i][j]->Fit("fgg","RQ", "",-0.3, 0.3);
	  float stdev = fgg->GetParameter(2)*sqrt(2);
	  float sedev = fgg->GetParError(2)*sqrt(2);
	  
	  h_sim_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_sim_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);

	  h_jer_sim_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_sim_im_3jet[j]->SetBinError(i+1, sedev);

	  h_aj_data[i][j]->Fit("fgg","RQ", "",-0.3, 0.3);
	  stdev = fgg->GetParameter(2)*sqrt(2);
	  sedev = fgg->GetParError(2)*sqrt(2);

	  //stdev = fgg->GetParameter(2)*sqrt(2);
	  //sedev = fgg->GetParError(2)*sqrt(2);

	  h_jer_data_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_data_im_3jet[j]->SetBinError(i+1, sedev);
	  
	  h_data_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_data_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);

	  h_aj_truth[i][j]->Fit("fgg","RQ", "",-0.15, 0.15);
	  stdev = fgg->GetParameter(2)*sqrt(2);
	  sedev = fgg->GetParError(2)*sqrt(2);

	  h_jer_truth_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_truth_im_3jet[j]->SetBinError(i+1, sedev);
	  
	  h_truth_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_truth_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);

	}
    }
  TF1 *fline = new TF1("fline","pol1", 0, 10);
  // Linear Fit And Extrapolation of soft corrections
  std::cout << "---------------- Estrapolation -------------------" <<std::endl; 
  for (int i = 0; i < 4 ; i++)
    {


      h_sim_3jet_width[i]->Fit("fline","R+","", 3, 10);
      float sim_at_zero = fline->Eval(0);
      for (int j = 0; j < 8; j++)
	{
	  
	  float sim_at_ten = fline->Eval(third_jet_cuts[j]);
	  float ksoft_sim = sim_at_zero/sim_at_ten;
	  h_sim_3jet_softcorr[j]->SetBinContent(i+1, ksoft_sim);
	  float res = h_jer_sim_im_3jet[j]->GetBinContent(i+1)*ksoft_sim;
	  float sigres = h_jer_sim_im_3jet[j]->GetBinError(i+1);
	  	  
	  std::cout << "sim " << i << " , " << j << " --> " << ksoft_sim << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_sim_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_sim_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      h_data_3jet_width[i]->Fit("fline","R+","", 3, 10);
      float data_at_zero = fline->Eval(0);
      //      float data_at_ten = fline->Eval(10);
      //float ksoft_data = data_at_zero/data_at_ten;
      for (int j = 0; j < 8; j++)
	{
	  float data_at_ten = fline->Eval(third_jet_cuts[j]);
		  //	  float data_at_ten = fline->Eval(h_jer_data_im_3jet[j]->GetBinCenter(i+1));
	  float ksoft_data = data_at_zero/data_at_ten;
	  h_data_3jet_softcorr[j]->SetBinContent(i+1, ksoft_data);
	  float res = h_jer_data_im_3jet[j]->GetBinContent(i+1)*ksoft_data;
	  float sigres = h_jer_data_im_3jet[j]->GetBinError(i+1)*ksoft_data;

	  std::cout << "data " << i << " , " << j << " --> " << ksoft_data << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_data_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_data_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      h_truth_3jet_width[i]->Fit("fline","R+","", 3, 10);
      float truth_at_zero = fline->Eval(0);
      //float truth_at_ten = fline->Eval(10);
      //float ksoft_truth = truth_at_zero/truth_at_ten;
      for (int j = 0; j < 8; j++)
	{
	  float truth_at_ten = fline->Eval(third_jet_cuts[j]);
	  float ksoft_truth = truth_at_zero/truth_at_ten;
	  h_truth_3jet_softcorr[j]->SetBinContent(i+1, ksoft_truth);
	  float res = h_jer_truth_im_3jet[j]->GetBinContent(i+1)*ksoft_truth;
	  float sigres = h_jer_truth_im_3jet[j]->GetBinError(i+1)*ksoft_truth;

	  std::cout << "truth " << i << " , " << j << " --> " << ksoft_truth << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_truth_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_truth_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      // No begins the subtraction of the truth
  
      for (int j = 0; j < 8; j++)
	{
	  float tres = h_jer_truth_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float tsigres = h_jer_truth_im_softcorr_3jet[j]->GetBinError(i+1);

	  float sres = h_jer_sim_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float ssigres = h_jer_sim_im_softcorr_3jet[j]->GetBinError(i+1);

	  float dres = h_jer_data_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float dsigres = h_jer_data_im_softcorr_3jet[j]->GetBinError(i+1);

	  float sim_new_res = sqrt(sres*sres - tres*tres);
	  float data_new_res = sqrt(dres*dres - tres*tres);
	  float data_new_sigres = sqrt(TMath::Power(tres*tsigres/data_new_res, 2) + TMath::Power(dres*dsigres/data_new_res, 2));
	  float sim_new_sigres = sqrt(TMath::Power(tres*tsigres/sim_new_res, 2) + TMath::Power(sres*ssigres/sim_new_res, 2));

	  std::cout << "data " << i << " , " << j << " --> " << dres << " - " << tres << " --> " << data_new_res << std::endl;
	  std::cout << "sim " << i << " , " << j << " --> " << sres << " - " << tres << " --> " << sim_new_res << std::endl;
	  h_jer_sim_im_softcorr_part_3jet[j]->SetBinContent(i+1, sim_new_res);
	  h_jer_sim_im_softcorr_part_3jet[j]->SetBinError(i+1, sim_new_sigres);
	  h_jer_data_im_softcorr_part_3jet[j]->SetBinContent(i+1, data_new_res);
	  h_jer_data_im_softcorr_part_3jet[j]->SetBinError(i+1, data_new_sigres);
	}

    }

  // 
  int color_sim = kRed;
  int color_truth = kBlack;
  int color_data = kBlue;

  TCanvas *ctt = new TCanvas("ctt","ctt", 500, 500);
  for (int i3 = 0; i3 < 8; i3++)
    {
      dlutility::SetMarkerAtt(h_data_3jet_softcorr[i3], color_data, 1, 8);
      dlutility::SetLineAtt(h_data_3jet_softcorr[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_truth_3jet_softcorr[i3], color_truth, 1, 8);
      dlutility::SetLineAtt(h_truth_3jet_softcorr[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_softcorr[i3], color_sim, 1, 8);
      dlutility::SetLineAtt(h_sim_3jet_softcorr[i3], color_sim, 1, 1);
      h_data_3jet_softcorr[i3]->SetMaximum(1.5);
      h_data_3jet_softcorr[i3]->SetMinimum(0);
      h_data_3jet_softcorr[i3]->Draw("same p");
      h_sim_3jet_softcorr[i3]->Draw("same p");
      h_truth_3jet_softcorr[i3]->Draw("same p");

      dlutility::DrawSPHENIXpp(0.22, 0.82);
      dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
      TLegend *legg = new TLegend(0.6, 0.77, 0.75, 0.87);
      legg->SetLineWidth(0);
      legg->AddEntry(h_data_3jet_softcorr[i3], "Data");
      legg->AddEntry(h_sim_3jet_softcorr[i3], "Sim");
      legg->AddEntry(h_truth_3jet_softcorr[i3], "Truth");

      legg->Draw("same");

      ctt->Print(Form("%s/jer_plots/softcorr_3jet%d_r%02d.pdf", rb.get_code_location().c_str(), i3, cone_size));
      ctt->Print(Form("%s/jer_plots/softcorr_3jet%d_r%02d.png", rb.get_code_location().c_str(), i3, cone_size));

    }
  
  
  TH1D *h_blank_res = new TH1D("hh", "; p_{T.3} [GeV]; #sigma(p_{T})/p_{T}", 1, 0, 10);
  h_blank_res->SetMinimum(0.0);
  h_blank_res->SetMaximum(0.5);
  
  TCanvas *ctg = new TCanvas("ctg", "ctg", 500, 500);
  for (int i = 0; i < 4; i++)
    {

      dlutility::SetLineAtt(h_sim_3jet_bi_width[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_bi_width[i], color_sim, 0.8, 21);
      dlutility::SetLineAtt(h_data_3jet_bi_width[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_bi_width[i], color_data, 0.8, 21);

      dlutility::SetLineAtt(h_sim_3jet_width[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_width[i], color_sim, 0.8, 20);
      dlutility::SetLineAtt(h_truth_3jet_width[i], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_truth_3jet_width[i], color_truth, 0.8, 20);
      dlutility::SetLineAtt(h_data_3jet_width[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_width[i], color_data, 0.8, 20);

      h_sim_3jet_width[i]->GetFunction("fline")->SetLineColor(color_sim);
      h_data_3jet_width[i]->GetFunction("fline")->SetLineColor(color_data);
      h_truth_3jet_width[i]->GetFunction("fline")->SetLineColor(color_truth);

      h_sim_3jet_width[i]->GetFunction("fline")->SetRange(0, 10);
      h_data_3jet_width[i]->GetFunction("fline")->SetRange(0, 10);
      h_truth_3jet_width[i]->GetFunction("fline")->SetRange(0, 10);

      h_blank_res->DrawCopy();
      h_data_3jet_width[i]->Draw("same p");
      h_sim_3jet_width[i]->Draw("same p");
      h_truth_3jet_width[i]->Draw("same p");
      h_data_3jet_bi_width[i]->Draw("same p");
      h_sim_3jet_bi_width[i]->Draw("same p");

      dlutility::DrawSPHENIXpp(0.2, 0.82);
      dlutility::drawText(Form("%d < #bar{p}_{T} < %d GeV", 20 + i*5, 25 + i*5), 0.2, 0.7); 

      TLegend *leg = new TLegend(0.55, 0.67, 0.8, 0.87);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_3jet_width[i], "Data - Imbalance");
      leg->AddEntry(h_sim_3jet_width[i], "Sim - Imbalance");
      leg->AddEntry(h_data_3jet_bi_width[i], "Data - Bisector");
      leg->AddEntry(h_sim_3jet_bi_width[i], "Sim - Bisector");

      leg->Draw("same");

      ctg->Print(Form("%s/jer_plots/width3jet_r%02d_pta%d.pdf", rb.get_code_location().c_str(), cone_size, i));
      ctg->Print(Form("%s/jer_plots/width3jet_r%02d_pta%d.png", rb.get_code_location().c_str(), cone_size, i));

    }


    for (int i = 0; i < 4; i++)
      {
	for (int j = 0; j < 8; j++)
	  {
	    ctg->SetLogy(1);
	    dlutility::SetLineAtt(h_aj_sim[i][j], color_sim, 1, 1);
	    dlutility::SetMarkerAtt(h_aj_sim[i][j], color_sim, 0.8, 1);
	    dlutility::SetLineAtt(h_aj_truth[i][j], color_truth, 1, 1);
	    dlutility::SetMarkerAtt(h_aj_truth[i][j], color_truth, 0.8, 1);
	    dlutility::SetLineAtt(h_aj_data[i][j], color_data, 1, 1);
	    dlutility::SetMarkerAtt(h_aj_data[i][j], color_data, 0.8, 8);

	    h_aj_sim[i][j]->SetFillColorAlpha(kRed - 3, 0.3);
	    h_aj_data[i][j]->GetFunction("fgg")->SetLineColor(kBlue + 2);//SetFillColorAlpha(kBlue - 3, 0.3);
	    h_aj_data[i][j]->GetFunction("fgg")->SetLineWidth(2);//SetFillColorAlpha(kBlue - 3, 0.3);
	    h_aj_sim[i][j]->GetFunction("fgg")->SetLineColor(kRed + 2);//SetFillColorAlpha(kBlue - 3, 0.3);
	    h_aj_sim[i][j]->GetFunction("fgg")->SetLineWidth(2);//SetFillColorAlpha(kBlue - 3, 0.3);

	    //	    h_aj_sim[i][j]->SetMaximum(100);
	    h_aj_truth[i][j]->Draw("p");
	    h_aj_sim[i][j]->Draw("same hist");
	    h_aj_sim[i][j]->Draw("same p");
	    h_aj_truth[i][j]->Draw("same p");
	    h_aj_data[i][j]->Draw("same p");
	    dlutility::DrawSPHENIXpp(0.22, 0.85);
	    dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = %0.1f}", cone_size*0.1), 0.22, 0.75);
	    dlutility::drawText(Form("%d < #bar{p}_{T} #leq %d GeV", 20+i*5, 25 + i*5), 0.22, 0.7);
	    dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[j]), 0.22, 0.65);
	    TLegend *leg = new TLegend(0.6, 0.77, 0.75, 0.87);
	    leg->SetLineWidth(0);
	    leg->AddEntry(h_aj_data[i][j], "Data");
	    leg->AddEntry(h_aj_sim[i][j], "Sim");
	    leg->Draw("same");

	    ctg->Print(Form("%s/jer_plots/aj_r%02d_3jet%d_pta%d.pdf", rb.get_code_location().c_str(), cone_size, j, i));
	    ctg->Print(Form("%s/jer_plots/aj_r%02d_3jet%d_pta%d.png", rb.get_code_location().c_str(), cone_size, j, i));
	  }
      }

  h_blank_res->SetMaximum(10);
  TCanvas *ctg2 = new TCanvas("ctg2", "ctg2", 500, 500);
  for (int i = 0; i < 4; i++)
    {

      dlutility::SetLineAtt(h_sim_3jet_sigp[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigp[i], color_sim, 0.8, 8);

      dlutility::SetLineAtt(h_data_3jet_sigp[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigp[i], color_data, 0.8, 8);

      dlutility::SetLineAtt(h_sim_3jet_sigv[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigv[i], color_sim, 0.8, 21);

      dlutility::SetLineAtt(h_data_3jet_sigv[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigv[i], color_data, 0.8, 21);

      dlutility::SetLineAtt(h_sim_3jet_sigma[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigma[i], color_sim, 0.8, 24);

      dlutility::SetLineAtt(h_data_3jet_sigma[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigma[i], color_data, 0.8, 24);

      h_blank_res->DrawCopy();
      h_data_3jet_sigp[i]->Draw("same p");
      h_sim_3jet_sigp[i]->Draw("same p");
      h_data_3jet_sigv[i]->Draw("same p");
      h_sim_3jet_sigv[i]->Draw("same p");
      h_data_3jet_sigma[i]->Draw("same p");
      h_sim_3jet_sigma[i]->Draw("same p");

      dlutility::DrawSPHENIXpp(0.22, 0.85);
      dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = %0.1f}", cone_size*0.1), 0.22, 0.75);
      dlutility::drawText(Form("%d < #bar{p}_{T} #leq %d GeV", 20+i*5, 25 + i*5), 0.22, 0.7);
      TLegend *legtg = new TLegend(0.6, 0.57, 0.75, 0.877);
      legtg->SetLineWidth(0);
      legtg->SetTextSize(0.04);
      legtg->SetTextFont(42);
      legtg->AddEntry(h_data_3jet_sigp[i], "Data #sigma(p_{T, #psi})");
      legtg->AddEntry(h_data_3jet_sigv[i],  "Data #sigma(p_{T, #eta})");
      legtg->AddEntry(h_data_3jet_sigma[i], "Data #sigma(p_{T, #psi} - p_{T, #eta} ");
      legtg->AddEntry(h_sim_3jet_sigp[i], "Sim #sigma(p_{T, #psi})");
      legtg->AddEntry(h_sim_3jet_sigv[i],  "Sim #sigma(p_{T, #eta})");
      legtg->AddEntry(h_sim_3jet_sigma[i], "Sim #sigma(p_{T, #psi} - p_{T, #eta} ");
      legtg->Draw("same");

      ctg2->Print(Form("%s/jer_plots/sigma3jet_r%02d_pta%d.pdf", rb.get_code_location().c_str(), cone_size, i));
      ctg2->Print(Form("%s/jer_plots/sigma3jet_r%02d_pta%d.png", rb.get_code_location().c_str(), cone_size, i));

    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *ctg1 = new TCanvas("ctg1", "ctg1", 600, 1200);
  ctg1->Divide(2,4);
  for (int i = 0; i < 4; i++)
    {
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      ctg1->cd(2*i+1);

      h2_aj_avgpt_data_3jet[i]->Draw("colz");
      gPad->SetLogz();
      ctg1->cd(2*i+2);
      h2_aj_avgpt_sim_3jet[i]->Draw("colz");
      gPad->SetLogz();

    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c1 = new TCanvas("c1","c1", 500, 500);

  h_ptavg_sim->Scale(1./h_ptavg_sim->Integral(), "width");
  h_ptavg_data->Scale(1./h_ptavg_data->Integral(), "width");

  dlutility::SetLineAtt(h_ptavg_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_ptavg_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_ptavg_sim, color_sim, 1, 1);
  dlutility::SetMarkerAtt(h_ptavg_data, color_data, 1, 1);

  c1->SetLogy();
  h_ptavg_data->Draw();
  h_ptavg_sim->Draw("same");
  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);

  dlutility::SetLineAtt(hp_trigger_ptv_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_trigger_ptv_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_trigger_ptv_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_trigger_ptv_pta_data, color_data, 1, 8);

  hp_trigger_ptv_pta_sim->Draw("p E1");
  hp_trigger_ptv_pta_data->Draw("same p E1");
  TCanvas *c3 = new TCanvas("c3","c3", 500, 500);

  dlutility::SetLineAtt(hp_trigger_ptp_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_trigger_ptp_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_trigger_ptp_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_trigger_ptp_pta_data, color_data, 1, 8);

  hp_trigger_ptp_pta_sim->Draw("p E1");
  hp_trigger_ptp_pta_data->Draw("same p E1");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c4 = new TCanvas("c4","c4", 500, 500);

  dlutility::SetLineAtt(hp_cosdphi_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_cosdphi_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_cosdphi_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_cosdphi_pta_data, color_data, 1, 8);

  hp_cosdphi_pta_sim->Draw("p E1");
  hp_cosdphi_pta_data->Draw("same p E1");

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < h_jer_sim->GetNbinsX(); i++)
    {
      TH1D *h_simv = (TH1D*) h_trigger_ptv_pta_sim->ProjectionY("hv", i+1, i+1);
      h_simv->Fit("fg", "R");
      float sigv = fg->GetParameter(2);
      float sig2v = fg->GetParError(2);
      TH1D *h_simp = (TH1D*) h_trigger_ptp_pta_sim->ProjectionY("hp", i+1, i+1);
      h_simp->Fit("fg", "R");
      float sigp = fg->GetParameter(2);
      float sig2p = fg->GetParError(2);

      h_sigv_sim->SetBinContent(i+1, sigv);
      h_sigv_sim->SetBinError(i+1, sig2v);
 
      h_sigp_sim->SetBinContent(i+1, sigp);
      h_sigp_sim->SetBinError(i+1, sig2p);

      float pta = hp_trigger_ptv_pta_sim->GetBinCenter(i+1);
      float sigpta = hp_trigger_ptv_pta_sim->GetBinWidth(i+1)/2.;
      float cosdphi = hp_cosdphi_pta_sim->GetBinContent(i+1);
      float sigcosdphi = hp_cosdphi_pta_sim->GetBinError(i+1);
      float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
      float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
      std::cout << "Sim " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
      if (sigp <= sigv) continue;
      h_jer_sim->SetBinContent(i+1, res);
      h_jer_sim->SetBinError(i+1, sigres);

    }


  for (int i = 0; i < h_jer_data->GetNbinsX(); i++)
    {
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      TH1D *h_datav = (TH1D*) h_trigger_ptv_pta_data->ProjectionY("hv", i+1, i+1);
      h_datav->Fit("fg", "R");
      float sigv = fg->GetParameter(2);
      float sig2v = fg->GetParError(2);
      TH1D *h_datap = (TH1D*) h_trigger_ptp_pta_data->ProjectionY("hp", i+1, i+1);
      h_datap->Fit("fg", "R");
      float sigp = fg->GetParameter(2);
      float sig2p = fg->GetParError(2);
 
      h_sigp_data->SetBinContent(i+1, sigp);
      h_sigp_data->SetBinError(i+1, sig2p);
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      float pta = hp_trigger_ptv_pta_data->GetBinCenter(i+1);
      float sigpta = hp_trigger_ptv_pta_data->GetBinWidth(i+1)/2.;
      float cosdphi = hp_cosdphi_pta_data->GetBinContent(i+1);
      float sigcosdphi = hp_cosdphi_pta_data->GetBinError(i+1);
      float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
      float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
      h_jer_data->SetBinContent(i+1, res);
      h_jer_data->SetBinError(i+1, sigres);
      if (sigp <= sigv) continue;
      std::cout << "Data " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int iem = 0; iem < 5; iem++)
    {
      for (int i = 0; i < h_em_jer_sim[iem]->GetNbinsX(); i++)
	{
	  h_em_trigger_ptv_pta_sim[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  h_em_trigger_ptp_pta_sim[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);

	  h_em_sigv_sim[iem]->SetBinContent(i+1, sigv);
	  h_em_sigv_sim[iem]->SetBinError(i+1, sig2v);
 
	  h_em_sigp_sim[iem]->SetBinContent(i+1, sigp);
	  h_em_sigp_sim[iem]->SetBinError(i+1, sig2p);
	  float pta = hp_em_trigger_ptv_pta_sim[iem]->GetBinCenter(i+1);
	  float sigpta = hp_em_trigger_ptv_pta_sim[iem]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_em_cosdphi_pta_sim[iem]->GetBinContent(i+1);
	  float sigcosdphi = hp_em_cosdphi_pta_sim[iem]->GetBinError(i+1);
	  float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + TMath::Power(sigcosdphi/cosdphi, 2)/4.));
	  h_em_jer_sim[iem]->SetBinContent(i+1, res);
	  h_em_jer_sim[iem]->SetBinError(i+1, sigres);
	  //std::cout << "Sim " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
	}

      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      for (int i = 0; i < h_em_jer_data[iem]->GetNbinsX(); i++)
	{
	  h_em_trigger_ptv_pta_data[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  h_em_trigger_ptp_pta_data[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);
	  h_em_sigv_data[iem]->SetBinContent(i+1, sigv);
	  h_em_sigv_data[iem]->SetBinError(i+1, sig2v);
 
	  h_em_sigp_data[iem]->SetBinContent(i+1, sigp);
	  h_em_sigp_data[iem]->SetBinError(i+1, sig2p);

	  float pta = hp_em_trigger_ptv_pta_data[iem]->GetBinCenter(i+1);
	  float sigpta = hp_em_trigger_ptv_pta_data[iem]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_em_cosdphi_pta_data[iem]->GetBinContent(i+1);
	  float sigcosdphi = hp_em_cosdphi_pta_data[iem]->GetBinError(i+1);
	  float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
	  h_em_jer_data[iem]->SetBinContent(i+1, res);
	  h_em_jer_data[iem]->SetBinError(i+1, sigres);
	  //	  std::cout << "Data " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
	}

    }
  TF1 *f = new TF1("f", "sqrt(TMath::Power([2]/x, 2) + TMath::Power([0]/sqrt(x), 2) + TMath::Power([1],2))", 20, 40);
  f->SetParameters(0.1, 0.14, 2);
  f->SetParLimits(0, 0.3, 2.0);
  f->SetParLimits(2, 1.2, 3.5);
  f->SetParLimits(1, 0.02, 0.3);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c5 = new TCanvas("c5","c5", 680, 900);
  dlutility::ratioPanelCanvas(c5, 0.4);

  c5->cd(1);
  gPad->SetBottomMargin(0.14);
  dlutility::SetLineAtt(h_jer_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_jer_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_jer_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(h_jer_data, color_data, 1, 8);
  dlutility::SetFont(h_jer_data, 42, 0.05, 0.06, 0.05, 0.05);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_data->GetYaxis()->SetTitleOffset(1.2);
  h_jer_data->SetMinimum(0);
  h_jer_data->SetMaximum(0.4);
  h_jer_data->Fit("f","R+");
  TF1 *f_data = (TF1*) f->Clone();
  f_data->SetName("f_data");
  f_data->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
  f->SetParameters(0.6313, 0.09508, 2.16);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->Fit("f","R+");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TF1 *f_sim = (TF1*) f->Clone();
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  f_sim->SetName("f_sim");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  f_sim->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TF1 *f_sim_calib = (TF1*) f->Clone();
  f_sim_calib->SetName("f_sim_calib");
  f_sim_calib->SetParameters(0.6313, 0.09508, 2.16);
  dlutility::SetLineAtt(f_sim_calib, kBlack, 2, 1);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  
  h_jer_data->GetFunction("f")->SetLineColor(kBlue);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->GetFunction("f")->SetLineColor(kRed);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_data->Draw("p");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->Draw("same p");
  f_sim_calib->Draw("same");
  dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
  dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.%d}", cone_size), 0.55, 0.7, 0, kBlack, 0.05);
  dlutility::drawText("Jet 8 GeV", 0.55, 0.65, 0, kBlack, 0.05);
  dlutility::drawText("N_{Jet} = 2", 0.55, 0.6, 0, kBlack, 0.05);
  TLegend *leg1 = new TLegend(0.20, 0.204, 0.455, 0.354);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.035);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_jer_data,Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data->GetFunction("f")->GetParameter(0), h_jer_data->GetFunction("f")->GetParameter(2), h_jer_data->GetFunction("f")->GetParameter(1)));
  leg1->AddEntry(h_jer_sim,Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim->GetFunction("f")->GetParameter(0), h_jer_sim->GetFunction("f")->GetParameter(2),h_jer_sim->GetFunction("f")->GetParameter(1)));//, h_jer_sim->GetFunction("f")->GetParameter(2)));
  leg1->Draw("same");  
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  c5->cd(2);
  gPad->SetTopMargin(0.05);

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TF1 *fr_compare = new TF1("fr_compare", "(f_data - f_sim)/f_sim", 20, 40);
  fr_compare->SetLineColor(kRed);
  fr_compare->SetParameters(f_data->GetParameter(0), f_data->GetParameter(1), f_data->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2));

  TF1 *fr_diff = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", 20, 40);
  fr_diff->SetLineColor(kRed);
  fr_diff->SetParameters(f_data->GetParameter(0), f_data->GetParameter(1), f_data->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2));

  TH1D *h_compare = (TH1D*) h_jer_data->Clone();
  h_compare->Reset();

  TH1D *h_dcompare = (TH1D*) h_jer_data->Clone();
  h_dcompare->Reset();
  h_dcompare->SetName("h_dcompare");

  for (int i = 0; i < h_compare->GetNbinsX(); i++)
    {
      float va = h_jer_data->GetBinContent(i+1);
      float ea = h_jer_data->GetBinError(i+1);
      float vb = h_jer_sim->GetBinContent(i+1);
      float eb = h_jer_sim->GetBinError(i+1);
      float rat = (va - vb)/vb;
            
      h_compare->SetBinContent(i+1, rat);
      float err = rat*sqrt(((va*va + vb*vb)*ea*ea + 2*(vb*vb - va*vb)*eb*eb)/(vb*vb*(va - vb)*(va - vb)));

      h_compare->SetBinError(i+1, err);
      
      if (va > vb)
	{
	  float sv = sqrt(va*va - vb*vb);
	  float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	  h_dcompare->SetBinContent(i+1, sv);
	  h_dcompare->SetBinError(i+1, se);
	}
    }
  dlutility::SetLineAtt(h_compare, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_compare, kRed, 1, 8);
  dlutility::SetLineAtt(h_dcompare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_dcompare, kBlack, 1, 8);

  dlutility::SetLineAtt(fr_diff, kBlack, 2, 1);

  dlutility::SetFont(h_compare, 42, 0.07, 0.077, 0.07, 0.07);
  h_compare->SetTitle("; <p_{T}> [ GeV ] ; (Data - MC)/MC");
  h_compare->GetYaxis()->SetTitleOffset(0.8);
  h_compare->SetMinimum(-0.3);
  h_compare->SetMaximum(0.6);
  h_compare->Draw("p E0");
  h_dcompare->Draw("p E0 same");
  fr_compare->Draw("same");
  fr_diff->Draw("same");
  TLine *linep = new TLine(20, 0.2, 40, 0.2);
  dlutility::SetLineAtt(linep, kBlack, 3, 4);
  linep->Draw("same");
  TLine *linen = new TLine(20, -0.2, 40, -0.2);
  dlutility::SetLineAtt(linen, kBlack, 3, 4);
  linen->Draw("same");
  TLine *lin0 = new TLine(20, 0, 40, 0);
  dlutility::SetLineAtt(lin0, kBlack, 2, 1);
  lin0->Draw("same");
  TLegend *legl = new TLegend(0.2, 0.78, 0.5, 0.89);
  legl->SetLineWidth(0);
  legl->SetTextFont(42);
  legl->SetTextSize(0.05);
  legl->AddEntry(fr_compare, "Fractional Difference");
  legl->AddEntry(fr_diff, "Quadrature Difference");
  legl->Draw("same");
  c5->Print(Form("%s/jer_plots/datasim_jer_r%02d.png", rb.get_code_location().c_str(), cone_size));
  c5->Print(Form("%s/jer_plots/datasim_jer_r%02d.pdf", rb.get_code_location().c_str(), cone_size));

  // JER with third jet cuts
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);
      dlutility::SetLineAtt(h_jer_sim_3jet[i3], color_sim, 1, 1);
      dlutility::SetLineAtt(h_jer_data_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_3jet[i3], color_sim, 1, 8);
      dlutility::SetMarkerAtt(h_jer_data_3jet[i3], color_data, 1, 8);
      dlutility::SetFont(h_jer_data_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);
      h_jer_data_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_3jet[i3]->SetMinimum(0);
      h_jer_data_3jet[i3]->SetMaximum(0.4);
      h_jer_data_3jet[i3]->Fit("f","R0");
      TF1 *f_data3 = (TF1*) f->Clone();
      f_data3->SetName("f_data3");
      f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f->SetParameters(0.6313, 0.09508, 2.16);
      h_jer_sim_3jet[i3]->Fit("f","R0");
      TF1 *f_sim3 = (TF1*) f->Clone();
      f_sim3->SetName("f_sim3");
      f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      h_jer_data_3jet[i3]->GetFunction("f")->SetLineColor(kBlue);
      h_jer_sim_3jet[i3]->GetFunction("f")->SetLineColor(kRed);
      h_jer_data_3jet[i3]->Draw("p");
      h_jer_sim_3jet[i3]->Draw("same p");
      f_sim_calib->Draw("same");
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.%d}", cone_size), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText("Jet 8 GeV", 0.55, 0.65, 0, kBlack, 0.05);
      dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.6, 0, kBlack, 0.05);
      TLegend *leg2 = new TLegend(0.20, 0.204, 0.455, 0.354);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.035);
      leg2->SetTextFont(42);
      leg2->AddEntry(h_jer_data_3jet[i3],Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(2), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(1)));
      leg2->AddEntry(h_jer_sim_3jet[i3],Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2),h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(1)));//, h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2)));
      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", 20, 40);
      fr_compare3->SetLineColor(kRed);
      fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", 20, 40);
      fr_diff3->SetLineColor(kRed);
      fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TH1D *h_compare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");

      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_dcompare3->SetBinContent(i+1, sv);
	      h_dcompare3->SetBinError(i+1, se);
	    }
	}
      dlutility::SetLineAtt(h_compare3, kRed, 1, 1);
      dlutility::SetMarkerAtt(h_compare3, kRed, 1, 8);
      dlutility::SetLineAtt(h_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_dcompare3, kBlack, 1, 8);

      dlutility::SetLineAtt(fr_diff3, kBlack, 2, 1);

      dlutility::SetFont(h_dcompare3, 42, 0.07, 0.077, 0.07, 0.07);
      h_dcompare3->SetTitle("; <p_{T}> [ GeV ] ; Quadrature Difference");
      h_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      h_dcompare3->SetMinimum(0);
      h_dcompare3->SetMaximum(0.18);
      h_dcompare3->Draw("p E0");
      fr_diff3->Draw("same");

      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_3jet%d.png", rb.get_code_location().c_str(), cone_size, i3));
      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_3jet%d.pdf", rb.get_code_location().c_str(), cone_size, i3));
    }
  c5 = new TCanvas("c5","c5", 680, 800);
  dlutility::ratioPanelCanvas(c5, 0.3);

  // JER with third jet cuts
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);
      dlutility::SetLineAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 8);

      dlutility::SetLineAtt(h_jer_truth_im_softcorr_3jet[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_jer_truth_im_softcorr_3jet[i3], color_truth, 1, 8);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 8);
      dlutility::SetFont(h_jer_data_im_softcorr_part_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);

      dlutility::SetLineAtt(h_jer_sim_im_softcorr_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_3jet[i3], color_sim, 1, 24);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_3jet[i3], color_data, 1, 24);

      dlutility::SetLineAtt(h_jer_sim_im_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_3jet[i3], color_sim, 1, 25);

      dlutility::SetLineAtt(h_jer_truth_im_3jet[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_jer_truth_im_3jet[i3], color_truth, 1, 25);

      dlutility::SetLineAtt(h_jer_data_im_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_3jet[i3], color_data, 1, 25);


      h_jer_data_im_softcorr_part_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMinimum(0);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMaximum(0.5);
      h_jer_data_im_softcorr_part_3jet[i3]->Fit("f","R0");
      TF1 *f_data3 = (TF1*) f->Clone();
      f_data3->SetName("f_data3");
      f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f->SetParameters(0.6313, 0.09508, 2.16);
      h_jer_sim_im_softcorr_part_3jet[i3]->Fit("f","R0");
      TF1 *f_sim3 = (TF1*) f->Clone();
      f_sim3->SetName("f_sim3");
      f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kBlue);
      h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kRed);
      h_jer_data_im_softcorr_part_3jet[i3]->Draw("p");
      h_jer_sim_im_softcorr_part_3jet[i3]->Draw("same p");
      h_jer_truth_im_softcorr_3jet[i3]->Draw("same p");
      h_jer_data_im_softcorr_3jet[i3]->Draw("same p"); 
      h_jer_sim_im_softcorr_3jet[i3]->Draw("same p"); 

      /* h_jer_truth_im_3jet[i3]->Draw("same p"); */
      //h_jer_data_im_3jet[i3]->Draw("same p");
      //h_jer_sim_im_3jet[i3]->Draw("same p");
      
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.%d}", cone_size), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText("Jet 8 GeV", 0.55, 0.65, 0, kBlack, 0.05);
      dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.6, 0, kBlack, 0.05);
      TLegend *leg2 = new TLegend(0.20, 0.55, 0.5, 0.8);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.035);
      leg2->SetTextFont(42);
      leg2->AddEntry(h_jer_truth_im_softcorr_3jet[i3],"Truth after LE");

      leg2->AddEntry(h_jer_data_im_softcorr_3jet[i3],"Data After LE");
      leg2->AddEntry(h_jer_sim_im_softcorr_3jet[i3],"Sim After LE");
      leg2->AddEntry(h_jer_data_im_softcorr_part_3jet[i3],"Data After LE and TS");
      leg2->AddEntry(h_jer_sim_im_softcorr_part_3jet[i3],"Sim After LE and TS");
      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", 20, 40);
      fr_compare3->SetLineColor(kRed);
      fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", 20, 40);
      fr_diff3->SetLineColor(kRed);
      fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TH1D *h_compare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");

      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_im_softcorr_part_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_dcompare3->SetBinContent(i+1, sv);
	      h_dcompare3->SetBinError(i+1, se);
	    }
	}
      dlutility::SetLineAtt(h_compare3, kRed, 1, 1);
      dlutility::SetMarkerAtt(h_compare3, kRed, 1, 8);
      dlutility::SetLineAtt(h_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_dcompare3, kBlack, 1, 8);

      dlutility::SetLineAtt(fr_diff3, kBlack, 2, 1);

      dlutility::SetFont(h_dcompare3, 42, 0.09);//, 0.09, 0.07, 0.07);
      h_dcompare3->SetTitle("; <p_{T}> [ GeV ] ; Quadrature Difference");
      h_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      h_dcompare3->SetMinimum(0);
      h_dcompare3->SetMaximum(0.2);
      //h_compare3->Draw("p E0");
      h_dcompare3->Draw("p E0");
      //fr_compare3->Draw("same");
      fr_diff3->Draw("same");
      TLegend *legl3 = new TLegend(0.2, 0.78, 0.5, 0.89);
      legl3->SetLineWidth(0);
      legl3->SetTextFont(42);
      legl3->SetTextSize(0.05);
      //      legl3->AddEntry(fr_compare3, "Fractional Difference");
      legl3->AddEntry(fr_diff3, "Quadrature Difference");
      //legl3->Draw("same");
      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_im_3jet%d.png", rb.get_code_location().c_str(), cone_size, i3));
      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_im_3jet%d.pdf", rb.get_code_location().c_str(), cone_size, i3));
    }

   // JER with third jet cuts by two methods
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);

      dlutility::SetLineAtt(h_jer_sim_3jet[i3], color_sim + 2, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_3jet[i3], color_sim + 2, 1, 21);

      dlutility::SetLineAtt(h_jer_data_3jet[i3], color_data + 2, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_3jet[i3], color_data + 2, 1, 21);

      dlutility::SetLineAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 20);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 20);

      dlutility::SetFont(h_jer_data_im_softcorr_part_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);


      h_jer_data_im_softcorr_part_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMinimum(0);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMaximum(0.4);
      h_jer_data_im_softcorr_part_3jet[i3]->Fit("f","R0");
      /* TF1 *f_data3 = (TF1*) f->Clone(); */
      /* f_data3->SetName("f_data3"); */
      /* f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2)); */
      /* f->SetParameters(0.6313, 0.09508, 2.16); */
      /* h_jer_sim_im_softcorr_part_3jet[i3]->Fit("f","R"); */
      /* TF1 *f_sim3 = (TF1*) f->Clone(); */
      /* f_sim3->SetName("f_sim3"); */
      /* f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2)); */
      /* h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kBlue); */
      /* h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kRed); */
      h_jer_data_im_softcorr_part_3jet[i3]->Draw("p");
      h_jer_sim_im_softcorr_part_3jet[i3]->Draw("same p");
      h_jer_data_3jet[i3]->Draw("same p");
      h_jer_sim_3jet[i3]->Draw("same p");
      f_sim_calib->Draw("same");
      //h_jer_truth_im_softcorr_3jet[i3]->Draw("same p");
 /*       h_jer_data_im_softcorr_3jet[i3]->Draw("same p"); */
 /*      h_jer_sim_im_softcorr_3jet[i3]->Draw("same p"); */

      /* h_jer_truth_im_3jet[i3]->Draw("same p"); */
      //h_jer_data_im_3jet[i3]->Draw("same p");
      //h_jer_sim_im_3jet[i3]->Draw("same p");
      
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.%d}", cone_size), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText("Jet 8 GeV", 0.55, 0.65, 0, kBlack, 0.05);
      dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.6, 0, kBlack, 0.05);
      TLegend *leg2 = new TLegend(0.20, 0.204, 0.455, 0.354);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.035);
      leg2->SetTextFont(42);
      leg2->AddEntry(h_jer_data_im_softcorr_part_3jet[i3],"Data - Imbalance");//Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(2), h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(1)));
      leg2->AddEntry(h_jer_sim_im_softcorr_part_3jet[i3], "Sim - Imbalance");//Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(2),h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(1)));//, h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->GetParameter(2)));
      leg2->AddEntry(h_jer_data_3jet[i3],"Data - Bisector");//Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(2), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(1)));
      leg2->AddEntry(h_jer_sim_3jet[i3], "Sim - Bisector");//Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2),h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(1)));//, h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2)));

      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      /* TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", 20, 40); */
      /* fr_compare3->SetLineColor(kRed); */
      /* fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2)); */

      /* TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", 20, 40); */
      /* fr_diff3->SetLineColor(kRed); */
      /* fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2)); */

      TH1D *h_im_compare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_im_compare3->Reset();
      h_im_compare3->SetName("h_im_compare3");
      TH1D *h_im_dcompare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_im_dcompare3->Reset();
      h_im_dcompare3->SetName("h_im_dcompare3");
      TH1D *h_compare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");

      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_dcompare3->SetBinContent(i+1, sv);
	      h_dcompare3->SetBinError(i+1, se);
	    }

	  va = h_jer_data_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  ea = h_jer_data_im_softcorr_part_3jet[i3]->GetBinError(i+1);
	  vb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  eb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinError(i+1);

	  rat = (va - vb)/vb;
            
	  h_im_compare3->SetBinContent(i+1, rat);
	  err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_im_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_im_dcompare3->SetBinContent(i+1, sv);
	      h_im_dcompare3->SetBinError(i+1, se);
	    }

	}

      dlutility::SetLineAtt(h_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_dcompare3, kBlack, 1, 21);

      dlutility::SetLineAtt(h_im_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_im_dcompare3, kBlack, 1, 20);

      dlutility::SetFont(h_dcompare3, 42, 0.07, 0.077, 0.07, 0.07);
      h_dcompare3->SetTitle("; <p_{T}> [ GeV ] ; Quadrature Difference");
      h_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      h_dcompare3->SetMinimum(0);
      h_dcompare3->SetMaximum(0.2);
      h_dcompare3->Draw("p E0");
      h_im_dcompare3->Draw("p E0 same");

      TLegend *legl3 = new TLegend(0.2, 0.78, 0.5, 0.89);
      legl3->SetLineWidth(0);
      legl3->SetTextFont(42);
      legl3->SetTextSize(0.05);
      legl3->AddEntry(h_dcompare3, "Bisector");
      legl3->AddEntry(h_im_dcompare3, "Imbalance");
      legl3->Draw("same");
      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_im_v_bi_3jet%d.png", rb.get_code_location().c_str(), cone_size, i3));
      c5->Print(Form("%s/jer_plots/datasim_jer_r%02d_im_v_bi_3jet%d.pdf", rb.get_code_location().c_str(), cone_size, i3));
    }

  /////////////
  TCanvas *cc = new TCanvas();
  h_dcompare->Draw();
  
  TFile *fint = new TFile(Form("%s/rootfiles/truth_sig.root", rb.get_code_location().c_str()), "r");
  TH1D *h_sigp_truth = (TH1D*) fint->Get("h_sigp");
  TH1D *h_sigv_truth = (TH1D*) fint->Get("h_sigv");
  TFile *fintf = new TFile(Form("%s/rootfiles/truth_sig_noFSR.root", rb.get_code_location().c_str()), "r");
  TH1D *h_sigp_truthfsr = (TH1D*) fintf->Get("h_sigp");
  TH1D *h_sigv_truthfsr = (TH1D*) fintf->Get("h_sigv");

  h_sigp_truth->SetTitle(";<p_{T}> [GeV]; #sigma(p_{T,#psi}) or #sigma(p_{T,#eta}) [GeV]");
  
  TCanvas *c10 = new TCanvas("c10","c10", 500, 500);

  dlutility::SetLineAtt(h_sigp_data, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_data, kGreen + 2, 1, 20);
  dlutility::SetLineAtt(h_sigv_data, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_data, kMagenta + 2, 1, 20);
  dlutility::SetLineAtt(h_sigp_sim, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_sim, kGreen + 2, 1, 24);
  dlutility::SetLineAtt(h_sigv_sim, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_sim, kMagenta + 2, 1, 24);
  dlutility::SetLineAtt(h_sigp_truth, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_truth, kGreen + 2, 1, 25);
  dlutility::SetLineAtt(h_sigv_truth, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_truth, kMagenta + 2, 1, 25);

  dlutility::SetLineAtt(h_sigp_truthfsr, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_truthfsr, kGreen + 2, 1, 21);
  dlutility::SetLineAtt(h_sigv_truthfsr, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_truthfsr, kMagenta + 2, 1, 21);

  h_sigp_truth->SetMaximum(10);
  h_sigp_truth->SetMinimum(0);
  dlutility::SetFont(h_sigp_truth, 42, 0.05, 0.05, 0.05, 0.05);
  h_sigp_truth->Draw("p");
  h_sigv_truth->Draw("same p");

  dlutility::DrawSPHENIXppInternalSpace(0.22, 0.85,0.1, 0, 1, 0, 1);
  dlutility::drawText("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.4}", 0.22, 0.70);
  dlutility::drawText("p_{T} > 10 GeV", 0.22, 0.65);
  TLegend *leg = new TLegend(0.67, 0.65, 0.83, 0.877);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_sigp_truth, "#sigma(p_{T, #psi}^{Truth})");
  leg->AddEntry(h_sigv_truth, "#sigma(p_{T, #eta}^{Truth})");

  leg->Draw("same");

  c10->Print(Form("%s/jer_plots/sigpsigv_truth_r%02d.pdf", rb.get_code_location().c_str(), cone_size));
  c10->Print(Form("%s/jer_plots/sigpsigv_truth_r%02d.png", rb.get_code_location().c_str(), cone_size));


  h_sigp_truth->SetMaximum(10);
  h_sigp_truth->SetMinimum(0);
  h_sigp_truth->Draw("p");
  h_sigv_truth->Draw("same p");
  h_sigp_truthfsr->Draw("same p");
  h_sigv_truthfsr->Draw("same p");

  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("p_{T} > 10 GeV", 0.22, 0.75);
  TLegend *leg4 = new TLegend(0.67, 0.65, 0.83, 0.877);
  leg4->SetLineWidth(0);
  leg4->SetTextFont(42);
  leg4->SetTextSize(0.04);
  leg4->AddEntry(h_sigp_truth, "#sigma(p_{T, #psi}^{Truth})");
  leg4->AddEntry(h_sigp_truthfsr, "no FSR #sigma(p_{T, #psi}^{Truth})");
  leg4->AddEntry(h_sigv_truth, "#sigma(p_{T, #eta}^{Truth})");

  leg4->Draw("same");

  c10->Print(Form("%s/jer_plots/sigpsigv_truthfsr.pdf", rb.get_code_location().c_str()));
  c10->Print(Form("%s/jer_plots/sigpsigv_truthfsr.png", rb.get_code_location().c_str()));

  h_sigp_truth->SetMaximum(20);
  h_sigp_truth->SetMinimum(0);
  //ydlutility::SetFont(h_sigp_truth, 42, 0.05);
  h_sigp_truth->Draw("p");
  h_sigv_sim->Draw("same p");
  h_sigp_sim->Draw("same p");
  h_sigv_truth->Draw("same p");

  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("p_{T} > 10 GeV", 0.22, 0.75);
  leg = new TLegend(0.67, 0.65, 0.83, 0.877);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_sigp_sim, "#sigma(p_{T, #psi}^{Reco})");
  leg->AddEntry(h_sigp_truth, "#sigma(p_{T, #psi}^{Truth})");
  leg->AddEntry(h_sigv_sim, "#sigma(p_{T, #eta}^{Reco})");

  leg->Draw("same");

  c10->Print(Form("%s/jer_plots/sigpsigv_sim_r%02d.pdf", rb.get_code_location().c_str(), cone_size));
  c10->Print(Form("%s/jer_plots/sigpsigv_sim_r%02d.png", rb.get_code_location().c_str(), cone_size));

  h_sigp_truth->SetMaximum(20);
  h_sigp_truth->SetMinimum(0);
  h_sigp_truth->Draw("p");
  h_sigv_sim->Draw("same p");
  h_sigp_data->Draw("same p");
  h_sigv_data->Draw("same p");
  h_sigp_sim->Draw("same p");
  h_sigv_truth->Draw("same p");

  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.4}", 0.22, 0.75);
  dlutility::drawText("p_{T} > 10 GeV", 0.22, 0.7);
  //  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg = new TLegend(0.6, 0.64, 0.754, 0.89);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_sigp_data, "#sigma(p_{T, #psi}^{Data})","p");
  leg->AddEntry(h_sigp_sim, "#sigma(p_{T, #psi}^{Reco})","p");
  leg->AddEntry(h_sigp_truth, "#sigma(p_{T, #psi}^{Truth})","p");
  leg->Draw("same");
  leg = new TLegend(0.754, 0.64, 0.908, 0.89);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_sigv_data, "#sigma(p_{T, #eta}^{Data})","p");
  leg->AddEntry(h_sigv_sim, "#sigma(p_{T, #eta}^{Reco})","p");
  leg->AddEntry(h_sigv_truth, "#sigma(p_{T, #eta}^{Truth})","p");

  leg->Draw("same");

  c10->Print(Form("%s/jer_plots/sigpsigv_data_r%02d.pdf", rb.get_code_location().c_str(), cone_size));
  c10->Print(Form("%s/jer_plots/sigpsigv_data_r%02d.png", rb.get_code_location().c_str(), cone_size));

  // Here is the average EM fraction result

  TCanvas *c_em = new TCanvas("c_em", "c_em", 800, 400);
  int color_em[5] = {kRed, kGreen, kRed, kBlack, kBlack};

  h_jer_data->SetMinimum(0);
  h_jer_data->SetMaximum(0.4);

  h_jer_sim->SetMinimum(0);
  h_jer_sim->SetMaximum(0.4);

  dlutility::SetLineAtt(h_jer_sim, color_data, 1, 1);
  dlutility::SetLineAtt(h_jer_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_jer_sim, color_data, 1, 20);
  dlutility::SetMarkerAtt(h_jer_data, color_data, 1, 20);

  dlutility::SetFont(h_jer_sim, 42, 0.04);
  dlutility::SetFont(h_jer_data, 42, 0.04);
  
  for (int iem = 0; iem < nbins_em; iem++)
    {
      dlutility::SetLineAtt(h_em_jer_data[iem], color_em[iem], 1, 1);
      dlutility::SetMarkerAtt(h_em_jer_data[iem], color_em[iem], 1, 20);
      dlutility::SetLineAtt(h_em_jer_sim[iem], color_em[iem], 1, 1);
      dlutility::SetMarkerAtt(h_em_jer_sim[iem], color_em[iem], 1, 20);
    }  

  c_em->Divide(2, 1);
  c_em->cd(1);
  h_jer_sim->Draw("p");
  h_em_jer_sim[0]->Draw("p same");
  h_em_jer_sim[1]->Draw("p same");
  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("p_{T} > 10 GeV", 0.22, 0.75);
  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg->SetLineWidth(0);  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_jer_sim, "ALL SIM","p");
  leg->AddEntry(h_em_jer_sim[0], "Sim <E_{T}^{EM}/E_{T}> < 0.5","p");
  leg->AddEntry(h_em_jer_sim[1], "Sim <E_{T}^{EM}/E_{T}> #geq 0.5","p");
  leg->Draw("same");
  c_em->cd(2);
  h_jer_data->Draw("p");
  h_em_jer_data[0]->Draw("p same");
  h_em_jer_data[1]->Draw("p same");
  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg->SetLineWidth(0);  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_jer_data, "ALL DATA","p");
  leg->AddEntry(h_em_jer_data[0], "Data <E_{T}^{EM}/E_{T}> < 0.5","p");
  leg->AddEntry(h_em_jer_data[1], "Data <E_{T}^{EM}/E_{T}> #geq 0.5","p");
  leg->Draw("same");

  c_em->Print(Form("%s/jer_plots/em_fraction_jer_r%02d.png", rb.get_code_location().c_str(), cone_size));
  c_em->Print(Form("%s/jer_plots/em_fraction_jer_r%02d.pdf", rb.get_code_location().c_str(), cone_size));  


  /// Widths fo sigeta and sigpsi
  TCanvas *cj3 = new TCanvas("cj3","cj3", 500, 500);
  cj3->SetLogy(1);
  for (int i3 = 0 ; i3 < 8 ; i3++)
    {
      for (int i = 0; i < 4; i++)
	{
	  TH1D *h_sigp_sim = (TH1D*) h_trigger_ptp_pta_sim_3jet[i3]->ProjectionY("h_sigp_sim", i+1, i+1);
	  TH1D *h_sigp_data = (TH1D*) h_trigger_ptp_pta_data_3jet[i3]->ProjectionY("h_sigp_data", i+1, i+1);

	  dlutility::SetLineAtt(h_sigp_sim, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_sigp_sim, kRed, 1, 8);

	  dlutility::SetLineAtt(h_sigp_data, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_sigp_data, kBlue, 1, 8);

	  h_sigp_data->Scale(1./h_sigp_data->Integral(), "width");
	  h_sigp_sim->Scale(1./h_sigp_sim->Integral(), "width");

	  h_sigp_data->SetTitle(";#sigma(p_{T, #psi}) [GeV]; Arb.");
	  h_sigp_data->SetFillColorAlpha(kBlue - 3, .3);
	  h_sigp_data->SetMaximum(10);
	  h_sigp_data->Draw("hist");
	  h_sigp_data->Draw("same p");
	  h_sigp_sim->Draw("same");
	  
	  dlutility::DrawSPHENIXpp(0.22, 0.85);
	  dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = %0.1f}", cone_size*0.1), 0.22, 0.75);
	  dlutility::drawText(Form("%d < #bar{p}_{T} #leq %d GeV", 20+i*5, 25 + i*5), 0.22, 0.7);
	  dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
	  TLegend *legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
	  legj3->SetLineWidth(0);
	  legj3->SetTextFont(42);
	  legj3->SetTextSize(0.04);
	  legj3->AddEntry(h_sigp_sim, "Pythia-8 Reco");
	  legj3->AddEntry(h_sigp_data, "Data");

	  legj3->Draw("same");

	  cj3->Print(Form("%s/jer_plots/sigp_r%02d_3jet%d_pta%d.pdf", rb.get_code_location().c_str(), cone_size, i3, i));
	  cj3->Print(Form("%s/jer_plots/sigp_r%02d_3jet%d_pta%d.png", rb.get_code_location().c_str(), cone_size, i3, i));

	  TH1D *h_sige_sim = (TH1D*) h_trigger_ptv_pta_sim_3jet[i3]->ProjectionY("h_sige_sim", i+1, i+1);
	  TH1D *h_sige_data = (TH1D*) h_trigger_ptv_pta_data_3jet[i3]->ProjectionY("h_sige_data", i+1, i+1);

	  dlutility::SetLineAtt(h_sige_sim, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_sige_sim, kRed, 1, 8);

	  dlutility::SetLineAtt(h_sige_data, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_sige_data, kBlue, 1, 8);

	  h_sige_data->Scale(1./h_sige_data->Integral(), "width");
	  h_sige_sim->Scale(1./h_sige_sim->Integral(), "width");

	  h_sige_data->SetTitle(";#sigma(p_{T, #eta}) [GeV]; Arb.");
	  h_sige_data->SetFillColorAlpha(kBlue - 3, .3);
	  h_sige_data->SetMaximum(10);
	  h_sige_data->Draw("hist");
	  h_sige_data->Draw("same p");
	  h_sige_sim->Draw("same");
	  
	  dlutility::DrawSPHENIXpp(0.22, 0.85);
	  dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.2]{#it{R} = %0.1f}", cone_size*0.1), 0.22, 0.75);
	  dlutility::drawText(Form("%d < #bar{p}_{T} #leq %d GeV", 20+i*5, 25 + i*5), 0.22, 0.7);
	  dlutility::drawText(Form("p_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
	  legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
	  legj3->SetLineWidth(0);
	  legj3->SetTextSize(0.04);
	  legj3->SetTextFont(42);
	  legj3->AddEntry(h_sige_sim, "Pythia-8 Reco");
	  legj3->AddEntry(h_sige_data, "Data");

	  legj3->Draw("same");

	  cj3->Print(Form("%s/jer_plots/sige_r%02d_3jet%d_pta%d.pdf", rb.get_code_location().c_str(), cone_size, i3, i));
	  cj3->Print(Form("%s/jer_plots/sige_r%02d_3jet%d_pta%d.png", rb.get_code_location().c_str(), cone_size, i3, i));

	  
	}
    }
  int color_sample[4] = {kBlue, kGreen, kRed, kOrange};
  int color_combined = kBlack;
  for (int i = 0; i < 4; i++)
    {
      dlutility::SetLineAtt(h_lead_sample[i], color_sample[i], 2, 1);
      dlutility::SetMarkerAtt(h_lead_sample[i], color_sample[i], 1.1, 20);
    }
  dlutility::SetLineAtt(h_lead_combined, color_combined, 2, 1);
  dlutility::SetMarkerAtt(h_lead_combined, color_combined, 1.1, 24);

  TCanvas *ccom = new TCanvas("ccom","ccom", 500, 500);
  ccom->SetLogy();
  h_lead_combined->Draw("p");
  h_lead_sample[0]->Draw("p same");
  h_lead_sample[1]->Draw("p same");
  h_lead_sample[2]->Draw("p same");
  h_lead_sample[3]->Draw("p same");
  h_lead_combined->Draw("p same");

  ccom->Print(Form("%s/jer_plots/sample_combine_r%02d.pdf", rb.get_code_location().c_str(), cone_size));
  return;
}
