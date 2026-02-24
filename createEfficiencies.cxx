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

int createEfficiencies(const std::string configfile = "binning.config" , const int cone_size = 4)
{

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
  
  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();
  //Int_t primer = rb.get_primer();

  //Double_t dphicut = rb.get_dphicut();

  Double_t vtx_cut = rb.get_vtx_cut();
  //Double_t njet_cut = rb.get_njet_cut();

  Int_t zyam_sys = rb.get_zyam_sys();
  Int_t inclusive_sys = rb.get_inclusive_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  Int_t prior_sys = rb.get_prior_sys();
  Int_t use_herwig = rb.get_herwig();
  
  TF1 *finjes = nullptr;

  TFile *infilejes = new TFile(Form("%s/r%02d/jes/jes_fits_r%02d_0JES.root", rb.get_jesr_location().c_str(), cone_size, cone_size), "r");
  finjes = (TF1*) infilejes->Get("fjesprime_2");

  TF1 *gjer = nullptr;
  TH1D *hjersmear = nullptr;

  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TFile *finjer = new TFile(Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)), "r");
   
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
  int nbins_smear = 0;
  //sWrite
  if (!hjersmear)
    {
      std::cout << "no histo for smear in " << Form("%s/r%02d/jer/jer_fits_r%02d_1JES_0_closure.root", rb.get_jesr_location().c_str(), (cone_size == 2 ? 3 : cone_size), (cone_size == 2 ? 3 : cone_size)) << std::endl;
      return 1;
    }

  nbins_smear = hjersmear->GetNbinsX();
  
  if (JES_sys != 0)
    {
      if (JES_sys < 0)
	sys_name = "negJES";
      else if (JES_sys > 0)
      	sys_name = "posJES";
    }

  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  Double_t dpt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      dpt_bins[i] = (Double_t) ipt_bins[i];
      if (verbosity > 1)
	{	  
	  std::cout << ipt_bins[i] << " -- " << " -- " << dpt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}
    }
  Int_t max_reco_bin = rb.get_maximum_reco_bin();

  std::string system_string = "pp";

  int sim_version = (use_herwig ? 10 : 11);
  
  std::string j_file[5];
  j_file[0] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_10_ana509_MDC2-00000028.root";
  j_file[1] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_15_ana509_MDC2-00000028.root";
  j_file[2] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
  j_file[3] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
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
  
      ttree[j]->SetBranchAddress(Form("jet_pt_%d", cone_size), &reco_jet_pt[j]);
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

  
  float cs_10 = (3.997e6);
  float cs_15 = (4.073e5);
  float cs_20 = 6.218e4;
  float cs_30 = (2.502e3);
  float cs_50 = (7.2695);

  float csh_10 = 1.5703e-9;
  float csh_30 = 1.473e-12;
  float scale_factor[5];
  scale_factor[0] = cs_10/cs_50;
  scale_factor[1] = cs_15/cs_50;
  scale_factor[2] = cs_20/cs_50;
  scale_factor[3] = cs_30/cs_50; 
  scale_factor[4] = 1;

  if (use_herwig)
    {
      scale_factor[0] = ((float) n_events[3]/ ((float)n_events[0])) * csh_10/(csh_30 * 4.);
      scale_factor[3] = 1;
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
      boundary_r4[1] = ipt_bins[measure_bins[1]];;
      boundary_r4[3] = ipt_bins[measure_bins[1]];;
      boundary_r4[4] = ipt_bins[nbins];;     
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

  TEfficiency *he_dijet_mbd_hit = new TEfficiency("he_dijet_mbd_hit",";Leading #it{p}^{Truth}_{T} [GeV];Subleading #it{p}^{Truth}_{T} [GeV]; MBD Hit Efficiency", nbins*nbins, 0, nbins*nbins);
  
  TEfficiency *he_dijet = new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_reco = new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_fake = new TEfficiency("he_dijet_fake",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_trigger = new TEfficiency("he_dijet_trigger",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_vertex = new TEfficiency("he_dijet_vertex",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  TEfficiency *he_dijet_vertex_truth = new TEfficiency("he_dijet_vertex_truth",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);


  TRandom *rng = new TRandom(0);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets2;
  std::vector<struct jet> myrecojets2;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;

  for (int isample = 0; isample < 5; isample++)
    {
      if (!use_sample[isample]) continue;
      std::cout << "Sample " << isample << std::endl;
      int entries0 = 0;
      int entries2 = ttree[isample]->GetEntries();
      if (verbosity > 5) entries2 = 100000;
      
      for (int i = entries0; i < entries2; i++)
	{
	  ttree[isample]->GetEntry(i);
	  double event_scale = scale_factor[isample];
	  if (verbosity > 5)
	    {
	      std::cout << "------------- Event " << i << " --------------" << std::endl;
	    }
	  
	  bool has_vertex = true;
	  bool has_mbd_hit = mbd_hit[isample];;
	  has_vertex = (fabs(mbd_vertex_z[isample]) < vtx_cut);

	  bool triggered = ( ( ( gl1_scaled[isample] >> 22 ) & 0x1) == 0x1);

	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (cone_size != 4)
	    {
	      maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
	    }

	  if (maxpttruth < boundary_r4[isample] || maxpttruth >= boundary_r4[isample+1]) continue;

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
	      float temppt = finjes->Eval(reco_jet_pt[isample]->at(j));

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
		  int ib = 2;
		  while (hjersmear->GetBinLowEdge(ib++) < temppt && ib <= nbins_smear)
		    {
		      continue;
		    }

		  float smear1 = hjersmear->GetBinContent(ib - 1);

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

	  bool reco_good = djf.check_dijet_reco(myrecojets2);

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
	  if (truth_good)
	    {
	      he_dijet_reco->Fill(reco_good, mytruthjets2.at(0).pt, event_scale);
	      he_dijet_vertex->Fill(has_vertex, mytruthjets2.at(0).pt, event_scale);
	      he_dijet_mbd_hit->Fill(has_mbd_hit, mytruthjets2.at(0).pt, event_scale);
	      he_dijet->Fill(matched && reco_good, max_truth, event_scale);
	    }
	  
	  if (reco_good)
	    {
	      he_dijet_fake->Fill(matched && truth_good, myrecojets2.at(0).pt, event_scale);
	    }
	}
    }
  TString responsepath = "unfolding_hists/efficiencies_" + system_string + "_r0" + std::to_string(cone_size);
	  
  responsepath += "_" + sys_name;

  responsepath += ".root";
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  
  he_dijet->Write();// new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_mbd_hit->Write();// new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_reco->Write();// new TEfficiency("he_dijet",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_fake->Write();// new TEfficiency("he_dijet_fake",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_trigger->Write();// new TEfficiency("he_dijet_trigger",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_vertex->Write();// new TEfficiency("he_dijet_vertex",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
  he_dijet_vertex_truth->Write();// new TEfficiency("he_dijet_vertex_truth",";Leading #it{p}^{Truth}_{T} [GeV]; Matching Efficiency", 50, 0, 100);
	      
  fr->Close();
 
}

int main(int argc, char *argv[])
{

  std::string config = "binning.config";
  int cone_size = 4;

  int set=0;

  for (int i = 1; i < argc; ++i)
    {
      std::string arg = argv[i];

      if (arg == "-c" && i + 1 < argc)
	{
	  set++;
	  config = argv[++i];  // Next argument as string
	}
      else if (arg == "-r" && i + 1 < argc)
	{
	  set++;
	  cone_size = std::stoi(argv[++i]);  // Convert next argument to double
	}
      else if (arg == "-v" && i+1 < argc)
	{
	  verbosity = std::stoi(argv[++i]);  // Convert next argument to double
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
      std::cout << "    ./createEfficiencies -c binning.config -r 4 " << std::endl;
      return 1;
    }
  
  return createEfficiencies(config, cone_size);
  
}
