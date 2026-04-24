#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "dijetfinder.h"

#include "TProfile.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"


static int verbosity = 0;

int makeHerwig_hist(const int cone_size = 4, const std::string configfile = "binning.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();


  read_binning rb(configfile.c_str());
  dijetfinder djf(cone_size);
  djf.SetVerbosity(verbosity);

  int sim_version = 10;
  std::string j_file[3];
  j_file[0] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
  j_file[1] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
  j_file[2] = rb.get_tntuple_location() + "/TREE_JET_SKIM_r0" + std::to_string(cone_size) + "_v" + std::to_string(sim_version) + "_40_ana509_MDC2-00000028.root";

  float n_events[3];
  TFile *finsim[3];


  TTree *ttree[3];
  ULong64_t gl1_scaled[3];

  std::vector<float> *truth_jet_pt_ref[3] = {0};
  std::vector<float> *truth_jet_pt[3] = {0};
  std::vector<float> *truth_jet_eta[3] = {0};
  std::vector<float> *truth_jet_phi[3] = {0};
  
  std::vector<float> *reco_jet_pt[3] = {0};
  std::vector<float> *reco_jet_emcal[3] = {0};
  std::vector<float> *reco_jet_e[3] = {0};
  std::vector<float> *reco_jet_eta[3] = {0};
  std::vector<float> *reco_jet_eta_det[3] = {0};
  std::vector<float> *reco_jet_phi[3] = {0};

  float truth_vertex_z[3];
  float mbd_vertex_z[3];
  int mbd_hit[3];

  for (int j = 0; j < 3; j++)
    {


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
  float cs[3] = {5.2613e4, 2.0694e3, 1.0510e2};

  float scale_factor[3];
  scale_factor[0] = cs[0]/cs[2];
  scale_factor[1] = cs[1]/cs[2];
  scale_factor[2] = 1;

  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();
  //Int_t primer = rb.get_primer();

  //Double_t dphicut = rb.get_dphicut();

  Double_t vtx_cut = rb.get_vtx_cut();
  //Double_t njet_cut = rb.get_njet_cut();

  std::string sys_name = "HERWIG";

  const int nbins = read_nbins;
  const int nbins_pt = nbins+1;
  float ipt_bins[nbins_pt+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  if (verbosity > 1)
    {
      for (int i = 0 ; i < nbins + 1; i++)
	{
	  std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}
    }

  ipt_bins[nbins_pt] = 100;
  
  Int_t max_reco_bin = rb.get_maximum_reco_bin();

  std::string system_string = "pp";

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  djf.setTruthCuts(truth_leading_cut, truth_subleading_cut);

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

  int exbin = 1;
  int measure_bins[10] = {0};
  int mbins = 3;  
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }

  float boundary_r[4];
  boundary_r[0] = ipt_bins[measure_bins[0]];
  boundary_r[1] = ipt_bins[measure_bins[1]];
  boundary_r[2] = ipt_bins[measure_bins[2]];
  boundary_r[3] = ipt_bins[measure_bins[3]];

  TH1D *h_truth_lead_sample[5];
  for (int i = 0; i < 5; i++)
    {
      h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }  

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2", ";p_{T,1}; p_{T,2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  
  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins_pt*nbins_pt, 0, nbins_pt*nbins_pt);

  TH1D *h_xj_range[3];
  TH1D *h_dphi_range[3];
  for (int i =0; i < 3; i++)
    {
      h_xj_range[i] = new TH1D(Form("h_xj_range_%d", i),"", 7, 0.3, 1.0);
      h_dphi_range[i] = new TH1D(Form("h_dphi_range_%d", i),"", 16, 7*TMath::Pi()/8., TMath::Pi());
    }
  const int nbinsdphi  = 16;
  float min_dphi = 3*TMath::Pi()/4.;
  float max_dphi = TMath::Pi();
  float stepdphi = (max_dphi - min_dphi)/(float)nbinsdphi;

  float idphi_bins[nbinsdphi+1];
  for (int i = 0; i < nbinsdphi+1; i++)
    {
      idphi_bins[i] = min_dphi + i*stepdphi;
    }


  TH3D *h_herwig_pt1pt2dphi = new TH3D("h_herwig_pt1pt2dphi",";#it{p}_{T,1, smear};#it{p}_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  int nbin_response = nbins_pt*nbins_pt;
  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> mytruthjets2;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;

  for (int isample = 0; isample < 3; isample++)
    {

      std::cout << "Sample " << isample << std::endl;
      int entries0 = 0;
      int entries2 = ttree[isample]->GetEntries();
      if (verbosity > 5) entries2 = 100000;
      
      for (int i = entries0; i < entries2; i++)
	{
	  ttree[isample]->GetEntry(i);
	  double event_scale = 1;
	  if (verbosity > 5)
	    {
	      std::cout << "------------- Event " << i << " --------------" << std::endl;
	    }
	  

	  bool has_vertex = (fabs(mbd_vertex_z[isample]) < vtx_cut);

	  bool has_mbd_hit = mbd_hit[isample];

	  bool mb_bad =  (!has_vertex || !has_mbd_hit);

	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (maxpttruth < boundary_r[isample] || maxpttruth >= boundary_r[isample+1]) continue;

	  mytruthjets.clear();
	  mytruthjets2.clear();

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

	  bool truth_good = djf.check_dijet_truth(mytruthjets);

	  if (truth_good)
	    {
	      mytruthjets2 = {mytruthjets.begin(), mytruthjets.begin()+2};
	    }
	  if (truth_good)
	    {
	      truth_good &= (mytruthjets2.at(0).pt < ipt_bins[nbins_pt]);
	    }

	  if (!truth_good)
	    {
	      continue;
	    }

	  if (verbosity > 5)
	    {
	      mytruthjets2.at(0).print();
	    }
	  float max_truth = mytruthjets.at(0).pt;
	  float min_truth = mytruthjets.at(1).pt;

	  float pt1_truth_bin = nbins_pt;
	  float pt2_truth_bin = nbins_pt;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  	  
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
	    }


	  

	  float dphi = fabs(mytruthjets.at(0).phi - mytruthjets.at(1).phi);
	  if (dphi > TMath::Pi())
	    {
	      dphi = 2*TMath::Pi() - dphi;
	    }
	  for (int im = 0 ; im < 3 ; im++)
	    {
	      if (e1 >= ipt_bins[measure_bins[im]] && e1 < ipt_bins[measure_bins[im+1]] && e2 >= measure_subleading_cut)
		{
		  h_xj_range[im]->Fill(e2/e1, scale_factor[isample]);
		  h_dphi_range[im]->Fill(dphi, scale_factor[isample]);
		}

	    }
	  if (mb_bad) continue;  
	
	  h_herwig_pt1pt2dphi->Fill(e1, e2, dphi); 
	  h_pt1pt2->Fill(maxit, minit, scale_factor[isample]);
	  h_pt1pt2->Fill(minit, maxit, scale_factor[isample]);
	  
	  h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, scale_factor[isample]);
	  h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, scale_factor[isample]);//, event_scale);
	
	  h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins_pt*pt2_truth_bin, scale_factor[isample]);
	  h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins_pt*pt1_truth_bin, scale_factor[isample]);

	}
    }
  
  h_flat_truth_pt1pt2->Scale(.5);

  h_flat_truth_to_response_pt1pt2->Scale(.5);

  TString responsepath = "truth_hists/herwig_hists_" + system_string + "_r0" + std::to_string(cone_size);
	  
  responsepath += "_" + sys_name;

  responsepath += ".root";
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  h_herwig_pt1pt2dphi->Write();
  h_flat_truth_pt1pt2->Write();
  h_flat_truth_to_response_pt1pt2->Write();	      
  h_pt1pt2->Write();
  for (int irange = 0; irange < 3; irange++)
    {
      h_xj_range[irange]->Write();
      h_dphi_range[irange]->Write();    
    }
fr->Close();
  

  return 0;
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
      std::cout << "    ./makeHerig_hist -c binning.config -r 4 " << std::endl;
      return 1;
    }

  makeHerwig_hist(cone_size, config);
  return 0;
  
}
