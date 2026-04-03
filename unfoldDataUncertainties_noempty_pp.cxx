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

#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TNtuple.h"

static int verbosity = 0;

int unfoldDataUncertainties_noempty_pp(const std::string configfile = "binning.config", const int niterations = 20, const int cone_size = 4, const int primer = 0)
{

  
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  std::string system_string = "pp";

  read_binning rb(configfile.c_str());

  Int_t prior_sys = rb.get_prior_sys();
  Int_t trigger_sys = rb.get_trigger_sys();

  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  Int_t herwig_sys = rb.get_herwig();
  
  std::string sys_name = "nominal";
  
  if (prior_sys)
    sys_name = "PRIOR";

  if (trigger_sys)
    sys_name = "TRIGGER";

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

  std::string sys_name_orig = sys_name;
  
  if (primer == 1)
    {
      sys_name = "PRIMER1_" + sys_name;
    }
  else if (primer == 2)
    {
      sys_name = "PRIMER2_" + sys_name;
    }

  Int_t read_nbins = rb.get_nbins();
  

  const int nbins = read_nbins;
  const int nbins_pt = read_nbins + 1;

  Double_t ipt_bins[nbins_pt+1];
  Double_t ixj_bins[nbins+1];

  float fipt_bins[nbins_pt+1];
  float fixj_bins[nbins+1];

  rb.get_pt_bins(fipt_bins);
  rb.get_xj_bins(fixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      ipt_bins[i] = fipt_bins[i];
      ixj_bins[i] = fixj_bins[i];
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

  ipt_bins[nbins_pt] = 100;
  Int_t max_reco_bin = rb.get_maximum_reco_bin();
  
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

  float sample_boundary[4] = {0};
  float low_trigger[3] = {0};

  const int mbins = rb.get_measure_bins();
  int measure_bins[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }

  
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ib = 0; ib < 3; ib++)
    {
      low_trigger[ib] = rb.get_low_trigger(ib);
      std::cout <<  low_trigger[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  
  std::string data_file = rb.get_tntuple_location() + "/TREE_DIJET_SKIM_r0" + std::to_string(cone_size) + "_v8_5_ana533_2025p007_v001_gl10-all.root";
  
  ULong64_t gl1_scaled;
  float mbd_vertex_z;
  
  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_emcal = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_eta_det = 0;
  std::vector<float> *reco_jet_phi = 0;

  TFile *fin = new TFile(data_file.c_str(), "r");
  TTree *ttree  = (TTree*) fin->Get("ttree");;

  if (!ttree)
    {
      std::cout << " no data "<< std::endl;
    }

  ttree->SetBranchAddress(Form("jet_pt_calib_%d", cone_size), &reco_jet_pt);
  ttree->SetBranchAddress(Form("jet_emcal_%d", cone_size), &reco_jet_emcal);
  ttree->SetBranchAddress(Form("jet_e_%d", cone_size), &reco_jet_e);
  ttree->SetBranchAddress(Form("jet_eta_%d", cone_size), &reco_jet_eta);
  ttree->SetBranchAddress(Form("jet_eta_det_%d", cone_size), &reco_jet_eta_det);
  ttree->SetBranchAddress(Form("jet_phi_%d", cone_size), &reco_jet_phi);
  ttree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  ttree->SetBranchAddress("gl1_scaled", &gl1_scaled);


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

  TFile *fresponse = new TFile(Form("%s/response_matrices/response_matrix_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"r");

  std::cout << Form("%s/response_matrices/response_matrix_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()) << std::endl;
  TH1D *h_flat_truth_pt1pt2_raw = (TH1D*) fresponse->Get("h_truth_flat_pt1pt2_raw"); 

  if (!h_flat_truth_pt1pt2_raw)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fresponse->Get("h_truth_flat_pt1pt2"); 

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

  
  RooUnfoldResponse *rooResponse = (RooUnfoldResponse*) fresponse->Get("response");

  if (!rooResponse)
    {
      std::cout << "no repsonse" << std::endl;
      return 1;
    }

  TH2D *h_flat_response_pt1pt2 = (TH2D*) fresponse->Get("h_flat_response_pt1pt2"); 

  if (!h_flat_response_pt1pt2)
    {
      std::cout << "no response" << std::endl;
      return 1;
    }



  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_" + system_string + "_r0" + std::to_string(cone_size);

  unfoldpath += "_" + sys_name;
  unfoldpath += ".root";
  
  TFile *funin = new TFile(unfoldpath.Data(),"r");
  
  TH1D *h_flat_data_pt1pt2 = (TH1D*) funin->Get("h_data_flat_pt1pt2");
  if (!h_flat_data_pt1pt2)
    {
      std::cout << " no h_flat_data_pt1pt2 " << std::endl;
      
    }

  TH1D *h_flat_unfold_pt1pt2[niterations];
  
  int niter = 1;

  TProfile2D *hp_pt1pt2_stats[niterations];
  TProfile *hp_xj[niterations];
  TProfile *hp_pt1pt2[niterations];
  TProfile *hp_xj_range[mbins][niterations];

  int nbins_pt1pt2 = nbins_pt*nbins_pt;
  
  for (int iter = 0; iter < niterations; iter++)
    {
      hp_pt1pt2_stats[iter] = new TProfile2D(Form("hp_pt1pt2_stats_%d", iter),";p_{T,1}; p_{T,2}; average +/- rms", nbins_pt, ipt_bins, nbins_pt, ipt_bins, "s");
      hp_pt1pt2[iter] = new TProfile(Form("hp_pt1pt2_%d", iter),";pt1pt2 bin; Avg +/- RMS", nbins_pt1pt2, 0, nbins_pt1pt2,"s");
      hp_xj[iter] = new TProfile(Form("hp_xj_%d", iter),";x_{J}; Avg +/- RMS", nbins, ixj_bins,"s");
      for (int i = 0; i < mbins; i++)
	{
	  hp_xj_range[i][iter] = new TProfile(Form("hp_xj_range_%d_%d", i, iter),";x_{J}; Avg +/- RMS", nbins, ixj_bins,"s");
	}
    }

  int ntoys = 100;
  
  TF1 *ferror_response  =  new TF1("fgaus_truth","gaus");

  ferror_response->SetParameters(1, 0, 1);
  ferror_response->SetRange(-5, 5);

  TH2D *h_pt1pt2_unfold[niterations];
  TH2D *h_pt1pt2_unfold_sym[niterations];
      
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  
  for (int itoy = 0; itoy < ntoys; ++itoy)
    {
      std::cout << "Toy: " << itoy << std::endl;

      TH2D *h_flat_response_clone = (TH2D*) h_flat_response_pt1pt2->Clone();
      for (int ibin = 0; ibin < h_flat_response_clone->GetXaxis()->GetNbins(); ++ibin)
	{
	  for (int jbin = 0; jbin < h_flat_response_clone->GetYaxis()->GetNbins(); ++jbin)
	    {
	      int ijbin = h_flat_response_pt1pt2->GetBin(ibin+1, jbin + 1);
	      float smear = ferror_response->GetRandom() * h_flat_response_pt1pt2->GetBinError(ijbin);
	      float value = h_flat_response_clone->GetBinContent(ijbin);
	      float newvalue = value + smear;
	      if (newvalue < 0) newvalue = 0;
	      h_flat_response_clone->SetBinContent(ijbin, newvalue);
	    }
	}
      
      RooUnfoldResponse rooResponse(h_flat_reco_pt1pt2, h_flat_truth_pt1pt2_raw, h_flat_response_clone);
      for (int iter = 0; iter < niterations; iter++ )
	{
	  RooUnfoldBayes   unfold (&rooResponse, h_flat_data_pt1pt2, iter + 1, false, true);    // OR
	  TH1D *h_flat_unfold_skim = (TH1D*) unfold.Hunfold();
	  h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
	  h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_flat_unfold_pt1pt2[iter]->Reset();
	  
	  histo_opps::fill_up_histo(h_flat_unfold_skim, h_flat_unfold_pt1pt2[iter], h_flat_truth_mapping_primer);


	}
    
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pt1pt2_unfold[iter]->Reset();
	  h_xj_unfold[iter]->Reset();
	  for (int ir = 0; ir < mbins; ir++)
	    {
	      h_xj_unfold_range[ir][iter]->Reset();
	    }
	}
      

      for (int ib = 0; ib < nbins_pt*nbins_pt; ib++)
	{
	  int xbin = ib/nbins_pt;
	  int ybin = ib%nbins_pt;
	  
	  int b = h_pt1pt2_unfold[0]->GetBin(xbin+1, ybin+1);
	  
	  for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinContent(b, h_flat_unfold_pt1pt2[iter]->GetBinContent(ib+1));
	  for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinError(b, h_flat_unfold_pt1pt2[iter]->GetBinError(ib+1));
	}

      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pt1pt2_unfold_sym[iter] = (TH2D*) h_pt1pt2_unfold[iter]->Clone();
	  histo_opps::project_xj(h_pt1pt2_unfold_sym[iter], h_xj_unfold[iter], nbins_pt, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	  for (int irange = 0; irange < mbins; irange++)
	    {
	      histo_opps::project_xj(h_pt1pt2_unfold_sym[iter], h_xj_unfold_range[irange][iter], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	    }
	  for (int ib = 0; ib < nbins; ib++)
	    {
	      hp_xj[iter]->Fill(hp_xj[iter]->GetBinCenter(ib+1), h_xj_unfold[iter]->GetBinContent(ib+1));
	      for (int irange = 0; irange < mbins; irange++)
		{
		  hp_xj_range[irange][iter]->Fill(hp_xj_range[irange][iter]->GetBinCenter(ib+1), h_xj_unfold_range[irange][iter]->GetBinContent(ib+1));
		}
	      for (int ib2 = 0; ib2 < nbins; ib2++)
		{
		  int ijbin = h_pt1pt2_unfold_sym[iter]->GetBin(ib+1, ib2+1);
		  hp_pt1pt2_stats[iter]->Fill(hp_pt1pt2_stats[iter]->GetXaxis()->GetBinCenter(ib+1), hp_pt1pt2_stats[iter]->GetYaxis()->GetBinCenter(ib2+1), h_pt1pt2_unfold_sym[iter]->GetBinContent(ijbin));
		}
	    }
	  for (int ib = 0; ib < nbins_pt1pt2; ib++)
	    {
	      hp_pt1pt2[iter]->Fill(hp_pt1pt2[iter]->GetBinCenter(ib+1), h_flat_unfold_pt1pt2[iter]->GetBinContent(ib+1));
	    }

	}

    }


  int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 
  TCanvas *c = new TCanvas("c","c", 500, 500);
  for (int iter = 0; iter < 5; iter++)
    {
      dlutility::SetLineAtt(hp_xj[iter], colors[iter], 1,1);
      dlutility::SetMarkerAtt(hp_xj[iter], colors[iter], 1,8);
    }

  hp_xj[0]->Draw("E1");
  hp_xj[1]->Draw("E1 same");
  hp_xj[2]->Draw("E1 same");
  hp_xj[3]->Draw("E1 same");
  hp_xj[4]->Draw("E1 same");
  

  
  TFile *fout = new TFile(Form("%s/uncertainties/uncertainties_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, sys_name.c_str()),"recreate");
  TEnv *penv = new TEnv("binning.config");
  penv->Write();
  for (int iter = 0; iter < niterations; iter++)
    {
      hp_xj[iter]->Write();
      hp_pt1pt2[iter]->Write();            
      for (int irange = 0; irange < mbins; irange++)
	{
	  hp_xj_range[irange][iter]->Write();            
	}
      hp_pt1pt2_stats[iter]->Write();      
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
  if (set < 3)
    {
      std::cout << "Not enough settings: " << std::endl;
      std::cout << "[usage] : " << std::endl;
      std::cout << "    ./unfoldDataUncertainties_noempty_pp -c binning.config -r 4 -n 10 -p 1 " << std::endl;
      return 1;
    }
  
  return unfoldDataUncertainties_noempty_pp(config, niterations, cone_size, primer);
  
}
