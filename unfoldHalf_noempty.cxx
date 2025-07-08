#include <iostream>

using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int unfoldHalf_noempty(const std::string configfile = "binning.config", const int niterations = 10, const int cone_size = 4)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());  
  std::string j10_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_10_new_ProdA_2024-00000021.root";
  std::string j20_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_20_new_ProdA_2024-00000021.root";
  std::string j30_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v6_30_new_ProdA_2024-00000021.root";

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float nrecojets[3];
  float dphi_reco[3];
  float match[3];
  float mbd_vertex[3];


  TFile *fin[3];
  fin[0] = new TFile(j10_file.c_str(), "r");
  fin[1] = new TFile(j20_file.c_str(), "r");
  fin[2] = new TFile(j30_file.c_str(), "r");

  TNtuple *tn[3];
  for (int i = 0 ; i < 3; i++)
    {
      if (!fin[i])
	{
	  std::cout << "no file: " << i << std::endl;
	  return 0;

	}

      tn[i] = (TNtuple*) fin[i]->Get("tn_match");
      tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
      tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
      tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
      tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
      tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
      tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
      tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
      tn[i]->SetBranchAddress("nrecojets", &nrecojets[i]);
      tn[i]->SetBranchAddress("matched", &match[i]);
      tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
    }
  float cs_10 = (2.889e-6);
  float cs_20 = 5.4067742e-8;
  float cs_30 = (2.505e-9);
  float scale_factor[3];
  scale_factor[0] = cs_10/cs_30;
  scale_factor[1] = cs_20/cs_30; 
  //scale_factor[1] = (3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  
  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();
  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  float width = 0.1 + JER_sys;
  fgaus->SetParameters(1, 0, width);

  if (JER_sys != 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
    }

  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
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

  Int_t njet_sys = rb.get_njet_sys();
  std::cout << "NJET scaling : " << njet_sys << std::endl;
  double njet_scale[10] = {0, 
			   0, 
			   1.22359, 
			   0.838055, 
			   0.669555, 
			   0.513105, 
			   0.398124, 
			   0.486207, 
			   0.263404, 
			   2.91882};

  Int_t prior_sys = rb.get_prior_sys();
  int prior_iteration = 0;
  
  // Vertex Reweight


  std::vector<std::pair<float, float>> vertex_scales;

  Int_t vtx_sys = rb.get_vtx_sys();
  if (primer != 1)
    {

      TFile *fvtx = new TFile(Form("%s/vertex/vertex_reweight_r%02d.root", rb.get_code_location().c_str(), cone_size),"r");
      TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");
      for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
	{
	  vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
	}
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

  float sample_boundary[4] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;

  TString responsepath = rb.get_code_location() + "/response_matrices/response_matrix_r0" + std::to_string(cone_size) + "_HALF.root";
  
  TFile *fresponse = new TFile(responsepath.Data(),"r");
  
  RooUnfoldResponse *rooResponse = (RooUnfoldResponse*) fresponse->Get("response_noempty");
  if (!rooResponse)
    {
      std::cout << "no repsonse" << std::endl;
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
  TH1D *h_flat_truth_mapping = (TH1D*) fresponse->Get("h_flat_truth_mapping"); 
  if (!h_flat_truth_mapping)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_reco_mapping = (TH1D*) fresponse->Get("h_flat_reco_mapping"); 
  if (!h_flat_reco_mapping)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }
  TH1D *h_flat_truth_skim = (TH1D*) fresponse->Get("h_flat_truth_skim"); 
  if (!h_flat_truth_skim)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_reco_skim = (TH1D*) fresponse->Get("h_flat_reco_skim"); 
  if (!h_flat_reco_skim)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }

  
  for (int isample = 0; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries1 = tn[isample]->GetEntries()/2;
      int entries2 = tn[isample]->GetEntries();
      for (int i = entries1; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  int inrecojets = nrecojets[isample];
	  double event_scale = scale_factor[isample];

	  float max_truth = 0;
	  float max_reco = 0;
	  float min_truth = 0;
	  float min_reco = 0;

	  if (pt1_truth[isample] >= pt2_truth[isample])
	    {
	      max_truth = pt1_truth[isample];
	      max_reco = pt1_reco[isample];
	      min_truth = pt2_truth[isample];
	      min_reco = pt2_reco[isample];
	    }
	  else
	    {
	      max_truth = pt2_truth[isample];
	      max_reco = pt2_reco[isample];
	      min_truth = pt1_truth[isample];
	      min_reco = pt1_reco[isample];
	    }
	    
	  float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;
	  float pt1_reco_bin = nbins;
	  float pt2_reco_bin = nbins;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;

	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;

	  double smear1 = fgaus->GetRandom();
	  double smear2 = fgaus->GetRandom();

	  if (JES_sys != 0)
	    {
	      es1 = es1 + (JES_sys + smear1)*e1;
	      es2 = es2 + (JES_sys + smear2)*e2; 
	    }
	  else
	    {
	      es1 = es1 + smear1*e1;
	      es2 = es2 + smear2*e2; 	  
	    }
	  
	  float maxi = std::max(es1, es2);
	  float mini = std::min(es1, es2);


	  for (int ib = 0; ib < nbins; ib++)
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
	  
	  bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] > dphicut);
	  if (reco_good)
	    {
	      h_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      h_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);

	      //h_count_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
	      //h_count_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
	      continue;
	    }
	}
    }

  h_flat_data_pt1pt2->Scale(.5);
  TH1D *h_flat_data_skim = (TH1D*) h_flat_reco_skim->Clone();
  h_flat_data_skim->SetName("h_flat_data_skim");
  h_flat_data_skim->Reset();
  histo_opps::skim_down_histo(h_flat_data_skim, h_flat_data_pt1pt2, h_flat_reco_mapping); 
  
  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];

  for (int iter = 0; iter < niterations; iter++ )
    {
      
      RooUnfoldBayes   unfold (rooResponse, h_flat_data_skim, iter + 1);    // OR
      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
      h_flat_unfold_pt1pt2[iter]->Reset();
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
      histo_opps::fill_up_histo(h_flat_unfold_skim[iter], h_flat_unfold_pt1pt2[iter], h_flat_truth_mapping);
    }

  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_r0" + std::to_string(cone_size) + "_HALF.root";
    
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");
  h_flat_data_skim->Write();
  h_flat_data_pt1pt2->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_unfold_pt1pt2[iter]->Write();
    }
  fout->Close();

  return 0;
}

