#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int makeTruth_dist(const std::string configfile = "binning.config", const int cone_size = 4)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  read_binning rb(configfile.c_str());  
  std::string j10_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" +std::to_string(cone_size) + "_v6_10_new_ProdA_2024-00000021.root";
  std::string j20_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" +std::to_string(cone_size) + "_v6_20_new_ProdA_2024-00000021.root";
  std::string j30_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" +std::to_string(cone_size) + "_v6_30_new_ProdA_2024-00000021.root";


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
  scale_factor[0] = (9987000./9988000)*cs_10/cs_30;
  scale_factor[1] = (9987000./9954000.)*cs_20/cs_30; 
  //scale_factor[1] = (3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  

  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

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

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
  std::cout << " measure_cut " << measure_subleading_cut << std::endl;
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  float measure_values[10] = {0};
  
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_values[ir] = ipt_bins[rb.get_measure_region(ir)];
    }

  TH1D *h_truth_lead_sample[3];
  for (int i = 0; i < 3; i++)
    {
      h_truth_lead_sample[i]= new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

  TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  // pure fills
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);

  TH1D *h_linear_truth_xj[3];
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i] = new TH1D(Form("h_linear_truth_xj_%d", i),";x_{J};;", 10, 0, 1.0);
    }
  
  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  for (int isample = 0; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries1 = tn[isample]->GetEntries()/2;
      int entries2 = tn[isample]->GetEntries();
      for (int i = 0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  int inrecojets = nrecojets[isample];
	  double event_scale = scale_factor[isample];

	  h_truth_lead_sample[isample]->Fill(pt1_truth[isample], event_scale);

	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;

	  float max_truth = 0;
	  float min_truth = 0;
	  
	  if (pt1_truth[isample] >= pt2_truth[isample])
	    {
	      max_truth = pt1_truth[isample];
	      min_truth = pt2_truth[isample];
	    }
	  else
	    {
	      max_truth = pt2_truth[isample];
	      min_truth = pt1_truth[isample];
	    }
	    

	  float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  
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
	    }
	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && dphi_truth[isample] > dphicuttruth);

	  
	  if (truth_good)
	    {
	      for (int im = 0; im < mbins; im++)
		{
		  if (e1 >= measure_values[im] && e1 < measure_values[im+1] && e2 >= measure_subleading_cut)
		    {
		      h_linear_truth_xj[im]->Fill(e2/e1, event_scale);
		    }
		}
	      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
	      h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
	    }
	}
    }


  TString responsepath = rb.get_code_location() + "/truth_hists/truth_hist_r0" + std::to_string(cone_size)+".root";
  
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i]->Write();
    }

  h_flat_truth_to_response_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_count_flat_truth_pt1pt2->Write();
  fr->Close();

  return 0;
}

