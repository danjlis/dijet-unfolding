#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif
#include "dlUtility.h"
#include "read_binning.h"


int createResponse_wfakes_full(const int niterations = 10)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  std::string mb_file = "TREE_MATCH_v4_8_new_ProdA_2024-00000021.root";
  std::string j10_file = "TREE_MATCH_v6_10_new_ProdA_2024-00000021.root";
  std::string j30_file = "TREE_MATCH_v6_30_new_ProdA_2024-00000021.root";

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float dphi_reco[3];
  float match[3];

  float sample_boundary_goal[4];
  sample_boundary_goal[0] = 3;
  sample_boundary_goal[1] = 13;
  sample_boundary_goal[2] = 29.99;
  sample_boundary_goal[3] = 100;

  float sample_boundary[4];
  sample_boundary[0] = 0;
  sample_boundary[1] = 0;
  sample_boundary[2] = 0;
  sample_boundary[3] = 100;


  TFile *fin[3];
  fin[0] = new TFile(mb_file.c_str(), "r");
  fin[1] = new TFile(j10_file.c_str(), "r");
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
      tn[i]->SetBranchAddress("matched", &match[i]);
    }
  float scale_factor[3];
  scale_factor[0] = (4197e-2)/(2.505e-9);//
  scale_factor[1] = (9987000./10000000.)*(3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  
  float resolution_cut = 4*0.16;
  //float bin_edges[11];


  read_binning rb("binning.config");

  Int_t read_nbins = rb.get_nbins();
  Int_t read_bbin = rb.get_bbins();
  Double_t read_fixed = rb.get_fixed();
  Double_t read_minimum = rb.get_minimum();

  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  float alpha = TMath::Power(read_fixed/read_minimum, 1/(float)read_bbin);
  float bin_min = read_minimum;
  float bin_max = read_minimum*TMath::Power(alpha, nbins);

  float bin_xj10 = 1.0;
  float bin_xj1 = 1.0*(bin_min/bin_max);

  std::cout << "Alpha = " << alpha <<std::endl;

  int measure_leading_bin = 0;
  int measure_subleading_bin = 0;
  int truth_leading_bin = 0;
  int truth_subleading_bin = 0;
  float truth_leading_cut = 0;
  float reco_leading_cut = 0;
  float reco_meas_leading_cut = 0;
  float truth_subleading_cut = 0;
  float reco_subleading_cut = 0;
  float reco_meas_subleading_cut = 0;

  float measure_leading_goal = 20;

  float truth_subleading_goal = 0;
  float truth_leading_goal = 13;

  float dphicut = 3*TMath::Pi()/4.;
  float dphicuttruth = 3*TMath::Pi()/4.;

  // Make binning
  for (int i = 0; i < nbins+1; i++)
    {
      float ipt = bin_min*TMath::Power(alpha, (float)i);
      float ixj = bin_xj1*TMath::Power(alpha, (float)i);
      ipt_bins[i] = ipt;
      ixj_bins[i] = ixj;
      for (int ib = 0; ib < 4; ib++)
	{
	  if (sample_boundary[ib] == 0 && ipt >= sample_boundary_goal[ib])
	    {
	      std::cout << "Found " << ib << std::endl;
	      sample_boundary[ib] = ipt;
	    }
	}
      if (measure_leading_bin == 0 && ipt >= measure_leading_goal) measure_leading_bin = i;
      if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
      if (truth_subleading_bin == 0 && ipt >= truth_subleading_goal) truth_subleading_bin = i;
      std::cout << i << " : " <<  ipt << " -- " << ixj <<  std::endl;
    }
  truth_subleading_bin = 0;
  truth_leading_cut=ipt_bins[truth_leading_bin];
  truth_subleading_cut=ipt_bins[truth_subleading_bin];

  reco_leading_cut=ipt_bins[truth_leading_bin + 2];
  reco_subleading_cut=ipt_bins[truth_subleading_bin + 2];

  reco_meas_leading_cut=ipt_bins[truth_leading_bin + 3];
  reco_meas_subleading_cut=ipt_bins[truth_subleading_bin + 3];
  measure_subleading_bin = truth_subleading_bin + 3;

  for (int ib = 0; ib < 4; ib++)
    {
      std::cout <<  sample_boundary[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  reco_meas_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  reco_meas_subleading_cut << std::endl;

  TH1D *h_truth_lead_sample[3];
  for (int i = 0; i < 3; i++)
    {
      h_truth_lead_sample[i]= new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

  TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_reco_lead = new TH1D("h_match_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_reco_sublead = new TH1D("h_match_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

  // pure fills
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

  TH1D *h_truth_dphi = new TH1D("h_truth_dphi",";#Delta#phi;1/N", 64, TMath::Pi()/2, TMath::Pi());
  TH1D *h_reco_dphi = new TH1D("h_reco_dphi",";#Delta#phi;1/N", 64, TMath::Pi()/2, TMath::Pi());
  TH1D *h_truth_match_dphi = new TH1D("h_truth_match_dphi",";#Delta#phi;1/N", 64, TMath::Pi()/2, TMath::Pi());
  TH1D *h_reco_match_dphi = new TH1D("h_reco_match_dphi",";#Delta#phi;1/N", 64, TMath::Pi()/2, TMath::Pi());

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  for (int isample = 1; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries1 = tn[isample]->GetEntries()/2;
      int entries2 = tn[isample]->GetEntries();
      for (int i = 0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  h_truth_lead_sample[isample]->Fill(pt1_truth[isample], scale_factor[isample]);
	  
	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;


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
	    
	  float maxi = std::max(pt1_reco[isample], pt2_reco[isample]);
	  float mini = std::min(pt1_reco[isample], pt2_reco[isample]);

	  float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;
	  float pt1_reco_bin = nbins;
	  float pt2_reco_bin = nbins;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;

	  if (e1 > truth_leading_cut) h_truth_lead->Fill(e1, scale_factor[isample]);
	  if (e2 >  truth_subleading_cut) h_truth_sublead->Fill(e2, scale_factor[isample]);

	  if (maxi >  reco_leading_cut) h_reco_lead->Fill(maxi, scale_factor[isample]);
	  if (mini >  reco_subleading_cut) h_reco_sublead->Fill(mini, scale_factor[isample]);


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
	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && dphi_truth[isample] > dphicuttruth);
	  //std::cout << 
	  bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] > dphicut);

	  h_truth_dphi->Fill(dphi_truth[isample], scale_factor[isample]);
	  h_reco_dphi->Fill(dphi_reco[isample], scale_factor[isample]);

	  if (!truth_good && !reco_good) continue;
	  //RooUnfoldResponse rooResponse_histo(*h_flat_reco_to_response_pt1pt2, *h_flat_truth_to_response, *flat_response_pt1pt2);
	  if (truth_good && !reco_good)
	    {

	      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      h_truth_xj->Fill(e2/e1, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      continue;
	    }
	  if (!truth_good && reco_good)
	    {
	      h_pt1pt2->Fill(es1, es2, scale_factor[isample]);
	      h_pt1pt2->Fill(es2, es1, scale_factor[isample]);

	      h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_reco_xj->Fill(mini/maxi, scale_factor[isample]);

	      h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      std::cout << " Fake 1: Truth bad Reco good with : " << dphicuttruth  << std::endl;

	      std::cout << e1 << " - " << e2 << " - " << dphi_truth[isample] << std::endl;

	      tn[isample]->Show();
	      rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      continue;
	    }
	  
	  if (reco_good && truth_good)
	    {
	      if (e1 > truth_leading_cut) h_match_truth_lead->Fill(e1, scale_factor[isample]);
	      if (e2 >  truth_subleading_cut) h_match_truth_sublead->Fill(e2, scale_factor[isample]);
	      if (maxi >  reco_leading_cut) h_match_reco_lead->Fill(maxi, scale_factor[isample]);
	      if (mini >  reco_subleading_cut) h_match_reco_sublead->Fill(mini, scale_factor[isample]);

	      h_truth_match_dphi->Fill(dphi_truth[isample], scale_factor[isample]);
	      h_reco_match_dphi->Fill(dphi_reco[isample], scale_factor[isample]);

	      
	      h_pt1pt2->Fill(es1, es2, scale_factor[isample]);
	      h_pt1pt2->Fill(es2, es1, scale_factor[isample]);
	      h_truth_xj->Fill(e2/e1, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);

	      h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_reco_xj->Fill(mini/maxi, scale_factor[isample]);

	      h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      continue;
	    }
	  else
	    {
	      h_truth_xj->Fill(e2/e1, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);

	      h_pt1pt2->Fill(es1, es2, scale_factor[isample]);
	      h_pt1pt2->Fill(es2, es1, scale_factor[isample]);

	      h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_reco_xj->Fill(mini/maxi, scale_factor[isample]);

	      rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      std::cout << " Fake 2: Truth good Reco good" << std::endl;
	      tn[isample]->Show();

	      h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, scale_factor[isample]);
	      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, scale_factor[isample]);
	      continue;
	    }

	  std::cout << "something not right" << std::endl;
	}
    }

  h_flat_reco_pt1pt2->Scale(.5, "width");
  h_flat_truth_pt1pt2->Scale(.5, "width");

  // int ntoys  = 100;
  // for (int itoy = 0; itoy < ntoys ; itoy++)
  //   {
  //     for (int ibin = 0; ibin < nbins; ibin ++)
  // 	{

  // 	  RooUnfoldResponse rooResponse_histo(h_flat_reco_to_response_pt1pt2, h_flat_truth_to_response_pt1pt2, h_flat_response_pt1pt2);

  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 3;
  for (int iter = 0; iter < niterations; iter++ )
    {
      
      RooUnfoldBayes   unfold (&rooResponse, h_flat_reco_pt1pt2, iter + 1);    // OR

      h_flat_unfold_pt1pt2[iter]= (TH1D*) unfold.Hunfold();
    }

  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  for (int iter = 0; iter < niterations; iter++)
    {
      dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
    }
  
  dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
  dlutility::SetLineAtt(h_flat_reco_pt1pt2, kBlue, 2, 1);

  h_flat_truth_pt1pt2->Draw("hist");
  h_flat_reco_pt1pt2->Draw("hist same");
  h_flat_unfold_pt1pt2[niter]->Draw("same p");

  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  TH1D *h_subleading_reco[nbins];
  TH1D *h_subleading_truth[nbins];
  TH1D *h_subleading_unfold[nbins];
  for (int ib = 0; ib; ib++)
    {
      h_subleading_unfold[ib]= new TH1D(Form("h_subleading_unfold_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
      h_subleading_truth[ib]= new TH1D(Form("h_subleading_truth_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
      h_subleading_reco[ib]= new TH1D(Form("h_subleading_reco_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
    }

  
  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_truth = new TH1D("h_xjunc_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xjunc_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xjunc_unfold[iter] = new TH1D(Form("h_xjunc_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int xbin = ib/nbins;
      int ybin = ib%nbins;
      
      int b = h_pt1pt2_reco->GetBin(xbin+1, ybin+1);

      h_pt1pt2_reco->SetBinContent(b, h_flat_reco_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_truth->SetBinContent(b, h_flat_truth_pt1pt2->GetBinContent(ib+1));
      for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinContent(b, h_flat_unfold_pt1pt2[iter]->GetBinContent(ib+1));
      h_pt1pt2_reco->SetBinError(b, h_flat_reco_pt1pt2->GetBinError(ib+1));
      h_pt1pt2_truth->SetBinError(b, h_flat_truth_pt1pt2->GetBinError(ib+1));
      for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinError(b, h_flat_unfold_pt1pt2[iter]->GetBinError(ib+1));
    }
  h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");
      
  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
  cpt1pt2->Divide(3, 1);
  cpt1pt2->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  h_pt1pt2_reco->Draw("colz");
  dlutility::DrawSPHENIXpp(0.22, 0.85);
  
  cpt1pt2->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_truth->Draw("colz");
  cpt1pt2->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_unfold[niter]->Draw("colz");


  cpt1pt2->Print("pt1pt2.pdf");
  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
	{
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);

	  if (ix > iy)
	    {

	      h_pt1pt2_reco->SetBinContent(bin, h_pt1pt2_reco->GetBinContent(bin)*2.);
	      h_pt1pt2_truth->SetBinContent(bin, h_pt1pt2_truth->GetBinContent(bin)*2.);
	      for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinContent(bin, h_pt1pt2_unfold[iter]->GetBinContent(bin)*2.);
	      
	      h_pt1pt2_reco->SetBinError(bin, h_pt1pt2_reco->GetBinError(bin)/sqrt(2));
	      h_pt1pt2_truth->SetBinError(bin, h_pt1pt2_truth->GetBinError(bin)/sqrt(2));
	      for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinError(bin, h_pt1pt2_unfold[iter]->GetBinError(bin)*sqrt(2));

	    }
	  else if (ix < iy)
	    {	      
	      h_pt1pt2_reco->SetBinContent(bin, 0);
	      h_pt1pt2_truth->SetBinContent(bin, 0);
	      for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinContent(bin, 0);
	    }

	}
    }

  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy <= ix; iy++)
	{
	  int low =  iy - ix - 1;
	  int high = iy - ix + 1;

	  int xjbin_low = nbins + low + 1;
	  int xjbin_high = nbins + low + 2;
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);
	  std::cout << ix << " -- " << iy << " --> " << low << " --> " << xjbin_high << "( " << h_xj_reco->GetBinCenter(xjbin_high)<<" ) " << " + " << xjbin_low << "( " << h_xj_reco->GetBinCenter(xjbin_low)<<" ) "<< std::endl;

	  if (ix < measure_leading_bin) continue;
	  if (iy < measure_subleading_bin) continue;
	  if (ix > nbins - 2) continue;
	  if (iy > nbins - 2) continue;

	  if (ix == iy)
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	    }

	  if (ix == iy)
	    {
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low), h_pt1pt2_truth->GetBinContent(bin));
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_truth->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_high), h_pt1pt2_truth->GetBinContent(bin)/2.);
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low), h_pt1pt2_truth->GetBinContent(bin)/2.);
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_truth->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_truth->GetBinError(bin)/sqrt(2), 2));
	    }
	  for (int iter = 0; iter < niterations; iter++ )
	    {
	      if (ix == iy)
		{
		  h_xj_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_low), h_pt1pt2_unfold[iter]->GetBinContent(bin));
		  h_xjunc_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_unfold[iter]->GetBinError(bin)/sqrt(2), 2));
		}
	      else
		{
		  h_xj_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_high), h_pt1pt2_unfold[iter]->GetBinContent(bin)/2.);
		  h_xj_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_low), h_pt1pt2_unfold[iter]->GetBinContent(bin)/2.);
		  h_xjunc_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_unfold[iter]->GetBinError(bin)/sqrt(2), 2));
		  h_xjunc_unfold[iter]->Fill(h_xj_unfold[iter]->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_unfold[iter]->GetBinError(bin)/sqrt(2), 2));
		}
	    }	  
	}
    }
  for (int i = 1; i <= nbins; i++)
    {
      h_xj_reco->SetBinError(i, sqrt(h_xjunc_reco->GetBinError(i)));
      h_xj_truth->SetBinError(i, sqrt(h_xjunc_truth->GetBinError(i)));
      for (int iter = 0; iter < niterations; iter++) h_xj_unfold[iter]->SetBinError(i, sqrt(h_xjunc_unfold[iter]->GetBinError(i)));
    }


  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_truth_xj, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_truth_xj, kRed, 1, 24);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 1, 24);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);
  dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter]->Scale(1./h_xj_unfold[iter]->Integral(0, -1), "width");
    }
  h_xj_truth->Scale(1./h_xj_truth->Integral(0, -1), "width");
  h_xj_reco->Scale(1./h_xj_reco->Integral(0, -1), "width");
  h_truth_xj->Scale(1./h_truth_xj->Integral(0, -1), "width");
  h_reco_xj->Scale(1./h_reco_xj->Integral(0, -1), "width");
  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_truth->Draw("p");
  h_xj_unfold[niter]->Draw("same p");
  h_xj_reco->Draw("hist same");
  h_xj_truth->Draw("same hist");
  h_xj_reco->Draw("p same");
  h_xj_unfold[niter]->Draw("same hist");
  h_xj_unfold[niter]->Draw("same p");
  dlutility::DrawSPHENIXpp(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
  TLegend *leg = new TLegend(0.22, 0.3, 0.4, 0.55);
  leg->SetLineWidth(0);
  leg->AddEntry(h_xj_reco, "Pythia Reco");
  leg->AddEntry(h_xj_truth, "Pythia8 Truth");
  leg->AddEntry(h_xj_unfold[niter], "Unfolded");
  leg->Draw("same");
  dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
    
  cxj->cd(2);

  TH1D *h_reco_compare = (TH1D*) h_xj_unfold[niter]->Clone();
  h_reco_compare->Divide(h_xj_truth);
  h_reco_compare->SetTitle(";x_{J}; Unfold / Truth");
  dlutility::SetFont(h_reco_compare, 42, 0.06);
  dlutility::SetLineAtt(h_reco_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_reco_compare, kBlack, 1,8);

 
  h_reco_compare->Draw("p");
  TLine *line = new TLine(0.1, 1, 1, 1);
  line->SetLineStyle(4);
  line->SetLineColor(kRed + 3);
  line->SetLineWidth(2);
  line->Draw("same");

  TCanvas *c_iter = new TCanvas("c_iter", "c_iter");
  TH1D *h_closure[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_closure[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
      h_closure[iter]->SetName(Form("h_closure_%d", iter));
		
      h_closure[iter]->Divide(h_xj_truth);
      h_closure[iter]->SetTitle(";x_{J}; Unfold / Truth");
      dlutility::SetFont(h_closure[iter], 42, 0.06);
    }
  int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 

  for (int iter = 0; iter < 5; iter++)
    {
      dlutility::SetLineAtt(h_closure[1 + 2*iter], colors[iter], 1,1);
      dlutility::SetMarkerAtt(h_closure[1 + 2*iter], colors[iter], 1,8);
    }
  h_closure[1]->Draw("p");
  h_closure[3]->Draw("p same");
  h_closure[5]->Draw("p same");
  h_closure[7]->Draw("p same");
  h_closure[9]->Draw("p same");
  TLine *line2 = new TLine(0.1, 1, 1, 1);
  line2->SetLineStyle(4);
  line2->SetLineColor(kRed + 3);
  line2->SetLineWidth(2);
  line2->Draw("same");


  TCanvas *cjet = new TCanvas("cjet","cjet", 700, 500);
  dlutility::createCutCanvas(cjet);
  cjet->cd(1);
  gPad->SetLogy();

  dlutility::SetLineAtt(h_truth_lead, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_truth_lead, kRed, 1, 24);

  dlutility::SetLineAtt(h_truth_sublead, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_truth_sublead, kBlue, 1, 24);

  dlutility::SetLineAtt(h_reco_lead, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_reco_lead, kRed, 1, 20);

  dlutility::SetLineAtt(h_reco_sublead, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_reco_sublead, kBlue, 1, 20);

  h_truth_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
    
  h_truth_lead->Draw();
  h_truth_sublead->Draw("same");
  h_reco_lead->Draw("same");
  h_reco_sublead->Draw("same");
  cjet->cd(2);
  dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.1, 0.74, 0, kBlack, 0.1);
  dlutility::drawText("All Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);

  TLegend *leg3 = new TLegend(0.01, 0.13, 0.7, 0.62);
  leg3->SetLineWidth(0);
  leg3->AddEntry(h_reco_lead, "Leading Reco");
  leg3->AddEntry(h_reco_sublead, "Subleading Reco");
  leg3->AddEntry(h_truth_lead, "Leading Truth");
  leg3->AddEntry(h_truth_sublead, "Subleading Truth");
  leg3->SetTextSize(0.1);
  leg3->Draw("same");

  TCanvas *cmjet = new TCanvas("cmjet","cmjet", 700, 500);
  dlutility::createCutCanvas(cmjet);
  cmjet->cd(1);
  gPad->SetLogy();

  dlutility::SetLineAtt(h_match_truth_lead, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_match_truth_lead, kRed, 1, 24);

  dlutility::SetLineAtt(h_match_truth_sublead, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_match_truth_sublead, kBlue, 1, 24);

  dlutility::SetLineAtt(h_match_reco_lead, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_match_reco_lead, kRed, 1, 20);

  dlutility::SetLineAtt(h_match_reco_sublead, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_match_reco_sublead, kBlue, 1, 20);

  h_match_reco_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
    
  h_match_reco_lead->Draw();
  h_match_truth_sublead->Draw("same");
  h_match_truth_lead->Draw("same");
  h_match_reco_sublead->Draw("same");
  cmjet->cd(2);
  dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.1, 0.74, 0, kBlack, 0.1);
  dlutility::drawText("Matched Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);
  TLegend *lef3 = new TLegend(0.01, 0.13, 0.7, 0.62);
  lef3->SetLineWidth(0);
  lef3->AddEntry(h_match_reco_lead, "Leading Reco");
  lef3->AddEntry(h_match_reco_sublead, "Subleading Reco");
  lef3->AddEntry(h_match_truth_lead, "Leading Truth");
  lef3->AddEntry(h_match_truth_sublead, "Subleading Truth");
  lef3->SetTextSize(0.1);
  lef3->Draw("same");


  TCanvas *cresponse = new TCanvas("fd","fd", 500, 500);
  h_flat_response_pt1pt2->Draw("colz");


  TFile *fr = new TFile("response_matrix.root","recreate");
  rooResponse.Write();

  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_flat_response_pt1pt2->Write();
  h_truth_dphi->Write();
  h_truth_match_dphi->Write();
  h_reco_dphi->Write();
  h_reco_match_dphi->Write();
  fr->Close();
  

  return 0;
  
}
