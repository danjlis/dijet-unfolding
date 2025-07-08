#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif
#include "../macros/dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int createResponse_half_noempty(const std::string configfile = "binning.config", const int niterations = 10)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  
  std::string j10_file = "TREE_MATCH_v6_10_new_ProdA_2024-00000021.root";
  std::string j20_file = "TREE_MATCH_v6_20_new_ProdA_2024-00000021.root";
  std::string j30_file = "TREE_MATCH_v6_30_new_ProdA_2024-00000021.root";

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float nrecojets[3];
  float dphi_reco[3];
  float match[3];


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
    }
  float cs_10 = (3.646e-6);
  float cs_20 = 5.4067742e-8;
  float cs_30 = (2.505e-9);
  float scale_factor[3];
  scale_factor[0] = (9987000./9988000)*cs_10/cs_30;
  scale_factor[1] = (9987000./9954000.)*cs_20/cs_30; 
  //scale_factor[1] = (3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  
  read_binning rb(configfile.c_str());

  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);
  fgaus->SetParameters(1, 0, 0.1);
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  if (JER_sys > 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
      fgaus->SetParameters(1, 0, JER_sys);
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

  TFile *fvtx = new TFile("vertex_reweight.root","r");
  TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");

  std::vector<std::pair<float, float>> vertex_scales;

  for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
    {
      vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
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

  TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  for (int isample = 0; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries_half[3] = {0, static_cast<int>(tn[isample]->GetEntries()/2), static_cast<int>(tn[isample]->GetEntries())};
      for (int ihalf = 0; ihalf < 2; ihalf++)
	{
	  for (int i = entries_half[ihalf]; i < entries_half[ihalf + 1]; i++)
	    {
	      tn[isample]->GetEntry(i);
	      int inrecojets = nrecojets[isample];
	      double event_scale = scale_factor[isample];
	      if (njet_sys)
		{
		  event_scale *= njet_scale[(inrecojets > 9 ? 0 : inrecojets)];
		}
	      // Vertex Rewieghting
	      for (int ib = 0; ib < vertex_scales.size(); ib++)
		{
		  if (mbd_vertex[isample] < vertex_scales.at(ib).first)
		    {
		      event_scale *= vertex_scales.at(ib).second;
		      if (i < 100) std:cout << "found z = " << vertex_scales.at(ib).first << " " <<vertex_scales.at(ib).second<<std::endl;
		      break;
		    }
		}

	      h_truth_lead_sample[isample]->Fill(pt1_truth[isample], event_scale);
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
	    

	      float pt1_truth_bin = nbins;
	      float pt2_truth_bin = nbins;
	      float pt1_reco_bin = nbins;
	      float pt2_reco_bin = nbins;

	      float e1 = max_truth;
	      float e2 = min_truth;
	      float es1 = max_reco;
	      float es2 = min_reco;

	      if (JES_sys != 0)
		{
		  es1 = es1 + JES_sys*e1;
		  es2 = es2 + JES_sys*e2; 
		}
	      else if (JER_sys > 0)
		{
		  double smear1 = fgaus->GetRandom();
		  double smear2 = fgaus->GetRandom();
		  es1 = es1 + smear1*e1;
		  es2 = es2 + smear2*e2; 
		}

	      float maxi = std::max(es1, es2);
	      float mini = std::min(es1, es2);
	      if (e1 > truth_leading_cut) h_truth_lead->Fill(e1, event_scale);
	      if (e2 >  truth_subleading_cut) h_truth_sublead->Fill(e2, event_scale);

	      if (maxi >  reco_leading_cut) h_reco_lead->Fill(maxi, event_scale);
	      if (mini >  reco_subleading_cut) h_reco_sublead->Fill(mini, event_scale);


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

	      bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] > dphicut);

	      if (!truth_good && !reco_good) continue;
	      //RooUnfoldResponse rooResponse_histo(*h_flat_reco_to_response_pt1pt2, *h_flat_truth_to_response, *flat_response_pt1pt2);
	      if (truth_good && !reco_good)
		{

		  if (e1 >= measure_leading_cut && e2 >= measure_subleading_cut)
		    {
		      h_truth_xj->Fill(e2/e1, event_scale);
		      h_linear_truth_xj->Fill(e2/e1, event_scale);
		    }

		  if (ihalf)
		    {
		      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		    }
		  else
		    {
		      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		  
		      rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		    }
		  continue;
		}
	      // if (!truth_good && reco_good)
	      //   {
	      //     h_pt1pt2->Fill(es1, es2, event_scale);
	      //     h_pt1pt2->Fill(es2, es1, event_scale);

	      //     h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
	      //     if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
	      // 	{
	      // 	  h_reco_xj->Fill(mini/maxi, event_scale);
	      // 	  h_linear_reco_xj->Fill(mini/maxi, event_scale);
	      // 	}
	      //     h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
	      //     rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
	      //     continue;
	      //   }
	  
	      if (match[isample] && reco_good && truth_good)
		{
		  if (e1 > truth_leading_cut) h_match_truth_lead->Fill(e1, event_scale);
		  if (e2 >  truth_subleading_cut) h_match_truth_sublead->Fill(e2, event_scale);
		  if (maxi >  reco_leading_cut) h_match_reco_lead->Fill(maxi, event_scale);
		  if (mini >  reco_subleading_cut) h_match_reco_sublead->Fill(mini, event_scale);


		  h_pt1pt2->Fill(es1, es2, event_scale);
		  h_pt1pt2->Fill(es2, es1, event_scale);
		  if (e1 >= measure_leading_cut && e2 >= measure_subleading_cut)
		    {
		      h_truth_xj->Fill(e2/e1, event_scale);
		      h_linear_truth_xj->Fill(e2/e1, event_scale);
		    }
		  if (ihalf)
		    {
		      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);

		      h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
		      h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
		      if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
			{
			  h_reco_xj->Fill(mini/maxi, event_scale);
			  h_linear_reco_xj->Fill(mini/maxi, event_scale);
			}
		    }
		  else
		    {
		      h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
		      h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
		      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		      h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		      rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
		      rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
		    }
		  continue;
		}
	      // else
	      //   {
	      //     h_truth_xj->Fill(e2/e1, event_scale);
	      //     h_linear_truth_xj->Fill(e2/e1, event_scale);
	      //     h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      //     h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);

	      //     h_pt1pt2->Fill(es1, es2, event_scale);
	      //     h_pt1pt2->Fill(es2, es1, event_scale);

	      //     h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
	      //     h_reco_xj->Fill(mini/maxi, event_scale);
	      //     h_linear_reco_xj->Fill(mini/maxi, event_scale);
	      //     rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      //     rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      //     rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);

	      //     h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
	      //     h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
	      //     h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      //     h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      //     continue;
	      //   }

	      //std::cout << "something not right" << std::endl;
	    }
	}
    }

  int nbinsx = h_flat_response_pt1pt2->GetXaxis()->GetNbins();
  int nbinsy = h_flat_response_pt1pt2->GetYaxis()->GetNbins();

  TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_pt1pt2->Clone();
  TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_pt1pt2->Clone();

  h_flat_truth_mapping->SetName("h_flat_truth_mapping");
  h_flat_reco_mapping->SetName("h_flat_reco_mapping");
  h_flat_reco_mapping->Reset();
  h_flat_truth_mapping->Reset();
  
  int nempty_reco = 0;
  std::vector<int> binnumbers_reco{};
  int nempty_truth = 0;
  std::vector<int> binnumbers_truth{};
  int nrecobins = 0;
  int ntruthbins = 0;
  for (int ib = 0; ib < nbins*nbins; ib++)
    {
      if (h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1) == 0)
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
      if (h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1) == 0)
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
  TH1D *h_flat_truth_skim = new TH1D("h_flat_truth_skim","", ntruthbins, 0, ntruthbins);
  TH1D *h_flat_truth_measure_skim = new TH1D("h_flat_truth_measure_skim","", ntruthbins, 0, ntruthbins);
  TH1D *h_flat_reco_skim = new TH1D("h_flat_reco_skim","", nrecobins, 0, nrecobins);
  TH1D *h_flat_reco_measure_skim = new TH1D("h_flat_reco_measure_skim","", nrecobins, 0, nrecobins);
  TH2D *h_flat_response_skim = new TH2D("h_flat_response_skim","", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

  histo_opps::skim_down_histo(h_flat_reco_measure_skim, h_flat_reco_pt1pt2, h_flat_reco_mapping);
  
  for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int bintruth = binnumbers_truth.at(ib);      
      int binreco = binnumbers_reco.at(ib);
      if (binreco)
	{
	  h_flat_reco_skim->SetBinContent(binreco, h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1));
	  h_flat_reco_skim->SetBinError(binreco, h_flat_reco_to_response_pt1pt2->GetBinError(ib+1));
	}
      if (bintruth)
	{
	 
	  h_flat_truth_skim->SetBinContent(bintruth, h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1));
	  h_flat_truth_skim->SetBinError(bintruth, h_flat_truth_to_response_pt1pt2->GetBinError(ib+1));

	  for (int ibr = 0; ibr < nbins*nbins; ibr++)
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
  std::cout <<" Nbins skim reco = "<<h_flat_truth_skim->GetNbinsX()<<std::endl;
  std::cout <<" Nbins skim reco = "<<h_flat_response_skim->GetXaxis()->GetNbins()<<"  --  " << h_flat_response_skim->GetYaxis()->GetNbins()<<std::endl; 
  RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);
  rooResponsehist.SetName("response_noempty");
  h_flat_reco_pt1pt2->Scale(.5);
  h_flat_truth_pt1pt2->Scale(.5);
  
  h_flat_reco_measure_skim->Scale(.5);  
  h_flat_reco_skim->Scale(.5);
  h_flat_truth_skim->Scale(.5);

  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 3;

  for (int iter = 0; iter < niterations; iter++ )
    {
      
      RooUnfoldBayes   unfold (&rooResponsehist, h_flat_reco_measure_skim, iter + 1);    // OR

      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
      h_flat_unfold_pt1pt2[iter]->Reset();
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
      for (int ib = 0; ib < nbins*nbins; ib++)
	{
	  int bin = binnumbers_truth.at(ib);
	  if (bin)
	    {
	      h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin));
	      h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin));
	    }
	}
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

  h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");

  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

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


  cpt1pt2->Print("pt1pt2_half.pdf");

  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_truth_xj, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_truth_xj, kRed, 2, 24);

  dlutility::SetLineAtt(h_linear_truth_xj, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_linear_truth_xj, kRed, 2, 25);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 2, 24);

  dlutility::SetLineAtt(h_linear_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_linear_reco_xj, kBlue, 2, 25);

  
  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

  dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
  
  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  histo_opps::normalize_histo(h_truth_xj, nbins);
  histo_opps::normalize_histo(h_reco_xj, nbins);
  h_linear_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
  h_linear_reco_xj->Scale(1./h_linear_reco_xj->Integral(), "width");
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }

  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_truth->SetMaximum(5);
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
  TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 500);
  cresponseskim->SetLogz();
  h_flat_response_skim->Draw("colz");
  dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);


  TString responsepath = "response_matrix_half.root";
  if (JES_sys > 0)
    {
      responsepath = "response_matrix_half_posJES.root";
    }
  else if (JES_sys < 0)
    {
      responsepath = "response_matrix_half_negJES.root";
    }
  else if (JER_sys > 0)
    {
      responsepath = "response_matrix_half_JER.root";
    }
  else if (njet_sys > 0)
    {
      responsepath = "response_matrix_half_NJET.root";
    }
  
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  rooResponsehist.Write();

  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_flat_response_pt1pt2->Write();
  h_flat_truth_mapping->Write();
  h_flat_reco_mapping->Write();
  h_flat_truth_skim->Write();
  h_flat_reco_skim->Write();
  h_flat_response_skim->Write();

  // Full closure drawing
  h_xj_truth->Write();
  for (int iter =0 ; iter < niterations; iter++)
    {
      h_xj_unfold[iter]->Write();
    }
  h_xj_reco->Write();

  fr->Close();

  std::cout << "Reco empty: " << nempty_reco << std::endl;
  std::cout << "Truth empty: " << nempty_truth << std::endl;
  return 0;
}

