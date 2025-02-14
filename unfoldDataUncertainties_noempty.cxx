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
int unfoldDataUncertainties_noempty(const int niterations = 20)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  std::string data_file = "TNTUPLE_DIJET_v6_2_ana462_2024p010_v001_gl10-00047289-00048291.root";//TNTUPLE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root";

  float pt1_reco;
  float pt2_reco;
  float dphi_reco;

  TFile *fresponse = new TFile("response_matrix.root","r");
  
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

  TH2D *h_flat_response_pt1pt2 = (TH2D*) fresponse->Get("h_flat_response_pt1pt2"); 

  if (!h_flat_response_pt1pt2)
    {
      std::cout << "no response" << std::endl;
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

  TH2D *h_flat_response_skim = (TH2D*) fresponse->Get("h_flat_response_skim"); 

  if (!h_flat_response_skim)
    {
      std::cout << "no response" << std::endl;
      return 1;
    }

  TFile *fin = new TFile(data_file.c_str(), "r");
  TNtuple *tn  = (TNtuple*) fin->Get("tn_dijet");;
  if (!tn)
    {
      std::cout << " no data "<< std::endl;
    }
  tn->SetBranchAddress("pt1_reco", &pt1_reco);
  tn->SetBranchAddress("pt2_reco", &pt2_reco);
  tn->SetBranchAddress("dphi_reco", &dphi_reco);

  read_binning rb("binning.config");

  Int_t read_nbins = rb.get_nbins();
  
  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  const int nbins2 = nbins*nbins;

  Double_t ipt_bins[nbins+1];
  Double_t ixj_bins[nbins+1];

  float fipt_bins[nbins+1];
  float fixj_bins[nbins+1];

  rb.get_pt_bins(fipt_bins);
  rb.get_xj_bins(fixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      ipt_bins[i] = fipt_bins[i];
      ixj_bins[i] = fixj_bins[i];
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


  TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_flat_data_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  int entries = tn->GetEntries();
  for (int i = 0; i < entries; i++)
    {
      tn->GetEntry(i);

      float maxi = std::max(pt1_reco, pt2_reco);
      float mini = std::min(pt1_reco, pt2_reco);

      float pt1_reco_bin = nbins;
      float pt2_reco_bin = nbins;

      float es1 = pt1_reco;
      float es2 = pt2_reco;
      if (es1 > 70 ) continue;
      for (int ib = 0; ib < nbins; ib++)
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
	  
	 
      bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco > dphicut);

      if (!reco_good) continue;
      
      h_pt1pt2->Fill(es1, es2);
      h_pt1pt2->Fill(es2, es1);
      h_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
      h_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
      h_reco_xj->Fill(mini/maxi);
    }

  std::cout << __LINE__ <<" :: " << "Done with going through data" << std::endl;
  h_flat_data_pt1pt2->Scale(.5);

  TH1D *h_flat_data_skim = (TH1D*) h_flat_reco_skim->Clone();
  h_flat_data_skim->SetName("h_flat_data_skim");
  h_flat_data_skim->Reset();
  histo_opps::skim_down_histo(h_flat_data_skim, h_flat_data_pt1pt2, h_flat_reco_mapping);
  TH1D *h_flat_unfold_pt1pt2[niterations];
  
  int niter = 1;

  TProfile2D *hp_pt1pt2_stats[niterations];
  TProfile *hp_xj[niterations];
  TProfile *hp_xj_range[mbins][niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      hp_pt1pt2_stats[iter] = new TProfile2D(Form("hp_pt1pt2_stats_%d", iter),";p_{T,1}; p_{T,2}; average +/- rms", nbins, ipt_bins, nbins, ipt_bins, "s");
      
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
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
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
      TH2D *h_flat_response_clone = (TH2D*) h_flat_response_skim->Clone();
      for (int ibin = 0; ibin < h_flat_response_clone->GetXaxis()->GetNbins(); ++ibin)
	{
	  for (int jbin = 0; jbin < h_flat_response_clone->GetYaxis()->GetNbins(); ++jbin)
	    {
	      int ijbin = h_flat_response_skim->GetBin(ibin+1, jbin + 1);
	      float smear = ferror_response->GetRandom() * h_flat_response_skim->GetBinError(ijbin);
	      float value = h_flat_response_clone->GetBinContent(ijbin);
	      float newvalue = value + smear;
	      if (newvalue < 0) newvalue = 0;
	      h_flat_response_clone->SetBinContent(ijbin, newvalue);
	    }
	}
      
      RooUnfoldResponse rooResponse(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_clone);
      for (int iter = 0; iter < niterations; iter++ )
	{
	  RooUnfoldBayes   unfold (&rooResponse, h_flat_data_skim, iter + 1);    // OR
	  TH1D *h_flat_unfold_skim = (TH1D*) unfold.Hunfold();
	  h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
	  h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_flat_unfold_pt1pt2[iter]->Reset();
	  
	  histo_opps::fill_up_histo(h_flat_unfold_skim, h_flat_unfold_pt1pt2[iter], h_flat_truth_mapping);
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
      

      for (int ib = 0; ib < nbins*nbins; ib++)
	{
	  int xbin = ib/nbins;
	  int ybin = ib%nbins;
	  
	  int b = h_pt1pt2_unfold[0]->GetBin(xbin+1, ybin+1);
	  
	  for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinContent(b, h_flat_unfold_pt1pt2[iter]->GetBinContent(ib+1));
	  for (int iter = 0; iter < niterations; iter++ ) h_pt1pt2_unfold[iter]->SetBinError(b, h_flat_unfold_pt1pt2[iter]->GetBinError(ib+1));
	}

      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pt1pt2_unfold_sym[iter] = (TH2D*) h_pt1pt2_unfold[iter]->Clone();
	  histo_opps::project_xj(h_pt1pt2_unfold_sym[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	  for (int irange = 0; irange < mbins; irange++)
	    {
	      histo_opps::project_xj(h_pt1pt2_unfold_sym[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
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
  

  
  TFile *fout = new TFile("uncertainties.root","recreate");
  TEnv *penv = new TEnv("binning.config");
  penv->Write();
  for (int iter = 0; iter < niterations; iter++)
    {
      hp_xj[iter]->Write();            
      for (int irange = 0; irange < mbins; irange++)
	{
	  hp_xj_range[irange][iter]->Write();            
	}
      hp_pt1pt2_stats[iter]->Write();      
    }

  fout->Close();

  return 0;
  
}
