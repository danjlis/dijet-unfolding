#include "../macros/dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

void makeIterationPlot(const int cone_size = 4)
{
  const int niterations = 10;

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

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

  int measured_bins[5] = { 7, 9, 10, 11, 14};
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

  float truth_leading_goal = 13;
  float truth_subleading_goal = 0;

  float dphicut = 3*TMath::Pi()/4.;

  // Make binning
  for (int i = 0; i < nbins+1; i++)
    {
      float ipt = bin_min*TMath::Power(alpha, (float)i);
      float ixj = bin_xj1*TMath::Power(alpha, (float)i);
      ipt_bins[i] = ipt;
      ixj_bins[i] = ixj;
      if (measure_leading_bin == 0 && ipt >= measure_leading_goal) measure_leading_bin = i;
      if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
      if (truth_subleading_bin == 0 && ipt >= truth_subleading_goal) truth_subleading_bin = i;
      std::cout << i << " : " <<  ipt << " -- " << ixj <<  std::endl;
    }
  std::cout << __LINE__ << std::endl;
  truth_leading_cut=ipt_bins[truth_leading_bin];
  truth_subleading_cut=ipt_bins[truth_subleading_bin];

  reco_leading_cut=ipt_bins[truth_leading_bin + 2];
  reco_subleading_cut=ipt_bins[truth_subleading_bin + 2];

  reco_meas_leading_cut=ipt_bins[truth_leading_bin + 3];
  reco_meas_subleading_cut=ipt_bins[truth_subleading_bin + 3];
  measure_subleading_bin = truth_subleading_bin + 3;

  
  std::cout << __LINE__ << std::endl;

  TFile *f_uncertainties = new TFile(Form("uncertainties/uncertainties_r%02d.root", cone_size),"r");

  if (!f_uncertainties)
    {
      std::cout << "no uncertainities" << std::endl;
      return;
    }
  TProfile *hp_xj[niterations];
    
  TFile *fin = new TFile(Form("unfolded_hists/unfolded_hists_r%02d.root", cone_size),"r");
  if (!fin)
    {
      std::cout << "no file" << std::endl;
      return;
    }

  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      hp_xj[iter] = (TProfile*) f_uncertainties->Get(Form("hp_xj_%d", iter));
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }
  std::cout << __LINE__ << std::endl;
  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }
  std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  TH1D *h_xj_profile_unfold[niterations];  
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
      h_xj_profile_unfold[iter] = new TH1D(Form("h_xj_profile_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  std::cout << __LINE__ << std::endl;
  TH1D *h_xjunc_data = new TH1D("h_xjunc_data", ";x_{J};", nbins, ixj_bins);
  std::cout << __LINE__ << std::endl;
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      std::cout << __LINE__ << std::endl;
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  std::cout << __LINE__ << std::endl;
  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      histo_opps::tprofile_to_histo(hp_xj[iter], h_xj_profile_unfold[iter], nbins);
    }
  histo_opps::normalize_histo(h_xj_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
      histo_opps::normalize_histo(h_xj_profile_unfold[iter], nbins);
    }

  TH1D *h_final_xj_data = (TH1D*) h_xj_data->Clone();
  h_final_xj_data->SetName("h_final_xj_data");
  h_final_xj_data->Reset();
  histo_opps::finalize_xj(h_xj_data,h_final_xj_data, nbins, 0.4);
  TH1D * h_final_xj_unfold[niterations];
  TH1D * h_final_xj_profile_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_final_xj_unfold[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
      h_final_xj_unfold[iter]->SetName(Form("h_final_xj_unfold_%d", iter));
      h_final_xj_unfold[iter]->Reset();
      h_final_xj_profile_unfold[iter] = (TH1D*) h_xj_profile_unfold[iter]->Clone();
      h_final_xj_profile_unfold[iter]->SetName(Form("h_final_xj_profile_unfold_%d", iter));
      h_final_xj_profile_unfold[iter]->Reset();

      histo_opps::finalize_xj(h_xj_unfold[iter], h_final_xj_unfold[iter], nbins, 0.4);
      histo_opps::finalize_xj(h_xj_profile_unfold[iter], h_final_xj_profile_unfold[iter], nbins, 0.4);
    }  
  
  TH1D *h_statistical_uncertainties = new TH1D("h_statistical_uncertainties","; N_{iter}; #sigma_{conv}",niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_unfold_uncertainties = new TH1D("h_unfold_uncertainties","; N_{iter}; #sigma_{conv}",niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_conv_uncertainties = new TH1D("h_conv_uncertainties","; N_{iter}; #sigma_{conv}",niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_binbybin_uncertainties = new TH1D("h_binbybin_uncertainties","; N_{iter}; #sigma_{conv}",niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_total_uncertainties = new TH1D("h_total_uncertainties","; N_{iter}; #sigma_{conv}",niterations + 1, -0.5, niterations + 0.5);
  std::cout << __LINE__ << std::endl;
  for (int iter = 1; iter < niterations; ++iter)
    {

      Double_t stat_unc = 0;
      Double_t unfold_unc = 0;
      Double_t binbybin_unc = 0;
      for (int ibin = 0; ibin < nbins; ibin++)
	{
	  if (iter == 0)
	    {
	      binbybin_unc += TMath::Power(h_final_xj_data->GetBinContent(ibin+1) - h_final_xj_unfold[iter]->GetBinContent(ibin+1),2);
	    }
	  else
	    {
	      binbybin_unc += TMath::Power(h_final_xj_unfold[iter - 1]->GetBinContent(ibin+1) - h_final_xj_unfold[iter]->GetBinContent(ibin+1),2);
	    }

	  stat_unc+= TMath::Power(h_final_xj_profile_unfold[iter]->GetBinError(ibin + 1), 2);
	  unfold_unc+= TMath::Power(h_final_xj_unfold[iter]->GetBinError(ibin + 1), 2);

	}

      std::cout << stat_unc << " + " << unfold_unc <<" + " << binbybin_unc << " = " <<  sqrt(stat_unc + unfold_unc + binbybin_unc) << std::endl;
      Double_t conv_unc = stat_unc + unfold_unc;

      stat_unc = sqrt(stat_unc);
      unfold_unc = sqrt(unfold_unc);

      Double_t total_unc = sqrt(binbybin_unc + conv_unc); 

      binbybin_unc = sqrt(binbybin_unc);
      conv_unc = sqrt(conv_unc);

      h_statistical_uncertainties->Fill(iter + 1, stat_unc);
      h_unfold_uncertainties->Fill(iter + 1, unfold_unc);
      h_conv_uncertainties->Fill(iter + 1, conv_unc);
      h_binbybin_uncertainties->Fill(iter + 1, binbybin_unc);
      h_total_uncertainties->Fill(iter + 1, total_unc);
    }


  TCanvas *c_unc = new TCanvas("c_unc","c_unc", 500, 500);
  dlutility::SetMarkerAtt(h_total_uncertainties, kBlack, 1, 8);
  dlutility::SetLineAtt(h_total_uncertainties, kBlack, 1, 1);

  dlutility::SetMarkerAtt(h_binbybin_uncertainties, kRed, 1, 8);
  dlutility::SetLineAtt(h_binbybin_uncertainties, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_conv_uncertainties, kBlue, 1, 8);
  dlutility::SetLineAtt(h_conv_uncertainties, kBlue, 1, 1);

  h_total_uncertainties->SetMaximum(0.2);
  h_total_uncertainties->Draw("p hist");
  h_conv_uncertainties->Draw("same p hist");
  h_binbybin_uncertainties->Draw("same p hist");

  dlutility::DrawSPHENIXpp(0.22, 0.87);
  dlutility::drawText("#sigma_{conv} = #sqrt{#sigma^{2}_{stat} + #sigma^{2}_{unfold} + #sigma^{2}_{bin-by-bin}}", 0.22, 0.77);
  TLegend *leg = new TLegend(0.218, 0.566, 0.397, 0.726);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_conv_uncertainties, "#sqrt{#sigma^{2}_{stat} + #sigma^{2}_{unfold}}","p");
  leg->AddEntry(h_binbybin_uncertainties, "#sigma_{bin-by-bin}","p");
  leg->AddEntry(h_total_uncertainties, "#sigma_{conv}","p");
  leg->Draw("same");

  c_unc->SaveAs(Form("unfolding_plots/iteration_tune_r%02d.pdf", cone_size));
  c_unc->SaveAs(Form("unfolding_plots/iteration_tune_r%02d.png", cone_size));
  return;
}
