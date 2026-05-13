#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
void drawRadialDPHI()
{
  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  const int nbins_pt = read_nbins+1;
  int first_bin = 0;
  float ipt_bins[nbins_pt+1];
  double  dxj_bins[nbins+1];

  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] > 0.3 && first_bin == 0) first_bin = i;

    }
  ipt_bins[nbins_pt] = 100;

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

  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int measure_bins[10] = {0};
  int subleading_measure_bins[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
      subleading_measure_bins[ir] = rb.get_subleading_measure_region(ir);
      std::cout << ipt_bins[measure_bins[ir]] << " -- " <<  ipt_bins[subleading_measure_bins[ir]] << std::endl;
    }


  TFile *fin[6];
  TGraph *h_average_herwig[3];
  TGraph *h_average_truth[3];
  TGraphAsymmErrors *h_average_data[3];
  TGraph *h_rms_herwig[3];
  TGraph *h_rms_truth[3];
  TGraphAsymmErrors *h_rms_data[3];
  TGraphAsymmErrors *h_rms_data_stat[3];
  for (int i = 0; i < 3; i++)
    {
      h_average_herwig[i] = new TGraph(6);
      h_average_truth[i] = new TGraph(6);
      h_average_data[i] = new TGraphAsymmErrors(6);
      h_rms_herwig[i] = new TGraph(6);
      h_rms_truth[i] = new TGraph(6);
      h_rms_data[i] = new TGraphAsymmErrors(6);
      h_rms_data_stat[i] = new TGraphAsymmErrors(6);
    }
  for (int i = 0; i < 6;i++)
    {
      fin[i] = new TFile(Form("dphihists/dphi_summary_pp_r0%d.root", i+3),"r");

      TH1D *h_a_h = (TH1D*) fin[i]->Get("h_average_dphi_herwig");
      TH1D *h_a_t = (TH1D*) fin[i]->Get("h_average_dphi_truth");
      TH1D *h_a_d = (TH1D*) fin[i]->Get("h_average_dphi");
      TH1D *h_r_h = (TH1D*) fin[i]->Get("h_sigma_dphi_herwig");
      TH1D *h_r_t = (TH1D*) fin[i]->Get("h_sigma_dphi_truth");
      TH1D *h_r_d = (TH1D*) fin[i]->Get("h_sigma_dphi");
      for (int j = 0; j < 3; j++)
	{
	  h_average_herwig[j]->SetPoint(i, 0.3 + (0.1*i), h_a_h->GetBinContent(j+1));
	  h_average_truth[j]->SetPoint(i, 0.3 + (0.1*i), h_a_t->GetBinContent(j+1));
	  h_average_data[j]->SetPoint(i, 0.3 + (0.1*i), h_a_d->GetBinContent(j+1, 1));
	  h_average_data[j]->SetPointError(i, 0.05, 0.05, h_a_d->GetBinContent(j+1, 2), h_a_d->GetBinContent(j+1, 3));

	  h_rms_herwig[j]->SetPoint(i, 0.3 + (0.1*i), h_r_h->GetBinContent(j+1));
	  h_rms_truth[j]->SetPoint(i, 0.3 + (0.1*i), h_r_t->GetBinContent(j+1));
	  h_rms_data[j]->SetPoint(i, 0.3 + (0.1*i), h_r_d->GetBinContent(j+1, 1));
	  h_rms_data[j]->SetPointError(i, 0.05, 0.05, h_r_d->GetBinContent(j+1, 2), h_r_d->GetBinContent(j+1, 3));
	  h_rms_data_stat[j]->SetPoint(i, 0.3 + (0.1*i), h_r_d->GetBinContent(j+1, 1));
	  h_rms_data_stat[j]->SetPointError(i, 0.05, 0.05, h_r_d->GetBinError(j+1, 1), h_r_d->GetBinError(j+1, 1));
	  std::cout << i << " / " << j << " : " << h_r_d->GetBinContent(j+1, 1) << " +- " << h_r_d->GetBinError(j+1, 1) <<std::endl;
	}
    }
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  TCanvas *c = new TCanvas("c","c",500, 500);
  int color_range[3] = {kAzure-2, kRed, kGreen+3};
  int marker_range[3] = {22, 21, 20};
  for (int j = 0; j < 3; j++)
    {
      dlutility::SetLineAtt(h_rms_truth[j], color_range[j], 3, 1);
      dlutility::SetLineAtt(h_rms_herwig[j], color_range[j], 3, 4);

      dlutility::SetLineAtt(h_rms_data_stat[j], color_range[j], 2, 1);
      dlutility::SetMarkerAtt(h_rms_data_stat[j], color_range[j], 1, marker_range[j]);

      dlutility::SetLineAtt(h_rms_data[j], color_range[j], 1, 1);
      dlutility::SetMarkerAtt(h_rms_data[j], color_range[j], 1, marker_range[j]);
      h_rms_data[j]->SetFillColorAlpha(color_range[j], 0.3);

    }
  TH1D *hb = new TH1D("hb","hb", 6, 0.25, 0.85);
  hb->GetYaxis()->SetTitleOffset(1.7);
  dlutility::SetFont(hb, 42, 0.05, 0.05, 0.05, 0.05);


  hb->SetTitle(";Jet Radius (#it{R}); #sigma(#Delta#phi)");
  hb->SetMaximum(0.4);
  hb->SetMinimum(0.1);

  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.17);

  hb->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_rms_data[i]->Draw("same p E2");
      //      h_rms_data_stat[i]->Draw("same p");
    }
  for (int i = 0; i < 3; i++)
    {
      h_rms_truth[i]->Draw("l same");
    }
  for (int i = 0; i < 3; i++)
    {
      h_rms_herwig[i]->Draw("l same");
    }

  dlutility::DrawSPHENIXpp(0.2, 0.87);

  dlutility::drawText("anti-#it{k}_{t}", 0.2, 0.87 - 2*0.05);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.2, 0.87 - 3*0.05);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.2, 0.87 - 4*0.05);

  TLegend *legrad = new TLegend(0.52, 0.72, 0.71, 0.92);
  legrad->SetHeader("Data");
  legrad->SetLineWidth(0);
  legrad->SetTextSize(0.03);
  legrad->SetTextFont(42);
  //legrad->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  int iir = 0;
  legrad->AddEntry(h_rms_data[0],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  iir = 1;
  legrad->AddEntry(h_rms_data[1],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV",ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]] ));
  iir = 2;
  legrad->AddEntry(h_rms_data[2],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  legrad->Draw("same");

  TLegend *legrad2 = new TLegend(0.52, 0.58, 0.71, 0.72);
  legrad2->SetHeader("MC");
  legrad2->SetLineWidth(0);
  legrad2->SetTextSize(0.03);
  legrad2->SetTextFont(42);
  //legrad2->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  legrad2->AddEntry(h_rms_truth[0],"PYTHIA-8");
  legrad2->AddEntry(h_rms_herwig[0],"HERWIG 7.3");
  legrad2->Draw("same");

  c->Print("h_final_dphi_radial_scan.pdf");  
}
