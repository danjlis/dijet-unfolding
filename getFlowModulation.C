
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int unfoldData_noempty_AA(const std::string configfile = "binning_AA.config", const int niterations = 10, const int cone_size = 3, const int centrality_bin = 0)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());

  std::string data_file = rb.get_tntuple_location() + "/TNTUPLE_DIJET_r0" + std::to_string(cone_size) + "_v9_0_492_2024p020_v007_gl10-all.root";

  float mbd_vertex;
  float pt1_reco;
  float pt2_reco;
  float dphi_reco;
  float nrecojets;
  float trigger;
  float centrality;
  TFile *fin = new TFile(data_file.c_str(), "r");
  TNtuple *tn  = (TNtuple*) fin->Get("tn_dijet");;
  if (!tn)
    {
      std::cout << " no data "<< std::endl;
    }
  tn->SetBranchAddress("pt1_reco", &pt1_reco);
  tn->SetBranchAddress("pt2_reco", &pt2_reco);
  tn->SetBranchAddress("dphi_reco", &dphi_reco);
  tn->SetBranchAddress("njets", &nrecojets);
  tn->SetBranchAddress("centrality", &centrality);
  tn->SetBranchAddress("trigger", &trigger);
  tn->SetBranchAddress("mbd_vertex", &mbd_vertex);

  Int_t read_nbins = rb.get_nbins();
  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;

  Int_t njet_sys = rb.get_njet_sys();
  Int_t prior_sys = rb.get_prior_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  if (JER_sys != 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
    }
  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
    }
  if (prior_sys != 0)
    {
      std::cout << "Calculating prior extra = " << prior_sys  << std::endl;
    }

  float dphicut_z_low = rb.get_zyam_low();
  float dphicut_z_high = rb.get_zyam_high();

  float ZYAM_scale = (TMath::Pi() - dphicut)/(dphicut_z_high - dphicut_z_low);

  const int n_centrality_bins = rb.get_number_centrality_bins();  
  float icentrality_bins[n_centrality_bins+1];

  rb.get_centrality_bins(icentrality_bins);
  
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

  TH1D *h_dphi_reco = new TH1D("h_reco_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  TH1D *h_dphi_reco_ZYAM = new TH1D("h_reco_ZYAM_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  TH1D *h_dphi_reco_ETACUT = new TH1D("h_reco_ETACUT_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  TH1D *h_dphi_reco_SIGNAL = new TH1D("h_reco_SIGNAL_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  
  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";#it{p}_{T,1, data};#it{p}_{T,2, data}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_ZYAM = new TH2D("h_pt1pt2_ZYAM",";#it{p}_{T,1, ZYAM};#it{p}_{T,2, ZYAM}", nbins, ipt_bins, nbins, ipt_bins);

  int nbin_response = nbins*nbins;
  
  int entries = tn->GetEntries();
  for (int i = 0; i < entries; i++)
    {
      tn->GetEntry(i);

      if (centrality < icentrality_bins[centrality_bin] || centrality >= icentrality_bins[centrality_bin + 1]) continue;

      float maxi = std::max(pt1_reco, pt2_reco);
      float mini = std::min(pt1_reco, pt2_reco);
 
      bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco >= dphicut);
      bool fill_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut);
      bool ZYAM_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco >= dphicut_z_low && dphi_reco < dphicut_z_high);

      if (fill_good)
	{
	  h_dphi_reco->Fill(dphi_reco);
	}
      
      if (reco_good)
	{
	  h_pt1pt2->Fill(es1, es2);
	  h_pt1pt2->Fill(es2, es1);
	  h_dphi_reco_SIGNAL->Fill(dphi_reco);
	}
      else if (ZYAM_good)
	{

	  h_pt1pt2_ZYAM->Fill(es1, es2);
	  h_pt1pt2_ZYAM->Fill(es2, es1);
	  h_dphi_reco_ZYAM->Fill(dphi_reco);
	}
    }

  TH2D *h_pt1pt2_raw = (TH2D*) h_pt1pt2->Clone();
  h_pt1pt2_raw->SetName("h_pt1pt2_raw");
  

  dlutility::SetLineAtt(h_dphi_reco, kBlack, 2, 1);
  dlutility::SetLineAtt(h_dphi_reco_ZYAM, kBlack, 2, 4);
  h_dphi_reco_ZYAM->SetFillColorAlpha(kYellow, 0.4);
  dlutility::SetLineAtt(h_dphi_reco_SIGNAL, kBlack, 2, 4);
  h_dphi_reco_SIGNAL->SetFillColorAlpha(kSpring+2, 0.4);

  TH1D *h_dphi_reco_subtract = (TH1D*) h_dphi_reco->Clone();
  h_dphi_reco_subtract->SetName("h_dphi_reco_subtract");
  dlutility::SetLineAtt(h_dphi_reco_subtract, kRed, 2, 1);

  TH1D *h_dphi_reco_subtract_flow = (TH1D*) h_dphi_reco->Clone();
  h_dphi_reco_subtract_flow->SetName("h_dphi_reco_subtract_flow");
  dlutility::SetLineAtt(h_dphi_reco_subtract_flow, kBlue, 2, 1);

  float nbins_x = 0;
  float average_y = 0;
  for (int i = 0; i < h_dphi_reco_ZYAM->GetNbinsX(); i++)
    {
      if (h_dphi_reco_ZYAM->GetBinCenter(i+1) > dphicut_z_low && h_dphi_reco_ZYAM->GetBinCenter(i+1) < dphicut_z_high)
	{
	  average_y += h_dphi_reco_ZYAM->GetBinContent(i+1);
	  nbins_x += 1.0;
	}
    }

  average_y /= nbins_x;
  
  float integral_ZYAM = average_y * (dphicut_z_high - dphicut_z_low);

  float integral_mod = f_hand_modulation_p->Integral(dphicut_z_low, dphicut_z_high);

  float a = (integral_ZYAM)/(dphicut_z_high - dphicut_z_low + integral_mod);

  f_hand_modulation->SetParameter(0, a);

  float integral_mod_sig = f_hand_modulation_p->Integral(dphicut, TMath::Pi());
  float integral_mod_zyam = f_hand_modulation_p->Integral(dphicut_z_low, dphicut_z_high);
  dlutility::SetLineAtt(f_hand_modulation, kBlue, 2, 2);

  for (int i = 0; i < h_dphi_reco_subtract->GetNbinsX(); i++)
    {
      float v = h_dphi_reco_subtract->GetBinContent(i+1);
      h_dphi_reco_subtract->SetBinContent(i+1, v - average_y);
      h_dphi_reco_subtract_flow->SetBinContent(i+1, v - f_hand_modulation->Eval(h_dphi_reco_subtract_flow->GetBinCenter(i+1)));
    }

  h_pt1pt2_ZYAM->Scale(ZYAM_scale);//integral_mod_sig/integral_mod_zyam);
  
  h_pt1pt2->Add(h_pt1pt2_ZYAM, -1);
  for (int ix = 0; ix < h_pt1pt2->GetXaxis()->GetNbins(); ix++)
    {
      for (int iy = 0; iy < h_pt1pt2->GetXaxis()->GetNbins(); iy++)
	{
	  int xybin = h_pt1pt2->GetBin(ix+1, iy+1);
	  if (h_pt1pt2->GetBinContent(xybin) < 0)
	    {
	      h_pt1pt2->SetBinContent(xybin, 0);
	    }
	}
    }
  
  TH1D *h_flat_data_pt1pt2_raw = (TH1D*) h_flat_data_pt1pt2->Clone();
  h_flat_data_pt1pt2_raw->SetName("h_flat_data_pt1pt2_raw");

  h_flat_data_pt1pt2->Add(h_flat_data_pt1pt2_ZYAM, -1);
  for (int ix = 0; ix < h_flat_data_pt1pt2->GetNbinsX(); ix++)
    {
      if (h_flat_data_pt1pt2->GetBinContent(ix+1) < 0)
	{
	  h_flat_data_pt1pt2->SetBinContent(ix+1, 0);
	}
    }
  
  h_flat_data_pt1pt2->Scale(.5);

  TH1D *h_flat_data_skim = (TH1D*) h_flat_reco_skim->Clone();
  h_flat_data_skim->SetName("h_flat_data_skim");
  h_flat_data_skim->Reset();
  histo_opps::skim_down_histo(h_flat_data_skim, h_flat_data_pt1pt2, h_flat_reco_mapping); 
  
  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 1;
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

  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  for (int iter = 0; iter < niterations; iter++)
    {
      dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
    }
  
  dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
  dlutility::SetLineAtt(h_flat_data_pt1pt2, kBlue, 2, 1);

  h_flat_truth_pt1pt2->Draw("hist");
  h_flat_data_pt1pt2->Draw("hist same");
  h_flat_unfold_pt1pt2[niter]->Draw("same p");

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  h_pt1pt2_data->SetTitle(";Data #it{p}_{T, 1} [GeV]; Data #it{p}_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth #it{p}_{T, 1} [GeV]; Truth #it{p}_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold #it{p}_{T, 1} [GeV]; Unfold #it{p}_{T, 2} [GeV]; Counts * lumi scale ");
      
  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
  cpt1pt2->Divide(3, 1);
  cpt1pt2->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  h_pt1pt2_data->Draw("colz");
  dlutility::DrawSPHENIX(0.22, 0.85);


  cpt1pt2->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_truth->Draw("colz");
  cpt1pt2->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_unfold[niter]->Draw("colz");


  cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2data_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin,  cone_size));

  cpt1pt2->cd(1);

  h_pt1pt2_raw->GetYaxis()->SetTitleSize(0.06);
  h_pt1pt2_raw->GetXaxis()->SetTitleSize(0.06);

  h_pt1pt2_raw->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2_raw->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  TLine *id_l = new TLine(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45, 45);
  dlutility::SetLineAtt(id_l, kBlack, 2, 2);
  h_pt1pt2_raw->Draw("colz");
  id_l->Draw("same");
  dlutility::DrawSPHENIX(0.22, 0.85);
  
  cpt1pt2->cd(2);
  h_pt1pt2_ZYAM->GetYaxis()->SetTitleSize(0.06);
  h_pt1pt2_ZYAM->GetXaxis()->SetTitleSize(0.06);
  h_pt1pt2_ZYAM->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2_ZYAM->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2_ZYAM->Draw("colz");
  id_l->Draw("same");
  dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.85);
  cpt1pt2->cd(3);
  h_pt1pt2->GetYaxis()->SetTitleSize(0.06);
  h_pt1pt2->GetXaxis()->SetTitleSize(0.06);

  h_pt1pt2->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetYaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2->Draw("colz");
  id_l->Draw("same");
  cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2dataZYAM_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin,  cone_size));

  TCanvas *ctext = new TCanvas("ctext","ctext", 600, 500);
  ctext->SetRightMargin(0.2);
  
  TH2D *h_pt1pt2_text = (TH2D*) h_pt1pt2_ZYAM->Clone();
  h_pt1pt2_text->Divide(h_pt1pt2_raw);

  h_pt1pt2_text->GetYaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);
  h_pt1pt2_text->GetXaxis()->SetRangeUser(h_pt1pt2_raw->GetXaxis()->GetBinLowEdge(1), 45);

  gStyle->SetPaintTextFormat("1.2f");
  h_pt1pt2_text->Draw("colz");
  h_pt1pt2_text->Draw("text,same");
  id_l->Draw("same");
  dlutility::DrawSPHENIX(0.22, 0.85);
  dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.75);
  dlutility::drawText("Background Fraction", 0.22, 0.7);
  ctext->Print(Form("%s/unfolding_plots/pt1pt2dataTEXT_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin,  cone_size));

  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }

  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_reco_xj, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }


  TCanvas *cdphi = new TCanvas("cphi","cdphi", 500, 500);

  
  
  h_dphi_reco->SetMinimum(-20);
  h_dphi_reco->SetMaximum(h_dphi_reco->GetBinContent(h_dphi_reco->GetMaximumBin())*1.7);
  h_dphi_reco->Draw("hist");
  h_dphi_reco_ZYAM->Draw("hist same");
  h_dphi_reco_SIGNAL->Draw("hist same");
  h_dphi_reco_subtract->Draw("hist same");
  //h_dphi_reco_subtract_flow->Draw("hist same");
  //f_hand_modulation->Draw("same");
  TLine *lii = new TLine(0, 0, TMath::Pi(), 0);
  lii->Draw("same");
  TLine *lis = new TLine(0, average_y, TMath::Pi(), average_y);
  dlutility::SetLineAtt(lis, kRed, 2, 2);
  lis->Draw("same");

  dlutility::DrawSPHENIX(0.22, 0.85);
  dlutility::drawText(Form("#it{p}_{T,1} #geq %2.1f GeV", reco_leading_cut), 0.22, 0.75);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", reco_subleading_cut), 0.22, 0.7);
  dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.65);
  TLegend *lc = new TLegend(0.5, 0.63, 0.75, 0.78);
  lc->SetTextSize(0.03);
  lc->SetLineWidth(0);
  lc->AddEntry(h_dphi_reco, "No ZYAM sub.");
  lc->AddEntry(h_dphi_reco_subtract, "w/ ZYAM sub.");
  //lc->AddEntry(h_dphi_reco_subtract_flow, "w/ ZYAM sub. + flow");
  lc->AddEntry(h_dphi_reco_SIGNAL, "Signal region");
  lc->AddEntry(h_dphi_reco_ZYAM, "ZYAM region");
  lc->Draw("same");
  cdphi->Print(Form("%s/unfolding_plots/dphi_ZYAM_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin, cone_size));

  for (int j = 0; j < 4; j++)
    {
      dlutility::SetLineAtt(h_dphi_reco_sub[j], kBlack, 2, 1);
      dlutility::SetLineAtt(h_dphi_reco_ZYAM_sub[j], kBlack, 2, 4);
      h_dphi_reco_ZYAM_sub[j]->SetFillColorAlpha(kYellow, 0.4);
      dlutility::SetLineAtt(h_dphi_reco_SIGNAL_sub[j], kBlack, 2, 4);
      h_dphi_reco_SIGNAL_sub[j]->SetFillColorAlpha(kSpring+2, 0.4);

      TH1D *h_dphi_reco_subtract_sub = (TH1D*) h_dphi_reco_sub[j]->Clone();
      h_dphi_reco_subtract_sub->SetName("h_dphi_reco_subtract_sub");
      dlutility::SetLineAtt(h_dphi_reco_subtract_sub, kRed, 2, 1);

      float nbins_x1 = 0;
      float average_y1 = 0;
      for (int i = 0; i < h_dphi_reco_ZYAM_sub[j]->GetNbinsX(); i++)
	{
	  if (h_dphi_reco_ZYAM_sub[j]->GetBinCenter(i+1) > dphicut_z_low && h_dphi_reco_ZYAM_sub[j]->GetBinCenter(i+1) < dphicut_z_high)
	    {
	      average_y1 += h_dphi_reco_ZYAM_sub[j]->GetBinContent(i+1);
	      nbins_x1 += 1.0;
	    }
	}

      average_y1 /= nbins_x1;
  
      for (int i = 0; i < h_dphi_reco_subtract_sub->GetNbinsX(); i++)
	{
	  float v = h_dphi_reco_subtract_sub->GetBinContent(i+1);
	  h_dphi_reco_subtract_sub->SetBinContent(i+1, v - average_y1);
	}

      h_dphi_reco_sub[j]->SetMinimum(-5);
      h_dphi_reco_sub[j]->SetMaximum(h_dphi_reco_sub[j]->GetBinContent(h_dphi_reco_sub[j]->GetMaximumBin())*1.7);
      h_dphi_reco_sub[j]->Draw("hist");
      h_dphi_reco_ZYAM_sub[j]->Draw("hist same");
      h_dphi_reco_SIGNAL_sub[j]->Draw("hist same");
      h_dphi_reco_subtract_sub->Draw("hist same");
      //h_dphi_reco_subtract_flow->Draw("hist same");
      //f_hand_modulation->Draw("same");
      //TLine *lii = new TLine(0, 0, TMath::Pi(), 0);
      lii->Draw("same");
      TLine *lis2 = new TLine(0, average_y1, TMath::Pi(), average_y1);
      dlutility::SetLineAtt(lis2, kRed, 2, 2);
      lis2->Draw("same");

      dlutility::DrawSPHENIX(0.22, 0.85);
      dlutility::drawText(Form("#it{p}_{T,1} #geq %2.1f GeV", reco_leading_cut), 0.22, 0.75);
      if (j == 0)
	dlutility::drawText(Form("%2.1f #leq #it{p}_{T,2}^{lead} < %2.1f GeV", ipt_bins[1], ipt_bins[2]), 0.22, 0.7);
      if (j == 1)
	dlutility::drawText(Form("%2.1f #leq #it{p}_{T,2}^{lead} < %2.1f GeV", ipt_bins[2], ipt_bins[3]), 0.22, 0.7);
      if (j == 2)
	dlutility::drawText(Form("%2.1f #leq #it{p}_{T,2}^{lead} < %2.1f GeV", ipt_bins[3], ipt_bins[5]), 0.22, 0.7);
      if (j == 3)
	dlutility::drawText(Form("%2.1f #leq #it{p}_{T,2}^{lead} < %2.1f GeV", ipt_bins[5], ipt_bins[7]), 0.22, 0.7);

      dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.22, 0.65);
      lc = new TLegend(0.65, 0.7, 0.75, 0.85);
      lc->SetTextSize(0.03);
      lc->SetLineWidth(0);
      lc->AddEntry(h_dphi_reco_sub[j], "No ZYAM sub.");
      lc->AddEntry(h_dphi_reco_subtract_sub, "w/ ZYAM sub.");
      //lc->AddEntry(h_dphi_reco_subtract_flow, "w/ ZYAM sub. + flow");
      lc->AddEntry(h_dphi_reco_SIGNAL_sub[j], "Signal region");
      lc->AddEntry(h_dphi_reco_ZYAM_sub[j], "ZYAM region");
      lc->Draw("same");
      cdphi->Print(Form("%s/unfolding_plots/dphi_ZYAM_AA_sub_%d_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), j, centrality_bin, cone_size));

    }

  
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 1, 24);


  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  //h_xj_truth->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_truth->SetMaximum(6);
  h_xj_truth->SetMinimum(0);
  h_xj_truth->Draw("p");
  h_xj_unfold[niter]->Draw("same p");
  h_xj_data->Draw("hist same");
  h_reco_xj->Draw("hist same");
  h_xj_truth->Draw("same hist");
  h_xj_data->Draw("p same");
  h_xj_unfold[niter]->Draw("same hist");
  h_xj_unfold[niter]->Draw("same p");
  dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("#it{p}_{T,2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/3", 0.22, 0.59);
  dlutility::drawText("\\mathscr{L} = 25.7 pb^{-1}", 0.22, 0.54);

  TLegend *leg = new TLegend(0.22, 0.35, 0.4, 0.50);
  leg->SetLineWidth(0);
  leg->AddEntry(h_xj_data, "Data");
  leg->AddEntry(h_xj_truth, "Pythia8");
  leg->AddEntry(h_xj_unfold[niter], "Unfolded");
  leg->Draw("same");
  dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
    
  cxj->cd(2);

  TH1D *h_data_compare = (TH1D*) h_xj_unfold[niter]->Clone();
  h_data_compare->Divide(h_xj_truth);
  //h_data_compare->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_data_compare->SetTitle(";x_{J}; Unfold / Truth");
  dlutility::SetFont(h_data_compare, 42, 0.06);
  dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);

 
  h_data_compare->Draw("p");
  TLine *line = new TLine(0.3, 1, 1, 1);
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



  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlack, 1, 24);
  
  TCanvas *cproj = new TCanvas("cproj","cproj", 500, 700);
  dlutility::ratioPanelCanvas(cproj);

  cproj->cd(1);
  h_xj_data->SetMaximum(3);
  h_xj_data->SetTitle(";x_{J}; #frac{1}{N_{pairs}}#frac{dN_{pair}}{dx_{J}}");
  dlutility::SetFont(h_xj_data, 42, 0.04);
  
  h_xj_data->Draw();
  h_reco_xj->Draw("same");
  dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("#it{p}_{T,2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
  TLegend *legp = new TLegend(0.22, 0.45, 0.4, 0.55);
  legp->SetLineWidth(0);
  legp->SetTextFont(42);
  legp->SetTextSize(0.04);
  legp->AddEntry(h_reco_xj, "Data Filled");
  legp->AddEntry(h_xj_data, "Data Projected");
  legp->Draw("same");

  cproj->cd(2);
  TH1D *h_fillproj_compare = (TH1D*) h_xj_data->Clone();

  dlutility::SetFont(h_fillproj_compare, 42, 0.06);
  h_fillproj_compare->Divide(h_reco_xj);

  h_fillproj_compare->SetMaximum(1.2);
  
  h_fillproj_compare->SetMinimum(0.8);
  h_fillproj_compare->SetTitle(";x_{J};Projected / Filled");
  dlutility::SetLineAtt(h_fillproj_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_fillproj_compare, kBlack, 1, 8);

  h_fillproj_compare->Draw();
  TLine *line3 = new TLine(0.1, 1, 1, 1);
  line3->SetLineStyle(4);
  line3->SetLineColor(kRed + 3);
  line3->SetLineWidth(2);
  line3->Draw("same");
  cproj->Print(Form("%s/unfolding_plots/proj_compar_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin, cone_size));
  
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_AA_cent_" + std::to_string(centrality_bin) + "_r0" + std::to_string(cone_size);
  if (JES_sys > 0)
    {
      unfoldpath += "_posJES";
    }
  if (JES_sys < 0)
    {
      unfoldpath += "_negJES";
    }

  else if (JER_sys > 0)
    {
      unfoldpath += "_posJER";
    }
  else if (JER_sys < 0)
    {
      unfoldpath += "_negJER";
    }
  else if (njet_sys > 0)
    {
      unfoldpath += "_NJET";
    }
  else if (prior_sys)
    {
      unfoldpath += "_PRIOR";
    }

  if (primer)
    {
      unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_AA_cent_" + std::to_string(centrality_bin) + "_r0" + std::to_string(cone_size) + "_PRIMER" + std::to_string(primer);
    }

  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");
  
  h_flat_data_skim->Write();
  h_flat_data_pt1pt2->Write();
  h_flat_data_pt1pt2_raw->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_mbd_vertex->Write();
  h_centrality->Write();
  
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_unfold_pt1pt2[iter]->Write();
    }
  fout->Close();
  

  return 0;
  
}
