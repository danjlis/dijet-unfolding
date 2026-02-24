#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

const bool NUCLEAR = true;
const int color_unfold_fill = kAzure - 4;
const int color_unfold = kAzure - 5;
const float marker_unfold = 20;
const float msize_unfold = 1.1;
const float lsize_unfold = 1.1;

const int color_pythia = kRed + 1;
const float msize_pythia = 1.1;
const float marker_pythia = 21;
const float lsize_pythia = 3;

const int color_herwig = kSpring - 1;
const float msize_herwig = 1.3;
const float marker_herwig = 22;
const float lsize_herwig = 3;

const int color_reco = kRed;
const float marker_reco = 21;
const float msize_reco = 0.9;
const float lsize_reco = 1.1;
const int color_data = kAzure - 6;
const float marker_data = 24;
const float msize_data = 0.9;
const float lsize_data = 1.1;
void drawFinalUnfold(const int cone_size = 4)
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  int first_bin = 0;
  float ipt_bins[nbins+1];
  double  dxj_bins[nbins+1];

  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] > 0.3 && first_bin == 0) first_bin = i-1;

    }
  ixj_bins[nbins] = 1;

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

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  

  // Get truth histograms without any reweighting

  TFile *fintrh = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_PRIMER1_HERWIG.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fintrh)
    {
      std::cout << "no herwig hists" << std::endl;
      return;
    }
  TH1D *h_flat_herwig_pt1pt2 = (TH1D*) fintrh->Get("h_truth_flat_pt1pt2");
  h_flat_herwig_pt1pt2->SetName("h_flat_herwig_pt1pt2");
  
  TFile *fintr = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_PRIMER1_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fintr)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fintr->Get("h_truth_flat_pt1pt2");


  
  const int niterations = 10;
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finu)
    {
      std::cout << " no unc " << std::endl;
      return;
    }

  // Get statistical uncertainties of the unfold
  TProfile *h_xj_rms[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  h_xj_rms[irange][iter] = (TProfile*) finu->Get(Form("hp_xj_range_%d_%d", irange, iter));
	}
    }

  // Get systematics
  TFile *fins = new TFile(Form("%s/uncertainties/systematics_pp_r%02d.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fins)
    {
      std::cout << " no sys " << std::endl;
      return;
    }
  TH1D *h_total_sys_range[mbins][niterations];
  TH1D *h_total_sys_neg_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  h_total_sys_range[irange][iter] = (TH1D*) fins->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	  h_total_sys_neg_range[irange][iter] = (TH1D*) fins->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));
	}
    }


  // Get the pt1pt2 histograms
  TFile *fin = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fin)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  
  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");

  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_herwig = new TH2D("h_pt1pt2_herwig", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }


  TH1D *h_xj_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_herwig_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_data_range[irange] = new TH1D(Form("h_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_herwig_range[irange] = new TH1D(Form("h_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  TH1D *h_xjunc_data_range[mbins];
  TH1D *h_xjunc_reco_range[mbins];
  TH1D *h_xjunc_truth_range[mbins];
  TH1D *h_xjunc_herwig_range[mbins];
  TH1D *h_xjunc_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xjunc_data_range[irange] = new TH1D(Form("h_xjunc_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_reco_range[irange] = new TH1D(Form("h_xjunc_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_truth_range[irange] = new TH1D(Form("h_xjunc_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_herwig_range[irange] = new TH1D(Form("h_xjunc_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
       for (int iter = 0; iter < niterations; iter++)
	{
	  h_xjunc_unfold_range[irange][iter] = new TH1D(Form("h_xjunc_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_herwig_pt1pt2, h_pt1pt2_herwig, nbins);  
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }


  TH1D *h_final_xj_data_range[mbins];
  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  TH1D *h_final_xj_herwig_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];

  TH1D *h_final_xj_rms[mbins][niterations];
  TH1D *h_final_xj_systematics[mbins][niterations];
  TH1D *h_final_xj_statistics[mbins][niterations];

  TGraphAsymmErrors *g_final_xj_systematics[mbins][niterations];
  TGraphAsymmErrors *g_final_xj_statistics[mbins][niterations];
  TGraph *g_final_xj_truth[mbins];
  TGraph *g_final_xj_herwig[mbins];

  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_data_range[irange] = new TH1D(Form("h_final_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_herwig_range[irange] = new TH1D(Form("h_final_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
   
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_final_xj_rms[irange][iter] = new TH1D(Form("h_final_xj_rms_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_final_xj_systematics[irange][iter] = new TH1D(Form("h_final_xj_systematics_%d_%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_final_xj_statistics[irange][iter] = new TH1D(Form("h_final_xj_statistics_%d_%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_data, h_xj_data_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_herwig_range[irange], h_final_xj_herwig_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_data_range[irange], h_final_xj_data_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);

      histo_opps::normalize_histo(h_final_xj_truth_range[irange], nbins);
      histo_opps::normalize_histo(h_final_xj_herwig_range[irange], nbins);
      histo_opps::normalize_histo(h_final_xj_data_range[irange], nbins);
      histo_opps::normalize_histo(h_final_xj_reco_range[irange], nbins);

      g_final_xj_truth[irange] = new TGraph(h_final_xj_truth_range[irange]);
      g_final_xj_herwig[irange] = new TGraph(h_final_xj_herwig_range[irange]);

      histo_opps::trim_tgraph(g_final_xj_truth[irange], nbins, first_xj);
      histo_opps::trim_tgraph(g_final_xj_herwig[irange], nbins, first_xj);

      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	  histo_opps::finalize_xj(h_xj_rms[irange][iter], h_final_xj_rms[irange][iter], nbins, first_xj);

	  histo_opps::normalize_histo(h_final_xj_unfold_range[irange][iter], nbins);
	  histo_opps::normalize_histo(h_final_xj_rms[irange][iter], nbins);

	  histo_opps::set_xj_errors(h_final_xj_unfold_range[irange][iter], h_final_xj_rms[irange][iter], nbins);

	  
	  h_final_xj_systematics[irange][iter] = (TH1D*) h_final_xj_unfold_range[irange][iter]->Clone();
	  h_final_xj_systematics[irange][iter]->SetName(Form("h_final_xj_systematics_%d_%d", irange, iter));

	  g_final_xj_systematics[irange][iter] = new TGraphAsymmErrors(h_final_xj_systematics[irange][iter]);
	  
	  histo_opps::trim_tgraph(g_final_xj_systematics[irange][iter], nbins, first_xj);
	  
	  histo_opps::get_xj_systematics(g_final_xj_systematics[irange][iter], h_total_sys_neg_range[irange][iter], h_total_sys_range[irange][iter], nbins);

	  h_final_xj_statistics[irange][iter] = (TH1D*) h_final_xj_unfold_range[irange][iter]->Clone();
	  h_final_xj_statistics[irange][iter]->SetName(Form("h_final_xj_statistics_%d_%d", irange, iter));
	  g_final_xj_statistics[irange][iter] = histo_opps::get_xj_statistics(h_final_xj_statistics[irange][iter], nbins);

	  histo_opps::trim_tgraph(g_final_xj_statistics[irange][iter], nbins, first_xj);
	}

    }

  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 500);
  dlutility::ratioPanelCanvas(cxj);
  int niter = 1;

  for (int irange = 0; irange < mbins; irange++)
    {
      cxj->cd(1);
      dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(g_final_xj_systematics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_systematics[irange][niter], color_unfold, msize_unfold, marker_unfold);
      g_final_xj_systematics[irange][niter]->SetFillColorAlpha(color_unfold_fill, 0.3); 
      dlutility::SetLineAtt(g_final_xj_statistics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_statistics[irange][niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
      dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

      dlutility::SetLineAtt(h_final_xj_data_range[irange], color_data, lsize_data, 1);
      dlutility::SetMarkerAtt(h_final_xj_data_range[irange], color_data, msize_data, marker_data);
      dlutility::SetLineAtt(h_final_xj_reco_range[irange], color_reco, lsize_reco, 1);
      dlutility::SetMarkerAtt(h_final_xj_reco_range[irange], color_reco, msize_reco, marker_reco);

      dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.05);
      h_final_xj_truth_range[irange]->SetMaximum(5);
      h_final_xj_truth_range[irange]->SetMinimum(0);
      h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      TH1D *ht = (TH1D*) h_final_xj_truth_range[irange]->Rebin(nbins - first_bin, Form("h_rebin_truth_%d", irange), &dxj_bins[first_bin]);
      TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);

      ht->Draw("p E1");
      g_final_xj_systematics[irange][niter]->Draw("same p E2");
      hu->Draw("same p");


      h_final_xj_data_range[irange]->Draw("same p");
      h_final_xj_reco_range[irange]->Draw("same p");
      
      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);

      TLegend *leg = new TLegend(0.2, 0.4, 0.4, 0.56);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_truth_range[irange], "PYTHIA-8","p");
      leg->AddEntry(h_final_xj_reco_range[irange], "PYTHIA-8 Reco","p");
      leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Data Unfold","p");
      leg->AddEntry(h_final_xj_data_range[irange], "Data Reco","p");
      leg->Draw("same");
   
      cxj->cd(2);

      TH1D *h_data_compare = (TH1D*) h_final_xj_truth_range[irange]->Clone();
      h_data_compare->Divide(h_final_xj_unfold_range[irange][niter]);
      h_data_compare->SetTitle(";x_{J}; PYTHIA-8 / Unfold");
      dlutility::SetFont(h_data_compare, 42, 0.1, 0.07, 0.07, 0.07);
      dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
      dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);
      h_data_compare->SetMaximum(2.3);
      h_data_compare->SetMinimum(0.0);
      TH1D *hd = (TH1D*) h_data_compare->Rebin(nbins - first_bin, "h_rebin_compare", &dxj_bins[first_bin]);

      hd->Draw("p");

      // no the systematics

      TGraphAsymmErrors *g_compare = new TGraphAsymmErrors(h_data_compare);
      for (int ib = 0; ib < h_data_compare->GetNbinsX(); ib++)
	{
	  g_compare->SetPointError(ib, 0, 0, g_final_xj_systematics[irange][niter]->GetErrorYlow(ib)/h_final_xj_truth_range[irange]->GetBinContent(ib+1), g_final_xj_systematics[irange][niter]->GetErrorYhigh(ib)/h_final_xj_truth_range[irange]->GetBinContent(ib+1));
	  g_compare->SetPointY(ib, 1);
	}
      dlutility::SetLineAtt(g_compare, kBlack, 1,1);
      dlutility::SetMarkerAtt(g_compare, kBlack, 1,1);
      g_compare->SetFillColorAlpha(kBlack, 0.3);
      g_compare->Draw("p E2 same");
      TLine *line = new TLine(hd->GetBinLowEdge(1), 1, 1, 1);
      line->SetLineStyle(4);
      line->SetLineColor(kRed + 3);
      line->SetLineWidth(2);
      line->Draw("same");
      cxj->Print(Form("%s/final_plots/h_xj_unfolded_pp_r%02d_range_%d.png",  rb.get_code_location().c_str(), cone_size, irange));
      cxj->Print(Form("%s/final_plots/h_xj_unfolded_pp_r%02d_range_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange));
    }

  
  TCanvas *cxj_money = new TCanvas("cxj_money","cxj_money", 500, 500);

  TH1D *hblank = new TH1D("hblank","", 1, 0.25, 1);
  
  for (int irange = 0; irange < mbins; irange++)
    {

      dlutility::SetLineAtt(g_final_xj_systematics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_systematics[irange][niter], color_unfold, msize_unfold, marker_unfold);
      g_final_xj_systematics[irange][niter]->SetFillColorAlpha(color_unfold_fill, 0.3); 

      dlutility::SetLineAtt(g_final_xj_statistics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_statistics[irange][niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(g_final_xj_truth[irange], color_pythia, lsize_pythia, 1);
      dlutility::SetMarkerAtt(g_final_xj_truth[irange], color_pythia, msize_pythia, marker_pythia);
      g_final_xj_truth[irange]->SetFillStyle(0);
      
      dlutility::SetLineAtt(g_final_xj_herwig[irange], color_herwig, lsize_herwig, 1);
      dlutility::SetMarkerAtt(g_final_xj_herwig[irange], color_herwig, msize_herwig, marker_herwig);
      g_final_xj_herwig[irange]->SetFillStyle(0);

      dlutility::SetLineAtt(h_final_xj_data_range[irange], color_data, lsize_data, 1);
      dlutility::SetMarkerAtt(h_final_xj_data_range[irange], color_data, msize_data, marker_data);
      dlutility::SetLineAtt(h_final_xj_reco_range[irange], color_reco, lsize_reco, 1);
      dlutility::SetMarkerAtt(h_final_xj_reco_range[irange], color_reco, msize_reco, marker_reco);

      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.17);

      
      hblank->GetYaxis()->SetTitleOffset(1.8);
      dlutility::SetFont(hblank, 42, 0.06, 0.04, 0.05, 0.05);
      hblank->SetMaximum(5);
      hblank->SetMinimum(0);

      hblank->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;
      hblank->Draw();
      g_final_xj_truth[irange]->Draw("l");
      g_final_xj_herwig[irange]->Draw("l");

      g_final_xj_systematics[irange][niter]->Draw("same p E2");
      g_final_xj_statistics[irange][niter]->Draw("same p E1");

      float top = 0.88;
      float ss = 0.05;
      dlutility::DrawSPHENIXpp(0.22, top);

      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, top - 2*ss);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, top - 3*ss);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, top - 4*ss);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, top - 5*ss);

      
      TLegend *leg = new TLegend(0.22, top - 8.5*ss, 0.4, top - 5.8*ss);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(g_final_xj_systematics[irange][niter], "Data");
      leg->AddEntry(g_final_xj_truth[irange], "PYTHIA-8","l");
      leg->AddEntry(g_final_xj_herwig[irange], "HERWIG","l");
      //leg->AddEntry(h_linear_herwig_xj[irange], "HERWIG 7.3","l");
      leg->Draw("same");

      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_%d.png",  rb.get_code_location().c_str(), cone_size, irange));
      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange));
    }

  TCanvas *cxj_moneyall = new TCanvas("cxj_moneyall","cxj_moneyall", 500, 500);

  TGraphAsymmErrors *g_final_shift_xj_unfold_range[mbins];
  TGraphAsymmErrors *g_final_shift_xj_sys_range[mbins];
  TGraphAsymmErrors *g_final_shift_xj_truth_range[mbins];

  //////////
  for (int i = 0; i < mbins; i++)
    {

      g_final_shift_xj_sys_range[i] = (TGraphAsymmErrors*) g_final_xj_systematics[i][niter]->Clone();
      g_final_shift_xj_unfold_range[i] = (TGraphAsymmErrors*) g_final_xj_statistics[i][niter]->Clone();
      g_final_shift_xj_truth_range[i] = (TGraphAsymmErrors*) g_final_xj_truth[i]->Clone();

      
      double binwidth = ixj_bins[first_bin + 1] - ixj_bins[first_bin];
      double fix_shift = binwidth / (float) mbins;

      for (int b = 0; b < g_final_shift_xj_sys_range[i]->GetN(); b++) {
	double ex = fix_shift/2.;
	
	double x = g_final_shift_xj_unfold_range[i]->GetPointX(b);
	double y = g_final_shift_xj_unfold_range[i]->GetPointY(b);
	double eyh = g_final_shift_xj_unfold_range[i]->GetErrorYhigh(b);
	double eyl = g_final_shift_xj_unfold_range[i]->GetErrorYlow(b);

	double xs = g_final_shift_xj_sys_range[i]->GetPointX(b);
	double ys = g_final_shift_xj_sys_range[i]->GetPointY(b);
	double eyhs = g_final_shift_xj_sys_range[i]->GetErrorYhigh(b);
	double eyls = g_final_shift_xj_sys_range[i]->GetErrorYlow(b);

	std::cout << eyhs << " / " << eyls << std::endl;
	double xt = g_final_shift_xj_truth_range[i]->GetPointX(b);
	double yt = g_final_shift_xj_truth_range[i]->GetPointY(b);
	double eyht = g_final_shift_xj_truth_range[i]->GetErrorYhigh(b);
	double eylt = g_final_shift_xj_truth_range[i]->GetErrorYlow(b);

	double sx=x + (i-2)*fix_shift;

	g_final_shift_xj_unfold_range[i]->SetPointError(b, 0,0, eyl, eyh);	    
	g_final_shift_xj_unfold_range[i]->SetPoint(b, sx, y);
	g_final_shift_xj_sys_range[i]->SetPointError(b, ex,ex, eyls, eyhs);	    
	g_final_shift_xj_sys_range[i]->SetPoint(b, sx, ys);
	g_final_shift_xj_truth_range[i]->SetPointError(b, 0,0, eylt, eyht);	    
	g_final_shift_xj_truth_range[i]->SetPoint(b, sx, yt);

      }
    }

  /////////
  
  //      h_linear_herwig_xj[irange]->Scale(1./h_linear_herwig_xj[irange]->Integral(0, -1, "width"));

  // dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
  // dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

  int marker_range[3] = {20, 21, 22};
  float msize_range[3] = {0.8, 0.7, 0.9};
    for (int irange = 0; irange < mbins; irange++)
    {
      dlutility::SetLineAtt(g_final_shift_xj_sys_range[irange], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_shift_xj_sys_range[irange], color_unfold, msize_range[irange], marker_range[irange]);
      g_final_shift_xj_sys_range[irange]->SetFillColorAlpha(color_unfold, 0.3); 

      dlutility::SetLineAtt(g_final_shift_xj_unfold_range[irange], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_shift_xj_unfold_range[irange], color_unfold, msize_range[irange], marker_range[irange]);

      dlutility::SetLineAtt(g_final_shift_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
      dlutility::SetMarkerAtt(g_final_shift_xj_truth_range[irange], color_pythia, msize_range[irange], marker_range[irange]);


    }
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.17);
  hblank->SetMaximum(5);
  hblank->Draw("");
  for (int irange = 0; irange < mbins; irange++)
    {
      g_final_shift_xj_truth_range[irange]->Draw("same p E1");
      g_final_shift_xj_sys_range[irange]->Draw("same p E2");
      g_final_shift_xj_unfold_range[irange]->Draw("same p E1");

    }
  float top = 0.88;
  float ss = 0.05;
  dlutility::DrawSPHENIXpp(0.2, top);

  dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.2, top - 2*ss);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.2, top - 3*ss);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.2, top - 4*ss);

      
  TLegend *legd = new TLegend(0.5, 0.78, 0.70, 0.91);
  legd->SetHeader("Data");
		      
  legd->SetLineWidth(0);
  legd->SetTextSize(0.03);
  legd->SetTextFont(42);
  legd->AddEntry(g_final_shift_xj_sys_range[0], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[0]], ipt_bins[measure_bins[0+1]]));
  legd->AddEntry(g_final_shift_xj_sys_range[1], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[1]], ipt_bins[measure_bins[1+1]]));
  legd->AddEntry(g_final_shift_xj_sys_range[2], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[2]], ipt_bins[measure_bins[2+1]]));
  legd->Draw("same");
  TLegend *legs = new TLegend(0.5, 0.63, 0.70, 0.76);
  legs->SetHeader("PYTHIA-8");
  legs->SetLineWidth(0);
  legs->SetTextSize(0.03);
  legs->SetTextFont(42);
      
  legs->AddEntry(g_final_shift_xj_truth_range[0], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[0]], ipt_bins[measure_bins[0+1]]));
  legs->AddEntry(g_final_shift_xj_truth_range[1], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[1]], ipt_bins[measure_bins[1+1]]));
  legs->AddEntry(g_final_shift_xj_truth_range[2], Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[2]], ipt_bins[measure_bins[2+1]]));
  legs->Draw("same");
      

  cxj_moneyall->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_all.png",  rb.get_code_location().c_str(), cone_size));
  cxj_moneyall->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_all.pdf",  rb.get_code_location().c_str(), cone_size));


  return;
}
