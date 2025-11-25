#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
const bool NUCLEAR = true;
const int color_unfold_fill = kAzure - 4;
const int color_unfold = kAzure - 5;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

const int color_pythia = kRed;
const float msize_pythia = 0.9;
const float marker_pythia = 20;
const float lsize_pythia = 1.1;

const int color_herwig = kViolet;
const float msize_herwig = 0.9;
const float marker_herwig = 20;
const float lsize_herwig = 1.1;

const int color_reco = kRed;
const float marker_reco = 24;
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
      if (dxj_bins[i] > 0.3 && first_bin == 0) first_bin = i;

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
  

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TFile *fintr = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_PRIMER1_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fintr)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fintr->Get("h_truth_flat_pt1pt2");


  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  /* TFile *finhe = new TFile("truth_hists/herwig_hist.root","r"); */

  /* TH1D *h_flat_herwig_pt1pt2 = (TH1D*) finhe->Get("h_truth_flat_pt1pt2"); */
  /* h_flat_herwig_pt1pt2->SetName("h_flat_herwig_pt1pt2"); */

  /* TH1D *h_linear_herwig_xj[3]; */
  /* for (int i = 0; i < 3; i++) */
  /*   { */
  /*     h_linear_herwig_xj[i] = (TH1D*) finhe->Get(Form("h_linear_truth_xj_%d", i)); */
  /*     h_linear_herwig_xj[i]->SetName(Form("h_linear_herwig_xj_%d", i)); */

  /*     if (!h_linear_herwig_xj[i]) */
  /* 	{ */
  /* 	  std::cout << "AHH"<< std::endl; */
  /* 	  return; */
  /* 	} */
  /*   } */
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  
  const int niterations = 10;
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finu)
    {
      std::cout << " no unc " << std::endl;
      return;
    }
  TProfile *h_xj_rms[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  h_xj_rms[irange][iter] = (TProfile*) finu->Get(Form("hp_xj_range_%d_%d", irange, iter));
	}
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;

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
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  //TH2D *h_pt1pt2_herwig = new TH2D("h_pt1pt2_herwig", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  //TH1D *h_xj_herwig = new TH1D("h_xj_herwig", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  //  TH1D *h_xj_herwig_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_data_range[irange] = new TH1D(Form("h_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //h_xj_herwig_range[irange] = new TH1D(Form("h_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
    if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data_subrange[mbins];
  TH1D *h_xj_reco_subrange[mbins];
  TH1D *h_xj_truth_subrange[mbins];
  //  TH1D *h_xj_herwig_subrange[mbins];
  TH1D *h_xj_unfold_subrange[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_data_subrange[irange] = new TH1D(Form("h_xj_data_subrange_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_subrange[irange] = new TH1D(Form("h_xj_reco_subrange_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_subrange[irange] = new TH1D(Form("h_xj_truth_subrange_%d", irange), ";x_{J};", nbins, ixj_bins);
      //      h_xj_herwig_subrange[irange] = new TH1D(Form("h_xj_herwig_subrange_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_subrange[irange][iter] = new TH1D(Form("h_xj_unfold_subrange_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  TH1D *h_xjunc_data = new TH1D("h_xjunc_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_truth = new TH1D("h_xjunc_truth", ";x_{J};",nbins, ixj_bins);
  //  TH1D *h_xjunc_herwig = new TH1D("h_xjunc_herwig", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xjunc_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_xjunc_unfold[iter] = new TH1D(Form("h_xjunc_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xjunc_data_range[mbins];
  TH1D *h_xjunc_reco_range[mbins];
  TH1D *h_xjunc_truth_range[mbins];
  //  TH1D *h_xjunc_herwig_range[mbins];
  TH1D *h_xjunc_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xjunc_data_range[irange] = new TH1D(Form("h_xjunc_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_reco_range[irange] = new TH1D(Form("h_xjunc_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_truth_range[irange] = new TH1D(Form("h_xjunc_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //      h_xjunc_herwig_range[irange] = new TH1D(Form("h_xjunc_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xjunc_unfold_range[irange][iter] = new TH1D(Form("h_xjunc_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  //  histo_opps::make_sym_pt1pt2(h_flat_herwig_pt1pt2, h_pt1pt2_herwig, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  //histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_truth, nbins);
  //  histo_opps::normalize_histo(h_xj_herwig, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_final_xj_data_range[mbins];
  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  //  TH1D *h_final_xj_herwig_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];
  TH1D *h_final_xj_systematics[mbins][niterations];
  TH1D *h_final_xj_statistics[mbins][niterations];

  TGraphAsymmErrors *g_final_xj_systematics[mbins][niterations];
  TGraphAsymmErrors *g_final_xj_statistics[mbins][niterations];
  TGraphAsymmErrors *g_final_xj_truth[mbins];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_data_range[irange] = new TH1D(Form("h_final_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //      h_final_xj_herwig_range[irange] = new TH1D(Form("h_final_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
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
      //      histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_truth_range[irange], nbins);
      //      histo_opps::normalize_histo(h_xj_herwig_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_data_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold_range[irange][iter], nbins);
	  histo_opps::normalize_histo(h_xj_rms[irange][iter], nbins);
	}


      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      g_final_xj_truth[irange] = histo_opps::get_xj_statistics(h_final_xj_truth_range[irange], nbins);
      //      histo_opps::finalize_xj(h_xj_herwig_range[irange], h_final_xj_herwig_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_data_range[irange], h_final_xj_data_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	  histo_opps::set_xj_errors(h_final_xj_unfold_range[irange][iter], h_xj_rms[irange][iter], nbins);

	  h_final_xj_systematics[irange][iter] = (TH1D*) h_final_xj_unfold_range[irange][iter]->Clone();
	  h_final_xj_systematics[irange][iter]->SetName(Form("h_final_xj_systematics_%d_%d", irange, iter));
	  g_final_xj_systematics[irange][iter] = histo_opps::get_xj_systematics(h_final_xj_systematics[irange][iter], h_total_sys_neg_range[irange][iter], h_total_sys_range[irange][iter], nbins);

	  h_final_xj_statistics[irange][iter] = (TH1D*) h_final_xj_unfold_range[irange][iter]->Clone();
	  h_final_xj_statistics[irange][iter]->SetName(Form("h_final_xj_statistics_%d_%d", irange, iter));
	  g_final_xj_statistics[irange][iter] = histo_opps::get_xj_statistics(h_final_xj_statistics[irange][iter], nbins);
	  


	}

    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  // Subleading
  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_data, h_xj_data_subrange[irange], nbins, measure_bins[1], measure_bins[2],subleading_measure_bins[irange], subleading_measure_bins[irange+1]);
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_subrange[irange], nbins, measure_bins[1], measure_bins[2], subleading_measure_bins[irange], subleading_measure_bins[irange+1]);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_subrange[irange], nbins, measure_bins[1], measure_bins[2], subleading_measure_bins[irange], subleading_measure_bins[irange+1]);
      //      histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig_subrange[irange], nbins, measure_bins[1], measure_bins[2], subleading_measure_bins[irange], subleading_measure_bins[irange+1]);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_subrange[irange][iter], nbins, measure_bins[1], measure_bins[2], subleading_measure_bins[irange], subleading_measure_bins[irange+1]);
	}


      histo_opps::normalize_histo(h_xj_truth_subrange[irange], nbins);
      //      histo_opps::normalize_histo(h_xj_herwig_subrange[irange], nbins);
      histo_opps::normalize_histo(h_xj_data_subrange[irange], nbins);
      histo_opps::normalize_histo(h_xj_reco_subrange[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold_subrange[irange][iter], nbins);
	}

    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
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

      //      h_final_xj_truth_range[irange]->SetFillColorAlpha(color_pythia, 0.3);

      dlutility::SetLineAtt(h_final_xj_data_range[irange], color_data, lsize_data, 1);
      dlutility::SetMarkerAtt(h_final_xj_data_range[irange], color_data, msize_data, marker_data);
      dlutility::SetLineAtt(h_final_xj_reco_range[irange], color_reco, lsize_reco, 1);
      dlutility::SetMarkerAtt(h_final_xj_reco_range[irange], color_reco, msize_reco, marker_reco);

      dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.05);
      h_final_xj_truth_range[irange]->SetMaximum(5);
      h_final_xj_truth_range[irange]->SetMinimum(0);
      h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      //h_final_xj_truth_range[irange]->Draw("E4 same");
      TH1D *ht = (TH1D*) h_final_xj_truth_range[irange]->Rebin(nbins - first_bin, Form("h_rebin_truth_%d", irange), &dxj_bins[first_bin]);
      //TH1D *hs = (TH1D*) h_final_xj_systematics[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_sys_%d", irange), &dxj_bins[first_bin]);
      TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);

      ht->Draw("p E1");
      g_final_xj_systematics[irange][niter]->Draw("same p E2");
      hu->Draw("same p");

      /* h_final_xj_truth_range[irange]->Draw("p E1"); */
      /* //h_final_xj_truth_range[irange]->Draw("E4 same"); */
      /* h_final_xj_systematics[irange][niter]->Draw("same p E2"); */
      /* h_final_xj_unfold_range[irange][niter]->Draw("same p E1"); */

      h_final_xj_data_range[irange]->Draw("same p");
      h_final_xj_reco_range[irange]->Draw("same p");

      
      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);
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
	  g_compare->SetPointY(ib, 1);//h_data_compare->GetBinWidth(ib+1)/2., h_data_compare->GetBinWidth(ib+1)/2., h_final_xj_systematics[irange][niter]->GetBinError(ib+1)/h_final_xj_truth_range[irange]->GetBinContent(ib+1), h_final_xj_systematics[irange][niter]->GetBinError(ib+1)/h_final_xj_truth_range[irange]->GetBinContent(ib+1));
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

  for (int irange = 0; irange < mbins; irange++)
    {

      //      h_linear_herwig_xj[irange]->Scale(1./h_linear_herwig_xj[irange]->Integral(0, -1, "width"));

      // dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
      // dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(g_final_xj_systematics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_systematics[irange][niter], color_unfold, msize_unfold, marker_unfold);
      g_final_xj_systematics[irange][niter]->SetFillColorAlpha(color_unfold_fill, 0.3); 

      dlutility::SetLineAtt(g_final_xj_statistics[irange][niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_statistics[irange][niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(g_final_xj_truth[irange], color_pythia, lsize_pythia, 1);
      dlutility::SetMarkerAtt(g_final_xj_truth[irange], color_pythia, msize_pythia, marker_pythia);




      //      dlutility::SetLineAtt(h_linear_herwig_xj[irange], color_herwig, 3, 1);
      //dlutility::SetMarkerAtt(h_linear_herwig_xj[irange], color_herwig, msize_herwig, marker_herwig);

      //      h_final_xj_truth_range[irange]->SetFillColorAlpha(color_pythia, 0.3);

      dlutility::SetLineAtt(h_final_xj_data_range[irange], color_data, lsize_data, 1);
      dlutility::SetMarkerAtt(h_final_xj_data_range[irange], color_data, msize_data, marker_data);
      dlutility::SetLineAtt(h_final_xj_reco_range[irange], color_reco, lsize_reco, 1);
      dlutility::SetMarkerAtt(h_final_xj_reco_range[irange], color_reco, msize_reco, marker_reco);

      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.17);

      TH1D *hblank = (TH1D*) h_final_xj_truth_range[irange]->Clone();

      hblank->Reset();

      TH1D *ht = (TH1D*) hblank->Rebin(nbins - first_bin, Form("h_rebin_blank_%d", irange), &dxj_bins[first_bin]);

      ht->GetYaxis()->SetTitleOffset(1.8);
      dlutility::SetFont(ht, 42, 0.06, 0.04, 0.05, 0.05);
      ht->SetMaximum(5);
      ht->SetMinimum(0);
      ht->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;


      ht->Draw("");

      g_final_xj_systematics[irange][niter]->Draw("same p E2");
      g_final_xj_statistics[irange][niter]->Draw("same p E1");
      g_final_xj_truth[irange]->Draw("same p E1");

      float top = 0.88;
      float ss = 0.05;
      dlutility::DrawSPHENIXppPrelim(0.22, top);

      dlutility::drawText(Form("anti-#it{k}_{t} #kern[-0.1]{#it{R} = %0.1f}", cone_size*0.1), 0.22, top - 2*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq #it{p}_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, top - 3*ss);
      dlutility::drawText(Form("#it{p}_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, top - 4*ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, top - 5*ss);

      
      TLegend *leg = new TLegend(0.22, top - 8.5*ss, 0.4, top - 5.8*ss);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(g_final_xj_systematics[irange][niter], "Data");
      leg->AddEntry(g_final_xj_truth[irange], "PYTHIA-8","l");
      //leg->AddEntry(h_linear_herwig_xj[irange], "HERWIG 7.3","l");
      leg->Draw("same");

      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_%d.png",  rb.get_code_location().c_str(), cone_size, irange));
      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_pp_r%02d_range_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange));
    }

  return;
}
