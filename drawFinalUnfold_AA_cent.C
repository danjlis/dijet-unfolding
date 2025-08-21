#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

const bool NUCLEAR = true;

const int color_pp_unfold_fill = kBlack;
const int color_pp_unfold = kBlack;
const float marker_pp_unfold = 20;
const float msize_pp_unfold = 0.9;
const float lsize_pp_unfold = 1.1;


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
const int color_pp_data = kAzure - 6;
const float marker_pp_data = 24;
const float msize_pp_data = 0.9;
const float lsize_pp_data = 1.1;

void drawFinalUnfold_AA_cent(const int cone_size = 3, const int centrality_bin = 0, const std::string configfile = "binning_AA.config")
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();


  int color_unfold_fill[4] = {kRed + 1, kOrange, kGreen + 2, kAzure - 4};
  int color_unfold[4] = {kRed + 1, kAzure + 2, kGreen+3, kViolet -1};
  float marker_unfold = 20;
  float msize_unfold = 0.9;
  float lsize_unfold = 2.2;

  
  read_binning rb(configfile.c_str());

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  std::string dphi_string = rb.get_dphi_string();

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins+1];
  rb.get_centrality_bins(icentrality_bins);

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
      if (dxj_bins[i] > 0.2 && first_bin == 0) first_bin = i;

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
  

  TFile *finjw = new TFile(Form("%s/truth_hists/JEWEL_for_Dan.root",  rb.get_code_location().c_str()),"r");
  if (!finjw)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }
  
  TH1D *h_jewel_xj_med = (TH1D*) finjw->Get("h1_medium");
  TH1D *h_jewel_xj_vac = (TH1D*) finjw->Get("h1_vacuum");
  dlutility::SetLineAtt(h_jewel_xj_med, kGreen + 2, 2, 1);
  dlutility::SetLineAtt(h_jewel_xj_vac, kBlack, 2, 1);
  
  TFile *fintr = new TFile(Form("%s/truth_hists/truth_hist_r%02d.root",  rb.get_code_location().c_str(),  cone_size),"r");
  if (!fintr)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }
  
  TH1D *h_linear_truth_xj[3];
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i] = (TH1D*) fintr->Get(Form("h_linear_truth_xj_%d", i));
      h_linear_truth_xj[i]->SetName(Form("h_linear_pythia_xj_%d", i));
      if (!h_linear_truth_xj[i])
	{
	  std::cout << "AHH"<< std::endl;
	  return;
	}
    }

  
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
  const int niter = 2;
  //pp results
  TProfile *hp_xj_rms_pp[mbins][niterations];
  TFile *finupp = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finupp)
    {
      std::cout << " no unc " << std::endl;
      return;
    }
  
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  hp_xj_rms_pp[irange][iter] = (TProfile*) finupp->Get(Form("hp_xj_range_%d_%d", irange, iter));
	}
    }

  TProfile *hp_xj_rms[cent_bins][mbins][niterations];

  TFile *finu[cent_bins];
  for (int i = 0; i < cent_bins; i++)
    {
      finu[i] = new TFile(Form("%s/uncertainties/uncertainties_AA_cent_%d_r%02d_nominal.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!finu[i])
	{
	  std::cout << " no unc " << std::endl;
	  return;
	}
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  for (int iter = 0; iter < niterations; ++iter)
	    {
	      hp_xj_rms[i][irange][iter] = (TProfile*) finu[i]->Get(Form("hp_xj_range_%d_%d", irange, iter));
	      hp_xj_rms[i][irange][iter]->SetName(Form("hp_xj_rms_%d_%d_%d", i, irange, iter));
	    }
	}
    }
  TH1D *h_total_sys_range[cent_bins][mbins][niterations];
  TH1D *h_total_sys_neg_range[cent_bins][mbins][niterations];
  
  TFile *fins[cent_bins];
  for (int i = 0; i < cent_bins; i++)
    {
      fins[i] = new TFile(Form("%s/uncertainties/systematics_AA_cent_%d_r%02d.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!fins[i])
	{
	  std::cout << " no sys fin " << i << std::endl;
	  return;
	}
      for (int irange = 0; irange < mbins; irange++)
	{
	  for (int iter = 0; iter < niterations; ++iter)
	    {
	      h_total_sys_range[i][irange][iter] = (TH1D*) fins[i]->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	      h_total_sys_neg_range[i][irange][iter] = (TH1D*) fins[i]->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));
	      
	    }    
	}
    }
  TH1D *h_total_pp_sys_range[mbins][niterations];
  TH1D *h_total_pp_sys_neg_range[mbins][niterations];
  
  TFile *fins_pp = new TFile(Form("%s/uncertainties/systematics_pp_r%02d.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fins_pp)
    {
      std::cout << " no sys pp"  << std::endl;
      return;
    }
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  h_total_pp_sys_range[irange][iter] = (TH1D*) fins_pp->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	  h_total_pp_sys_neg_range[irange][iter] = (TH1D*) fins_pp->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));
	  h_total_pp_sys_range[irange][iter]->SetName(Form("h_total_pp_sys_range_%d_%d", irange, iter));
	  h_total_pp_sys_neg_range[irange][iter]->SetName(Form("h_total_pp_sys_neg_range_%d_%d", irange, iter));
	}    
    }


  TFile *fin[cent_bins];
  TH1D *h_flat_data_pt1pt2[cent_bins];
  TH1D *h_flat_unfold_pt1pt2[cent_bins][niterations];

  for (int i = 0; i <  cent_bins; i++)
    {
      fin[i]= new TFile(Form("%s/unfolding_hists/unfolding_hists_AA_cent_%d_r%02d_nominal.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!fin[i])
	{
	  std::cout << " no file " << std::endl;
	  return;
	}
      h_flat_data_pt1pt2[i] = (TH1D*) fin[i]->Get("h_data_flat_pt1pt2");
      h_flat_data_pt1pt2[i]->SetName(Form("h_flat_data_pt1pt2_%d", i));
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_flat_unfold_pt1pt2[i][iter] = (TH1D*) fin[i]->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_flat_unfold_pt1pt2[i][iter]->SetName(Form("h_flat_unfold_pt1pt2_%d_%d", i, iter));
	}
    }


  TFile *finpp = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finpp)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  TH1D *h_flat_pp_data_pt1pt2 = (TH1D*) finpp->Get("h_data_flat_pt1pt2");
  h_flat_pp_data_pt1pt2->SetName("h_flat_pp_data_pt1pt2");
  
  TH1D *h_flat_pp_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_pp_unfold_pt1pt2[iter] = (TH1D*) finpp->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_pp_unfold_pt1pt2[iter]->SetName(Form("h_flat_pp_unfold_pt1pt2_%d", iter));
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  TH2D *h_pt1pt2_data[cent_bins];
  TH2D *h_pt1pt2_pp_data = new TH2D("h_pt1pt2_pp_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  //TH2D *h_pt1pt2_herwig = new TH2D("h_pt1pt2_herwig", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[cent_bins][niterations];
  TH2D *h_pt1pt2_pp_unfold[niterations];
  for (int i = 0; i < cent_bins; i++)
    {
      h_pt1pt2_data[i] = new TH2D(Form("h_pt1pt2_data_%d", i), ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int i = 0; i < cent_bins; i++)
	{
	  h_pt1pt2_unfold[i][iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
	  h_pt1pt2_unfold[i][iter]->SetName(Form("h_pt1pt2_unfold_%d_iter%d", i, iter));
	}
      h_pt1pt2_pp_unfold[iter] = new TH2D("h_pt1pt2_pp_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_pp_unfold[iter]->SetName(Form("h_pt1pt2_pp_unfold_iter%d", iter));

    }

  TH1D *h_xj_data[cent_bins];
  TH1D *h_xj_pp_data = new TH1D("h_xj_pp_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  //TH1D *h_xj_herwig = new TH1D("h_xj_herwig", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[cent_bins][niterations];
  TH1D *h_xj_pp_unfold[niterations];
  for (int i = 0; i < cent_bins; i++)
    {
      h_xj_data[i]= new TH1D(Form("h_xj_data_%d", i), ";x_{J};", nbins, ixj_bins);
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int i = 0; i < cent_bins; i++)
	{
	  h_xj_unfold[i][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", i, iter), ";x_{J};",nbins, ixj_bins);
	}
      h_xj_pp_unfold[iter] = new TH1D(Form("h_xj_pp_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data_range[cent_bins][mbins];
  TH1D *h_xj_pp_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  //  TH1D *h_xj_herwig_range[mbins];
  TH1D *h_xj_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_xj_pp_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {

      h_xj_pp_data_range[irange] = new TH1D(Form("h_xj_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //h_xj_herwig_range[irange] = new TH1D(Form("h_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0 ; i < cent_bins; i++)
	{
	  h_xj_data_range[i][irange] = new TH1D(Form("h_xj_data_range_%d_%d",i, irange), ";x_{J};", nbins, ixj_bins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_xj_unfold_range[i][irange][iter] = new TH1D(Form("h_xj_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_pp_unfold_range[irange][iter] = new TH1D(Form("h_xj_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  TH1D *h_xjunc_data[cent_bins];
  TH1D *h_xjunc_pp_data = new TH1D("h_xjunc_pp_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_truth = new TH1D("h_xjunc_truth", ";x_{J};",nbins, ixj_bins);
  //  TH1D *h_xjunc_herwig = new TH1D("h_xjunc_herwig", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xjunc_unfold[cent_bins][niterations];
  TH1D *h_xjunc_pp_unfold[niterations];
  for (int i = 0 ; i < cent_bins; i++)
    {
      h_xjunc_data[i]= new TH1D(Form("h_xjunc_data_%d", i), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xjunc_unfold[i][iter] = new TH1D(Form("h_xjunc_unfold_%d_iter%d", i, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xjunc_pp_unfold[iter] = new TH1D(Form("h_xjunc_pp_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xjunc_data_range[cent_bins][mbins];
  TH1D *h_xjunc_pp_data_range[mbins];
  TH1D *h_xjunc_reco_range[mbins];
  TH1D *h_xjunc_truth_range[mbins];
  //  TH1D *h_xjunc_herwig_range[mbins];
  TH1D *h_xjunc_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_xjunc_pp_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {

      h_xjunc_pp_data_range[irange] = new TH1D(Form("h_xjunc_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_reco_range[irange] = new TH1D(Form("h_xjunc_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_truth_range[irange] = new TH1D(Form("h_xjunc_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //      h_xjunc_herwig_range[irange] = new TH1D(Form("h_xjunc_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0; i < cent_bins; i++)
	{
	  h_xjunc_data_range[i][irange] = new TH1D(Form("h_xjunc_data_range_%d_%d", i, irange), ";x_{J};", nbins, ixj_bins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_xjunc_unfold_range[i][irange][iter] = new TH1D(Form("h_xjunc_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	  for (int iter = 0; iter < niterations; iter++)
	    h_xjunc_pp_unfold_range[irange][iter] = new TH1D(Form("h_xjunc_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
    
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  //  histo_opps::make_sym_pt1pt2(h_flat_herwig_pt1pt2, h_pt1pt2_herwig, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_pp_data_pt1pt2, h_pt1pt2_pp_data, nbins);
  for (int i = 0; i < cent_bins; i++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2[i], h_pt1pt2_data[i], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[i][iter], h_pt1pt2_unfold[i][iter], nbins);
	}

    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_pp_unfold_pt1pt2[iter], h_pt1pt2_pp_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  histo_opps::project_xj(h_pt1pt2_pp_data, h_xj_pp_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  //histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);

  for (int i = 0; i < cent_bins; i++)
    {
      histo_opps::project_xj(h_pt1pt2_data[i], h_xj_data[i], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[i][iter], h_xj_unfold[i][iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_pp_unfold[iter], h_xj_pp_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_pp_data, nbins);
  for (int i = 0 ; i < cent_bins; i++)
    {
      histo_opps::normalize_histo(h_xj_data[i], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold[i][iter], nbins);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_pp_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  TH1D *h_final_xj_pp_data_range[mbins];
  TH1D *h_final_xj_pp_unfold_range[mbins][niterations];

  TH1D *h_final_xj_pp_systematics[mbins][niterations];
  TGraphAsymmErrors *g_final_xj_pp_systematics[mbins][niterations];

  TH1D *h_final_xj_data_range[cent_bins][mbins];
  TH1D *h_final_xj_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_final_xj_systematics[cent_bins][mbins][niterations];
  TGraphAsymmErrors *g_final_xj_systematics[cent_bins][mbins][niterations];
  
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_pp_data_range[irange] = new TH1D(Form("h_final_xj_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      //      h_final_xj_herwig_range[irange] = new TH1D(Form("h_final_xj_herwig_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0; i < cent_bins; i++)
	{
	  h_final_xj_data_range[i][irange] = new TH1D(Form("h_final_xj_data_range_%d_%d", i, irange), ";x_{J};", nbins, ixj_bins);

	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_final_xj_unfold_range[i][irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	      h_final_xj_systematics[i][irange][iter] = new TH1D(Form("h_final_xj_systematics_%d_%d_%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  
	  h_final_xj_pp_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_final_xj_pp_systematics[irange][iter] = new TH1D(Form("h_final_xj_pp_systematics_%d_%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int irange = 0; irange < mbins; irange++)
    {

      histo_opps::project_xj(h_pt1pt2_pp_data, h_xj_pp_data_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);

      //      histo_opps::project_xj(h_pt1pt2_herwig, h_xj_herwig_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int i = 0; i < cent_bins; i++)
	{
	  histo_opps::project_xj(h_pt1pt2_data[i], h_xj_data_range[i][irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::project_xj(h_pt1pt2_unfold[i][iter], h_xj_unfold_range[i][irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_pp_unfold[iter], h_xj_pp_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_pp_data_range[irange], nbins);
      for (int i = 0 ; i < cent_bins; i++)
	{
	  histo_opps::normalize_histo(h_xj_data_range[i][irange], nbins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::normalize_histo(h_xj_unfold_range[i][irange][iter], nbins);
	      histo_opps::normalize_histo(hp_xj_rms[i][irange][iter], nbins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_pp_unfold_range[irange][iter], nbins);
	  histo_opps::normalize_histo(hp_xj_rms_pp[irange][iter], nbins);
	}


      histo_opps::finalize_xj(h_xj_pp_data_range[irange], h_final_xj_pp_data_range[irange], nbins, first_xj);
      for (int i = 0; i < cent_bins; i++)
	{
	  histo_opps::finalize_xj(h_xj_data_range[i][irange], h_final_xj_data_range[i][irange], nbins, first_xj);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::finalize_xj(h_xj_unfold_range[i][irange][iter], h_final_xj_unfold_range[i][irange][iter], nbins, first_xj);
	      histo_opps::set_xj_errors(h_final_xj_unfold_range[i][irange][iter], hp_xj_rms[i][irange][iter], nbins);


	      h_final_xj_systematics[i][irange][iter] = (TH1D*) h_final_xj_unfold_range[i][irange][iter]->Clone();
	      h_final_xj_systematics[i][irange][iter]->SetName(Form("h_final_xj_systematics_%d_%d_%d", i, irange, iter));
	      g_final_xj_systematics[i][irange][iter] = histo_opps::get_xj_systematics(h_final_xj_systematics[i][irange][iter], h_total_sys_neg_range[i][irange][iter], h_total_sys_range[i][irange][iter], nbins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_pp_unfold_range[irange][iter], h_final_xj_pp_unfold_range[irange][iter], nbins, first_xj);
	  histo_opps::set_xj_errors(h_final_xj_pp_unfold_range[irange][iter], hp_xj_rms_pp[irange][iter], nbins);
	      
	  h_final_xj_pp_systematics[irange][iter] = (TH1D*) h_final_xj_pp_unfold_range[irange][iter]->Clone();
	  h_final_xj_pp_systematics[irange][iter]->SetName(Form("h_final_xj_pp_systematics_%d_%d", irange, iter));
	  g_final_xj_pp_systematics[irange][iter] = histo_opps::get_xj_systematics(h_final_xj_pp_systematics[irange][iter], h_total_pp_sys_neg_range[irange][iter], h_total_pp_sys_range[irange][iter], nbins);
	      
	}

    }

  TCanvas *cxj_money = new TCanvas("cxj_money","cxj_money", 500, 500);

  for (int irange = 0; irange < mbins; irange++)
    {

      h_linear_truth_xj[irange]->Scale(1./h_linear_truth_xj[irange]->Integral(0, -1, "width"));
      for (int ic = 0; ic <cent_bins; ic++)
	{
	  dlutility::SetLineAtt(h_final_xj_unfold_range[ic][irange][niter], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[ic][irange][niter], color_unfold[ic], msize_unfold, marker_unfold);
	  dlutility::SetLineAtt(g_final_xj_systematics[ic][irange][niter], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_xj_systematics[ic][irange][niter], color_unfold[ic], msize_unfold, marker_unfold);
	  g_final_xj_systematics[ic][irange][niter]->SetFillColorAlpha(color_unfold_fill[ic], 0.3); 
	  g_final_xj_systematics[ic][irange][niter]->SetFillStyle(0);//ColorAlpha(color_unfold_fill[ic], 0.3);


	}
      
      dlutility::SetLineAtt(h_final_xj_pp_unfold_range[irange][niter], color_pp_unfold, lsize_pp_unfold, 1);
      dlutility::SetMarkerAtt(h_final_xj_pp_unfold_range[irange][niter], color_pp_unfold, msize_pp_unfold, marker_pp_unfold);

      
      dlutility::SetLineAtt(g_final_xj_pp_systematics[irange][niter], color_pp_unfold, lsize_pp_unfold, 1);
      dlutility::SetMarkerAtt(g_final_xj_pp_systematics[irange][niter], color_pp_unfold, msize_pp_unfold, marker_pp_unfold);
      g_final_xj_pp_systematics[irange][niter]->SetFillColorAlpha(color_pp_unfold_fill, 0.3); 

      dlutility::SetLineAtt(h_linear_truth_xj[irange], color_pythia, 3, 1);
      dlutility::SetMarkerAtt(h_linear_truth_xj[irange], color_pythia, msize_pythia, marker_pythia);



      h_linear_truth_xj[irange]->SetMaximum(4.5);
      h_linear_truth_xj[irange]->SetMinimum(0);
      h_linear_truth_xj[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

            //h_final_xj_truth_range[irange]->Draw("E4 same");
      //TH1D *hh = (TH1D*) h_final_xj_herwig_range[irange]->Rebin(nbins - first_bin, Form("h_rebin_herwig_%d", irange), &dxj_bins[first_bin]);
      TH1D *hs[cent_bins];
      TH1D *hs_pp = (TH1D*) h_final_xj_pp_systematics[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_pp_sys_%d", irange), &dxj_bins[first_bin]);
      TH1D *hu[cent_bins];
      for (int ic = 0; ic < cent_bins; ic++)
	{
	  hs[ic] = (TH1D*) h_final_xj_systematics[ic][irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_sys_%d_%d", ic, irange), &dxj_bins[first_bin]);
	  hu[ic] = (TH1D*) h_final_xj_unfold_range[ic][irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d_%d", ic, irange), &dxj_bins[first_bin]);
	}
      TH1D *hpp = (TH1D*) h_final_xj_pp_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_pp_unf_%d", irange), &dxj_bins[first_bin]);

      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.17);
      hu[0]->GetYaxis()->SetTitleOffset(1.8);
      dlutility::SetFont(hu[0], 42, 0.06, 0.04, 0.05, 0.05);
      hu[0]->SetMaximum(5);
      hu[0]->SetMinimum(0);
      hu[0]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      //ht->Draw("E4 same");
      hu[0]->Draw("p E1");
      for (int ic = 0; ic < cent_bins; ic ++)
	{
	  g_final_xj_systematics[ic][irange][niter]->Draw("same p E2");
	}
      g_final_xj_pp_systematics[irange][niter]->Draw("same p E2");
      //hs->Draw("same p E2");
      //h_linear_truth_xj[irange]->Draw("same hist C");
      //      h_linear_herwig_xj[irange]->Draw("same hist C");
      hpp->Draw("same p E1");
      for (int i = 0; i < cent_bins; i++)
	{
	  hu[i]->Draw("same p E1");
	}


      float top = 0.88;
      float ss = 0.05;
      dlutility::DrawSPHENIX(0.22, top);

      dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.1]{#it{R = %0.1f}}", cone_size*0.1), 0.22, top - 2*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, top - 3*ss);
      dlutility::drawText(Form("p_{T,2}^{lead} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, top - 4*ss);
      dlutility::drawText(Form("#Delta#phi #geq %s", dphi_string.c_str()), 0.22, top - 5*ss);
      
      
      TLegend *leg = new TLegend(0.65, top - 3*ss, 0.8, top);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      
      leg->SetTextFont(42);
      leg->AddEntry(hu[0], "0-10%");
      leg->AddEntry(hu[1], "10-30%");
      leg->AddEntry(hu[2], "30-50%");
      leg->AddEntry(hu[3], "50-90%");
      leg->AddEntry(hpp, "pp");
      leg->Draw("same");

      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_AA_cent_all_r%02d_range_%d.png",  rb.get_code_location().c_str(), cone_size, irange));
      cxj_money->Print(Form("%s/final_plots/h_final_xj_unfolded_AA_cent_all_r%02d_range_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange));
    }

  return;
}
