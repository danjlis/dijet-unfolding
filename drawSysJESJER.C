#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

const int color_unfold = kAzure - 6;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

const int color_njet = kOrange + 7;
const float marker_njet = 20;
const float msize_njet = 0.9;
const float lsize_njet = 1.1;

const int color_half = kOrange + 7;
const float marker_half = 20;
const float msize_half = 0.9;
const float lsize_half = 1.1;

const int color_herwig = kBlue + 1;
const float marker_herwig = 20;
const float msize_herwig = 0.9;
const float lsize_herwig = 1.1;

const int color_prior = kRed - 2;
const float marker_prior = 20;
const float msize_prior = 0.9;
const float lsize_prior = 1.1;

const int color_njer = kPink + 2;
const float marker_njer = 20;
const float msize_njer = 0.9;
const float lsize_njer = 1.1;

const int color_pjer = kMagenta + 2;
const float marker_pjer = 20;
const float msize_pjer = 0.9;
const float lsize_pjer = 1.1;

const int color_pjes = kCyan + 1;
const float marker_pjes = 20;
const float msize_pjes = 0.9;
const float lsize_pjes = 1.1;

const int color_njes = kGreen  - 2;
const float marker_njes = 20;
const float msize_njes = 0.9;
const float lsize_njes = 1.1;

const int color_pythia = kBlack;
const float msize_pythia = 0.7;
const float marker_pythia = 21;
const float lsize_pythia = 1.1;
const int color_reco = kRed - 2;
const float marker_reco = 21;
const float msize_reco = 0.7;
const float lsize_reco = 1.1;
const int color_data = kBlue - 9;
const float marker_data = 20;
const float msize_data = 0.9;
const float lsize_data = 1.1;

const bool doJER = true;
const bool doJES = true;
const bool doPRIOR = true;
const bool doNJET = false;
void drawSysJESJER(const int cone_size = 4)
{

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  int first_bin = 0;
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];
  double  dxj_bins[nbins+1];
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
      std::cout << measure_bins[ir] << " -- " <<  subleading_measure_bins[ir] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  
  const int niterations = 10;
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root", rb.get_code_location().c_str(), cone_size),"r");

  TProfile *h_xj_rms[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_xj_rms[iter] = (TProfile*) finu->Get(Form("hp_xj_%d", iter));
    }
  TFile *fin = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_nominal.root", rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }
  // Closure
  TH1D *h_sys_half[mbins][niterations];
  TFile *finhalf = new TFile(Form("%s/uncertainties/half_closure_pp_r%02d_nominal.root", rb.get_code_location().c_str(), cone_size),"r");
  std::cout << "Half" << std::endl;
  TH1D *h_closure_test[mbins][10];
  std::cout << "Half" << std::endl;
  for (int irange = 0; irange < mbins; irange++)
    {
      std::cout << "Half" << std::endl;
      for (int iter = 0; iter < 10; iter++)
	{

	  h_closure_test[irange][iter] = (TH1D*) finhalf->Get(Form("h_closure_test_%d_%d", irange, iter));

	  if (!h_closure_test[irange][iter])
	    {
	      std::cout << "nothing in half " << irange << " / " << iter << std::endl;
	    }
	  h_sys_half[irange][iter] = (TH1D*) h_closure_test[irange][iter]->Clone();
	  h_sys_half[irange][iter]->SetName(Form("h_sys_half_%d_%d", irange, iter));
	  h_sys_half[irange][iter]->SetMinimum(-5);
	  h_sys_half[irange][iter]->SetMaximum(5);
	  for (int ib = 0; ib < h_sys_half[irange][iter]->GetNbinsX(); ib++)
	    {
	      h_sys_half[irange][iter]->SetBinContent(ib+1, h_sys_half[irange][iter]->GetBinContent(ib+1));
	    }
	  h_sys_half[irange][iter] = (TH1D*) h_sys_half[irange][iter]->Rebin(nbins - first_bin, Form("h_rebin_half_%d_%d", irange, iter), &dxj_bins[first_bin]);
	  
	}
    }
  std::cout << "Half" << std::endl;
	
  // HERWIG
  TFile *finherwig = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_HERWIG.root", rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_herwig_flat_data_pt1pt2 = (TH1D*) finherwig->Get("h_data_flat_pt1pt2");
  h_herwig_flat_data_pt1pt2->SetName("h_herwig_flat_data_pt1pt2");
  TH1D *h_herwig_flat_reco_pt1pt2 = (TH1D*) finherwig->Get("h_reco_flat_pt1pt2");
  h_herwig_flat_reco_pt1pt2->SetName("h_herwig_flat_reco_pt1pt2");
  TH1D *h_herwig_flat_truth_pt1pt2 = (TH1D*) finherwig->Get("h_truth_flat_pt1pt2");
  h_herwig_flat_truth_pt1pt2->SetName("h_herwig_flat_truth_pt1pt2");
  TH1D *h_herwig_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_herwig_flat_unfold_pt1pt2[iter] = (TH1D*) finherwig->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_herwig_flat_unfold_pt1pt2[iter]->SetName(Form("h_herwig_flat_unfold_pt1pt2_%d", iter));
    }
  // JES +
  TFile *finpjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_posJES.root", rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_pjes_flat_data_pt1pt2 = (TH1D*) finpjes->Get("h_data_flat_pt1pt2");
  h_pjes_flat_data_pt1pt2->SetName("h_pjes_flat_data_pt1pt2");
  TH1D *h_pjes_flat_reco_pt1pt2 = (TH1D*) finpjes->Get("h_reco_flat_pt1pt2");
  h_pjes_flat_reco_pt1pt2->SetName("h_pjes_flat_reco_pt1pt2");
  TH1D *h_pjes_flat_truth_pt1pt2 = (TH1D*) finpjes->Get("h_truth_flat_pt1pt2");
  h_pjes_flat_truth_pt1pt2->SetName("h_pjes_flat_truth_pt1pt2");
  TH1D *h_pjes_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_flat_unfold_pt1pt2[iter] = (TH1D*) finpjes->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_pjes_flat_unfold_pt1pt2[iter]->SetName(Form("h_pjes_flat_unfold_pt1pt2_%d", iter));
    }

  // JES -
  TFile *finnjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_negJES.root", rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_njes_flat_data_pt1pt2 = (TH1D*) finnjes->Get("h_data_flat_pt1pt2");
  h_njes_flat_data_pt1pt2->SetName("h_njes_flat_data_pt1pt2");
  TH1D *h_njes_flat_reco_pt1pt2 = (TH1D*) finnjes->Get("h_reco_flat_pt1pt2");
  h_njes_flat_reco_pt1pt2->SetName("h_njes_flat_reco_pt1pt2");
  TH1D *h_njes_flat_truth_pt1pt2 = (TH1D*) finnjes->Get("h_truth_flat_pt1pt2");
  h_njes_flat_truth_pt1pt2->SetName("h_njes_flat_truth_pt1pt2");
  TH1D *h_njes_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_flat_unfold_pt1pt2[iter] = (TH1D*) finnjes->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_njes_flat_unfold_pt1pt2[iter]->SetName(Form("h_njes_flat_unfold_pt1pt2_%d", iter));
    }
  // posJER
  TFile *finpjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_posJER.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_pjer_flat_data_pt1pt2 = (TH1D*) finpjer->Get("h_data_flat_pt1pt2");
  h_pjer_flat_data_pt1pt2->SetName("h_pjer_flat_data_pt1pt2");
  TH1D *h_pjer_flat_reco_pt1pt2 = (TH1D*) finpjer->Get("h_reco_flat_pt1pt2");
  h_pjer_flat_reco_pt1pt2->SetName("h_pjer_flat_reco_pt1pt2");
  TH1D *h_pjer_flat_truth_pt1pt2 = (TH1D*) finpjer->Get("h_truth_flat_pt1pt2");
  h_pjer_flat_truth_pt1pt2->SetName("h_pjer_flat_truth_pt1pt2");
  TH1D *h_pjer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjer_flat_unfold_pt1pt2[iter] = (TH1D*) finpjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_pjer_flat_unfold_pt1pt2[iter]->SetName(Form("h_pjer_flat_unfold_pt1pt2_%d", iter));
    }
  // njer
  TFile *finnjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_negJER.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_njer_flat_data_pt1pt2 = (TH1D*) finnjer->Get("h_data_flat_pt1pt2");
  h_njer_flat_data_pt1pt2->SetName("h_njer_flat_data_pt1pt2");
  TH1D *h_njer_flat_reco_pt1pt2 = (TH1D*) finnjer->Get("h_reco_flat_pt1pt2");
  h_njer_flat_reco_pt1pt2->SetName("h_njer_flat_reco_pt1pt2");
  TH1D *h_njer_flat_truth_pt1pt2 = (TH1D*) finnjer->Get("h_truth_flat_pt1pt2");
  h_njer_flat_truth_pt1pt2->SetName("h_njer_flat_truth_pt1pt2");
  TH1D *h_njer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njer_flat_unfold_pt1pt2[iter] = (TH1D*) finnjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_njer_flat_unfold_pt1pt2[iter]->SetName(Form("h_njer_flat_unfold_pt1pt2_%d", iter));
    }

  /* // VTX */
  /* TFile *finvtx = new TFile("%s/unfolding_hists_VTX.root","r"); */

  /* TH1D *h_vtx_flat_data_pt1pt2 = (TH1D*) finvtx->Get("h_data_flat_pt1pt2"); */
  /* h_vtx_flat_data_pt1pt2->SetName("h_vtx_flat_data_pt1pt2"); */
  /* TH1D *h_vtx_flat_reco_pt1pt2 = (TH1D*) finvtx->Get("h_reco_flat_pt1pt2"); */
  /* h_vtx_flat_reco_pt1pt2->SetName("h_vtx_flat_reco_pt1pt2"); */
  /* TH1D *h_vtx_flat_truth_pt1pt2 = (TH1D*) finvtx->Get("h_truth_flat_pt1pt2"); */
  /* h_vtx_flat_truth_pt1pt2->SetName("h_vtx_flat_truth_pt1pt2"); */
  /* TH1D *h_vtx_flat_unfold_pt1pt2[niterations]; */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_flat_unfold_pt1pt2[iter] = (TH1D*) finvtx->Get(Form("h_flat_unfold_pt1pt2_%d", iter)); */
  /*     h_vtx_flat_unfold_pt1pt2[iter]->SetName(Form("h_vtx_flat_unfold_pt1pt2_%d", iter)); */
  /*   } */

  // NJET 

  // TFile *finnjet = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_NJET.root",  rb.get_code_location().c_str(), cone_size),"r");

  // TH1D *h_njet_flat_data_pt1pt2 = (TH1D*) finnjet->Get("h_data_flat_pt1pt2");
  // h_njet_flat_data_pt1pt2->SetName("h_njet_flat_data_pt1pt2");
  // TH1D *h_njet_flat_reco_pt1pt2 = (TH1D*) finnjet->Get("h_reco_flat_pt1pt2");
  // h_njet_flat_reco_pt1pt2->SetName("h_njet_flat_reco_pt1pt2");
  // TH1D *h_njet_flat_truth_pt1pt2 = (TH1D*) finnjet->Get("h_truth_flat_pt1pt2");
  // h_njet_flat_truth_pt1pt2->SetName("h_njet_flat_truth_pt1pt2");
  // TH1D *h_njet_flat_unfold_pt1pt2[niterations];
  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     h_njet_flat_unfold_pt1pt2[iter] = (TH1D*) finnjet->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
  //     h_njet_flat_unfold_pt1pt2[iter]->SetName(Form("h_njet_flat_unfold_pt1pt2_%d", iter));
  //   }

  // PRIOR 
  TFile *finprior = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_PRIMER2_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_prior_flat_data_pt1pt2 = (TH1D*) finprior->Get("h_data_flat_pt1pt2");
  h_prior_flat_data_pt1pt2->SetName("h_prior_flat_data_pt1pt2");
  TH1D *h_prior_flat_reco_pt1pt2 = (TH1D*) finprior->Get("h_reco_flat_pt1pt2");
  h_prior_flat_reco_pt1pt2->SetName("h_prior_flat_reco_pt1pt2");
  TH1D *h_prior_flat_truth_pt1pt2 = (TH1D*) finprior->Get("h_truth_flat_pt1pt2");
  h_prior_flat_truth_pt1pt2->SetName("h_prior_flat_truth_pt1pt2");
  TH1D *h_prior_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_flat_unfold_pt1pt2[iter] = (TH1D*) finprior->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_prior_flat_unfold_pt1pt2[iter]->SetName(Form("h_prior_flat_unfold_pt1pt2_%d", iter));
    }


  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_xj_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_data_range[irange] = new TH1D(Form("h_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  // herwig
  TH2D *h_herwig_pt1pt2_reco = new TH2D("h_herwig_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_herwig_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_herwig_pt1pt2_unfold[iter] = new TH2D("h_herwig_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_herwig_pt1pt2_unfold[iter]->SetName(Form("h_herwig_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_herwig_xj_reco = new TH1D("h_herwig_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_herwig_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_herwig_xj_unfold[iter] = new TH1D(Form("h_herwig_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  
  TH1D *h_herwig_xj_reco_range[mbins];
  TH1D *h_herwig_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_herwig_xj_reco_range[irange] = new TH1D(Form("h_herwig_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_herwig_xj_unfold_range[irange][iter] = new TH1D(Form("h_herwig_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  // pJES
  TH2D *h_pjes_pt1pt2_reco = new TH2D("h_pjes_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pjes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_pt1pt2_unfold[iter] = new TH2D("h_pjes_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pjes_pt1pt2_unfold[iter]->SetName(Form("h_pjes_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_pjes_xj_reco = new TH1D("h_pjes_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_pjes_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_xj_unfold[iter] = new TH1D(Form("h_pjes_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  
  TH1D *h_pjes_xj_reco_range[mbins];
  TH1D *h_pjes_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjes_xj_reco_range[irange] = new TH1D(Form("h_pjes_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjes_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjes_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  // nJES

  TH2D *h_njes_pt1pt2_reco = new TH2D("h_njes_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_njes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_pt1pt2_unfold[iter] = new TH2D("h_njes_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_njes_pt1pt2_unfold[iter]->SetName(Form("h_njes_pt1pt2_unfold_iter%d", iter));
    }
  TH1D *h_njes_xj_reco = new TH1D("h_njes_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_njes_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_xj_unfold[iter] = new TH1D(Form("h_njes_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  
  TH1D *h_njes_xj_reco_range[mbins];
  TH1D *h_njes_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njes_xj_reco_range[irange] = new TH1D(Form("h_njes_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njes_xj_unfold_range[irange][iter] = new TH1D(Form("h_njes_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  //POSJER
  TH2D *h_pjer_pt1pt2_reco = new TH2D("h_pjer_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pjer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjer_pt1pt2_unfold[iter] = new TH2D("h_pjer_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pjer_pt1pt2_unfold[iter]->SetName(Form("h_pjer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_pjer_xj_reco = new TH1D("h_pjer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_pjer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjer_xj_unfold[iter] = new TH1D(Form("h_pjer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_pjer_xj_reco_range[mbins];
  TH1D *h_pjer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjer_xj_reco_range[irange] = new TH1D(Form("h_pjer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjer_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  //NJER
  TH2D *h_njer_pt1pt2_reco = new TH2D("h_njer_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_njer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_njer_pt1pt2_unfold[iter] = new TH2D("h_njer_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_njer_pt1pt2_unfold[iter]->SetName(Form("h_njer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_njer_xj_reco = new TH1D("h_njer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_njer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njer_xj_unfold[iter] = new TH1D(Form("h_njer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_njer_xj_reco_range[mbins];
  TH1D *h_njer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njer_xj_reco_range[irange] = new TH1D(Form("h_njer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njer_xj_unfold_range[irange][iter] = new TH1D(Form("h_njer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  /* //VTX */
  /* TH2D *h_vtx_pt1pt2_reco = new TH2D("h_vtx_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins); */
  /* TH2D *h_vtx_pt1pt2_unfold[niterations]; */

  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_pt1pt2_unfold[iter] = new TH2D("h_vtx_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins); */
  /*     h_vtx_pt1pt2_unfold[iter]->SetName(Form("h_vtx_pt1pt2_unfold_iter%d", iter)); */
  /*   } */

  /* TH1D *h_vtx_xj_reco = new TH1D("h_vtx_xj_reco", ";x_{J};", nbins, ixj_bins); */
  /* TH1D *h_vtx_xj_unfold[niterations]; */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_xj_unfold[iter] = new TH1D(Form("h_vtx_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins); */
  /*   } */

  /* TH1D *h_vtx_xj_reco_range[mbins]; */
  /* TH1D *h_vtx_xj_unfold_range[mbins][niterations]; */
  /* for (int irange = 0; irange < mbins; irange++) */
  /*   { */
  /*     h_vtx_xj_reco_range[irange] = new TH1D(Form("h_vtx_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins); */
  /*     for (int iter = 0; iter < niterations; iter++) */
  /* 	{ */
  /* 	  h_vtx_xj_unfold_range[irange][iter] = new TH1D(Form("h_vtx_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins); */
  /* 	} */
  /*   } */
  //NJET
  // TH2D *h_njet_pt1pt2_reco = new TH2D("h_njet_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  // TH2D *h_njet_pt1pt2_unfold[niterations];

  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     h_njet_pt1pt2_unfold[iter] = new TH2D("h_njet_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
  //     h_njet_pt1pt2_unfold[iter]->SetName(Form("h_njet_pt1pt2_unfold_iter%d", iter));
  //   }

  // TH1D *h_njet_xj_reco = new TH1D("h_njet_xj_reco", ";x_{J};", nbins, ixj_bins);
  // TH1D *h_njet_xj_unfold[niterations];
  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     h_njet_xj_unfold[iter] = new TH1D(Form("h_njet_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
  //   }

  // TH1D *h_njet_xj_reco_range[mbins];
  // TH1D *h_njet_xj_unfold_range[mbins][niterations];
  // for (int irange = 0; irange < mbins; irange++)
  //   {
  //     h_njet_xj_reco_range[irange] = new TH1D(Form("h_njet_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
  //     for (int iter = 0; iter < niterations; iter++)
  // 	{
  // 	  h_njet_xj_unfold_range[irange][iter] = new TH1D(Form("h_njet_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
  // 	}
  //   }
  //PRIOR
  TH2D *h_prior_pt1pt2_reco = new TH2D("h_prior_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_prior_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_pt1pt2_unfold[iter] = new TH2D("h_prior_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_prior_pt1pt2_unfold[iter]->SetName(Form("h_prior_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_prior_xj_reco = new TH1D("h_prior_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_prior_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_xj_unfold[iter] = new TH1D(Form("h_prior_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_prior_xj_reco_range[mbins];
  TH1D *h_prior_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_prior_xj_reco_range[irange] = new TH1D(Form("h_prior_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_prior_xj_unfold_range[irange][iter] = new TH1D(Form("h_prior_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }

  //herwig
  histo_opps::make_sym_pt1pt2(h_herwig_flat_reco_pt1pt2, h_herwig_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_herwig_flat_unfold_pt1pt2[iter], h_herwig_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_herwig_pt1pt2_reco, h_herwig_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_herwig_pt1pt2_unfold[iter], h_herwig_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_herwig_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_herwig_xj_unfold[iter], nbins);
    }
  //pjes
  histo_opps::make_sym_pt1pt2(h_pjes_flat_reco_pt1pt2, h_pjes_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_pjes_flat_unfold_pt1pt2[iter], h_pjes_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pjes_pt1pt2_reco, h_pjes_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pjes_pt1pt2_unfold[iter], h_pjes_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_pjes_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_pjes_xj_unfold[iter], nbins);
    }

  //njes
  histo_opps::make_sym_pt1pt2(h_njes_flat_reco_pt1pt2, h_njes_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_njes_flat_unfold_pt1pt2[iter], h_njes_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_njes_pt1pt2_reco, h_njes_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_njes_pt1pt2_unfold[iter], h_njes_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_njes_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_njes_xj_unfold[iter], nbins);
    }
  //pjer
  histo_opps::make_sym_pt1pt2(h_pjer_flat_reco_pt1pt2, h_pjer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_pjer_flat_unfold_pt1pt2[iter], h_pjer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pjer_pt1pt2_reco, h_pjer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pjer_pt1pt2_unfold[iter], h_pjer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_pjer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_pjer_xj_unfold[iter], nbins);
    }
  //njer
  histo_opps::make_sym_pt1pt2(h_njer_flat_reco_pt1pt2, h_njer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_njer_flat_unfold_pt1pt2[iter], h_njer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_njer_pt1pt2_reco, h_njer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_njer_pt1pt2_unfold[iter], h_njer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_njer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_njer_xj_unfold[iter], nbins);
    }
  /* //vtx */
  /* histo_opps::make_sym_pt1pt2(h_vtx_flat_reco_pt1pt2, h_vtx_pt1pt2_reco, nbins); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::make_sym_pt1pt2(h_vtx_flat_unfold_pt1pt2[iter], h_vtx_pt1pt2_unfold[iter], nbins); */
  /*   } */

  /* histo_opps::project_xj(h_vtx_pt1pt2_reco, h_vtx_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::project_xj(h_vtx_pt1pt2_unfold[iter], h_vtx_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2); */
  /*   } */
  
  /* histo_opps::normalize_histo(h_vtx_xj_reco, nbins); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::normalize_histo(h_vtx_xj_unfold[iter], nbins); */
  /*   } */

  //njet
  // histo_opps::make_sym_pt1pt2(h_njet_flat_reco_pt1pt2, h_njet_pt1pt2_reco, nbins);
  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     histo_opps::make_sym_pt1pt2(h_njet_flat_unfold_pt1pt2[iter], h_njet_pt1pt2_unfold[iter], nbins);
  //   }

  // histo_opps::project_xj(h_njet_pt1pt2_reco, h_njet_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     histo_opps::project_xj(h_njet_pt1pt2_unfold[iter], h_njet_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  //   }
  
  // histo_opps::normalize_histo(h_njet_xj_reco, nbins);
  // for (int iter = 0; iter < niterations; iter++)
  //   {
  //     histo_opps::normalize_histo(h_njet_xj_unfold[iter], nbins);
  //   }

  //prior
  histo_opps::make_sym_pt1pt2(h_prior_flat_reco_pt1pt2, h_prior_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_prior_flat_unfold_pt1pt2[iter], h_prior_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_prior_pt1pt2_reco, h_prior_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_prior_pt1pt2_unfold[iter], h_prior_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_prior_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_prior_xj_unfold[iter], nbins);
    }
  
  TH1D *h_final_xj_data_range[mbins];
  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_data_range[irange] = new TH1D(Form("h_final_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //herwig
  TH1D *h_herwig_final_xj_reco_range[mbins];
  TH1D *h_herwig_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_herwig_final_xj_reco_range[irange] = new TH1D(Form("h_herwig_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_herwig_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_herwig_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //pjes
  TH1D *h_pjes_final_xj_reco_range[mbins];
  TH1D *h_pjes_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjes_final_xj_reco_range[irange] = new TH1D(Form("h_pjes_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjes_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjes_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //njes
  TH1D *h_njes_final_xj_reco_range[mbins];
  TH1D *h_njes_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njes_final_xj_reco_range[irange] = new TH1D(Form("h_njes_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njes_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_njes_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //pjer
  TH1D *h_pjer_final_xj_reco_range[mbins];
  TH1D *h_pjer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjer_final_xj_reco_range[irange] = new TH1D(Form("h_pjer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //njer
  TH1D *h_njer_final_xj_reco_range[mbins];
  TH1D *h_njer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njer_final_xj_reco_range[irange] = new TH1D(Form("h_njer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_njer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  /* //vtx */
  /* TH1D *h_vtx_final_xj_reco_range[mbins]; */
  /* TH1D *h_vtx_final_xj_unfold_range[mbins][niterations]; */
  /* for (int irange = 0; irange < mbins; irange++) */
  /*   { */
  /*     h_vtx_final_xj_reco_range[irange] = new TH1D(Form("h_vtx_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins); */
  /*     for (int iter = 0; iter < niterations; iter++) */
  /* 	{ */
  /* 	  h_vtx_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_vtx_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins); */
  /* 	} */
  /*   } */

  //njet
  // TH1D *h_njet_final_xj_reco_range[mbins];
  // TH1D *h_njet_final_xj_unfold_range[mbins][niterations];
  // for (int irange = 0; irange < mbins; irange++)
  //   {
  //     h_njet_final_xj_reco_range[irange] = new TH1D(Form("h_njet_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
  //     for (int iter = 0; iter < niterations; iter++)
  // 	{
  // 	  h_njet_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_njet_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
  // 	}
  //   }

  //prior
  TH1D *h_prior_final_xj_reco_range[mbins];
  TH1D *h_prior_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_prior_final_xj_reco_range[irange] = new TH1D(Form("h_prior_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_prior_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_prior_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }


  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_data, h_xj_data_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_truth_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_data_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_data_range[irange], h_final_xj_data_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

    }

  for (int irange = 0; irange < mbins; irange++)
    {
      //herwig
      histo_opps::project_xj(h_herwig_pt1pt2_reco, h_herwig_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_herwig_pt1pt2_unfold[iter], h_herwig_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}

      histo_opps::normalize_histo(h_herwig_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_herwig_xj_unfold_range[irange][iter], nbins);
	}

      
      histo_opps::finalize_xj(h_herwig_xj_reco_range[irange], h_herwig_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_herwig_xj_unfold_range[irange][iter], h_herwig_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

      //pjes
      histo_opps::project_xj(h_pjes_pt1pt2_reco, h_pjes_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pjes_pt1pt2_unfold[iter], h_pjes_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}

      histo_opps::normalize_histo(h_pjes_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_pjes_xj_unfold_range[irange][iter], nbins);
	}

      
      histo_opps::finalize_xj(h_pjes_xj_reco_range[irange], h_pjes_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_pjes_xj_unfold_range[irange][iter], h_pjes_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}


      //njes
      histo_opps::project_xj(h_njes_pt1pt2_reco, h_njes_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_njes_pt1pt2_unfold[iter], h_njes_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_njes_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_njes_xj_unfold_range[irange][iter], nbins);
	}


      histo_opps::finalize_xj(h_njes_xj_reco_range[irange], h_njes_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_njes_xj_unfold_range[irange][iter], h_njes_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

      
      //pjer
      histo_opps::project_xj(h_pjer_pt1pt2_reco, h_pjer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pjer_pt1pt2_unfold[iter], h_pjer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}

      histo_opps::normalize_histo(h_pjer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_pjer_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_pjer_xj_reco_range[irange], h_pjer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_pjer_xj_unfold_range[irange][iter], h_pjer_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

      

      //njer
      histo_opps::project_xj(h_njer_pt1pt2_reco, h_njer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_njer_pt1pt2_unfold[iter], h_njer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_njer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_njer_xj_unfold_range[irange][iter], nbins);
	}

      
      histo_opps::finalize_xj(h_njer_xj_reco_range[irange], h_njer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_njer_xj_unfold_range[irange][iter], h_njer_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}


      /* //vtx */
      /* histo_opps::project_xj(h_vtx_pt1pt2_reco, h_vtx_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::project_xj(h_vtx_pt1pt2_unfold[iter], h_vtx_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2); */
      /* 	} */


      /* histo_opps::normalize_histo(h_vtx_xj_reco_range[irange], nbins); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::normalize_histo(h_vtx_xj_unfold_range[irange][iter], nbins); */
      /* 	} */

      /* histo_opps::finalize_xj(h_vtx_xj_reco_range[irange], h_vtx_final_xj_reco_range[irange], nbins, first_xj); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::finalize_xj(h_vtx_xj_unfold_range[irange][iter], h_vtx_final_xj_unfold_range[irange][iter], nbins, first_xj); */
      /* 	} */
      //njet

      // histo_opps::project_xj(h_njet_pt1pt2_reco, h_njet_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      // for (int iter = 0; iter < niterations; iter++)
      // 	{
      // 	  histo_opps::project_xj(h_njet_pt1pt2_unfold[iter], h_njet_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      // 	}


      // histo_opps::normalize_histo(h_njet_xj_reco_range[irange], nbins);
      // for (int iter = 0; iter < niterations; iter++)
      // 	{
      // 	  histo_opps::normalize_histo(h_njet_xj_unfold_range[irange][iter], nbins);
      // 	}

      // histo_opps::finalize_xj(h_njet_xj_reco_range[irange], h_njet_final_xj_reco_range[irange], nbins, first_xj);
      // for (int iter = 0; iter < niterations; iter++)
      // 	{
      // 	  histo_opps::finalize_xj(h_njet_xj_unfold_range[irange][iter], h_njet_final_xj_unfold_range[irange][iter], nbins, first_xj);
      // 	}

      //prior
      histo_opps::project_xj(h_prior_pt1pt2_reco, h_prior_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_prior_pt1pt2_unfold[iter], h_prior_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}

      histo_opps::normalize_histo(h_prior_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_prior_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_prior_xj_reco_range[irange], h_prior_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_prior_xj_unfold_range[irange][iter], h_prior_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}



    }
    

  // JES/JER systematics
  TH1D *h_sys_pjer[mbins][niterations];
  TH1D *h_sys_njer[mbins][niterations];
  TH1D *h_sys_pjes[mbins][niterations];  
  TH1D *h_sys_njes[mbins][niterations];

  // NJET systematics
  //TH1D *h_sys_njet[mbins][niterations];
  // PRIOR systematics
  TH1D *h_sys_prior[mbins][niterations];
  TH1D *h_sys_herwig[mbins][niterations];

  // VTX reweight systeamtic
  //  TH1D *h_sys_vtx[mbins][niterations];
  // total systematics
  TH1D *h_total_sys_range[mbins][niterations];
  TH1D *h_total_sys_neg_range[mbins][niterations];

  TCanvas *cxj_split = new TCanvas("cxj_split","cxj_split", 1800, 600);

  for (int niter = 0; niter < niterations; niter++)
    {
      cxj_split->Clear();
      
      dlutility::systematic_split_canvas(cxj_split, 3);
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxj_split->cd(1 + irange*2);

	  dlutility::SetLineAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, msize_pjes, marker_pjes);
	  dlutility::SetLineAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, lsize_njes, 1);
	  dlutility::SetMarkerAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, msize_njes, marker_njes);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetFont(h_final_xj_unfold_range[irange][niter], 42, 0.08);

	  h_final_xj_unfold_range[irange][niter]->SetMaximum(3.5);
	  h_final_xj_unfold_range[irange][niter]->SetMinimum(0);
	  h_final_xj_unfold_range[irange][niter]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hpjes = (TH1D*) h_pjes_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_pjes_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hnjes = (TH1D*) h_njes_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_njes_%d", irange), &dxj_bins[first_bin]);
	  
	  hu->Draw("p");
	  hpjes->Draw("same p");
	  hnjes->Draw("same p");
	  if (irange == 0)
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.4, 0.8 ,0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.4, 0.7, 0, kBlack, 0.08);
	    }
	  else
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.1, 0.8, 0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.1, 0.7, 0, kBlack, 0.08);
	    }

    
	  cxj_split->cd(2 + irange*2);

	  TH1D *h_pjes_compare = (TH1D*) hpjes->Clone();
	  h_pjes_compare->Divide(hu);
	  h_pjes_compare->SetTitle(";x_{J}; #frac{Var. - Nom.}{Nom}");
	  dlutility::SetFont(h_pjes_compare, 42, 0.12);
	  dlutility::SetLineAtt(h_pjes_compare, color_pjes, 1,1);
	  dlutility::SetMarkerAtt(h_pjes_compare, color_pjes, 1,8);

	  h_pjes_compare->SetMaximum(0.2);
	  h_pjes_compare->SetMinimum(-0.2);

	  TH1D *h_njes_compare = (TH1D*) hnjes->Clone();
	  h_njes_compare->Divide(hu);
	  dlutility::SetLineAtt(h_njes_compare, color_njes, 1,1);
	  dlutility::SetMarkerAtt(h_njes_compare, color_njes, 1,8);


	  for (int ib = 0; ib < h_njes_compare->GetNbinsX(); ib++)
	    {

	      if (h_njes_compare->GetBinCenter(ib+1) < 0.3) continue;
	      h_njes_compare->SetBinContent(ib+1, h_njes_compare->GetBinContent(ib+1) - 1);
	      h_pjes_compare->SetBinContent(ib+1, h_pjes_compare->GetBinContent(ib+1) - 1);
	    }

	  h_sys_njes[irange][niter] = (TH1D*) h_njes_compare->Clone();
	  h_sys_pjes[irange][niter] = (TH1D*) h_pjes_compare->Clone();

	  h_pjes_compare->SetFillColor(color_pjes);
	  h_njes_compare->SetFillColor(color_njes);
	  
	  h_pjes_compare->SetLineColor(kBlack);
	  h_njes_compare->SetLineColor(kBlack);

	  h_pjes_compare->Draw("hist");
	  h_njes_compare->Draw("hist same");

	  
	  TLine *line = new TLine(h_njes_compare->GetBinLowEdge(1), 0, 1, 0);
	  line->SetLineStyle(1);
	  line->SetLineColor(kBlack);
	  line->SetLineWidth(1);
	  line->Draw("same");

	}
      cxj_split->cd(7);
      dlutility::DrawSPHENIXpp(0.22, 0.84, 0.08);
      dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74, 0, kBlack, 0.08);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59, 0, kBlack, 0.08);
      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49, 0, kBlack, 0.08);
      TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.07);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_unfold_range[0][niter], "Nominal");
      leg->AddEntry(h_pjes_final_xj_unfold_range[0][niter], "+5% JES");
      leg->AddEntry(h_njes_final_xj_unfold_range[0][niter], "-5% JES");
      leg->Draw("same");

      cxj_split->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, niter));
      cxj_split->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, niter));
    }


  for (int niter = 0; niter < niterations; niter++)
    {
      cxj_split->Clear();
      
      dlutility::systematic_split_canvas(cxj_split, 3);
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxj_split->cd(1 + irange*2);

	  dlutility::SetLineAtt(h_pjer_final_xj_unfold_range[irange][niter], color_pjer, lsize_pjer, 1);
	  dlutility::SetMarkerAtt(h_pjer_final_xj_unfold_range[irange][niter], color_pjer, msize_pjer, marker_pjer);
	  dlutility::SetLineAtt(h_njer_final_xj_unfold_range[irange][niter], color_njer, lsize_njer, 1);
	  dlutility::SetMarkerAtt(h_njer_final_xj_unfold_range[irange][niter], color_njer, msize_njer, marker_njer);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetFont(h_final_xj_unfold_range[irange][niter], 42, 0.08);

	  h_final_xj_unfold_range[irange][niter]->SetMaximum(3.5);
	  h_final_xj_unfold_range[irange][niter]->SetMinimum(0);
	  h_final_xj_unfold_range[irange][niter]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hpjer = (TH1D*) h_pjer_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_pjer_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hnjer = (TH1D*) h_njer_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_njer_%d", irange), &dxj_bins[first_bin]);
	  
	  hu->Draw("p");
	  hpjer->Draw("same p");
	  hnjer->Draw("same p");
	  if (irange == 0)
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.4, 0.8 ,0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.4, 0.7, 0, kBlack, 0.08);
	    }
	  else
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.1, 0.8, 0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.1, 0.7, 0, kBlack, 0.08);
	    }

    
	  cxj_split->cd(2 + irange*2);

	  TH1D *h_pjer_compare = (TH1D*) hpjer->Clone();
	  h_pjer_compare->Divide(hu);
	  h_pjer_compare->SetTitle(";x_{J}; #frac{Var. - Nom.}{Nom}");
	  dlutility::SetFont(h_pjer_compare, 42, 0.12);
	  dlutility::SetLineAtt(h_pjer_compare, color_pjer, 1,1);
	  dlutility::SetMarkerAtt(h_pjer_compare, color_pjer, 1,8);

	  h_pjer_compare->SetMaximum(0.2);
	  h_pjer_compare->SetMinimum(-0.2);

	  TH1D *h_njer_compare = (TH1D*) hnjer->Clone();
	  h_njer_compare->Divide(hu);
	  dlutility::SetLineAtt(h_njer_compare, color_njer, 1,1);
	  dlutility::SetMarkerAtt(h_njer_compare, color_njer, 1,8);


	  for (int ib = 0; ib < h_njer_compare->GetNbinsX(); ib++)
	    {

	      if (h_njer_compare->GetBinCenter(ib+1) < 0.3) continue;
	      h_njer_compare->SetBinContent(ib+1, h_njer_compare->GetBinContent(ib+1) - 1);
	      h_pjer_compare->SetBinContent(ib+1, h_pjer_compare->GetBinContent(ib+1) - 1);
	    }

	  h_sys_njer[irange][niter] = (TH1D*) h_njer_compare->Clone();
	  h_sys_pjer[irange][niter] = (TH1D*) h_pjer_compare->Clone();

	  
	  h_pjer_compare->SetFillColor(color_pjer);
	  h_njer_compare->SetFillColor(color_njer);
	  
	  h_pjer_compare->SetLineColor(kBlack);
	  h_njer_compare->SetLineColor(kBlack);

	  h_pjer_compare->Draw("hist");
	  h_njer_compare->Draw("hist same");

	  
	  TLine *line = new TLine(h_njer_compare->GetBinLowEdge(1), 0, 1, 0);
	  line->SetLineStyle(1);
	  line->SetLineColor(kBlack);
	  line->SetLineWidth(1);
	  line->Draw("same");

	}
      cxj_split->cd(7);
      dlutility::DrawSPHENIXpp(0.22, 0.84, 0.08);
      dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74, 0, kBlack, 0.08);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59, 0, kBlack, 0.08);
      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49, 0, kBlack, 0.08);
      TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.07);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_unfold_range[0][niter], "Nominal");
      leg->AddEntry(h_pjer_final_xj_unfold_range[0][niter], "+5% JER");
      leg->AddEntry(h_njer_final_xj_unfold_range[0][niter], "-5% JER");
      leg->Draw("same");

      cxj_split->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, niter));
      cxj_split->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, niter));
    }

  int niter_prior = 2;
  TH1D *hprior[3];
  
  for (int ir = 0; ir < 3; ir++)
    {
      hprior[ir] = (TH1D*) h_prior_final_xj_unfold_range[ir][niter_prior]->Rebin(nbins - first_bin, Form("h_rebin_prior_%d", ir), &dxj_bins[first_bin]);
    }
  
  for (int niter = 0; niter < niterations; niter++)
    {
      cxj_split->Clear();
      
      dlutility::systematic_split_canvas(cxj_split, 3);
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxj_split->cd(1 + irange*2);

	  dlutility::SetLineAtt(hprior[irange], color_prior, lsize_prior, 1);
	  dlutility::SetMarkerAtt(hprior[niter], color_prior, msize_prior, marker_prior);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetFont(h_final_xj_unfold_range[irange][niter], 42, 0.08);

	  h_final_xj_unfold_range[irange][niter]->SetMaximum(3.5);
	  h_final_xj_unfold_range[irange][niter]->SetMinimum(0);
	  h_final_xj_unfold_range[irange][niter]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);
	  
	  hu->Draw("p");
	  hprior[irange]->Draw("same p");
	  if (irange == 0)
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.4, 0.8 ,0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.4, 0.7, 0, kBlack, 0.08);
	    }
	  else
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.1, 0.8, 0, kBlack, 0.08);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.1, 0.7, 0, kBlack, 0.08);
	    }

    
	  cxj_split->cd(2 + irange*2);

	  TH1D *h_prior_compare = (TH1D*) hprior[irange]->Clone();
	  h_prior_compare->Divide(hu);
	  h_prior_compare->SetTitle(";x_{J}; #frac{Var. - Nom.}{Nom}");
	  dlutility::SetFont(h_prior_compare, 42, 0.12);
	  dlutility::SetLineAtt(h_prior_compare, color_prior, 1,1);
	  dlutility::SetMarkerAtt(h_prior_compare, color_prior, 1,8);

	  h_prior_compare->SetMaximum(0.2);
	  h_prior_compare->SetMinimum(-0.2);


	  for (int ib = 0; ib < h_prior_compare->GetNbinsX(); ib++)
	    {

	      if (h_prior_compare->GetBinCenter(ib+1) < 0.3) continue;
	      h_prior_compare->SetBinContent(ib+1, h_prior_compare->GetBinContent(ib+1) - 1);
	    }

	  h_sys_prior[irange][niter] = (TH1D*) h_prior_compare->Clone();
		  
	  h_prior_compare->SetFillColor(color_prior);
	  
	  h_prior_compare->SetLineColor(kBlack);

	  h_prior_compare->Draw("hist");
	  
	  TLine *line = new TLine(h_prior_compare->GetBinLowEdge(1), 0, 1, 0);
	  line->SetLineStyle(1);
	  line->SetLineColor(kBlack);
	  line->SetLineWidth(1);
	  line->Draw("same");

	}
      cxj_split->cd(7);
      dlutility::DrawSPHENIXpp(0.22, 0.84, 0.08);
      dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74, 0, kBlack, 0.08);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59, 0, kBlack, 0.08);
      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49, 0, kBlack, 0.08);
      TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.07);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_unfold_range[0][niter], "Nominal");
      leg->AddEntry(h_prior_final_xj_unfold_range[0][niter], " No Prior Weighting");
      leg->Draw("same");

      cxj_split->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, niter));
      cxj_split->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, niter));
    }

  for (int niter = 0; niter < niterations; niter++)
    {
      cxj_split->Clear();
      
      dlutility::systematic_split_canvas(cxj_split, 3);
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxj_split->cd(1 + irange*2);

	  dlutility::SetLineAtt(h_herwig_final_xj_unfold_range[irange][niter], color_herwig, lsize_herwig, 1);
	  dlutility::SetMarkerAtt(h_herwig_final_xj_unfold_range[irange][niter], color_herwig, msize_herwig, marker_herwig);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetFont(h_final_xj_unfold_range[irange][niter], 42, 0.08);

	  h_final_xj_unfold_range[irange][niter]->SetMaximum(3.5);
	  h_final_xj_unfold_range[irange][niter]->SetMinimum(0);
	  h_final_xj_unfold_range[irange][niter]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hherwig = (TH1D*) h_herwig_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_herwig_%d", irange), &dxj_bins[first_bin]);
	  
	  hu->Draw("p");
	  hherwig->Draw("same p");

    
	  cxj_split->cd(2 + irange*2);

	  TH1D *h_herwig_compare = (TH1D*) hherwig->Clone();
	  h_herwig_compare->Divide(hu);
	  h_herwig_compare->SetTitle(";x_{J}; #frac{Var. - Nom.}{Nom}");
	  dlutility::SetFont(h_herwig_compare, 42, 0.12);
	  dlutility::SetLineAtt(h_herwig_compare, color_herwig, 1,1);
	  dlutility::SetMarkerAtt(h_herwig_compare, color_herwig, 1,8);

	  h_herwig_compare->SetMaximum(0.2);
	  h_herwig_compare->SetMinimum(-0.2);


	  for (int ib = 0; ib < h_herwig_compare->GetNbinsX(); ib++)
	    {

	      if (h_herwig_compare->GetBinCenter(ib+1) < 0.3) continue;
	      h_herwig_compare->SetBinContent(ib+1, h_herwig_compare->GetBinContent(ib+1) - 1);
	    }

	  h_sys_herwig[irange][niter] = (TH1D*) h_herwig_compare->Clone();
	  h_herwig_compare->SetFillColor(color_herwig);
	  
	  h_herwig_compare->SetLineColor(kBlack);

	  h_herwig_compare->Draw("hist");
	  
	  TLine *line = new TLine(h_herwig_compare->GetBinLowEdge(1), 0, 1, 0);
	  line->SetLineStyle(1);
	  line->SetLineColor(kBlack);
	  line->SetLineWidth(1);
	  line->Draw("same");

	}
      cxj_split->cd(7);
      dlutility::DrawSPHENIXpp(0.22, 0.84, 0.08);
      dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74, 0, kBlack, 0.08);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59, 0, kBlack, 0.08);
      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49, 0, kBlack, 0.08);
      TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.07);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_unfold_range[0][niter], "Nominal");
      leg->AddEntry(h_herwig_final_xj_unfold_range[0][niter], " Herwig Unfold");
      leg->Draw("same");

      cxj_split->SaveAs(Form("%s/systematic_plots/h_HERWIG_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, niter));
      cxj_split->SaveAs(Form("%s/systematic_plots/h_HERWIG_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, niter));
    }

  TCanvas *cxjsys = new TCanvas("cxjsys","cxjsys", 1500, 400);
  for (int niter = 0; niter < 10; niter++)
    {
      dlutility::systematic_split_canvas(cxjsys, 3, 0);
      TLegend *leg4 = new TLegend(0.02, 0.30, 0.9, 0.58);
      leg4->SetTextSize(0.08);
      leg4->SetTextFont(42);
      leg4->SetLineWidth(0);

      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxjsys->cd(irange+1);
	  
	  TH1D *h_sys_jer = (TH1D*) h_sys_pjer[irange][niter]->Clone();
	  TH1D *h_sys_jer_flip = (TH1D*) h_sys_pjer[irange][niter]->Clone();

	  TH1D *h_sys_half_flip = (TH1D*) h_sys_half[irange][niter]->Clone();
	  TH1D *h_sys_prior_flip = (TH1D*) h_sys_prior[irange][niter]->Clone();
	  TH1D *h_sys_herwig_flip = (TH1D*) h_sys_herwig[irange][niter]->Clone();

	  TH1D *h_sys_jes = (TH1D*) h_sys_njes[irange][niter]->Clone();
	  TH1D *h_sys_jes_flip = (TH1D*) h_sys_njes[irange][niter]->Clone();	  

	  TH1D *h_total_sys = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys->Reset();
	  h_total_sys->SetName("h_total_sys");
	  h_total_sys->SetTitle(";x_{J}; #frac{Var - Nom}{Nom}");

	  TH1D *h_total_sys_flip = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_flip->Reset();
	  h_total_sys_flip->SetName("h_total_sys_flip");
	  h_total_sys_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_range[irange][niter]->Reset();
	  h_total_sys_range[irange][niter]->SetName(Form("h_total_sys_range_%d_iter_%d", irange, niter));
	  h_total_sys_neg_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_neg_range[irange][niter]->Reset();
	  h_total_sys_neg_range[irange][niter]->SetName(Form("h_total_sys_neg_range_%d_iter_%d", irange, niter));

      
	  for (int i = 0; i < h_sys_prior_flip->GetNbinsX(); i++)
	    {
	      
	      //h_sys_njet_flip->SetBinContent(i+1, 2 - h_sys_njet[irange][niter]->GetBinContent(i+1));
	      h_sys_prior_flip->SetBinContent(i+1, -1*h_sys_prior[irange][niter]->GetBinContent(i+1));
	      h_sys_herwig_flip->SetBinContent(i+1, -1*h_sys_herwig[irange][niter]->GetBinContent(i+1));
	      h_sys_half_flip->SetBinContent(i+1, -1*h_sys_half[irange][niter]->GetBinContent(i+1));

	      float maxjer = h_sys_pjer[irange][niter]->GetBinContent(i+1);
	      float minjer = h_sys_njer[irange][niter]->GetBinContent(i+1);;
	      if (maxjer < 0 && minjer > 0)
		{
		  maxjer = h_sys_njer[irange][niter]->GetBinContent(i+1);
		  minjer = h_sys_pjer[irange][niter]->GetBinContent(i+1);
		}
	      else if(maxjer > 0 && minjer > 0)
		{
		  maxjer = std::max(maxjer, minjer);
		  minjer = 0;
		}
	      else if(maxjer < 0 && minjer < 0)
		{
		  minjer = std::min(maxjer, minjer);
		  maxjer = 0;
		}
	      float maxjes = h_sys_pjes[irange][niter]->GetBinContent(i+1);
	      float minjes = h_sys_njes[irange][niter]->GetBinContent(i+1);
	      if (maxjes < 0 && minjes > 0)
		{
		  maxjes = h_sys_njes[irange][niter]->GetBinContent(i+1);
		  minjes = h_sys_pjes[irange][niter]->GetBinContent(i+1);
		}
	      else if(maxjes > 0 && minjes > 0)
		{
		  maxjes = std::max(maxjes, minjes);
		  minjes = 0;
		}
	      else if(maxjes < 0 && minjes < 0)
		{
		  minjes = std::min(maxjes, minjes);
		  maxjes = 0;
		}

	      std::cout << "iter: " << niter << " / range "<< irange << " / JES " << maxjes << " /  JER " << maxjer << " / HALF " << h_sys_half[irange][niter]->GetBinContent(i+1) << std::endl; 

	      h_sys_jer->SetBinContent(i+1, maxjer);
	      h_sys_jer_flip->SetBinContent(i+1, -1*fabs(minjer));

	      h_sys_jes->SetBinContent(i+1, maxjes);
	      h_sys_jes_flip->SetBinContent(i+1, -1*fabs(minjes));


	      float toal_neg = sqrt(TMath::Power(minjer, 2) +
				    TMath::Power(minjes, 2) +
				    TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1), 2) +
				    //TMath::Power(h_sys_herwig[irange][niter]->GetBinContent(i+1) - 1, 2) +
				    TMath::Power(h_sys_half[irange][niter]->GetBinContent(i+1), 2));
      
	      float toal_pos = sqrt(TMath::Power(maxjer, 2) +
				    TMath::Power(maxjes, 2) +
				    TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1), 2) +
				    //TMath::Power(h_sys_herwig[irange][niter]->GetBinContent(i+1), 2) +
				    TMath::Power(h_sys_half[irange][niter]->GetBinContent(i+1), 2));
      
	      h_total_sys_range[irange][niter]->SetBinContent(i+1, toal_pos);
	      h_total_sys_neg_range[irange][niter]->SetBinContent(i+1, -1*toal_neg);
	      h_total_sys->SetBinContent(i+1, toal_pos);
	      h_total_sys_flip->SetBinContent(i+1, -1* toal_neg);

	    }


 	  dlutility::SetLineAtt(h_sys_half[irange][niter], color_half, 1 + lsize_half, 1);
	  dlutility::SetMarkerAtt(h_sys_half[irange][niter], color_half, msize_half, 1);

	  dlutility::SetLineAtt(h_sys_half_flip, color_half, 1 + lsize_half, 1);
	  dlutility::SetMarkerAtt(h_sys_half_flip, color_half, msize_half, 1);

	  dlutility::SetLineAtt(h_sys_jer, color_pjer, 1 + lsize_pjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer, color_pjer, msize_pjer, 1);

	  dlutility::SetLineAtt(h_sys_jer_flip, color_pjer, 1 + lsize_pjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer_flip, color_pjer, msize_pjer, 1);

	  dlutility::SetLineAtt(h_sys_jes, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes, color_pjes, msize_pjes, 1);
	  dlutility::SetLineAtt(h_sys_jes_flip, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes_flip, color_pjes, msize_pjes, 1);

	  // dlutility::SetLineAtt(h_sys_njet[irange][niter], color_njet, 1 + lsize_njet, 1);
	  // dlutility::SetMarkerAtt(h_sys_njet[irange][niter], color_njet, msize_njet, 1);
	  // dlutility::SetLineAtt(h_sys_njet_flip, color_njet, 1 + lsize_njet, 1);
	  // dlutility::SetMarkerAtt(h_sys_njet_flip, color_njet, msize_njet, 1);

	  dlutility::SetLineAtt(h_sys_prior[irange][niter], color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior[irange][niter], color_prior, msize_prior, 1);
	  dlutility::SetLineAtt(h_sys_prior_flip, color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior_flip, color_prior, msize_prior, 1);

	  dlutility::SetLineAtt(h_sys_herwig[irange][niter], color_herwig, 1 + lsize_herwig, 1);
	  dlutility::SetMarkerAtt(h_sys_herwig[irange][niter], color_herwig, msize_herwig, 1);
	  dlutility::SetLineAtt(h_sys_herwig_flip, color_herwig, 1 + lsize_herwig, 1);
	  dlutility::SetMarkerAtt(h_sys_herwig_flip, color_herwig, msize_herwig, 1);

	  dlutility::SetLineAtt(h_total_sys, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys, kBlack, 1, 8);
	  dlutility::SetLineAtt(h_total_sys_flip, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys_flip, kBlack, 1, 8);


	  h_total_sys->SetMinimum(-0.2);
	  h_total_sys->SetMaximum(0.3);

	  /* if (irange == 0) */
	  /*   { */
	  /*     dlutility::SetFont(h_total_sys, 42, 0.06, 0.05, 0.05, 0.05); */
	  /*     h_total_sys->GetXaxis()->SetTitleOffset(1.3); */
	  /*     h_total_sys->GetYaxis()->SetTitleOffset(2.0); */
	  /*     h_total_sys->GetXaxis()->SetLabelOffset(0.02); */
	  /*     h_total_sys->GetYaxis()->SetLabelOffset(0.03); */

	  /*   } */
	  
	  dlutility::SetFont(h_total_sys, 42, 0.08, 0.07, 0.07, 0.07);
	
	  
	  h_total_sys->GetXaxis()->SetRangeUser(dxj_bins[first_bin], dxj_bins[nbins]);
	  h_total_sys->Draw("p");
	  h_total_sys_flip->Draw("p same");

	  /* h_sys_njet[irange][niter]->Draw("hist same"); */
	  /* h_sys_njet_flip->Draw("hist same"); */

	  h_sys_jer->Draw("hist same");
	  h_sys_jer_flip->Draw("hist same");


	  h_sys_jes->Draw("hist same");
	  h_sys_jes_flip->Draw("hist same");

	  //h_sys_njet[irange][niter]->Draw("hist same");
	  //h_sys_njet_flip->Draw("hist same");
	  h_sys_prior[irange][niter]->Draw("hist same");
	  //h_sys_herwig[irange][niter]->Draw("hist same");
	  h_sys_half[irange][niter]->Draw("hist same");

	  h_sys_half_flip->Draw("hist same");
	  h_sys_prior_flip->Draw("hist same");
	  //h_sys_herwig_flip->Draw("hist same")
	  if (irange == 0)
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.4, 0.85 ,0, kBlack, 0.07);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.4, 0.79, 0, kBlack, 0.07);
	    }
	  else
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.1, 0.85, 0, kBlack, 0.07);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.1, 0.79, 0, kBlack, 0.07);
	    }
	  if (irange == 0)
	    {
	      leg4->AddEntry(h_total_sys, "Total Systematics","p");
	      leg4->AddEntry(h_sys_jes_flip, "JES Systematics");
	      leg4->AddEntry(h_sys_jer_flip, "JER Systematics");
	      leg4->AddEntry(h_sys_half_flip, "Half Closure");
	      //leg4->AddEntry(h_sys_njet_flip, "Event Jet Multiplicity");
	      leg4->AddEntry(h_sys_prior_flip, "Sensitivity to Prior");
	      //leg4->AddEntry(h_sys_herwig_flip, "Unfold w/ Herwig");

	    }
	}
      cxjsys->cd(4);
      dlutility::DrawSPHENIXppsize(0.02, 0.9, 0.06);
      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", 0.1*cone_size), 0.02, 0.74, 0, kBlack, 0.06);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.02,0.63, 0, kBlack, 0.06);
      leg4->Draw();
      cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, niter));
      cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, niter));
      
    
    }
  TFile *fsys = new TFile(Form("%s/uncertainties/systematics_pp_r%02d.root",  rb.get_code_location().c_str(), cone_size), "recreate");
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  h_sys_half[irange][iter]->Write();
	  h_total_sys_range[irange][iter]->Write();
	  h_total_sys_neg_range[irange][iter]->Write();
	}
    }
  fsys->Close();
  return;
}
