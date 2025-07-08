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

const int color_prior = kRed - 2;
const float marker_prior = 20;
const float msize_prior = 0.9;
const float lsize_prior = 1.1;

const int color_negjer = kPink + 2;
const float marker_negjer = 20;
const float msize_negjer = 0.9;
const float lsize_negjer = 1.1;

const int color_posjer = kMagenta + 2;
const float marker_posjer = 20;
const float msize_posjer = 0.9;
const float lsize_posjer = 1.1;

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
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_r%02d.root", rb.get_code_location().c_str(), cone_size),"r");

  TProfile *h_xj_rms[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_xj_rms[iter] = (TProfile*) finu->Get(Form("hp_xj_%d", iter));
    }
  TFile *fin = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d.root", rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }

  // JES +
  TFile *finpjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_posJES.root", rb.get_code_location().c_str(), cone_size),"r");

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
  // JES +
  TFile *finnjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_negJES.root", rb.get_code_location().c_str(), cone_size),"r");

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
  TFile *finposjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_posJER.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_posjer_flat_data_pt1pt2 = (TH1D*) finposjer->Get("h_data_flat_pt1pt2");
  h_posjer_flat_data_pt1pt2->SetName("h_posjer_flat_data_pt1pt2");
  TH1D *h_posjer_flat_reco_pt1pt2 = (TH1D*) finposjer->Get("h_reco_flat_pt1pt2");
  h_posjer_flat_reco_pt1pt2->SetName("h_posjer_flat_reco_pt1pt2");
  TH1D *h_posjer_flat_truth_pt1pt2 = (TH1D*) finposjer->Get("h_truth_flat_pt1pt2");
  h_posjer_flat_truth_pt1pt2->SetName("h_posjer_flat_truth_pt1pt2");
  TH1D *h_posjer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_flat_unfold_pt1pt2[iter] = (TH1D*) finposjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_posjer_flat_unfold_pt1pt2[iter]->SetName(Form("h_posjer_flat_unfold_pt1pt2_%d", iter));
    }
  // negJER
  TFile *finnegjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_negJER.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_negjer_flat_data_pt1pt2 = (TH1D*) finnegjer->Get("h_data_flat_pt1pt2");
  h_negjer_flat_data_pt1pt2->SetName("h_negjer_flat_data_pt1pt2");
  TH1D *h_negjer_flat_reco_pt1pt2 = (TH1D*) finnegjer->Get("h_reco_flat_pt1pt2");
  h_negjer_flat_reco_pt1pt2->SetName("h_negjer_flat_reco_pt1pt2");
  TH1D *h_negjer_flat_truth_pt1pt2 = (TH1D*) finnegjer->Get("h_truth_flat_pt1pt2");
  h_negjer_flat_truth_pt1pt2->SetName("h_negjer_flat_truth_pt1pt2");
  TH1D *h_negjer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_flat_unfold_pt1pt2[iter] = (TH1D*) finnegjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_negjer_flat_unfold_pt1pt2[iter]->SetName(Form("h_negjer_flat_unfold_pt1pt2_%d", iter));
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
  TFile *finnjet = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_NJET.root",  rb.get_code_location().c_str(), cone_size),"r");

  TH1D *h_njet_flat_data_pt1pt2 = (TH1D*) finnjet->Get("h_data_flat_pt1pt2");
  h_njet_flat_data_pt1pt2->SetName("h_njet_flat_data_pt1pt2");
  TH1D *h_njet_flat_reco_pt1pt2 = (TH1D*) finnjet->Get("h_reco_flat_pt1pt2");
  h_njet_flat_reco_pt1pt2->SetName("h_njet_flat_reco_pt1pt2");
  TH1D *h_njet_flat_truth_pt1pt2 = (TH1D*) finnjet->Get("h_truth_flat_pt1pt2");
  h_njet_flat_truth_pt1pt2->SetName("h_njet_flat_truth_pt1pt2");
  TH1D *h_njet_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njet_flat_unfold_pt1pt2[iter] = (TH1D*) finnjet->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_njet_flat_unfold_pt1pt2[iter]->SetName(Form("h_njet_flat_unfold_pt1pt2_%d", iter));
    }

  // PRIOR 
  TFile *finprior = new TFile(Form("%s/unfolding_hists/unfolding_hists_r%02d_PRIOR.root",  rb.get_code_location().c_str(), cone_size),"r");

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

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
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

  // pJES
  TH2D *h_pjes_pt1pt2_reco = new TH2D("h_pjes_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pjes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_pt1pt2_unfold[iter] = new TH2D("h_pjes_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
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

  TH2D *h_njes_pt1pt2_reco = new TH2D("h_njes_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_njes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_pt1pt2_unfold[iter] = new TH2D("h_njes_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
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
  TH2D *h_posjer_pt1pt2_reco = new TH2D("h_posjer_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_posjer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_pt1pt2_unfold[iter] = new TH2D("h_posjer_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_posjer_pt1pt2_unfold[iter]->SetName(Form("h_posjer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_posjer_xj_reco = new TH1D("h_posjer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_posjer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_xj_unfold[iter] = new TH1D(Form("h_posjer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_posjer_xj_reco_range[mbins];
  TH1D *h_posjer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_posjer_xj_reco_range[irange] = new TH1D(Form("h_posjer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_posjer_xj_unfold_range[irange][iter] = new TH1D(Form("h_posjer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  //NEGJER
  TH2D *h_negjer_pt1pt2_reco = new TH2D("h_negjer_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_negjer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_pt1pt2_unfold[iter] = new TH2D("h_negjer_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_negjer_pt1pt2_unfold[iter]->SetName(Form("h_negjer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_negjer_xj_reco = new TH1D("h_negjer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_negjer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_xj_unfold[iter] = new TH1D(Form("h_negjer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_negjer_xj_reco_range[mbins];
  TH1D *h_negjer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_negjer_xj_reco_range[irange] = new TH1D(Form("h_negjer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_negjer_xj_unfold_range[irange][iter] = new TH1D(Form("h_negjer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  /* //VTX */
  /* TH2D *h_vtx_pt1pt2_reco = new TH2D("h_vtx_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins); */
  /* TH2D *h_vtx_pt1pt2_unfold[niterations]; */

  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_pt1pt2_unfold[iter] = new TH2D("h_vtx_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins); */
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
  TH2D *h_njet_pt1pt2_reco = new TH2D("h_njet_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_njet_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_njet_pt1pt2_unfold[iter] = new TH2D("h_njet_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_njet_pt1pt2_unfold[iter]->SetName(Form("h_njet_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_njet_xj_reco = new TH1D("h_njet_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_njet_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njet_xj_unfold[iter] = new TH1D(Form("h_njet_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_njet_xj_reco_range[mbins];
  TH1D *h_njet_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njet_xj_reco_range[irange] = new TH1D(Form("h_njet_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njet_xj_unfold_range[irange][iter] = new TH1D(Form("h_njet_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //PRIOR
  TH2D *h_prior_pt1pt2_reco = new TH2D("h_prior_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_prior_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_pt1pt2_unfold[iter] = new TH2D("h_prior_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
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
  //posjer
  histo_opps::make_sym_pt1pt2(h_posjer_flat_reco_pt1pt2, h_posjer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_posjer_flat_unfold_pt1pt2[iter], h_posjer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_posjer_pt1pt2_reco, h_posjer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_posjer_pt1pt2_unfold[iter], h_posjer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_posjer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_posjer_xj_unfold[iter], nbins);
    }
  //negjer
  histo_opps::make_sym_pt1pt2(h_negjer_flat_reco_pt1pt2, h_negjer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_negjer_flat_unfold_pt1pt2[iter], h_negjer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_negjer_pt1pt2_reco, h_negjer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_negjer_pt1pt2_unfold[iter], h_negjer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_negjer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_negjer_xj_unfold[iter], nbins);
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
  histo_opps::make_sym_pt1pt2(h_njet_flat_reco_pt1pt2, h_njet_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_njet_flat_unfold_pt1pt2[iter], h_njet_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_njet_pt1pt2_reco, h_njet_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_njet_pt1pt2_unfold[iter], h_njet_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_njet_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_njet_xj_unfold[iter], nbins);
    }

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
  //posjer
  TH1D *h_posjer_final_xj_reco_range[mbins];
  TH1D *h_posjer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_posjer_final_xj_reco_range[irange] = new TH1D(Form("h_posjer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_posjer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_posjer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //negjer
  TH1D *h_negjer_final_xj_reco_range[mbins];
  TH1D *h_negjer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_negjer_final_xj_reco_range[irange] = new TH1D(Form("h_negjer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_negjer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_negjer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
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
  TH1D *h_njet_final_xj_reco_range[mbins];
  TH1D *h_njet_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njet_final_xj_reco_range[irange] = new TH1D(Form("h_njet_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njet_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_njet_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

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
      //posjer
      histo_opps::project_xj(h_posjer_pt1pt2_reco, h_posjer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_posjer_pt1pt2_unfold[iter], h_posjer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_posjer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_posjer_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_posjer_xj_reco_range[irange], h_posjer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_posjer_xj_unfold_range[irange][iter], h_posjer_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}
      //negjer
      histo_opps::project_xj(h_negjer_pt1pt2_reco, h_negjer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_negjer_pt1pt2_unfold[iter], h_negjer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_negjer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_negjer_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_negjer_xj_reco_range[irange], h_negjer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_negjer_xj_unfold_range[irange][iter], h_negjer_final_xj_unfold_range[irange][iter], nbins, first_xj);
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
      histo_opps::project_xj(h_njet_pt1pt2_reco, h_njet_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_njet_pt1pt2_unfold[iter], h_njet_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_njet_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_njet_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_njet_xj_reco_range[irange], h_njet_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_njet_xj_unfold_range[irange][iter], h_njet_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

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
  TH1D *h_sys_posjer[mbins][niterations];
  TH1D *h_sys_negjer[mbins][niterations];
  TH1D *h_sys_pjes[mbins][niterations];
  TH1D *h_sys_njes[mbins][niterations];

  // NJET systematics
  TH1D *h_sys_njet[mbins][niterations];
  // PRIOR systematics
  TH1D *h_sys_prior[mbins][niterations];
  // VTX reweight systeamtic
  //  TH1D *h_sys_vtx[mbins][niterations];
  // total systematics
  TH1D *h_total_sys_range[mbins][niterations];
  TH1D *h_total_sys_neg_range[mbins][niterations];

  
  TCanvas *cxjjes = new TCanvas("cxjjes","cxjjes", 500, 700);
  dlutility::ratioPanelCanvas(cxjjes);

  for (int niter = 0; niter < niterations; niter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjjes->cd(1);
	  dlutility::SetLineAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, msize_pjes, marker_pjes);
	  dlutility::SetLineAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, lsize_njes, 1);
	  dlutility::SetMarkerAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, msize_njes, marker_njes);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_pjes_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_njes_final_xj_unfold_range[irange][niter]->Draw("same p");

	  dlutility::DrawSPHENIXpp(0.22, 0.84);
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 Truth");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_pjes_final_xj_unfold_range[irange][niter], "Unfolded + 6% JES");
	  leg->AddEntry(h_njes_final_xj_unfold_range[irange][niter], "Unfolded - 6% JES");
	  leg->Draw("same");

    
	  cxjjes->cd(2);

	  TH1D *h_pjes_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_pjes_compare->Divide(h_pjes_final_xj_unfold_range[irange][niter]);
	  h_pjes_compare->SetTitle(";x_{J}; Unfold / Unfold #pm 6% JES");
	  dlutility::SetFont(h_pjes_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_pjes_compare, color_pjes, 1,1);
	  dlutility::SetMarkerAtt(h_pjes_compare, color_pjes, 1,8);

	  h_pjes_compare->SetMaximum(3.0);
	  h_pjes_compare->SetMinimum(0);
	  h_pjes_compare->Draw("p");

	  TH1D *h_njes_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_njes_compare->Divide(h_njes_final_xj_unfold_range[irange][niter]);
	  dlutility::SetLineAtt(h_njes_compare, color_njes, 1,1);
	  dlutility::SetMarkerAtt(h_njes_compare, color_njes, 1,8);
      
	  h_njes_compare->Draw("p same");

	  h_sys_njes[irange][niter] = (TH1D*) h_njes_compare->Clone();
	  h_sys_pjes[irange][niter] = (TH1D*) h_pjes_compare->Clone();
      
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjjes->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), cone_size, irange, niter));
	  cxjjes->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange, niter));
	}

  
      TCanvas *cxjjer = new TCanvas("cxjjer","cxjjer", 500, 700);
      dlutility::ratioPanelCanvas(cxjjer);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjjer->cd(1);
	  dlutility::SetLineAtt(h_posjer_final_xj_unfold_range[irange][niter], color_posjer, lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_posjer_final_xj_unfold_range[irange][niter], color_posjer, msize_posjer, marker_posjer);

	  dlutility::SetLineAtt(h_negjer_final_xj_unfold_range[irange][niter], color_negjer, lsize_negjer, 1);
	  dlutility::SetMarkerAtt(h_negjer_final_xj_unfold_range[irange][niter], color_negjer, msize_negjer, marker_negjer);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  
	  h_posjer_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_negjer_final_xj_unfold_range[irange][niter]->Draw("same p");

	  dlutility::DrawSPHENIXpp(0.22, 0.84);
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.22, 0.24, 0.4, 0.44);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 Truth");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_posjer_final_xj_unfold_range[irange][niter], "Unfolded + 5% JER");
	  leg->AddEntry(h_negjer_final_xj_unfold_range[irange][niter], "Unfolded - 5% JER");
	  leg->Draw("same");

    
	  cxjjer->cd(2);

	  TH1D *h_posjer_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_posjer_compare->Divide(h_posjer_final_xj_unfold_range[irange][niter]);
	  h_posjer_compare->SetTitle(";x_{J}; Unfold / Unfold + 5% JER");
	  dlutility::SetFont(h_posjer_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_posjer_compare, color_posjer, 1,1);
	  dlutility::SetMarkerAtt(h_posjer_compare, color_posjer, 1,8);
	  TH1D *h_negjer_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_negjer_compare->Divide(h_negjer_final_xj_unfold_range[irange][niter]);
	  h_negjer_compare->SetTitle(";x_{J}; Unfold / Unfold - 5% JER");
	  dlutility::SetFont(h_negjer_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_negjer_compare, color_negjer, 1,1);
	  dlutility::SetMarkerAtt(h_negjer_compare, color_negjer, 1,8);

	  h_posjer_compare->SetMaximum(3.0);
	  h_posjer_compare->SetMinimum(0.0);
	  h_posjer_compare->Draw("p");
	  h_negjer_compare->Draw("p same");
	  h_sys_posjer[irange][niter] = (TH1D*) h_posjer_compare->Clone();
	  h_sys_negjer[irange][niter] = (TH1D*) h_negjer_compare->Clone();
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjjer->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), cone_size, irange, niter));
	  cxjjer->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange, niter));

	}
      /* TCanvas *cxjvtx = new TCanvas("cxjvtx","cxjvtx", 500, 700); */
      /* dlutility::ratioPanelCanvas(cxjvtx); */

      /* for (int irange = 0; irange < mbins; irange++) */
      /* 	{ */
      /* 	  cxjvtx->cd(1); */
      /* 	  dlutility::SetLineAtt(h_vtx_final_xj_unfold_range[irange][niter], color_vtx, lsize_vtx, 1); */
      /* 	  dlutility::SetMarkerAtt(h_vtx_final_xj_unfold_range[irange][niter], color_vtx, msize_vtx, marker_vtx); */

      /* 	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1); */
      /* 	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold); */

      /* 	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1); */
      /* 	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia); */

      /* 	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04); */
      /* 	  h_final_xj_truth_range[irange]->SetMaximum(4.5); */
      /* 	  h_final_xj_truth_range[irange]->SetMinimum(0); */
      /* 	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}"); */

      /* 	  h_final_xj_truth_range[irange]->Draw("p"); */
      /* 	  h_final_xj_unfold_range[irange][niter]->Draw("same p"); */
      /* 	  h_vtx_final_xj_unfold_range[irange][niter]->Draw("same p"); */

      /* 	  dlutility::DrawSPHENIXpp(0.22, 0.84); */
      /* 	  dlutility::drawText("anti-k_{t} R = 0.4", 0.22, 0.74); */
      /* 	  dlutility::drawText(Form("%2.1f #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69); */
      /* 	  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64); */
      /* 	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59); */
      /* 	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54); */
      /* 	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49); */
      /* 	  TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44); */
      /* 	  leg->SetLineWidth(0); */
      /* 	  leg->SetTextSize(0.04); */
      /* 	  leg->SetTextFont(42); */
      /* 	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 Truth"); */
      /* 	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded"); */
      /* 	  leg->AddEntry(h_vtx_final_xj_unfold_range[irange][niter], "Unfolded + z_{vtx} reweight"); */
      /* 	  leg->Draw("same"); */

    
      /* 	  cxjvtx->cd(2); */

      /* 	  TH1D *h_vtx_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone(); */
      /* 	  h_vtx_compare->Divide(h_vtx_final_xj_unfold_range[irange][niter]); */
      /* 	  h_vtx_compare->SetTitle(";x_{J}; Unfold / Unfold + z_{vtx} reweight"); */
      /* 	  dlutility::SetFont(h_vtx_compare, 42, 0.06); */
      /* 	  dlutility::SetLineAtt(h_vtx_compare, color_vtx, 1,1); */
      /* 	  dlutility::SetMarkerAtt(h_vtx_compare, color_vtx, 1,8); */

      /* 	  h_vtx_compare->SetMaximum(2.0); */
      /* 	  h_vtx_compare->SetMinimum(0.0); */
      /* 	  h_vtx_compare->Draw("p"); */
      /* 	  h_sys_vtx[irange][niter] = (TH1D*) h_vtx_compare->Clone(); */
      /* 	  TLine *line = new TLine(0.1, 1, 1, 1); */
      /* 	  line->SetLineStyle(4); */
      /* 	  line->SetLineColor(kRed + 3); */
      /* 	  line->SetLineWidth(2); */
      /* 	  line->Draw("same"); */
      /* 	  cxjvtx->SaveAs(Form("systematic_plots/h_VTX_xj_unfolded_range_%d_iter_%d.png", irange, niter)); */
      /* 	  cxjvtx->SaveAs(Form("systematic_plots/h_VTX_xj_unfolded_range_%d_iter_%d.pdf", irange, niter)); */
      /* 	} */

      TCanvas *cxjnjet = new TCanvas("cxjnjet","cxjnjet", 500, 700);
      dlutility::ratioPanelCanvas(cxjnjet);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjnjet->cd(1);
	  dlutility::SetLineAtt(h_njet_final_xj_unfold_range[irange][niter], color_njet, lsize_njet, 1);
	  dlutility::SetMarkerAtt(h_njet_final_xj_unfold_range[irange][niter], color_njet, msize_njet, marker_njet);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_njet_final_xj_unfold_range[irange][niter]->Draw("same p");

	  dlutility::DrawSPHENIXpp(0.22, 0.84);
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 Truth");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_njet_final_xj_unfold_range[irange][niter], "Unfolded w/ N_{Jet} weighting");
	  leg->Draw("same");

    
	  cxjnjet->cd(2);

	  TH1D *h_njet_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_njet_compare->Divide(h_njet_final_xj_unfold_range[irange][niter]);
	  h_njet_compare->SetTitle(";x_{J}; Unfold / Unfold w/ N_{Jet} Weighting");
	  dlutility::SetFont(h_njet_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_njet_compare, color_njet, 1,1);
	  dlutility::SetMarkerAtt(h_njet_compare, color_njet, 1,8);

	  h_njet_compare->SetMaximum(3.0);
	  h_njet_compare->SetMinimum(0.0);
	  h_njet_compare->Draw("p");

	  h_sys_njet[irange][niter] = (TH1D*) h_njet_compare->Clone();
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjnjet->SaveAs(Form("%s/systematic_plots/h_NJET_xj_unfolded_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), cone_size, irange, niter));
	  cxjnjet->SaveAs(Form("%s/systematic_plots/h_NJET_xj_unfolded_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange, niter));
	}


      TCanvas *cxjprior = new TCanvas("cxjprior","cxjprior", 500, 700);
      dlutility::ratioPanelCanvas(cxjprior);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjprior->cd(1);
	  dlutility::SetLineAtt(h_prior_final_xj_unfold_range[irange][niter], color_prior, lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_prior_final_xj_unfold_range[irange][niter], color_prior, msize_prior, marker_prior);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_prior_final_xj_unfold_range[irange][niter]->Draw("same p");

	  dlutility::DrawSPHENIXpp(0.22, 0.84);
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.22, 0.30, 0.4, 0.44);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 Truth");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_prior_final_xj_unfold_range[irange][niter], "Unfolded w/o Prior");
	  leg->Draw("same");
    
	  cxjprior->cd(2);

	  TH1D *h_prior_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();
	  h_prior_compare->Divide(h_prior_final_xj_unfold_range[irange][niter]);
	  h_prior_compare->SetTitle(";x_{J}; Unfold / Unfold w/ Prior");
	  dlutility::SetFont(h_prior_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_prior_compare, color_prior, 1,1);
	  dlutility::SetMarkerAtt(h_prior_compare, color_prior, 1,8);

	  h_prior_compare->SetMaximum(3.0);
	  h_prior_compare->SetMinimum(0.0);
	  h_prior_compare->Draw("p");

	  h_sys_prior[irange][niter] = (TH1D*) h_prior_compare->Clone();
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjprior->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), cone_size, irange, niter));
	  cxjprior->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange, niter));
	}

      // all systematics on one plot
      TCanvas *cxjsys = new TCanvas("cxjsys","cxjsys", 700, 500);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjsys->Clear();
	  dlutility::createCutCanvas(cxjsys);
	  cxjsys->cd(1);
	  
	  TH1D *h_sys_jer = (TH1D*) h_sys_posjer[irange][niter]->Clone();
	  TH1D *h_sys_jer_flip = (TH1D*) h_sys_posjer[irange][niter]->Clone();

	  //TH1D *h_sys_vtx_flip = (TH1D*) h_sys_vtx[irange][niter]->Clone();
	  TH1D *h_sys_njet_flip = (TH1D*) h_sys_njet[irange][niter]->Clone();
	  TH1D *h_sys_prior_flip = (TH1D*) h_sys_prior[irange][niter]->Clone();

	  TH1D *h_sys_jes = (TH1D*) h_sys_njes[irange][niter]->Clone();
	  TH1D *h_sys_jes_flip = (TH1D*) h_sys_njes[irange][niter]->Clone();	  

	  TH1D *h_total_sys = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys->Reset();
	  h_total_sys->SetName("h_total_sys");
	  h_total_sys->SetTitle(";x_{J}; Nominal / Variation");
	  TH1D *h_total_sys_flip = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_flip->Reset();
	  h_total_sys_flip->SetName("h_total_sys_flip");
	  h_total_sys_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_range[irange][niter]->Reset();
	  h_total_sys_range[irange][niter]->SetName(Form("h_total_sys_range_%d_iter_%d", irange, niter));
	  h_total_sys_neg_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_neg_range[irange][niter]->Reset();
	  h_total_sys_neg_range[irange][niter]->SetName(Form("h_total_sys_neg_range_%d_iter_%d", irange, niter));

      
	  for (int i = 0; i < h_sys_njet_flip->GetNbinsX(); i++)
	    {
	  
	      h_sys_njet_flip->SetBinContent(i+1, 2 - h_sys_njet[irange][niter]->GetBinContent(i+1));
	      h_sys_prior_flip->SetBinContent(i+1, 2 - h_sys_prior[irange][niter]->GetBinContent(i+1));

	      float maxjer = h_sys_posjer[irange][niter]->GetBinContent(i+1) - 1;
	      float minjer = h_sys_negjer[irange][niter]->GetBinContent(i+1) - 1;
	      if (maxjer < 0 && minjer > 0)
		{
		  maxjer = h_sys_negjer[irange][niter]->GetBinContent(i+1) - 1;
		  minjer = h_sys_posjer[irange][niter]->GetBinContent(i+1) - 1;

		}
	      float maxjes = h_sys_pjes[irange][niter]->GetBinContent(i+1) - 1;
	      float minjes = h_sys_njes[irange][niter]->GetBinContent(i+1) - 1;
	      if (maxjes < 0 && minjes > 0)
		{
		  maxjes = h_sys_njes[irange][niter]->GetBinContent(i+1) - 1;
		  minjes = h_sys_pjes[irange][niter]->GetBinContent(i+1) - 1;

		}
	      std::cout << "iter: " << niter << " / range "<< irange << " / JES " << maxjes << " /  JER " << maxjer << std::endl; 

	      h_sys_jer->SetBinContent(i+1, maxjer + 1);
	      h_sys_jer_flip->SetBinContent(i+1, 1-fabs(minjer));

	      h_sys_jes->SetBinContent(i+1, maxjes + 1);
	      h_sys_jes_flip->SetBinContent(i+1, 1 - fabs(minjes));

	      float toal_neg = sqrt(TMath::Power(minjer, 2) +
				    TMath::Power(minjes, 2) +
				    TMath::Power(h_sys_njet[irange][niter]->GetBinContent(i+1) - 1, 2) +
				    TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1) - 1, 2));
      
	      float toal_pos = sqrt(TMath::Power(maxjer, 2) +
				    TMath::Power(maxjes, 2) +
				    TMath::Power(h_sys_njet[irange][niter]->GetBinContent(i+1) - 1, 2) +
				    TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1) - 1, 2));
      
	      h_total_sys_range[irange][niter]->SetBinContent(i+1, toal_pos);
	      h_total_sys_neg_range[irange][niter]->SetBinContent(i+1, toal_neg);
	      h_total_sys->SetBinContent(i+1, 1 + toal_pos);
	      h_total_sys_flip->SetBinContent(i+1, 1 - toal_neg);

	    }

 	  dlutility::SetLineAtt(h_sys_jer, color_posjer, 1 + lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer, color_posjer, msize_posjer, 1);

	  dlutility::SetLineAtt(h_sys_jer_flip, color_posjer, 1 + lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer_flip, color_posjer, msize_posjer, 1);

	  dlutility::SetLineAtt(h_sys_jes, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes, color_pjes, msize_pjes, 1);
	  dlutility::SetLineAtt(h_sys_jes_flip, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes_flip, color_pjes, msize_pjes, 1);

	  dlutility::SetLineAtt(h_sys_njet[irange][niter], color_njet, 1 + lsize_njet, 1);
	  dlutility::SetMarkerAtt(h_sys_njet[irange][niter], color_njet, msize_njet, 1);
	  dlutility::SetLineAtt(h_sys_njet_flip, color_njet, 1 + lsize_njet, 1);
	  dlutility::SetMarkerAtt(h_sys_njet_flip, color_njet, msize_njet, 1);

	  dlutility::SetLineAtt(h_sys_prior[irange][niter], color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior[irange][niter], color_prior, msize_prior, 1);
	  dlutility::SetLineAtt(h_sys_prior_flip, color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior_flip, color_prior, msize_prior, 1);
	  
	  dlutility::SetLineAtt(h_total_sys, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys, kBlack, 1, 8);
	  dlutility::SetLineAtt(h_total_sys_flip, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys_flip, kBlack, 1, 8);


	  h_total_sys->SetMinimum(0.0);
	  h_total_sys->SetMaximum(3.0);
	  dlutility::SetFont(h_total_sys, 42, 0.06, 0.05, 0.05, 0.05);
	  h_total_sys->Rebin(nbins - first_bin, Form("h_rebin_sys_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_total_sys_flip->Rebin(nbins - first_bin, Form("h_rebin_sys_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jer->Rebin(nbins - first_bin, Form("h_rebin_jer_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jer_flip->Rebin(nbins - first_bin, Form("h_rebin_jer_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jes->Rebin(nbins - first_bin, Form("h_rebin_jes_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jes_flip->Rebin(nbins - first_bin, Form("h_rebin_jes_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);

	  h_sys_prior[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_prior_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_prior_flip->Rebin(nbins - first_bin, Form("h_rebin_prior_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);

	  h_sys_njet[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_njet_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_njet_flip->Rebin(nbins - first_bin, Form("h_rebin_njet_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_total_sys->GetXaxis()->SetRangeUser(dxj_bins[first_bin], dxj_bins[nbins]);
	  h_total_sys->Draw("p");
	  h_total_sys_flip->Draw("p same");

	  /* h_sys_njet[irange][niter]->Draw("hist same"); */
	  /* h_sys_njet_flip->Draw("hist same"); */

	  h_sys_jer->Draw("hist same");
	  h_sys_jer_flip->Draw("hist same");
	  h_sys_jes->Draw("hist same");
	  h_sys_jes_flip->Draw("hist same");

	  h_sys_njet[irange][niter]->Draw("hist same");
	  h_sys_njet_flip->Draw("hist same");
	  h_sys_prior[irange][niter]->Draw("hist same");
	  h_sys_prior_flip->Draw("hist same");



	  cxjsys->cd(2);
	  dlutility::DrawSPHENIXppPrelimsize(0.02, 0.9, 0.1);
	  dlutility::drawText(Form("anti-#it{k_{t}} #kern[-0.3]{#it{R}} = %0.1f", 0.1*cone_size), 0.02, 0.76, 0, kBlack, 0.1);
	  dlutility::drawText(Form("%2.1f #kern[-0.1]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.02, 0.69, 0, kBlack, 0.1);
	  dlutility::drawText(Form("p_{T,2} #kern[-0.15]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.02,0.62, 0, kBlack, 0.1);
	  dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.02,0.55, 0, kBlack, 0.1);
	  TLegend *leg4 = new TLegend(0.02, 0.20, 0.9, 0.48);
	  leg4->SetTextSize(0.08);
	  leg4->SetTextFont(42);
	  leg4->SetLineWidth(0);
	  leg4->AddEntry(h_total_sys, "Total Systematics","p");
	  leg4->AddEntry(h_sys_jes_flip, "JES Systematics");
	  leg4->AddEntry(h_sys_jer_flip, "JER Systematics");
	  leg4->AddEntry(h_sys_njet_flip, "Event Jet Multiplicity");
	  leg4->AddEntry(h_sys_prior_flip, "Sensitivity to Prior");
	  leg4->Draw();
	  cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), cone_size, irange, niter));
	  cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, irange, niter));

	}
    }
  TFile *fsys = new TFile(Form("%s/uncertainties/systematics_r%02d.root",  rb.get_code_location().c_str(), cone_size), "recreate");
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  h_total_sys_range[irange][iter]->Write();
	  h_total_sys_neg_range[irange][iter]->Write();
	}
    }
  fsys->Close();
  return;
}
