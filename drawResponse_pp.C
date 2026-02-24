#include "dlUtility.h"
#include "read_binning.h"

void drawResponse_pp(const int cone_size = 4, const int primer = 0, const std::string configfile = "binning.config")
{

  read_binning rb(configfile.c_str());

  Int_t prior_sys = rb.get_prior_sys();

  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  Int_t herwig_sys = rb.get_herwig();
  
  TString sys_name = "nominal";
  TString generator = "PYTHIA-8";
  if (prior_sys)
    sys_name = "PRIOR";
  if (herwig_sys)
    {
      sys_name = "HERWIG";
      generator = "HERWIG 7.8";
    }
  if (JER_sys < 0)
    sys_name = "negJER";

  if (JER_sys > 0)
    sys_name = "posJER";

  if (JES_sys < 0)
    sys_name = "negJES";

  if (JES_sys > 0)
    sys_name = "posJES";

  if (primer)
    {
      sys_name = "PRIMER" + std::to_string(primer) + "_" + sys_name.Data();
    }
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;


  const int nbins = read_nbins;

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


  TString responsepath = Form("response_matrices/response_matrix_pp_r%02d_%s.root",  cone_size, sys_name.Data());

  TFile *fr = new TFile(responsepath.Data(),"r");

  TH2D *h_flat_response_skim = (TH2D*) fr->Get("h_flat_response_skim");

  TH1D *h_xj_reco = (TH1D*) fr->Get("h_xj_reco");
  TH1D *h_xj_truth = (TH1D*) fr->Get("h_xj_truth");
  TH1D *h_xj_unfold[10];
  for (int iter =0; iter < 10; iter++)
    {
      h_xj_unfold[iter] = (TH1D*) fr->Get(Form("h_xj_unfold_iter%d", iter));
    }

  
  dlutility::SetyjPadStyle();
  TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 600);
  cresponseskim->SetLogz();

  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.25);
  h_flat_response_skim->SetTitle(";#it{p}_{T,1}^{reco, bin} #times nbins + #it{p}_{T,2}^{reco, bin};#it{p}_{T,1}^{truth, bin} #times nbins + #it{p}_{T,2}^{truth, bin}");
  h_flat_response_skim->Draw("col");

  dlutility::DrawSPHENIXpp(0.13, 0.95 , 0.04, 0, 1, 1, 1, "PYTHIA-8");
  dlutility::drawText(Form("%s", generator.Data()), 0.13, 0.9);
  dlutility::drawText(Form("PRIMER/SYS: %s", sys_name.Data()), 0.13, 0.85);
  dlutility::drawText(Form("anti-#it{k}_{t} R = 0.%d", cone_size), 0.87, 0.9, 1);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.87, 0.85, 1);

  cresponseskim->Print(Form("%s/unfolding_plots/response_matrix_pp_r%02d_%s.pdf", rb.get_code_location().c_str(), cone_size, sys_name.Data()));
  
  std::cout << "bins: " << h_flat_response_skim->GetXaxis()->GetNbins() << " by " << h_flat_response_skim->GetYaxis()->GetNbins() << std::endl;
  TCanvas *cresponseskimzoom = new TCanvas("fdsz","fdsz", 500, 500);
  cresponseskimzoom->SetLogz();
  
  gPad->SetRightMargin(0.05);
  h_flat_response_skim->GetXaxis()->SetRangeUser(46, 60);
  h_flat_response_skim->GetYaxis()->SetRangeUser(121, 140);
  h_flat_response_skim->Draw("col");
  dlutility::DrawSPHENIX_Prelim(0.93, 0.4);
  dlutility::drawText("Response matrix", 0.93, 0.3, 1);

  cresponseskimzoom->Print(Form("%s/unfolding_plots/response_matrix_zoom_pp_r%02d_%s.pdf", rb.get_code_location().c_str(), cone_size, sys_name.Data()));



  return;
}
