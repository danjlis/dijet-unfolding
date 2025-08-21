

#include <iostream>
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

using std::cout;
using std::endl;

struct jet
{
  int id;
  float pt = 0;
  int ptbin = 0;
  float eta = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
};

const float etacut = 1.1;
const float eta_cut_dijet = 0.8;

const float vertex_cut = 60;

void getProbability(const int cone_size = 3,  const std::string configfile = "binning_AA.config")
{

  std::cout << "Staring" << std::endl;
  read_binning rb(configfile.c_str());
  std::cout << "Read" << std::endl;
  Int_t read_nbins = rb.get_nbins();
  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();

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


  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
    
  std::string infile = "../../trees/TREE_DIJET_v10_2_502_2024p022_v001_gl10-all.root";

  TFile *f = new TFile(infile.c_str(), "r");

  TH1D *h_pt2_correction[100];
  TH1D *h_pt2_bin_correction[100];

  TH1D *h_centrality = (TH1D*) f->Get("h_centrality");

  TH1D *h_jet_spectra[100];
  for ( int i = 0; i < 100; i++)
    {
      h_pt2_correction[i] = new TH1D(Form("h_pt2_correction_%d", i),";p_{T};Efficiency", 45, 5, 50);
      h_pt2_bin_correction[i] = new TH1D(Form("h_pt2_bin_correction_%d", i),";p_{T};Efficiency", 45, 5, 50);
      h_jet_spectra[i] = (TH1D*) f->Get(Form("h_jet_spectra_%d", i));
      if (h_jet_spectra[i]->Integral() == 0) continue;

      h_jet_spectra[i]->Scale(1./h_centrality->GetBinContent(i+1),"width");
      for (int j = 0 ; j < h_pt2_correction[i]->GetNbinsX(); j++)
	{
	  double bin_low_edge = h_pt2_correction[i]->GetBinLowEdge(j+1);
	  int bin_number = h_jet_spectra[i]->FindBin(bin_low_edge);
	  double bin_integral = h_jet_spectra[i]->Integral(bin_number, -1, "width")*dphicut/TMath::Pi();
	  double prob = TMath::Exp(-1*bin_integral);
	  h_pt2_correction[i]->SetBinContent(j+1, prob);
	}
    }

  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/probability_hists_AA_r0" + std::to_string(cone_size);
  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");

  for (int i = 0 ; i < 100; i++)
    {
      h_pt2_correction[i]->Write();
    }
  fout->Close();

}
