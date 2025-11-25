#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
struct jet
{

  double pt;
  double eta;
  double phi;
};

int makeHerwig_hist(const std::string configfile = "binning.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  
  std::string filenames[3] = {"../herwig/herwig_10.root", "../herwig/herwig_20.root", "../herwig/herwig_30.root"};

  std::vector<double> *pt[3] = {0};
  std::vector<double> *phi[3] = {0};
  std::vector<double> *eta[3] = {0};
  TFile *fin[3];
  
  TTree *tn[3];
  for (int i = 0; i < 3; i++)
    {
      fin[i] = new TFile(filenames[i].c_str(), "r");
      tn[i]= (TTree*) fin[i]->Get("JetTree");

      tn[i]->SetBranchAddress("pt", &pt[i]);
      tn[i]->SetBranchAddress("eta", &eta[i]);
      tn[i]->SetBranchAddress("phi", &phi[i]);
    }
  
  float scale_factor[3] = {1, 1, 1};

  
  read_binning rb(configfile.c_str());

  Int_t minentries = rb.get_minentries();
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

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
  std::cout << " measure_cut " << measure_subleading_cut << std::endl;
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  float measure_values[10] = {0};
  float measure_sub_values[10] = {0};
  
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_values[ir] = ipt_bins[rb.get_measure_region(ir)];
      measure_sub_values[ir] = ipt_bins[rb.get_subleading_measure_region(ir)];
    }
  float herwig_sample_boundary[4] = {14, measure_values[1], measure_values[3],  70};
  TH1D *h_truth_lead_sample[3];
  for (int i = 0; i < 3; i++)
    {
      h_truth_lead_sample[i]= new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet #it{p}_{T} [GeV]; counts", 100, 0, 100);
    }
  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet #it{p}_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet #it{p}_{T} [GeV]; counts", 100, 0, 100);

  TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet #it{p}_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet #it{p}_{T} [GeV]; counts", 100, 0, 100);

  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  // pure fills
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);

  TH1D *h_linear_truth_xj[3];
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i] = new TH1D(Form("h_linear_truth_xj_%d", i),";x_{J};;", 10, 0, 1.0);
    }

  TH1D *h_truth_leading_dphi[3];
  TH1D *h_truth_subleading_dphi[3];

  for (int i = 0; i < 3; i++)
    {
      h_truth_leading_dphi[i] = new TH1D(Form("h_truth_leading_dphi_%d", i),";x_{J};;", 32, TMath::Pi()*3./4., TMath::Pi());
      h_truth_subleading_dphi[i] = new TH1D(Form("h_truth_subleading_dphi_%d", i),";x_{J};;", 32, TMath::Pi()*3./4., TMath::Pi());
    }

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);
  for (int isample = 0; isample < 3; isample++)
    {

 
  
      int entries2 = tn[isample]->GetEntries();

      std::vector<struct jet> myjets;
      for (int i = 0; i < entries2; i++)
	{
	  myjets.clear();
	  
	  tn[isample]->GetEntry(i);
	  int inrecojets = pt[isample]->size();

	  for (int ij = 0; ij < inrecojets; ij++)
	    {
	      if (pt[isample]->at(ij) < 3) continue;
	      if (fabs(eta[isample]->at(ij)) > 0.7) continue;
	      
	      struct jet tempjet;
	      tempjet.pt = pt[isample]->at(ij);
	      tempjet.eta = eta[isample]->at(ij);
	      tempjet.phi = phi[isample]->at(ij);
	      
	      myjets.push_back(tempjet);
	    }
	  if (myjets.size() < 2) continue;

      
	  std::sort(myjets.begin(), myjets.end(), [] (auto a, auto b) { return a.pt > b.pt;});
	  if (myjets.at(0).pt < herwig_sample_boundary[isample] || myjets.at(0).pt >= herwig_sample_boundary[isample+1]) continue;
	  auto leading_iter = myjets.begin();
	  auto subleading_iter = myjets.begin() + 1;

	  float maxpttruth = leading_iter->pt;
	  float pt1_truth = leading_iter->pt;
	  float pt2_truth = subleading_iter->pt;
	  float dphi_truth = fabs(leading_iter->phi - subleading_iter->phi);
	  if (dphi_truth > TMath::Pi()) dphi_truth = 2*TMath::Pi() - dphi_truth;
      
	  double event_scale = scale_factor[isample];

	  h_truth_lead_sample[0]->Fill(pt1_truth, event_scale);

	  float max_truth = 0;
	  float min_truth = 0;
	  
	  if (pt1_truth >= pt2_truth)
	    {
	      max_truth = pt1_truth;
	      min_truth = pt2_truth;
	    }
	  else
	    {
	      max_truth = pt2_truth;
	      min_truth = pt1_truth;
	    }
	    

	  float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  
	  for (int ib = 0; ib < nbins; ib++)
	    {
	      if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
		{
		  pt1_truth_bin = ib;
		}
	      if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
		{
		  pt2_truth_bin = ib;
		}
	    }
	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && dphi_truth > dphicuttruth);

	  
	  if (truth_good)
	    {
	      for (int im = 0; im < mbins; im++)
		{
		  if (e1 >= measure_values[im] && e1 < measure_values[im+1] && e2 >= measure_subleading_cut)
		    {
		      h_linear_truth_xj[im]->Fill(e2/e1, event_scale);
		      h_truth_leading_dphi[im]->Fill(dphi_truth, event_scale);
		      if (im == 1)
			{
			  for (int jm = 0; jm < mbins; jm++)
			    {
			      if (e2 >= measure_sub_values[jm] && e2 < measure_sub_values[jm+1])
				{

				  h_truth_subleading_dphi[jm]->Fill(dphi_truth, event_scale);
				}
			    }
			}
		    }
		}
	      h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
	      h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
	      h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
	      h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
	    }
	}
    }


  TString responsepath = "herwig_hist.root";
  
  TFile *fr = new TFile(responsepath.Data(),"recreate");
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i]->Write();
      h_truth_leading_dphi[i]->Write();
      h_truth_subleading_dphi[i]->Write();
    }

  h_flat_truth_to_response_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_count_flat_truth_pt1pt2->Write();
  fr->Close();

  return 0;
}

