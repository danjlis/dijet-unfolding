#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "dlUtility.h"
int play()
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  TF1 *fgaus = new TF1("fgaus", "gaus", 0, 2);
  TF1 *fexp = new TF1("fexp", "[0]*exp([1]*x)", 5, 70);

  TF1 *faj = new TF1("faj", "gaus", 0, 1);

  faj->SetParameters(1./sqrt(TMath::Pi()), 0, 0.10);

  fexp->SetParameters(10000000, -0.3);
  faj->Draw();

  fgaus->SetParameters(1, 1.0, 0.25);

  fexp->SetRange(5, 70);

  //float bin_edges[11];
  const int nbins = 15;
  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];
  float bin_min = 5;
  float bin_max = 70;

  float bin_xj10 = 1.0;
  float bin_xj1 = 1.0*(bin_min/bin_max);
  float alpha = TMath::Power(bin_max/bin_min, 1/(float)nbins);
  std::cout << "Alpha = " << alpha <<std::endl;
  int truth_leading_bin = 0;
  int measure_leading_bin = 0;
  int measure_subleading_bin = 0;
  int truth_subleading_bin = 0;
  float truth_leading_cut = 0;
  float reco_leading_cut = 0;
  float reco_meas_leading_cut = 0;
  float truth_subleading_cut = 0;
  float reco_subleading_cut = 0;
  float reco_meas_subleading_cut = 0;

  float measure_leading_goal = 25;
  float truth_leading_goal = 12;
  float truth_subleading_goal = 5;

  float dphicut = 3*TMath::Pi()/4.;

  // Make binning
  for (int i = 0; i < nbins+1; i++)
    {
      float ipt = bin_min*TMath::Power(alpha, (float)i);
      float ixj = bin_xj1*TMath::Power(alpha, (float)i);
      ipt_bins[i] = ipt;
      ixj_bins[i] = ixj;
      if (measure_leading_bin == 0 && ipt >= measure_leading_goal) measure_leading_bin = i;
      if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
      if (truth_subleading_bin == 0 && ipt >= truth_subleading_goal) truth_subleading_bin = i;
      std::cout << i << " : " <<  ipt << " -- " << ixj <<  std::endl;
    }

  truth_leading_cut=ipt_bins[truth_leading_bin];
  truth_subleading_cut=ipt_bins[truth_subleading_bin];

  reco_leading_cut=ipt_bins[truth_leading_bin + 3];
  reco_subleading_cut=ipt_bins[truth_subleading_bin + 3];

  reco_meas_leading_cut=ipt_bins[truth_leading_bin + 3];
  reco_meas_subleading_cut=ipt_bins[truth_subleading_bin + 3];

  measure_subleading_bin = truth_subleading_bin + 3;

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  reco_meas_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  reco_meas_subleading_cut << std::endl;

  TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
  TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
  // pure fills
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  float truth_cut = 5;
  float reco_cut = 8;
  
  for (int i = 0; i < 20000000; i++)
    {

      float smear1 = fgaus->GetRandom();
      float smear2 = fgaus->GetRandom();

      float e1 = fexp->GetRandom();

      h_truth_lead->Fill(e1);

      float aj = faj->GetRandom();

      float e2 = e1*(1. - aj)/(1 + aj);
      bool match = true;
      //      if (gRandom->Uniform() > 0.7) match = false;
      
      h_truth_lead->Fill(e1);
      h_truth_sublead->Fill(e2);

      float es1 = e1*smear1;
      float es2 = e2*smear2;

      h_reco_lead->Fill(max(es1, es2));
      h_reco_sublead->Fill(min(es1, es2));

      float maxi = std::max(es1, es2);
      float mini = std::min(es1, es2);
      float pt1_truth_bin = nbins;
      float pt2_truth_bin = nbins;
      float pt1_reco_bin = nbins;
      float pt2_reco_bin = nbins;

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
	  if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
	    {
	      pt1_reco_bin = ib;
	    }
	  if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
	    {
	      pt2_reco_bin = ib;
	    }
	}
      
      bool truth_good = (e1 >= truth_subleading_cut && e2 >= truth_subleading_cut);// && dphi_truth[isample] > dphicut);
      bool reco_good = (maxi >= reco_subleading_cut && mini >= reco_subleading_cut);// && dphi_reco[isample] > dphicut);
      if (!reco_good && !truth_good) continue;
      if (reco_good)
	{

	  h_pt1pt2->Fill(es1, es2);
	  h_pt1pt2->Fill(es2, es1);

	  h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
	  h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
	  h_reco_xj->Fill(mini/maxi);
	}
      if (truth_good)
	{
	  h_truth_xj->Fill(e2/e1);
	  h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
	  h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
	}


      //RooUnfoldResponse rooResponse_histo(*h_flat_reco_to_response_pt1pt2, *h_flat_truth_to_response, *flat_response_pt1pt2);
      if (truth_good && !reco_good)
	{
	  std::cout << "Miss" << std::endl;
	  rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin);
	  rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin);
	  continue;
	}
      if (!truth_good && reco_good)
	{
	  std::cout << "Fake" << std::endl;
	  rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin);
	  rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin);
	  continue;
	}
	  
      if (match && reco_good && truth_good)
	{
	  rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin);
	  rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin);
	  continue;
	}
      else
	{
	  // rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin);
	  // rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin);
	  // rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin);
	  // rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin);

	  continue;
	}

      std::cout << "something not right" << std::endl;

      
    }

  h_flat_reco_pt1pt2->Scale(.5);
  h_flat_truth_pt1pt2->Scale(.5);

  RooUnfoldBayes   unfold (&rooResponse, h_flat_reco_pt1pt2, 4);    // OR

  TH1D* h_flat_unfold_pt1pt2= (TH1D*) unfold.Hunfold();

  TCanvas *c = new TCanvas("c","c", 500, 500);

  dlutility::SetLineAtt(h_flat_unfold_pt1pt2, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2, kBlack, 1, 8);

  dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
  dlutility::SetLineAtt(h_flat_reco_pt1pt2, kBlue, 2, 1);

  h_flat_truth_pt1pt2->Draw("hist");
  h_flat_reco_pt1pt2->Draw("hist same");
  h_flat_unfold_pt1pt2->Draw("same p");

  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold = new TH1D("h_xj_unfold", ";x_{J};",nbins, ixj_bins);

  TH1D *h_subleading_reco[nbins];
  TH1D *h_subleading_truth[nbins];
  TH1D *h_subleading_unfold[nbins];
  for (int ib = 0; ib; ib++)
    {
      h_subleading_unfold[ib]= new TH1D(Form("h_subleading_unfold_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
      h_subleading_truth[ib]= new TH1D(Form("h_subleading_truth_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
      h_subleading_reco[ib]= new TH1D(Form("h_subleading_reco_%d", ib), ";Subleading Jet p_{T} [GeV];",nbins, ixj_bins);
    }
  
  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_truth = new TH1D("h_xjunc_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xjunc_unfold = new TH1D("h_xjunc_unfold", ";x_{J};",nbins, ixj_bins);

  for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int xbin = ib/nbins;
      int ybin = ib%nbins;
      
      int b = h_pt1pt2_reco->GetBin(xbin+1, ybin+1);

      h_pt1pt2_reco->SetBinContent(b, h_flat_reco_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_truth->SetBinContent(b, h_flat_truth_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_unfold->SetBinContent(b, h_flat_unfold_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_reco->SetBinError(b, h_flat_reco_pt1pt2->GetBinError(ib+1));
      h_pt1pt2_truth->SetBinError(b, h_flat_truth_pt1pt2->GetBinError(ib+1));
      h_pt1pt2_unfold->SetBinError(b, h_flat_unfold_pt1pt2->GetBinError(ib+1));

    }

  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 900, 300);
  cpt1pt2->Divide(3, 1);
  cpt1pt2->cd(1);
  gPad->SetLogz();
  h_pt1pt2_reco->Draw("colz");
  cpt1pt2->cd(2);
  gPad->SetLogz();
  h_pt1pt2_truth->Draw("colz");
  cpt1pt2->cd(3);
  gPad->SetLogz();
  h_pt1pt2_unfold->Draw("colz");

  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
	{
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);

	  if (ix > iy)
	    {

	      h_pt1pt2_reco->SetBinContent(bin, h_pt1pt2_reco->GetBinContent(bin)*2.);
	      h_pt1pt2_truth->SetBinContent(bin, h_pt1pt2_truth->GetBinContent(bin)*2.);
	      h_pt1pt2_unfold->SetBinContent(bin, h_pt1pt2_unfold->GetBinContent(bin)*2.);
	      h_pt1pt2_reco->SetBinError(bin, h_pt1pt2_reco->GetBinError(bin)/sqrt(2));
	      h_pt1pt2_truth->SetBinError(bin, h_pt1pt2_truth->GetBinContent(bin)/sqrt(2));
	      h_pt1pt2_unfold->SetBinError(bin, h_pt1pt2_unfold->GetBinContent(bin)/sqrt(2));

	    }
	  else if (ix < iy)
	    {	      
	      h_pt1pt2_reco->SetBinContent(bin, 0);
	      h_pt1pt2_truth->SetBinContent(bin, 0);
	      h_pt1pt2_unfold->SetBinContent(bin, 0);
	    }
	}
    }
  for (int ib = 0; ib < nbins; ib++)
    {
      h_subleading_unfold[ib] = (TH1D*) h_pt1pt2_unfold->ProjectionY(Form("h_subleading_unfold_%d", ib), ib + 1, ib+1);
      h_subleading_reco[ib] = (TH1D*) h_pt1pt2_reco->ProjectionY(Form("h_subleading_reco_%d", ib), ib + 1, ib+1);
      h_subleading_truth[ib] = (TH1D*) h_pt1pt2_truth->ProjectionY(Form("h_subleading_truth_%d", ib), ib + 1, ib+1);
    }
  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy <= ix; iy++)
	{
	  int low =  iy - ix - 1;
	  int high = iy - ix + 1;

	  int xjbin_low = nbins + low + 1;
	  int xjbin_high = nbins + low + 2;
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);
	  std::cout << ix << " -- " << iy << " --> " << low << " --> " << xjbin_high << "( " << h_xj_reco->GetBinCenter(xjbin_high)<<" ) " << " + " << xjbin_low << "( " << h_xj_reco->GetBinCenter(xjbin_low)<<" ) "<< std::endl;
	  if (ix < measure_leading_bin) continue;
	  if (iy < measure_subleading_bin) continue;
	  if (ix > nbins - 2) continue;
	  if (iy > nbins - 2) continue;

	  if (ix == iy)
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	    }

	  if (ix == iy)
	    {
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low), h_pt1pt2_truth->GetBinContent(bin));
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_truth->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_high), h_pt1pt2_truth->GetBinContent(bin)/2.);
	      h_xj_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low), h_pt1pt2_truth->GetBinContent(bin)/2.);
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_truth->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_truth->Fill(h_xj_truth->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_truth->GetBinError(bin)/sqrt(2), 2));
	    }
	
	  if (ix == iy)
	    {
	      h_xj_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_low), h_pt1pt2_unfold->GetBinContent(bin));
	      h_xjunc_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_unfold->GetBinError(bin)/sqrt(2), 2));
	    }
	  else
	    {
	      h_xj_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_high), h_pt1pt2_unfold->GetBinContent(bin)/2.);
	      h_xj_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_low), h_pt1pt2_unfold->GetBinContent(bin)/2.);
	      h_xjunc_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_unfold->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_unfold->Fill(h_xj_unfold->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_unfold->GetBinError(bin)/sqrt(2), 2));
	    }
	  
	}
    }
  for (int i = 1; i <= nbins; i++)
    {
      h_xj_reco->SetBinError(i, sqrt(h_xjunc_reco->GetBinError(i)));
      h_xj_truth->SetBinError(i, sqrt(h_xjunc_truth->GetBinError(i)));
      h_xj_unfold->SetBinError(i, sqrt(h_xjunc_unfold->GetBinError(i)));
    }


    TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold, kBlack, 1, 8);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
  // float scale = h_xj_reco->Integral();
  // h_xj_truth->Scale(scale/h_xj_truth->Integral());
  // h_xj_unfold->Scale(scale/h_xj_unfold->Integral());
  
  h_xj_truth->Draw("hist");
  h_xj_reco->Draw("hist same");
  h_xj_truth->Draw("same p");
  h_xj_reco->Draw("same p");

  h_xj_unfold->Draw("same p");

  cxj->cd(2);

  TH1D *h_truth_compare = (TH1D*) h_xj_truth->Clone();
  h_truth_compare->Divide(h_xj_unfold);
  dlutility::SetLineAtt(h_truth_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_truth_compare, kBlack, 1,8);

  h_truth_compare->SetMinimum(0.5);
  h_truth_compare->SetMaximum(1.5);

  h_truth_compare->Draw("p");

  TCanvas *c_close = new TCanvas("c_close","c_close", 500, 500);

  TH2D *h_close_2 = (TH2D*) h_pt1pt2_unfold->Clone();
  h_close_2->Divide(h_pt1pt2_truth);
  h_close_2->Draw("colz");
  
  return 0;
  
}
