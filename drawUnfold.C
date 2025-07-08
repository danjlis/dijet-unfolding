#include "../macros/dlUtility.h"
void drawUnfold()
{
  TFile *fin = new TFile("UnfoldHists.root","r");

  TH1D *h_matched_reco_xj = (TH1D*) fin->Get("h_matched_reco_xj");
  TH1D *h_xj_reco = (TH1D*) fin->Get("h_xj_reco");
  TH1D *h_matched_truth_xj = (TH1D*) fin->Get("h_matched_truth_xj");
  TH1D *h_xj_truth = (TH1D*) fin->Get("h_xj_truth");
  TH1D *h_xj_unfold = (TH1D*) fin->Get("h_xj_unfold");

  for (int i = 0; i < 10; i++)
    {
      h_matched_reco_xj->SetBinError(i+1, sqrt(h_matched_reco_xj->GetBinContent(i+1)));
      h_matched_truth_xj->SetBinError(i+1, sqrt(h_matched_truth_xj->GetBinContent(i+1)));
      h_xj_truth->SetBinError(i+1, sqrt(h_xj_truth->GetBinContent(i+1)));
      h_xj_reco->SetBinError(i+1, sqrt(h_xj_reco->GetBinContent(i+1)));
      h_xj_unfold->SetBinError(i+1, sqrt(h_xj_unfold->GetBinContent(i+1)));
    }
  // Now divide
  TH1D *h_reco_divide = (TH1D*) h_matched_reco_xj->Clone();
  h_reco_divide->SetName("h_reco_divide");
  h_reco_divide->Divide(h_xj_reco);
  TH1D *h_truth_divide = (TH1D*) h_matched_truth_xj->Clone();
  h_truth_divide->SetName("h_truth_divide");
  h_truth_divide->Divide(h_xj_truth);
  TH1D *h_unfold_divide = (TH1D*) h_xj_unfold->Clone();
  h_unfold_divide->SetName("h_unfold_divide");
  h_unfold_divide->Divide(h_xj_truth);
  // colors

  int color_truth = kBlack;
  int color_truth_raw = kSpring + 2;
  int color_reco = kBlue;
  int color_reco_raw = kBlue - 2;
  int color_unfold = kRed;
  int style_truth = 4;
  int style_truth_raw = 8;
  int style_reco = 25;
  int style_reco_raw = 21;
  int style_unfold = 8;
  float size_truth = 1.5;
  float size_truth_raw = 1.3;
  float size_reco = 1.6;
  float size_reco_raw = 1.4;
  float size_unfold = 1.5;

  dlutility::SetMarkerAtt(h_matched_reco_xj, color_reco_raw, size_reco_raw, style_reco_raw);
  dlutility::SetMarkerAtt(h_matched_truth_xj, color_truth_raw, size_truth_raw, style_truth_raw);
  dlutility::SetMarkerAtt(h_xj_reco, color_reco, size_reco, style_reco);
  dlutility::SetMarkerAtt(h_xj_truth, color_truth, size_truth, style_truth);
  dlutility::SetMarkerAtt(h_xj_unfold, color_unfold, size_unfold, style_unfold);
  dlutility::SetMarkerAtt(h_reco_divide, color_reco, size_reco, style_reco);
  dlutility::SetMarkerAtt(h_truth_divide, color_truth, size_truth, style_truth);

  dlutility::SetLineAtt(h_matched_reco_xj, color_reco_raw, 1, 1);
  dlutility::SetLineAtt(h_matched_truth_xj, color_truth_raw, 1, 1);
  dlutility::SetLineAtt(h_xj_reco, color_reco, 1, 1);
  dlutility::SetLineAtt(h_xj_truth, color_truth, 1, 1);
  dlutility::SetLineAtt(h_xj_unfold, color_unfold, 1, 1);
  dlutility::SetLineAtt(h_reco_divide, color_reco, 1, 1);
  dlutility::SetLineAtt(h_truth_divide, color_truth, 1, 1);

  h_reco_divide->SetTitle(";x_{J}; Raw / Proj");
  h_truth_divide->SetTitle(";x_{J}; Raw / Proj");
  TH2D *h_pt1pt2_truth = (TH2D*) fin->Get("h_pt1pt2_truth");
  TH2D *h_pt1pt2_reco = (TH2D*) fin->Get("h_pt1pt2_reco");
  TH2D *h_pt1pt2_unfold = (TH2D*) fin->Get("h_pt1pt2_unfold");

  TH2D *h_pt1pt2_truth_symm = (TH2D*) h_pt1pt2_truth->Clone();
  TH2D *h_pt1pt2_reco_symm = (TH2D*) h_pt1pt2_reco->Clone();
  TH2D *h_pt1pt2_unfold_symm = (TH2D*) h_pt1pt2_unfold->Clone();

  for (int i = 0; i < h_pt1pt2_truth->GetNbinsX(); i++)
    {
      for (int j = 0; j < i; j++)
	{
	  int bin = h_pt1pt2_truth->GetBin(i+1, j+1);
	  int binsym = h_pt1pt2_truth->GetBin(j+1, i+1);
	  float value = h_pt1pt2_truth_symm->GetBinContent(bin)/2.;
	  h_pt1pt2_truth_symm->SetBinContent(bin,value);
	  h_pt1pt2_truth_symm->SetBinContent(binsym,value);
	  value = h_pt1pt2_reco_symm->GetBinContent(bin)/2.;
	  h_pt1pt2_reco_symm->SetBinContent(bin,value);
	  h_pt1pt2_reco_symm->SetBinContent(binsym,value);
	  value = h_pt1pt2_unfold_symm->GetBinContent(bin)/2.;
	  h_pt1pt2_unfold_symm->SetBinContent(bin,value);
	  h_pt1pt2_unfold_symm->SetBinContent(binsym,value);	  
	}
    }
  dlutility::SetyjPadStyle();
  gStyle->SetOptStat(0);

  TCanvas *c_pt = new TCanvas("c_pt","c_pt", 800, 400);
  c_pt->Divide(3,1);
  c_pt->cd(1);
  gPad->SetLogz();
  h_pt1pt2_reco_symm->Draw("col");
  c_pt->cd(2);
  gPad->SetLogz();
  h_pt1pt2_truth_symm->Draw("col");
  c_pt->cd(3);
  gPad->SetLogz();
  h_pt1pt2_unfold_symm->Draw("col");


  c_pt->Print("pt1pt2.pdf");
  TCanvas *c_reco = new TCanvas("c_reco","c_reco", 500, 700);
  dlutility::ratioPanelCanvas(c_reco);
  c_reco->cd(1);
  h_matched_reco_xj->SetMaximum(h_matched_reco_xj->GetBinContent(h_matched_reco_xj->GetMaximumBin())*1.5);;
  h_matched_reco_xj->Draw();
  h_xj_reco->Draw("same");
  dlutility::DrawSPHENIXpp(0.2, 0.8, 1, 0, 1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.2, 0.68);
  TLegend *l_reco = new TLegend(0.6, 0.67, 0.78, 0.82);
  l_reco->SetLineWidth(0);
  l_reco->AddEntry(h_matched_reco_xj,"Reco raw");
  l_reco->AddEntry(h_xj_reco,"Reco proj");

  l_reco->Draw("same");

  c_reco->cd(2);
  h_reco_divide->Draw();
  c_reco->Print("reco_xj.pdf");

  TCanvas *c_truth = new TCanvas("c_truth","c_truth", 500, 700);
  dlutility::ratioPanelCanvas(c_truth);
  c_truth->cd(1);
  h_matched_truth_xj->SetMaximum(h_matched_truth_xj->GetBinContent(h_matched_truth_xj->GetMaximumBin())*1.5);;
  h_matched_truth_xj->Draw();
  h_xj_truth->Draw("same");
  dlutility::DrawSPHENIXpp(0.2, 0.8, 1, 0, 1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.2, 0.68);
  TLegend *l_truth = new TLegend(0.6, 0.67, 0.78, 0.82);
  l_truth->SetLineWidth(0);
  l_truth->AddEntry(h_matched_truth_xj,"Pythia 8 raw");
  l_truth->AddEntry(h_xj_truth,"Pythia 8 proj");

  l_truth->Draw("same");

  c_truth->cd(2);
  h_truth_divide->Draw();
  c_truth->Print("truth_xj.pdf");

  TCanvas *c_unfold = new TCanvas("c_unfold","c_unfold", 500, 700);
  dlutility::ratioPanelCanvas(c_unfold);
  c_unfold->cd(1);

  h_xj_truth->SetMaximum(h_xj_truth->GetBinContent(h_xj_truth->GetMaximumBin())*1.5);;
  h_xj_truth->Draw();
  h_xj_reco->Draw("same");
  h_xj_unfold->Draw("same");
  h_xj_truth->Draw("same");
  dlutility::DrawSPHENIXpp(0.2, 0.8, 1, 0, 1);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.2, 0.68);
  TLegend *l_unfold = new TLegend(0.6, 0.67, 0.78, 0.82);
  l_unfold->SetLineWidth(0);
  l_unfold->AddEntry(h_xj_truth,"Pythia 8");
  l_unfold->AddEntry(h_xj_reco,"Reco");
  l_unfold->AddEntry(h_xj_unfold,"Unfolded");

  l_unfold->Draw("same");
  c_unfold->cd(2);
  h_unfold_divide->Draw(); 
  c_unfold->Print("unfold_xj.pdf");
}
