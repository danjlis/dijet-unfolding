#include "dlUtility.h"
void drawMatching()
{
  std::string infile3 = "TREE_MATCH_v6_30_new_ProdA_2024-00000021.root";
  TFile *f3 = new TFile(infile3.c_str(), "r");

  TEfficiency *h_fake_rate_v_reco30 = (TEfficiency*) f3->Get("h_fake_rate_v_reco");
  h_fake_rate_v_reco30->SetName("h_fake_rate_v_reco30");
  
  std::string infile = "TREE_MATCH_v6_10_new_ProdA_2024-00000021.root";
  TFile *f = new TFile(infile.c_str(), "r");

  TEfficiency *h_fake_rate_v_reco = (TEfficiency*) f->Get("h_fake_rate_v_reco");
  h_fake_rate_v_reco->Add(*h_fake_rate_v_reco30);
  TH1D *h  = new TH1D("h","; Reconstructed Jet p_{T} [GeV]; 1 - Fake Rate", 1, 0, 50);
  h->SetMaximum(1.15);
  h->SetMinimum(0.7);
  dlutility::SetyjPadStyle();
  TCanvas *c = new TCanvas("c","c", 500, 500);
  gStyle->SetOptStat(0);
  h->Draw();
  h_fake_rate_v_reco->Draw("same p");
  dlutility::DrawSPHENIXpp(0.2, 0.85, 0, 1, 0, 1);
  dlutility::drawText("Pythia8", 0.2, 0.75);
  dlutility::drawText("anti-k_{t} R = 0.4", 0.9, 0.85, 1);
  dlutility::drawText("Matching #Delta R #leq 0.3", 0.9, 0.80, 1);
  TLine *line = new TLine(0, 1, 50, 1);
  dlutility::SetLineAtt(line, kBlack, 1.5, 4);

  line->Draw("same");

  c->SaveAs("fakerate.pdf");
}
