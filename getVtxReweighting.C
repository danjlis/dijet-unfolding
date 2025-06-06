#include "../macros/dlUtility.h"
int color_sim = kRed - 2;
int color_data = kAzure - 6;

void getVtxReweighting(const int cone_size = 4)
{
  dlutility::SetyjPadStyle();
  TFile *f_sim = new TFile(Form("response_matrices/response_matrix_r%02d_PRIMER1.root", cone_size),"r");

  TH1D *h_mbd_sim = (TH1D*) f_sim->Get("h_mbd_vertex");
  h_mbd_sim->SetName("h_mbd_sim");
  h_mbd_sim->Rebin(5);
  TH1D *h_njet_sim = (TH1D*) f_sim->Get("h_njets");
  h_njet_sim->SetName("h_njet_sim");

  TFile *f_data = new TFile(Form("unfolded_hists/unfolded_hists_r%02d_PRIMER1.root", cone_size),"r");

  TH1D *h_mbd_data = (TH1D*) f_data->Get("h_mbd_vertex");
  h_mbd_data->SetName("h_mbd_data");
  h_mbd_data->Rebin(5);
  TH1D *h_njet_data = (TH1D*) f_data->Get("h_njets");
  h_njet_data->SetName("h_njet_data");

  h_mbd_data->Scale(1./h_mbd_data->Integral(),"width");
  h_mbd_sim->Scale(1./h_mbd_sim->Integral(),"width");
  h_njet_data->Scale(1./h_njet_data->Integral(),"width");
  h_njet_sim->Scale(1./h_njet_sim->Integral(),"width");

  TCanvas *c5 = new TCanvas("c5","c5", 500, 500);
  dlutility::ratioPanelCanvas(c5, 0.4);
  c5->cd(1);
  gPad->SetBottomMargin(0.1);
  dlutility::SetLineAtt(h_mbd_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_mbd_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_mbd_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_mbd_data, color_data, 0.5, 8);
  dlutility::SetFont(h_mbd_data, 42, 0.05);
  h_mbd_data->SetMinimum(0);
  h_mbd_data->SetMaximum(0.02);
  h_mbd_data->SetTitle("; z_{vtx}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dz_{vtx}} ");
  h_mbd_data->Draw("p");
  h_mbd_sim->Draw("same p");

  dlutility::DrawSPHENIXppsize(0.7, 0.8, 0.05);
  dlutility::drawText("Jet 10 GeV", 0.7, 0.64, 0, kBlack, 0.05);
  dlutility::drawText("|z_{vtx}| < 60", 0.7, 0.56, 0, kBlack, 0.05);
  TLegend *leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_mbd_data,"Data");
  leg1->AddEntry(h_mbd_sim,"Sim Reco");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_compare = (TH1D*) h_mbd_data->Clone();
  h_compare->SetName("h_mbd_reweight");
  h_compare->SetTitle("; z_{vtx} [cm]; Data/MC");
  h_compare->Divide(h_mbd_sim);

  dlutility::SetLineAtt(h_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_compare, kBlack, 1, 8);

  dlutility::SetFont(h_compare, 42, 0.07);

  h_compare->SetMinimum(0.5);
  h_compare->SetMaximum(1.5);
  h_compare->Draw();
  TLine *linl = new TLine(-60, 1, 60 ,1);
  linl->SetLineColor(kBlack);
  linl->SetLineWidth(2);
  linl->SetLineStyle(4);
  linl->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print("unfolding_plots/datasim_mbd.png");
  c5->Print("unfolding_plots/datasim_mbd.pdf");

  c5->cd(1);

  dlutility::SetLineAtt(h_njet_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_njet_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_njet_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_njet_data, color_data, 0.5, 8);
  dlutility::SetFont(h_njet_data, 42, 0.05);
  h_njet_data->SetMinimum(0);
  h_njet_data->SetMaximum(1.0);
  h_njet_data->SetTitle("; N_{jet}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dN_{jet}} ");
  h_njet_data->Draw("p");
  h_njet_sim->Draw("same p");

  dlutility::DrawSPHENIXppsize(0.7, 0.8, 0.05);
  dlutility::drawText("Jet 10 GeV", 0.7, 0.64, 0, kBlack, 0.05);
  dlutility::drawText("|z_{vtx}| < 60", 0.7, 0.56, 0, kBlack, 0.05);
   leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_njet_data,"Data");
  leg1->AddEntry(h_njet_sim,"Sim Reco");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_njet_compare = (TH1D*) h_njet_data->Clone();
  h_njet_compare->SetName("h_njet_reweight");
  h_njet_compare->SetTitle("; z_{vtx} [cm]; Data/MC");
  h_njet_compare->Divide(h_njet_sim);

  dlutility::SetLineAtt(h_njet_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_njet_compare, kBlack, 1, 8);

  dlutility::SetFont(h_njet_compare, 42, 0.07);

  h_njet_compare->SetMinimum(0.5);
  h_njet_compare->SetMaximum(1.5);
  h_njet_compare->Draw();
  linl = new TLine(-60, 1, 60 ,1);
  linl->SetLineColor(kBlack);
  linl->SetLineWidth(2);
  linl->SetLineStyle(4);
  linl->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print("unfolding_plots/datasim_njet.png");
  c5->Print("unfolding_plots/datasim_njet.pdf");

  TFile *fout = new TFile(Form("vertex/vertex_reweight_r%02d.root", cone_size),"recreate");
  h_compare->Write();
  fout->Write();
  fout->Close();

  fout = new TFile(Form("njet/njet_reweight_r%02d.root", cone_size),"recreate");
  h_njet_compare->Write();
  fout->Write();
  fout->Close();
  
}
