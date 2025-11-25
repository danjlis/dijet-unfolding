#include "dlUtility.h"
#include "read_binning.h"

int color_sim = kRed - 2;
int color_data = kAzure - 6;

void getCentralityReweighting(const int cone_size = 4, const int centrality_bin = 9,  const std::string configfile = "binning_AA.config")
{

  bool ispp = ( centrality_bin < 0 );
  std::string system_string = (ispp?"pp":"AA_cent_" + std::to_string(centrality_bin));
  dlutility::SetyjPadStyle();
  read_binning rb(configfile.c_str());

  
  Int_t zyam_sys = rb.get_zyam_sys();
  Int_t inclusive_sys = rb.get_inclusive_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  Int_t prior_sys = rb.get_prior_sys();
  
  std::string sys_name = "nominal";
  
  if (prior_sys)
    sys_name = "PRIOR";
  
  if (zyam_sys)
    sys_name = "ZYAM";

  if (inclusive_sys)
    sys_name = "INCLUSIVE";
    
  if (JER_sys < 0)
    sys_name = "negJER";

  if (JER_sys > 0)
    sys_name = "posJER";

  if (JES_sys < 0)
    sys_name = "negJES";

  if (JES_sys > 0)
    sys_name = "posJES";
    

  TFile *f_sim = new TFile(Form("%s/response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"r");

  TH1D *h_mbd_sim = (TH1D*) f_sim->Get("h_mbd_vertex");
  h_mbd_sim->SetName("h_mbd_sim");
  h_mbd_sim->Rebin(5);
  TH1D *h_centrality_sim = (TH1D*) f_sim->Get("h_centrality");
  h_centrality_sim->SetName("h_centrality_sim");

  TFile *f_data = new TFile(Form("%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root", rb.get_code_location().c_str(), system_string.c_str(),  cone_size),"r");

  TH1D *h_mbd_data = (TH1D*) f_data->Get("h_mbd_vertex");
  h_mbd_data->SetName("h_mbd_data");
  h_mbd_data->Rebin(5);
  TH1D *h_centrality_data = (TH1D*) f_data->Get("h_centrality");
  h_centrality_data->SetName("h_centrality_data");

  h_mbd_data->Scale(1./h_mbd_data->Integral(),"width");
  h_mbd_sim->Scale(1./h_mbd_sim->Integral(),"width");

  h_centrality_data->Scale(1./h_centrality_data->Integral(),"width");
  h_centrality_sim->Scale(1./h_centrality_sim->Integral(),"width");

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

  dlutility::DrawSPHENIX(0.7, 0.8);
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

  c5->Print(Form("%s/unfolding_plots/datasim_mbd_%s_r%02d.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_mbd_%s_r%02d.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size));

  c5->cd(1);

  dlutility::SetLineAtt(h_centrality_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_centrality_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_centrality_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_centrality_data, color_data, 0.5, 8);
  dlutility::SetFont(h_centrality_data, 42, 0.05);
  h_centrality_data->SetMinimum(0);
  h_centrality_data->SetMaximum(1.0);
  h_centrality_data->SetTitle("; N_{jet}; #frac{1}{N_{pair}}#frac{dN_{pair}}{d Centrality} ");
  h_centrality_data->Draw("p");
  h_centrality_sim->Draw("same p");

  dlutility::DrawSPHENIX(0.7, 0.8);
  leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_centrality_data,"Data");
  leg1->AddEntry(h_centrality_sim,"Sim Reco");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_centrality_compare = (TH1D*) h_centrality_data->Clone();
  h_centrality_compare->SetName("h_centrality_reweight");
  h_centrality_compare->SetTitle("; Centrality ; Data/MC");
  h_centrality_compare->Divide(h_centrality_sim);

  dlutility::SetLineAtt(h_centrality_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_centrality_compare, kBlack, 1, 8);

  dlutility::SetFont(h_centrality_compare, 42, 0.07);

  h_centrality_compare->SetMinimum(0.5);
  h_centrality_compare->SetMaximum(1.5);
  h_centrality_compare->Draw();
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

  c5->Print(Form("%s/unfolding_plots/datasim_centrality_%s_r%02d.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_centrality_%s_r%02d.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size));

  TFile *fout = new TFile(Form("%s/vertex/vertex_reweight_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"recreate");
  h_compare->Write();
  fout->Write();
  fout->Close();

  if (!ispp)
    {
      fout = new TFile(Form("%s/centrality/centrality_reweight_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"recreate");
      h_centrality_compare->Write();
      fout->Write();
      fout->Close();
    }
}
