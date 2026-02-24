#include "dlUtility.h"
#include "read_binning.h"

int color_sim = kRed - 2;
int color_data = kAzure - 6;

void getVtxReweighting(const int cone_size = 4, const std::string configfile = "binning.config", const int input_generator = 0, const int unfold_generator = 0, const int full_or_half = 0)
{
  int primer = 1;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  read_binning rb(configfile.c_str());

  Int_t prior_sys = rb.get_prior_sys();

  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  Int_t herwig_sys = rb.get_herwig();
  
  std::string sys_name = "nominal";
  
  if (prior_sys)
    sys_name = "PRIOR";
  if (herwig_sys)
    sys_name = "HERWIG";
  
  if (JER_sys < 0)
    sys_name = "negJER";

  if (JER_sys > 0)
    sys_name = "posJER";

  if (JES_sys < 0)
    sys_name = "negJES";

  if (JES_sys > 0)
    sys_name = "posJES";

  if (unfold_generator == 2 && input_generator == 1)
    {
      sys_name = "PuH";
    }
  if (unfold_generator == 1 && input_generator == 1)
    {
      sys_name = "PuP";
    }
  if (unfold_generator == 2 && input_generator == 2)
    {
      sys_name = "HuH";
    }
  if (unfold_generator == 1 && input_generator == 2)
    {
      sys_name = "HuP";
    }
  if (full_or_half)
    {
      sys_name = "HALF_" + sys_name;
    }
  
  TFile *f_sim = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_PRIMER1_%s.root", rb.get_code_location().c_str(), cone_size, sys_name.c_str()),"r");

  TH1D *h_mbd_sim = (TH1D*) f_sim->Get("h_mbd_vertex");
  h_mbd_sim->SetName("h_mbd_sim");
  h_mbd_sim->Rebin(5);
  TH1D *h_njet_sim = (TH1D*) f_sim->Get("h_njets");
  h_njet_sim->SetName("h_njet_sim");
  TH2D *h_eta_sim = (TH2D*) f_sim->Get("h_eta_lead_sublead");
  h_eta_sim->SetName("h_eta_sim");
  h_eta_sim->Rebin2D(10, 10);
  
  TFile *f_data = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_PRIMER1_%s.root", rb.get_code_location().c_str(), cone_size, sys_name.c_str()),"r");

  TH1D *h_mbd_data = (TH1D*) f_data->Get("h_mbd_vertex");
  h_mbd_data->SetName("h_mbd_data");
  h_mbd_data->Rebin(5);
  TH1D *h_njet_data = (TH1D*) f_data->Get("h_njets");
  h_njet_data->SetName("h_njet_data");
  TH2D *h_eta_data = (TH2D*) f_data->Get("h_eta_lead_sublead");
  h_eta_data->SetName("h_eta_data");
  h_eta_data->Rebin2D(10, 10);
  
  h_mbd_data->Scale(1./h_mbd_data->Integral(),"width");
  h_mbd_sim->Scale(1./h_mbd_sim->Integral(),"width");
  h_njet_data->Scale(1./h_njet_data->Integral(),"width");
  h_njet_sim->Scale(1./h_njet_sim->Integral(),"width");
  h_eta_data->Scale(1./h_eta_data->Integral(0, -1, 0, -1),"width");
  h_eta_sim->Scale(1./h_eta_sim->Integral(0, -1, 0, -1),"width");

  TH1D *h_eta_lead_data = (TH1D*) h_eta_data->ProjectionX("h_eta_lead_data");
  TH1D *h_eta_lead_sim = (TH1D*) h_eta_sim->ProjectionX("h_eta_lead_sim");
  TH1D *h_eta_sublead_data = (TH1D*) h_eta_data->ProjectionY("h_eta_sublead_data");
  TH1D *h_eta_sublead_sim = (TH1D*) h_eta_sim->ProjectionY("h_eta_sublead_sim");

  TH2D *h_eta_compare = (TH2D*) h_eta_data->Clone();
  h_eta_compare->SetName("h_eta_reweight");
  h_eta_compare->Divide(h_eta_sim);
  
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
  TLegend *leg1 = new TLegend(0.22, 0.65, 0.35, 0.85);
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

  c5->Print(Form("%s/unfolding_plots/datasim_mbd_pp_r%02d_%s.png", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));
  c5->Print(Form("%s/unfolding_plots/datasim_mbd_pp_r%02d_%s.pdf", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));

  c5->cd(1);

  dlutility::SetLineAtt(h_njet_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_njet_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_njet_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_njet_data, color_data, 0.5, 8);
  dlutility::SetFont(h_njet_data, 42, 0.05);
  gPad->SetLogy(1);
  h_njet_data->SetMinimum(0.001);
  h_njet_data->SetMaximum(10);
  h_njet_data->SetTitle("; N_{jet}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dN_{jet}} ");
  h_njet_data->Draw("p");
  h_njet_sim->Draw("same p");

  dlutility::DrawSPHENIXppsize(0.7, 0.8, 0.05);
  dlutility::drawText("Jet 10 GeV", 0.7, 0.64, 0, kBlack, 0.05);
  dlutility::drawText("|z_{vtx}| < 60", 0.7, 0.56, 0, kBlack, 0.05);
  leg1 = new TLegend(0.7, 0.4, 0.85, 0.5);
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
  h_njet_compare->SetTitle("; N_{jets}; Data/MC");
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

  c5->Print(Form("%s/unfolding_plots/datasim_njet_pp_r%02d_%s.png", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));
  c5->Print(Form("%s/unfolding_plots/datasim_njet_pp_r%02d_%s.pdf", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));

  c5->cd(1);

  dlutility::SetLineAtt(h_eta_lead_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_eta_lead_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_eta_lead_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_eta_lead_data, color_data, 0.5, 8);

  dlutility::SetLineAtt(h_eta_sublead_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_eta_sublead_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_eta_sublead_sim, color_sim, 0.5, 24);
  dlutility::SetMarkerAtt(h_eta_sublead_data, color_data, 0.5, 24);

  dlutility::SetFont(h_eta_lead_data, 42, 0.05);
  h_eta_lead_data->SetTitle("; #eta ; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#eta} ");

  h_eta_lead_data->SetMaximum(18);
  h_eta_lead_data->Draw("p");
  h_eta_lead_sim->Draw("same p");
  h_eta_sublead_data->Draw("same p");
  h_eta_sublead_sim->Draw("same p");

  dlutility::DrawSPHENIXppsize(0.7, 0.8, 0.05);
  dlutility::drawText("Jet 10 GeV", 0.7, 0.64, 0, kBlack, 0.05);
  dlutility::drawText("|z_{vtx}| < 60", 0.7, 0.56, 0, kBlack, 0.05);
  leg1 = new TLegend(0.22, 0.65, 0.35, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);

  leg1->AddEntry(h_eta_lead_data,"Data L");
  leg1->AddEntry(h_eta_lead_sim,"Sim Reco L");
  leg1->AddEntry(h_eta_sublead_data,"Data S");
  leg1->AddEntry(h_eta_sublead_sim,"Sim Reco S");

  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_eta_lead_compare = (TH1D*) h_eta_lead_data->Clone();
  h_eta_lead_compare->SetName("h_eta_lead_reweight");
  h_eta_lead_compare->SetTitle("; #eta; Data/MC");
  h_eta_lead_compare->Divide(h_eta_lead_sim);

  dlutility::SetLineAtt(h_eta_lead_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_eta_lead_compare, kBlack, 1, 8);

  TH1D *h_eta_sublead_compare = (TH1D*) h_eta_sublead_data->Clone();
  h_eta_sublead_compare->SetName("h_eta_sublead_reweight");
  h_eta_sublead_compare->SetTitle("; z_{vtx} [cm]; Data/MC");
  h_eta_sublead_compare->Divide(h_eta_sublead_sim);

  dlutility::SetLineAtt(h_eta_sublead_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_eta_sublead_compare, kBlack, 1, 24);

  dlutility::SetFont(h_eta_lead_compare, 42, 0.07);

  h_eta_lead_compare->SetMinimum(0.7);
  h_eta_lead_compare->SetMaximum(1.3);
  h_eta_lead_compare->Draw();
  h_eta_sublead_compare->Draw("same");
    
  linl = new TLine(-1.1, 1, 1.1 ,1);
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

  c5->Print(Form("%s/unfolding_plots/datasim_eta_lead_pp_r%02d_%s.png", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));
  c5->Print(Form("%s/unfolding_plots/datasim_eta_lead_pp_r%02d_%s.pdf", rb.get_code_location().c_str(), cone_size, sys_name.c_str()));

  
  TFile *fout = new TFile(Form("%s/vertex/vertex_reweight_pp_r%02d_%s.root", rb.get_code_location().c_str(), cone_size, sys_name.c_str()),"recreate");
  h_compare->Write();
  h_eta_compare->Write();
  fout->Write();
  fout->Close();

  fout = new TFile(Form("%s/njet/njet_reweight_pp_r%02d_%s.root", rb.get_code_location().c_str(), cone_size, sys_name.c_str()),"recreate");
  h_njet_compare->Write();
  fout->Write();
  fout->Close();
  
}
