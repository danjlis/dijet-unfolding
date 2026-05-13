#include <iostream>
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

using std::cout;
using std::endl;

bool default1 = false;

struct jet
{
  int id;
  float pt = 0;
  float t = 0;
  float eta = 0;
  float eta_det = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;

  void print()
  {
    std::cout << "Jet " << id << std::endl;
    std::cout << "    pt/eta/phi : " << pt << " / " << eta << " / " << phi << std::endl;
  };

};

const bool Debug = false;
const bool Nuclear = false;
const float dRcut = 0.75;
static float cone_size = .3;
const float truth_cut = 3;
const float reco_cut = 3;
const float etacut = 0.7;
const float isocut = 1.0;
float etacut_cone = 1.1;
const float dphicut = 3*TMath::Pi()/4.;
const float dphicutloose = TMath::Pi()/2.;
const double deltatimecut=6.0;
const double leadtimecut=9.0;
const float vertex_cut = 60;
const float pt2cut = 8.20847;
const float pt1cut = 18.2836;

float getDPHI(float phi1, float phi2);
bool check_dijet_reco(std::vector<struct jet> myrecojets, int cone_size);
bool passtime(std::vector<struct jet> myrecojets);
double getDR(struct jet j1, struct jet j2)
{ 
  double dphi = fabs(j1.phi - j2.phi);
  if (dphi > TMath::Pi())
    {
      dphi = 2*TMath::Pi() - dphi;
    }

  double dR = sqrt(TMath::Power(j1.eta - j2.eta, 2) + TMath::Power(dphi, 2));
  return dR;
}
void getPileup(const int version = 39)
{

  read_binning rb("binning.config");
  Int_t read_nbins = rb.get_nbins();
  const int nbins = read_nbins;
  const int nbins_pt = read_nbins + 1;
  std::cout << nbins_pt << std::endl;
  float ipt_bins[nbins_pt+1];
  double dpt_bins[nbins_pt+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      dpt_bins[i] = (double) ipt_bins[i];
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }
  ipt_bins[nbins_pt] = 100;
  dpt_bins[nbins_pt] = 100;

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float reco_leading_bin = rb.get_reco_leading_bin();
  float reco_subleading_bin = rb.get_reco_subleading_bin();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
  
  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  int exbin = 1;
  int measure_bins[10] = {0};
  int mbins = 3;  
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }


  std::string infile = rb.get_tntuple_location() + "/TREE_JET_SKIM_r04_v11_" + std::to_string(version) + "_ana509_MDC2-00000028.root";
  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "HIST_JET_PILEUP";
  std::string oldname = "TREE_JET";
  mycopy.replace(mycopy.find(oldname), oldname.length(), newname);
  std::string outfile_name = rb.get_code_location() + "/efficiencies/" + mycopy;

  TFile *f = new TFile(infile.c_str(),"r");

  TTree *ttree = (TTree *) f->Get("ttree");
  
  
  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_calib  = 0;
  std::vector<float> *reco_jet_emcal  = 0;
  std::vector<float> *reco_jet_e  = 0;
  std::vector<float> *reco_jet_eta  = 0;
  std::vector<float> *reco_jet_eta_det  = 0;
  std::vector<float> *reco_jet_phi  = 0;

  float mbd_vertex_z;
  int mbd_hit;

  int conesize = 4;
  ttree->SetBranchAddress(Form("jet_pt_%d", conesize), &reco_jet_pt);
  ttree->SetBranchAddress(Form("jet_pt_calib_%d", conesize), &reco_jet_pt_calib);
  ttree->SetBranchAddress(Form("jet_emcal_%d", conesize), &reco_jet_emcal);
  ttree->SetBranchAddress(Form("jet_e_%d", conesize), &reco_jet_e);
  ttree->SetBranchAddress(Form("jet_eta_%d", conesize), &reco_jet_eta);
  ttree->SetBranchAddress(Form("jet_eta_det_%d", conesize), &reco_jet_eta_det);
  ttree->SetBranchAddress(Form("jet_phi_%d", conesize), &reco_jet_phi);


  ttree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  ttree->SetBranchAddress("mbd_hit", &mbd_hit);
 
  int nevents = ttree->GetEntries();

  TH2D *h_pt1pt2 = new TH2D(Form("h_pt1pt2_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);  
  TEfficiency *he_mbd = new TEfficiency("he_mbd",";pt;MBD Efficiency",20, 0, 100); 

  short reco_good = 0;
  std::vector<struct jet> myrecojets;
  for (int i = 0; i < nevents; i++)
    {
      ttree->GetEntry(i);
      std::cout << "Event " << i << " \r "<< std::flush;
      
      int nrecojets = reco_jet_pt->size();
	  
      // all reco jets
      myrecojets.clear();
      reco_good = 0;
      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt_calib->at(j) < reco_cut) continue;
	  if (reco_jet_e->at(j) < 0) continue;
	  
	  struct jet tempjet;
	  tempjet.pt = reco_jet_pt_calib->at(j);
	  
	  tempjet.emcal = reco_jet_emcal->at(j);
	  tempjet.eta = reco_jet_eta->at(j);
	  tempjet.eta_det = reco_jet_eta_det->at(j);
	  tempjet.phi = reco_jet_phi->at(j);
	  tempjet.id = j;
	  myrecojets.push_back(tempjet);
	}
      
      
      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      
      reco_good = check_dijet_reco(myrecojets, 4);
      if (!reco_good) continue;
      he_mbd->Fill(mbd_hit, myrecojets.at(0).pt);
      if (std::isnan(mbd_vertex_z) || fabs(mbd_vertex_z) > vertex_cut) continue;
      if (!mbd_hit) continue;
      
      h_pt1pt2->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt);
      h_pt1pt2->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt);
    }
  
  TFile *fout = new TFile(outfile_name.c_str(), "recreate");
  h_pt1pt2->Write();
  he_mbd->Write();
  fout->Close();
}

bool check_dijet_reco(std::vector<struct jet> myrecojets, int conesize)
{

  float etacut_cone = 1.1 - 0.1*conesize;
  
  if (myrecojets.size() < 2) return false;
    
  auto leading_iter = myrecojets.begin();
  
  auto subleading_iter = myrecojets.begin() + 1;
  
  if (fabs(leading_iter->eta) > etacut_cone || fabs(subleading_iter->eta) > etacut_cone) return false;
  if (fabs(leading_iter->eta_det) > etacut_cone || fabs(subleading_iter->eta_det) > etacut_cone) return false;

  float dphir = getDPHI(leading_iter->phi, subleading_iter->phi);

  if (!(leading_iter->pt >= pt1cut && subleading_iter->pt >= pt2cut && dphir >= dphicut)) return false;

  return true;
  
}
bool passtime(std::vector<struct jet> myrecojets)
{

  if (myrecojets.size() < 2) return false;
    
  auto leading_iter = myrecojets.begin();
  
  auto subleading_iter = myrecojets.begin() + 1;

  double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
  double jetleadtime = 17.6*(leading_iter->t);
  bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < leadtimecut );
  bool passdijettime = (TMath::Abs(jetdeltatime) < deltatimecut);	  
  
  bool passbothtime = (passdijettime) && (passleadtime);
  
  if (!passbothtime) return false;
  return true;
}

float getDPHI(float phi1, float phi2)
{
  float dphitr = phi1 - phi2;
  if (dphitr > TMath::Pi())
    {
      dphitr -= 2*TMath::Pi();
    }
  if (dphitr < -1*TMath::Pi())
    {
      dphitr += 2*TMath::Pi();
    }
  return dphitr;
  
}

void calculate_and_draw()
{

  TFile *f_nop = new TFile("efficiencies/HIST_JET_PILEUP_SKIM_r04_v11_12_ana509_MDC2-00000028.root","r");

  TH2D *h_pt1pt2_nop = (TH2D*) f_nop->Get("h_pt1pt2_r04");
  h_pt1pt2_nop->SetName("h_pt1pt2_nop");

  TFile *f_pile = new TFile("efficiencies/HIST_JET_PILEUP_SKIM_r04_v11_39_ana509_MDC2-00000028.root","r");

  TH2D *h_pt1pt2_pile = (TH2D*) f_pile->Get("h_pt1pt2_r04");
  h_pt1pt2_pile->SetName("h_pt1pt2_pile");

  TH2D *h_pt1pt2_mix = (TH2D*) h_pt1pt2_nop->Clone();

  float entries = h_pt1pt2_nop->GetEntries();
  float norm_entries = h_pt1pt2_pile->GetEntries();

  TH2D *h_pt1pt2_pile_norm = (TH2D*) h_pt1pt2_pile->Clone();
  float mix = 0.22;
  h_pt1pt2_pile_norm->Scale(mix*(entries)/((1 - mix)*norm_entries));

  h_pt1pt2_mix->Add(h_pt1pt2_pile_norm);

  h_pt1pt2_nop->Scale(1./h_pt1pt2_nop->Integral(),"width");
  h_pt1pt2_pile->Scale(1./h_pt1pt2_pile->Integral(),"width");
  h_pt1pt2_mix->Scale(1./h_pt1pt2_mix->Integral(),"width");
  TH2D *h_pt1pt2_ratio = (TH2D*) h_pt1pt2_mix->Clone();
  h_pt1pt2_ratio->SetName("h_pt1pt2_ratio");
  h_pt1pt2_ratio->Divide(h_pt1pt2_nop);

  dlutility::SetyjPadStyle();
  gStyle->SetOptStat(0);
  
  TCanvas *c = new TCanvas("c","c", 1200, 1000);
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLogz();
  h_pt1pt2_nop->SetTitle("No Pileup;p_{T1};p_{T2};arb.");
  gPad->SetRightMargin(0.2);
  h_pt1pt2_nop->Draw("colz");
  c->cd(2);
  gPad->SetLogz();
  h_pt1pt2_pile->SetTitle("Pileup;p_{T1};p_{T2};arb.");
  gPad->SetRightMargin(0.2);
  h_pt1pt2_pile->Draw("colz");
  c->cd(3);
  gPad->SetLogz();
  h_pt1pt2_mix->SetTitle("Mix Pileup;p_{T1};p_{T2};arb.");
  gPad->SetRightMargin(0.2);
  h_pt1pt2_mix->Draw("colz");
  c->cd(4);
  h_pt1pt2_ratio->SetTitle("Mix Ratio;p_{T1};p_{T2};Mix/NOP");
  gPad->SetRightMargin(0.2);
  h_pt1pt2_ratio->Draw("colz");
}
