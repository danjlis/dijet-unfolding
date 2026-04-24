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
  float pt_uncalib = 0;
  float e = 0;
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
const double deltatimecut=3.0;
const double leadtimecut=6.0;

const double loose_deltatimecut=6.0;
const double loose_leadtimecut=9.0;

const float vertex_cut = 60;
const float pt2cut = 8.20847;
const float pt1cut = 18.2836;
float getDPHI(float phi1, float phi2);
bool check_dijet_reco(std::vector<struct jet> myrecojets, int cone_size);
bool pass_time_ellipse(double lead_time, double delta_time, double r1, double r2);

bool passtime(std::vector<struct jet> myrecojets, int loose);
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
void getPurity(int conesize = 4)
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

  
  std::string infile = rb.get_tntuple_location() + "/TREE_DIJET_SKIM_r04_v8_7_ana533_2025p007_v001_gl10-alltime.root";
  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "HIST_DIJET_PURITY";
  std::string oldname = "TREE_DIJET";
  mycopy.replace(mycopy.find(oldname), oldname.length(), newname);

  std::string outfile_name = rb.get_code_location() + "/efficiencies/purity_correction_r0" + std::to_string(conesize) + ".root";

  std::string corrfile_name = rb.get_code_location() + "/efficiencies/mbd_time_correction.root";

  TFile *fcalib = new TFile(corrfile_name.c_str(), "r");
  TProfile *hp_mbd_t_runnumber = (TProfile*) fcalib->Get("hp_mbd_t_runnumber");

  int lowrunbin = 47000;
  

  TFile *f = new TFile(infile.c_str(),"r");

  TTree *ttree = (TTree *) f->Get("ttree");
  
  ULong64_t gl1_scaled;
  ULong64_t gl1_live;
  int runnumber;
  
  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_calib  = 0;
  std::vector<float> *reco_jet_t  = 0;
  std::vector<float> *reco_jet_emcal  = 0;
  std::vector<float> *reco_jet_e  = 0;
  std::vector<float> *reco_jet_eta  = 0;
  std::vector<float> *reco_jet_eta_det  = 0;
  std::vector<float> *reco_jet_phi  = 0;

  float time_zero;  
  float mbd_vertex_z;
  int mbd_hit;
  double mbd_maxdelta;
  double mbd_avgdelta;
  double mbd_proddelta;
  double mbd_maxsigma;
  double mbd_avgsigma;
  double mbd_prodsigma;
  double calib_lead_time;
  double calib_sublead_time;
  double calib_mbd_time;
  
  ttree->SetBranchAddress(Form("jet_pt_%d", conesize), &reco_jet_pt);
  ttree->SetBranchAddress(Form("jet_pt_calib_%d", conesize), &reco_jet_pt_calib);
  ttree->SetBranchAddress(Form("jet_t_%d", conesize), &reco_jet_t);
  ttree->SetBranchAddress(Form("jet_emcal_%d", conesize), &reco_jet_emcal);
  ttree->SetBranchAddress(Form("jet_e_%d", conesize), &reco_jet_e);
  ttree->SetBranchAddress(Form("jet_eta_%d", conesize), &reco_jet_eta);
  ttree->SetBranchAddress(Form("jet_eta_det_%d", conesize), &reco_jet_eta_det);
  ttree->SetBranchAddress(Form("jet_phi_%d", conesize), &reco_jet_phi);


  ttree->SetBranchAddress("runnumber", &runnumber);
  
  ttree->SetBranchAddress("calib_lead_time", &calib_lead_time);
  ttree->SetBranchAddress("calib_delta_time", &calib_sublead_time);
  ttree->SetBranchAddress("calib_mbd_time", &calib_mbd_time);
      
  ttree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  ttree->SetBranchAddress("mbd_hit", &mbd_hit);
  ttree->SetBranchAddress("time_zero", &time_zero);
  ttree->SetBranchAddress("gl1_scaled", &gl1_scaled);
  ttree->SetBranchAddress("gl1_live", &gl1_live);
  ttree->SetBranchAddress("mbd_avgsigma", &mbd_avgsigma);
  ttree->SetBranchAddress("mbd_prodsigma", &mbd_prodsigma);
  ttree->SetBranchAddress("mbd_maxsigma", &mbd_maxsigma);
  ttree->SetBranchAddress("mbd_avgdelta", &mbd_avgdelta);
  ttree->SetBranchAddress("mbd_proddelta", &mbd_proddelta);
  ttree->SetBranchAddress("mbd_maxdelta", &mbd_maxdelta);
 
  int nevents = ttree->GetEntries();

  // histograms
  int nbinst = 80;
  float lowt = -20;
  float hight = 20;
  TH1D *h_mbd_t = new TH1D(Form("h_mbd_t_r0%d", conesize), ";MBD Time;counts", nbinst, lowt, hight);

  TH2D *h_t_dt = new TH2D(Form("h_t_dt_r0%d", conesize), ";lead jet t; #Delta t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_mbd_dt = new TH2D(Form("h_mbd_dt_r0%d", conesize), ";MBD t; #Delta t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_t_mbd = new TH2D(Form("h_t_mbd_r0%d", conesize), ";lead jet t; MBD Time;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_dmbd_dt = new TH2D(Form("h_dmbd_dt_r0%d", conesize), ";DMBD t; lead t - sublead t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_t_dmbd = new TH2D(Form("h_t_dmbd_r0%d", conesize), ";lead t; DMBD t;counts", nbinst, lowt, hight, nbinst, lowt, hight);

  TH1D *h_emfrac = new TH1D(Form("h_emfrac_r%02d", conesize), ";emfrac;", 120, -0.1, 1.1);
  TH1D *h_emfrac_pass = new TH1D(Form("h_emfrac_pass_r%02d", conesize), ";emfrac;", 120, -0.1, 1.1);
  TH1D *h_emfrac_bad = new TH1D(Form("h_emfrac_bad_r%02d", conesize), ";emfrac;", 200, -0.1, 1.1);
  TH1D *h_emfracdiff = new TH1D(Form("h_emfracdiff_r%02d", conesize), ";emfracdiff;", 100, 0, 1);
  TH1D *h_emfracdiff_pass = new TH1D(Form("h_emfracdiff_pass_r%02d", conesize), ";emfracdiff;", 100, 0, 1);
  TH1D *h_emfracdiff_bad = new TH1D(Form("h_emfracdiff_bad_r%02d", conesize), ";emfracdiff;", 100, 0, 1);
      
  TH2D *h_calib_t_t = new TH2D(Form("h_calib_t_t_r0%d", conesize), ";lead t; lead t - sublead t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_calib_dt_dt = new TH2D(Form("h_calib_dt_dt_r0%d", conesize), ";lead t; lead t - sublead t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  
  TH2D *h_calib_t_dt = new TH2D(Form("h_calib_t_dt_r0%d", conesize), ";lead jet t; #Delta t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_calib_mbd_dt = new TH2D(Form("h_calib_mbd_dt_r0%d", conesize), ";MBD t; lead t - sublead t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_calib_t_mbd = new TH2D(Form("h_calib_t_mbd_r0%d", conesize), ";lead t; MBD t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_calib_dmbd_dt = new TH2D(Form("h_calib_dmbd_dt_r0%d", conesize), ";DMBD t; lead t - sublead t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  TH2D *h_calib_t_dmbd = new TH2D(Form("h_calib_t_dmbd_r0%d", conesize), ";lead t; DMBD t;counts", nbinst, lowt, hight, nbinst, lowt, hight);
  
  TH1D *h_purity = new TH1D(Form("h_purity_r%02d", conesize), ";pt12;counts", nbins_pt*nbins_pt, 0, nbins_pt);

  TH2D *h_pt1pt2 = new TH2D(Form("h_pt1pt2_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_pass = new TH2D(Form("h_pt1pt2_pass_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  TProfile2D *hp_emfrac_pt1pt2_pass = new TProfile2D(Form("hp_emfrac_pt1pt2_pass_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TProfile2D *hp_emfrac_pt1pt2_bad = new TProfile2D(Form("hp_emfrac_pt1pt2_bad_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, dpt_bins, nbins_pt, dpt_bins);

  TProfile2D *hp_emfracdiff_pt1pt2_pass = new TProfile2D(Form("hp_emfracdiff_pt1pt2_pass_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, dpt_bins, nbins_pt, dpt_bins);
  TProfile2D *hp_emfracdiff_pt1pt2_bad = new TProfile2D(Form("hp_emfracdiff_pt1pt2_bad_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, dpt_bins, nbins_pt, dpt_bins);

  TH2D *h_pt1pt2_pass_e = new TH2D(Form("h_pt1pt2_pass_e_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  TH2D *h_pt1pt2_mbd_good = new TH2D(Form("h_pt1pt2_mbd_good_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_mbd_bad = new TH2D(Form("h_pt1pt2_mbd_bad_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  TH2D *h_pt1pt2_bad = new TH2D(Form("h_pt1pt2_bad_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_bad_loose = new TH2D(Form("h_pt1pt2_bad_loose_r0%d", conesize), ";pt1;pt2;counts", nbins_pt, ipt_bins, nbins_pt, ipt_bins);

  
  double bad_r1 = 15*sqrt(2);
  double bad_r2 = 60*sqrt(2);
  double pass_r1 = 2*sqrt(2);
  double pass_r2 = 6*sqrt(2);
 
  short reco_good = 0;
  std::vector<struct jet> myrecojets;
  for (int i = 0; i < nevents; i++)
    {
      ttree->GetEntry(i);
      std::cout << "Event " << i << " \r "<< std::flush;
      if (std::isnan(mbd_vertex_z) || fabs(mbd_vertex_z) > vertex_cut) continue;
      if (!mbd_hit) continue;
      
      bool trig = ((gl1_scaled >> 18) & 0x1) == 0x1 || ((gl1_scaled >> 22) & 0x1) == 0x1 || ((gl1_scaled >> 34) & 0x1) == 0x1;
      if (!trig) continue;
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
	  tempjet.pt_uncalib = reco_jet_pt->at(j);
	  tempjet.e = reco_jet_e->at(j);
	  
	  tempjet.emcal = reco_jet_emcal->at(j);
	  tempjet.eta = reco_jet_eta->at(j);
	  tempjet.eta_det = reco_jet_eta_det->at(j);
	  tempjet.phi = reco_jet_phi->at(j);
	  tempjet.id = j;
	  tempjet.t = reco_jet_t->at(j);
	  myrecojets.push_back(tempjet);
	}
      
      
      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      
      reco_good = check_dijet_reco(myrecojets, 4);
      if (!reco_good) continue;
      //if (myrecojets.at(0).pt < 30) continue;
      float pt1_reco_bin = nbins_pt;
      float pt2_reco_bin = nbins_pt;

      float maxi = myrecojets.at(0).pt;
      float mini = myrecojets.at(1).pt;
      float es1 = maxi;
      float es2 = mini;

      //if () continue;

      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
	    {
	      pt1_reco_bin = ib;
	    }
	  if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
	    {
	      pt2_reco_bin = ib;
	    }
	}

      
      h_pt1pt2->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt);
      h_pt1pt2->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt);
      auto leading_iter = myrecojets.begin();
  
      auto subleading_iter = myrecojets.begin() + 1;

      double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
      double jetleadtime = 17.6*(leading_iter->t)+2.0;

      float time_corr = hp_mbd_t_runnumber->GetBinContent(1 + runnumber - lowrunbin);
      float corr_mbd_time = calib_mbd_time - time_corr;
      h_mbd_t->Fill(corr_mbd_time);
      h_t_dt->Fill(jetleadtime, jetdeltatime);
      h_mbd_dt->Fill(time_zero, jetdeltatime);
      h_t_mbd->Fill(jetleadtime, time_zero);
      h_dmbd_dt->Fill(jetleadtime - time_zero, jetdeltatime);
      h_t_dmbd->Fill(jetleadtime, jetleadtime - time_zero);
      
      double calib_delta_time = calib_lead_time - calib_sublead_time;
      h_calib_t_dt->Fill(calib_lead_time, calib_delta_time);
      h_calib_mbd_dt->Fill(corr_mbd_time, calib_delta_time);
      h_calib_t_mbd->Fill(calib_lead_time, corr_mbd_time);
      h_calib_dmbd_dt->Fill(calib_lead_time - corr_mbd_time, calib_delta_time);
      h_calib_t_dmbd->Fill(calib_lead_time, calib_lead_time - corr_mbd_time);

      bool passes_time_cut = (fabs(calib_lead_time) < 6 && fabs(calib_delta_time) < 3);
      h_emfrac->Fill(myrecojets.at(0).emcal);
      h_emfracdiff->Fill(fabs(myrecojets.at(0).emcal - myrecojets.at(1).emcal));
      //bool passes_time_cut = pass_time_ellipse(calib_lead_time, calib_delta_time, pass_r1, pass_r2);
      //bool bad_time_cut = pass_time_ellipse(calib_lead_time, calib_delta_time, bad_r1, bad_r2);
	
      if (passes_time_cut)
	{
	  hp_emfracdiff_pt1pt2_pass->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, fabs( myrecojets.at(0).emcal-myrecojets.at(1).emcal));
	  hp_emfracdiff_pt1pt2_pass->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt, fabs( myrecojets.at(0).emcal-myrecojets.at(1).emcal));
	  hp_emfrac_pt1pt2_pass->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, myrecojets.at(0).emcal);
	  hp_emfrac_pt1pt2_pass->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt, myrecojets.at(0).emcal);

	  h_pt1pt2_pass->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt);
	  h_pt1pt2_pass->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt);
	  h_pt1pt2_pass_e->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt);
	  h_pt1pt2_pass_e->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt);	  	  
	  h_emfrac_pass->Fill(myrecojets.at(0).emcal);
	  h_emfracdiff_pass->Fill(fabs(myrecojets.at(0).emcal - myrecojets.at(1).emcal));
	}
      else
	{
	  hp_emfracdiff_pt1pt2_bad->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, fabs( myrecojets.at(0).emcal-myrecojets.at(1).emcal));
	  hp_emfracdiff_pt1pt2_bad->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt, fabs( myrecojets.at(0).emcal-myrecojets.at(1).emcal));
	  hp_emfrac_pt1pt2_bad->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, myrecojets.at(0).emcal);
	  hp_emfrac_pt1pt2_bad->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt, myrecojets.at(0).emcal);
	  h_emfrac_bad->Fill(myrecojets.at(0).emcal);
	  h_emfracdiff_bad->Fill(fabs(myrecojets.at(0).emcal - myrecojets.at(1).emcal));
	  h_pt1pt2_bad->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt);
	  h_pt1pt2_bad->Fill(myrecojets.at(1).pt, myrecojets.at(0).pt);	  
	}

    }


  
  int nx = h_pt1pt2_bad->GetNbinsX();
  int ny = h_pt1pt2_bad->GetNbinsY();
  
  TMatrixD M(nx, ny);
  
  for (int i = 1; i <= nx; i++) {
    for (int j = 1; j <= ny; j++) {
     M(i-1, j-1) = h_pt1pt2_bad->GetBinContent(i,j);
    }
  }

  TDecompSVD svd(M);

  TMatrixD U = svd.GetU();
  TMatrixD V = svd.GetV();
  TVectorD S = svd.GetSig();
  for(int i=0;i<S.GetNrows();i++)
    std::cout << i << "  " << S[i] << std::endl;

  TMatrixD M_smooth(nx, ny);
  M_smooth.Zero();
  int k = 2;
  for(int m = 0; m < k; m++) {

    for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {

	M_smooth(i,j) += S[m] * U(i,m) * V(j,m);

      }
    }

  }
  TH2D* h_pt1pt2_badSmooth = (TH2D*)h_pt1pt2_bad->Clone("h_pt1pt2_badSmooth");
  h_pt1pt2_badSmooth->Reset();

  for(int i=1;i<=nx;i++){
    for(int j=1;j<=ny;j++){

      double val = M_smooth(i-1,j-1);
      if(val < 0) val = 0;

      h_pt1pt2_badSmooth->SetBinContent(i,j,val);

    }
  }

  
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  gStyle->SetPalette(kRainBow);
  TCanvas *c2 = new TCanvas("c2","c2", 300, 700);
  c2->Divide(2,3);
  c2->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  dlutility::SetFont(h_pt1pt2, 42, 0.06);
  h_pt1pt2->DrawCopy("colz");
  c2->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  dlutility::SetFont(h_pt1pt2_pass, 42, 0.06);
  h_pt1pt2_pass->DrawCopy("colz");
  c2->cd(3);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  dlutility::SetFont(h_pt1pt2_bad, 42, 0.06);  
  h_pt1pt2_bad->DrawCopy("colz");
  c2->cd(4);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  dlutility::SetFont(h_pt1pt2_bad_loose, 42, 0.06);
  h_pt1pt2_bad_loose->DrawCopy("colz");
  c2->cd(5);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  h_pt1pt2_bad->Smooth();
  h_pt1pt2_bad->DrawCopy("colz");
  c2->cd(6);
  gPad->SetLogz();

  
  gPad->SetRightMargin(0.2);

  h_pt1pt2_badSmooth->SetMinimum(0.001);
  h_pt1pt2_badSmooth->DrawCopy("colz");

  std::cout << h_pt1pt2_badSmooth->Integral() << " / " << h_pt1pt2_bad->Integral() << std::endl;

  float a = 40*40;
  float b = 6*12;
  float A = a  - b;
  h_pt1pt2_bad->Scale(b/A);
  h_pt1pt2_badSmooth->Scale(b/A);
  TCanvas *c3 = new TCanvas("c3","c3", 500, 600);
  c3->SetRightMargin(0.2);
  //c3->SetLogz();

  TH2D *h2_purity = (TH2D*) h_pt1pt2_pass->Clone("h2_purity");
  TH2D *h2_base = (TH2D*) h_pt1pt2_pass->Clone("h2_base");
  h2_purity->Add(h_pt1pt2_badSmooth, -1);
  h2_purity->Divide(h2_base);
  h2_purity->Draw("colz");

  for (int i = 0; i < nbins_pt; i++)
    {
      for (int j = 0; j < nbins_pt; j++)
	{
	  int gbin = h2_purity->GetBin(i+1, j+1);
	  h_purity->SetBinContent(1 + nbins_pt*j + i, h2_purity->GetBinContent(gbin));
	}
    }

  TCanvas *c4 = new TCanvas("c4","c4", 500, 600);
  c4->SetRightMargin(0.2);
  c4->SetLogz();
  TH2D *h_corrected_pt1pt2 = (TH2D*) h_pt1pt2_pass->Clone();
  h_corrected_pt1pt2->Multiply(h2_purity);
  h_corrected_pt1pt2->Draw("colz");

  TCanvas *c5 = new TCanvas("c5","c5", 500, 600);
  c5->SetRightMargin(0.2);
  int irange = 2;
  TH1D *h_xj_pass_e = new TH1D("h_xj_pass_e", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_pass = new TH1D("h_xj_pass", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_corr = new TH1D("h_xj_corr", ";x_{J};",nbins, ixj_bins);

  histo_opps::project_xj(h_pt1pt2_pass, h_xj_pass, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
  histo_opps::project_xj(h_pt1pt2_pass_e, h_xj_pass_e, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
  histo_opps::project_xj(h_corrected_pt1pt2, h_xj_corr, nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 1);
  histo_opps::normalize_histo(h_xj_corr, nbins);
  histo_opps::normalize_histo(h_xj_pass, nbins);
  histo_opps::normalize_histo(h_xj_pass_e, nbins); 
  dlutility::SetLineAtt(h_xj_pass, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_pass, kBlue, 1, 8);
  dlutility::SetLineAtt(h_xj_pass_e, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_pass_e, kRed, 1, 8);
  dlutility::SetLineAtt(h_xj_corr, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_xj_corr, kBlack, 1, 8);

  h_xj_pass->Draw();
  h_xj_pass_e->Draw("same");
  h_xj_corr->Draw("same");
  
  TCanvas *c1 = new TCanvas("c1","c1", 700, 700);
  
  // Layout fractions
  double mainFrac = 0.55;
  
  // Pads
  TPad *padMain   = new TPad("padMain","",0.0,1-mainFrac,mainFrac,1.0);
  TPad *padBottom = new TPad("padBottom","",0.0,0.0,mainFrac,1-mainFrac);
  TPad *padMBD = new TPad("padMBD","",mainFrac,0.0,1.0,1-mainFrac);
  TPad *padRight  = new TPad("padRight","",mainFrac,1-mainFrac,1.0,1.0);
    
  padMain->Draw();
  padBottom->Draw();
  padRight->Draw();
  padMBD->Draw();
  padMBD->cd();
  //padMBD->SetLogy();
  h_mbd_t->Draw();
  // Draw main TH2
  padMain->cd();
  padMain->SetLogz();
  dlutility::SetFont(h_t_dt, 42, 0.06);
  h_calib_t_dt->Draw("COLZ");

  // X projection
  padBottom->cd();
  padBottom->SetLogz();

  TH1D *projX = h_t_dt->ProjectionX(Form("%s_projX",h_t_dt->GetName()));
  dlutility::SetFont(projX, 42, 0.06);
  dlutility::SetLineAtt(projX, kBlue+1, 2, 1);
  float maxt = projX->GetBinContent(projX->GetMaximumBin())*1.3;
  projX->SetMaximum(maxt);
  TH1D *projY = h_calib_t_dt->ProjectionY(Form("%s_projY",h_t_dt->GetName()));
  dlutility::SetFont(projY, 42, 0.06);
  dlutility::SetLineAtt(projY, kBlue+1, 2, 1);
  float maxdt = projY->GetBinContent(projY->GetMaximumBin())*1.3;
  projY->SetMaximum(maxdt);

  TF1 *fitgaus_t = new TF1("fitgaus_t", "gaus", -20, 20);
  //TF1 *fitgaus_dt = new TF1("fitgaus_dt", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)", -20, 20);
  // Alternatively, using the built-in "gaus" function syntax:
  TF1 *fitgaus_dt = new TF1("fitgaus_dt", "gaus(0) + gaus(3)", -20, 20);
  fitgaus_dt->SetParameters(23000, 0, 1.27, 230, 0, 10);
  //  TF1 *fitgaus_dt = new TF1("fitgaus_dt", "gaus+gaus", -20, 20, 6);
  projX->Fit("fitgaus_t");
  projY->Fit("fitgaus_dt");
  padBottom->SetLogy();
  projX->Draw("hist");
  fitgaus_t->Draw("same");
  padRight->cd();
  padRight->SetLogy();

  projY->Draw("hist");
  fitgaus_dt->Draw("same"); 

  TCanvas *ce2 = new TCanvas("ce2","ce2", 700, 600);
  hp_emfrac_pt1pt2_pass->SetTitle(";Symmetrized #it{p}_{T1,2};Symmetrized #it{p}_{T1,2}; Average EM Fraction");
  ce2->SetRightMargin(0.2);
  ce2->SetTopMargin(0.05);
  hp_emfrac_pt1pt2_pass->GetZaxis()->SetTitleOffset(1.3);
  hp_emfrac_pt1pt2_pass->Draw("colz");
  dlutility::DrawSPHENIXpp(0.52, 0.87);

  dlutility::drawText("anti-#it{k}_{t} #it{R} = 0.4", 0.52, 0.87 - 2*0.05);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.52, 0.87 - 3*0.05);
  ce2->Print("./efficiencies/hp_emfrac_r04.pdf");
  
  TCanvas *ce3 = new TCanvas("ce3","ce3", 1400, 700);
  ce3->Divide(2,1);
  ce3->cd(1);
  hp_emfracdiff_pt1pt2_pass->SetMaximum(1);
  hp_emfracdiff_pt1pt2_pass->SetMinimum(0);
  hp_emfracdiff_pt1pt2_pass->Draw("colz");
  ce3->cd(2);
  hp_emfracdiff_pt1pt2_bad->SetMaximum(1);
  hp_emfracdiff_pt1pt2_bad->SetMinimum(0);
  hp_emfracdiff_pt1pt2_bad->Draw("colz");

  TCanvas *ce = new TCanvas("ce","ce", 1400, 700);
  ce->Divide(2,1);

  ce->cd(1);
  h_emfrac_bad->Rebin(4);

  h_emfrac_pass->Scale(1./h_emfrac_pass->Integral(),"width");
  h_emfrac_bad->Scale(1./h_emfrac_bad->Integral(),"width");
  h_emfrac->Scale(1./h_emfrac->Integral(),"width");
  dlutility::SetLineAtt(h_emfrac, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_emfrac, kBlack, 1, 20);
  dlutility::SetLineAtt(h_emfrac_bad, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_emfrac_bad, kRed, 1, 20);
  dlutility::SetLineAtt(h_emfrac_pass, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_emfrac_pass, kBlue, 1, 20);
  h_emfrac->SetTitle(";Lead Jet EM Fraction; Self Norm.");
  h_emfrac->Draw();
  h_emfrac_pass->Draw("same");
  h_emfrac_bad->Draw("same");

  TLegend *l1 = new TLegend(0.25, 0.5, 0.5, 0.8);
  l1->SetLineWidth(0);
  l1->SetTextSize(0.04);
  l1->AddEntry(h_emfrac,"No timing cut");
  l1->AddEntry(h_emfrac_pass,"Passes timing cut");
  l1->AddEntry(h_emfrac_bad,"Doesn't pass timing cut");
  l1->Draw("same");
  ce->cd(2);
  h_emfracdiff_bad->Rebin(4);

  h_emfracdiff_pass->Scale(1./h_emfracdiff_pass->Integral(),"width");
  h_emfracdiff_bad->Scale(1./h_emfracdiff_bad->Integral(),"width");
  h_emfracdiff->Scale(1./h_emfracdiff->Integral(),"width");
  dlutility::SetLineAtt(h_emfracdiff, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_emfracdiff, kBlack, 1, 20);
  dlutility::SetLineAtt(h_emfracdiff_bad, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_emfracdiff_bad, kRed, 1, 20);
  dlutility::SetLineAtt(h_emfracdiff_pass, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_emfracdiff_pass, kBlue, 1, 20);
  h_emfracdiff->SetTitle(";|EM_{1} - EM_{2}|; Self Norm.");
  h_emfracdiff->Draw();
  h_emfracdiff_pass->Draw("same");
  h_emfracdiff_bad->Draw("same");
  TLegend *l2 = new TLegend(0.5, 0.6, 0.8, 0.7);
  l2->SetLineWidth(0);
  l2->SetTextSize(0.04);
  l2->AddEntry(h_emfracdiff,"No timing cut");
  l2->AddEntry(h_emfracdiff_pass,"Passes timing cut");
  l2->AddEntry(h_emfracdiff_bad,"Doesn't pass timing cut");
  l2->Draw("same");

  TFile *fout = new TFile(outfile_name.c_str(), "recreate");
  h_t_dt->Write();
  h_mbd_dt->Write();
  h_t_mbd->Write();
  h_dmbd_dt->Write();
  h_t_dmbd->Write();      
  h_calib_t_dt->Write();
  h_calib_mbd_dt->Write();
  h_calib_t_mbd->Write();
  h_calib_dmbd_dt->Write();
  h_calib_t_dmbd->Write();      
  h_pt1pt2->Write();
  h_pt1pt2_pass->Write();
  h_pt1pt2_bad->Write();
  h_purity->Write();
  h2_purity->Write();

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
bool passtime(std::vector<struct jet> myrecojets, int loose)
{

  float dtc = deltatimecut;
  float tc = leadtimecut;
  if (loose)
    {
      dtc = loose_deltatimecut;
      tc = loose_leadtimecut;
    }
  if (myrecojets.size() < 2) return false;
    
  auto leading_iter = myrecojets.begin();
  
  auto subleading_iter = myrecojets.begin() + 1;

  double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
  double jetleadtime = 17.6*(leading_iter->t);
  bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < tc );
  bool passdijettime = (TMath::Abs(jetdeltatime) < dtc);
  
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
bool pass_time_ellipse(double lead_time, double delta_time, double r1, double r2)
{
  double x = lead_time;
  double y = delta_time;
  double xr = x*cos(45) + y*sin(45);
  double yr = -x*sin(45) + y*cos(45);
  double dr = sqrt(TMath::Power(xr/r1, 2) + TMath::Power(yr/r2, 2));
  if (dr <= 1) return true;
  return false;
}
