#include <iostream>
#include "read_binning.h"


using std::cout;
using std::endl;

struct jet
{
  int id;
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
};

const bool Debug = false;

const float dRcut = 1.3;
float cone_size = 4;
const float truth_cut = 3;
const float reco_cut = 8;
const float etacut = 1.1;

const float dphicut = 0;//3*TMath::Pi()/4.;
const float dphicutloose = 0;//3*TMath::Pi()/4.;

const float vertex_cut = 60;

void makeDataTreeAuAu(const std::string infile, const int cone_size_int = 4, const int isAuAu = 0)
{

  read_binning rb("binning_AA.config");

  std::string infile = rb.get_data_location() + "/TREE_DIJET_v10_2_502_2024p022_v001_gl10-all.root";
  cone_size = (float) cone_size_int;
  std::cout << cone_size << std::endl;
  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "TNTUPLE_DIJET_r0" + std::to_string(cone_size_int);
  std::string oldname = "TREE_DIJET";
  mycopy.replace(mycopy.find(oldname), oldname.length(), newname);
  std::string newfile = rb.get_tntuple_location() + "/" + mycopy;
  std::cout << newfile << std::endl;

  TFile *f = new TFile(infile.c_str(), "r");
  TTree *t = (TTree*) f->Get("ttree");

  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_unsub = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;

  int centrality;
  int minbias;
  float mbd_vertex_z;
  ULong64_t gl1_scaled;
  float mbd_time_zero;
  t->SetBranchAddress("gl1_scaled", &gl1_scaled);
  t->SetBranchAddress("centrality", &centrality);
  t->SetBranchAddress("minbias", &minbias);
  t->SetBranchAddress("mbd_time_zero", &mbd_time_zero);
  t->SetBranchAddress(Form("jet_pt_%d_sub", cone_size_int), &reco_jet_pt);
  t->SetBranchAddress(Form("jet_pt_unsub_%d_sub", cone_size_int), &reco_jet_pt_unsub);
  t->SetBranchAddress(Form("jet_e_%d_sub", cone_size_int), &reco_jet_e);
  t->SetBranchAddress(Form("jet_e_unsub_%d_sub", cone_size_int), &reco_jet_e_unsub);
  t->SetBranchAddress(Form("jet_eta_%d_sub", cone_size_int), &reco_jet_eta);
  t->SetBranchAddress(Form("jet_phi_%d_sub", cone_size_int), &reco_jet_phi);
  t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  
  int entries = t->GetEntries();

  TFile *fout = new TFile(newfile.c_str(), "recreate");
  TNtuple *tn_dijet = new TNtuple("tn_dijet","matched truth and reco","pt1_reco:pt2_reco:dphi_reco:deta_reco:trigger:njets:centrality:mbd_vertex");

  std::pair<int, float> id_leaders[2];
  TF1 *fcut = new TF1("fcut","[0]+[1]*TMath::Exp(-[2]*x)",0.0,100.0);
  //  fcut->SetParameters(2.5,36.2,0.035);
  fcut->SetParameters(0.0,40,0.038);
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      std::cout << "Event: " << i << " \r" << std::flush;

      if (fabs(mbd_vertex_z) > vertex_cut) continue;
      if (!minbias) continue;

      int trigger_fired = 0;

      if ((gl1_scaled >> 10) & 0x1U) trigger_fired = 10;

      if (trigger_fired == 0) continue;

      int nrecojets = reco_jet_pt->size();
      // this is the reco index for the leading and subleading jet
      myrecojets.clear();
      if (Debug) std::cout << __LINE__ << std::endl;
      bool found_reco_dijet = false;
      id_leaders[0] = std::make_pair(0, 0);
      id_leaders[1] = std::make_pair(0, 0);
      int njet_good = 0;
      int bad_event = 0;
      float cut_value = fcut->Eval(centrality);
      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt->at(j) < reco_cut) continue;

	  if (reco_jet_e->at(j) < 0) continue;
	  if (reco_jet_e_unsub->at(j) < 0) continue;
	  if (fabs(reco_jet_eta->at(j)) > etacut) continue;
	  
	  float pt_unsub = reco_jet_pt_unsub->at(j) - reco_jet_pt->at(j);

	  if (pt_unsub > cut_value)
	    {
	      bad_event = 1;
	      continue;
	    }


	  njet_good++;
	  if (reco_jet_pt->at(j) > id_leaders[0].second)
	    {
	      id_leaders[1] = id_leaders[0];
	      id_leaders[0] = std::make_pair(j, reco_jet_pt->at(j));
	    }
	  else if (reco_jet_pt->at(j) > id_leaders[1].second)
	    {
	      id_leaders[1] = std::make_pair(j, reco_jet_pt->at(j));
	    } 	  

	}
      if (bad_event) continue;
      if (njet_good < 2) continue;
      
      struct jet tempjet1;
      tempjet1.pt = reco_jet_pt->at(id_leaders[0].first);
      tempjet1.eta = reco_jet_eta->at(id_leaders[0].first);
      tempjet1.phi = reco_jet_phi->at(id_leaders[0].first);
      tempjet1.id = id_leaders[0].first;
      myrecojets.push_back(tempjet1);
      struct jet tempjet2;
      tempjet2.pt = reco_jet_pt->at(id_leaders[1].first);
      tempjet2.eta = reco_jet_eta->at(id_leaders[1].first);
      tempjet2.phi = reco_jet_phi->at(id_leaders[1].first);
      tempjet2.id = id_leaders[1].first;;
      myrecojets.push_back(tempjet2);

      double dphir = 0;
      double detar = 0;
      auto leading_iter = myrecojets.begin();
      // Check if the reco jet satisfies the cuts
      // check dphi
      auto subleading_iter = myrecojets.begin() + 1;
      // std::find_if(leading_iter + 1, myrecojets.end(), [=] (auto a)
      // 			    {
      // 			      double temp_dphir = fabs(leading_iter->phi - a.phi);
      // 			      if (temp_dphir > TMath::Pi())
      // 				{
      // 				  temp_dphir = 2*TMath::Pi() - temp_dphir;
      // 				}
      // 			      if (temp_dphir > dphicut) return true;
      // 			      return false;
      // 			    });
      if (subleading_iter != myrecojets.end())
	{
	  detar = fabs(leading_iter->eta - subleading_iter->eta);
	  dphir = fabs(leading_iter->phi - subleading_iter->phi);
	  if (dphir > TMath::Pi())
	    {
	      dphir = 2*TMath::Pi() - dphir;
	    }
	  if (dphir >= dphicut && leading_iter->pt >= 15 && subleading_iter->pt >= 5) found_reco_dijet = true;
	  
	}
      if (found_reco_dijet) tn_dijet->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, dphir, detar, trigger_fired, nrecojets,centrality, mbd_vertex_z);
    
    }

  tn_dijet->Write();

  fout->Close();

}
