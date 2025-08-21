#include <iostream>
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

const float dRcut = 0.75;
float cone_size = 4;
const float truth_cut = 3;
const float reco_cut = 3;
const float etacut = 0.7;

const float dphicut = 3*TMath::Pi()/4.;
const float dphicutloose = 3*TMath::Pi()/4.;
const float auau_cut = 25;
const float vertex_cut = 60;

void makeDataTree(const std::string infile, const int cone_size_int = 4, const int isAuAu = 0)
{

  cone_size = (float) cone_size_int;
  std::cout << cone_size << std::endl;
  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "TNTUPLE_DIJET_r0" + std::to_string(cone_size_int);
  std::string oldname = "TREE_DIJET";
  mycopy.replace(mycopy.find(oldname), oldname.length(), newname);

  std::cout << mycopy << std::endl;

  TFile *f = new TFile(infile.c_str(), "r");
  TTree *t = (TTree*) f->Get("ttree");

  std::vector<float> *reco_jet_pt_4 = 0;
  std::vector<float> *reco_jet_emcal_4 = 0;
  std::vector<float> *reco_jet_eta_4 = 0;
  std::vector<float> *reco_jet_eta_det_4 = 0;
  std::vector<float> *reco_jet_phi_4 = 0;
  float mbd_vertex_z;
  float centrality;
  int minbias;
  ULong64_t gl1_scaled;
  if (isAuAu)
    {
      t->SetBranchAddress("minbias", &minbias);
      t->SetBranchAddress("centrality", &centrality);
      t->SetBranchAddress(Form("jet_pt_%d_sub", cone_size_int), &reco_jet_pt_4);
      t->SetBranchAddress(Form("jet_emcal_%d_sub", cone_size_int), &reco_jet_emcal_4);
      t->SetBranchAddress(Form("jet_eta_%d_sub", cone_size_int), &reco_jet_eta_4);
      t->SetBranchAddress(Form("jet_phi_%d_sub", cone_size_int), &reco_jet_phi_4);
    }
  else
    {
      t->SetBranchAddress(Form("jet_pt_%d", cone_size_int), &reco_jet_pt_4);
      t->SetBranchAddress(Form("jet_emcal_%d", cone_size_int), &reco_jet_emcal_4);
      t->SetBranchAddress(Form("jet_eta_%d", cone_size_int), &reco_jet_eta_4);
      t->SetBranchAddress(Form("jet_eta_det_%d", cone_size_int), &reco_jet_eta_det_4);
      t->SetBranchAddress(Form("jet_phi_%d", cone_size_int), &reco_jet_phi_4);
    }
  t->SetBranchAddress("gl1_scaled", &gl1_scaled);
  t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);

 std:cout << "entering" << std::endl;
  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  
  int entries = t->GetEntries();

  TFile *fout = new TFile(mycopy.c_str(), "recreate");
  TNtuple *tn_dijet = 0;
  if (isAuAu)
    {
      tn_dijet = new TNtuple("tn_dijet","matched truth and reco","pt1_reco:pt2_reco:em1_reco:em2_reco:dphi_reco:trigger:njets:thirdjetpt:centrality:mbd_vertex");
    }
  else
    {
      tn_dijet = new TNtuple("tn_dijet","matched truth and reco","pt1_reco:pt2_reco:em1_reco:em2_reco:dphi_reco:trigger:njets:thirdjetpt:mbd_vertex");
    }

  //10 GeV sample
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      std::cout << "Event: " << i << " \r" << std::flush;
      if (isAuAu && !minbias) continue;
      if (fabs(mbd_vertex_z) > vertex_cut) continue;

      int trigger_fired = 0;

      if (isAuAu)
	{

	  if (((gl1_scaled >> 10) & 0x1) == 0x1) trigger_fired = 10;;
	  if (((gl1_scaled >> 12) & 0x1) == 0x1) trigger_fired = 12;;
	  if (((gl1_scaled >> 13) & 0x1) == 0x1) trigger_fired = 13;;
	  if (((gl1_scaled >> 14) & 0x1) == 0x1) trigger_fired = 14;;
	}	
      else	
	{
	  // if (((gl1_scaled >> 17) & 0x1) == 0x1) trigger_fired = 17;;
	  if (((gl1_scaled >> 18) & 0x1) == 0x1) trigger_fired = 18;;
	  // if (((gl1_scaled >> 19) & 0x1) == 0x1) trigger_fired = 19;;
	  // if (((gl1_scaled >> 21) & 0x1) == 0x1) trigger_fired = 21;;
	  if (((gl1_scaled >> 22) & 0x1) == 0x1) trigger_fired = 22;;
	  // if (((gl1_scaled >> 23) & 0x1) == 0x1) trigger_fired = 23;;
	  // if (((gl1_scaled >> 33) & 0x1) == 0x1) trigger_fired = 33;;
	  if (((gl1_scaled >> 34) & 0x1) == 0x1) trigger_fired = 34;;
	  // if (((gl1_scaled >> 35) & 0x1) == 0x1) trigger_fired = 35;;
	}

      if (trigger_fired == 0) continue;

      int nrecojets = reco_jet_pt_4->size();
      // this is the reco index for the leading and subleading jet
      myrecojets.clear();
      if (Debug) std::cout << __LINE__ << std::endl;
      bool found_reco_dijet = false;

      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt_4->at(j) < reco_cut) continue;
	  if (fabs(reco_jet_eta_det_4->at(j)) > etacut) continue;
	  if (fabs(reco_jet_eta_4->at(j)) > etacut) continue;

	  struct jet tempjet;
	  tempjet.pt = reco_jet_pt_4->at(j);
	  tempjet.eta = reco_jet_eta_4->at(j);
	  tempjet.phi = reco_jet_phi_4->at(j);
	  tempjet.emcal = reco_jet_emcal_4->at(j);
	  tempjet.id = j;
	  myrecojets.push_back(tempjet);	  
	}

      if (myrecojets.size() == 0) continue;

      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      int njet = myrecojets.size();
      if (isAuAu && myrecojets.at(0).pt < auau_cut) continue;
      double dphir = 0;
      double thirdjetpt = 0;

      if (myrecojets.size() >= 2 && !isAuAu) 
	{
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
	      dphir = fabs(leading_iter->phi - subleading_iter->phi);
	      if (dphir > TMath::Pi())
		{
		  dphir = 2*TMath::Pi() - dphir;
		}
	      if (leading_iter->pt >= 10 && subleading_iter->pt >= 5 && dphir >= dphicut) found_reco_dijet = true;

	    }
	  if (njet >= 3) thirdjetpt = myrecojets.at(2).pt;
	  if (found_reco_dijet) tn_dijet->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, myrecojets.at(0).emcal, myrecojets.at(1).emcal, dphir, trigger_fired, njet,thirdjetpt, mbd_vertex_z);
	}
      else if (isAuAu)
	{

	  auto leading_iter = myrecojets.begin();
	  if (njet > 1)
	    {
	      auto subleading_iter = myrecojets.begin() + 1;

	      if (subleading_iter != myrecojets.end())
		{
		  dphir = fabs(leading_iter->phi - subleading_iter->phi);
		  if (dphir > TMath::Pi())
		    {
		      dphir = 2*TMath::Pi() - dphir;
		    }
		}
	    }
	  if (njet >= 3) thirdjetpt = myrecojets.at(2).pt;
	  if (njet > 1)
	    tn_dijet->Fill(myrecojets.at(0).pt, myrecojets.at(1).pt, myrecojets.at(0).emcal, myrecojets.at(1).emcal, dphir, trigger_fired, njet, thirdjetpt, centrality, mbd_vertex_z);
	  else if (njet == 1)
	    tn_dijet->Fill(myrecojets.at(0).pt, 0, myrecojets.at(0).emcal, 0, dphir, trigger_fired, njet,thirdjetpt, centrality, mbd_vertex_z);
	}
      // if (myrecojets.size() >= 2) 
      // 	{

      // 	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });


      // 	  // check dphi
      // 	  dphir = fabs(myrecojets.at(0).phi - myrecojets.at(1).phi);
      // 	  if (dphir > TMath::Pi())
      // 	    {
      // 	      dphir = 2*TMath::Pi() - dphir;
      // 	    }
      // 	  if (dphir >= dphicut) 
      // 	    {

      // 	    }
      // 	}
    }

  tn_dijet->Write();

  fout->Close();

}
