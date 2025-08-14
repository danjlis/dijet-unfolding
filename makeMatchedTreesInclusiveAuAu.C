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

  void print()
  {
    std::cout << "Jet " << id << std::endl;
    std::cout << "    pt/eta/phi : " << pt << " / " << eta << " / " << phi << std::endl;
  };
};

const bool Debug = false;
const bool Nuclear = false;
const float dRcut = 1.0;
static float cone_size = 3;
const float truth_cut = 3;
const float reco_cut = 3;
const float etacut = 1.1;

const float dphicut = 3*TMath::Pi()/4.;
const float dphicutloose = TMath::Pi()/2.;

const float vertex_cut = 60;

float getDPHI(float phi1, float phi2);
std::vector<std::pair<struct jet, struct jet>>  match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets);

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
void makeMatchedTreesInclusiveAuAu(const int sim_version = 20, const int cone_size_int = 3)
{

  std::string infile = Form("/sphenix/tg/tg01/jets/dlis/sim/hijing/v14/TREE_JET_SIM_v14_%d_new_ProdA_2024-00000030.root", sim_version);

  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "TREE_MATCH_r0" + std::to_string(cone_size_int);
  std::string oldname = "TREE_JET_SIM";
  mycopy.replace(mycopy.find(oldname), oldname.length(), newname);

  std::cout << mycopy << std::endl;


  cone_size = (float) cone_size_int;
  TFile *f = new TFile(infile.c_str(), "r");
  if (!f)
    {
      std::cout << " no file" << std::endl;
      return;
    }
  TTree *t = (TTree*) f->Get("ttree");
  if (!t)
    {
      std::cout << " no tree" << std::endl;
      return;
    }
  
  std::vector<float> *truth_jet_pt_ref = 0;
  std::vector<float> *truth_jet_pt = 0;
  std::vector<float> *truth_jet_eta = 0;
  std::vector<float> *truth_jet_phi = 0;

  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;
  float mbd_vertex_z;
  //float truth_vertex_z;
  int centrality = 0;
  int minbias = 0;
  int flow_fail = 0;
  float psi = 0;

  t->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_ref);
  t->SetBranchAddress(Form("truth_jet_pt_%d", cone_size_int), &truth_jet_pt);
  t->SetBranchAddress(Form("truth_jet_eta_%d", cone_size_int), &truth_jet_eta);
  t->SetBranchAddress(Form("truth_jet_phi_%d", cone_size_int), &truth_jet_phi);
  t->SetBranchAddress(Form("jet_pt_%d_sub", cone_size_int), &reco_jet_pt);
  t->SetBranchAddress(Form("jet_e_%d_sub", cone_size_int), &reco_jet_e);
  t->SetBranchAddress(Form("jet_e_unsub_%d_sub", cone_size_int), &reco_jet_e_unsub);
  t->SetBranchAddress(Form("jet_eta_%d_sub", cone_size_int), &reco_jet_eta);
  t->SetBranchAddress(Form("jet_phi_%d_sub", cone_size_int), &reco_jet_phi);
  t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  t->SetBranchAddress("centrality", &centrality);
  t->SetBranchAddress("minbias", &minbias);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;

  int entries = t->GetEntries();
  if (Debug) entries = 100000;
  std::cout << entries << std::endl;
  TFile *fout = new TFile(mycopy.c_str(), "recreate");
  TNtuple *tn_stats = new TNtuple("tn_stats","run stats","nevents");
  tn_stats->Fill(entries);
  
  TNtuple *tn_match = new TNtuple("tn_match","matched truth and reco","maxpttruth:pt1_truth:pt2_truth:dphi_truth:pt1_reco:pt2_reco:dphi_reco:matched:nrecojets:centrality:mbd_vertex");
  TEfficiency *he_pt_truth[10];
  TEfficiency *he_pt_dijet[10];
  TH2D *h2_ptt_ptrptt[10];
  for (int i = 0; i < 10; i++)
    {
      h2_ptt_ptrptt[i] = new TH2D(Form("h2_ptt_ptrptt_%d", i),"", 25, 0, 50, 120, 0, 1.2);
      he_pt_truth[i] = new TEfficiency(Form("he_pt_truth_%d", i),"", 25, 0, 50);
      he_pt_dijet[i] = new TEfficiency(Form("he_pt_dijet_%d", i),"", 25, 0, 50);
    }

  //10 GeV sample
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      if (fabs(mbd_vertex_z) > vertex_cut) continue;
      if (!minbias) continue;
      if (!Debug) std::cout << "Event: " << i  << "\r"<<std::flush;
      int cent_bin = centrality/10;
      int ntruthjets = truth_jet_pt->size();      
      int nrecojets = reco_jet_pt->size();
      if (Debug) std::cout << __LINE__ <<  " " << nrecojets << std::endl;
      if (ntruthjets == 0) continue;
      if (Debug) std::cout << __LINE__ << std::endl;
      float max_pt_dijet_r04 = *std::max_element(truth_jet_pt_ref->begin(), truth_jet_pt_ref->end());;

      // this is the reco index for the leading and subleading jet
      std::pair<int,int> dijet_ids = std::make_pair(-1,-1);

      matched_dijets.clear();
      mytruthjets.clear();
      myrecojets.clear();

      bool found_truth_dijet = false;
      bool found_reco_dijet = false;
      
      for (int j = 0; j < ntruthjets;j++)
	{
	  struct jet tempjet;
	  if (truth_jet_pt->at(j) < truth_cut) continue;
	  //	  if (fabs(truth_jet_eta->at(j)) > etacut) continue;	  

	  tempjet.pt = truth_jet_pt->at(j);
	  tempjet.eta = truth_jet_eta->at(j);
	  tempjet.phi = truth_jet_phi->at(j);
	  tempjet.id = j;
	  mytruthjets.push_back(tempjet);	  
	  
	}
      if (Debug) std::cout << __LINE__ << std::endl;
      int nntruthjets = mytruthjets.size();
      double dphitr = 0;      
      if (mytruthjets.size() >= 2) 
	{
	  std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  auto leading_iter = mytruthjets.begin();
	  auto subleading_iter = mytruthjets.begin() + 1;
	  if (subleading_iter != mytruthjets.end())
	    {
	      dphitr = leading_iter->phi - subleading_iter->phi;
	      if (dphitr > TMath::Pi())
		{
		  dphitr -= 2*TMath::Pi();
		}
	      if (dphitr < -1*TMath::Pi())
		{
		  dphitr += 2*TMath::Pi();
		}
	      if (fabs(dphitr) >= dphicut) found_truth_dijet = true;
	    }

	}
            
      // Getting the reco jets
      for (int j = 0; j < nrecojets;j++)
	{

	  if (reco_jet_pt->at(j) < reco_cut) continue;
	  if (reco_jet_e->at(j) < 0) continue;
	  if (reco_jet_e_unsub->at(j) < 0) continue;
	  //	  if (fabs(reco_jet_eta->at(j)) > etacut) continue;

	  struct jet tempjet;
	  tempjet.pt = reco_jet_pt->at(j);
	  tempjet.eta = reco_jet_eta->at(j);
	  tempjet.phi = reco_jet_phi->at(j);
	  tempjet.id = j;
	  myrecojets.push_back(tempjet);

	}

      if (Debug) std::cout << __LINE__ << std::endl;
      int nnrecojets = myrecojets.size();
      double dphir = 0;

      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      matched_dijets = match_dijets(myrecojets, mytruthjets);
      if (Debug) std::cout << __LINE__ << std::endl;
      for (auto tjet : mytruthjets)
	{
	  he_pt_truth[cent_bin]->Fill(std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) {return tjet.id == a.first.id; }) != matched_dijets.end(), tjet.pt);	  
	}
      for (auto mjets : matched_dijets)
	{
	  float ptt = mjets.first.pt;
	  float ptr = mjets.second.pt;
	  h2_ptt_ptrptt[cent_bin]->Fill(ptt, ptr/ptt);
	}

      // if top reco jets match a dijet, fill, if not, a fake
      if (matched_dijets.size() > 0)
	{

	  std::sort(matched_dijets.begin(), matched_dijets.end(), [] (auto a, auto b) { return a.first.pt > b.first.pt; });
	  if (Debug)
	    {
	      std::cout << "Matched Jets" << std::endl;
	      for (auto mjets : matched_dijets)
		{
		  std::cout << "truth: " << std::endl;
		  mjets.first.print();
		  std::cout << "reco: " << std::endl;
		  mjets.second.print();
		}
	    }
	}
      
      if (Debug) std::cout << __LINE__ << std::endl;
      // for (auto tjet : mytruthjets)
      // 	{
      // 	  h_miss_rate_v_truth->Fill(std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id; }) != matched_dijets.end(), tjet.pt);
      // 	}
      // for (auto rjet : mytruthjets)
      // 	{
      // 	  h_fake_rate_v_reco->Fill(std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == rjet.id; }) != matched_dijets.end(), rjet.pt);
      // 	}

      // check if leading and subleading are matched
      if (mytruthjets.size() >= 2)
	{
	  auto leading_truth_iter = mytruthjets.begin();
	  auto subleading_truth_iter = mytruthjets.begin() + 1;
	  auto found_leading = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == leading_truth_iter->id; });
	  auto found_subleading = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == subleading_truth_iter->id; });
	  if (Debug) std::cout << __LINE__ << std::endl;
	  if (found_leading != matched_dijets.end() && found_subleading != matched_dijets.end())
	    {
	      he_pt_dijet[cent_bin]->Fill(1, leading_truth_iter->pt);
	      if (Debug)
		{
		  std::cout << "Found Lead: "<<std::endl;
		  found_leading->first.print();
		  std::cout << "      Reco: "<<std::endl;
		  found_leading->second.print();
		  std::cout << "      Sublead: "<<std::endl;
		  found_subleading->first.print();
		  std::cout << "      Reco: "<<std::endl;
		  found_subleading->second.print();
		}
	      float pt1truth = found_leading->first.pt;
	      float eta1truth = found_leading->first.eta;
	      float phi1truth = found_leading->first.phi;
	      
	      float pt2truth = found_subleading->first.pt;
	      float eta2truth = found_subleading->first.eta;
	      float phi2truth = found_subleading->first.phi;

	      float pt1reco = found_leading->second.pt;
	      float eta1reco = found_leading->second.eta;
	      float phi1reco = found_leading->second.phi;
	      
	      float pt2reco = found_subleading->second.pt;
	      float eta2reco = found_subleading->second.eta;
	      float phi2reco = found_subleading->second.phi;

	      float dphi_reco_match = fabs(getDPHI(phi1reco, phi2reco));
	      float dphi_truth_match = fabs(getDPHI(phi1truth, phi2truth));

	      tn_match->Fill(max_pt_dijet_r04, pt1truth, pt2truth, dphi_truth_match, pt1reco, pt2reco, dphi_reco_match, 1, myrecojets.size(), centrality, mbd_vertex_z);	      
	    }
	  else
	    {
	      he_pt_dijet[cent_bin]->Fill(0, leading_truth_iter->pt);
	      float pt1truth = mytruthjets.at(0).pt;
	      float eta1truth = mytruthjets.at(0).eta;
	      float phi1truth = mytruthjets.at(0).phi;
	      
	      float pt2truth = mytruthjets.at(1).pt;
	      float eta2truth = mytruthjets.at(1).eta;
	      float phi2truth = mytruthjets.at(1).phi;

	      float dphi_truth_match = fabs(getDPHI(phi1truth, phi2truth));

	      tn_match->Fill(max_pt_dijet_r04, pt1truth, pt2truth, dphi_truth_match, 0, 0, 0, 0, 0, centrality, mbd_vertex_z);	      
	    }
	} 
    }
  for (int i = 0; i < 10; i++)
    {
      he_pt_truth[i]->Write();
      he_pt_dijet[i]->Write();
      h2_ptt_ptrptt[i]->Write();
    }
  tn_match->Write();
  tn_stats->Write();
  fout->Close();

}
std::vector<std::pair<struct jet, struct jet>>  match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto tjet : mytruthjets)
    {
      for (auto jet : myrecojets)
	{
	  float dR = getDR(jet, tjet);
	      
	  if (dR < dRcut*cone_size)
	    {
	      if (std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;}) != matched_dijets.end()) continue;

	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;

	      matched_dijets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }

  return matched_dijets;
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
