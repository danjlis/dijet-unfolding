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
const float dRcut = 0.75;
static float cone_size = 4;
const float truth_cut = 3;
const float reco_cut = 3;
const float etacut = 0.7;

const float dphicut = 1*TMath::Pi()/2.;
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
void makeMatchedTreesInclusive(const int sim_version = 10, const int cone_size_int = 4, const int is_auau = 0)
{

  int version = 6;
  if (is_auau) cersion = 7;

  std::string infile = Form("/sphenix/user/dlis/Projects/jet/simtrees/TREE_JET_v%d_%d_new_ProdA_2024-00000021.root", version, sim_version);

  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "TREE_MATCH_r0" + std::to_string(cone_size_int);
  std::string oldname = "TREE_JET";
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
  
  std::vector<float> *truth_jet_pt_4 = 0;
  //std::vector<float> *truth_jet_pt_ref = 0;
  std::vector<float> *truth_jet_eta_4 = 0;
  std::vector<float> *truth_jet_phi_4 = 0;
  std::vector<float> *reco_jet_pt_4 = 0;
  std::vector<float> *reco_jet_emcal_4 = 0;
  std::vector<float> *reco_jet_eta_4 = 0;
  std::vector<float> *reco_jet_eta_det_4 = 0;
  std::vector<float> *reco_jet_phi_4 = 0;
  float mbd_vertex_z;
  float truth_vertex_z;

  t->SetBranchAddress(Form("truth_jet_pt_%d", cone_size_int), &truth_jet_pt_4);
  t->SetBranchAddress(Form("truth_jet_eta_%d", cone_size_int), &truth_jet_eta_4);
  t->SetBranchAddress(Form("truth_jet_phi_%d", cone_size_int), &truth_jet_phi_4);
  if (!is_auau)
    {
      t->SetBranchAddress(Form("jet_pt_%d", cone_size_int), &reco_jet_pt_4);
      t->SetBranchAddress(Form("jet_emcal_%d", cone_size_int), &reco_jet_emcal_4);
      t->SetBranchAddress(Form("jet_eta_%d", cone_size_int), &reco_jet_eta_4);
      t->SetBranchAddress(Form("jet_eta_det_%d", cone_size_int), &reco_jet_eta_det_4);
      t->SetBranchAddress(Form("jet_phi_%d", cone_size_int), &reco_jet_phi_4);
    }
  else
    {
      t->SetBranchAddress(Form("jet_pt_%d_sub", cone_size_int), &reco_jet_pt_4);
      t->SetBranchAddress(Form("jet_emcal_%d_sub", cone_size_int), &reco_jet_emcal_4);
      t->SetBranchAddress(Form("jet_eta_%d_sub", cone_size_int), &reco_jet_eta_4);
      t->SetBranchAddress(Form("jet_eta_det_%d_sub", cone_size_int), &reco_jet_eta_det_4);
      t->SetBranchAddress(Form("jet_phi_%d_sub", cone_size_int), &reco_jet_phi_4);
    }
  t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  t->SetBranchAddress("truth_vertex_z", &truth_vertex_z);

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;

  int entries = t->GetEntries();
  if (Debug) entries = 100000;
  std::cout << entries << std::endl;
  TFile *fout = new TFile(mycopy.c_str(), "recreate");
  TNtuple *tn_match = new TNtuple("tn_match","matched truth and reco","maxpttruth:pt1_truth:pt2_truth:dphi_truth:pt1_reco:pt2_reco:em1_reco:em2_reco:dphi_reco:matched:thirdjet_pt:thirdjet_pt_truth:nrecojets:mbd_vertex");

  TH1D *h_dR[2];
  TH2D *h_ptrptt_v_ptt[2];
  TH2D *h_ptrptt_deta_fake[2];
  TH2D *h_ptrptt_dphi_fake[2];
  TH2D *h_ptrptt_deta_miss[2];
  TH2D *h_ptrptt_dphi_miss[2];
  TH2D *h_dvtx_deta_fake[2];
  TH2D *h_dvtx_dphi_fake[2];
  TH2D *h_dvtx_deta_miss[2];
  TH2D *h_dvtx_dphi_miss[2];
  TH1D *h_ptrptt[2];
  TEfficiency *he_truth_match[2];
  TEfficiency *he_reco_fake[2];

  TH1D *h_deta_miss[2];
  TH1D *h_deta_fake[2];
  TH1D *h_dphi_miss[2];
  TH1D *h_dphi_fake[2];

  TH1D *h_dR_miss[2];
  TH1D *h_dR_fake[2];
  TH1D *h_ptrptt_miss[2];
  TH1D *h_ptrptt_fake[2];

  for (int i = 0; i < 2; i++)
    {
      h_dR[i] = new TH1D(Form("h_dR_%d", i), ";dR;", 100, 0.0, 1.0);

      he_truth_match[i] = new TEfficiency(Form("he_truth_math_%d", i), "; p_{T, truth};", 100, 0, 100);
      he_reco_fake[i] = new TEfficiency(Form("he_reco_fake_%d", i), "; p_{T, truth};", 100, 0, 100);

      h_ptrptt_v_ptt[i] = new TH2D(Form("h_ptrptt_v_ptt_%d", i), ";p_{t, truth}; p_{T, reco}/p_{T, truth};", 100, 0.0, 100.0, 100, 0, 2.);
      h_ptrptt[i] = new TH1D(Form("h_ptrptt_%d", i), ";p_{T, reco}/p_{T, truth};", 100, 0.0, 2.0);
      h_dR_miss[i] = new TH1D(Form("h_dR_miss_%d", i), ";dR_miss;", 100, -1.0, 1.0);
      h_dR_fake[i] = new TH1D(Form("h_dR_fake_%d", i), ";dR_fake;", 100, -1.0, 1.0);
      h_deta_miss[i] = new TH1D(Form("h_deta_miss_%d", i), ";deta_miss;", 100, -1.0, 1.0);
      h_deta_fake[i] = new TH1D(Form("h_deta_fake_%d", i), ";deta_fake;", 100, -1.0, 1.0);
      h_dphi_miss[i] = new TH1D(Form("h_dphi_miss_%d", i), ";dphi_miss;", 100, -1.0, 1.0);    
      h_dphi_fake[i] = new TH1D(Form("h_dphi_fake_%d", i), ";dphi_fake;", 100, -1.0, 1.0);    
      h_ptrptt_deta_miss[i] = new TH2D(Form("h_ptrptt_deta_miss_%d", i), ";deta_miss;", 100, 0, 2, 100, -1.0, 1.0);
      h_ptrptt_deta_fake[i] = new TH2D(Form("h_ptrptt_deta_fake_%d", i), ";deta_fake;", 100, 0, 2, 100, -1.0, 1.0);
      h_ptrptt_dphi_miss[i] = new TH2D(Form("h_ptrptt_dphi_miss_%d", i), ";dphi_miss;", 100, 0, 2,  100, -1.0, 1.0);    
      h_ptrptt_dphi_fake[i] = new TH2D(Form("h_ptrptt_dphi_fake_%d", i), ";dphi_fake;", 100, 0, 2,  100, -1.0, 1.0);    

      h_dvtx_deta_miss[i] = new TH2D(Form("h_dvtx_deta_miss_%d", i), ";deta_miss;", 100, -100, 100, 100, -1.0, 1.0);
      h_dvtx_deta_fake[i] = new TH2D(Form("h_dvtx_deta_fake_%d", i), ";deta_fake;", 100, -100, 100, 100, -1.0, 1.0);
      h_dvtx_dphi_miss[i] = new TH2D(Form("h_dvtx_dphi_miss_%d", i), ";dphi_miss;", 100, -100, 100,  100, -1.0, 1.0);    
      h_dvtx_dphi_fake[i] = new TH2D(Form("h_dvtx_dphi_fake_%d", i), ";dphi_fake;", 100, -100, 100,  100, -1.0, 1.0);    

      h_ptrptt_miss[i] = new TH1D(Form("h_ptrptt_miss_%d", i), ";ptrptt_miss;", 100, 0.0, 2.0);    
      h_ptrptt_fake[i] = new TH1D(Form("h_ptrptt_fake_%d", i), ";ptrptt_fake;", 100, 0.0, 2.0);      
   }

  TEfficiency *h_fake_rate_v_reco = new TEfficiency("h_fake_rate_v_reco", ";Reco p_{T} [GeV]; Fake Rate;", 100, 0, 100);
  TEfficiency *h_miss_rate_v_truth = new TEfficiency("h_miss_rate_v_truth", ";Truth p_{T} [GeV]; Miss Rate;", 100, 0, 100);
  TEfficiency *h_miss_rate_v_vertex = new TEfficiency("h_miss_rate_v_vertex", ";Vertex p_{T} [GeV]; Miss Rate;", 100, -100, 100);

  TEfficiency *h2_fake_rate_v_pos = new TEfficiency("h2_fake_rate_v_pos", ";#eta; #phi; Fake Rate", 24, -1.1, 1.1, 64, -1*TMath::Pi(), TMath::Pi());
  TEfficiency *h2_miss_rate_v_pos = new TEfficiency("h2_miss_rate_v_pos", ";#eta; #phi; Fake Rate;", 24, -1.1, 1.1, 64, -1*TMath::Pi(), TMath::Pi());

  //10 GeV sample
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      if (Nuclear) std::cout << __LINE__ << std::endl;
      if (fabs(mbd_vertex_z) > vertex_cut) continue;
      if (!Debug) std::cout << "Event: " << i  << "\r"<<std::flush;
      if (Nuclear) std::cout << __LINE__ << std::endl;
      int ntruthjets = truth_jet_pt_4->size();      
      int nrecojets = reco_jet_pt_4->size();
      if (Nuclear) std::cout << __LINE__ << std::endl;
      if (ntruthjets == 0) continue;

      float max_pt_dijet_r04 = *std::max_element(truth_jet_pt_4->begin(), truth_jet_pt_4->end());;

      if (Nuclear) std::cout << __LINE__ << std::endl;
      // this is the reco index for the leading and subleading jet
      std::pair<int,int> dijet_ids = std::make_pair(-1,-1);
      matched_dijets.clear();
      mytruthjets.clear();
      myrecojets.clear();

      bool found_truth_dijet = false;
      bool found_reco_dijet = false;
      if (Nuclear) std::cout << __LINE__ << std::endl;
      float third_jet_pt_truth = 0;
      for (int j = 0; j < ntruthjets;j++)
	{
	  struct jet tempjet;
	  if (truth_jet_pt_4->at(j) < truth_cut) continue;
	  if (fabs(truth_jet_eta_4->at(j)) > etacut) continue;	  

	  tempjet.pt = truth_jet_pt_4->at(j);
	  tempjet.eta = truth_jet_eta_4->at(j);
	  tempjet.phi = truth_jet_phi_4->at(j);
	  tempjet.id = j;
	  mytruthjets.push_back(tempjet);	  
	  
	}


      int nntruthjets = mytruthjets.size();
      double dphitr = 0;      
      if (mytruthjets.size() >= 2) 
	{

	  std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  auto leading_iter = mytruthjets.begin();
	  auto subleading_iter = mytruthjets.begin() + 1;
	  if (mytruthjets.size() >= 3)
	    {
	      third_jet_pt_truth = mytruthjets.at(2).pt;
	    }

	  // Check if the truth jet satisfies the cuts
	  // check dphi
	  // auto subleading_iter = std::find_if(leading_iter + 1, mytruthjets.end(), [=] (auto a)
	  // 			    {
	  // 			      double temp_dphitr = fabs(leading_iter->phi - a.phi);
	  // 			      if (temp_dphitr > TMath::Pi())
	  // 				{
	  // 				  temp_dphitr = 2*TMath::Pi() - temp_dphitr;
	  // 				}
	  // 			      if (temp_dphitr > dphicut) return true;
	  // 			      return false;
	  // 			    });

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

	  if (reco_jet_pt_4->at(j) < reco_cut) continue;
	  if (fabs(reco_jet_eta_4->at(j)) < etacut) continue;

	  struct jet tempjet;
	  tempjet.pt = reco_jet_pt_4->at(j);
	  tempjet.eta = reco_jet_eta_4->at(j);
	  tempjet.emcal = reco_jet_emcal_4->at(j);
	  tempjet.phi = reco_jet_phi_4->at(j);
	  tempjet.id = j;
	  myrecojets.push_back(tempjet);

	}
      float third_jet_pt = 0;

      int nnrecojets = myrecojets.size();
      double dphir = 0;
      if (Nuclear) std::cout << __LINE__ << std::endl;
      if (myrecojets.size() >= 2) 
	{

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  auto leading_iter = myrecojets.begin();
	  // Check if the reco jet satisfies the cuts
	  // check dphi
	  auto subleading_iter = myrecojets.begin() + 1;
	  if (nnrecojets >= 3)
	    {
	      third_jet_pt = myrecojets.at(2).pt;
	    }
	  // std::find_if(leading_iter + 1, myrecojets.end(), [=] (auto a)
	  // 				    {
	  // 				      double temp_dphir = fabs(leading_iter->phi - a.phi);
	  // 				      if (temp_dphir > TMath::Pi())
	  // 					{
	  // 					  temp_dphir = 2*TMath::Pi() - temp_dphir;
	  // 					}
	  // 				      if (temp_dphir > dphicut) return true;
	  // 				      return false;
	  // 				    });

	  if (subleading_iter != mytruthjets.end())
	    {
	      dphir = leading_iter->phi - subleading_iter->phi;
	      if (dphir > TMath::Pi())
		{
		  dphir -= 2*TMath::Pi();
		}
	      if (dphir < -1*TMath::Pi())
		{
		  dphir += 2*TMath::Pi();
		}
	      if (fabs(dphir) >= dphicut) found_reco_dijet = true;

	      myrecojets = {*leading_iter, *subleading_iter};
	    }

	}

      // if (myrecojets.size() >= 2) 
      // 	{

      // 	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      // 	  dphir = fabs(myrecojets.at(0).phi - myrecojets.at(1).phi);
      // 	  if (dphir > TMath::Pi())
      // 	    {
      // 	      dphir = 2*TMath::Pi() - dphir;
      // 	    }
      // 	  if (dphir > dphicut) found_reco_dijet = true;


      // 	}
      if (Nuclear) std::cout << __LINE__ << std::endl;
      // if (found_reco_dijet && !found_truth_dijet) 
      // 	{
      // 	  //if (Debug) std::cout << " FAKE " <<std::endl;
      // 	  h_fake_rate_v_reco->Fill(1, myrecojets.at(0).pt);
      // 	  h2_fake_rate_v_pos->Fill(1, tjet.eta, tjet.phi);  
      // 	  //h2_fake_rate_v_pos->Fill(1, myrecojets.at(0).eta, myrecojets.at(0).phi);
      // 	  tn_match->Fill(0,0,0, myrecojets.at(0).pt, myrecojets.at(1).pt, dphir, 0);
      // 	  continue;
      // 	}
      // else if (!found_reco_dijet && found_truth_dijet) 
      // 	{
      // 	  //if (Debug) std::cout << " MISS" <<std::endl;
      // 	  h_miss_rate_v_truth->Fill(0, mytruthjets.at(0).pt);
      // 	  h2_miss_rate_v_pos->Fill(1, tjet.eta, tjet.phi);  
      // 	  //	  h2_miss_rate_v_pos->Fill(0, mytruthjets.at(0).eta, mytruthjets.at(0).phi);
      // 	  tn_match->Fill(mytruthjets.at(0).pt, mytruthjets.at(1).pt, dphitr, 0,0,0, 0);
      // 	  continue;
      // 	}
      if (Nuclear) std::cout << __LINE__ << std::endl;
      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });
      if (Nuclear) std::cout << __LINE__ << std::endl;
      matched_dijets = match_dijets(myrecojets, mytruthjets);
      if (Nuclear) std::cout << __LINE__ << std::endl;
      // if top reco jets match a dijet, fill, if not, a fake
      if (matched_dijets.size() > 0)
	{
	  std::sort(matched_dijets.begin(), matched_dijets.end(), [] (auto a, auto b) { return a.second.pt > b.second.pt; });
	}
      if (Nuclear) std::cout << __LINE__ << std::endl;
      for (int itruth = 0; itruth < mytruthjets.size(); itruth++)
	{
	  if (Nuclear) std::cout << __LINE__ << std::endl;
	  auto tjet = mytruthjets.at(itruth);
	  if (fabs(tjet.eta) > etacut) continue;
	  auto titer = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id; });
	  if (titer == matched_dijets.end())
	    {
	      struct jet minjet;
	      for (auto jet : myrecojets)
		{
		  float dR = getDR(jet, tjet);
		      
		  if (dR < minjet.dR)
		    {
		      minjet = jet;
		      minjet.dR = dR;
		    }
		}
	      if (minjet.dR < dRcut*cone_size && minjet.matched)
		{
		  continue;
		}
	      
	      float deta = tjet.eta - minjet.eta;		  
	      h_deta_miss[0]->Fill(deta);
	      float dphitm = tjet.phi - minjet.phi;
	      if (dphitm > TMath::Pi()) dphitm = -2*TMath::Pi() + dphitm ;
	      if (dphitm < -1*TMath::Pi()) dphitm = 2*TMath::Pi() + dphitm ;
	      h_dphi_miss[0]->Fill(dphitm);

	      h_dR_miss[0]->Fill(sqrt(dphitm*dphitm + deta*deta));

	      h_ptrptt_miss[0]->Fill(minjet.pt/tjet.pt);
	      h_ptrptt_dphi_miss[0]->Fill(minjet.pt/tjet.pt, dphitm);
	      h_ptrptt_deta_miss[0]->Fill(minjet.pt/tjet.pt, deta);
	      if (Nuclear) std::cout << __LINE__ << std::endl;
	      if (minjet.pt/tjet.pt > 0.2)
		{
		  h_dvtx_dphi_miss[0]->Fill(truth_vertex_z - mbd_vertex_z, dphitm);
		  h_dvtx_deta_miss[0]->Fill(truth_vertex_z - mbd_vertex_z, deta);
		}
	      h_miss_rate_v_truth->Fill(0, tjet.pt);
	    }
	  else
	    {	
	      h_miss_rate_v_truth->Fill(1, tjet.pt);
	    }
	}

      for (int ireco = 0; ireco < myrecojets.size(); ireco++)
	{
	  if (Nuclear) std::cout << __LINE__ << std::endl;
	  auto tjet = myrecojets.at(ireco);
	  if (fabs(tjet.eta) > etacut) continue;
	  auto titer = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == tjet.id; });
	  if (titer == matched_dijets.end())
	    {
	      struct jet minjet;
	      for (auto jet : mytruthjets)
		{
		  float dR = getDR(jet, tjet);
		      
		  if (dR < minjet.dR)
		    {
		      minjet = jet;
		      minjet.dR = dR;
		    }
		}
	      if (minjet.matched && minjet.dR < dRcut*cone_size)
		{
		  continue;
		}
	      h_fake_rate_v_reco->Fill(0, tjet.pt);
	    }
	  else
	    {	
	      h_fake_rate_v_reco->Fill(1, tjet.pt);
	    }
	}

      if (Nuclear) std::cout << __LINE__ << std::endl; 
      float pt1truth = 0;
      float pt2truth = 0;
      float eta1truth = 0;
      float eta2truth = 0;
      float phi1truth = 0;
      float phi2truth = 0;

      if (mytruthjets.size() >= 1)
	{
	  pt1truth = mytruthjets.at(0).pt;
	  eta1truth = mytruthjets.at(0).eta;
	  phi1truth = mytruthjets.at(0).phi;
	}
      if (mytruthjets.size() >= 2)
	{
	  pt2truth = mytruthjets.at(1).pt;
	  eta2truth = mytruthjets.at(1).eta;
	  phi2truth = mytruthjets.at(1).phi;
	}

      float pt1reco = 0;
      float pt2reco = 0;
      float eta1reco = 0;
      float eta2reco = 0;
      float phi1reco = 0;
      float phi2reco = 0;

      if (myrecojets.size() >= 1)
	{
	  pt1reco = myrecojets.at(0).pt;
	  eta1reco = myrecojets.at(0).eta;
	  phi1reco = myrecojets.at(0).phi;
	}
      if (myrecojets.size() >= 2)
	{
	  pt2reco = myrecojets.at(1).pt;
	  eta2reco = myrecojets.at(1).eta;
	  phi2reco = myrecojets.at(1).phi;
	}
      // if top reco jets do not make a dijet...

      // std::cout << "Truth Good? " << (found_truth_dijet ? "Yes" : "No") << "  Reco Good? " << (found_reco_dijet?"Yes":"No") << std::endl;
      // // check if there is one an thats a miss
      // std::cout << "Matched Jets: " << std::endl;
      // for (auto dijet : matched_dijets)
      // 	{
      // 	  std::cout << dijet.first.pt << "\t"<< dijet.second.pt << std::endl;
      // 	  std::cout << dijet.first.eta << "\t"<< dijet.second.eta << std::endl;
      // 	  std::cout << dijet.first.phi << "\t"<< dijet.second.phi << std::endl;
      // 	} 
      // {
      // 	int im = 0;
      // 	int it = 0;
      // 	for (auto match : mytruthjets)
      // 	  {
      // 	    std::cout << "Truth" << it++ << std::endl;
      // 	    std::cout << "    " << match.pt << std::endl;
      // 	    std::cout << "    " << match.eta << std::endl;
      // 	    std::cout << "    " << match.phi << std::endl;
      // 	  }
      // 	int ir = 0;
      // 	for (auto match : myrecojets)
      // 	  {
      // 	    std::cout << "Reco " << ir++ << std::endl;
      // 	    std::cout << "    " << match.pt << std::endl;
      // 	    std::cout << "    " << match.eta << std::endl;
      // 	    std::cout << "    " << match.phi << std::endl;	      
      // 	  }
      // }
      if (Nuclear) std::cout << __LINE__ << std::endl;
      for (auto dijet : matched_dijets)
	{
	  h_ptrptt[0]->Fill(dijet.second.pt/dijet.first.pt);
	}
      if (Nuclear) std::cout << __LINE__ << std::endl;
      if (matched_dijets.size() >= 2)
	{

	  float dphi_reco_match = getDPHI(matched_dijets.at(0).second.phi, matched_dijets.at(1).second.phi);
	  float dphi_truth_match = getDPHI(matched_dijets.at(0).first.phi, matched_dijets.at(1).first.phi);

	  tn_match->Fill(max_pt_dijet_r04, matched_dijets.at(0).first.pt, matched_dijets.at(1).first.pt, dphi_truth_match, matched_dijets.at(0).second.pt, matched_dijets.at(1).second.pt, matched_dijets.at(0).second.emcal, matched_dijets.at(1).second.emcal, dphi_reco_match, 1, third_jet_pt_truth, third_jet_pt, nnrecojets, mbd_vertex_z);
	  continue;
	}

      if (found_reco_dijet && !found_truth_dijet)
	{
	  tn_match->Fill(max_pt_dijet_r04, pt1truth, pt2truth, dphitr, myrecojets.at(0).pt, myrecojets.at(1).pt,myrecojets.at(0).emcal, myrecojets.at(1).emcal, dphir, 0, third_jet_pt_truth, third_jet_pt, nnrecojets, mbd_vertex_z);
	  continue;
	}
      else if (found_truth_dijet && !found_reco_dijet)
	{
	  tn_match->Fill(max_pt_dijet_r04, mytruthjets.at(0).pt, mytruthjets.at(1).pt, dphitr, pt1reco, pt2reco, 0, 0, dphir, 0, third_jet_pt_truth, third_jet_pt, nnrecojets, mbd_vertex_z);
	  continue;
	}
      else if (found_truth_dijet && found_reco_dijet)
	{
	  tn_match->Fill(max_pt_dijet_r04, mytruthjets.at(0).pt, mytruthjets.at(1).pt, dphitr, myrecojets.at(0).pt, myrecojets.at(1).pt,myrecojets.at(0).emcal, myrecojets.at(1).emcal, dphir, 0, third_jet_pt_truth, third_jet_pt, nnrecojets, mbd_vertex_z);
	  continue;
	}
      else
	{
	  continue;
	}
    }

  tn_match->Write();
  h_fake_rate_v_reco->Write();
  h_miss_rate_v_truth->Write();
  h_miss_rate_v_vertex->Write();
  h2_fake_rate_v_pos->Write();
  h2_miss_rate_v_pos->Write();
  //h_leading_dR->Write();
  //h_subleading_dR->Write();

  for (int i = 0; i <2; i++)
    {
      h_dR[i]->Write();
      h_ptrptt_v_ptt[i]->Write();
      h_ptrptt[i]->Write();
      he_truth_match[i]->Write();
      he_reco_fake[i]->Write();

      h_deta_miss[i]->Write();
      h_deta_fake[i]->Write();
      h_dphi_miss[i]->Write();
      h_dphi_fake[i]->Write();

      h_dR_miss[i]->Write();
      h_dR_fake[i]->Write();
      h_ptrptt_miss[i]->Write();
      h_ptrptt_fake[i]->Write();
      h_ptrptt_deta_miss[i]->Write();
      h_ptrptt_deta_fake[i]->Write();
      h_ptrptt_dphi_miss[i]->Write();
      h_ptrptt_dphi_fake[i]->Write();

      h_dvtx_deta_miss[i]->Write();
      h_dvtx_deta_fake[i]->Write();
      h_dvtx_dphi_miss[i]->Write();
      h_dvtx_dphi_fake[i]->Write();

    }
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
	      
	  if (dR < dRcut*cone_size && dR < jet.dR)
	    {
		  
	      auto riter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id; });
	      auto titer = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id; });
	      // if both aren't matched match them
	      if (titer == matched_dijets.end() && riter == matched_dijets.end())
		{
		  tjet.matched = 1;
		  tjet.dR = dR;
		  jet.matched = 1;
		  jet.dR = dR;

		  matched_dijets.push_back(std::make_pair(tjet, jet));
		}
	      else if (titer == matched_dijets.end())
		{
		  if (tjet.pt < riter->first.pt) continue;

		  jet.matched = 1;
		  jet.dR = dR;
		  tjet.matched = 1;
		  tjet.dR = dR;

		  matched_dijets.erase(riter);
		  matched_dijets.push_back(std::make_pair(tjet, jet));
		}
	      else if (riter == matched_dijets.end())
		{
		  if (titer->second.pt < jet.pt)
		    {
		      jet.matched = 1;
		      jet.dR = dR;
		      tjet.matched = 1;
		      tjet.dR = dR;

		      matched_dijets.erase(titer);
		      matched_dijets.push_back(std::make_pair(tjet, jet));
		    }
		}
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
