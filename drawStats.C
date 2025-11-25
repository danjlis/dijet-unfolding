#include "dlUtility.h"
#include "read_binning.h"
float vertex_cut = 60;

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

const float truth_cut = 3;
const float reco_cut = 3;
const float etacut = 0.8;

void drawStats(const std::string configfile = "binning_AA.config")
{
  const int cone_size_int = 3;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());
    
  std::string data_file = rb.get_tntuple_location() + "TNTUPLE_DIJET_r03_v10_1_492_2024p020_v007_gl10-all.root";
  std::string j10_file[3];
  int sim_version = 10;
  std::string infile = Form("%s/TREE_JET_SIM_v14_%d_new_ProdA_2024-00000030.root", rb.get_sim_location().c_str(), sim_version);
  j10_file[0] = infile;
  sim_version = 20;
  infile = Form("%s/TREE_JET_SIM_v14_%d_new_ProdA_2024-00000030.root", rb.get_sim_location().c_str(), sim_version);
  j10_file[1] =  infile;
  sim_version = 30;
  infile = Form("%s/TREE_JET_SIM_v14_%d_new_ProdA_2024-00000030.root", rb.get_sim_location().c_str(), sim_version);
  j10_file[2] =  infile;

  float cs_10 = (2.889e-6);
  float cs_20 = 5.4067742e-8;
  float cs_30 = (2.505e-9);
  
  float scale_factor[3];
  scale_factor[0] = cs_10/cs_30;
  scale_factor[1] = cs_20/cs_30; 
  scale_factor[2] = 1;

  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = 3.*TMath::Pi()/4.;//rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;
  
  TFile *find = new TFile(data_file.c_str(), "r");
  TNtuple *tnd  = (TNtuple*) find->Get("tn_dijet");;
  float pt1_data = 0;
  float pt2_data = 0;
  float dphi_data = 0;
  if (!tnd)
    {
      std::cout << " no data "<< std::endl;
    }
  tnd->SetBranchAddress("pt1_reco", &pt1_data);
  tnd->SetBranchAddress("pt2_reco", &pt2_data);
  tnd->SetBranchAddress("dphi_reco", &dphi_data);

  const int nbinsdphi  = 16;
  float min_dphi = 3*TMath::Pi()/4.;
  float max_dphi = TMath::Pi();
  float stepdphi = (max_dphi - min_dphi)/(float)nbinsdphi;

  float idphi_bins[nbinsdphi+1];
  for (int i = 0; i < nbinsdphi+1; i++)
    {
      idphi_bins[i] = min_dphi + i*stepdphi;
    }
  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

  float truth_leading_cut = 15;//rb.get_truth_leading_cut();
  float truth_subleading_cut = 8;//rb.get_truth_subleading_cut();

  float reco_leading_cut = 20;//rb.get_reco_leading_cut();
  float reco_subleading_cut = 8;//rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();
  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int binranges[10] = {0};
  int binrangesmin[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      binranges[ir] = rb.get_measure_region(ir);
      binrangesmin[ir] = rb.get_subleading_measure_region(ir);
      std::cout << binranges[ir] << " -- " <<  binrangesmin[ir] << std::endl;
    }


  TH1D *h_data_leading = new TH1D("h_data_leading",";#Delta#phi;1/N", 42, 8, 50);
  TH1D *h_truth_leading = new TH1D("h_truth_leading",";#Delta#phi;1/N", 42, 8, 50);
  TH1D *h_reco_leading = new TH1D("h_reco_leading",";#Delta#phi;1/N", 42, 8, 50);

  TH1D *h_data_subleading = new TH1D("h_data_subleading",";#Delta#phi;1/N", 42, 8, 50);
  TH1D *h_truth_subleading = new TH1D("h_truth_subleading",";#Delta#phi;1/N", 42, 8, 50);
  TH1D *h_reco_subleading = new TH1D("h_reco_subleading",";#Delta#phi;1/N", 42, 8, 50);

  for (int i = 0 ; i < 3; i++)
    {
      std::cout <<" sample : " << i << " " << j10_file[i] << std::endl;      
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
      TFile *fin = new TFile(j10_file[i].c_str(), "r");
      
      TTree *t = (TTree*) fin->Get("ttree");
      if (!t) return;
      
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

      int n_events = t->GetEntries();

      scale_factor[i] /= (float) n_events;
      std::vector<struct jet> mytruthjets;
      std::vector<struct jet> myrecojets;

      for (int ie = 0; ie < n_events; ie++)
	{
	  t->GetEntry(ie);
	  std::cout << "Event " << ie << " \r" << std::flush;
	  if (fabs(mbd_vertex_z) > vertex_cut) continue;
	  if (!minbias) continue;
	  int cent_bin = centrality/10;
	  int ntruthjets = truth_jet_pt->size();      
	  int nrecojets = reco_jet_pt->size();
	  if (ntruthjets == 0) continue;
	  float max_pt_dijet_r04 = *std::max_element(truth_jet_pt_ref->begin(), truth_jet_pt_ref->end());;
	  if (max_pt_dijet_r04 < sample_boundary[i] || max_pt_dijet_r04 >= sample_boundary[i+1]) continue;
	  // this is the reco index for the leading and subleading jet
	  std::pair<int,int> dijet_ids = std::make_pair(-1,-1);

	  mytruthjets.clear();
	  myrecojets.clear();

	  bool found_truth_dijet = false;
	  bool found_reco_dijet = false;
      
	  for (int j = 0; j < ntruthjets;j++)
	    {
	      struct jet tempjet;
	      if (truth_jet_pt->at(j) < truth_cut) continue;
	      if (fabs(truth_jet_eta->at(j)) > etacut) continue;	  

	      tempjet.pt = truth_jet_pt->at(j);
	      tempjet.eta = truth_jet_eta->at(j);
	      tempjet.phi = truth_jet_phi->at(j);
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
		  if (fabs(dphitr) >= dphicut && leading_iter->pt >= 20 && subleading_iter->pt >= 8)
		    {
		      h_truth_leading->Fill(leading_iter->pt, scale_factor[i]);
		      h_truth_subleading->Fill(subleading_iter->pt, scale_factor[i]);		      
		    }
		}

	    }
            
	  // Getting the reco jets
	  for (int j = 0; j < nrecojets;j++)
	    {

	      if (reco_jet_pt->at(j) < reco_cut) continue;
	      if (reco_jet_e->at(j) < 0) continue;
	      if (reco_jet_e_unsub->at(j) < 0) continue;
	      if (fabs(reco_jet_eta->at(j)) > etacut) continue;

	      struct jet tempjet;
	      tempjet.pt = reco_jet_pt->at(j);
	      tempjet.eta = reco_jet_eta->at(j);
	      tempjet.phi = reco_jet_phi->at(j);
	      tempjet.id = j;
	      myrecojets.push_back(tempjet);

	    }

	  int nnrecojets = myrecojets.size();
	  double dphir = 0;
	  if (myrecojets.size() >= 2) 
	    {
	      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	      auto leading_iter = myrecojets.begin();
	      auto subleading_iter = myrecojets.begin() + 1;
	      if (subleading_iter != myrecojets.end())
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
		  if (fabs(dphir) >= dphicut && leading_iter->pt >= 20 && subleading_iter->pt >= 8)
		    {
		      h_reco_leading->Fill(leading_iter->pt, scale_factor[i]);
		      h_reco_subleading->Fill(subleading_iter->pt, scale_factor[i]);		      
		    }
		}

	    }
	}
    }
  

  int entriesd = tnd->GetEntries();
 
  for (int i = 0; i < entriesd; i++)
    {
      tnd->GetEntry(i);
      float maxi = std::max(pt1_data, pt2_data);
      float mini = std::min(pt1_data, pt2_data);
	 
      bool data_good = (maxi >= 20 && mini >= 8 && dphi_data > dphicut);

      if (!data_good) continue;
      
      h_data_leading->Fill(maxi);
      h_data_subleading->Fill(mini);

    }


  TCanvas *c = new TCanvas("c","c", 700, 500);
  c->SetLogy();

  h_reco_leading->Scale(1./h_reco_leading->Integral(),"width");
  h_reco_subleading->Scale(1./h_reco_subleading->Integral(),"width");

  h_truth_leading->Scale(1./h_truth_leading->Integral(),"width");
  h_truth_subleading->Scale(1./h_truth_subleading->Integral(),"width");

  h_data_leading->Scale(1./h_data_leading->Integral(),"width");
  h_data_subleading->Scale(1./h_data_subleading->Integral(),"width"); 


  dlutility::SetLineAtt(h_truth_leading, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_truth_leading, kBlack, 1, 20);
  dlutility::SetLineAtt(h_reco_leading, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_reco_leading, kRed, 1, 20);
  dlutility::SetLineAtt(h_data_leading, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_data_leading, kBlue, 1, 20);
  dlutility::SetLineAtt(h_truth_subleading, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_truth_subleading, kBlack, 1, 21);
  dlutility::SetLineAtt(h_reco_subleading, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_reco_subleading, kRed, 1, 21);
  dlutility::SetLineAtt(h_data_subleading, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_data_subleading, kBlue, 1, 21);


  h_truth_leading->SetMaximum(3);
  h_truth_leading->SetMinimum(0.0001);
  h_truth_leading->SetTitle(";#it{p}_{T} [GeV]; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#it{p}_{T}}");
  h_truth_leading->Draw("p");
  h_truth_subleading->Draw("same p");
  
  h_reco_leading->Draw("same p");
  h_reco_subleading->Draw("same p");
  h_data_leading->Draw("same p");
  h_data_subleading->Draw("same p");
  dlutility::DrawSPHENIX(0.6, 0.85);
  TLegend *l = new TLegend(0.67, 0.5, 0.85, 0.75);
  l->SetLineWidth(0);
  l->SetTextSize(0.03);
  l->AddEntry(h_truth_leading,"Truth leading");
  l->AddEntry(h_truth_subleading,"Truth subleading");
  l->AddEntry(h_reco_leading,"Reco-sim leading");
  l->AddEntry(h_reco_subleading,"reco-sim subleading");
  l->AddEntry(h_data_leading,"Data leading");
  l->AddEntry(h_data_subleading,"Data subleading");
  l->Draw("same")  ;
  c->Print("ptspectra.pdf");
}
