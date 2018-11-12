#define CAF_cxx
#ifdef CAF_cxx

#include "CAF.h"

CAF::CAF( std::string filename )
{
  cafFile = new TFile( filename.c_str(), "RECREATE" );
  cafMVA = new TTree( "caf", "caf" );
  cafPOT = new TTree( "meta", "meta" );

  cafMVA->Branch( "run", &run, "run/I" );
  cafMVA->Branch( "subrun", &subrun, "subrun/I" );
  cafMVA->Branch( "event", &event, "event/I" );
  cafMVA->Branch( "isFD", &isFD, "isFD/I" );
  cafMVA->Branch( "isFHC", &isFHC, "isFHC/I" );
  cafMVA->Branch( "isCC", &isCC, "isCC/I" );

  cafMVA->Branch( "nuPDG", &neutrinoPDG, "nuPDG/I" );
  cafMVA->Branch( "nuPDGunosc", &neutrinoPDGunosc, "nuPDGunosc/I");
  cafMVA->Branch( "NuMomX", &NuMomX, "NuMomX/D" );
  cafMVA->Branch( "NuMomY", &NuMomY, "NuMomY/D" );
  cafMVA->Branch( "NuMomZ", &NuMomZ, "NuMomZ/D" );
  cafMVA->Branch( "Ev", &Ev, "Ev/D" );
  cafMVA->Branch( "mode", &mode, "mode/I" );
  cafMVA->Branch( "LepPDG", &LepPDG, "LepPDG/I" );
  cafMVA->Branch( "LepMomX", &LepMomX, "LepMomX/D" );
  cafMVA->Branch( "LepMomY", &LepMomY, "LepMomY/D" );
  cafMVA->Branch( "LepMomZ", &LepMomZ, "LepMomZ/D" );
  cafMVA->Branch( "LepE", &LepE, "LepE/D" );
  cafMVA->Branch( "LepNuAngle", &LepNuAngle, "LepNuAngle/D" );
  cafMVA->Branch( "Q2", &Q2, "Ev/D" );
  cafMVA->Branch( "W", &W, "Ev/D" );
  cafMVA->Branch( "X", &X, "Ev/D" );
  cafMVA->Branch( "Y", &Y, "Ev/D" );

  cafMVA->Branch( "nP", &nP, "nP/I" );
  cafMVA->Branch( "nN", &nN, "nN/I" );
  cafMVA->Branch( "nipip", &nipip, "nipip/I" );
  cafMVA->Branch( "nipim", &nipim, "nipim/I" );
  cafMVA->Branch( "nipi0", &nipi0, "nipi0/I" );
  cafMVA->Branch( "nikp", &nikp, "nikp/I" );
  cafMVA->Branch( "nikm", &nikm, "nikm/I" );
  cafMVA->Branch( "nik0", &nik0, "nik0/I" );
  cafMVA->Branch( "niem", &niem, "niem/I" );
  cafMVA->Branch( "niother", &niother, "niother/I" );
  cafMVA->Branch( "nNucleus", &nNucleus, "nNucleus/I" );
  cafMVA->Branch( "nUNKNOWN", &nUNKNOWN, "nUNKNOWN/I" );

  cafMVA->Branch("eP",        &eP,         "eP/D");
  cafMVA->Branch("eN",        &eN,         "eN/D");
  cafMVA->Branch("ePip",      &ePip,       "ePip/D");
  cafMVA->Branch("ePim",      &ePim,       "ePim/D");
  cafMVA->Branch("ePi0",      &ePi0,       "ePi0/D");
  cafMVA->Branch("eOther",    &eOther,     "eOther/D");

  cafMVA->Branch( "det_x", &det_x, "det_x/D" );
  cafMVA->Branch( "vtx_x", &vtx_x, "vtx_x/D" );
  cafMVA->Branch( "vtx_y", &vtx_y, "vtx_y/D" );
  cafMVA->Branch( "vtx_z", &vtx_z, "vtx_z/D" );

  cafMVA->Branch( "Ev_reco", &Ev_reco, "Ev_reco/D" );
  cafMVA->Branch( "Elep_reco", &Elep_reco, "Elep_reco/D" );
  cafMVA->Branch( "reco_numu", &reco_numu, "reco_numu/I" );
  cafMVA->Branch( "reco_nue", &reco_nue, "reco_nue/I" );
  cafMVA->Branch( "reco_nc", &reco_nc, "reco_nc/I" );
  cafMVA->Branch( "reco_q", &reco_q, "reco_q/I" );
  cafMVA->Branch( "muon_contained", &muon_contained, "muon_contained/I" );
  cafMVA->Branch( "muon_tracker", &muon_tracker, "muon_tracker/I" );
  cafMVA->Branch( "muon_ecal", &muon_ecal, "muon_ecal/I" );
  cafMVA->Branch( "muon_exit", &muon_exit, "muon_exit/I" );
  cafMVA->Branch( "Ehad_veto", &Ehad_veto, "Ehad_veto/D" );

  cafPOT->Branch( "pot", &pot, "pot/D" );
  cafPOT->Branch( "run", &meta_run, "run/I" );
  cafPOT->Branch( "subrun", &meta_subrun, "subrun/I" );
  cafPOT->Branch( "version", &version, "version/I" );
}

CAF::~CAF() {}

void CAF::fill()
{
  cafMVA->Fill();
}

void CAF::Print()
{
  printf( "Event %d:\n", event );
  printf( "   Truth: Ev = %3.3f Elep = %3.3f Q2 = %3.3f W = %3.3f x = %3.3f y = %3.3f lepton %d mode %d\n", Ev, LepE, Q2, W, X, Y, LepPDG, mode );
  printf( "    Reco: Ev = %3.3f Elep = %3.3f q %d mu/e/nc %d%d%d cont/trk/ecal/exit %d%d%d%d had veto %2.1f\n\n", Ev_reco, Elep_reco, reco_q, reco_numu, reco_nue, reco_nc, muon_contained, muon_tracker, muon_ecal, muon_exit, Ehad_veto );
}

void CAF::fillPOT()
{
  printf( "Filling metadata\n" );
  cafPOT->Fill();
}

void CAF::write()
{
  cafFile->cd();
  cafMVA->Write();
  cafPOT->Write();
}

void CAF::addRWbranch( int parId, std::string name, std::string wgt_var, std::vector<double> &vars )
{
  cafMVA->Branch( Form("%s_nshifts", name.c_str()), &nwgt[parId], Form("%s_nshifts/I", name.c_str()) );
  cafMVA->Branch( Form("%s_cv%s", name.c_str(), wgt_var.c_str()), &cvwgt[parId], Form("%s_cv%s/D", name.c_str(), wgt_var.c_str()) );
  cafMVA->Branch( Form("%s_%s", wgt_var.c_str(), name.c_str()), wgt[parId], Form("%s_%s[%s_nshifts]/D", wgt_var.c_str(), name.c_str(), name.c_str()) );
}

#endif
