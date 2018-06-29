#define CAF_cxx
#ifdef CAF_cxx

#include "CAF.h"

CAF::CAF( std::string filename )
{
  cafFile = new TFile( filename.c_str(), "RECREATE" );
  mvaselect = cafFile->mkdir( "mvaselect" );
  mvaselect->cd();
  cafMVA = new TTree( "MVASelection", "MVASelection" );
  cafPOT = new TTree( "pottree", "pottree" );

  cafMVA->Branch( "Ev_reco", &Ev_reco, "Ev_reco/D" );
  cafMVA->Branch( "Ev", &Ev, "Ev/D" );
  cafMVA->Branch( "Elep_reco", &Elep_reco, "Elep_reco/D" );
  cafMVA->Branch( "Elep", &Elep, "Elep/D" );
  cafMVA->Branch( "Q2", &Q2, "Q2/D" );
  cafMVA->Branch( "mvaresult", &mvaresult, "mvaresult/D" );
  cafMVA->Branch( "numu_pid", &numu_pid, "numu_pid/D" );
  cafMVA->Branch( "nue_pid", &nue_pid, "nue_pid/D" );
  cafMVA->Branch( "ccnc", &ccnc, "ccnc/I" );
  cafMVA->Branch( "LepPDG", &LepPDG, "LepPDG/I" );
  cafMVA->Branch( "beamPdg", &beamPdg, "beamPdg/I" );
  cafMVA->Branch( "neu", &neu, "neu/I" );
  cafMVA->Branch( "nipip", &nipip, "nipip/I" );
  cafMVA->Branch( "nipim", &nipim, "nipim/I" );
  cafMVA->Branch( "nipi0", &nipi0, "nipi0/I" );
  cafMVA->Branch( "mode", &mode, "mode/I" );
  cafMVA->Branch( "isFD", &isFD, "isFD/I" );
  cafMVA->Branch( "isFHC", &isFHC, "isFHC/I" );
  cafMVA->Branch( "reco_q", &reco_q, "reco_q/I" );
  cafMVA->Branch( "run", &run, "run/I" );

  cafPOT->Branch( "pot", &pot, "pot/D" );
  cafPOT->Branch( "run", &meta_run, "run/I" );
}

CAF::~CAF() {}

void CAF::fill()
{
  cafMVA->Fill();
}

void CAF::Print()
{
  printf( "Event FD %d FHC %d:\n", isFD, isFHC );
  printf( "   Truth: Ev = %3.3f Elep = %3.3f Q2 = %3.3f, lepton %d mode %d\n", Ev, Elep, Q2, LepPDG, mode );
  printf( "    Reco: Ev = %3.3f Elep = %3.3f q %d numu %1.1f nue %1.1f\n\n", Ev_reco, Elep_reco, reco_q, numu_pid, nue_pid );
  printf( "    weights: %d\n", nwgt[0] );
  for( int i = 0; i < nwgt[0]; ++i ) printf( "%f ", wgt[0][i] );
  printf( "\n" );
}

void CAF::fillPOT()
{
  cafPOT->Fill();
}

void CAF::write()
{
  mvaselect->cd();
  cafMVA->Write();
  cafPOT->Write();
}

void CAF::addRWbranch( int parId, std::string name, std::vector<double> &vars )
{
  cafMVA->Branch( Form("%s_nshifts", name.c_str()), &nwgt[parId], Form("%s_nshifts/I", name.c_str()) );
  cafMVA->Branch( Form("wgt_%s", name.c_str()), wgt[parId], Form("wgt_%s[%s_nshifts]/D", name.c_str(), name.c_str()) );
}

#endif
