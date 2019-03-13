#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TVectorT.h"
#include "TCanvas.h"

const int n_Ebins = 22;
const int n_ybins = 7;
double Ebins[23] = { 0., 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 5., 5.5, 6., 7., 8., 10. };
double ybins[8] = { 0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0 };

double plbins[29] = { -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 9., 10., 10.5 };
double ptbins[17] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.25 };
double hbins[22] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 5.25 };

const int nu = 100;

int get1Dbin( int bx, int by )
{
    return (bx-1) * n_ybins + by;
}

void get2Dbins( int b1D, int &eb, int &yb )
{
  eb = ((b1D-1) / n_ybins) + 1;
  yb = ((b1D-1) % n_ybins) + 1;
}

double getP( double e, int pdg )
{
  if( abs(pdg) == 11 ) return sqrt(e*e - 0.000511*0.000511);
  else if( abs(pdg) == 13 ) return sqrt(e*e - 0.105658*0.105658);
  else if( abs(pdg) == 211 ) return sqrt(e*e - 0.13957*0.13957);
  else if( abs(pdg) == 321 ) return sqrt(e*e - 0.49366*0.49366);
  else if( abs(pdg) == 2212 ) return sqrt((e+0.93827)*(e+0.93827) - 0.93827*0.93827);
  else return e;
}

double getE( double p, int pdg )
{
  if( abs(pdg) == 11 ) return sqrt(p*p + 0.000511*0.000511);
  else if( abs(pdg) == 13 ) return sqrt(p*p + 0.105658*0.105658);
  else if( abs(pdg) == 211 ) return sqrt(p*p + 0.13957*0.13957);
  else if( abs(pdg) == 321 ) return sqrt(p*p + 0.49366*0.49366);
  else if( abs(pdg) == 2212 ) return sqrt(p*p + 0.93827*0.93827) - 0.93827;
  else return p;
}

void fix( TMatrixD &cov )
{
  // get rid of super-tiny negative eigenvalues
  TVectorD evals;
  const TMatrixD evecs = cov.EigenVectors( evals );
  TMatrixD evalmat(evals.GetNrows(), evals.GetNrows());
  for( int i = 0; i < evals.GetNrows(); ++i ) {
    evalmat(i, i) = std::max(1e-14, evals[i]);
  }

  TMatrixD evecs_inv(TMatrixD::kTransposed, evecs);
  cov = evecs*evalmat*evecs_inv;
}

void makeCov()
{

  TRandom3 * rando = new TRandom3(12345);

  // Get acceptance uncertainty histograms
  TFile * tf_AccUnc = new TFile( "/dune/data/users/marshalc/CAFs/mcc11_v3/ND_eff_syst.root" );
  TH2D * hMuUnc = (TH2D*) tf_AccUnc->Get( "unc" );
  TH1D * hHadUnc = (TH1D*) tf_AccUnc->Get( "hunc" );

  TChain * cafTree = new TChain( "cafTree", "cafTree" );
  cafTree->Add( "/pnfs/dune/data/users/LBL_TDR/v4/ND_FHC_*.root" );

  TFile * fdFileMu = new TFile( "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_nonswap.root" );
  TFile * fdFileE = new TFile( "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_nueswap.root" );
  TTree * cafFDmu = (TTree*) fdFileMu->Get( "cafTree" );
  TTree * cafFDe  = (TTree*) fdFileE->Get( "cafTree" );

  // gas TPC files
  TFile * gasFile = new TFile( "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/NDgas_FHC.root" );
  TTree * gasCaf = (TTree*) gasFile->Get( "cafTree" );

  // CAF variables for ND & FD
  double vtx_x, vtx_y, vtx_z, LepE, LepNuAngle, Ev, Ev_reco, Elep_reco;
  double eRecoP, eRecoN, eRecoPip, eRecoPim, eRecoPi0;
  double eP, eN, ePip, ePim, ePi0;
  int LepPDG;

  // Additional gas TPC variables for pion counting
  int gastpc_pi_min_mult, gastpc_pi_pl_mult, nFSP;
  int pdg[100];
  double trkLen[100], partEvReco[100];

  // CAF variables for ND only
  int reco_numu, muon_contained, muon_tracker, reco_q;
  double Ehad_veto;

  // CAF variables for FD only
  double cvnnumu, cvnnue;

  // Universe histograms in analysis bins
  TH2D * histCV = new TH2D( "histCV", ";Reconstructed E_{#nu};Reconstructed y", n_Ebins, Ebins, n_ybins, ybins );
  TH1D * histCV_FDmu = new TH1D( "histCV_FDmu", ";Reconstructed E_{#nu}", n_Ebins, Ebins );
  TH1D * histCV_FDe = new TH1D( "histCV_FDe", ";Reconstructed E_{#nu}", n_Ebins, Ebins );
  TH2D * histCV_gas = new TH2D( "histCV_gas", ";Number of charged pions;Reconstructed E_{#nu}", 3, 0., 3., n_Ebins, Ebins );
  TH2D * hists[nu];
  TH1D * hists_FDmu[nu];
  TH1D * hists_FDe[nu];
  TH2D * histsAccOnly[nu];
  TH2D * histsEscaleOnly[nu];
  TH2D * hists_gas[nu];

  // Uncertainties for each universe -- ND
  TH2D * muAccThrow[nu];
  TH1D * hAccThrow[nu];
  TF1 * EtotThrow[nu];
  TF1 * EmuThrowLAr[nu];
  TF1 * EmuThrowGAr[nu];
  TF1 * EhadThrow[nu];
  TF1 * EEMThrow[nu];
  TF1 * EneutThrow[nu];
  double EmuRes[nu];
  double EhadRes[nu];
  double EEMRes[nu];
  double EneutRes[nu];

  // Gas TPC uncertainties
  TF1 * Pscale[nu];
  double trkThreshold[nu];

  // FD
  TF1 * EtotThrowFD[nu];
  TF1 * EmuThrowFD[nu];
  TF1 * EhadThrowFD[nu];
  TF1 * EEMThrowFD[nu];
  TF1 * EneutThrowFD[nu];
  double EmuResFD[nu];
  double EhadResFD[nu];
  double EEMResFD[nu];
  double EneutResFD[nu];

  for( int u = 0; u < nu; ++u ) {
    hists[u] = new TH2D( Form("h%03d", u), ";Reco E_{#nu} (GeV);Reco y", n_Ebins, Ebins, n_ybins, ybins );
    histsAccOnly[u] = new TH2D( Form("hAO%03d", u), ";Reco E_{#nu} (GeV);Reco y", n_Ebins, Ebins, n_ybins, ybins );
    histsEscaleOnly[u] = new TH2D( Form("hEO%03d", u), ";Reco E_{#nu} (GeV);Reco y", n_Ebins, Ebins, n_ybins, ybins );

    hists_FDmu[u] = new TH1D( Form("hFDmu%03d", u), ";Reco E_{#nu} (GeV)", n_Ebins, Ebins );
    hists_FDe[u] = new TH1D( Form("hFDe%03d", u), ";Reco E_{#nu} (GeV)", n_Ebins, Ebins );

    muAccThrow[u] = new TH2D( Form("muAccThrow%03d", u), ";Muon p_{L};Muon p_{T}", 28, plbins, 16, ptbins );
    hAccThrow[u] = new TH1D( Form("hAccThrow%03d", u), ";Hadronic energy", 21, hbins );

    hists_gas[u] = new TH2D( Form("hGas%03d",u), ";Number of charged pions;Reconstructed E_{#nu}", 3, 0., 3., n_Ebins, Ebins );

    // All the energy systematics as TF1 vs. energy
    EtotThrow[u]   = new TF1( Form("Etot%03d", u),   "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EmuThrowLAr[u] = new TF1( Form("EmuLAr%03d", u), "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EmuThrowGAr[u] = new TF1( Form("EmuGAr%03d", u), "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EhadThrow[u]   = new TF1( Form("Ehad%03d", u),   "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EEMThrow[u]    = new TF1( Form("EEM%03d", u),    "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EneutThrow[u]  = new TF1( Form("Eneut%03d", u),  "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    // gas TPC
    Pscale[u] = new TF1( Form("Pscale%03d", u), "[0] + [1]*x + [2]*pow(x,2.)", 0., 100. );

    // FD
    EtotThrowFD[u]   = new TF1( Form("EtotFD%03d", u),   "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EmuThrowFD[u] = new TF1( Form("EmuFD%03d", u), "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EhadThrowFD[u]   = new TF1( Form("EhadFD%03d", u),   "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EEMThrowFD[u]    = new TF1( Form("EEMFD%03d", u),    "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );
    EneutThrowFD[u]  = new TF1( Form("EneutFD%03d", u),  "[0] + [1]*x + [2]*pow(x+0.1,-0.5)", 0., 100. );

    EmuRes[u] = rando->Gaus(0., 0.1);
    EhadRes[u] = rando->Gaus(0., 0.1);
    EEMRes[u] = rando->Gaus(0., 0.1);
    EneutRes[u] = rando->Gaus(0., 0.3);
    EmuResFD[u] = rando->Gaus(0., 0.1);
    EhadResFD[u] = rando->Gaus(0., 0.1);
    EEMResFD[u] = rando->Gaus(0., 0.1);
    EneutResFD[u] = rando->Gaus(0., 0.3);

    trkThreshold[u] = rando->Gaus( 6., 3. ); // 5 MeV threshold, 2.5 MeV width
    if( trkThreshold[u] < 1. ) trkThreshold[u] = 1.; // truncate gaussian at 1 MeV threshold

    // set the parameters
    EtotThrow[u]->SetParameter( 0, rando->Gaus(0., 0.02) );
    EtotThrow[u]->SetParameter( 1, rando->Gaus(0., 0.01) );
    EtotThrow[u]->SetParameter( 2, rando->Gaus(0., 0.02) );

    EmuThrowLAr[u]->SetParameter( 0, rando->Gaus(0., 0.02) );
    EmuThrowLAr[u]->SetParameter( 1, rando->Gaus(0., 0.005) );
    EmuThrowLAr[u]->SetParameter( 2, rando->Gaus(0., 0.02) );

    EmuThrowGAr[u]->SetParameter( 0, rando->Gaus(0., 0.01) );
    EmuThrowGAr[u]->SetParameter( 1, rando->Gaus(0., 0.00001) );
    EmuThrowGAr[u]->SetParameter( 2, rando->Gaus(0., 0.01) );

    EhadThrow[u]->SetParameter( 0, rando->Gaus(0., 0.05) );
    EhadThrow[u]->SetParameter( 1, rando->Gaus(0., 0.05) );
    EhadThrow[u]->SetParameter( 2, rando->Gaus(0., 0.05) );

    EEMThrow[u]->SetParameter( 0, rando->Gaus(0., 0.05) );
    EEMThrow[u]->SetParameter( 1, rando->Gaus(0., 0.05) );
    EEMThrow[u]->SetParameter( 2, rando->Gaus(0., 0.05) );

    EneutThrow[u]->SetParameter( 0, rando->Gaus(0., 0.2) );
    EneutThrow[u]->SetParameter( 1, rando->Gaus(0., 0.3) );
    EneutThrow[u]->SetParameter( 2, rando->Gaus(0., 0.3) );

    // gas TPC
    Pscale[u]->SetParameter( 0, rando->Gaus(0., 0.01) );
    Pscale[u]->SetParameter( 1, rando->Gaus(0., 0.002) );
    Pscale[u]->SetParameter( 2, rando->Gaus(0., 0.001) );

    // FD
    EtotThrowFD[u]->SetParameter( 0, rando->Gaus(0., 0.02) );
    EtotThrowFD[u]->SetParameter( 1, rando->Gaus(0., 0.01) );
    EtotThrowFD[u]->SetParameter( 2, rando->Gaus(0., 0.02) );

    EmuThrowFD[u]->SetParameter( 0, rando->Gaus(0., 0.02) );
    EmuThrowFD[u]->SetParameter( 1, rando->Gaus(0., 0.005) );
    EmuThrowFD[u]->SetParameter( 2, rando->Gaus(0., 0.02) );

    EhadThrowFD[u]->SetParameter( 0, rando->Gaus(0., 0.05) );
    EhadThrowFD[u]->SetParameter( 1, rando->Gaus(0., 0.05) );
    EhadThrowFD[u]->SetParameter( 2, rando->Gaus(0., 0.05) );

    EEMThrowFD[u]->SetParameter( 0, rando->Gaus(0., 0.05) );
    EEMThrowFD[u]->SetParameter( 1, rando->Gaus(0., 0.05) );
    EEMThrowFD[u]->SetParameter( 2, rando->Gaus(0., 0.05) );

    EneutThrowFD[u]->SetParameter( 0, rando->Gaus(0., 0.2) );
    EneutThrowFD[u]->SetParameter( 1, rando->Gaus(0., 0.3) );
    EneutThrowFD[u]->SetParameter( 2, rando->Gaus(0., 0.3) );

  }

  // Build throw histograms for acceptance uncertainties        
  // for each bin, throw the uncertainty, as if totally uncorrelated bin to bin
  for( int u = 0; u < nu; ++u ) {
    for( int b = 1; b <= hHadUnc->GetNbinsX(); ++b ) {
      if( hHadUnc->GetBinContent(b) > 0. ) {
        hAccThrow[u]->SetBinContent( b, rando->Gaus(0., hHadUnc->GetBinContent(b)) );
      }
    }
    for( int bx = 1; bx <= hMuUnc->GetNbinsX(); ++bx ) {
      for( int by = 1; by <= hMuUnc->GetNbinsY(); ++by ) {
        if( hMuUnc->GetBinContent(bx, by) > 1.E-6 ) {
          muAccThrow[u]->SetBinContent( bx, by, rando->Gaus(0., hMuUnc->GetBinContent(bx,by)) );
        }
      }
    }
    // for smoothing, set unfilled bins to the value of their neighbors
    for( int bx = 1; bx <= hMuUnc->GetNbinsX(); ++bx ) {
      for( int by = 1; by <= hMuUnc->GetNbinsY(); ++by ) {
        if( hMuUnc->GetBinContent(bx,by) < 1.E-6 ) { // bin not filled
          int near_by = by - 1;
          while( near_by && hMuUnc->GetBinContent(bx, near_by) < 1.E-6 ) --near_by;
          if( !near_by ) {
            int near_bx = bx+1;
            while( hMuUnc->GetBinContent(near_bx, by) < 1.E-6 ) ++near_bx;
            muAccThrow[u]->SetBinContent( bx, by, muAccThrow[u]->GetBinContent(near_bx, by) );
          } else {
            muAccThrow[u]->SetBinContent( bx, by, muAccThrow[u]->GetBinContent(bx, near_by) );
          }
        }
      }
    }
    // Now smooth it, so that it allows any smooth function in the envelope of the uncertainty          
    hAccThrow[u]->Smooth(2);
    muAccThrow[u]->Smooth(1);
  }

  // some validation plots
  TH2D * val_Ev[nu];
  TH2D * val_y[nu];

  TH2D * val_npi_gas[nu];
  TH2D * val_Ev_gas[nu];
  // add gas tpc validation plots
  for( int u = 0; u < nu; ++u ) {
    val_Ev[u] = new TH2D( Form("val_Ev_%03d",u), ";Reco E_{#nu};Shifted E_{#nu}", 100, 0., 10., 100, 0., 10. );
    val_y[u] = new TH2D( Form("val_y_%03d",u), ";Reco y;Shifted y", 100, 0., 1., 100, 0., 1. );

    val_npi_gas[nu] = new TH2D( Form("val_gas_npi_%03d",u), ";CV N_{#pi};Shifted N_{#pi}", 3, 0., 3. );
    val_Ev_gas[nu] = new TH2D( Form("val_gas_Ev_%03d",u), ";CV Ev;Shifted Ev", 100, 0., 10., 100, 0., 10. );
  }

    
  // Loop over ND events and fill the analysis bin histograms
  cafTree->SetBranchAddress( "vtx_x", &vtx_x );
  cafTree->SetBranchAddress( "vtx_y", &vtx_y );
  cafTree->SetBranchAddress( "vtx_z", &vtx_z );
  cafTree->SetBranchAddress( "LepE", &LepE );
  cafTree->SetBranchAddress( "LepNuAngle", &LepNuAngle );
  cafTree->SetBranchAddress( "Ehad_veto", &Ehad_veto );
  cafTree->SetBranchAddress( "Ev", &Ev );
  cafTree->SetBranchAddress( "Ev_reco", &Ev_reco );
  cafTree->SetBranchAddress( "Elep_reco", &Elep_reco );
  cafTree->SetBranchAddress( "LepPDG", &LepPDG );
  cafTree->SetBranchAddress( "reco_numu", &reco_numu );
  cafTree->SetBranchAddress( "muon_contained", &muon_contained );
  cafTree->SetBranchAddress( "muon_tracker", &muon_tracker );
  cafTree->SetBranchAddress( "reco_q", &reco_q );
  cafTree->SetBranchAddress( "eRecoP", &eRecoP );
  cafTree->SetBranchAddress( "eRecoN", &eRecoN );
  cafTree->SetBranchAddress( "eRecoPip", &eRecoPip );
  cafTree->SetBranchAddress( "eRecoPim", &eRecoPim );
  cafTree->SetBranchAddress( "eRecoPi0", &eRecoPi0 );
  cafTree->SetBranchAddress( "eP", &eP );
  cafTree->SetBranchAddress( "eN", &eN );
  cafTree->SetBranchAddress( "ePip", &ePip );
  cafTree->SetBranchAddress( "ePim", &ePim );
  cafTree->SetBranchAddress( "ePi0", &ePi0 );

  cafTree->SetBranchStatus( "*", 0 );
  cafTree->SetBranchStatus( "vtx_x", 1 );
  cafTree->SetBranchStatus( "vtx_y", 1 );
  cafTree->SetBranchStatus( "vtx_z", 1 );
  cafTree->SetBranchStatus( "LepPDG", 1 );
  cafTree->SetBranchStatus( "LepE", 1 );
  cafTree->SetBranchStatus( "LepNuAngle", 1 );
  cafTree->SetBranchStatus( "reco_numu", 1 );
  cafTree->SetBranchStatus( "muon_contained", 1 );
  cafTree->SetBranchStatus( "muon_tracker", 1 );
  cafTree->SetBranchStatus( "reco_q", 1 );
  cafTree->SetBranchStatus( "Ehad_veto", 1 );
  cafTree->SetBranchStatus( "Ev", 1 );
  cafTree->SetBranchStatus( "Ev_reco", 1 );
  cafTree->SetBranchStatus( "Elep_reco", 1 );
  cafTree->SetBranchStatus( "eRecoP", 1 );
  cafTree->SetBranchStatus( "eRecoN", 1 );
  cafTree->SetBranchStatus( "eRecoPip", 1 );
  cafTree->SetBranchStatus( "eRecoPim", 1 );
  cafTree->SetBranchStatus( "eRecoPi0", 1 );
  cafTree->SetBranchStatus( "eP", 1 );
  cafTree->SetBranchStatus( "eN", 1 );
  cafTree->SetBranchStatus( "ePip", 1 );
  cafTree->SetBranchStatus( "ePim", 1 );
  cafTree->SetBranchStatus( "ePi0", 1 );
  int N = cafTree->GetEntries();
  for( int ii = 0; ii < N; ++ii ) {
    cafTree->GetEntry(ii);

    if( ii % 100000 == 0 ) printf( "Event %d of %d...\n", ii, N );

    // FV cut
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z < 350. ) continue;
    // true numu CC cut
    if( LepPDG != 13 ) continue;
    // reco numu CC cut
    bool numuCC = (reco_numu && reco_q == -1 && (muon_contained || muon_tracker));
    if( !numuCC ) continue;

    // determine quantities for acceptance uncertainties, including overflow bins
    double p = sqrt(LepE*LepE - 0.105658*0.105658);
    double pl = p*cos(LepNuAngle);
    if( pl > 10.25 ) pl = 10.25; // overflow
    double pt = p*sin(LepNuAngle);
    if( pt > 3.2 ) pt = 3.2; // overflow
    double ehad = Ev_reco - Elep_reco;
    if( ehad > 5.1 ) ehad = 5.1;

    histCV->Fill( Ev_reco, (Ev_reco-Elep_reco)/Ev_reco, 1. );
    for( int u = 0; u < nu; ++u ) {
      double wgt_mu = 1. + muAccThrow[u]->GetBinContent( muAccThrow[u]->FindBin(pl,pt) );
      double wgt_had = 1. + hAccThrow[u]->GetBinContent( hAccThrow[u]->FindBin(ehad) );

      // determine the shifted energies
      double shiftChargedHad = EhadThrow[u]->Eval(eRecoP + eRecoPip + eRecoPim);
      double shiftEM = EEMThrow[u]->Eval(eRecoPi0);
      double shiftN = EneutThrow[u]->Eval(eRecoN);
      double shiftMu = ( muon_contained ? EmuThrowLAr[u]->Eval(Elep_reco) : EmuThrowGAr[u]->Eval(Elep_reco) );
      double shiftTot = EtotThrow[u]->Eval(Ev_reco - Elep_reco);

      double Ehad_reco_shift = Ev_reco - Elep_reco;
      Ehad_reco_shift += shiftChargedHad*(eRecoP + eRecoPip + eRecoPim);
      Ehad_reco_shift += shiftEM*eRecoPi0;
      Ehad_reco_shift += shiftN*eRecoN;

      double Elep_reco_shift = Elep_reco*(1.+shiftMu);

      Ehad_reco_shift *= (1.+shiftTot);
      if( muon_contained ) Elep_reco_shift *= (1.+shiftTot);

      // resolution uncertainties
      Elep_reco_shift += (LepE - Elep_reco)*EmuRes[u];
      Ehad_reco_shift += ((eP + ePip + ePim) - (eRecoP + eRecoPip + eRecoPim))*EhadRes[u];
      Ehad_reco_shift += (ePi0 - eRecoPi0)*EEMRes[u];
      Ehad_reco_shift += (eN - eRecoN)*EneutRes[u];

      double Ev_reco_shift = Elep_reco_shift + Ehad_reco_shift;

      hists[u]->Fill( Ev_reco_shift, Ehad_reco_shift/Ev_reco_shift, wgt_mu*wgt_had );
      histsAccOnly[u]->Fill( Ev_reco, (Ev_reco - Elep_reco)/Ev_reco, wgt_mu*wgt_had );
      histsEscaleOnly[u]->Fill( Ev_reco_shift, Ehad_reco_shift/Ev_reco_shift, 1. );

      val_Ev[u]->Fill( Ev_reco, Ev_reco_shift );
      val_y[u]->Fill( (Ev_reco - Elep_reco)/Ev_reco, Ehad_reco_shift/Ev_reco_shift );
    }
  }

  // Gas TPC loop
  // Loop over ND events and fill the analysis bin histograms
  gasCaf->SetBranchAddress( "vtx_x", &vtx_x );
  gasCaf->SetBranchAddress( "vtx_y", &vtx_y );
  gasCaf->SetBranchAddress( "vtx_z", &vtx_z );
  gasCaf->SetBranchAddress( "Ev", &Ev );
  gasCaf->SetBranchAddress( "Ev_reco", &Ev_reco );
  gasCaf->SetBranchAddress( "LepPDG", &LepPDG );
  gasCaf->SetBranchAddress( "reco_numu", &reco_numu );
  gasCaf->SetBranchAddress( "reco_q", &reco_q );
  gasCaf->SetBranchAddress( "gastpc_pi_min_mult", &gastpc_pi_min_mult );
  gasCaf->SetBranchAddress( "gastpc_pi_pl_mult", &gastpc_pi_pl_mult );
  gasCaf->SetBranchAddress( "nFSP", &nFSP );
  gasCaf->SetBranchAddress( "pdg", pdg );
  gasCaf->SetBranchAddress( "trkLen", trkLen );
  gasCaf->SetBranchAddress( "partEvReco", partEvReco );

  gasCaf->SetBranchStatus( "*", 0 );
  gasCaf->SetBranchStatus( "vtx_x", 1 );
  gasCaf->SetBranchStatus( "vtx_y", 1 );
  gasCaf->SetBranchStatus( "vtx_z", 1 );
  gasCaf->SetBranchStatus( "LepPDG", 1 );
  gasCaf->SetBranchStatus( "reco_numu", 1 );
  gasCaf->SetBranchStatus( "reco_q", 1 );
  gasCaf->SetBranchStatus( "Ev", 1 );
  gasCaf->SetBranchStatus( "Ev_reco", 1 );
  gasCaf->SetBranchStatus( "gastpc_pi_min_mult", 1 );
  gasCaf->SetBranchStatus( "gastpc_pi_pl_mult", 1 );
  gasCaf->SetBranchStatus( "nFSP", 1 );
  gasCaf->SetBranchStatus( "pdg", 1 );
  gasCaf->SetBranchStatus( "trkLen", 1 );
  gasCaf->SetBranchStatus( "partEvReco", 1 );

  int N = gasCaf->GetEntries();
  for( int ii = 0; ii < N; ++ii ) {
    gasCaf->GetEntry(ii);

    if( ii % 100000 == 0 ) printf( "Event %d of %d...\n", ii, N );

    // FV cut
    if( abs(vtx_x) > 200. ) continue; // endcap cut
    double r = sqrt((vtx_z-952.5)*(vtx_z-952.5) + (vtx_y+72.5)*(vtx_y+72.5));
    if( r > 200. ) continue; // circle cut

    // true numu CC cut
    if( LepPDG != 13 ) continue;
    // reco numu CC cut
    bool numuCC = (reco_numu && reco_q == -1);
    if( !numuCC ) continue;

    int cvpimult = gastpc_pi_pl_mult+gastpc_pi_min_mult;
    if( cvpimult > 2 ) cvpimult = 2;
    histCV_gas->Fill( gastpc_pi_pl_mult+gastpc_pi_min_mult, Ev_reco, 1. );
    for( int u = 0; u < nu; ++u ) {

      double shift_Ev_reco = 0.;
      int pimult = 0;
      for( int i = 0; i < nFSP; ++i ) {
        double preco = getP( partEvReco[i], pdg[i] );
        double pshiftfrac = Pscale[u]->Eval(preco);
        double precoshift = preco*(1.+pshiftfrac);

        if( trkLen[i] > trkThreshold[u] && (pdg[i]== 211 || pdg[i]==-211) ) pimult++;
        if( trkLen[i] > trkThreshold[u] ) shift_Ev_reco += getE( precoshift, pdg[i] );
      }

      if( pimult > 2 ) pimult = 2;
      hists_gas[u]->Fill( pimult, shift_Ev_reco );

      val_npi_gas[u]->Fill( cvpimult, pimult );
      val_Ev_gas[u]->Fill( Ev_reco, shift_Ev_reco );

    }
  }

  // FD numu
  cafFDmu->SetBranchAddress( "vtx_x", &vtx_x );
  cafFDmu->SetBranchAddress( "vtx_y", &vtx_y );
  cafFDmu->SetBranchAddress( "vtx_z", &vtx_z );
  cafFDmu->SetBranchAddress( "LepE", &LepE );
  cafFDmu->SetBranchAddress( "Ev", &Ev );
  cafFDmu->SetBranchAddress( "Ev_reco_numu", &Ev_reco );
  cafFDmu->SetBranchAddress( "RecoLepEnNumu", &Elep_reco );
  cafFDmu->SetBranchAddress( "LepPDG", &LepPDG );
  cafFDmu->SetBranchAddress( "cvnnumu", &cvnnumu );
  cafFDmu->SetBranchAddress( "cvnnue", &cvnnue );
  cafFDmu->SetBranchAddress( "eRecoP", &eRecoP );
  cafFDmu->SetBranchAddress( "eRecoN", &eRecoN );
  cafFDmu->SetBranchAddress( "eRecoPip", &eRecoPip );
  cafFDmu->SetBranchAddress( "eRecoPim", &eRecoPim );
  cafFDmu->SetBranchAddress( "eRecoPi0", &eRecoPi0 );
  cafFDmu->SetBranchAddress( "eP", &eP );
  cafFDmu->SetBranchAddress( "eN", &eN );
  cafFDmu->SetBranchAddress( "ePip", &ePip );
  cafFDmu->SetBranchAddress( "ePim", &ePim );
  cafFDmu->SetBranchAddress( "ePi0", &ePi0 );

  cafFDmu->SetBranchStatus( "*", 0 );
  cafFDmu->SetBranchStatus( "vtx_x", 1 );
  cafFDmu->SetBranchStatus( "vtx_y", 1 );
  cafFDmu->SetBranchStatus( "vtx_z", 1 );
  cafFDmu->SetBranchStatus( "LepPDG", 1 );
  cafFDmu->SetBranchStatus( "LepE", 1 );
  cafFDmu->SetBranchStatus( "Ev", 1 );
  cafFDmu->SetBranchStatus( "Ev_reco_numu", 1 );
  cafFDmu->SetBranchStatus( "RecoLepEnNumu", 1 );
  cafFDmu->SetBranchStatus( "cvnnumu", 1 );
  cafFDmu->SetBranchStatus( "cvnnue", 1 );
  cafFDmu->SetBranchStatus( "eRecoP", 1 );
  cafFDmu->SetBranchStatus( "eRecoN", 1 );
  cafFDmu->SetBranchStatus( "eRecoPip", 1 );
  cafFDmu->SetBranchStatus( "eRecoPim", 1 );
  cafFDmu->SetBranchStatus( "eRecoPi0", 1 );
  cafFDmu->SetBranchStatus( "eP", 1 );
  cafFDmu->SetBranchStatus( "eN", 1 );
  cafFDmu->SetBranchStatus( "ePip", 1 );
  cafFDmu->SetBranchStatus( "ePim", 1 );
  cafFDmu->SetBranchStatus( "ePi0", 1 );

  N = cafFDmu->GetEntries();
  for( int ii = 0; ii < N; ++ii ) {
    cafFDmu->GetEntry(ii);

    if( ii % 100000 == 0 ) printf( "FD mu event %d of %d...\n", ii, N );

    //printf( "vtx (%1.1f, %1.1f, %1.1f) LepPDG %d cvn mu %f e %f\n", vtx_x, vtx_y, vtx_z, LepPDG, cvnnumu, cvnnue );

    // FV cut
    if( abs(vtx_x) > 310. || abs(vtx_y) > 550. || vtx_z < 50. || vtx_z > 1244. ) continue; 
    // true numu CC cut
    if( abs(LepPDG) != 13 ) continue;
    // reco numu CC cut
    bool numuCC = (cvnnumu > 0.5 && cvnnue < 0.85);
    if( !numuCC ) continue;

    histCV_FDmu->Fill( Ev_reco, 1. );
    for( int u = 0; u < nu; ++u ) {

      // determine the shifted energies
      double shiftChargedHad = EhadThrowFD[u]->Eval(eRecoP + eRecoPip + eRecoPim);
      double shiftEM = EEMThrowFD[u]->Eval(eRecoPi0);
      double shiftN = EneutThrowFD[u]->Eval(eRecoN);
      double shiftMu = EmuThrowFD[u]->Eval(Elep_reco);
      double shiftTot = EtotThrowFD[u]->Eval(Ev_reco - Elep_reco);

      double Ehad_reco_shift = Ev_reco - Elep_reco;
      Ehad_reco_shift += shiftChargedHad*(eRecoP + eRecoPip + eRecoPim);
      Ehad_reco_shift += shiftEM*eRecoPi0;
      Ehad_reco_shift += shiftN*eRecoN;

      double Elep_reco_shift = Elep_reco*(1.+shiftMu);

      Ehad_reco_shift *= (1.+shiftTot);

      // resolution uncertainties
      Elep_reco_shift += (LepE - Elep_reco)*EmuResFD[u];
      Ehad_reco_shift += ((eP + ePip + ePim) - (eRecoP + eRecoPip + eRecoPim))*EhadResFD[u];
      Ehad_reco_shift += (ePi0 - eRecoPi0)*EEMResFD[u];
      Ehad_reco_shift += (eN - eRecoN)*EneutResFD[u];

      double Ev_reco_shift = Elep_reco_shift + Ehad_reco_shift;

      hists_FDmu[u]->Fill( Ev_reco_shift, 1. );
      //if( u == 7 ) printf( "True %3.3f Reco %3.3f shift %3.3f\n", Ev, Ev_reco, Ev_reco_shift );
    }
  }

  // FD nue
  cafFDe->SetBranchAddress( "vtx_x", &vtx_x );
  cafFDe->SetBranchAddress( "vtx_y", &vtx_y );
  cafFDe->SetBranchAddress( "vtx_z", &vtx_z );
  cafFDe->SetBranchAddress( "LepE", &LepE );
  cafFDe->SetBranchAddress( "Ev", &Ev );
  cafFDe->SetBranchAddress( "Ev_reco_nue", &Ev_reco );
  cafFDe->SetBranchAddress( "RecoLepEnNue", &Elep_reco );
  cafFDe->SetBranchAddress( "LepPDG", &LepPDG );
  cafFDe->SetBranchAddress( "cvnnumu", &cvnnumu );
  cafFDe->SetBranchAddress( "cvnnue", &cvnnue );
  cafFDe->SetBranchAddress( "eRecoP", &eRecoP );
  cafFDe->SetBranchAddress( "eRecoN", &eRecoN );
  cafFDe->SetBranchAddress( "eRecoPip", &eRecoPip );
  cafFDe->SetBranchAddress( "eRecoPim", &eRecoPim );
  cafFDe->SetBranchAddress( "eRecoPi0", &eRecoPi0 );
  cafFDe->SetBranchAddress( "eP", &eP );
  cafFDe->SetBranchAddress( "eN", &eN );
  cafFDe->SetBranchAddress( "ePip", &ePip );
  cafFDe->SetBranchAddress( "ePim", &ePim );
  cafFDe->SetBranchAddress( "ePi0", &ePi0 );

  cafFDe->SetBranchStatus( "*", 0 );
  cafFDe->SetBranchStatus( "vtx_x", 1 );
  cafFDe->SetBranchStatus( "vtx_y", 1 );
  cafFDe->SetBranchStatus( "vtx_z", 1 );
  cafFDe->SetBranchStatus( "LepPDG", 1 );
  cafFDe->SetBranchStatus( "LepE", 1 );
  cafFDe->SetBranchStatus( "Ev", 1 );
  cafFDe->SetBranchStatus( "Ev_reco_nue", 1 );
  cafFDe->SetBranchStatus( "RecoLepEnNue", 1 );
  cafFDe->SetBranchStatus( "cvnnumu", 1 );
  cafFDe->SetBranchStatus( "cvnnue", 1 );
  cafFDe->SetBranchStatus( "eRecoP", 1 );
  cafFDe->SetBranchStatus( "eRecoN", 1 );
  cafFDe->SetBranchStatus( "eRecoPip", 1 );
  cafFDe->SetBranchStatus( "eRecoPim", 1 );
  cafFDe->SetBranchStatus( "eRecoPi0", 1 );
  cafFDe->SetBranchStatus( "eP", 1 );
  cafFDe->SetBranchStatus( "eN", 1 );
  cafFDe->SetBranchStatus( "ePip", 1 );
  cafFDe->SetBranchStatus( "ePim", 1 );
  cafFDe->SetBranchStatus( "ePi0", 1 );

  N = cafFDe->GetEntries();
  for( int ii = 0; ii < N; ++ii ) {
    cafFDe->GetEntry(ii);

    if( ii % 100000 == 0 ) printf( "FD e event %d of %d...\n", ii, N );

    // FV cut
    if( abs(vtx_x) > 310. || abs(vtx_y) > 550. || vtx_z < 50. || vtx_z > 1244. ) continue; 
    // true numu CC cut
    if( abs(LepPDG) != 11 ) continue;
    // reco numu CC cut
    bool nueCC = (cvnnue > 0.85 && cvnnumu < 0.5);
    if( !nueCC ) continue;

    histCV_FDe->Fill( Ev_reco, 1. );
    for( int u = 0; u < nu; ++u ) {

      // determine the shifted energies
      double shiftChargedHad = EhadThrowFD[u]->Eval(eRecoP + eRecoPip + eRecoPim);
      double shiftEM = EEMThrowFD[u]->Eval(eRecoPi0);
      double shiftN = EneutThrowFD[u]->Eval(eRecoN);
      double shiftMu = EmuThrowFD[u]->Eval(Elep_reco);
      double shiftTot = EtotThrowFD[u]->Eval(Ev_reco - Elep_reco);

      double Ehad_reco_shift = Ev_reco - Elep_reco;
      Ehad_reco_shift += shiftChargedHad*(eRecoP + eRecoPip + eRecoPim);
      Ehad_reco_shift += shiftEM*eRecoPi0;
      Ehad_reco_shift += shiftN*eRecoN;

      double Elep_reco_shift = Elep_reco*(1.+shiftMu);

      Ehad_reco_shift *= (1.+shiftTot);

      // resolution uncertainties
      Elep_reco_shift += (LepE - Elep_reco)*EmuResFD[u];
      Ehad_reco_shift += ((eP + ePip + ePim) - (eRecoP + eRecoPip + eRecoPim))*EhadResFD[u];
      Ehad_reco_shift += (ePi0 - eRecoPi0)*EEMResFD[u];
      Ehad_reco_shift += (eN - eRecoN)*EneutResFD[u];

      double Ev_reco_shift = Elep_reco_shift + Ehad_reco_shift;

      hists_FDe[u]->Fill( Ev_reco_shift, 1. );
    }
  }

  // Now determine the actual covariance
  int n_bins = n_Ebins * n_ybins;

  TMatrixD cov( n_bins, n_bins );
  TMatrixD covAcc( n_bins, n_bins );
  TMatrixD covEscale( n_bins, n_bins );
  for( int b0 = 0; b0 < n_bins; ++b0 ) {
    for( int b1 = 0; b1 < n_bins; ++b1 ) {
      cov[b0][b1] = 0.;
      covAcc[b0][b1] = 0.;
      covEscale[b0][b1] = 0.;
    }
  }
  for( int b0 = 0; b0 < n_bins; ++b0 ) {
    int bE0, by0;
    get2Dbins( b0+1, bE0, by0 );
    double cv0 = histCV->GetBinContent( bE0, by0 );

    for( int b1 = 0; b1 < n_bins; ++b1 ) {
      int bE1, by1;
      get2Dbins( b1+1, bE1, by1 );
      double cv1 = histCV->GetBinContent( bE1, by1 );

      for( int u = 0; u < nu; ++u ) {
        double u0 = hists[u]->GetBinContent( bE0, by0 );
        double u1 = hists[u]->GetBinContent( bE1, by1 );

        double u0a = histsAccOnly[u]->GetBinContent( bE0, by0 );
        double u1a = histsAccOnly[u]->GetBinContent( bE1, by1 );

        double u0e = histsEscaleOnly[u]->GetBinContent( bE0, by0 );
        double u1e = histsEscaleOnly[u]->GetBinContent( bE1, by1 );

        // fractional covariance, dividing out number of universes at the same time
        if( cv0*cv1 > 1.E-12 ) {
          cov[b0][b1] += (u0-cv0)*(u1-cv1)/(cv0*cv1*nu);
          covAcc[b0][b1] += (u0a-cv0)*(u1a-cv1)/(cv0*cv1*nu);
          covEscale[b0][b1] += (u0e-cv0)*(u1e-cv1)/(cv0*cv1*nu);
        }
      }
    }
  }

  // matrices are not positive definite due to numerical precision; make them positive definite
  fix( cov );
  fix( covAcc );
  fix( covEscale );

  // test that they are invertible, this will barf if they are singular
  TMatrixD covInv( TMatrixD::kInverted, cov );
  TMatrixD covAccInv( TMatrixD::kInverted, covAcc );
  TMatrixD covEscaleInv( TMatrixD::kInverted, covEscale );

  // Gas TPC
  int n_binsgas = n_Ebins * 3;

  TMatrixD covGas( n_binsgas, n_binsgas );
  for( int b0 = 0; b0 < n_binsgas; ++b0 ) {
    for( int b1 = 0; b1 < n_binsgas; ++b1 ) {
      covGas[b0][b1] = 0.;
    }
  }

  for( int b0 = 0; b0 < n_binsgas; ++b0 ) {
    int bx0 = (b0 % 3) + 1;
    int by0 = (b0 / 3) + 1;
    
    double cv0 = histCV_gas->GetBinContent( bx0, by0 );

    for( int b1 = 0; b1 < n_binsgas; ++b1 ) {
      int bx1 = (b1 % 3) + 1;
      int by1 = (b1 / 3) + 1;
      double cv1 = histCV_gas->GetBinContent( bx1, by1 );

      for( int u = 0; u < nu; ++u ) {
        double u0 = hists_gas[u]->GetBinContent( bx0, by0 );
        double u1 = hists_gas[u]->GetBinContent( bx1, by1 );

        // fractional covariance, dividing out number of universes at the same time
        if( cv0*cv1 > 1.E-12 ) {
          covGas[b0][b1] += (u0-cv0)*(u1-cv1)/(cv0*cv1*nu);
        }
      }
    }
  }

  // matrices are not positive definite due to numerical precision; make them positive definite
  fix( covGas );

  // test that they are invertible, this will barf if they are singular
  TMatrixD covGasInv( TMatrixD::kInverted, covGas );


  // FD matrices
  // Now determine the actual covariance
  n_bins = n_Ebins;

  TMatrixD covMu( n_bins, n_bins );
  TMatrixD covE( n_bins, n_bins );
  for( int b0 = 0; b0 < n_bins; ++b0 ) {
    for( int b1 = 0; b1 < n_bins; ++b1 ) {
      covMu[b0][b1] = 0.;
      covE[b0][b1] = 0.;
    }
  }
  for( int b0 = 0; b0 < n_bins; ++b0 ) {
    double cv0mu = histCV_FDmu->GetBinContent( b0+1 );
    double cv0e = histCV_FDe->GetBinContent( b0+1 );

    for( int b1 = 0; b1 < n_bins; ++b1 ) {
      double cv1mu = histCV_FDmu->GetBinContent( b1+1 );
      double cv1e = histCV_FDe->GetBinContent( b1+1 );

      for( int u = 0; u < nu; ++u ) {
        double u0mu = hists_FDmu[u]->GetBinContent( b0+1 );
        double u1mu = hists_FDmu[u]->GetBinContent( b1+1 );

        double u0e = hists_FDe[u]->GetBinContent( b0+1 );
        double u1e = hists_FDe[u]->GetBinContent( b1+1 );

        // fractional covariance, dividing out number of universes at the same time
        if( cv0mu*cv1mu > 1.E-12 ) {
          covMu[b0][b1] += (u0mu-cv0mu)*(u1mu-cv1mu)/(cv0mu*cv1mu*nu);
        }
        if( cv0e*cv1e > 1.E-12 ) {
          covE[b0][b1] += (u0e-cv0e)*(u1e-cv1e)/(cv0e*cv1e*nu);
        }
      }
    }
  }

  // matrices are not positive definite due to numerical precision; make them positive definite
  fix( covMu );
  fix( covE );

  // test that they are invertible, this will barf if they are singular
  TMatrixD covMuInv( TMatrixD::kInverted, covMu );
  TMatrixD covEInv( TMatrixD::kInverted, covE );

  TFile * outfile = new TFile( "ND_syst_cov.root", "RECREATE" );
  cov.Write("nd_frac_cov");
  covAcc.Write("nd_frac_cov_accOnly");
  covEscale.Write("nd_frac_cov_EscaleOnly");
  covMu.Write("fd_numu_frac_cov");
  covE.Write("fd_nue_frac_cov");
  covGas.Write("nd_frac_cov_gasTPC");

  TCanvas * c = new TCanvas();
  covAcc.Draw("colz");
  c->Print( "ND_syst_cov_acc.png" );
  covEscale.Draw("colz");
  c->Print( "ND_syst_cov_Escale.png" );
  cov.Draw("colz");
  c->Print( "ND_syst_cov.png" );
  c->SetLogz(1);
  c->Print( "ND_syst_cov_log.png" );    
  c->SetLogz(0);


  // Energy projection, for display purposes only
  TH2D * covE_acc = new TH2D( "covE_acc", ";Neutrino energy (GeV);Neutrino energy (GeV)", n_Ebins, Ebins, n_Ebins, Ebins );
  TH2D * covE_scale = new TH2D( "covE_scale", ";Neutrino energy (GeV);Neutrino energy (GeV)", n_Ebins, Ebins, n_Ebins, Ebins );

  for( int b0 = 1; b0 <= n_Ebins; ++b0 ) {
    double cv0 = histCV->ProjectionX()->GetBinContent( b0 );
    for( int b1 = 1; b1 <= n_Ebins; ++b1 ) {
      double cv1 = histCV->ProjectionX()->GetBinContent( b1 );

      for( int u = 0; u < nu; ++u ) {
        double u0a = histsAccOnly[u]->ProjectionX()->GetBinContent( b0 );
        double u1a = histsAccOnly[u]->ProjectionX()->GetBinContent( b1 );

        double u0e = histsEscaleOnly[u]->ProjectionX()->GetBinContent( b0 );
        double u1e = histsEscaleOnly[u]->ProjectionX()->GetBinContent( b1 );

        // fractional covariance, dividing out number of universes at the same time
        if( cv0*cv1 > 1.E-12 ) {
          covE_acc->Fill( covE_acc->GetXaxis()->GetBinCenter(b0), covE_acc->GetXaxis()->GetBinCenter(b1), (u0a-cv0)*(u1a-cv1)/(cv0*cv1*nu) );
          covE_scale->Fill( covE_acc->GetXaxis()->GetBinCenter(b0), covE_acc->GetXaxis()->GetBinCenter(b1), (u0e-cv0)*(u1e-cv1)/(cv0*cv1*nu) );
        }
      }
    }
  }

  covE_acc->SetMaximum(0.0005);
  covE_acc->Draw("colz");
  c->Print( "ND_syst_cov_projE_acc.png" );
  covE_scale->SetMaximum(0.025);
  covE_scale->Draw("colz");
  c->Print( "ND_syst_cov_projE_scale.png" );


  TFile * val = new TFile( "out.root", "RECREATE" );
  histCV_gas->Write();
  for( int u = 0; u < nu; ++u ) {
    hists_gas[u]->Write();
    val_npi_gas[u]->Write();
    val_Ev_gas[u]->Write();
  }
  covGas->Write( "gas_cov" );
}
















