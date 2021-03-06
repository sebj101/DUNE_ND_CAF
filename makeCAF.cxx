#include "CAF.C"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>
#include <vector>

TRandom3 * rando;
const double mmu = 0.1056583745;
const double e   = 2.7182818285;
TF1 * tsmear; // angular resolution function
TSpline3 * sp; // pion interaction in ECAL probability

// params will be extracted from command line, and passed to the reconstruction
struct params {
  double OA_xcoord;
  bool fhc, grid, IsGasTPC;
  int seed, run, subrun, first, n, nfiles;
  double trk_muRes, LAr_muRes, ECAL_muRes;
  double em_const, em_sqrtE;
  double michelEff;
  double CC_trk_length;
  double pileup_frac, pileup_max;
  double gastpc_len, gastpc_B, gastpc_padPitch, gastpc_X0;
};

// Fill reco variables for muon reconstructed in magnetized tracker
void recoMuonTracker( CAF &caf, params &par )
{
  // smear momentum by resolution
  double p = sqrt(caf.LepE*caf.LepE - mmu*mmu);
  double reco_p = rando->Gaus( p, p*par.trk_muRes );
  caf.Elep_reco = sqrt(reco_p*reco_p + mmu*mmu);

  double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
  double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
  double evalTsmear = tsmear->Eval(caf.Elep_reco - mmu);
  if( evalTsmear < 0. ) evalTsmear = 0.;
  double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
  double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
  caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

  // assume perfect charge reconstruction
  caf.reco_q = (caf.LepPDG > 0 ? -1 : 1);

  // assume always muon for tracker-matched
  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 0; caf.muon_tracker = 1; caf.muon_ecal = 0; caf.muon_exit = 0;
  caf.Ev_reco = caf.Elep_reco;
}

// Fill reco muon variables for muon contained in LAr
void recoMuonLAr( CAF &caf, params &par )
{
  // range-based, smear kinetic energy
  double ke = caf.LepE - mmu;
  double reco_ke = rando->Gaus( ke, ke*par.LAr_muRes );
  caf.Elep_reco = reco_ke + mmu;

  double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
  double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
  double evalTsmear = tsmear->Eval(caf.Elep_reco - mmu);
  if( evalTsmear < 0. ) evalTsmear = 0.;
  double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
  double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
  caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );


  // assume negative for FHC, require Michel for RHC
  if( par.fhc ) caf.reco_q = -1;
  else {
    double michel = rando->Rndm();
    if( caf.LepPDG == -13 && michel < par.michelEff ) caf.reco_q = 1; // correct mu+
    else if( caf.LepPDG == 13 && michel < par.michelEff*0.25 ) caf.reco_q = 1; // incorrect mu-
    else caf.reco_q = -1; // no reco Michel
  }

  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 1; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
  caf.Ev_reco = caf.Elep_reco;
}

// Fill reco variables for muon reconstructed in magnetized tracker
void recoMuonECAL( CAF &caf, params &par )
{
  // range-based KE
  double ke = caf.LepE - mmu;
  double reco_ke = rando->Gaus( ke, ke*par.ECAL_muRes );
  caf.Elep_reco = reco_ke + mmu;

  double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
  double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
  double evalTsmear = tsmear->Eval(caf.Elep_reco - mmu);
  if( evalTsmear < 0. ) evalTsmear = 0.;
  double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
  double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
  caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

  // assume perfect charge reconstruction -- these are fairly soft and should curve a lot in short distance
  caf.reco_q = (caf.LepPDG > 0 ? -1 : 1);

  // assume always muon for ecal-matched
  caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
  caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 1; caf.muon_exit = 0;
  caf.Ev_reco = caf.Elep_reco;
}

// Fill reco variables for true electron
void recoElectron( CAF &caf, params &par )
{
  caf.reco_q = 0; // never know charge
  caf.reco_numu = 0;
  caf.muon_contained = 0; caf.muon_tracker = 1; caf.muon_ecal = 0; caf.muon_exit = 0;

  // fake efficiency...threshold of 300 MeV, eff rising to 100% by 700 MeV
  if( rando->Rndm() > (caf.LepE-0.3)*2.5 ) { // reco as NC
    caf.Elep_reco = 0.;
    caf.reco_nue = 0; caf.reco_nc = 1;
    caf.Ev_reco = caf.LepE; // include electron energy in Ev anyway, since it won't show up in reco hadronic energy
  } else { // reco as CC
    caf.Elep_reco = rando->Gaus( caf.LepE, caf.LepE*(par.em_const + par.em_sqrtE/sqrt(caf.LepE)) );
    caf.reco_nue = 1; caf.reco_nc = 0;
    caf.Ev_reco = caf.Elep_reco;
  }

  double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
  double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
  double evalTsmear = 3. + tsmear->Eval(caf.Elep_reco - mmu);
  if( evalTsmear < 0. ) evalTsmear = 0.;
  double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
  double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
  caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

}

void decayPi0( TLorentzVector pi0, TVector3 &gamma1, TVector3 &gamma2 )
{
  double e = pi0.E();
  double mp = 134.9766; // pi0 mass

  double beta = sqrt( 1. - (mp*mp)/(e*e) ); // velocity of pi0
  double theta = 3.1416*rando->Rndm(); // theta of gamma1 w.r.t. pi0 direction
  double phi = 2.*3.1416*rando->Rndm(); // phi of gamma1 w.r.t. pi0 direction

  double p = mp/2.; // photon momentum in pi0 rest frame
  TLorentzVector g1( 0., 0., p, p ); // pre-rotation photon 1
  TLorentzVector g2( 0., 0., -p, p ); // pre-rotation photon 2 is opposite

  // rotate to the random decay axis in pi0 rest frame. choice of rotation about x instead of y is arbitrary
  g1.RotateX( theta );
  g2.RotateX( theta );
  g1.RotateZ( phi );
  g2.RotateZ( phi );

  // boost to lab frame with pi0 velocity. pi0 direction is z axis for this
  g1.Boost( 0., 0., beta );
  g2.Boost( 0., 0., beta );

  // make gamma1 the more energetic one
  if( g1.E() > g2.E() ) {
    gamma1 = g1.Vect();
    gamma2 = g2.Vect();
  } else {
    gamma1 = g2.Vect();
    gamma2 = g1.Vect();
  }

  // rotate from frame where pi0 is z' direction into neutrino frame
  TVector3 pi0dir = pi0.Vect().Unit(); // actually w.r.t. neutrino direction
  gamma1.RotateUz( pi0dir );
  gamma2.RotateUz( pi0dir );
}

// main loop function
void loop( CAF &caf, params &par, TTree * tree, std::string ghepdir, std::string fhicl_filename )
{
  // read in edep-sim output file
  int ifileNo, ievt, lepPdg, muonReco, nFS;
  float lepKE, muGArLen, hadTot, hadCollar;
  float hadP, hadN, hadPip, hadPim, hadPi0, hadOther;
  float p3lep[3], vtx[3], muonExitPt[3], muonExitMom[3];
  int fsPdg[100];
  float fsPx[100], fsPy[100], fsPz[100], fsE[100], fsTrkLen[100], fsTrkLenPerp[100];
  tree->SetBranchAddress( "ifileNo", &ifileNo );
  tree->SetBranchAddress( "ievt", &ievt );
  tree->SetBranchAddress( "lepPdg", &lepPdg );
  tree->SetBranchAddress( "muonReco", &muonReco );
  tree->SetBranchAddress( "lepKE", &lepKE );
  tree->SetBranchAddress( "muGArLen", &muGArLen );
  tree->SetBranchAddress( "hadTot", &hadTot );
  tree->SetBranchAddress( "hadCollar", &hadCollar );
  tree->SetBranchAddress( "hadP", &hadP );
  tree->SetBranchAddress( "hadN", &hadN );
  tree->SetBranchAddress( "hadPip", &hadPip );
  tree->SetBranchAddress( "hadPim", &hadPim );
  tree->SetBranchAddress( "hadPi0", &hadPi0 );
  tree->SetBranchAddress( "hadOther", &hadOther );
  tree->SetBranchAddress( "p3lep", p3lep );
  tree->SetBranchAddress( "vtx", vtx );
  tree->SetBranchAddress( "muonExitPt", muonExitPt );
  tree->SetBranchAddress( "muonExitMom", muonExitMom );
  tree->SetBranchAddress( "nFS", &nFS );
  tree->SetBranchAddress( "fsPdg", fsPdg );
  tree->SetBranchAddress( "fsPx", fsPx );
  tree->SetBranchAddress( "fsPy", fsPy );
  tree->SetBranchAddress( "fsPz", fsPz );
  tree->SetBranchAddress( "fsE", fsE );
  tree->SetBranchAddress( "fsTrkLen", fsTrkLen );
  tree->SetBranchAddress( "fsTrkLenPerp", fsTrkLenPerp );

  // Get GHEP file for genie::EventRecord from other file
  int current_file = -1;
  TFile * ghep_file = NULL;
  TTree * gtree = NULL;

  std::string mode = ( par.fhc ? "neutrino" : "antineutrino" );

  // DUNE reweight getter
  nusyst::response_helper rh( fhicl_filename );
  // Get list of variations, and make CAF branch for each one
  std::vector<unsigned int> parIds = rh.GetParameters();
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations.size() );
    bool is_wgt = head.isWeightSystematicVariation;
    std::string wgt_var = ( is_wgt ? "wgt" : "var" );
    caf.addRWbranch( parIds[i], head.prettyName, wgt_var, head.paramVariations );
    caf.iswgt[parIds[i]] = is_wgt;
  }

  caf.pot = 0.;

  // Main event loop
  int N = tree->GetEntries();
  if( par.n > 0 && par.n < N ) N = par.n + par.first;
  for( int ii = par.first; ii < N; ++ii ) {

    tree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d of %d...\n", ii, N );

    caf.setToBS();

    // set defaults
    for( int j = 0; j < 100; ++j ) {
      caf.nwgt[j] = 7;
      caf.cvwgt[j] = ( caf.iswgt[j] ? 1. : 0. );
      for( unsigned int k = 0; k < 100; ++k ) {
        caf.wgt[j][k] = ( caf.iswgt[j] ? 1. : 0. );
      }
    }

    // make sure ghep file matches the current one, otherwise update to the current ghep file
    if( ifileNo != current_file ) {
      // close the previous file
      if( ghep_file ) ghep_file->Close();

      if( par.grid ) ghep_file = new TFile( Form("genie.%d.root", ifileNo) );
      else if( !par.IsGasTPC ) ghep_file = new TFile( Form("%s/%02d/LAr.%s.%d.ghep.root", ghepdir.c_str(), ifileNo/1000, mode.c_str(), ifileNo) );
      else ghep_file = new TFile( Form("%s/%02d/GAr.%s.%d.ghep.root", ghepdir.c_str(), ifileNo/1000, mode.c_str(), ifileNo) );
      
      gtree = (TTree*) ghep_file->Get( "gtree" );
      caf.pot += gtree->GetWeight();
      printf( "New GHEP file with %g POT, total = %g\n", gtree->GetWeight(), caf.pot );

      // can't find GHepRecord
      if( gtree == NULL ) {
        printf( "Can't find ghep event record for file %d!!!\n", ifileNo );
        continue;
      }

      gtree->SetBranchAddress( "gmcrec", &caf.mcrec );
      current_file = ifileNo;
    }

    caf.vtx_x = vtx[0];
    caf.vtx_y = vtx[1];
    caf.vtx_z = vtx[2]; 
    caf.det_x = par.OA_xcoord;

    // configuration variables in CAF file; we don't use mvaresult so just set it to zero
    caf.run = par.run;
    caf.subrun = par.subrun;
    caf.event = ii;
    caf.isFD = 0;
    caf.isFHC = par.fhc;

    // get GENIE event record
    gtree->GetEntry( ievt );
    genie::EventRecord * event = caf.mcrec->event;
    genie::Interaction * in = event->Summary();

    // Get truth stuff out of GENIE ghep record
    caf.neutrinoPDG = in->InitState().ProbePdg();
    caf.neutrinoPDGunosc = in->InitState().ProbePdg(); // fill this for similarity with FD, but no oscillations
    caf.mode = in->ProcInfo().ScatteringTypeId();
    caf.Ev = in->InitState().ProbeE(genie::kRfLab);
    caf.LepPDG = in->FSPrimLeptonPdg();
    caf.isCC = (abs(caf.LepPDG) == 13 || abs(caf.LepPDG) == 11);
    
    TLorentzVector lepP4;
    TLorentzVector nuP4nuc = *(in->InitState().GetProbeP4(genie::kRfHitNucRest));
    TLorentzVector nuP4 = *(in->InitState().GetProbeP4(genie::kRfLab));

    caf.nP = 0;
    caf.nN = 0;
    caf.nipip = 0;
    caf.nipim = 0;
    caf.nipi0 = 0;
    caf.nikp = 0;
    caf.nikm = 0;
    caf.nik0 = 0;
    caf.niem = 0;
    caf.niother = 0;
    caf.nNucleus = 0;
    caf.nUNKNOWN = 0; // there is an "other" category so this never gets used
    caf.eP = 0.;
    caf.eN = 0.;
    caf.ePip = 0.;
    caf.ePim = 0.;
    caf.ePi0 = 0.;
    caf.eOther = 0.;
    caf.eRecoP = 0.;
    caf.eRecoN = 0.;
    caf.eRecoPip = 0.;
    caf.eRecoPim = 0.;
    caf.eRecoPi0 = 0.;
    caf.eOther = 0.;
    for( int i = 0; i < nFS; ++i ) {
      double ke = 0.001*(fsE[i] - sqrt(fsE[i]*fsE[i] - fsPx[i]*fsPx[i] - fsPy[i]*fsPy[i] - fsPz[i]*fsPz[i]));
      if( fsPdg[i] == caf.LepPDG ) {
        lepP4.SetPxPyPzE( fsPx[i]*0.001, fsPy[i]*0.001, fsPz[i]*0.001, fsE[i]*0.001 );
        caf.LepE = fsE[i]*0.001;
      }
      else if( fsPdg[i] == 2212 ) {caf.nP++; caf.eP += ke;}
      else if( fsPdg[i] == 2112 ) {caf.nN++; caf.eN += ke;}
      else if( fsPdg[i] ==  211 ) {caf.nipip++; caf.ePip += ke;}
      else if( fsPdg[i] == -211 ) {caf.nipim++; caf.ePim += ke;}
      else if( fsPdg[i] ==  111 ) {caf.nipi0++; caf.ePi0 += ke;}
      else if( fsPdg[i] ==  321 ) {caf.nikp++; caf.eOther += ke;}
      else if( fsPdg[i] == -321 ) {caf.nikm++; caf.eOther += ke;}
      else if( fsPdg[i] == 311 || fsPdg[i] == -311 || fsPdg[i] == 130 || fsPdg[i] == 310 ) {caf.nik0++; caf.eOther += ke;}
      else if( fsPdg[i] ==   22 ) {caf.niem++; caf.eOther += ke;}
      else if( fsPdg[i] > 1000000000 ) caf.nNucleus++;
      else {caf.niother++; caf.eOther += ke;}
    }

    // true 4-momentum transfer
    TLorentzVector q = nuP4-lepP4;

    // Q2, W, x, y frequently do not get filled in GENIE Kinematics object, so calculate manually
    caf.Q2 = -q.Mag2();
    caf.W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2()); // "Wexp"
    caf.X = -q.Mag2()/(2*0.939*q.E());
    caf.Y = q.E()/caf.Ev;

    caf.theta_reco = -1.; // default value

    caf.NuMomX = nuP4.X();
    caf.NuMomY = nuP4.Y();
    caf.NuMomZ = nuP4.Z();
    caf.LepMomX = lepP4.X();
    caf.LepMomY = lepP4.Y();
    caf.LepMomZ = lepP4.Z();
    caf.LepE = lepP4.E();
    caf.LepNuAngle = nuP4.Angle( lepP4.Vect() );

    // Add DUNErw weights to the CAF

    systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
    for( systtools::event_unit_response_w_cv_t::iterator it = resp.begin(); it != resp.end(); ++it ) {
      caf.nwgt[(*it).pid] = (*it).responses.size();
      caf.cvwgt[(*it).pid] = (*it).CV_response;
      for( unsigned int i = 0; i < (*it).responses.size(); ++i ) {
        caf.wgt[(*it).pid][i] = (*it).responses[i];
      }
    }

    //--------------------------------------------------------------------------
    // Parameterized reconstruction
    //--------------------------------------------------------------------------
    if( !par.IsGasTPC ) {
      // Loop over final-state particles
      double longest_mip = 0.;
      double longest_mip_KE = 0.;
      int longest_mip_charge = 0;
      caf.reco_lepton_pdg = 0;
      int electrons = 0;
      double electron_energy = 0.;
      int reco_electron_pdg = 0;
      for( int i = 0; i < nFS; ++i ) {
        int pdg = fsPdg[i];
        double p = sqrt(fsPx[i]*fsPx[i] + fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]);
        double KE = fsE[i] - sqrt(fsE[i]*fsE[i] - p*p);

        if( (abs(pdg) == 13 || abs(pdg) == 211) && fsTrkLen[i] > longest_mip ) {
          longest_mip = fsTrkLen[i];
          longest_mip_KE = KE;
          caf.reco_lepton_pdg = pdg;
          if( pdg == 13 || pdg == -211 ) longest_mip_charge = -1;
          else longest_mip_charge = 1;
        }

        // pi0 as nu_e
        if( pdg == 111 ) {
          TVector3 g1, g2;
          TLorentzVector pi0( fsPx[i], fsPy[i], fsPz[i], fsE[i] );
          decayPi0( pi0, g1, g2 );
          double g1conv = rando->Exp( 14. ); // conversion distance
          bool compton = (rando->Rndm() < 0.15); // dE/dX misID probability for photon
          // if energetic gamma converts in first wire, and other gamma is either too soft or too colinear
          if( g1conv < 2.0 && compton && (g2.Mag() < 50. || g1.Angle(g2) < 0.01) ) electrons++;
          electron_energy = g1.Mag();
          reco_electron_pdg = 111;
        }
      }

      // True CC reconstruction
      if( abs(lepPdg) == 11 ) { // true nu_e
        recoElectron( caf, par );
        electrons++;
        reco_electron_pdg = lepPdg;
      } else if( abs(lepPdg) == 13 ) { // true nu_mu
        if     ( muonReco == 2 ) recoMuonTracker( caf, par ); // gas TPC match
        else if( muonReco == 1 ) recoMuonLAr( caf, par ); // LAr-contained muon, this might get updated to NC...
        else if( muonReco == 3 ) recoMuonECAL( caf, par ); // ECAL-stopper
        else { // exiting but poorly-reconstructed muon
          caf.Elep_reco = longest_mip * 0.0022;
          caf.reco_q = 0;
          caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
          caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 1;

          double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
          double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
          double evalTsmear = tsmear->Eval(caf.Elep_reco - mmu);
          if( evalTsmear < 0. ) evalTsmear = 0.;
          double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
          double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
          caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
        }
      } else { // NC -- set PID variables, will get updated later if fake CC
        caf.Elep_reco = 0.;
        caf.reco_q = 0;
        caf.reco_numu = 0; caf.reco_nue = 0; caf.reco_nc = 1;
        caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
      }

      // CC/NC confusion
      if( electrons == 1 && muonReco <= 1 ) { // NC or numuCC reco as nueCC
        caf.Elep_reco = electron_energy*0.001;
        caf.reco_q = 0;
        caf.reco_numu = 0; caf.reco_nue = 1; caf.reco_nc = 0;
        caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
        caf.reco_lepton_pdg = reco_electron_pdg;
      } else if( muonReco <= 1 && !(abs(lepPdg) == 11 && caf.Elep_reco > 0.) && (longest_mip < par.CC_trk_length || longest_mip_KE/longest_mip > 3.) ) { 
        // reco as NC
        caf.Elep_reco = 0.;
        caf.reco_q = 0;
        caf.reco_numu = 0; caf.reco_nue = 0; caf.reco_nc = 1;
        caf.muon_contained = 0; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
        caf.reco_lepton_pdg = 0;
      } else if( (abs(lepPdg) == 12 || abs(lepPdg) == 14) && longest_mip > par.CC_trk_length && longest_mip_KE/longest_mip < 3. ) { // true NC reco as CC numu
        caf.Elep_reco = longest_mip_KE*0.001 + mmu;
        if( par.fhc ) caf.reco_q = -1;
        else {
          double michel = rando->Rndm();
          if( longest_mip_charge == 1 && michel < par.michelEff ) caf.reco_q = 1; // correct mu+
          else if( michel < par.michelEff*0.25 ) caf.reco_q = 1; // incorrect mu-
          else caf.reco_q = -1; // no reco Michel
        }
        caf.reco_numu = 1; caf.reco_nue = 0; caf.reco_nc = 0;
        caf.muon_contained = 1; caf.muon_tracker = 0; caf.muon_ecal = 0; caf.muon_exit = 0;
      }

      // Hadronic energy calorimetrically
      caf.Ev_reco = caf.Elep_reco + hadTot*0.001;
      caf.Ehad_veto = hadCollar;
      caf.eRecoP = hadP*0.001;
      caf.eRecoN = hadN*0.001;
      caf.eRecoPip = hadPip*0.001;
      caf.eRecoPim = hadPim*0.001;
      caf.eRecoPi0 = hadPi0*0.001;
      caf.eRecoOther = hadOther*0.001;

      caf.pileup_energy = 0;
      if( rando->Rndm() < par.pileup_frac ) caf.pileup_energy = rando->Rndm() * par.pileup_max;
      caf.Ev_reco += caf.pileup_energy;
    } else {
      std::vector<double> piRecoE;
      // gas TPC: FS particle loop look for long enough tracks and smear momenta
      caf.Ev_reco = 0.;
      caf.nFSP = nFS;
      // Need to determine which of the tracks is reco'd as the muon
      double highestMipPdg, highestMipMom = 0.;

      for( int i = 0; i < nFS; ++i ) {
        double ptrue = 0.001*sqrt(fsPx[i]*fsPx[i] + fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]);
        double mass = 0.001*sqrt(fsE[i]*fsE[i] - fsPx[i]*fsPx[i] - fsPy[i]*fsPy[i] - fsPz[i]*fsPz[i]);
        caf.pdg[i] = fsPdg[i];
        caf.ptrue[i] = ptrue;
        caf.trkLen[i] = fsTrkLen[i];
        caf.trkLenPerp[i] = fsTrkLenPerp[i];
	double pT = 0.001*sqrt(fsPy[i]*fsPy[i] + fsPz[i]*fsPz[i]); // transverse to B field, in GeV

        // track length cut 6cm according to T Junk
        if( fsTrkLen[i] > 0. && fsPdg[i] != 2112) {
	  // basically select charged particles; somehow neutrons ocasionally get nonzero track length
	  double nHits = fsTrkLen[i] / par.gastpc_padPitch; // doesn't matter if not integer as only used in eq
	  // Gluckstern formula, sigmapT/pT, with sigmaX and L in meters
	  double fracSig_meas = sqrt(720./(nHits+4)) * (0.01*par.gastpc_padPitch/sqrt(12.)) * pT / (0.3 * par.gastpc_B * 0.0001 * fsTrkLenPerp[i]*fsTrkLenPerp[i]);
	  // multiple scattering term
	  double fracSig_MCS = 0.052 / (par.gastpc_B * sqrt(par.gastpc_X0*fsTrkLenPerp[i]*0.0001));

	  double sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
	  double preco = rando->Gaus( ptrue, sigmaP );
	  double ereco = sqrt( preco*preco + mass*mass ) - mass; // kinetic energy

	  if(fsTrkLen[i] > par.gastpc_len && abs(fsPdg[i]) != 211 && abs(fsPdg[i]) != 13) {
	    caf.gastpc_nRecoFS++;
	    if (fsPdg[i] == 2212 && preco > 1.5) {
	      // Energy reconstruction for particles in the ECAL
	      double keTrue = sqrt(ptrue*ptrue + mass*mass) - mass;
	      double keCalReco = keTrue * (1.+rando->Gaus(0., 0.2)+rando->Gaus(0., 0.3)/sqrt(keTrue));
	      // Check if energy from momentum reco is closer to proton or pion expectation
	      double eIfPi  = sqrt(preco*preco + 0.1395*0.1395) - 0.1395;
	      double eIfPro = sqrt(preco*preco + 0.9383*0.9383) - 0.9383;
	      // Is reco'd as a pion
	      if (abs(eIfPi-keCalReco) < abs(eIfPro-keCalReco)) {
		// Add in pion mass
		ereco += 0.1395;
		piRecoE.push_back(ereco);
		caf.gastpc_pi_pl_mult++;
	      }
	      else {
		caf.gastpc_pro_mult++;
	      }
	    }
	    else if (fsPdg[i] == 2212 && preco < 1.5) {
	      // For lower momentum protons, we can get PID from de/dx
	      caf.gastpc_pro_mult++;
	    }
	    else { // Unidentified hadron
	      caf.gastpc_other_had_mult++;
	    }
	    caf.partEvReco[i] = ereco;
	    caf.Ev_reco += ereco;
	  }
	  else if (fsTrkLen[i] > par.gastpc_len && (abs(fsPdg[i]) == 211 || abs(fsPdg[i]) == 13)) {
	    caf.gastpc_nRecoFS++;
	    // Is a MIP
	    caf.Ev_reco += (ereco+mass);
	    if (abs(fsPdg[i])==13) {
	      // If highest momentum MIP is muon then we're all good
	      if (highestMipMom > 0.) {
		piRecoE.push_back(sqrt(highestMipMom*highestMipMom+0.1395*0.1395));
		highestMipPdg > 0 ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
	      }
	      highestMipMom = preco;
	      highestMipPdg = fsPdg[i];
	      caf.Elep_reco = sqrt(preco*preco + mass*mass);
	      // angle reconstruction
	      double true_tx = 1000.*atan(caf.LepMomX / caf.LepMomZ);
	      double true_ty = 1000.*atan(caf.LepMomY / caf.LepMomZ);
	      double evalTsmear = tsmear->Eval(caf.Elep_reco - mmu);
	      if( evalTsmear < 0. ) evalTsmear = 0.;
	      double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
	      double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
	      caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
	      caf.gastpc_RecoLepMomZ = preco / (sqrt(tan(reco_tx*0.001)*tan(reco_tx*0.001)+tan(reco_ty*0.001)*tan(reco_ty*0.001)+1));
	      caf.gastpc_RecoLepMomX = caf.gastpc_RecoLepMomZ * tan(reco_tx*0.001); 
	      caf.gastpc_RecoLepMomY = caf.gastpc_RecoLepMomZ * tan(reco_ty*0.001);
	      // To reco as numu CC require track length perpendicular to B field > 1m 
	      // or muon exits and ranges out in ECAL
	      bool muonExits = (muonExitPt[2] > 0.) ? true : false;
	      double muonMom = sqrt( muonExitMom[0]*muonExitMom[0] + muonExitMom[1]*muonExitMom[1] + muonExitMom[2]*muonExitMom[2] );
	      if (fsTrkLenPerp[i] > 100. || (muonExits && muonMom < 380.)) { // Reco as numu CC
		// assume perfect charge reconstruction
		caf.reco_q = (fsPdg[i] > 0 ? -1 : 1);
		caf.reco_numu=1; caf.reco_nue=0; caf.reco_nc=0;
		caf.muon_tracker = 1;
	      }
	    }
	    else if (abs(fsPdg[i])==211 && preco < 0.4) {
	      // We have good PID here
	      piRecoE.push_back(ereco+0.1395);
	      fsPdg[i] > 0 ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
	    }
	    else if (abs(fsPdg[i])==211 && preco > 0.4 && preco < 1.) {
	      // Use a TSpline3 to check probability of this pion interacting in ECAL
	      if (preco > highestMipMom && rando->Rndm() < sp->Eval(ptrue*1000.)) {
		// We now think this pion is the lepton
		if (highestMipMom > 0.) {
		  piRecoE.push_back(sqrt(highestMipMom*highestMipMom+0.1395*0.1395));
		  (highestMipPdg==-13 || highestMipPdg==211) ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
		}
		highestMipMom = preco;
		highestMipPdg = fsPdg[i];
		caf.Elep_reco = sqrt(preco*preco + mmu*mmu);
		// angle reconstruction
		double true_tx = 1000.*atan(fsPx[i] / fsPz[i]);
		double true_ty = 1000.*atan(fsPy[i] / fsPz[i]);
		double evalTsmear = tsmear->Eval(ereco);
		if( evalTsmear < 0. ) evalTsmear = 0.;
		double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
		double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
		caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
		caf.gastpc_RecoLepMomZ = preco / (sqrt(tan(reco_tx*0.001)*tan(reco_tx*0.001)+tan(reco_ty*0.001)*tan(reco_ty*0.001)+1));
		caf.gastpc_RecoLepMomX = caf.gastpc_RecoLepMomZ * tan(reco_tx*0.001); 
		caf.gastpc_RecoLepMomY = caf.gastpc_RecoLepMomZ * tan(reco_ty*0.001);
		if (fsTrkLenPerp[i] > 100.) {
		  // assume perfect charge reconstruction
		  caf.reco_q = ((fsPdg[i]==13 || fsPdg[i]==-211) > 0 ? -1 : 1);
		  caf.reco_numu=1; caf.reco_nue=0; caf.reco_nc=0;
		  caf.muon_tracker = 1; 
		}
	      }
	      else {
		piRecoE.push_back(ereco+0.1395);
		fsPdg[i] > 0 ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
	      }
	    } // Longest track is a pion
	    else if (abs(fsPdg[i])==211 && preco > 1.) {
	      // If very high energy and pi plus, there is a chance to reco as a proton
	      if (preco > 1.5 && fsPdg[i]==211) {
		// Energy reconstruction for particles in the ECAL
		double keTrue = sqrt(ptrue*ptrue + mass*mass) - mass;
		double keCalReco = keTrue * (1.+rando->Gaus(0., 0.2)+rando->Gaus(0., 0.3)/sqrt(keTrue));
		// Check if energy from momentum reco is closer to proton or pion expectation
		double eIfPi  = sqrt(preco*preco + 0.1395*0.1395) - 0.1395;
		double eIfPro = sqrt(preco*preco + 0.9383*0.9383) - 0.9383;
		// Is reco'd as a pion
		if (abs(eIfPi-keCalReco) < abs(eIfPro-keCalReco)) {
		  piRecoE.push_back(ereco+0.1395);
		  caf.gastpc_pi_pl_mult++;
		}
		else {
		  caf.Ev_reco -= 0.1395; // Remove pion mass
		  caf.gastpc_pro_mult++;
		} // pion mis-ID'd as a proton
	      } // Very high energy pion
	      else {
		if (preco > highestMipMom && rando->Rndm() < 1/(e*e*e)) {
		  if (highestMipMom > 0.) {
		    piRecoE.push_back(sqrt(highestMipMom*highestMipMom+0.1395*0.1395));
		    (highestMipPdg==-13 || highestMipPdg==211) ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
		  }
		  highestMipMom = preco;
		  highestMipPdg = fsPdg[i];
		  caf.Elep_reco = sqrt(preco*preco + mmu*mmu);
		  // angle reconstruction
		  double true_tx = 1000.*atan(fsPx[i] / fsPz[i]);
		  double true_ty = 1000.*atan(fsPy[i] / fsPz[i]);
		  double evalTsmear = tsmear->Eval(ereco);
		  if( evalTsmear < 0. ) evalTsmear = 0.;
		  double reco_tx = true_tx + rando->Gaus(0., evalTsmear/sqrt(2.));
		  double reco_ty = true_ty + rando->Gaus(0., evalTsmear/sqrt(2.));
		  caf.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
		  caf.gastpc_RecoLepMomZ = preco / (sqrt(tan(reco_tx*0.001)*tan(reco_tx*0.001)+tan(reco_ty*0.001)*tan(reco_ty*0.001)+1));
		  caf.gastpc_RecoLepMomX = caf.gastpc_RecoLepMomZ * tan(reco_tx*0.001); 
		  caf.gastpc_RecoLepMomY = caf.gastpc_RecoLepMomZ * tan(reco_ty*0.001);
		  if (fsTrkLen[i] > 100.) {
		    // assume perfect charge reconstruction
		    caf.reco_q = ((fsPdg[i]==13 || fsPdg[i]==-211) > 0 ? -1 : 1);
		    caf.reco_numu=1; caf.reco_nue=0; caf.reco_nc=0;
		    caf.muon_tracker = 1; 
		  }
		}
		else {
		  piRecoE.push_back(ereco+0.1395);
		  fsPdg[i] > 0 ? caf.gastpc_pi_pl_mult++ : caf.gastpc_pi_min_mult++;
		}
	      }
	    }
	  } // Is a MIP
	}
	else if( fsPdg[i] == 111 ){
	  // Decay pi0
	  TVector3 g1, g2;
	  TLorentzVector pi0( fsPx[i], fsPy[i], fsPz[i], fsE[i] );
	  decayPi0( pi0, g1, g2 );
	  double g1conv = rando->Exp( 14. ); // conversion distance
	  bool compton = (rando->Rndm() < 0.15); // dE/dX misID probability for photon
	  // Need to check if two photons are reco'd as separate (to check if pi0 is reco'd)
	  double ereco = 0.001 * (rando->Gaus(g1.Mag(), 0.1*g1.Mag()) + rando->Gaus(g2.Mag(), 0.1*g2.Mag()));
	  caf.Ev_reco += ereco;
	  // Photons are not overlaid on one another and low energy photon is not too low
	  if( !(g1conv < 2.0 && compton && (g2.Mag() < 50. || g1.Angle(g2) < 0.01)) ) {
	    caf.gastpc_pi_0_mult++;
	    piRecoE.push_back(ereco);
	    caf.gastpc_nRecoFS++;
	  }
	}
	else if (fsPdg[i] == 22 ) {
	  double ereco = 0.001 * rando->Gaus( fsE[i], 0.1*fsE[i] );
	  caf.Ev_reco += ereco;
	}
      }
         
      // Now loop over reconstructed pion energies and find the leading and sub-leading ones
      double leadpi    = 0.;
      double subleadpi = 0.;
      for (unsigned int pi=0; pi<piRecoE.size(); pi++) {
	if (piRecoE[pi] > leadpi) {
	  subleadpi = leadpi;
	  leadpi = piRecoE[pi];
	}
	else if (piRecoE[pi] > subleadpi && piRecoE[pi] < leadpi) {
	  subleadpi = piRecoE[pi];
	}
      }
      caf.gastpc_lead_pi_E = leadpi;
      caf.gastpc_sublead_pi_E = subleadpi;
    }

    //printf( "Ev reco %f pion mult %d %d Elep reco %f reco numu %d reco q %d Ehad_veto %f muon_tracker %d\n", caf.Ev_reco, caf.gastpc_pi_pl_mult, caf.gastpc_pi_min_mult, caf.Elep_reco, caf.reco_numu, caf.reco_q, caf.Ehad_veto, caf.muon_tracker );

    caf.fill();
  }

  // set POT
  caf.meta_run = par.run;
  caf.meta_subrun = par.subrun;

}

int main( int argc, char const *argv[] ) 
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  // get command line options
  std::string ghepdir;
  std::string outfile;
  std::string edepfile;
  std::string fhicl_filename;

  // Make parameter object and set defaults
  params par;
  par.IsGasTPC = false;
  par.OA_xcoord = 0.; // on-axis by default
  par.fhc = true;
  par.grid = false;
  par.seed = 7; // a very random number
  par.run = 1; // CAFAna doesn't like run number 0
  par.subrun = 0;
  par.n = -1;
  par.nfiles = 1;
  par.first = 0;
  par.trk_muRes = 0.02; // fractional muon energy resolution of HP GAr TPC
  par.LAr_muRes = 0.05; // fractional muon energy resolution of muons contained in LAr
  par.ECAL_muRes = 0.1; // fractional muon energy resolution of muons ending in ECAL
  par.em_const = 0.03; // EM energy resolution constant term: A + B/sqrt(E) (GeV)
  par.em_sqrtE = 0.1; // EM energy resolution 1/sqrt(E) term: A + B/sqrt(E) (GeV)
  par.michelEff = 0.75; // Michel finder efficiency
  par.CC_trk_length = 100.; // minimum track length for CC in cm
  par.pileup_frac = 0.1; // fraction of events with non-zero pile-up
  par.pileup_max = 0.5; // GeV
  par.gastpc_len = 6.; // track length cut in cm
  par.gastpc_B = 0.4; // B field strength in Tesla
  par.gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
  par.gastpc_X0 = 1300.; // cm = 13m radiation length

  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--edepfile") ) {
      edepfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--ghepdir") ) {
      ghepdir = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--outfile") ) {
      outfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--fhicl") ) {
      fhicl_filename = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--seed") ) {
      par.seed = atoi(argv[i+1]);
      par.run = par.seed;
      i += 2;
    } else if( argv[i] == std::string("--nevents") ) {
      par.n = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--nfiles") ) {
      par.nfiles = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--first") ) {
      par.first = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--oa") ) {
      par.OA_xcoord = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--grid") ) {
      par.grid = true;
      i += 1;
    } else if( argv[i] == std::string("--rhc") ) {
      par.fhc = false;
      i += 1;
    } else if( argv[i] == std::string("--gastpc") ) {
      par.IsGasTPC = true;
      i += 1;
    } else i += 1; // look for next thing
  }

  printf( "Making CAF from edep-sim tree dump: %s\n", edepfile.c_str() );
  printf( "Searching for GENIE ghep files here: %s\n", ghepdir.c_str() );
  if( par.fhc ) printf( "Running neutrino mode (FHC)\n" );
  else printf( "Running antineutrino mode (RHC)\n" );
  if( par.grid ) printf( "Running grid mode\n" );
  else printf( "Running test mode\n" );
  printf( "Output CAF file: %s\n", outfile.c_str() );
  if( par.IsGasTPC ) printf( "Running gas TPC\n" );

  rando = new TRandom3( par.seed );
  CAF caf( outfile, par.IsGasTPC );

  // LAr driven smearing, maybe we want to change for gas?
  tsmear = new TF1( "tsmear", "0.162 + 3.407*pow(x,-1.) + 3.129*pow(x,-0.5)", 0., 999.9 );
  // For checking pion ECAL interaction
  TFile *fMip = new TFile("p_range.root","read");
  TGraph *grMip = (TGraph*)fMip->Get("non_int_prob");
  sp = new TSpline3("sp", grMip);
  fMip->Close();

  TFile * tf = new TFile( edepfile.c_str() );
  TTree * tree = (TTree*) tf->Get( "tree" );

  loop( caf, par, tree, ghepdir, fhicl_filename );

  caf.version = 4;
  printf( "Run %d POT %g\n", caf.meta_run, caf.pot );
  caf.fillPOT();
  caf.write();

  caf.cafFile->Close();

  printf( "-30-\n" );


}
