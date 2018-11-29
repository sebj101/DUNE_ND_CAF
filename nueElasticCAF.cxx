#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <stdio.h>
#include <math.h>
#include "CAF.C"

const double me  = 0.511E-3; // GeV
const double mmu = 0.10537; // GeV

double getPOT( TTree * meta )
{

  double totpot = 0.;

  double pot;
  meta->SetBranchAddress( "pot", &pot );
  for( int i = 0; i < meta->GetEntries(); ++i ) {
    meta->GetEntry(i);
    totpot += pot;
  }

  return totpot;

}

void decayPi0( TLorentzVector pi0, TVector3 &gamma1, TVector3 &gamma2, TRandom3 * rand )
{
  double e = pi0.E();
  double mp = 0.1349766; // pi0 mass

  double beta = sqrt( 1. - (mp*mp)/(e*e) ); // velocity of pi0
  double theta = 3.1416*rand->Rndm(); // theta of gamma1 w.r.t. pi0 direction
  double phi = 2.*3.1416*rand->Rndm(); // phi of gamma1 w.r.t. pi0 direction

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

void loop( TTree * tree, int cat, CAF &caf )
{

  TRandom3 * rando = new TRandom3(12345);

  // initialize variable resolutions
  TF1 * esmear = new TF1( "esmear", "0.03 + 0.05*pow(x,-0.5)", 0., 999.9 );
  TF1 * tsmear = new TF1( "tsmear", "3. + 0.162 + 3.407*pow(x,-1.) + 3.129*pow(x,-0.5)", 0., 999.9 );

  // Signal
  if( cat == 0 ) {
    // event variables
    int evtNo, nfsp;
    double Enu, vtx_x, vtx_y, vtx_z, nu_thetaX, nu_thetaY;

    // particle variables
    int pdg[100];
    double E[100];
    //double px[100], py[100], pz[100];
    double perf_px[100], perf_py[100], perf_pz[100];
    double best_px[100], best_py[100], best_pz[100];

    tree->SetBranchAddress( "evt", &evtNo );
    tree->SetBranchAddress( "Enu", &Enu );
    tree->SetBranchAddress( "vtx_x", &vtx_x );
    tree->SetBranchAddress( "vtx_y", &vtx_y );
    tree->SetBranchAddress( "vtx_z", &vtx_z );
    tree->SetBranchAddress( "nu_thetaX", &nu_thetaX );
    tree->SetBranchAddress( "nu_thetaY", &nu_thetaY );
    tree->SetBranchAddress( "nfsp", &nfsp );
    tree->SetBranchAddress( "pdg", pdg );
    //tree->SetBranchAddress( "px", px );
    //tree->SetBranchAddress( "py", py );
    //tree->SetBranchAddress( "pz", pz );
    tree->SetBranchAddress( "perf_px", perf_px );
    tree->SetBranchAddress( "perf_py", perf_py );
    tree->SetBranchAddress( "perf_pz", perf_pz );
    tree->SetBranchAddress( "best_px", best_px );
    tree->SetBranchAddress( "best_py", best_py );
    tree->SetBranchAddress( "best_pz", best_pz );
    tree->SetBranchAddress( "E", E );
  
    int N = tree->GetEntries();

    // loop over events
    for( int ii = 0; ii < N; ++ii ) {
      tree->GetEntry( ii );

      if( ii % 1000000 == 0 ) printf( "Event %d of %d...\n", ii, N );

      // Set basic CAF variables
      caf.run = 10;
      caf.subrun = 10;
      caf.event = ii;
      caf.isCC = 0;
      caf.mode = 7; // nu+e
      caf.LepPDG = 11; 
      caf.Ev = Enu;
      caf.Q2 = 0.;
      caf.W = 0.;
      caf.X = 0.;
      caf.Y = 0.;
      caf.NuMomX = Enu*sin(nu_thetaX);
      caf.NuMomY = Enu*sin(nu_thetaY);
      caf.NuMomZ = sqrt( Enu*Enu - caf.NuMomX*caf.NuMomX - caf.NuMomY-caf.NuMomY );

      caf.nP = 0; caf.nN = 0; caf.nipip = 0; caf.nipi0 = 0; caf.nikm = 0; caf.nik0 = 0; 
      caf.niem = 0; caf.niother = 0; caf.nNucleus = 0; caf.nUNKNOWN = 0;
      caf.eP = 0.; caf.eN = 0.; caf.ePip = 0.; caf.ePi0 = 0.; caf.eOther = 0.;

      caf.vtx_x = vtx_x;
      caf.vtx_y = vtx_y;
      caf.vtx_z = vtx_z;

      // nonsense variables
      caf.reco_numu = 0;
      caf.reco_nue = 1;
      caf.reco_nc = 0;
      caf.reco_q = 0;
      caf.muon_contained = -1; caf.muon_tracker = -1; caf.muon_ecal = -1; caf.muon_exit = -1;
      caf.reco_lepton_pdg = 11;
      caf.pileup_energy = 0.;

      // loop over final-state particles
      for( int i = 0; i < nfsp; ++i ) {

        if( pdg[i] == 11 ) {
          double Ttrue = E[i] - me;

          double best_thetaX = 1000.*atan( best_px[i] / best_pz[i] );
          double best_thetaY = 1000.*atan( best_py[i] / best_pz[i] );

          double evalEsmear = esmear->Eval(Ttrue);
          if( evalEsmear < 0. ) evalEsmear = 0.;

          double evalTsmear = tsmear->Eval(Ttrue);
          if( evalTsmear < 0. ) evalTsmear = 0.;

          double ereco = Ttrue * ( 1. + rando->Gaus(0., evalEsmear) );
          double smearx = best_thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
          double smeary = best_thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
          double reco_theta = sqrt( smearx*smearx + smeary*smeary );

          // Lepton truth info
          caf.LepMomX = perf_px[i]; // wrt true neutrino
          caf.LepMomY = perf_py[i];
          caf.LepMomZ = perf_pz[i];
          caf.LepE = E[i];
          double true_tx = atan( perf_px[i] / perf_pz[i] );
          double true_ty = atan( perf_py[i] / perf_pz[i] );
          caf.LepNuAngle = sqrt( true_tx*true_tx + true_ty*true_ty );

          double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
          double reco_enu = ereco / reco_y;

          // fill CAF
          caf.Elep_reco = ereco;
          caf.theta_reco = reco_theta;
          caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative
          caf.Ehad_veto = 0.;
          
        } else if( abs(pdg[i] == 12) || abs(pdg[i]) == 14 ) {
          caf.neutrinoPDG = pdg[i]; // actually don't know
          caf.neutrinoPDGunosc = pdg[i];
        }
      } // fsp
    } // evt loop
  } // if is signal
  else { // bkg
    // event variables
    int evtNo, nfsp, scat;
    double Enu, Q2, W, x, y, nu_thetaX, nu_thetaY;
    double vtx_x, vtx_y, vtx_z;

    // particle variables
    int pdg[100];
    double E[100];
    double perf_px[100], perf_py[100], perf_pz[100];
    double best_px[100], best_py[100], best_pz[100];

    tree->SetBranchAddress( "evt", &evtNo );
    tree->SetBranchAddress( "Enu", &Enu );
    tree->SetBranchAddress( "Q2", &Q2 );
    tree->SetBranchAddress( "W", &W );
    tree->SetBranchAddress( "x", &x );
    tree->SetBranchAddress( "y", &y );
    tree->SetBranchAddress( "nfsp", &nfsp );
    tree->SetBranchAddress( "scat", &scat );
    tree->SetBranchAddress( "nu_thetaX", &nu_thetaX );
    tree->SetBranchAddress( "nu_thetaY", &nu_thetaY );
    tree->SetBranchAddress( "pdg", pdg );
    tree->SetBranchAddress( "best_px", best_px );
    tree->SetBranchAddress( "best_py", best_py );
    tree->SetBranchAddress( "best_pz", best_pz );
    tree->SetBranchAddress( "perf_px", perf_px );
    tree->SetBranchAddress( "perf_py", perf_py );
    tree->SetBranchAddress( "perf_pz", perf_pz );
    tree->SetBranchAddress( "E", E );
    tree->SetBranchAddress( "vtx_x", &vtx_x );
    tree->SetBranchAddress( "vtx_y", &vtx_y );
    tree->SetBranchAddress( "vtx_z", &vtx_z );

    int N = tree->GetEntries();
    for( int ii = 0; ii < 10000; ++ii ) {
      tree->GetEntry( ii );

      if( ii % 100000 == 0 ) printf( "Event %d of %d...\n", ii, N );

      // "background" sample is the full physics list, which includes nu+e signal, which we don't want to count as background
      if( scat == 7 ) continue;

      // Set basic CAF variables
      caf.run = 10 + cat;
      caf.subrun = 10 + cat;
      caf.event = ii;
      caf.isCC = 0; // change to 1 if we find charged lepton
      caf.mode = scat;
      caf.Ev = Enu;
      caf.Q2 = Q2;
      caf.W = W;
      caf.X = x;
      caf.Y = y;
      caf.NuMomX = Enu*sin(nu_thetaX);
      caf.NuMomY = Enu*sin(nu_thetaY);
      caf.NuMomZ = sqrt( Enu*Enu - caf.NuMomX*caf.NuMomX - caf.NuMomY-caf.NuMomY );

      caf.nP = 0; caf.nN = 0; caf.nipip = 0; caf.nipi0 = 0; caf.nikm = 0; caf.nik0 = 0; 
      caf.niem = 0; caf.niother = 0; caf.nNucleus = 0; caf.nUNKNOWN = 0;
      caf.eP = 0.; caf.eN = 0.; caf.ePip = 0.; caf.ePi0 = 0.; caf.eOther = 0.;

      caf.vtx_x = vtx_x;
      caf.vtx_y = vtx_y;
      caf.vtx_z = vtx_z;

      // nonsense variables
      caf.reco_numu = 0;
      caf.reco_nue = 1;
      caf.reco_nc = 0;
      caf.reco_q = 0;
      caf.muon_contained = -1; caf.muon_tracker = -1; caf.muon_ecal = -1; caf.muon_exit = -1;
      caf.reco_lepton_pdg = 11;

      double extraE = 0.;
      int electron_candidates = 0;
      int photon_candidates = 0;

      // loop over final-state particles
      for( int i = 0; i < nfsp; ++i ) {

        double ke = 0.001*(E[i]*E[i] - sqrt(E[i]*E[i] - best_px[i]*best_px[i] - best_py[i]*best_py[i] - best_pz[i]*best_pz[i]));
        if( abs(pdg[i]) == 14 || abs(pdg[i]) == 12 ) {
          caf.LepPDG = pdg[i];
          caf.isCC = 0;
          caf.neutrinoPDG = pdg[i];
          caf.neutrinoPDGunosc = pdg[i];

          caf.LepMomX = perf_px[i]; // wrt true neutrino
          caf.LepMomY = perf_py[i];
          caf.LepMomZ = perf_pz[i];
          caf.LepE = E[i];
          double true_tx = atan( perf_px[i] / perf_pz[i] );
          double true_ty = atan( perf_py[i] / perf_pz[i] );
          caf.LepNuAngle = sqrt( true_tx*true_tx + true_ty*true_ty );
        } else if( abs(pdg[i]) == 13 || abs(pdg[i]) == 11 ) {
          caf.LepPDG = pdg[i];
          caf.isCC = 1;
          caf.neutrinoPDG = (pdg[i] > 0 ? pdg[i]+1 : pdg[i]-1);
          caf.neutrinoPDGunosc = caf.neutrinoPDG;

          caf.LepMomX = perf_px[i]; // wrt true neutrino
          caf.LepMomY = perf_py[i];
          caf.LepMomZ = perf_pz[i];
          caf.LepE = E[i];
          double true_tx = atan( perf_px[i] / perf_pz[i] );
          double true_ty = atan( perf_py[i] / perf_pz[i] );
          caf.LepNuAngle = sqrt( true_tx*true_tx + true_ty*true_ty );
        }
        else if( pdg[i] == 2212 ) {caf.nP++; caf.eP += ke;}
        else if( pdg[i] == 2112 ) {caf.nN++; caf.eN += ke;}
        else if( pdg[i] ==  211 ) {caf.nipip++; caf.ePip += ke;}
        else if( pdg[i] == -211 ) {caf.nipim++; caf.ePim += ke;}
        else if( pdg[i] ==  111 ) {caf.nipi0++; caf.ePi0 += ke;}
        else if( pdg[i] ==  321 ) {caf.nikp++; caf.eOther += ke;}
        else if( pdg[i] == -321 ) {caf.nikm++; caf.eOther += ke;}
        else if( pdg[i] == 311 || pdg[i] == -311 || pdg[i] == 130 || pdg[i] == 310 ) {caf.nik0++; caf.eOther += ke;}
        else if( pdg[i] ==   22 ) {caf.niem++; caf.eOther += ke;}
        else if( pdg[i] > 1000000000 ) caf.nNucleus++;
        else {caf.niother++; caf.eOther += ke;}

        // final-state electron
        if( abs(pdg[i]) == 11 ) {
          ++electron_candidates;

          double Ttrue = E[i] - me;

          double thetaX = atan( best_px[i] / best_pz[i] ); // true angle, including beam divergence
          double thetaY = atan( best_py[i] / best_pz[i] );

          double evalEsmear = esmear->Eval(Ttrue);
          if( evalEsmear < 0. ) evalEsmear = 0.;

          double evalTsmear = tsmear->Eval(Ttrue);
          if( evalTsmear < 0. ) evalTsmear = 0.;

          double ereco = Ttrue * ( 1. + rando->Gaus(0., evalEsmear) );
          double smearx = 1000*thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
          double smeary = 1000*thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
          double reco_theta = sqrt( smearx*smearx + smeary*smeary );

          // Lepton truth info
          caf.LepMomX = perf_px[i]; // wrt true neutrino
          caf.LepMomY = perf_py[i];
          caf.LepMomZ = perf_pz[i];
          caf.LepE = E[i];
          double true_tx = atan( perf_px[i] / perf_pz[i] );
          double true_ty = atan( perf_py[i] / perf_pz[i] );
          caf.LepNuAngle = sqrt( true_tx*true_tx + true_ty*true_ty );

          double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
          double reco_enu = ereco / reco_y;

          // fill CAF
          caf.Elep_reco = ereco;
          caf.theta_reco = reco_theta;
          caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative
          caf.Ehad_veto = 0.;

        }
        else if( pdg[i] == 111 ) { // pi0 production

          TVector3 gamma1, gamma2;
          TLorentzVector pi0( best_px[i], best_py[i], best_pz[i], E[i] );
          decayPi0( pi0, gamma1, gamma2, rando ); // sets photon vectors

          double evalEsmear = esmear->Eval(gamma1.Mag());
          if( evalEsmear < 0. ) evalEsmear = 0.;
          double evalTsmear = tsmear->Eval(gamma1.Mag());
          if( evalTsmear < 0. ) evalTsmear = 0.;

          double reco_e_g1 = gamma1.Mag() * ( 1. + rando->Gaus(0., evalEsmear) );
          double reco_e_g2 = gamma2.Mag() * ( 1. + rando->Gaus(0., evalEsmear) );

          double ereco = 0.;
          double reco_theta = 0.;

          // plausible to reconstruct if a) gamma2 is < 50 MeV, b) angle is < resolution
          if( reco_e_g2 < 0.05 ) {
            ++photon_candidates;
            ereco = reco_e_g1;
            double thetaX = atan( gamma1.x() / gamma1.z() );
            double thetaY = atan( gamma1.y() / gamma1.z() ); // convert to mrad for smearing
            double smearx = 1000*thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
            double smeary = 1000*thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
            reco_theta = sqrt( smearx*smearx + smeary*smeary );
          } else if( 1000.*gamma1.Angle(gamma2) < evalTsmear ) {
            ++photon_candidates;
            ereco = reco_e_g1 + reco_e_g2;
            double thetaX = atan( gamma1.x() / gamma1.z() );
            double thetaY = atan( gamma1.y() / gamma1.z() ); // convert to mrad for smearing
            double smearx = 1000*thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
            double smeary = 1000*thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
            reco_theta = sqrt( smearx*smearx + smeary*smeary );
          } else {
            extraE += (reco_e_g1 + reco_e_g2);
          }

          double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
          double reco_enu = ereco / reco_y;

          // fill CAF
          caf.Elep_reco = ereco;
          caf.theta_reco = reco_theta;
          caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative

        } else if( abs(pdg[i]) == 12 || abs(pdg[i]) == 14 || pdg[i] == 2112 || pdg[i] > 9999 ) { // neutrinos, neutrons, nuclear fragments
          continue; // skip these; they contribute nothing to extra energy
        } else if( abs(pdg[i]) == 211 ) { // charged pion
          extraE += E[i] - 0.1395;
        } else if( pdg[i] == 2212 ) { // proton -- should there be a threshold?
          extraE += E[i] - 0.9382;
        } else { // anything exotic (kaons?) assume it will cause event to be rejected; it's mass alone will put it over this cut
          extraE += E[i];
        }
      } // fsp loop

      caf.Ehad_veto = extraE;
      caf.pileup_energy = 0.;

      if( cat == 1 && electron_candidates == 1 && photon_candidates == 0 ) caf.fill();
      else if( cat == 2 && electron_candidates == 0 && photon_candidates == 1 ) caf.fill();
    }
  }
}

int main()
{

  // Signal
  TChain * tree = new TChain( "tree", "tree" );
  tree->Add( "/dune/data/users/marshalc/nue/FHC/*.root" );
  
  TChain * meta = new TChain( "meta", "meta" );
  meta->Add( "/dune/data/users/marshalc/nue/FHC/*.root" );
  double pot = getPOT( meta );

  CAF signal( "/dune/data/users/marshalc/CAFs/mcc11_v2/ND_nue_signal.root" );
  //loop( tree, pot, signal );

  printf( "Got %g POT for %lld events\n", pot, tree->GetEntries() );

}

