#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <stdio.h>
#include <math.h>
#include "nusystematics/artless/response_helper.hh"
#include "CAF.C"

// genie includes
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Conventions/Units.h"
#include "FluxDrivers/GSimpleNtpFlux.h"

using namespace genie;

const double me  = 0.511E-3; // GeV
const double mmu = 0.10537; // GeV

// initialize variable resolutions
TF1 * esmear;
TF1 * tsmear;

TRandom3 * rando;

nusyst::response_helper rh( "./fhicl.fcl" );

void init()
{
  esmear = new TF1( "esmear", "0.03 + 0.05*pow(x,-0.5)", 0., 999.9 );
  tsmear = new TF1( "tsmear", "3. + 0.162 + 3.407*pow(x,-1.) + 3.129*pow(x,-0.5)", 0., 999.9 );

  rando = new TRandom3(12345);
}

void decayPi0( TLorentzVector &pi0, TVector3 &gamma1, TVector3 &gamma2 )
{
  double e = pi0.E();
  double mp = 0.1349766; // pi0 mass

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

void RotateZu( TVector3 &v, TVector3 uu )
{

  // new z axis
  TVector3 u = uu.Unit();
  double u1 = u.x();
  double u2 = u.y();
  double u3 = u.z();
  double up = u1*u1 + u2*u2;

  // new components
  double fX = v.x();
  double fY = v.y();
  double fZ = v.z();
  if( up ) {
    up = TMath::Sqrt(up);
    // old vector components
    double px = v.x();
    double py = v.y();
    double pz = v.z();

    // for RotateUz
    //fX = (u1*u3*px - u2*py + u1*up*pz)/up;
    //fY = (u2*u3*px + u1*py + u2*up*pz)/up;
    //fZ = (u3*u3*px -    px + u3*up*pz)/up;

    fX = (-u2*px   + u1*py)/up;
    fY = -(u1*u3*px + u2*u3*py - up*up*pz)/up;
    fZ = (u1*(1.-u3*u3)*px + u2*(1.-u3*u3)*py + u3*up*up*pz)/(up*up);

  } else if( u3 < 0. ) { // theta == pi, phi == 0
    fX = -fX;
    fZ = -fZ;
  } else {
    printf( "Can't rotate the vector\n" );
  }

  v.SetXYZ( fX, fY, fZ );
}

void loop( TTree * tree, int cat, CAF &caf )
{

  //nusyst::response_helper rh( "./fhicl.fcl" );

  // midpoint of the decay pipe, relative to detector center at (0,0,0)
  TVector3 origin(0., 4823.6, -46048.);

  // Signal
  int N = tree->GetEntries();

  NtpMCEventRecord * mcrec = NULL;
  tree->SetBranchAddress("gmcrec", &mcrec);

  if( cat == 0 ) {

    // loop over events
    for( int ii = 0; ii < N; ++ii ) {
      tree->GetEntry( ii );

      caf.setToBS();

      // Set basic CAF variables
      caf.run = 10;
      caf.subrun = 10;
      caf.event = ii;
      caf.isCC = 0;
      caf.mode = 7; // nu+e
      caf.LepPDG = 11; 

      caf.nP = 0; caf.nN = 0; caf.nipip = 0; caf.nipi0 = 0; caf.nikm = 0; caf.nik0 = 0; 
      caf.niem = 0; caf.niother = 0; caf.nNucleus = 0; caf.nUNKNOWN = 0;
      caf.eP = 0.; caf.eN = 0.; caf.ePip = 0.; caf.ePi0 = 0.; caf.eOther = 0.;

      // nonsense variables
      caf.reco_numu = 0;
      caf.reco_nue = 1;
      caf.reco_nc = 0;
      caf.reco_q = 0;
      caf.muon_contained = -1; caf.muon_tracker = -1; caf.muon_ecal = -1; caf.muon_exit = -1;
      caf.reco_lepton_pdg = 11;
      caf.pileup_energy = 0.;

      // get the GENIE event
      EventRecord *event = mcrec->event;
      Interaction *in = event->Summary();

      TLorentzVector lep = in->Kine().FSLeptonP4();
      TLorentzVector nu = *(in->InitState().GetProbeP4(kRfLab));
      TVector3 nudir = nu.Vect().Unit();
      TLorentzVector q = nu-lep;

      TVector3 vtxO = event->Vertex()->Vect();
      TVector3 vtx( vtxO.x()*100., vtxO.y()*100. - 305., vtxO.z()*100. - 5. );
      caf.vtx_x = vtx.x();
      caf.vtx_y = vtx.y();
      caf.vtx_z = vtx.z();

      caf.NuMomX = nu.X();
      caf.NuMomY = nu.Y();
      caf.NuMomZ = nu.Z();

      // interaction-level variables
      caf.X = -q.Mag2()/(2*0.939*q.E());
      caf.Y = q.E() / nu.E();
      caf.Q2 = -q.Mag2();
      caf.W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2());
      caf.mode = in->ProcInfo().ScatteringTypeId();
      caf.Ev = in->InitState().ProbeE(kRfLab);

      // Loop over all particles in this event, fill particle variables
      GHepParticle * p = 0;
      TIter event_iter(event);

      while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

        if( p->Status() == kIStStableFinalState ) {
          if( p->Pdg() == 11 ) {
            TLorentzVector mom = *(p->P4());
            double Ttrue = mom.E() - me;

            TVector3 bestnudir = (vtx - origin).Unit();
            TVector3 best = mom.Vect();
            RotateZu( best, bestnudir );

            TVector3 perf = mom.Vect();
            RotateZu( perf, nudir );

            double best_thetaX = 1000.*atan( best.x() / best.z() );
            double best_thetaY = 1000.*atan( best.y() / best.z() );

            double evalEsmear = esmear->Eval(Ttrue);
            if( evalEsmear < 0. ) evalEsmear = 0.;

            double evalTsmear = tsmear->Eval(Ttrue);
            if( evalTsmear < 0. ) evalTsmear = 0.;

            double ereco = Ttrue * ( 1. + rando->Gaus(0., evalEsmear) );
            double smearx = best_thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
            double smeary = best_thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
            double reco_theta = sqrt( smearx*smearx + smeary*smeary );

            // Lepton truth info
            caf.LepMomX = perf.x(); // wrt true neutrino
            caf.LepMomY = perf.y();
            caf.LepMomZ = perf.z();
            caf.LepE = Ttrue;
            caf.LepNuAngle = perf.Angle( nudir );

            double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
            double reco_enu = ereco / reco_y;

            // fill CAF
            caf.Elep_reco = ereco;
            caf.theta_reco = 0.001*reco_theta;
            caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative
            caf.Ehad_veto = 0.;
          } // if electron
          else if( abs(p->Pdg()) == 12 || abs(p->Pdg()) == 14 ) {
            caf.neutrinoPDG = p->Pdg();
            caf.neutrinoPDGunosc = p->Pdg();
          }
        }
      }// end loop over particles 

      caf.fill();

      // clear current mc event record
      mcrec->Clear();
    } // event loop
  } // if signal
  else { // bkg

    // loop over events
    for( int ii = 0; ii < N; ++ii ) {
      tree->GetEntry( ii );

      caf.setToBS();

      // Set basic CAF variables
      caf.run = 10 + cat;
      caf.subrun = 10 + cat;
      caf.event = ii;
      caf.isCC = 0; // overwritten in particle loop
      caf.LepPDG = 0; // overwritten in particle loop

      caf.nP = 0; caf.nN = 0; caf.nipip = 0; caf.nipi0 = 0; caf.nikm = 0; caf.nik0 = 0; 
      caf.niem = 0; caf.niother = 0; caf.nNucleus = 0; caf.nUNKNOWN = 0;
      caf.eP = 0.; caf.eN = 0.; caf.ePip = 0.; caf.ePi0 = 0.; caf.eOther = 0.;

      // nonsense variables
      caf.reco_numu = 0;
      caf.reco_nue = 1;
      caf.reco_nc = 0;
      caf.reco_q = 0;
      caf.muon_contained = -1; caf.muon_tracker = -1; caf.muon_ecal = -1; caf.muon_exit = -1;
      caf.reco_lepton_pdg = 11;
      caf.pileup_energy = 0.;

      // get the GENIE event
      EventRecord *event = mcrec->event;
      Interaction *in = event->Summary();

/*
      systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
      for( systtools::event_unit_response_w_cv_t::iterator it = resp.begin(); it != resp.end(); ++it ) {
        caf.nwgt[(*it).pid] = (*it).responses.size();
        caf.cvwgt[(*it).pid] = (*it).CV_response;
        for( unsigned int i = 0; i < (*it).responses.size(); ++i ) {
          caf.wgt[(*it).pid][i] = (*it).responses[i];
        }
      }
*/
      TVector3 vtxO = event->Vertex()->Vect();
      TVector3 vtx( vtxO.x()*100., vtxO.y()*100. - 305., vtxO.z()*100. - 5. );
      caf.vtx_x = vtx.x();
      caf.vtx_y = vtx.y();
      caf.vtx_z = vtx.z();

      TLorentzVector lep = in->Kine().FSLeptonP4();
      TLorentzVector nu = *(in->InitState().GetProbeP4(kRfLab));
      TVector3 nudir = nu.Vect().Unit();
      TLorentzVector q = nu-lep;

      caf.NuMomX = nu.X();
      caf.NuMomY = nu.Y();
      caf.NuMomZ = nu.Z();

      // interaction-level variables
      caf.X = -q.Mag2()/(2*0.939*q.E());
      caf.Y = q.E() / nu.E();
      caf.Q2 = -q.Mag2();
      caf.W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2());
      caf.mode = in->ProcInfo().ScatteringTypeId();
      caf.Ev = in->InitState().ProbeE(kRfLab);

      if( caf.mode == 7 ) {
        mcrec->Clear();
        continue; // nu+e in background file
      }

      // Loop over all particles in this event, fill particle variables
      GHepParticle * p = 0;
      TIter event_iter(event);

      double extraE = 0.;
      int electron_candidates = 0;
      int photon_candidates = 0;

      while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

        if( p->Status() == kIStStableFinalState ) {

          TLorentzVector mom = *(p->P4());
          double ke = mom.E() - mom.M();
          int pdg = p->Pdg();

          if( abs(pdg) >= 11 && abs(pdg) <= 16 ) {
            caf.LepPDG = pdg;

            TVector3 perf = mom.Vect();
            RotateZu( perf, nudir );

            caf.LepMomX = perf.x(); // wrt true neutrino
            caf.LepMomY = perf.y();
            caf.LepMomZ = perf.z();
            caf.LepE = mom.E();
            caf.LepNuAngle = perf.Angle( nudir );

            if( abs(pdg) == 11 || abs(pdg) == 13 ) {
              caf.isCC = 1;
              caf.neutrinoPDG = (pdg > 0 ? pdg+1 : pdg-1);
              caf.neutrinoPDGunosc = caf.neutrinoPDG;
            } else {
              caf.isCC = 0;
              caf.neutrinoPDG = pdg;
              caf.neutrinoPDGunosc = caf.neutrinoPDG;
            }
          } // if lepton
          else if( pdg == 2212 ) {caf.nP++; caf.eP += ke;}
          else if( pdg == 2112 ) {caf.nN++; caf.eN += ke;}
          else if( pdg ==  211 ) {caf.nipip++; caf.ePip += ke;}
          else if( pdg == -211 ) {caf.nipim++; caf.ePim += ke;}
          else if( pdg ==  111 ) {caf.nipi0++; caf.ePi0 += ke;}
          else if( pdg ==  321 ) {caf.nikp++; caf.eOther += ke;}
          else if( pdg == -321 ) {caf.nikm++; caf.eOther += ke;}
          else if( pdg == 311 || pdg == -311 || pdg == 130 || pdg == 310 ) {caf.nik0++; caf.eOther += ke;}
          else if( pdg ==   22 ) {caf.niem++; caf.eOther += ke;}
          else if( pdg > 1000000000 ) caf.nNucleus++;
          else {caf.niother++; caf.eOther += ke;}


          // background reco stuff now
          if( abs(pdg) == 11 ) {
            ++electron_candidates;

            TVector3 bestnudir = (vtx - origin).Unit();
            TVector3 best = mom.Vect();
            RotateZu( best, bestnudir );

            TVector3 perf = mom.Vect();
            RotateZu( perf, nudir );

            double thetaX = 1000.*atan( best.x() / best.z() );
            double thetaY = 1000.*atan( best.y() / best.z() );
            double Ttrue = mom.E() - me;

            double evalEsmear = esmear->Eval(Ttrue);
            if( evalEsmear < 0. ) evalEsmear = 0.;

            double evalTsmear = tsmear->Eval(Ttrue);
            if( evalTsmear < 0. ) evalTsmear = 0.;

            double ereco = Ttrue * ( 1. + rando->Gaus(0., evalEsmear) );
            double smearx = thetaX + rando->Gaus(0., evalTsmear/sqrt(2.));
            double smeary = thetaY + rando->Gaus(0., evalTsmear/sqrt(2.));
            double reco_theta = sqrt( smearx*smearx + smeary*smeary );

            double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
            double reco_enu = ereco / reco_y;

            // fill CAF
            caf.Elep_reco = ereco;
            caf.theta_reco = 0.001*reco_theta;
            caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative
          }
          else if( pdg == 111 ) { // pi0 production

            TVector3 gamma1, gamma2;
            TLorentzVector pi0( mom.X(), mom.Y(), mom.Z(), mom.E() );
            decayPi0( pi0, gamma1, gamma2 ); // sets photon vectors

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

            if( cat == 2 ) {
              double reco_y = 1. - (ereco * (1. - cos(reco_theta/1000.)))/me;
              double reco_enu = ereco / reco_y;

              caf.Elep_reco = ereco;
              caf.theta_reco = 0.001*reco_theta;
              caf.Ev_reco = reco_enu; // 2D neutrino energy reco, can sometimes be negative
            }  

          } else if( abs(pdg) == 12 || abs(pdg) == 14 || pdg == 2112 || pdg > 9999 ) { // neutrinos, neutrons, nuclear fragments
            continue; // skip these; they contribute nothing to extra energy
          } else if( abs(pdg) == 211 || pdg == 2212 ) { // charged pion
            extraE += mom.E() - mom.M();
            if( pdg == 211 ) caf.pileup_energy = 1.; // Michel veto?
          } else {
            extraE += mom.E();
          }
        } // fsp if stable fs
      } // fsp loop

      if( electron_candidates + photon_candidates == 1 ) {

        caf.Ehad_veto = extraE*1000.; // MeV

        if( cat == 1 && electron_candidates == 1 && photon_candidates == 0 ) caf.fill();
        else if( cat == 2 && electron_candidates == 0 && photon_candidates == 1 ) caf.fill();
      }

      // clear current mc event record
      mcrec->Clear();

    } // event loop

  } // if bkg
}

int main()
{

  init();
  std::vector<unsigned int> parIds = rh.GetParameters();
/*
  // nu+e signal
  CAF signal( "/dune/data/users/marshalc/CAFs/mcc11_v2/ND_nue_signal.root" );
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations.size() );
    bool is_wgt = head.isWeightSystematicVariation;
    std::string wgt_var = ( is_wgt ? "wgt" : "var" );
    signal.addRWbranch( parIds[i], head.prettyName, wgt_var, head.paramVariations );
    signal.iswgt[parIds[i]] = is_wgt;
  }
  signal.pot = 0.;
  signal.meta_run = 0;
  signal.meta_subrun = 0;
  signal.version = 2;
  for( int i = 0; i <= 999; ++i ) {
    printf( "nu+e signal file %d, so far %4.4g POT\n", i, signal.pot );
    TFile * tf = new TFile( Form("/pnfs/dune/persistent/users/marshalc/CAF/genieNuESignal/FHC/LAr.neutrino.%d.ghep.root",i) );
    if( tf == NULL ) delete tf;
    else {
      TTree * tree = (TTree*) tf->Get("gtree");
      if( tree && !tf->TestBit(TFile::kRecovered) ) {
        signal.pot += 1.0E21;
        loop( tree, 0, signal );
      }
      delete tree;
      tf->Close();
      delete tf;
    }
  }
  signal.fillPOT();
  signal.write();

  // nu_e CC background
  CAF bkg1( "/dune/data/users/marshalc/CAFs/mcc11_v2/ND_nue_CCbkg.root" );
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations.size() );
    bool is_wgt = head.isWeightSystematicVariation;
    std::string wgt_var = ( is_wgt ? "wgt" : "var" );
    bkg1.addRWbranch( parIds[i], head.prettyName, wgt_var, head.paramVariations );
    bkg1.iswgt[parIds[i]] = is_wgt;
  }
  bkg1.pot = 0.;
  bkg1.meta_run = 0;
  bkg1.meta_subrun = 0;
  bkg1.version = 2;
  for( int i = 0; i <= 999; ++i ) {
    printf( "nue CC background file %d, so far %4.4g POT\n", i, bkg1.pot );
    TFile * tf = new TFile( Form("/pnfs/dune/persistent/users/marshalc/CAF/genieNuEBkg/FHC/LAr.neutrino.%d.ghep.root",i) );
    if( tf != NULL ) {
      TTree * tree = (TTree*) tf->Get("gtree");
      if( tree && !tf->TestBit(TFile::kRecovered) ) {
        bkg1.pot += 1.0E18;
        loop( tree, 1, bkg1 );
      }
      tf->Close();
    }
    delete tf;
  }
  bkg1.fillPOT();
  bkg1.write();
*/

  // NC background
  CAF bkg2( "/dune/data/users/marshalc/CAFs/mcc11_v2/ND_nue_NCbkg.root" );
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations.size() );
    bool is_wgt = head.isWeightSystematicVariation;
    std::string wgt_var = ( is_wgt ? "wgt" : "var" );
    bkg2.addRWbranch( parIds[i], head.prettyName, wgt_var, head.paramVariations );
    bkg2.iswgt[parIds[i]] = is_wgt;
  }
  bkg2.pot = 0.;
  bkg2.meta_run = 0;
  bkg2.meta_subrun = 0;
  bkg2.version = 2;
  int nevt = 0;
  for( int i = 0; i <= 999; ++i ) {
    printf( "NC background file %d, so far %4.4g POT for %d events\n", i, bkg2.pot, nevt );
    TFile * tf = new TFile( Form("/pnfs/dune/persistent/users/marshalc/CAF/genieNewFluxv2/LAr/FHC/00/LAr.neutrino.%d.ghep.root",i) );
    if( tf != NULL ) {
      TTree * tree = (TTree*) tf->Get("gtree");
      if( tree && !tf->TestBit(TFile::kRecovered) ) {
        bkg2.pot += 5.0E16;
        loop( tree, 2, bkg2 );
        nevt += tree->GetEntries();
      }
      tf->Close();
    }
    delete tf;
  }
  bkg2.fillPOT();
  bkg2.write();
  printf( "Got %g POT for %d events\n", bkg2.pot, nevt );

}

