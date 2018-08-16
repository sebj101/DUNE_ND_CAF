#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"

class CAF {

public:
  CAF( std::string filename );
  ~CAF();
  void fill();
  void fillPOT();
  void write();
  void addRWbranch( int parId, std::string name, std::vector<double> &vars );
  void Print();

  // Make ntuple variables public so they can be set from other file

  // configuration variables
  int isFD, isFHC;
  // event accounting
  int run, subrun, event;
  // Truth information
  int isCC, neutrinoPDG, mode, LepPDG; 
  double Ev, Q2, W, X, Y, NuMomX, NuMomY, NuMomZ, LepMomX, LepMomY, LepMomZ, LepE, LepNuAngle;
  // True particle counts
  int nP, nN, nipip, nipim, nipi0, nikp, nikm, nik0, niem, niother, nNucleus, nUNKNOWN;

  // Reco information CV
  double Ev_reco, Elep_reco;
  int reco_numu, reco_nue, reco_nc, reco_q;
  int muon_contained, muon_tracker, muon_ecal, muon_exit;
  double Ehad_veto;

  // reweights -- make sure big enough to hold all the variations for each knob, and all the knobs
  // the names, and what they actually mean, are determined automatically from the fhicl input file
  int nwgt[100];
  double wgt[100][100];

  // meta
  double pot;
  int meta_run, meta_subrun;

  TFile * cafFile;
  TDirectory * mvaselect;
  TTree * cafMVA;
  TTree * cafPOT;
};

#endif

