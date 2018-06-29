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
  int run, isFD, isFHC;

  // Truth information that does not get reweighted
  int ccnc, beamPdg, neu, LepPDG, mode, nipi0, nipip, nipim;
  double Ev, Elep, Q2;

  // Reco information CV
  double Ev_reco, Elep_reco, numu_pid, nue_pid, mvaresult;
  int reco_q;

  // reweights -- make sure big enough to hold all the variations for each knob, and all the knobs
  int nwgt[100];
  double wgt[100][20];

  // meta
  double pot;
  int meta_run;

  TFile * cafFile;
  TDirectory * mvaselect;
  TTree * cafMVA;
  TTree * cafPOT;
};

#endif

