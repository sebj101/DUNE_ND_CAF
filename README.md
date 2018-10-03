# DUNE_ND_CAF
Make DUNE near detector CAFs

% git clone https://github.com/cmmarshall/DUNE_ND_CAF.git
% cd DUNE_ND_CAF
% export NDCAF=`pwd`
% git clone https://github.com/luketpickering/nusystematics.git

Set up the environment using setup_nd_cafmaker.sh
% source setup_nd_cafmaker.sh

Build nusystematics with cmake (which will build dependent systematicstools)
% cd nusystematics
% mkdir build; cd build
% cmake ../ -DUSEART=0 -DLIBXML2_LIB=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/lib/ -DLIBXML2_INC=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/include/libxml2 -DPYTHIA6=/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28i/Linux64bit+2.6-2.12-gcc640-prof/lib
% make
% make install

Add nusystematics stuff to LD_LIBRARY_PATH
Note: this should be done every time, so you probably want to add these lines to a setup script
% export NUSYST=$NDCAF/nusystematics
% export LD_LIBRARY_PATH=$NUSYST/build/Linux/lib:$LD_LIBRARY_PATH
% export LD_LIBRARY_PATH=$NUSYST/build/nusystematics/artless:$LD_LIBRARY_PATH

You also need edep-sim to run the Geant4 stage, and also to make the intermediate TTree that is used by makeCAF.
Install edep-sim by following the instructions in the README found here:
https://github.com/ClarkMcGrew/edep-sim
There is a build script that just works in my experience.

