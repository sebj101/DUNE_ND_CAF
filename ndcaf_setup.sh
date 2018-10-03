source /grid/fermiapp/products/dune/setup_dune.sh
# cmake v3 is required for edep-sim
setup cmake v3_9_0
# gcc v6 is required for c++14, used in nusyst
setup gcc v6_4_0
# needed for grid submission
setup jobsub_client
setup pycurl

# flux file format
setup dk2nu        v01_05_01b -q e15:prof

# GENIE
setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau

# Geant4 for edep-sim
setup geant4 v4_10_3_p01b -q e15:prof

# copy files to and from grid nodes
setup ifdhc

# we need something called TBB for ROOT
export TBB_DIR=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1
export TBB_LIB=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/lib
export TBB_INC=/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_1/Linux64bit+2.6-2.12-e15-prof/include

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# After building nusystematics, add these lines
# export NUSYST=`pwd`/nusystematics
# export LD_LIBRARY_PATH=$NUSYST/build/Linux/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=$NUSYST/build/nusystematics/artless:$LD_LIBRARY_PATH
