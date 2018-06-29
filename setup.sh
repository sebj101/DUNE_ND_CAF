source /grid/fermiapp/products/dune/setup_dune.sh
# cmake v3 is required for edep-sim
setup cmake v3_9_0
# gcc v6 is required for c++14, used in nusyst
setup gcc v6_4_0

# setup ROOT
source /cvmfs/larsoft.opensciencegrid.org/products/root/v6_10_08b/Linux64bit+2.6-2.12-e15-nu-prof/bin/thisroot.sh

setup genie        v2_12_8c   -q e15:prof
setup genie_xsec   v2_12_8    -q DefaultPlusMECWithNC
setup genie_phyopt v2_12_8    -q dkcharmtau
setup dk2nu        v01_05_01b -q e15:prof
setup geant4 v4_10_3_p01b -q e15:prof
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

# change tab title in this window
echo -ne "\033]0;nusyst  --  ${HOSTNAME%%.*}"; echo -ne "\007"

cd edep-sim
source setup.sh
cd ..

export NUSYST=`pwd`/nusyst
export LD_LIBRARY_PATH=$NUSYST/build/Linux/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NUSYST/build/nusyst/artless:$LD_LIBRARY_PATH
