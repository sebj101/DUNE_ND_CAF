#! /usr/bin/env bash

##################################################

FIRSTRUN=$((${PROCESS} * 100))
LASTRUN=$((${PROCESS} * 100 + 99))
RDIR=0$((${PROCESS} / 10))

RNDSEED=${PROCESS}
HORN="FHC"
NEUTRINO="neutrino"
RHC=""
CP="ifdh cp"

INPUTTOP="/pnfs/dune/persistent/users/jmalbos/GArTPC_ND/data"
DUMPDIR="/pnfs/dune/persistent/users/marshalc/CAF/dump"
CAFDIR="/pnfs/dune/persistent/users/marshalc/CAF/CAF"
NUSYST="/pnfs/dune/persistent/users/marshalc/CAF/nusyst.tar.gz"
EDEPSIM="/pnfs/dune/persistent/users/marshalc/CAF/edep-sim.tar.gz"
EDEPSIMEVENTS="/pnfs/dune/persistent/users/marshalc/CAF/EDepSimEvents.tar.gz"
DUMPTREE="/pnfs/dune/persistent/users/marshalc/CAF/dumpTree.py"
MAKECAF="/pnfs/dune/persistent/users/marshalc/CAF/makeCAF"
FHICL="/pnfs/dune/persistent/users/marshalc/CAF/fhicl.fcl"

##################################################

## Setup UPS and required products
echo "Setting up software..."

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

source /cvmfs/larsoft.opensciencegrid.org/products/root/v6_10_08b/Linux64bit+2.6-2.12-e15-nu-prof/bin/thisroot.sh
setup genie        v2_12_8c   -q e15:prof
setup genie_xsec   v2_12_8    -q DefaultPlusMECWithNC
setup genie_phyopt v2_12_8    -q dkcharmtau
setup dk2nu        v01_05_01b -q e15:prof
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

##################################################
## Fetch the genie and edep-sim output files
echo "Copying input files..."
for RUN in $(seq ${FIRSTRUN} ${LASTRUN})
do
  echo "${INPUTTOP}/sim/LArDipole/${RDIR}/LArDipole.${NEUTRINO}.${RUN}.edepsim.root"
  echo "${INPUTTOP}/genie/LAr/${RDIR}/LAr.${NEUTRINO}.${RUN}.ghep.root"
  ${CP} ${INPUTTOP}/sim/LArDipole/${RDIR}/LArDipole.${NEUTRINO}.${RUN}.edepsim.root edep.${RUN}.root
  ${CP} ${INPUTTOP}/genie/LAr/${RDIR}/LAr.${NEUTRINO}.${RUN}.ghep.root genie.${RUN}.root
done

##################################################
## Copy edep-sim and nusyst binaries and untar them

echo "Getting edep-sim and nusyst code"
${CP} ${EDEPSIM} edep-sim.tar.gz
tar xzf edep-sim.tar.gz
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/edep-sim/lib

${CP} ${NUSYST} nusyst.tar.gz
tar xzf nusyst.tar.gz
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/nusyst/build/Linux/lib:${PWD}/nusyst/artless

${CP} ${EDEPSIMEVENTS} EDepSimEvents.tar.gz
tar xzf EDepSimEvents.tar.gz

## copy dumpTree.py, makeCAF binary, fhicl file
${CP} ${DUMPTREE} dumpTree.py
${CP} ${MAKECAF} makeCAF
${CP} ${FHICL} fhicl.fcl

## Run dumpTree
echo "Running dumpTree.py"
python dumpTree.py --topdir ${PWD} --first_run ${FIRSTRUN} --last_run ${LASTRUN} ${RHC} --grid --outfile dump.root

## Run makeCAF
echo "Running makeCAF"
./makeCAF --edepfile dump.root --ghepdir ${PWD} --outfile CAF.root --fhicl fhicl.fcl --seed ${JOBNO} ${RHC} --grid

## copy outputs
${CP} dump.root ${DUMPDIR}/${HORN}_${JOBNO}.root
${CP} CAF.root ${CAFDIR}/CAF_${HORN}_${JOBNO}.root

##################################################


