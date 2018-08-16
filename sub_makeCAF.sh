#! /usr/bin/env bash

# This script is intended for submitting ND CAF-maker to the Fermilab grid
# It first runs dumpTree (extract edep-sim output to flat tree)
# Then runs makeCAF (create CAF file from extracted output)
# The intermediate flat tree is also saved
# Syntax for the jobsub_submit command is:
# jobsub_submit --group dune --role=Analysis -N 100 --OS=SL6 --expected-lifetime=12h --memory=4000MB --group=dune file://`pwd`/sub_makeCAF.sh
##################################################

HORN=$1
NPER=$2
TEST=$3
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

if [ "${NPER}" = "" ]; then
echo "Number of runs per file not specified, using 50"
NPER=50
fi

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "In TEST mode, assuming interactive running"
CP="cp"
PROCESS=0
fi


FIRSTRUN=$((PROCESS * NPER))
LASTRUN=$((PROCESS * NPER + NPER))

NEUTRINO="neutrino"
RHC=""
if [ "${HORN}" = "RHC" ]; then
echo "Using RHC beam mode"
NEUTRINO="antineutrino"
RHC=" --rhc"
fi

INPUTTOP="/pnfs/dune/persistent/users/jmalbos/GArTPC_ND/data"
DUMPDIR="/pnfs/dune/persistent/users/marshalc/CAF/dump"
CAFDIR="/pnfs/dune/persistent/users/marshalc/CAF/CAF"
#NUSYST="/pnfs/dune/persistent/users/marshalc/CAF/nusyst.tar.gz"
#EDEPSIM="/pnfs/dune/persistent/users/marshalc/CAF/edep-sim.tar.gz"
#EDEPSIMEVENTS="/pnfs/dune/persistent/users/marshalc/CAF/EDepSimEvents.tar.gz"
#DUMPTREE="/pnfs/dune/persistent/users/marshalc/CAF/dumpTree.py"
#MAKECAF="/pnfs/dune/persistent/users/marshalc/CAF/makeCAF"
#FHICL="/pnfs/dune/persistent/users/marshalc/CAF/fhicl.fcl"
STUFF="/pnfs/dune/persistent/users/marshalc/CAF/DUNE_ND_CAF.tar.gz"

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
echo "Copying input files for runs ${FIRSTRUN} to ${LASTRUN}..."
for RUN in $(seq ${FIRSTRUN} ${LASTRUN})
do
  RDIR=0$((${RUN} / 1000))
  echo "Copying: ${INPUTTOP}/sim/LArDipole/${RDIR}/LArDipole.${NEUTRINO}.${RUN}.edepsim.root"
  ${CP} ${INPUTTOP}/sim/LArDipole/${RDIR}/LArDipole.${NEUTRINO}.${RUN}.edepsim.root edep.${RUN}.root
  echo "Copying: ${INPUTTOP}/genie/LAr/${RDIR}/LAr.${NEUTRINO}.${RUN}.ghep.root"
  ${CP} ${INPUTTOP}/genie/LAr/${RDIR}/LAr.${NEUTRINO}.${RUN}.ghep.root genie.${RUN}.root
done

##################################################
## Copy edep-sim and nusyst binaries and untar them

echo "Getting edep-sim and nusyst code"
${CP} ${STUFF} DUNE_ND_CAF.tar.gz
tar xzf DUNE_ND_CAF.tar.gz
mv DUNE_ND_CAF/* .
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/edep-sim/lib

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/nusyst/build/Linux/lib:${PWD}/nusyst/artless

## Run dumpTree
echo "Running dumpTree.py..."
echo "python dumpTree.py --topdir ${PWD} --first_run ${FIRSTRUN} --last_run $((LASTRUN-1)) ${RHC} --grid --outfile dump.root"
python dumpTree.py --topdir ${PWD} --first_run ${FIRSTRUN} --last_run $((LASTRUN-1)) ${RHC} --grid --outfile dump.root

## Run makeCAF
echo "Running makeCAF..."
echo "./makeCAF --edepfile dump.root --ghepdir ${PWD} --outfile CAF.root --fhicl fhicl.fcl --seed ${PROCESS} --grid ${RHC}"
./makeCAF --edepfile dump.root --ghepdir ${PWD} --outfile CAF.root --fhicl fhicl.fcl --seed ${PROCESS} ${RHC} --grid

## copy outputs
echo "Copying outputs..."
${CP} dump.root ${DUMPDIR}/${HORN}_${PROCESS}.root
${CP} CAF.root ${CAFDIR}/CAF_${HORN}_${PROCESS}.root

##################################################

