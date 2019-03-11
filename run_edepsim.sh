#! /usr/bin/env bash

##################################################

HORN=$1
NPER=$2
FIRST=$3
TEST=$4
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

MODE="neutrino"
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
fi

if [ "${NPER}" = "" ]; then
echo "Number of events per job not specified, will run entire input file"
NPER=-1
fi

if [ "${FIRST}" = "" ]; then
echo "First run number not specified, using 0"
FIRST=0
fi

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "Test mode"
#CP="cp"
PROCESS=0
fi

echo "Running edepsim for ${HORN} mode, ${NPER} events"

RNDSEED=$((${PROCESS}+${FIRST}))

GEOMETRY="lar_mpt"

OUTFLAG="LAr"

RDIR=0$((${RNDSEED} / 1000))

USERDIR="/pnfs/dune/persistent/users/marshalc/CAF"
INDIR="/pnfs/dune/persistent/users/marshalc/CAF/genieNewFluxv2/${OUTFLAG}/${HORN}/${RDIR}"
OUTDIR="/pnfs/dune/persistent/users/marshalc/CAF/edepNewFluxv2"

##################################################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof
setup ifdhc

##################################################

## Fetch the input GENIE file and convert it to the rootracker format

echo "${CP} ${INDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root input_file.ghep.root"
${CP} ${INDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root input_file.ghep.root

## Copy a tarball of all the code we need
echo "${CP} ${USERDIR}/edep.tar.gz ${PWD}/edep.tar.gz"
${CP} ${USERDIR}/edep.tar.gz ${PWD}/edep.tar.gz
tar -xzf edep.tar.gz
mv edep/* .
${CP} ${USERDIR}/edep-sim.tar.gz ${PWD}/edep-sim.tar.gz
tar -xzf edep-sim.tar.gz

export LD_LIBRARY_PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}

echo "About to run gntpc, ls:"
ls
echo "LD_LIBRARY_PATH:"
echo ${LD_LIBRARY_PATH}
gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# run all the events in a file
if [ "${NPER}" = -1 ]; then
echo "Specified all events, determining how many events there are"
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)
echo "There are ${NPER} events"
fi

##################################################

## Run edep-sim
ls
echo "edep-sim -C -g ${GEOMETRY}.gdml -o ${PWD}/output_file.root -u -e ${NPER} dune-nd.mac"
edep-sim \
    -C \
    -g ${GEOMETRY}.gdml \
    -o ${PWD}/output_file.root \
    -u \
    -e ${NPER} \
    dune-nd.mac

## Copy output file to dCache persistent
echo "${CP} output_file.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.edepsim.root"
${CP} output_file.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.edepsim.root

##################################################


