#! /usr/bin/env bash

##################################################

HORN=$1
FIRST=$2
TEST=$3
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

MODE="neutrino"
RHC=""
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
RHC=" --rhc"
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

OUTFLAG="GAr"

RDIR=0$((${RNDSEED} / 1000))

USERDIR="/pnfs/dune/persistent/users/marshalc/CAF"
INDIR="/pnfs/dune/persistent/users/sbjones/CAF/genie/${OUTFLAG}/${HORN}/${RDIR}"
DUMPDIR="/pnfs/dune/persistent/users/sbjones/CAF/dump"
CAFDIR="/pnfs/dune/persistent/users/sbjones/CAF/outCAFs/v10"
STUFF="/pnfs/dune/persistent/users/sbjones/CAF/DUNE_ND_CAF.tar.gz"


##################################################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof
setup ifdhc

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

##################################################

## Fetch the input GENIE file and convert it to the rootracker format

echo "${CP} ${INDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root input_file.ghep.root"
${CP} ${INDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root input_file.ghep.root

## Copy a tarball of all the code we need
echo "${CP} ${USERDIR}/edep.tar.gz ${PWD}/edep.tar.gz"
${CP} ${USERDIR}/edep.tar.gz ${PWD}/edep.tar.gz
tar -xzf edep.tar.gz
mv edep/* .

${CP} ${STUFF} DUNE_ND_CAF.tar.gz
tar xzf DUNE_ND_CAF.tar.gz
mv DUNE_ND_CAF/* .

#${CP} ${USERDIR}/edep-sim.tar.gz ${PWD}/edep-sim.tar.gz
#tar -xzf edep-sim.tar.gz

export LD_LIBRARY_PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}

gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# run all the events in a file
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)
echo "There are ${NPER} events"

##################################################

## Run edep-sim
ls
echo "edep-sim -C -g ${GEOMETRY}.gdml -o ${PWD}/edep.${RNDSEED}.root -u -e ${NPER} dune-nd.mac"
edep-sim \
    -C \
    -g ${GEOMETRY}.gdml \
    -o ${PWD}/edep.${RNDSEED}.root \
    -u \
    -e ${NPER} \
    dune-nd.mac

## Don't copy output file to dCache persistent
#echo "${CP} output_file.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.edepsim.root"
#${CP} output_file.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.edepsim.root

##################################################

# makeCAF expects ghep file to be genie.RUN.root
mv input_file.ghep.root genie.${RNDSEED}.root

# hack
unset LD_LIBRARY_PATH
setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof
setup ifdhc
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/nusystematics/build/Linux/lib:${PWD}/nusyst/artless

## Run dumpTree
echo "Running dumpTree.py..."
echo "python dumpTree.py --topdir ${PWD} --first_run ${RNDSEED} --last_run ${RNDSEED} ${RHC} --grid --outfile dump.root"
python dumpTree.py --topdir ${PWD} --first_run ${RNDSEED} --last_run ${RNDSEED} ${RHC} --grid --outfile dump.root

## Run makeCAF
echo "Running makeCAF..."
echo "./makeCAF --edepfile dump.root --ghepdir ${PWD} --outfile CAF.root --fhicl fhicl.fcl --seed ${RNDSEED} --grid ${RHC}"
./makeCAF --edepfile dump.root --ghepdir ${PWD} --outfile CAF.root --fhicl ./fhicl.fcl --seed ${RNDSEED} ${RHC} --grid

## copy outputs
echo "Copying outputs..."
echo "${CP} dump.root ${DUMPDIR}/${HORN}_${RNDSEED}.root"
${CP} dump.root ${DUMPDIR}/${HORN}_${RNDSEED}.root
echo "${CP} CAF.root ${CAFDIR}/CAF_${HORN}_${RNDSEED}.root"
${CP} CAF.root ${CAFDIR}/CAF_${HORN}_${RNDSEED}.root




