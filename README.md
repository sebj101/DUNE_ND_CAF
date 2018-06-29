# DUNE_ND_CAF
Make DUNE near detector CAFs

Requires edep-sim Geant4 package (https://github.com/ClarkMcGrew/edep-sim), and DUNE reweighting framework nusyst (https://github.com/luketpickering/nusyst). Both are easy to build at Fermilab. Setup script uses Fermilab ups, version for LBL NERSC forthcoming.

Usage:

dumpTree.py
Reads in edep-sim output file and produces flat ROOT tree with information necessary to make CAFs.

python dumpTree.py --topdir /path/to/edep-sim/output/files \
                   --outfile /path/to/output/filename.root \
                   --first_run X \
                   --last_run Y
Note it assumes a particular format of edep-sim files, which you can see from looking at the code. Will look for typeinfo for edep-sim classes and create it in a directory called EDepSimEvents using TFile::MakeProject. There is an optional flag --rhc to run in RHC mode

makeCAF
Compile with the Makefile provided, which requires some environment variables that are created by the setup script. 

./makeCAF --edepfile /path/to/output/from/dumpTree/file.root \
          --ghepdir /path/to/genie/ghep/files \
          --outfile /path/to/CAF/file/CAF.root \
          --fhicl /path/to/fhicl/file/with/DUNErw
Optionally you can set the random number seed with --seed for reproduceability, and --rhc for RHC mode. The fhicl file is what configures DUNErw shifts, using nusyst. Branches are generated automatically based on the configuration in that file.
