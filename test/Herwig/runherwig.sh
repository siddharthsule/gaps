#!/bin/bash

# Number of events and cores, and preamble for mac docker
nevents=${1:-100000}
ncores=${2:-5}
option=${3:-"all"}

# Build Analysis
rivet-build RivetGAPS_LEP.so GAPS_LEP.cc
rivet-build RivetGAPS_LHC.so GAPS_LHC.cc
export RIVET_ANALYSIS_PATH=$PWD

# Array of configurations
if [ "$option" == "all" ]; then
  configs=("LEP" "LHC" "LEPNLO" "LHCNLO")
else
  configs=("$option")
fi

# Loop through each configuration
for config in "${configs[@]}"; do

  # Remove Herwig cache folder if exists
  if [ -d "Herwig-cache" ]; then
    rm -rf Herwig-cache
  fi

  # Run Herwig
  Herwig read ${config}.in
  Herwig run ${config}.run -N $nevents -j $ncores
  rivet-merge ${config}-*.yoda -o ${config}.yoda -e

  # Rename Rivet analyses - lowercase
  sed -i "s/GAPS_LEP/gaps_lep/g" ${config}.yoda
  sed -i "s/GAPS_LHC/gaps_lhc/g" ${config}.yoda

  # Fix any small values like 1e-16 to 0
  python yoda-fixer.py ${config}.yoda

  # Rename Yoda files - lowercase
  mv ${config}.yoda herwig-${config,,}.yoda

  # Cleanup
  rm ${config}-*.yoda ${config}-*.log ${config}-*.out ${config}-*.tex
  rm ${config}.run
done

# Make Plots with Rivet
# rivet-mkhtml --mc-errs -c plots.conf *.yoda