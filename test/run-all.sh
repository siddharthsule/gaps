#!/bin/bash

# Number of Events and Number of CPU Cores
nevents=1000000
ncores=4

# Head to main gaps dir
cd ..

# Clean up any leftover cluster output from previous runs so it doesn't
# get swept into the first yodamerge below
rm -f cpu-*.yoda

# LEP LO
./rungaps -n $nevents
mv gpu.yoda gpu_lep_lo.yoda
./rungaps -n $nevents -r cpu-cluster -ncpu $ncores
yodamerge -o cpu_lep_lo.yoda cpu-*.yoda; rm cpu-*.yoda

# LEP NLO
./rungaps -nlo -n $nevents
mv gpu.yoda gpu_lep_nlo.yoda
./rungaps -nlo -n $nevents -r cpu-cluster -ncpu $ncores
yodamerge -o cpu_lep_nlo.yoda cpu-*.yoda; rm cpu-*.yoda

# LEP NLO Hadronised
./rungaps -nlo -cmw -hadronise -n $nevents
mv gpu.yoda gpu_lep_nlo_had.yoda
./rungaps -nlo -cmw -hadronise -n $nevents -r cpu-cluster -ncpu $ncores
yodamerge -o cpu_lep_nlo_had.yoda cpu-*.yoda; rm cpu-*.yoda

# LHC LO
./rungaps -p LHC -n $nevents
mv gpu.yoda gpu_lhc_lo.yoda
./rungaps -p LHC -n $nevents -r cpu-cluster -ncpu $ncores
yodamerge -o cpu_lhc_lo.yoda cpu-*.yoda; rm cpu-*.yoda

# LHC NLO
./rungaps -p LHC -nlo -n $nevents
mv gpu.yoda gpu_lhc_nlo.yoda
./rungaps -p LHC -nlo -n $nevents -r cpu-cluster -ncpu $ncores
yodamerge -o cpu_lhc_nlo.yoda cpu-*.yoda; rm cpu-*.yoda

# return to testdir
cd test
