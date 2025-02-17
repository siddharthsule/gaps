#------------------------------------------------------------------------------
#!/bin/bash
# ------------------------------------------------------------------------------

# GAPS - Run Script
# -----------------
# This script is used to compile and run the GAPS and C++ Shower codes. It
# provides a number of options to control the number of events, the number of
# cores to use, and the type of run to perform. The run types are:
#   - gaps: Run the GAPS simulation
#   - cpp: Run the C++ Shower simulation
#   - compare: Run both the GAPS and C++ Shower and compare the results
#   - full: Run both the GAPS and C++ Shower for a range of event numbers


# ------------------------------------------------------------------------------
# Default Parameters

events=10000
energy=91.2
ncores=1

# Defaults to just the CUDA simulation
runtype="gaps"

# ------------------------------------------------------------------------------
# Use optargs to adjust default values and define run

while getopts "n:e:c:r:h" opt; do
    case $opt in
        n) events=$OPTARG ;;
        e) energy=$OPTARG ;;
        c) ncores=$OPTARG ;;
        r) runtype=$OPTARG ;;
        h) echo "Usage: $0 [-n nevents] [-e energy] [-c cores] [-r runtype] [-h help]"
        echo "  -n: set the number of events (default: 10000)"
        echo "  -e: set the CoM energy of the system (default: 91.2)"
        echo "  -c: set the number of cores (default: 1)"
        echo "  -r: set the run type (default: gaps, options: gaps, cpp, compare, full)"
        echo "  -h: display this help message"
        exit 0
        ;;
    esac
done

# ------------------------------------------------------------------------------
# Compile Code

compile() {
    dir=$1
    echo "Compiling $dir"
    (cd $dir && mkdir -p build && cd build && cmake .. && make -j $ncores)
}

if [ "$runtype" = "gaps" ] || [ "$runtype" = "compare" ] || [ "$runtype" = "full" ]; then
    compile "gaps"
fi

if [ "$runtype" = "cpp" ] || [ "$runtype" = "compare" ] || [ "$runtype" = "full" ]; then
    compile "cpp-shower"
fi 

# ------------------------------------------------------------------------------
# Default: Just run GAPS

if [ "$runtype" = "gaps" ]; then
    echo "Running GAPS"
    ./gaps/bin/gaps $events $energy
fi

# ------------------------------------------------------------------------------
# Run C++ Shower

if [ "$runtype" = "cpp" ]; then
    echo "Running C++ Shower"
    ./cpp-shower/bin/cpp-shower $events $energy
fi

# ------------------------------------------------------------------------------
# Compare GAPS and C++ Shower

if [ "$runtype" = "compare" ]; then
    echo "Running GAPS"
    ./gaps/bin/gaps $events $energy
    echo "Running C++ Shower"
    ./cpp-shower/bin/cpp-shower $events $energy
fi

# ------------------------------------------------------------------------------
# Do a full analysis, comparing GAPS and C++ Shower at many values of N

if [ "$runtype" = "full" ]; then

    # Remove previous results and make folder
    rm -rf results
    mkdir -p results

    # Clear previous log files
    rm cpp-time.dat
    rm gaps-time.dat

    # Run the comparison 100 times, for different number of events
    neventsarray=(1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000)
    for n in "${neventsarray[@]}"
    do
        # Run and store the output in a log file
        for i in {1..100}
        do
            echo "Running GAPS with $n events"
            ./gaps/bin/gaps $n $energy
            echo "Running C++ Shower with $n events"
            ./cpp-shower/bin/cpp-shower $n $energy
        done
    done

    # Move the log files to the results directory
    mv cpp-time.dat results/
    mv gaps-time.dat results/
    mv cpp.yoda results/
    mv gaps.yoda results/

fi

