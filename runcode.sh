#!/bin/bash
# ------------------------------------------------------------------------------
nevents=${2:-10000} # Number of Events to Simulate
ncores=${3:-1} # Number of Cores to build both Generators
# ------------------------------------------------------------------------------
# Function to build the C++ code
build_cpp() {
    local ncores=$1
    echo "Building C++ code..."
    cd cpp-shower
    mkdir -p build
    cd build
    cmake ..
    make -j $ncores
    cd ../..
}
# ------------------------------------------------------------------------------
# Function to build the CUDA code
build_cuda() {
    local ncores=$1
    echo "Building CUDA code..."
    mkdir -p build
    cd build
    cmake ..
    make -j $ncores
    cd ..
}
# ------------------------------------------------------------------------------
# Function to run C++ simulations - (Single Core Only as of 16-02-2024)
run_cpp() {
    local nevents=$1
    local ncores=$2
    echo "Running C++ simulations..."
    build_cpp $ncores
    ./cpp-shower/bin/main $nevents
}
# ------------------------------------------------------------------------------
# Function to run CUDA simulations
run_cuda() {
    local nevents=$1
    local ncores=$2
    echo "Running CUDA simulations..."
    build_cuda $ncores
    ./bin/gaps $nevents
}
# ------------------------------------------------------------------------------
# Function to compare CUDA and C++ for a set of events
compare() {
    local nevents=$1
    local ncores=$2
    echo "Comparing CUDA and C++..."
    run_cpp $nevents $ncores
    run_cuda $nevents $ncores
}
# ------------------------------------------------------------------------------
# Function to do a full analysis
full_analysis() {
    local ncores=$1
    echo "Running full analysis..."
    rm -rf results
    mkdir -p results
    mkdir -p results/cycles
    # Clear previous log files
    rm cpp-time.dat
    rm cuda-time.dat
    # Run the comparison 10 times, for different number of events
    neventsarray=(1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000)
        for n in "${neventsarray[@]}"
        do
            # Run and store the output in a log file
            for i in {1..10}
            do
                compare $n $ncores
                mv cuda-cycles.dat cuda-cycles-$n-$i.dat
            done
        done
    mv cpp-time.dat results/
    mv cuda-time.dat results/
    mv cpp-shower.yoda results/
    mv cuda-shower.yoda results/
    mv cuda-cycles-*.dat results/cycles/
}
# ------------------------------------------------------------------------------
# Function to display help
help() {
    echo "Usage: $0 {cuda|cpp|compare|full} [nevents] [ncores]"
    echo "Defaults: nevents=10000, ncores=1"
    echo ""
    echo "  cuda: Run CUDA simulations (builds CUDA code before running)"
    echo "    Usage: $0 cuda [nevents]"
    echo "    Example: $0 cuda 5000"
    echo ""
    echo "  cpp: Run C++ simulations (builds C++ code before running)"
    echo "    Usage: $0 cpp [nevents] [ncores]"
    echo "    Example: $0 cpp 5000 4"
    echo ""
    echo "  compare: Compare CUDA and C++ (builds both CUDA and C++ code before running)"
    echo "    Usage: $0 compare [nevents] [ncores]"
    echo "    Example: $0 compare 5000 4"
    echo ""
    echo "  full: Run full analysis (builds both CUDA and C++ code before running)"
    echo "    Usage: $0 full [ncores]"
    echo "    Example: $0 full 4"
    echo ""
    echo "  help: Display this help message"
}

# Check the first command-line argument
case $1 in
    cuda)
        run_cuda $nevents $ncores
        ;;
    cpp)
        run_cpp $nevents $ncores
        ;;
    compare)
        compare $nevents $ncores
        ;;
    full)
        # no nevents so ncores is the first argument
        full_analysis ${2:-1}
        ;;
    help)
        help
        ;;
    *)
        echo "Usage: $0 {cuda|cpp|compare|full} <nevents> <ncores>"
        echo "Run $0 help for more information"
esac
