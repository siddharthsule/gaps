#!/bin/bash

# Run rungaps once, incase compilation is needed
../../rungaps -r gpu
rm gpu.yoda gpu-time.dat

# Define number of samples for tuning
num_samples=20

# Start a timer
start_time=$(date +%s)

# Sample
prof2-sample limits.txt -t template.sh -n $num_samples --seed 1234

# Run the tuning jobs
for i in $(seq 0 $((num_samples-1))); do
    cd scan/$(printf "%04d" $i)
    source template.sh
    cd ../../
done

# Interpolate
prof2-ipol scan/ --wfile weights.txt --order 3

# Tune
prof2-tune ipol.dat -d data/ --wfile weights.txt -r scan/ --filter

# End timer and print elapsed time
end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total tuning time: $elapsed seconds"