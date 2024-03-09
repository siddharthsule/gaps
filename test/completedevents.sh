# Algorithm to get all the completed events data
# Not really worth putting in gaps.sh, so in Test

# Move to the root directory
cd ..

# Remove the old and make the new results directory
rm -rf results
mkdir -p results

# Run for 1000, 10000, 100000 events and get results
for i in 1000 10000 100000
do
    echo "Running for $i events"
    # Run the algorithm
    ./bin/gaps $i
    mv gaps-cycles.dat results/gaps-cycles-$i.dat
done

# For 1000000 events, profile and store the results
echo "Running for 1000000 events"
# Run the algorithm
nsys profile --stats=true ./bin/gaps 1000000 > profile.log
mv gaps-cycles.dat results/gaps-cycles-1000000.dat
mv profile.log results/profile.log