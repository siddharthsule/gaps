# Noether Cluster Only
source /gluster/data/theory/ssule/setup-sid.sh

# Copy Git Repo Here and Enter
cp -r ~/gaps .
cd gaps

echo "Running full analysis... (C++ ONLY)"

# Run cuda once so that building is done ;)
source runcode.sh full 64

# Copy Results to Local Machine
cp -r results ~

# Remove Copied Folder (Just in Case!)
cd ..
rm -rf gaps