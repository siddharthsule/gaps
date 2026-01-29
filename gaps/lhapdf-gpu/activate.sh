#!/bin/bash

# Set up install location
export INSTALL_LOC=$PWD/install
echo "INSTALL_LOC is set to: $INSTALL_LOC"

# Set up paths
export PATH=$PATH:$INSTALL_LOC/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_LOC/lib
export PYTHONPATH=$PYTHONPATH:$INSTALL_LOC/lib/python3.9/site-packages

# Install CT14lo set
lhapdf install CT14lo
echo "LHAPDF ready to use!"
