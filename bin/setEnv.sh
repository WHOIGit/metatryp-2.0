#!/bin/bash


# Get the proteomics base dir.
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
echo $BASE_DIR
# Activate virtualenv
#conda activate $BASE_DIR/metatryp2_env/

# Set python path to include proteomics
export PYTHONPATH="$PYTHONPATH:$BASE_DIR/lib"
echo $PYTHONPATH


