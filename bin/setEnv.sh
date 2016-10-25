#!/bin/bash


# Get the proteomics base dir.
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
echo $BASE_DIR
# Activate virtualenv
source $BASE_DIR/py2.7/bin/activate

# Set python path to include proteomics
export PYTHONPATH="$PYTHONPATH:$BASE_DIR/lib"
echo $PYTHONPATH
# Set path to database
export PROTEOMICS_DB="postgresql://metatryp:Pr0t3om1cs@localhost:5432/proteomics"
export PROTEOMICS_DB_NAME="proteomics"

