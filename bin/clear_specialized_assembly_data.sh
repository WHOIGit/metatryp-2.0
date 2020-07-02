#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. $DIR/setEnv.sh

python -m proteomics.scripts.clear_specialized_assembly_data $@
