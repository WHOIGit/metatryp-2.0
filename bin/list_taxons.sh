#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. $DIR/setEnv.sh

psql -U protein -d $PROTEOMICS_DB_NAME -c "SELECT * FROM taxon;"

