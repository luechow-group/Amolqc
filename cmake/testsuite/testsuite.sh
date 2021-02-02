#! /bin/bash

# to make exit status non-zero if tests.sh returns non-zero
set -e
set -o pipefail

export AMOLQC_TEST=$1
echo AMOLQC = $AMOLQC_TEST
export AMOLQCRUN=$2
echo AMOLQCRUN = $AMOLQCRUN
cd testsuite/run;bash ./tests.sh $3 $4 | tee test.log
