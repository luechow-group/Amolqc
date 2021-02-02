#! /bin/bash

# to make exit status non-zero if tests.sh returns non-zero
set -e
set -o pipefail

cd testsuite/run;bash ./tests.sh $1 $2 | tee test.log
