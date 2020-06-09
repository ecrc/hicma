#!/bin/bash -le
module purge

if [ "$HOSTNAME" == "thana" ]; then
	. ./scripts/power8.modules
else
	. ./scripts/intel.modules
fi

export MORSE_TESTING_VERBOSE=1
./build/testing/testing_zposv 1 0 posv 1000 1000 200 1000 250 1e-7 0 250 250 1 1 1
