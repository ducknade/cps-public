#!/bin/bash

# This script is used to automatically install APPROPIATE version of CPS and all
# its dependencies without bfm.
#
# This script would install CPS in $cps and all the dependencies in $prefix.
# Change $cps if you need to install in other directories. These variables are 
# specified in conf.sh.
#
# CAUTION! This script could remove files in $cps silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin. Edited by Jiqun Tu for the KPI for DBW2 action.

. conf.sh

if $build_libs ; then
    if [ -e $prefix ] ; then
        echo "$prefix already exist, continue to build will erase all its contents."
        echo "Use ./scripts/cps.sh to build CPS only."
        echo "Ctrl-C to stop."
        for i in {10..0} ; do
            echo -n "$i "
            sleep 1;
        done
    fi
	rm -rf $prefix || true
	mkdir -p $prefix
fi

if $build_libs ; then
        ./scripts/gsl.sh
        ./scripts/qmp.sh
        # ./scripts/false.sh
        # ./scripts/noweb.sh
        # ./scripts/qa0.sh
#      	./scripts/mdwf.sh
#        ./scripts/fftw.sh
#        ./scripts/cmake.sh
#        ./scripts/eigen.sh
        ./scripts/qio.sh
        # ./scripts/libxml2.sh
        # ./scripts/qdp_pp.sh
        # ./scripts/chroma.sh
        # ./scripts/bagel.sh
fi

if $build_cps ; then
#    ./scripts/timer.sh
#    ./scripts/gaugefield.sh
    time ./scripts/cps.sh
fi

rm -rf $temp_dir
