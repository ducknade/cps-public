#!/bin/bash

. conf.sh

name=cps

echo "!!!! build CPS !!!!"

cps="$(readlink -m $cps)"
echo "Real location is : ""$cps"

rm -rf $cps/build* >/dev/null 2>&1 || true
rm $wd/log >/dev/null 2>&1 || true

mkdir -p $cps || true
cd "$cps"

export TIMEFORMAT="%R secs"

profile () {
    echo
    echo "$@"
    time { "$@" >>$wd/log 2>&1 ; }
}

if [ $arch = bgq ] ; then

    profile rsync -avz --delete $wd/cps_pp/ ./cps_pp/
    cd cps_pp
    profile autoconf
    cd "$cps"
    mkdir build
    cd build
    profile ../cps_pp/configure \
    	--build=powerpc64-unknown-linux-gnu --host=powerpc64-bgq-linux \
	--enable-target=bgq \
        --enable-debug=no \
        --enable-fftw=$prefix \
        --enable-qmp=$prefix \
        --enable-xml=$prefix \
        --enable-qio \
        --enable-gmp \
	--enable-bfm=$prefix/bfm \
        --enable-gsl=$prefix \
	--enable-eigen=$prefix
	profile make -j60
	# profile make
    cd "$cps"

elif [ $arch = amd64 ] ; then

##    profile rsync -avz --delete $wd/cps_pp/ ./cps_pp/
##    cd cps_pp
##    profile autoconf
##    cd "$cps"
##    mkdir build
##    cd build
#    profile ../cps_pp/configure \
#        --enable-fftw=$prefix \
#        --enable-eigen=$prefix \
#        --enable-debug=no
##    profile ../cps_pp/configure --enable-debug=no
#    profile make -j30
##	profile make

	cd "$cps"

    profile rsync -avz --delete $wd/cps_pp/ ./cps_pp-qmp/
    cd cps_pp-qmp
    profile autoconf
    cd "$cps"
    mkdir build-qmp
    cd build-qmp
    profile ../cps_pp-qmp/configure \
        --enable-qmp=$prefix \
        --enable-qio \
        --enable-sse=no \
        --enable-debug=no
	profile make -j30
    cd "$cps"

fi

echo
cd $wd
echo "!!!! CPS build !!!!"

rm -rf $temp_dir
