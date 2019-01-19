#!/bin/bash

. conf.sh

name=cps

echo "!!!! building CPS !!!!"

cps="$(cd $cps; pwd)"
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
cp cps.a libcps.a
cd "$cps"

echo
cd $wd
echo "!!!! CPS built !!!!"

rm -rf $temp_dir
