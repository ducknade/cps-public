#!/bin/bash

. conf.sh

name=qio

echo "!!!! build qio !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xf $distfiles/qio-2.3.8.tar.gz

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir

if [ $arch = bgq ] ; then
    $src_dir/qio-2.3.8/configure \
        --prefix=$prefix \
        --host=powerpc64-bgq-linux \
        --build=none \
        --with-qmp=$prefix \
        --enable-largefile CFLAGS=-O2
else
    $src_dir/qio-2.3.8/configure \
        --prefix=$prefix \
        --build=none \
        --with-qmp=$prefix \
        --enable-largefile CFLAGS=-O2
fi

make -j60
make install

cd $wd
echo "!!!! qio build !!!!"

rm -rf $temp_dir
