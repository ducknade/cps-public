#!/bin/bash

. conf.sh

name=qmp

echo "!!!! build qmp !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
rsync -avz --delete $distfiles/$name .

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir
cd $src_dir/qmp*/
autoconf
cd $build_dir

if [ $arch = bgq ] ; then
    $src_dir/qmp*/configure \
        --prefix=$prefix \
        CXX=mpicxx CC=mpicc \
        --host=powerpc64-bgq-linux \
        --build=none \
        --with-qmp-comms-type=mpi \
        --enable-bgq
else
    $src_dir/qmp*/configure \
        --prefix=$prefix \
        --with-qmp-comms-type=mpi \
        CXX=mpicxx CC=mpicc
fi

make
make install

cd $wd
echo "!!!! qmp build !!!!"

rm -rf $temp_dir
