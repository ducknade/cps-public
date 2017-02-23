#!/bin/bash

if [ -z "$cps" ] ; then
    cps=$HOME/cps-jtu-build/public/
fi

prefix=$cps/local

if [ -z "$build_cps" ] ; then
    build_cps=true
fi
if [ -z "$build_libs" ] ; then
    build_libs=true
fi

set -e

if [ $(uname -m) = ppc64 ] ; then
    arch=bgq
elif [ $(uname -m) = x86_64 ] ; then
    arch=amd64
else
    arch=unknown
fi

if [ $arch = bgq ] ; then
    if [ "$host" = "DD1" ] ; then
        # export PATH=/bgsys/drivers/DRV2012_0130_2243/ppc64-rhel60/gnu-linux/bin:$PATH
        # export PATH=/bgsys/drivers/DRV2012_0130_2243/ppc64-rhel60/comm/gcc/bin:$PATH
        export PATH=/bgsys/drivers/ppcfloor/gnu-linux/bin:$PATH
        export PATH=/bgsys/drivers/ppcfloor/comm/gcc/bin:$PATH
    else
        export PATH=/bgsys/drivers/ppcfloor/gnu-linux/bin:$PATH
        export PATH=/bgsys/drivers/ppcfloor/comm/gcc/bin:$PATH
    fi
fi

add-to-colon-list () {
    local name="$1"
    local new_value="$2"
    local value="${!name}"
    if [ -z "$value" ] ; then
        export "$name"="$new_value"
    else
        export "$name"="$new_value":"$value"
    fi
}

add-to-colon-list PATH "$prefix/bin"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib:$prefix/lib64"
add-to-colon-list LIBRARY_PATH "$prefix/lib:$prefix/lib64"
add-to-colon-list C_INCLUDE_PATH "$prefix/include:$prefix/include/ncurses"
add-to-colon-list CPLUS_INCLUDE_PATH "$prefix/include:$prefix/include/ncurses"
add-to-colon-list PKG_CONFIG_PATH "$prefix/lib/pkgconfig"

if [ -z "$MANPATH" ] ; then
    export MANPATH="$prefix/share/man:$(man --path)"
else
    export MANPATH="$prefix/share/man:$MANPATH"
fi

if [ -z "$TEXINPUTS" ] ; then
    export TEXINPUTS=".:$prefix/share/texmf/tex/latex/local//:"
else
    export TEXINPUTS=".:$prefix/share/texmf/tex/latex/local//:${TEXINPUTS#.:}"
fi

wd=$(pwd)
distfiles=$wd/distfiles
temp_dir=/dev/shm/$(whoami)/temp/cps-build
src_dir=$temp_dir/src
build_dir=$temp_dir/build

rm -rf $temp_dir || true
