#!/bin/bash

rm -f ../results/*

runjob --block $1 --args -qmp-geom native --exe ../binaries/BGQ.x
