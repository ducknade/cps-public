#!/bin/bash


exec >log.out 2>&1

runjob --block $1 --label --args -qmp-geom native --exe ../binaries/BGQ.x
