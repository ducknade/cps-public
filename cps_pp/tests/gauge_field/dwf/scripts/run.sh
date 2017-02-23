#!/bin/bash

export BG_ALLOW_CACHELINE_LOCKING=1
export BG_L2LOCK_L1_ONLY=1
export MUSPI_NUMINJFIFOS=64
export MUSPI_NUMRECFIFOS=64
export MUSPI_NUMBATIDS=256
export MUSPI_NUMCLASSROUTES=1

runjob --block $1 --env-all --args -qmp-geom native --exe ../binaries/BGQ.x
