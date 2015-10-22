#!/bin/bash

./cern-ldb -vs "ATLAS:LUMI_REGION_SIZE_X" -fn 1400 -bm1 STABLE -bm2 RAMPDOWN   -C ldb.conf -N Fill1400/Fill1400Atlasx -F CSV
