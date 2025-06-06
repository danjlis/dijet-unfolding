#!/usr/bin/env bash

i=$1
source ./roounfold/RooUnfold-master/build/setup.sh

root -l -q -b "makeTruth_dist.C(\"binning.config\", ${i})"

root -l -q -b "makeIterationPlot.C(${i})"

root -l -q -b "drawFullClosure.C(${i})"

root -l -q -b "drawHalfClosure.C(${i})"

root -l -q -b "drawSysJESJER.C(${i})"

root -l -q -b "drawFinalUnfold.C(${i})"
