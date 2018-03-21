#!/bin/bash

make distclean && make -j12 omp promd && ./promd_omp -nt 4 -imd ../senseanalysis/examples/salome/methane/methane.imd -cnf ../senseanalysis/examples/salome/methane/methane.cnf -top ../senseanalysis/examples/salome/methane/methane.top -omd methane-run.omd -trc methane-run.g96 -trf methane-run.trf

#make distclean && make -j12 omp promd && ./promd  -imd ../senseanalysis/examples/salome/methane/methane.imd -cnf ../senseanalysis/examples/salome/methane/methane.cnf -top ../senseanalysis/examples/salome/methane/methane.top -omd methane-run.omd -trc methane-run.g96 -trf methane-run.trf

./analyses.bash methane-run.omd
