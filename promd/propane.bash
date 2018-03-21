#!/bin/bash

make distclean && make -j12 omp promd && ./promd_omp -nt 4 -imd ../senseanalysis/examples/salome/propane/propane.imd -cnf ../senseanalysis/examples/salome/propane/propane.cnf -top ../senseanalysis/examples/salome/propane/propane.top -omd propane-run.omd -trc propane-run.g96 -trf propane-run.trf

./analyses.bash propane-run.omd
