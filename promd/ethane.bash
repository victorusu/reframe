#!/bin/bash

make distclean && make -j12 omp promd && ./promd_omp -nt 4 -imd ../senseanalysis/examples/salome/ethane/ethane.imd -cnf ../senseanalysis/examples/salome/ethane/ethane.cnf -top ../senseanalysis/examples/salome/ethane/ethane.top -omd ethane-run.omd -trc ethane-run.g96 -trf ethane-run.trf

./analyses.bash ethane-run.omd
