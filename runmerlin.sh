#!/bin/sh
#$ -cwd
#$ -S /bin/bash

merlin -d seizures.dat -p seizures.ped -m seizures.map -f seizures.freq --npl --exp --tabulate --markerNames --pdf
mv merlin.pdf seizures.npl-markerNames.pdf
mv merlin-nonparametric.tbl seizures.npl-markerNames-table.txt
