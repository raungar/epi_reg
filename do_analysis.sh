#!/bin/bash
./scripts/make_chip_matrix.sh hotair chr12 53962308 53974956
./scripts/average_rna_celltypes.sh hotair
Rscript scripts/correlation.R hotair
