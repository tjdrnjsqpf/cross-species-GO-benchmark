#!/bin/bash
set -e
RLIB=${CONDA_PREFIX}/lib/R/library
cd /tmp
URL="https://bioconductor.org/packages/3.20/data/annotation/src/contrib/GenomeInfoDbData_1.2.13.tar.gz"
ALT="https://depot.galaxyproject.org/software/bioconductor-genomeinfodbdata/bioconductor-genomeinfodbdata_1.2.13_src_all.tar.gz"
echo "downloading..."
curl -fSL --retry 3 -o gid.tar.gz "$URL" || curl -fSL --retry 3 -o gid.tar.gz "$ALT"
ls -l gid.tar.gz; file gid.tar.gz
R CMD INSTALL -l "$RLIB" gid.tar.gz
echo "INSTALL_DONE"
