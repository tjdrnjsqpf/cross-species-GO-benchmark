#!/bin/bash
# install a placeholder bioconductor DATA package by key (e.g. go.db-3.20.0)
set -e
ENV=/var2/lsg/miniforge3/envs/gotx-r
KEY=$1
JSON=$ENV/share/bioconductor-data-packages/dataURLs.json
URL=$($ENV/bin/yq -r ".\"$KEY\".urls[0]" $JSON)
ALT=$($ENV/bin/yq -r ".\"$KEY\".urls[2] // .\"$KEY\".urls[1]" $JSON)
cd /tmp
echo "KEY=$KEY URL=$URL"
curl -fSL --retry 3 -o pkg.tar.gz "$URL" || curl -fSL --retry 3 -o pkg.tar.gz "$ALT"
ls -l pkg.tar.gz; file pkg.tar.gz
$ENV/bin/R CMD INSTALL -l $ENV/lib/R/library pkg.tar.gz
rm -f pkg.tar.gz
echo "DONE_$KEY"
