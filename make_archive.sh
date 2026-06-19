#!/usr/bin/env bash
# Build the Zenodo deposit archive from this folder.
set -e
cd "$(dirname "$0")"
rm -f cross-species-GO-benchmark_deposit.zip
zip -r -q cross-species-GO-benchmark_deposit.zip code data README.md LICENSE DATA_SOURCES.md RESULTS_SUMMARY.md
echo "built cross-species-GO-benchmark_deposit.zip"
