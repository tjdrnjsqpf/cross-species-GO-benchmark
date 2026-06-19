#!/bin/bash
set -u
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"; cd $ROOT
GOTX=python3
S=logs/s17_run.log; : > $S
say(){ echo "[$(date +%H:%M:%S)] $*" >> $S; }
say "waiting for STRING download"
while ! grep -q "STRING DOWNLOAD DONE" logs/s17_string_dl.log 2>/dev/null; do sleep 30; done
say "download done ($(du -sh data/string|cut -f1)); running S17b"
$GOTX scripts/S17b_string_conservation.py >> $S 2>&1; say "S17b rc=$?"
$GOTX scripts/S17c_integrate.py >> $S 2>&1; say "S17c rc=$?"
say "S17 STRING ALL DONE"
