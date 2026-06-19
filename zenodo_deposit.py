#!/usr/bin/env python3
"""Deposit the code+data archive to Zenodo and obtain a DOI (for the Data-availability statement).

Setup once:
  - Zenodo (or https://sandbox.zenodo.org for testing) > top-right menu > Applications >
    Personal access tokens > New token, scopes: deposit:write deposit:actions
  - export ZENODO_TOKEN=<the token>
  - EDIT the METADATA below (especially `creators` = all authors).

Usage:
  python zenodo_deposit.py cross-species-GO-benchmark_deposit.zip            # create DRAFT (review on web, then publish)
  python zenodo_deposit.py cross-species-GO-benchmark_deposit.zip --sandbox  # test on sandbox first
  python zenodo_deposit.py cross-species-GO-benchmark_deposit.zip --publish  # create AND publish -> DOI

Recommended: run on --sandbox once, confirm the record looks right, then run real (DRAFT), review on
the Zenodo web UI, and click Publish there. Put the resulting DOI in the manuscript Data-availability.
"""
import os, sys, json, requests

SANDBOX = "--sandbox" in sys.argv
PUBLISH = "--publish" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
ARCHIVE = args[0] if args else "cross-species-GO-benchmark_deposit.zip"
BASE = "https://sandbox.zenodo.org/api" if SANDBOX else "https://zenodo.org/api"
TOKEN = os.environ.get("ZENODO_TOKEN")
assert TOKEN, "Set ZENODO_TOKEN (Zenodo > Applications > Personal access tokens; scope deposit:write)."
assert os.path.exists(ARCHIVE), f"archive not found: {ARCHIVE} (run: bash make_archive.sh)"

# ---- EDIT before publishing ----
METADATA = {"metadata": {
    "title": "How far can you borrow? Divergence limits of cross-species Gene Ontology enrichment — code and processed data",
    "upload_type": "dataset",
    "description": ("Pipeline code and processed result tables for the cross-species GO-transfer "
        "divergence-tolerance benchmark. Raw data derive from public resources (UniProt, GO Annotation "
        "database, Bgee, EBI Expression Atlas, TimeTree); see DATA_SOURCES.md. Reproduce with code/."),
    "creators": [
        {"name": "[Last, First]", "affiliation": "[affiliation]"},
        # add every author, e.g. {"name": "Lee, Sung-Gwon", "affiliation": "..."},
    ],
    "license": "cc-by-4.0",
    "keywords": ["Gene Ontology", "enrichment analysis", "orthology", "non-model organisms",
                 "functional annotation", "cross-species", "benchmark"],
}}
# --------------------------------

P = {"access_token": TOKEN}
s = requests.Session(); s.headers.update({"Content-Type": "application/json"})

print(f"[1/3] creating deposition on {BASE} ...")
r = s.post(f"{BASE}/deposit/depositions", params=P, json={}); r.raise_for_status()
dep = r.json(); dep_id = dep["id"]; bucket = dep["links"]["bucket"]
print(f"      deposition id = {dep_id}")

print(f"[2/3] uploading {ARCHIVE} ({os.path.getsize(ARCHIVE)/1e6:.1f} MB) ...")
with open(ARCHIVE, "rb") as fp:
    r = requests.put(f"{bucket}/{os.path.basename(ARCHIVE)}", data=fp, params=P); r.raise_for_status()
print(f"      uploaded, checksum {r.json().get('checksum')}")

print("[3/3] setting metadata ...")
r = s.put(f"{BASE}/deposit/depositions/{dep_id}", params=P, data=json.dumps(METADATA)); r.raise_for_status()

if PUBLISH:
    r = s.post(f"{BASE}/deposit/depositions/{dep_id}/actions/publish", params=P); r.raise_for_status()
    o = r.json()
    print("PUBLISHED.  DOI:", o.get("doi"), "| record:", o["links"].get("record_html"))
else:
    print("DRAFT created (NOT published). Review/edit metadata and publish here:")
    print("   ", dep["links"]["html"])
