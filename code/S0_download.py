#!/usr/bin/env python3
"""S0 (generic) — config-driven data acquisition. Downloads proteome FASTA + GAF
(GOA dir, or QuickGO by taxid) for every species in a track config. Skips existing.
Usage: S0_download.py config/track_mammal.yaml
"""
import os, sys, gzip, time, urllib.request, yaml

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
PBASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota"
GBASE = "https://ftp.ebi.ac.uk/pub/databases/GO/goa"
QGO = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
QLIMIT = 2000000

def fetch(url, out, tries=3):
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        return "skip"
    for i in range(tries):
        try:
            urllib.request.urlretrieve(url, out); return "ok"
        except Exception as e:
            last = e; time.sleep(5)
    print(f"   FAIL {url}: {last}"); return "fail"

def main():
    cfg = yaml.safe_load(open(sys.argv[1]))
    pdir = f"{ROOT}/data/proteomes"; gdir = f"{ROOT}/data/gaf"
    os.makedirs(pdir, exist_ok=True); os.makedirs(gdir, exist_ok=True)
    for sp, v in cfg["species"].items():
        up = v["proteome"]; updir = up.split("_")[0]; taxid = v["taxid"]
        # proteome
        fa = f"{pdir}/{sp}.fasta.gz"
        r = fetch(f"{PBASE}/{updir}/{up}.fasta.gz", fa)
        print(f"[{sp}] proteome {r} ({os.path.getsize(fa) if os.path.exists(fa) else 0} B)")
        # GAF
        gaf = f"{gdir}/{sp}.gaf.gz"
        if os.path.exists(gaf) and os.path.getsize(gaf) > 1000:
            print(f"[{sp}] gaf skip"); continue
        if v["go_source"] == "goa":
            r = fetch(f"{GBASE}/{v['goa_path']}", gaf)
            print(f"[{sp}] gaf(goa) {r}")
        else:  # quickgo by taxid -> gzip
            tmp = f"{gdir}/{sp}.gaf"
            url = f"{QGO}?taxonId={taxid}&taxonUsage=exact&downloadLimit={QLIMIT}"
            try:
                req = urllib.request.Request(url, headers={"Accept": "text/gaf"})
                with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fo:
                    fo.write(resp.read())
                with open(tmp, "rb") as fi, gzip.open(gaf, "wb") as fo:
                    fo.writelines(fi)
                os.remove(tmp)
                n = sum(1 for l in gzip.open(gaf, "rt") if not l.startswith("!"))
                flag = "  <<< HIT LIMIT, may be truncated" if n >= QLIMIT else ""
                print(f"[{sp}] gaf(quickgo) ok  {n} lines{flag}")
            except Exception as e:
                print(f"[{sp}] gaf(quickgo) FAIL: {e}")
    print("[S0] download done for", cfg["track"])

if __name__ == "__main__":
    main()
