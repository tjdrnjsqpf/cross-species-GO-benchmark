#!/usr/bin/env python3
"""gen_config.py — resolve UniProt reference-proteome IDs (API) for a species registry
and emit a track YAML (same schema as track_fish.yaml). GO source auto: GOA dir if the
species has one, else QuickGO by taxid. Logs species with no reference proteome (skipped).
Registries are defined inline below; run:  gen_config.py <trackname>
"""
import sys, time, urllib.request, urllib.parse

GOA_DIRS = {  # species with a dedicated UniProt-GOA species GAF (keeps MOD experimental codes)
 "human":"HUMAN/goa_human.gaf.gz","mouse":"MOUSE/goa_mouse.gaf.gz","rat":"RAT/goa_rat.gaf.gz",
 "cow":"COW/goa_cow.gaf.gz","dog":"DOG/goa_dog.gaf.gz","pig":"PIG/goa_pig.gaf.gz",
 "chicken":"CHICKEN/goa_chicken.gaf.gz","zebrafish":"ZEBRAFISH/goa_zebrafish.gaf.gz",
 "fruitfly":"FLY/goa_fly.gaf.gz","yeast":"YEAST/goa_yeast.gaf.gz","celegans":"WORM/goa_worm.gaf.gz",
 "arabidopsis":"ARABIDOPSIS/goa_arabidopsis.gaf.gz",
}

def resolve_upid(taxid):
    url=f"https://rest.uniprot.org/proteomes/search?query=organism_id:{taxid}+AND+reference:true&format=tsv&fields=upid,organism_id"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                lines=r.read().decode().splitlines()
            if len(lines)>=2:
                up,org=lines[1].split("\t")[:2]; return up, org
            return None,None
        except Exception:
            time.sleep(4)
    return None,None

# ---- registries: (species, taxid, My, role, richness, flags) ; focal first ----
REG = {
 "fungi": dict(focal="yeast", eqd=["pombe","aspergillus","neurospora"], species=[
    ("yeast",559292,0,"focal","high","wgd=post"),
    ("s_paradoxus",226125,5,"reference","mid","wgd=post"),
    ("s_uvarum",230603,20,"reference","low","wgd=post"),
    ("c_glabrata",284593,80,"reference","mid","wgd=post"),
    ("n_castellii",1064592,90,"reference","low","wgd=post"),
    ("z_rouxii",559307,100,"reference","low","wgd=pre"),
    ("k_lactis",284590,110,"reference","mid","wgd=pre"),
    ("torulaspora",4950,110,"reference","low","wgd=pre"),
    ("e_gossypii",284811,120,"reference","mid","wgd=pre"),
    ("k_marxianus",4911,110,"reference","low","wgd=pre"),
    ("c_albicans",237561,200,"reference","mid","wgd=none"),
    ("yarrowia",284591,300,"reference","mid","wgd=none"),
    ("aspergillus",227321,450,"reference","mid","wgd=none"),
    ("neurospora",367110,450,"reference","mid","wgd=none"),
    ("magnaporthe",242507,450,"reference","low","wgd=none"),
    ("pombe",284812,550,"reference","high","wgd=none"),
    ("cryptococcus",235443,800,"reference","mid","wgd=none"),
    ("ustilago",237631,800,"reference","low","wgd=none"),
    ("human",9606,1500,"reference","highest","wgd=none"),
 ]),
 "insect": dict(focal="fruitfly", eqd=["apis","tribolium","anopheles"], species=[
    ("fruitfly",7227,0,"focal","high","wgd=none"),
    ("d_simulans",7240,5,"reference","low","wgd=none"),
    ("d_yakuba",7245,10,"reference","low","wgd=none"),
    ("d_ananassae",7217,25,"reference","low","wgd=none"),
    ("d_pseudoobscura",7237,30,"reference","low","wgd=none"),
    ("d_virilis",7244,40,"reference","low","wgd=none"),
    ("d_grimshawi",7222,50,"reference","low","wgd=none"),
    ("anopheles",7165,250,"reference","mid","wgd=none"),
    ("aedes",7159,250,"reference","mid","wgd=none"),
    ("apis",7460,330,"reference","mid","wgd=none"),
    ("nasonia",7425,330,"reference","low","wgd=none"),
    ("tribolium",7070,350,"reference","mid","wgd=none"),
    ("bombyx",7091,350,"reference","mid","wgd=none"),
    ("aphid",7029,400,"reference","low","wgd=none"),
    ("daphnia",6669,500,"reference","low","wgd=none"),
    ("celegans",6239,700,"reference","mid","wgd=none"),
    ("mouse",10090,700,"reference","highest","wgd=none"),
    ("yeast",559292,1500,"reference","high","wgd=none"),
 ]),
}

# ---- high-%id plateau anchors (addendum §1/§4: REQUIRED for ID50; circularity-filtered in S3) ----
ANCHORS = {
 "fish": [  # cyprinids/characins close to zebrafish -> 73-88% plateau
    ("astyanax",7994,150,"mid","none"),("carp",7962,50,"mid","post"),       # carp WGD
    ("goldfish",7957,50,"mid","post"),("sinocyclocheilus",75366,40,"low","post"),
    ("grasscarp",79684,50,"low","none"),
 ],
 "plant_rice": [  # wild Oryza + sister genus -> 80-99% plateau
    ("o_rufipogon",4529,1,"low","none"),("o_nivara",4536,2,"low","none"),
    ("o_glaberrima",4538,3,"low","none"),("o_barthii",65489,3,"low","none"),
    ("o_punctata",4537,8,"low","none"),("o_brachyantha",4533,15,"low","none"),
    ("leersia",77586,20,"low","none"),
 ],
 "plant_arabidopsis": [  # Brassicaceae relatives -> 80-92% plateau
    ("a_lyrata",59689,10,"low","none"),("a_halleri",81970,8,"low","none"),
    ("capsella",81985,15,"mid","none"),("camelina",90675,20,"mid","none"),
    ("brassica_rapa",3711,25,"mid","none"),("brassica_oleracea",109376,25,"mid","none"),
 ],
}

# ---- densification: extra reference species to fill the %identity 30-65% gap ----
EXTRA = {
 "fish": [  # focal zebrafish
    ("stickleback",69293,230,"low","none"),("fugu",31033,230,"low","none"),
    ("tilapia",8128,230,"mid","none"),("salmon",8030,230,"mid","post"),  # salmonid WGD
    ("cod",8049,230,"low","none"),("gar",7918,430,"low","none"),         # slow-evolving, high %id
    ("coelacanth",7897,430,"low","none"),("elephantshark",7868,460,"low","none"),
    ("lamprey",7757,600,"low","none"),("amphioxus",7739,650,"low","none"),
    ("seaurchin",7668,700,"low","none"),("celegans",6239,700,"mid","none"),
 ],
 "mammal": [  # focal mouse
    ("rabbit",9986,90,"mid","none"),("horse",9796,95,"mid","none"),
    ("opossum",13616,160,"mid","none"),("platypus",9258,180,"low","none"),
    ("anole",28377,320,"low","none"),("gar",7918,435,"low","none"),
    ("lamprey",7757,600,"low","none"),("amphioxus",7739,650,"low","none"),
    ("seaurchin",7668,700,"low","none"),("celegans",6239,700,"mid","none"),
 ],
 "plant_rice": [  # focal rice
    ("barley",4513,50,"mid","none"),("grape",29760,160,"mid","none"),
    ("potato",4113,160,"mid","none"),("selaginella",88036,420,"low","none"),
    ("marchantia",3197,480,"low","none"),("klebsormidium",327967,800,"low","none"),
 ],
 "plant_arabidopsis": [  # focal arabidopsis
    ("grape",29760,110,"mid","none"),("potato",4113,110,"mid","none"),
    ("barley",4513,160,"mid","none"),("selaginella",88036,420,"low","none"),
    ("marchantia",3197,480,"low","none"),("klebsormidium",327967,800,"low","none"),
 ],
}

def extend(track, which="EXTRA"):
    import yaml
    path=f"/var2/lsg/Claude_Code/Cross-species-GeneOntology/config/track_{track}.yaml"
    cfg=yaml.safe_load(open(path)); added=0
    reg = (ANCHORS if which=="ANCHORS" else EXTRA).get(track, [])
    for name,taxid,My,rich,wgd in reg:
        if name in cfg["species"]:
            print(f"  have {name}"); continue
        up,org=resolve_upid(taxid)
        if not up: print(f"  SKIP {name} ({taxid}): no reference proteome"); continue
        gs="goa" if name in GOA_DIRS else "quickgo"
        e={"role":"reference","proteome":f"{up}_{org}","taxid":int(org),
           "go_source":gs,"timetree_My":My,"richness_class":rich,"wgd":wgd}
        if gs=="goa": e["goa_path"]=GOA_DIRS[name]
        cfg["species"][name]=e; added+=1
        print(f"  + {name:14s} {up}_{org} goa={gs} My={My}")
    yaml.safe_dump(cfg,open(path,"w"),sort_keys=False,default_flow_style=False)
    print(f"[extend] {track}: +{added} species -> {path}")

def emit(track):
    r=REG[track]; lines=[]
    lines.append(f"# Track {track} (auto-generated by gen_config.py). focal={r['focal']}.")
    lines.append(f"track: {track}")
    lines.append(f"focal: {r['focal']}")
    lines.append(f"equidistant_group_names: [{', '.join(r['eqd'])}]")
    lines.append("id_space: uniprot")
    lines.append("experimental_evidence: [EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP]")
    lines.append("transfer_methods: [besthit, rbh]")
    lines.append("diamond: {evalue: 1.0e-5, query_cover: 50, subject_cover: 50, max_target_seqs: 5, sensitivity: more-sensitive}")
    lines.append("species:")
    skipped=[]
    for name,taxid,My,role,rich,flags in r["species"]:
        up,org=resolve_upid(taxid)
        if not up:
            skipped.append((name,taxid)); print(f"  SKIP {name} ({taxid}): no reference proteome"); continue
        gs = "goa" if name in GOA_DIRS else "quickgo"
        goa = f", goa_path: {GOA_DIRS[name]}" if gs=="goa" else ""
        wgd = dict(kv.split("=") for kv in flags.split(",")).get("wgd","none")
        print(f"  {name:16s} {up}_{org}  goa={gs}  My={My} wgd={wgd}")
        lines.append(f"  {name}: {{role: {role}, proteome: {up}_{org}, taxid: {org}, go_source: {gs}{goa}, timetree_My: {My}, richness_class: {rich}, wgd: {wgd}}}")
    path=f"/var2/lsg/Claude_Code/Cross-species-GeneOntology/config/track_{track}.yaml"
    open(path,"w").write("\n".join(lines)+"\n")
    print(f"[gen] wrote {path} ({len(r['species'])-len(skipped)} species, {len(skipped)} skipped)")

if __name__=="__main__":
    if sys.argv[1]=="extend":
        extend(sys.argv[2], sys.argv[3] if len(sys.argv)>3 else "EXTRA")
    else:
        emit(sys.argv[1])
