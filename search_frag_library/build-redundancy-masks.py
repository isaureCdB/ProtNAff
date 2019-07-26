import json
import os
import numpy as np

seqids = 30, 40, 50, 70, 90, 95, 100
seqids_str = [str(seqid) for seqid in seqids]

# For now, we consider every nucleotide to be in contact with every protein chain
#  So we start to make a mapping from PDB code to seqclust IDs / PFAM IDs / UniProt IDs

pdb_to_seqclust = { seqid:{} for seqid in seqids}
pdb_to_pfam = {}
pdb_to_uniprot = {}

structures = json.load(open("structures.json"))
for seqid, seqid_str in zip(seqids, seqids_str):
    curr_pdb_to_seqclust = pdb_to_seqclust[seqid]
    for pdb, struc in structures.items():
        seq_clusts = set([clusts.get(seqid_str) for clusts in struc["seqclust"].values()])
        seq_clusts.discard(None)
        curr_pdb_to_seqclust[pdb.encode()] = sorted(seq_clusts)
for pdb, struc in structures.items():
    pfams = sorted(struc["pfam"].values())
    uniprots = sorted(struc["uniprot"].values())
    pdb_to_pfam[pdb.encode()] = pfams
    pdb_to_uniprot[pdb.encode()] = uniprots

try:
    os.mkdir("redundancy-masks")
except:
    pass
fragments = np.load("fragments_clust.npy")

def build_mask(clust, clust_center, pdb_to_ids):
    used_ids = {}
    mask = np.zeros(len(fragments))
    for fragnr, fragment in enumerate(fragments):
        if not clust_center[fragnr]:
            continue        
        clustid = fragment["motif"], clust[fragnr]
        assert clustid not in used_ids
        pdb = fragment["structure"]
        used_ids[clustid] = set(pdb_to_ids.get(pdb,[]))    
        mask[fragnr] = True

    for fragnr, fragment in enumerate(fragments):
        if clust_center[fragnr]:
            continue        
        clustid = fragment["motif"], clust[fragnr]
        if clustid not in used_ids: ### BUG in pipeline; center is not calculated for cluster 0...
            used_ids[clustid] = set()            
        curr_used_ids = used_ids[clustid]
        pdb = fragment["structure"]        
        new_used_ids = set(pdb_to_ids.get(pdb,[]))
        if curr_used_ids.intersection(new_used_ids):
            continue
        mask[fragnr] = True
        curr_used_ids.update(new_used_ids)
    return mask

print("%d fragments" % len(fragments))
print()
for clust_cutoff in "2.0", "1.0", "0.2":
    clust = fragments["clust" + clust_cutoff]
    clust_center = fragments["clust" + clust_cutoff + "_center"]
    print("Fragment clustering at %s A: %d clusters" % (clust_cutoff,clust_center.sum()))
    
    mask = build_mask(clust, clust_center, pdb_to_pfam)
    filename = "redundancy-masks/%s-pfam.npy" % clust_cutoff
    print("%s: every PFAM ID only once per fragment cluster, %d fragments" % (filename, mask.sum()))
    np.save(filename, mask)

    mask = build_mask(clust, clust_center, pdb_to_uniprot)
    filename = "redundancy-masks/%s-uniprot.npy" % clust_cutoff
    print("%s: every Uniprot ID only once per fragment cluster, %d fragments" % (filename, mask.sum()))
    np.save(filename, mask)    

    for seqid, seqid_str in zip(seqids, seqids_str):
        mask = build_mask(clust, clust_center, pdb_to_seqclust[seqid])
        filename = "redundancy-masks/%s-seqclust-%s.npy" % (clust_cutoff, seqid)
        print("%s: every %s %% seqid cluster only once per fragment cluster, %d fragments" % (filename, seqid, mask.sum()))
        np.save(filename, mask)    

    print()