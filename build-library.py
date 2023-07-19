""" See build-library.md for documentation

Also print out fragment stats (PDB code, PDB chain, resid) + rmsd for every red and quasi-unique singleton
"""
import sys
import numpy as np
import json
from tqdm import tqdm
from nefertiti.functions.superimpose import superimpose, superimpose_array

from collections import namedtuple

Cluster = namedtuple("Cluster", ("pdb", "fragment_index", "rmsd", "rmsd2", "size", "npdbs"))

conformer_list = [] # all fragment indices that are "conformers" (primary and quaternary list, 
                    # quasi-unique list,
                    # and those tertiary that come from singletons)
primary_list = []   # PDB code + fragment index
secondary_list = [] # (primary) PDB code + (replacement) fragment index
tertiary_list = []  # Cluster instances. "pdb" refers to primary PDB code. "size"/"npdbs" refer to replacement. 
                    # "rmsd" is the RMSD between replacement and primary, "rmsd2" is the RMSD if the replacement is not used
quaternary_list = [] # (primary) PDB code + (primary) fragment index + (replacement) fragment index + rmsd
eliminated_list = []  # Cluster instances. "pdb" refers to primary PDB code. "size"/"npdbs" refer to replacement
quasi_unique_list = []  # PDB code + primary fragment index + quasi-replacement fragment index + RMSD
intra_cluster_list = [] # non-"conformers", coming from the secondary list and the tertiary list from non-singletons

frag_json = sys.argv[1] # fragments.json
motif = sys.argv[2] # e.g. AAA
cluster_threshold = float(sys.argv[3]) # e.g. 0.5
greenred_threshold = float(sys.argv[4]) # e.g. 1.0
clust02 = sys.argv[5] # e.g trilib/AAA-aa-fit-clust0.2
clust02npy = sys.argv[6] # e.g trilib/AAA-aa-fit-clust0.2.npy
final_clust = sys.argv[7] # e.g trilib/AAA-recluster-0.5-all.list
closest_cluster_rmsd = sys.argv[8] # e.g trilib/AAA-recluster-0.5-closest-cluster.rmsd
outfile_pattern = sys.argv[9] # e.g trilib/AAA-lib . 
                               # Files created: 
                               # X.list, X-secondary.list, X-tertiary.list, 
                               # X-quaternary.list, X-quaternary.stats,
                               # X-quasi-unique.list, X-quasi-unique.stats,
                               # X-conformer.list, X-conformer.npy, 
                               # X-intracluster.list, X-intracluster.npy


def read_clusterfile(clusterfile):
    clusters = []
    for lnr, l in enumerate(open(clusterfile)):
        ll = l.split()
        assert ll[0].lower() == "cluster"
        assert ll[1] == str(lnr+1)
        assert ll[2] == "->" or ll[2] == "=>"
        c = [int(lll) for lll in ll[3:]]
        clusters.append(c)
    return clusters

clusters_02 = read_clusterfile(clust02)
cluster_coordinates = np.load(clust02npy)
clusters_final = read_clusterfile(final_clust)

with open(frag_json) as f:
    frag = json.load(f)[motif]

sets_02 = []
for cluster_02 in clusters_02:
    set_02 = set()
    for fr in cluster_02:
        d = frag[str(fr)]
        struc = d["structure"][:4]
        set_02.add(struc)
    sets_02.append(set_02)


closest = []

membership = {}
with open(closest_cluster_rmsd) as f:
    for cluster_nr, l in enumerate(f):
        ll = l.split()
        cl = int(ll[1]) - 1
        rmsd = float(ll[0])
        closest.append((cl, rmsd))
        
        if rmsd >= cluster_threshold:
            continue
        clf = clusters_final[cl][0]
        if clf not in membership:
            membership[clf] = set()
        membership[clf].add(cluster_nr)
        clus_h = clusters_final[cluster_nr][0]
        if clus_h not in membership:
            membership[clus_h] = set()
        membership[clus_h].add(cl)
        rmsd = float(ll[0])
        if rmsd >= cluster_threshold:
            continue

assert cluster_nr == len(clusters_final) - 1
for cluster_nr, cluster in enumerate(clusters_final):
    for fr in cluster:          
        if fr not in membership:
            membership[fr] = set()
        membership[fr].add(cluster_nr)

# Sort clusters from smallest to largest
clus = np.argsort([len(c) for c in clusters_final])

is_singleton = {} # bool, PDB code
for cluster_nr in clus:
    cluster = clusters_final[cluster_nr]
    s = set()
    for fr in cluster:
        s.update(sets_02[fr-1])
    is_singleton[cluster_nr] = (len(s) == 1), (list(s)[0] if len(s) == 1 else None)


# First, remove all clusters where *all* members are part of another cluster
clus2 = []
for cluster_nr in clus:
    cluster = clusters_final[cluster_nr]
    for fr in cluster:
        if len(membership[fr]) == 1:
            break
    else:
        print(f"Eliminate fully redundant cluster {cluster_nr+1}")
        for fr in cluster:
            membership[fr].remove(cluster_nr)
        continue
    clus2.append(cluster_nr)
clus = clus2

cluster_struc = cluster_coordinates[[c[0]-1 for c in clusters_final]]
intra_cluster = set()

_rmsd = {}
def get_rmsd(cluster_nr):
    rmsd = _rmsd.get(cluster_nr)
    if rmsd is None:
        _, rmsd = superimpose_array(cluster_struc, cluster_struc[cluster_nr])
        _rmsd[cluster_nr] = rmsd
    return rmsd

# Determine quasi-unique clusters

for cluster_nr in tqdm(clus):
    if not is_singleton[cluster_nr][0]:
        continue
    pdb_code = is_singleton[cluster_nr][1]
    cl, r0 = closest[cluster_nr]
    if not is_singleton[cl][0]:
        continue
    pdb_code2 = is_singleton[cl][1]
    if pdb_code2 == pdb_code:   
        rmsd = get_rmsd(cluster_nr)
        for c_ind, close_cluster in enumerate(np.argsort(rmsd)):
            if c_ind == 0:
                assert rmsd[close_cluster] < 0.01 and close_cluster == cluster_nr
            elif c_ind == 1:
                assert abs(rmsd[close_cluster] - r0) < 0.01 and close_cluster == cl
            is_sing, pdb_code2 = is_singleton[close_cluster]
            if is_sing and pdb_code2 == pdb_code:
                continue
            r = rmsd[close_cluster]
            break
        else:
            raise Exception
        if (r0 <= greenred_threshold) and (r > greenred_threshold):
            quasi_unique_list.append((pdb_code, cluster_nr, cl, r0))
            #print(cluster_nr+1, cl+1, pdb_code, pdb_code2, r0, r)
        closest[cluster_nr] = (close_cluster, r)

print(f"{len(quasi_unique_list)} quasi-unique clusters")

# Red non-singletons

for cluster_nr in clus:
    if is_singleton[cluster_nr][0]:
        continue
    if closest[cluster_nr][1] <= greenred_threshold:
        continue
    heart = clusters_final[cluster_nr][0]
    pdb_codes = sets_02[heart-1]
    if len(pdb_codes) == 1:
        pdb_code = list(pdb_codes)[0]
    else:
        pdb_code = None
    primary_list.append((pdb_code, heart))
    if pdb_code is not None:
        for fr in clusters_final[cluster_nr]:
            pdb_codes = sets_02[fr-1]
            if len(pdb_codes) > 1 or list(pdb_codes)[0] != pdb_code:
                secondary_list.append((pdb_code, fr)) #primary pdb code!
                intra_cluster.add(fr)                
                break
        else:
            raise Exception # can't happen for non-singleton
            
print(f"{len(primary_list)} red non-singletons added to the primary list")
print(f"{len(secondary_list)} fragments in the secondary list")
checkpoint = len(primary_list)

# Green non-singletons
for cluster_nr in clus:
    if is_singleton[cluster_nr][0]:
        continue
    if closest[cluster_nr][1] > greenred_threshold:
        continue
    cluster = clusters_final[cluster_nr]
    heart = cluster[0]
    pdb_codes = sets_02[heart-1]
    if len(pdb_codes) == 1:
        pdb_code = list(pdb_codes)[0]
    else:
        pdb_code = None
    primary_list.append((pdb_code, heart))
    if pdb_code is not None:
        s = set()
        for fr in cluster:
            s.update(sets_02[fr-1])
        for fr in cluster:
            pdb_codes = sets_02[fr-1]
            if len(pdb_codes) > 1 or list(pdb_codes)[0] != pdb_code:
                _, r = superimpose(cluster_coordinates[heart-1], cluster_coordinates[fr-1])
                clust = Cluster(pdb=pdb_code, fragment_index=fr, size=len(cluster), npdbs=len(s), rmsd=r, rmsd2=closest[cluster_nr][1])
                intra_cluster.add(fr)
                tertiary_list.append(clust)
                break
        else:
            raise Exception # can't happen for non-singleton


print(f"{len(primary_list)-checkpoint} green non-singletons added to the primary list")
print(f"{len(tertiary_list)} fragments in the tertiary list")
checkpoint = len(primary_list)
checkpoint3 = len(tertiary_list)

# Type I green singletons

for cluster_nr in clus:
    if not is_singleton[cluster_nr][0]:
        continue
    nb, r  = closest[cluster_nr]
    if r > greenred_threshold:
        continue
    if is_singleton[nb][0]:
        continue
    pdb_code = is_singleton[cluster_nr][1]
    heart = clusters_final[cluster_nr][0]
    nb_clust = clusters_final[nb]
    s = set()
    for fr in nb_clust:
        s.update(sets_02[fr-1])
    _, rmsd_repl = superimpose(cluster_coordinates[heart-1], cluster_coordinates[nb_clust[0]-1])
    clust = Cluster(pdb=pdb_code, fragment_index=heart, rmsd=r, size=len(nb_clust), npdbs=len(s), rmsd2=None)
    eliminated_list.append(clust)

print(f"{len(eliminated_list)} green singletons eliminated")    
checkpoint_elim = len(eliminated_list)

# Type II green singletons

for cluster_nr in tqdm(clus):
    if not is_singleton[cluster_nr][0]:
        continue
    nb, r  = closest[cluster_nr] 
    if r > greenred_threshold:
        continue
    if not is_singleton[nb][0]:
        continue
    pdb_code = is_singleton[cluster_nr][1]
    heart = clusters_final[cluster_nr][0]
    primary_list.append((pdb_code, heart))
    
    nb_clust = clusters_final[nb]
    s = set()
    for fr in nb_clust:
        s.update(sets_02[fr-1])
    clust = Cluster(pdb=pdb_code, fragment_index=nb_clust[0], rmsd=r, size=len(nb_clust), npdbs=1, rmsd2=None)
    if nb_clust[0] in [c[1] for c in primary_list]:
        eliminated_list.append(clust)
    else:
        rmsds = get_rmsd(cluster_nr)
        rmsds[cluster_nr] = 99999
        rmsds[nb] = 99999
        for c_ind, c_nr in enumerate(np.argsort(rmsds)):
            is_sing, pdb_code2 = is_singleton[c_nr]
            if is_sing and pdb_code2 == pdb_code:
                continue
            rmsd2 = rmsds[c_nr]
            break
        clust = clust._replace(rmsd2 = rmsd2)
        tertiary_list.append(clust)

print(f"{len(primary_list)-checkpoint} green singletons added to the primary list")
if len(eliminated_list) > checkpoint_elim:
    print(f"{len(eliminated_list)-checkpoint_elim} fragments eliminated (replacement already in primary list)")
print(f"{len(tertiary_list)-checkpoint3} fragments added to the tertiary list")

checkpoint = len(primary_list)
checkpoint3 = len(tertiary_list)

# Red singletons

for cluster_nr in clus:
    if not is_singleton[cluster_nr][0]:
        continue
    nb, r  = closest[cluster_nr] 
    if r <= greenred_threshold:
        continue

    pdb_code = is_singleton[cluster_nr][1]    
    nb_clust = clusters_final[nb]
    item = (pdb_code, clusters_final[cluster_nr][0], nb_clust[0], r)
    quaternary_list.append(item)
print(f"{len(quaternary_list)} fragments added to the quaternary list")

conformer_list = set([c[1] for c in primary_list])
for pdb_code, fragment_index in secondary_list:
    intra_cluster_list.append(fragment_index)
for clust in tertiary_list:
    fragment_index = clust.fragment_index
    if fragment_index not in intra_cluster:
        conformer_list.add(fragment_index)
    else:
        intra_cluster_list.append(fragment_index)
for pdb_code, fragment_index0, fragment_index, r in quaternary_list:
    conformer_list.add(fragment_index)
for pdb_code, fragment_index0, fragment_index, r in quasi_unique_list:
    conformer_list.add(fragment_index)

conformer_list = sorted(list(set(conformer_list)))
intra_cluster_list = sorted(list(set(intra_cluster_list)))

print(f"Final list sizes: {len(primary_list)} (primary), {len(secondary_list)} (secondary), {len(tertiary_list)} (tertiary), {len(quaternary_list)} (quaternary), {len(quasi_unique_list)} (quasi-unique), {len(conformer_list)} (conformers), {len(intra_cluster_list)} (intra-cluster)")

# Write out lists, stats, and coordinate arrays

with open(outfile_pattern + ".list", "w") as f:
    for pdb_code, fragment_index in primary_list:
        if pdb_code is None:
            print(fragment_index, file=f)
        else:
            print(fragment_index, pdb_code, file=f)

with open(outfile_pattern + "-secondary.list", "w") as f:
    for pdb_code, fragment_index in secondary_list:
        print(pdb_code, fragment_index, file=f)

with open(outfile_pattern + "-tertiary.list", "w") as f:
    for cluster in sorted(tertiary_list, key= lambda c: c.pdb):
        print(cluster.pdb, cluster.fragment_index, "{:.3f}".format(cluster.rmsd), "{:.3f}".format(cluster.rmsd2), cluster.size, cluster.npdbs, file=f)

with open(outfile_pattern + "-eliminated.list", "w") as f:
    for cluster in sorted(eliminated_list, key= lambda c: -c.rmsd):
        print(cluster.pdb, cluster.fragment_index, "{:.3f}".format(cluster.rmsd), cluster.size, cluster.npdbs, file=f)

with open(outfile_pattern + "-quaternary.list", "w") as f:
    for pdb_code, fragment_index1, fragment_index2, r in sorted(quaternary_list, key = lambda c: -c[3]):
        print(pdb_code, fragment_index2, "{:.3f}".format(r), file=f)

with open(outfile_pattern + "-quaternary.stats", "w") as f:
    print("#PDB chain model resid rmsd", file=f)
    for pdb_code, fragment_index1, fragment_index2, r in sorted(quaternary_list, key = lambda c: -c[3]):
        for fr in clusters_02[fragment_index1-1]:
            curr_frag = frag[str(fr)]
            assert pdb_code == curr_frag["structure"][:4]            
            print(pdb_code, curr_frag["chain"], curr_frag["model"], f'{curr_frag["resid"][0]}-{curr_frag["resid"][2]}', "{:.3f}".format(r), file=f)

with open(outfile_pattern + "-quasi-unique.list", "w") as f:
    for pdb_code, cluster_index1, cluster_index2, r in sorted(quasi_unique_list, key = lambda c: c[0]):
        fragment_index1 = clusters_final[cluster_index1][0]
        fragment_index2 = clusters_final[cluster_index2][0]
        print(pdb_code, fragment_index1, fragment_index2, file=f)

with open(outfile_pattern + "-quasi-unique.stats", "w") as f:
    print("#PDB chain model resid chain2 model2 resid2 rmsd", file=f)
    for pdb_code, cluster_index1, cluster_index2, r in sorted(quasi_unique_list, key = lambda c: c[0]):
        fragment_index1 = clusters_final[cluster_index1][0]
        fragment_index2 = clusters_final[cluster_index2][0]
        for fr in clusters_02[fragment_index1-1]:
            curr_frag = frag[str(fr)]
            assert pdb_code == curr_frag["structure"][:4]
            for fr2 in clusters_02[fragment_index2-1]:
                curr_frag2 = frag[str(fr2)]
                assert pdb_code == curr_frag2["structure"][:4]

                print(pdb_code, 
                      curr_frag["chain"], curr_frag["model"], f'{curr_frag["resid"][0]}-{curr_frag["resid"][2]}',
                      curr_frag2["chain"], curr_frag2["model"], f'{curr_frag2["resid"][0]}-{curr_frag2["resid"][2]}'
                      , "{:.3f}".format(r), file=f)

with open(outfile_pattern + "-conformer.list", "w") as f:
    for fragment_index in sorted(conformer_list):
        print(fragment_index, file=f)
np.save(outfile_pattern + "-conformer.npy", cluster_coordinates[np.array(sorted(conformer_list))-1])

with open(outfile_pattern + "-intracluster.list", "w") as f:
    for fragment_index in sorted(intra_cluster_list):
        print(fragment_index, file=f)
np.save(outfile_pattern + "-intracluster.npy", cluster_coordinates[np.array(sorted(intra_cluster_list))-1])