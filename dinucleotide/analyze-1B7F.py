import numpy as np
import json


def read_clusterfile(clusterfile):
    clusters = []
    for lnr, l in enumerate(open(clusterfile)):
        ll = l.split()
        assert ll[0].lower() == "cluster"
        assert ll[1] == str(lnr + 1)
        assert ll[2] == "->" or ll[2] == "=>"
        c = [int(lll) for lll in ll[3:]]
        clusters.append(c)
    return clusters


motifs = []
for n1 in ("A", "C"):
    for n2 in ("A", "C"):
        for n3 in ("A", "C"):
            motifs.append(n1 + n2 + n3)

rmsd = {}
closest = {}
for motif in motifs:
    for pos in 0, 1:
        spos = str(pos)
        f = f"database/dinucleotide-completeness-{motif}-{pos}.txt"
        with open(f) as fp:
            curr_rmsd = []
            curr_closest = []
            for lnr, l in enumerate(fp):
                ll = l.split()
                assert ll[0] == motif
                assert ll[1] == spos
                assert ll[2] == str(lnr), (lnr, l)
                curr_rmsd.append(float(ll[3]))
                curr_closest.append(ll[4:7])
            curr_rmsd = np.array(curr_rmsd)
            rmsd[motif, pos] = curr_rmsd
            closest[motif, pos] = curr_closest

with open("database/fragments.json") as f:
    frags = json.load(f)

b7f = {"P": {}, "Q": {}}
fragj = {}
for motif in motifs:
    ffrags = frags[motif]
    ffragj = []
    fragj[motif] = ffragj

    f = f"database/trilib/{motif}-aa-fit-clust0.2"
    clusters = read_clusterfile(f)

    for clusnr, clus in enumerate(clusters):
        d0 = ffrags[str(clus[0])]
        ffragj.append((d0["structure"][:4], d0["chain"], d0["indices"][0]))
        for fr in clus:
            d = ffrags[str(fr)]
            struc = d["structure"][:4]
            if struc == "1B7F":
                break
        else:
            continue
        ind = d["indices"][0]
        chain = d["chain"]
        if ind == 1:
            b7f[chain][ind] = (
                motif[:2],
                clusnr,
                rmsd[motif, 0][clusnr],
                closest[motif, 0][clusnr],
            )
        b7f[chain][ind + 1] = (
            motif[1:],
            clusnr,
            rmsd[motif, 1][clusnr],
            closest[motif, 1][clusnr],
        )

for chain in ("P", "Q"):
    b = b7f[chain]
    for k in sorted(b.keys()):
        f1, f2, f3, f4 = b[k]
        other_motif = f4[0]
        other_pos = int(f4[1])
        other_clust = int(f4[2])
        print(
            chain,
            k,
            f1,
            f2,
            f3,
            other_motif,
            other_pos,
            fragj[other_motif][other_clust],
        )
    print()
