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
    f = f"database/trinucleotide-completeness-{motif}.txt"
    with open(f) as fp:
        curr_rmsd = []
        curr_closest = []
        for lnr, l in enumerate(fp):
            ll = l.split()
            assert ll[0] == motif
            assert len(ll) in (3, 4), lnr
            assert ll[1] == str(lnr), (lnr, l)
            curr_rmsd.append(float(ll[2]))
            cc = None
            if len(ll) == 4:
                cc = int(ll[3])
            curr_closest.append(cc)
        curr_rmsd = np.array(curr_rmsd)
        rmsd[motif] = curr_rmsd
        closest[motif] = curr_closest

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
    assert len(clusters) == len(closest[motif]), (
        motif,
        len(clusters),
        len(closest[motif]),
    )

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
        b7f[chain][ind] = (
            motif,
            clusnr,
            rmsd[motif][clusnr],
            closest[motif][clusnr],
        )

for chain in ("P", "Q"):
    b = b7f[chain]
    for k in sorted(b.keys()):
        f1, f2, f3, f4 = b[k]
        motif = f1
        other_clust = f4
        print(
            chain,
            k,
            f1,
            f2,
            f3,
            fragj[motif][other_clust],
        )
    print()
