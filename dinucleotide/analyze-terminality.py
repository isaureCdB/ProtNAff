import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

frags = np.load("database/fragments_clust.npy")

data = {}
terminality = {}
fragdata = {}
for motif3 in ("AAA", "AAC", "ACA", "ACC", "CAA", "CAC", "CCA", "CCC"):
    fr = frags[frags["motif"] == motif3.encode()]
    fr2 = fr[fr["clust0.2A_center"]]
    fr_order = np.argsort(fr2["clust0.2A"])
    fr3 = fr2[fr_order]
    fragdata[motif3] = fr3
    for pos in 0, 1:
        f = f"database/dinucleotide-completeness-{motif3}-{pos}.txt"
        d = np.loadtxt(f, usecols=[2, 3])
        ind = d[:, 0]
        assert np.allclose(ind, np.arange(len(d))), f
        data[motif3, pos] = d[:, 1]
    assert len(data[motif3, 0]) == len(data[motif3, 1])
    nclust = len(data[motif3, 0])

    f = f"database/trilib/{motif3}-clust0.2A-terminality.npy"
    term0 = np.load(f)
    assert len(term0) == nclust, motif3

    # pos = 0
    term = np.zeros(
        (nclust, 3), bool
    )  # 5terminal (first 2), middle, 3terminal (last 2)
    term[:, 0] = term0[:, 0] | term0[:, 1]
    term[:, 1] = term0[:, 2] | term0[:, 3]
    term[:, 2] = term0[:, 4]
    terminality[motif3, 0] = term

    # pos = 1
    term = np.zeros(
        (nclust, 3), bool
    )  # 5terminal (first 2), middle, 3terminal (last 2)
    term[:, 0] = term0[:, 0]
    term[:, 1] = term0[:, 1] | term0[:, 2]
    term[:, 2] = term0[:, 3] | term0[:, 4]
    terminality[motif3, 1] = term

for motif in ("AA", "AC", "CA", "CC"):
    cdata = []
    cterm = []
    cfrag = []
    for pos in 0, 1:
        for motif3 in ("AAA", "AAC", "ACA", "ACC", "CAA", "CAC", "CCA", "CCC"):
            if motif3[pos : pos + 2] != motif:
                continue
            cdata.append(data[motif3, pos])
            cterm.append(terminality[motif3, pos])
            cfrag.append(fragdata[motif3])
    cdata = np.concatenate(cdata)
    cterm = np.concatenate(cterm).astype(bool)
    cfrag = np.concatenate(cfrag)

    assert len(cdata) == len(cterm) == len(cfrag)

    mask_all = np.ones(len(cdata), bool)
    mask_terminal = cterm[:, 0] | cterm[:, 2]
    mask_middle_exclusive = ~mask_terminal

    for maskname, mask in (
        ("all", mask_all),
        ("exclusive middle position", mask_middle_exclusive),
        ("terminal position", mask_terminal),
    ):
        nearest = cdata[mask]
        pct_mask = mask.mean() * 100
        if maskname == "all":
            pct_mask = len(mask)
        pct_05 = (nearest > 0.5).mean() * 100
        pct_1 = (nearest > 1).mean() * 100
        print(
            f"{motif}, {maskname} fragments ({pct_mask:.1f} %), {pct_05:.1f} % above 0.5 A, {pct_1:.1f} % above 1.0 A"
        )

        plt.cla()
        sns.histplot(nearest, bins=300)
        plt.yscale("log")
        plt.savefig(
            f"database/dinucleotide-completeness-{motif}-{maskname.replace(' ', '_')}.png"
        )
    bad = cdata > 1
    pct_bad = (bad & mask_terminal).sum() / bad.sum() * 100
    print(f"{motif}: terminal position fragments: {pct_bad:.1f} % of >1A fragments")
    print()
    npdb = len(np.unique(cfrag["structure"]))
    pct1 = int(npdb / 100 * 1)
    pct2 = int(npdb / 100 * 2)
    pct5 = int(npdb / 100 * 5)
    bad_per_pdb = np.unique(cfrag[bad]["structure"], return_counts=True)
    bad_per_pdb_sorted = np.sort(bad_per_pdb[1])[::-1]
    nbad = bad.sum()
    pct_bad1 = bad_per_pdb_sorted[:pct1].sum() / nbad.sum() * 100
    pct_bad2 = bad_per_pdb_sorted[:pct2].sum() / nbad.sum() * 100
    pct_bad5 = bad_per_pdb_sorted[:pct5].sum() / nbad.sum() * 100
    print(f"{motif}: Frags >1A in the worst 1% of PDB structures: {pct_bad1:.1f} %")
    print(f"{motif}: Frags >1A in the worst 2% of PDB structures: {pct_bad2:.1f} %")
    print(f"{motif}: Frags >1A in the worst 5% of PDB structures: {pct_bad5:.1f} %")
    print()
