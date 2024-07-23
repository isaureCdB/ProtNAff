import itertools
import os
import numpy as np

sequences = ["".join(seq) for seq in itertools.product(('A','C'), repeat=3)]

old_f1, old_f2 = None, None
d1, d2 = None, None
for seq in sequences:
    for mode in "forward", "backward":
        if mode == "forward":
            ###pairs = [f'{seq}{c}-compatibility-rmsd-full.npy' for c in ("A", "C")]
            pairs = [f'{seq}{c}-compatibility-rmsd.npy' for c in ("A", "C")] ###
        else:
            ###pairs = [f'{c}{seq}-compatibility-rmsd-full.npy' for c in ("A", "C")]
            pairs = [f'{c}{seq}-compatibility-rmsd.npy' for c in ("A", "C")] ###
        f1, f2 = pairs
        if not os.path.exists(f1) or not os.path.exists(f2):
            continue
        print(f"{mode} analysis for {seq[:2]} in {seq}")
        old_d1, old_d2 = d1, d2
        if f1 == old_f1:
            d1 = old_d1
        elif f1 == old_f2:
            d1 = old_d2
        else:
            d1 = np.load(f1)
        if f2 == old_f1:
            d2 = old_d1
        elif f2 == old_f2:
            d2 = old_d2
        else:
            d2 = np.load(f2)
        del old_d1, old_d2
        
        print(d1.shape, d2.shape, mode)
        if mode == "forward":
            ax1, ax2 = 1, 1
        else:
            ax1, ax2 = 0, 0
        bestind1 = d1.argmin(axis=ax1)
        bestind2 = d2.argmin(axis=ax2)
        
        if ax1 == 0:
            best1 = d1[bestind1]
            d1[bestind1] = 999
        else:
            best1 = d1[:, bestind1]
            d1[:, bestind1] = 999
        if ax2 == 0:
            best2 = d2[bestind2]
            ###d2[bestind2] = 999
        else:
            best2 = d2[:, bestind2]
            ###d2[:, bestind2] = 999

        m1 = d1.min(axis=ax1)
        m2 = d2.min(axis=ax2)
        m = np.minimum(m1, m2)
        
        if ax1 == 0:
            d1[bestind1] = best1
        else:
            d1[:, bestind1] = best1
        if ax2 == 0:
            d2[bestind2] = best2
        else:
            d2[:, bestind2] = best2

        print(np.histogram(m))
        print("Percentage above 0.5: {:.2f}".format((m>0.5).mean() * 100))
        print("Percentage above 0.75: {:.2f}".format((m>0.75).mean() * 100))
        print("Percentage above 1: {:.2f}".format((m>1).mean() * 100))
        print("Percentage above 1.25: {:.2f}".format((m>1.25).mean() * 100))
        print()