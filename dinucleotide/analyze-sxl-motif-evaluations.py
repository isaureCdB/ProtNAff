import numpy as np


def readf(f):
    seqsigs = []
    scores = []
    with open(f) as fp:
        for l in fp.readlines():
            seqsig, score = l.split()
            seqsigs.append(seqsig)
            scores.append(float(score))
    scores = np.array(scores)
    seqsigs = np.array(seqsigs)
    return seqsigs, scores


tp_seqsigs, tp_scores = readf("sxl-prot-motif-evaluate-predictions-truepos.txt")
fp_seqsigs, fp_scores = readf("sxl-prot-motif-evaluate-predictions-falsepos.txt")
rrm_doublestack_seqsigs = set(
    [l.strip().split()[1] for l in open("rrm-doublestack-seqsigs.list")]
)
rrm_nonstack_seqsigs = set(
    [l.strip().split()[1] for l in open("rrm-nonstack-seqsigs.list")]
)


def compare_seqsig(seqsig1, seqsig2):
    identities = 0
    for c1, c2 in zip(seqsig1, seqsig2):
        if c1 == c2:
            identities += 1
        elif c1 == "-" or c2 == "-":
            identities += 0.5
    return (identities - 3) / 30.0


def in_seqsigs(query_seqsig, seqsigs):
    max_score = 0
    max_seqsig = None
    for seqsig in seqsigs:
        score = compare_seqsig(query_seqsig, seqsig)
        if score >= 0.8 and score > max_score:
            max_score = score
            max_seqsig = seqsig
    return max_seqsig


print("True positives against RRM stacking")

n = 0
unique_tp_seqsigs = set()
for seqsig, score in zip(tp_seqsigs, tp_scores):
    for useqsig in unique_tp_seqsigs:
        if compare_seqsig(useqsig, seqsig) >= 0.8:
            break
    else:
        unique_tp_seqsigs.add(seqsig)
        closest_rrm_doublestack_seqsig = in_seqsigs(seqsig, rrm_doublestack_seqsigs)
        n += 1
        print(n, seqsig, "%.3f" % score, closest_rrm_doublestack_seqsig)
print()

print("Unique false positives <= 0.5A")
unique_fp_seqsigs = set()
for seqsig, score in zip(fp_seqsigs, fp_scores):
    if score >= 0.5:
        continue
    for useqsig in unique_fp_seqsigs:
        if compare_seqsig(useqsig, seqsig) >= 0.8:
            break
    else:
        unique_fp_seqsigs.add(seqsig)
        print(len(unique_fp_seqsigs), seqsig, "%.3f" % score)
print()

print("False positives against RRM non-stacking")
nonstack_seqsigs = rrm_nonstack_seqsigs.copy()
for seqsig, score in zip(fp_seqsigs, fp_scores):
    closest_rrm_nonstack_seqsig = in_seqsigs(seqsig, nonstack_seqsigs)
    if closest_rrm_nonstack_seqsig:
        print(seqsig, "%.3f" % score, closest_rrm_nonstack_seqsig)
        nonstack_seqsigs.remove(closest_rrm_nonstack_seqsig)
