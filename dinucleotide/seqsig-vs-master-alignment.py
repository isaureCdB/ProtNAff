# Master alignment (master.fasta) and leader seqs (rrm-leadseqs.fasta) from RRM GitHub project
seqsigs = []
with open("double-stacking-sorted-seqsig.txt") as fp:
    for l in fp.readlines():
        seqsigs.append(l.strip())

sequences = []
sequence_titles = []

with open("master.fasta") as fp:
    lines = fp.readlines()
    for n in range(0, len(lines), 2):
        title = lines[n].strip()
        seq = lines[n + 1].strip()
        sequences.append(seq)
        sequence_titles.append(title[1:].strip())

with open("rrm-leadseqs.fasta") as fp:
    lines = fp.readlines()
    for n in range(len(sequences)):
        title = lines[2 * n].strip()
        leadseq = lines[2 * n + 1].strip()
        assert len(leadseq) == 7, leadseq
        title = title[1:].strip()
        assert title == sequence_titles[n], (title, sequence_titles[n])
        sequences[n] = leadseq + sequences[n]


def compare_seqsig(seqsig1, seqsig2):
    identities = 0
    for c1, c2 in zip(seqsig1, seqsig2):
        if c1 == c2:
            identities += 1
        elif c1 == "-" or c2 == "-":
            identities += 0.5
    if identities - 3 >= 24:  # 80 %
        return True
    return False


unique_stack_seqsigs = {}
unique_nonstack_seqsigs = {}
nstack = 0
nnonstack = 0
for n, seq0 in enumerate(sequences):
    seq = seq0[7:]
    if not (seq[1] in "FYHW" and seq[298] in "FYHW"):
        print("Not two aromatics", sequence_titles[n])
        continue
    offset1 = seq[:2].count("-")
    offset2 = seq[:299].count("-")
    seqr = seq0.replace("-", "")
    leaderlen = len(seq0[:7].replace("-", ""))
    pos1, pos2 = 1 - offset1 + leaderlen, 298 - offset2 + leaderlen
    assert seqr[pos1] == seq[1], seq0
    assert seqr[pos2] == seq[298], seq0
    rrm_seqsig = seqr[max(pos1 - 7, 0) : pos1 + 8] + "..." + seqr[pos2 - 7 : pos2 + 8]
    if pos1 - 7 < 0:
        rrm_seqsig = (7 - pos1) * "-" + rrm_seqsig
    rrm_seqsig = rrm_seqsig.upper()
    # print(sequence_titles[n], rrm_seqsig)

    for seqsig in seqsigs:
        if compare_seqsig(rrm_seqsig, seqsig):
            print("Has (n - n+1) double stacking,", sequence_titles[n])
            unique_seqsigs = unique_stack_seqsigs
            nstack += 1
            break
    else:
        print("No (n - n+1) double stacking", sequence_titles[n])
        unique_seqsigs = unique_nonstack_seqsigs
        nnonstack += 1

    for unique_seqsig in unique_seqsigs:
        if compare_seqsig(rrm_seqsig, unique_seqsig):
            break
    else:
        unique_seqsigs[rrm_seqsig] = sequence_titles[n]

    # print(motif)
    # break

print()
f1 = "rrm-doublestack-seqsigs.list"
print(
    f"{nstack} double stackings, {len(unique_stack_seqsigs)} unique seqsigs, write to {f1}"
)
with open(f1, "w") as fp:
    for seqsig in unique_stack_seqsigs:
        print(unique_stack_seqsigs[seqsig], seqsig, file=fp)
f2 = "rrm-nonstack-seqsigs.list"
print(
    f"{nnonstack} double aromatics without double stackings, {len(unique_nonstack_seqsigs)} unique seqsigs, write to {f2}"
)
with open(f2, "w") as fp:
    for seqsig in unique_nonstack_seqsigs:
        print(unique_nonstack_seqsigs[seqsig], seqsig, file=fp)
