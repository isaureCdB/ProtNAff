motifs = []
for n1 in ("A", "C"):
    for n2 in ("A", "C"):
        for n3 in ("A", "C"):
            motifs.append(n1 + n2 + n3)

for motif in motifs:
    print(motif)
    data = {}
    for hexind in range(32):
        fname = f"database/trinucleotide-completeness-{motif}-{hexind+1}"
        with open(fname) as f:
            for l in f:
                ll = l.split()
                assert ll[0] == motif, (fname, l)
                pos = int(ll[1])
                if pos in data:
                    assert data[pos] == l, (data[pos], l)
                else:
                    data[pos] = l
    fname2 = f"database/trinucleotide-completeness-{motif}.txt"
    with open(fname2, "w") as f:
        print(len(data))
        for n in range(len(data)):
            l = data[n]
            print(l.rstrip(), file=f)
