import json
import glob


def find_double_stacking(code, f):
    with open(f) as fp:
        d = json.load(fp)
    nuc = {}
    for stacking in d["stacks"]:
        nts0 = stacking["nts_long"].split(",")
        all_nuc = True
        nts = []
        for nt in nts0:
            l = nt.split(".")
            chain, resname = l[2:4]
            resid = int(l[4])
            if resname not in ("RU", "RG", "RC", "RA"):
                all_nuc = False
            nts.append((chain, resname, resid))
        if all_nuc:
            continue
        for chain, resname, resid in nts:
            if resname in ("RU", "RG", "RC", "RA"):
                mol = nuc
            else:
                continue
            if chain not in mol:
                mol[chain] = {}
            if resid not in mol[chain]:
                mol[chain][resid] = [resname]
        for n in range(len(nts)):
            chain1, resname1, resid1 = nts[n]
            if resname1 not in ("RU", "RG", "RC", "RA"):
                continue
            c1 = nuc[chain1][resid1]
            for nn in range(len(nts)):
                if n == nn:
                    continue
                chain2, resname2, resid2 = nts[nn]
                if resname2 in ("RU", "RG", "RC", "RA"):
                    continue
                if resname2 not in ("TYR", "PHE", "TRP", "ARG", "HIS"):
                    continue
                c1.append((chain2, resname2, resid2))

    for chain, stacks in nuc.items():
        for n in sorted(stacks):
            n2 = n + 1
            if n2 not in stacks:
                continue
            resname1 = stacks[n][0]
            for res1 in stacks[n][1:]:
                prot_chain, prot_resname1, prot_resid1 = res1
                resname2 = stacks[n2][0]
                for res2 in stacks[n2][1:]:
                    prot_chain2, prot_resname2, prot_resid2 = res2
                    if prot_chain2 != prot_chain:
                        continue
                    if prot_resid1 == prot_resid2:
                        continue
                    print(
                        code,
                        chain,
                        resname1,
                        n,
                        resname2,
                        n2,
                        prot_chain,
                        prot_resname1,
                        prot_resid1,
                        prot_resname2,
                        prot_resid2,
                    )


for file in glob.glob("../database/3dnaj/*-1-dssr.json"):
    pos1 = file.rindex("/") + 1
    pos2 = len(file) - len("-1-dssr.json")
    code = file[pos1:pos2]
    find_double_stacking(code, file)
