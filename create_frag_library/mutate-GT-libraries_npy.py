#!/usr/bin/python3
# Copyright 2017 - 2018 Isaure Chauvot de Beauchene (CNRS)

import sys, os, json
from subprocess import check_call

na = sys.argv[1]
pattern = sys.argv[2]

print((na, pattern), file=sys.stderr)

dirname = os.path.dirname(os.path.abspath(__file__))

def execute(command):
    print(command)
    check_call(command, shell = True)

x = "D"
dict_mutations = {'G':'A', 'T':'C'}

s = ['G','T']
sequences = [a+b+c for a in s for b in s for c in s ]
m = [0,1]
mutations = [(a, b, c) for a in m for b in m for c in m ]
mutations = mutations[1:]

diffatoms={
 'A':["N6 "],
 'G':["O6 ","N2 "],
 'U':["O4 "],
 'C':["N4 "],
 'T':["O4", "C7"]
 }

for sequence in sequences:
    seq = [ s for s in sequence ]
    print(sequence)
    os.system('mkdir -p /tmp/%s'%sequence)
    execute('%s/npy2pdb.py %s-%s-aa.npy templates/%s.pdb > /tmp/%s.pdb' \
             %(dirname, sequence, pattern, sequence, sequence))
    for mut in mutations:
        mutant = ""
        for s, m in zip(seq ,mut):
            if m:
                c = dict_mutations[s]
            else:
                c = s
            mutant = mutant + c
        name = "%s-%s"%(mutant, pattern)
        os.system('mkdir /tmp/%s; mkdir %s'%(mutant, mutant))
        #mutate
        execute('cp /tmp/%s.pdb /tmp/%s.pdb'%(sequence, mutant))
        for i in range(3):
            execute("sed -i 's/%s%s     %i/%s%s     %i/' /tmp/%s.pdb" \
                    %(x, seq[i], i+1, x, mutant[i], i+1, mutant))
        execute('%s/splitmodel.py /tmp/%s.pdb /tmp/%s/conf-aa > /tmp/%s-aa.list' \
                 %(dirname, mutant, mutant, name))
        execute("sed 's/\/tmp\///' /tmp/%s-aa.list > %s-aa.list" \
                %(name,  name))
        execute('%s/pdbcompletion.py /tmp/%s-aa.list %s-aa.list \
                --batch --heavy --%s --patch None 1 None > /dev/null 2> /dev/null' \
                %(dirname, name, name, na))
        #reduce
        execute("sed 's/\-aa/r/' %s-aa.list > %sr.list "\
                %(name, name))
        execute('%s/pdb2npy.py %s-aa.list --list --outp %s-aa.npy'\
                %(dirname, name, mutant))
        execute('%s/reduce-npy.py %s-aa.npy templates/%s.pdb %sr.npy' \
                 %(dirname, mutant, mutant, mutant))
        #cleanup
        execute("rm -rf /tmp/%s*"%mutant)
