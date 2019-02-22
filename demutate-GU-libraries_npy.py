#!/usr/bin/python2.7

import sys, os, json
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
from subprocess import check_call
import rmsdlib

def execute(command):
    print(command)
    check_call(command, shell = True)

def create_template(seq):
    execute('cat /dev/null > templates/%s.pdb'%seq)
    for i in range(3):
        execute("sed 's/D%s     1/D%s     %i/' $ATTRACTTOOLS/../allatom/%slib/%s.pdb >> templates/%s.pdb"%(seq[i], seq[i], i+1, na, seq[i], seq))

na = sys.argv[1]

m = [0,1]
s = ['G','U']
if na == 'dna':
    s = ['G','T']

mutations = [(a, b, c) for a in m for b in m for c in m ]
mutations = mutations[1:]
sequences = [a+b+c for a in s for b in s for c in s ]

###mutations = [(1, 0, 0)] ###
dict_mutations={
  'G':'A',
  'U':'C',
  'T':'C'
  }

diffatoms={
 'A':["N6 "],
 'G':["O6 ","N2 "],
 'U':["O4 "],
 'C':["N4 "],
 'T':["O4", "C7"]
 }

for sequence in sequences:
    #if sequence == 'GUU' : continue #TODO: remove this line!
#for sequence in ['GUU']:
    seq = [ s for s in sequence ]
    print(sequence)
    os.system('mkdir -p /tmp/%s'%sequence)
    #create_template(sequence)
    execute('$SCRIPTS/npy2pdb.py %s-2nd.npy templates/%s.pdb > /tmp/%s.pdb'%(sequence, sequence, sequence))
    for mut in mutations:
    #for mut in [(0,1,0)]:
        mutant = ""
        for s, m in zip(seq ,mut):
            if m:
                c = dict_mutations[s]
            else:
                c = s
            mutant = mutant + c
        print mutant
        create_template(mutant)
        execute('cp /tmp/%s.pdb /tmp/%s.pdb'%(sequence, mutant))
        for i in range(3):
            if na == 'dna':
                execute("sed -i 's/D%s     %i/D%s     %i/' /tmp/%s.pdb"%(seq[i], i+1, mutant[i], i+1, mutant))
            else:
                execute("sed -i 's/R%s     %i/R%s     %i/' /tmp/%s.pdb"%(seq[i], i+1, mutant[i], i+1, mutant))
        #mutate
        os.system('mkdir /tmp/%s; mkdir %s'%(mutant, mutant))
        execute('$ATTRACTTOOLS/splitmodel /tmp/%s.pdb /tmp/%s/conf-aa > /tmp/%s-2nd.list'%(mutant, mutant, mutant))
        execute("sed 's/\/tmp\///' /tmp/%s-2nd.list > %s-2nd.list"%(mutant, mutant))
        execute('python2 $ATTRACTTOOLS/../allatom/aareduce.py /tmp/%s-2nd.list %s-2nd.list --batch --heavy --%s --nalib --patch 1 None > /dev/null 2> /dev/null'%(mutant, mutant, na))
        #reduce
        execute("sed 's/\-aa/r/' %s-2nd.list > %s-2ndr.list "%(mutant, mutant))
        execute('$SCRIPTS/pdb2npy.py %s-2nd.list --list --outp %s-2nd.npy'%(mutant, mutant))
        execute('python2 $ATTRACTTOOLS/reduce-npy.py %s-2nd.npy templates/%s.pdb %s-2ndr.npy'%(mutant, mutant, mutant))
        #cleanup
        execute("rm -rf /tmp/%s*"%mutant)
