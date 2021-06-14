#!/usr/bin/env python

import sys
f = sys.argv[1]
template = f[:-4]
if len(sys.argv) > 2: template = sys.argv[2]
size = None
if len(sys.argv) > 3: size = int(sys.argv[3])

header = []
count = 0
for l in open(f):
  l = l.rstrip("\n")
  if l.startswith("MODEL"):
    if count > 0: ff.close()
    count += 1
    target = template+"-"+str(count)+".pdb"
    print(target)
    ff = open(target, "w")
    for h in header: ff.write(h + "\n")
    continue
  if l.startswith("ENDMDL"): continue
  if count == 0: header.append(l)
  else: ff.write(l + "\n")

if size is not None and count != size:
  raise Exception("ERROR in splitmodel: model count should be %d, is %d" % (size, count))
