import json, sys, os
import argparse

'''

'''

parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("file", type=str, help="File to check and modify.")
args = parser.parse_args()

# Put input in variable names
checked_file=args.file

list_to_remove = []

f = open(checked_file, 'r')
lines = f.readlines()
f.close()

f = open(checked_file, 'w')

for line in lines:
    if line[22:26] not in list_to_remove:
        for coord in line[30:54]:
            if coord == 'X':
                list_to_remove.append(line[22:26])
                break

for line in lines:
    if line[22:26] not in list_to_remove:
        f.write(line)

f.close()
