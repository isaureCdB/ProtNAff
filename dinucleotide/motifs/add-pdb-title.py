import requests, sys
import glob

for f in glob.glob("*.list"):
    f2 = f[:-5] + "-titles.list"
    print(f)
    with open(f2, "w") as fp2:
        for code in open(f):
            code = code.strip()
            r = requests.get(f'http://rcsb.org/structure/{code}')
            al = r.text
            title = al[al.find('<title>') + 7 : al.find('</title>')][17:]
            print(code, title, file=fp2)