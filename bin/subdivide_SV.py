#!/usr/bin/env python3
import re
import sys
if len(sys.argv) < 3:
    sys.exit("Usage: python3 {} syri.sv.out outpre".format(sys.argv[0]))

def shorter(line):
    new = []
    tmp = line.strip().split("\t")
    new = tmp[0:3] + tmp[5:9]
    new.append(tmp[10])
    return "\t".join(new) + "\n"

def main():
    with open(sys.argv[1], 'r') as fh:
        cg = open(sys.argv[2] + ".copygain.tab", 'w')
        cl = open(sys.argv[2] + ".copyloss.tab", 'w')
        inv = open(sys.argv[2] + ".inv.tab", 'w')
        tra = open(sys.argv[2] + ".trans.tab", 'w')
        invtr = open(sys.argv[2] + ".invtr.tab", 'w')
        for i in fh:
            if i.startswith("#"): continue
            tmp = i.strip().split("\t")
            al = re.compile(r'AL$')
            if al.search(tmp[10]) or tmp[10] == "SYN":
                continue
            if tmp[11] == "copygain":
                cg.write(shorter(i))
            elif tmp[11] == "copyloss":
                cl.write(shorter(i))
            elif tmp[10] == "INV":
                inv.write(shorter(i))
            elif tmp[10] == "INVTR":
                invtr.write(shorter(i))
            elif tmp[10] == "TRANS":
                tra.write(shorter(i))
            else:
                print(shorter(i))
    cg.close()
    cl.close()
    inv.close()
    tra.close()
    invtr.close()

if __name__ == '__main__':
    main()

