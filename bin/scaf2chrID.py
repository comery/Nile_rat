#!/usr/bin/env python3
import sys
if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {}".format(sys.argv[0], "scaf2chr.txt", "genome.fsa", "prefix"))
    exit()


accs = {}
prefix = sys.argv[3]
with open(sys.argv[1], "r") as ac:
    for i in ac:
        if i.startswith("#") or i.startswith("-"):
            continue
        tmp = i.strip().split(",")
        accs[tmp[0]] = tmp[1]

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, '\n'.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, '\n'.join(seq))


with open(sys.argv[2], 'r') as fp:
    for name, seq in read_fasta(fp):
        name = name.replace(">", "")
        seq_id = name.split(" ")[0]
        if seq_id in accs.keys():
            if "X" in accs[seq_id] or 'Y' in accs[seq_id]:
                continue
            else:
                print(">" + prefix + "_Chr" + accs[seq_id])
                print(seq)
