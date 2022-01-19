#!/usr/bin/env python3
import os
import sys
import subprocess
if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {} ".format(sys.argv[0], "Assemblytics_structural_variants.bed",
                                                   "mat.fa", "pat.fa"))
    exit()


def read_one_fa(fa):
    name = ""
    seq = ""
    with open(fa, 'r') as fh:
        for i in fh:
            if i.startswith(">"):
                name = i.strip().replace(">", "")
            else:
                seq += i.strip()
    return name, seq

def out2fa(f, string):
    with open(f, 'w') as fw:
        fw.write(string)

def Ncontent(sequence):
    sequence = sequence.upper()
    seqlen = len(sequence)
    n = sequence.count("N")
    n_content = n / seqlen
    return n_content

def get_pos(info):
    tmp = info.split(":")
    tmp1 = tmp[1].split("-")
    return tmp[0], int(tmp1[0]), int(tmp1[1])

def main():
    ref_name, ref_seq = read_one_fa(sys.argv[2])
    ref_name, que_seq = read_one_fa(sys.argv[3])
    ref_len = len(ref_seq)
    que_len = len(que_seq)
    with open(sys.argv[1], 'r') as fh:
        for r in fh:
            if r.startswith("#"):
                continue
            tmp = r.strip().split("\t")
            svtype = tmp[6]
            ref_id = tmp[0]
            ref_start = int(tmp[1])
            ref_end = int(tmp[2])
            ref_varlen = abs(ref_end - ref_start)
            que_id, que_start, que_end = get_pos(tmp[9])
            que_varlen = abs(que_end - que_start)
            varseq_ref = ref_seq[ref_start-1:ref_end]
            varseq_que = que_seq[que_start-1:que_end]
            ref_n_content = Ncontent(varseq_ref)
            que_n_content = Ncontent(varseq_que)
            if ref_n_content <= 0.5 and que_n_content <= 0.5:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(ref_id, ref_start, ref_end,
                                                                 que_id, que_start, que_end,
                                                                 svtype, tmp[4], tmp[10]))


if __name__ == '__main__':
    main()
