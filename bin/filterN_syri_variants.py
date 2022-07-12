#!/usr/bin/env python3
import sys
import re
if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {} ".format(sys.argv[0], "syri.out",
                                                   "ref.fa", "querry.fa"))
    exit()


def read_fasta(fa):
    fasta_dict = {}
    name = ""
    seq = ""
    with open(fa, 'r') as fh:
        for i in fh:
            if i.startswith(">"):
                if name:
                    fasta_dict[name] = seq
                else:
                    name = i.strip().replace(">", "")
            else:
                seq += i.strip()
        if name:
            fasta_dict[name] = seq

    return fasta_dict


def Ncontent(sequence):
    sequence = sequence.upper()
    seqlen = len(sequence)
    n = sequence.count("N")
    n_content = n / seqlen
    return n_content


def shorter(line):
    new = []
    tmp = line.strip().split("\t")
    new = tmp[0:3] + tmp[5:9]
    new.append(tmp[10])
    return "\t".join(new)


def main():
    ref_seqs = read_fasta(sys.argv[2])
    que_seqs = read_fasta(sys.argv[3])

    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split("\t")
            ref_id = tmp[0]
            ref_start = int(tmp[1])
            ref_end = int(tmp[2])
            que_id = tmp[5]
            que_start = int(tmp[6])
            que_end = int(tmp[7])
            varseq_ref = ref_seqs[ref_id][ref_start-1:ref_end]
            varseq_que = que_seqs[que_id][que_start-1:que_end]
            ref_n_content = Ncontent(varseq_ref)
            que_n_content = Ncontent(varseq_que)
            if ref_n_content <= 0.5 and que_n_content <= 0.5:
                print(shorter(i))


if __name__ == '__main__':
    main()

