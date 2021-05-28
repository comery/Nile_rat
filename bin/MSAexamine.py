#!/usr/bin/env python3

import sys
import json
import argparse

def getLen(length_file):
    lens = {}
    with open(length_file, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            lens[tmp[0]] = int(tmp[1])
    return lens


def Add_records(dict, k, v):
    if k in dict.keys():
        dict[k].append(v)
    else:
        dict[k] = [v,]

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "\n".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "\n".join(seq))

def alignrate(s1, s2):
    slen = len(s1)
    aligned = 0
    for i in range(slen):
        if s1[i] != "-" and s2[i] != "-":
            aligned += 1
    return aligned/slen


def get_alignrate(ref_species, ref_seq, sequence):
    rates = {}
    for s in sequence.keys():
        rates[s] = alignrate(ref_seq, sequence[s])
    return rates

def get_stop_condon(ref_species, ref_seq, sequence):
    stop_codon = {}
    stop_codon[ref_species] = ref_seq.count("U")
    for s in sequence.keys():
        stop_codon[s] = sequence[s].count("U")
    return stop_codon

def generate_report(gene, align_rate, stop_condon_number, length, out_format):
    min_align_rate = sorted(align_rate.values())[0]
    all_stop = sum(stop_condon_number.values())
    if min_align_rate >= args.rate and all_stop == 0 and length[0] >= args.minL and length[0] <= args.maxL:
        sign = 'PASS'
    else:
        sign = 'Filtered'

    if out_format == 'tab':
        t1 = []
        for a in align_rate.keys():
            t1.append(a + ":{:.4f}".format(align_rate[a]))
        t2 = []
        for a in stop_condon_number.keys():
            t2.append(a + ":{}".format(stop_condon_number[a]))
        print(gene + "\t" + sign + "\t" + "\t".join(t1) + "\t" + "\t".join(t2) + "\t" + str(length[0]))

    if out_format == 'json':
        results = {
            'gene': gene,
            'status': sign,
            'align_rate_vs_ref': align_rate,
            'stop_codon_number': stop_condon_number,
            'alignment_length': length[0],
        }
        json_output(results)

def json_output(data):
        #json.dump(data, sys.stdout)
        print(json.dumps(data, indent=4))

def main(args):
    sequence = {}
    with open(args.alignment, 'r') as fh:
        length = []
        for name, seq in read_fasta(fh):
            info = name.replace(">", "").split("_")
            species = info[0]
            gene = info[1]
            length.append(len(seq))
            if species == args.rid:
                ref_species = species
                ref_seq = seq
            else:
                sequence[species] = seq

        align_rate = {}
        stop_condon_number = {}
        length.sort()
        if length[0] != length[-1]:
            sys.exit("uneven length, this is not a alignment")
        align_rate = get_alignrate(ref_species, ref_seq, sequence)
        stop_condon_number = get_stop_condon(ref_species, ref_seq, sequence)
        generate_report(gene, align_rate, stop_condon_number, length, args.out_format)

if __name__ == '__main__':
    des = """
    Check protein alignment status by:
    print(all_stop)
        - align rate [0.8]
        - stop codon [none stop]
        - alignment length [minL, maxL]
    and give a short report [tab, json]

    yangchentao at genomics.cn, BGI.
    """
    parser = argparse.ArgumentParser(description=des,
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-a", type=str, required=True, metavar='<FILE>',
                        dest='alignment', help="alignment file")
    parser.add_argument("-f", type=str, required=False, metavar='<STR>',
                        choices=['tab', 'json'], default='json',
                        dest='out_format', help="output format, tab, json")
    parser.add_argument("-rid", type=str, required=True, metavar='<STR>',
                        help="species ID for reference ")
    parser.add_argument("-minL", type=int, required=False, metavar='<INT>',
                        help="min alignment length required, default=50",
                        default=50)
    parser.add_argument("-maxL", type=int, required=False, metavar='<INT>',
                        default=5000, help="max alignment length allowed, default=5000")
    parser.add_argument("-rate", type=float, required=False,
                        metavar='<FLOAT>', default=0.8,
                        help="min aligning rate allow")
    args =  parser.parse_args()
    main(args)
