#!/usr/bin/env python3

import os
import sys

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} all.Assemblytics_structural_variants.filterN.bed pat.sv.mrna.overlap.bed mat.sv.mrna.overlap.bed ")


def read_bed_annot(infile):
    records = {}
    with open(infile, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            bed = "\t".join(tmp[0:3])
            records[bed] = i.strip()

    return records


def main():
    original = sys.argv[1]
    ref_annot = read_bed_annot(sys.argv[2])
    que_annot = read_bed_annot(sys.argv[3])
    with open(original, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            ref_bed = "\t".join(tmp[0:3])
            que_bed = "\t".join(tmp[3:6])
            non = f".\t.\t.\t."
            info = "\t".join(tmp[6:])
            new_record = []
            if ref_bed in ref_annot.keys():
                new_record = [ref_annot[ref_bed], ]
            else:
                new_record = [non, ]
            if que_bed in que_annot.keys():
                new_record.append(que_annot[que_bed])
            else:
                new_record.append(non)
            new_record.append(info)
            print("\t".join(new_record))


if __name__ == '__main__':
    main()

