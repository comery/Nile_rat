#!/usr/bin/env python3
import os, sys, re
import numpy as np
from collections import defaultdict
if len(sys.argv) < 2:
    sys.exit("python3 {} <codeml.out>".format(sys.argv[0]))

codeml = open(sys.argv[1], "r") # codeml pairwise ML output here

status = 0
genome_dnds = defaultdict(list)
#print "first\tsecond\tdnds\tdn\tds"
output = []
for i in codeml.readlines():
    if i.startswith("pairwise comparison, codon frequencies"):
        status = 1

    if status == 1:
        if i[0].isdigit():
            line = i.rstrip()
            line2 = re.sub("\(", "", line)
            line3 = re.sub("\)", "", line2)
            spaces = line3.split(" ")
            first = spaces[1]
            second = spaces[4]

            first_split = first.split("..")
            second_split = second.split("..")
            g1 = first_split[0]
            g2 = second_split[0]

        if i.startswith("t="):
            line = i.rstrip()
            line1 = re.sub("=", "", line)
            line2 = re.sub("\s+", "\t", line1)
            tabs = line2.split("\t")
            dnds = tabs[7]
            dn = tabs[9]
            ds = tabs[11]
            #if float(ds) < 2 and float(ds) > 0.01 and float(dnds) < 10:
            tmp = f"{first}\t{second}\t{dnds}\t{dn}\t{ds}"
            output.append(tmp)
if len(output) > 1:
    print("Gene_1\tGene_2\tdnds\tdN\tdS")
    for i in output:
        print(i)
