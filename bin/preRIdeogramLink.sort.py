#!/usr/bin/env python3

import sys

# Import palettable
try:
    from palettable.tableau import GreenOrange_12
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)

try:
    import palettable.cartocolors.qualitative
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)


def getLen(length_file):
    lens = {}
    with open(length_file, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            lens[tmp[0]] = int(tmp[1])
    return lens

def give_color():
    """ colors theme  """""
    colors = []
    print(
            "colors theme is designed by palettable." +\
            " \nsee more:http://colorbrewer2.org\n"
    )

    colors = palettable.cartocolors.qualitative.Antique_10.colors
    colors += palettable.cartocolors.qualitative.Bold_10.colors
    colors += palettable.cartocolors.qualitative.Pastel_10.colors
    colors += palettable.cartocolors.qualitative.Prism_10.colors
    colors += palettable.cartocolors.qualitative.Safe_10.colors
    colors += palettable.cartocolors.qualitative.Vivid_10.colors
    colors += palettable.cartocolors.diverging.ArmyRose_7.colors
    colors += palettable.cartocolors.diverging.Earth_7.colors
    colors += palettable.cartocolors.diverging.Fall_7.colors
    colors += palettable.cartocolors.diverging.Geyser_7.colors
    colors += palettable.cartocolors.diverging.TealRose_7.colors
    colors += palettable.cartocolors.diverging.Temps_7.colors
    colors += palettable.cartocolors.diverging.Tropic_7.colors
    colors += palettable.cartocolors.diverging.Geyser_7.colors
    colors += palettable.tableau.BlueRed_12.colors
    colors += palettable.tableau.ColorBlind_10.colors
    colors += palettable.tableau.PurpleGray_12.colors
    return colors

def Add_records(dict, k, v):
    if k in dict.keys():
        dict[k].append(v)
    else:
        dict[k] = [v,]

def main():
    records = []
    lens = getLen(sys.argv[2])
    colors = give_color()
    out = open("synteny_comparison.txt", 'w')
    species1_kary = {}
    species2_kary = {}
    with open(sys.argv[1], 'r') as fh:
        for b in fh:
            tmp = b.strip().split()
            if len(species1_kary.keys()) == 0:
                species1_kary[tmp[0]] = 1
                species1_current = 1
            elif tmp[0] not in species1_kary.keys():
                species1_current += 1
                species1_kary[tmp[0]] = species1_current

            if len(species2_kary.keys()) == 0:
                species2_kary[tmp[3]] = 1
                species2_current = 1
            elif tmp[3] not in species2_kary.keys():
                species2_current += 1
                species2_kary[tmp[3]] = species2_current
            out.write("{}\t{}\t{}\t{}\t{}\t{}\tcccccc\n".format(species1_kary[tmp[0]], tmp[1], tmp[2],
                                                                species2_kary[tmp[3]], tmp[4], tmp[5]))
    out.close()
    with open("RIdeogram.karyotype.txt", 'w') as ka:
        ka.write("Chr\tStart\tEnd\tfill\tspecies\tsize\tcolor\n")
        if len(species1_kary.keys()) + len(species2_kary.keys()) > 150:
            print(
                "Sorry, It just can handle less than 150 chrs because " +\
                "too many chrs will make a mass!, however, you can change" +\
                "this code for suit\n"
                )
            exit(0)

        sorted_s1 = sorted(species1_kary.keys(), key=lambda k:species1_kary[k])
        for s1 in sorted_s1:
            order = species1_kary[s1]
            r, g, b = colors[order]
            color = '{:02x}{:02x}{:02x}'.format(r, g, b)
            ka.write("{}\t1\t{}\t{}\tmaternal\t12\t969696\n".format(order, lens[s1], color))
        sorted_s2 = sorted(species2_kary.keys(), key=lambda k:species2_kary[k])
        for s2 in sorted_s2:
            order = species2_kary[s2]
            r, g, b = colors[-order]
            color = '{:02x}{:02x}{:02x}'.format(r, g, b)
            ka.write("{}\t1\t{}\t{}\tpaternal\t12\t969696\n".format(order, lens[s2], color))

    with open("RIdeogram.plot.R", 'w') as rh:
        plot = """library(RIdeogram)
karyotype <- read.table("RIdeogram.karyotype.txt", sep = "\\t", header = T, stringsAsFactors = F)
comparison <- read.table("synteny_comparison.txt", sep = "\\t", header = F, stringsAsFactors = F)
names(comparison) <- c('Species_1', 'Start_1', 'End_1', 'Species_2', 'Start_2', 'End_2', 'fill')
ideogram(karyotype = karyotype, synteny = comparison, output = "synteny_comparison.svg")
library(rsvg)
rsvg_pdf("synteny_comparison.svg", "synteny_comparison.pdf")
        """
        rh.write(plot)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("python3 {} <*.link> <all.lens>".format(sys.argv[0]))
    main()
