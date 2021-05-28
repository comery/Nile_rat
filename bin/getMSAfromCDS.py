import sys
import os
usage = """
deal wiht all.cds file, generate a species_gene matrix and
split single gene into subdir by gene name.
"""
if len(sys.argv) < 3:
    print(usage)
    sys.exit("python3 {} <all.cds> <outdir>".format(sys.argv[0]))

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

def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})

def main():
    outdir = sys.argv[2]
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    DATA = {}
    sequences = {}
    species = set()
    with open(sys.argv[1], 'r') as fp:
        for name, seq in read_fasta(fp):
            name = name.replace(">", "")
            name = name.replace("-D0", "")
            [species_tag, gene_id] = name.split("_")
            species.add(species_tag)
            addtwodimdict(DATA, gene_id, species_tag, name)
            sequences[name] = seq
    species_list = list(species)

    with open("ortholog.matrix.txt", 'w') as mh:
        gene_count = 0 # total gene number which should be saved
        loop = 200  # how many genes in a subdir
        dir_count = 0 # how many subdirs
        print("#Gene\t" + "Path\t" + "\t".join(species_list), file=mh)
        for gene in DATA.keys():
            tmp_out = [gene,]
            absent = 0 # to check whether all species have this gene
            for sp in species_list:
                if sp in DATA[gene].keys():
                    tmp_out.append(DATA[gene][sp])
                else:
                    tmp_out.append('-')
                    absent += 1
            # if all species have this gene, then ouput their sequence into subdir
            if absent == 0:
                gene_count += 1
                if gene_count % loop == 1: # for example: 200 to 201, 001 to 002
                    dir_count += 1
                    threeLetterDir = "{:03d}".format(dir_count)

                subdir = outdir + "/" + threeLetterDir + "/" + gene
                tmp_out.insert(1, subdir)
                if os.path.exists(subdir) == False:
                    os.makedirs(subdir)
                with open(subdir + "/" + gene + ".cds", 'w') as fh:
                    for s in tmp_out[2:]:
                        print(">" + s + "\n" + sequences[s], file=fh)
            else:
                tmp_out.insert(1, 'na')
            # output specie_genename into a matrix
            print("\t".join(tmp_out), file=mh)


if __name__ == '__main__':
    main()
