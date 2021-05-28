#!/usr/bin/env python3
import sys
import os
import subprocess

def get_seqid(aligned_seq):
    seqid = {}
    with open(aligned_seq, 'r') as fh:
        for i in fh:
            if i.startswith(">"):
                long_id = i.replace(">", "").strip()
                tmp = long_id.split("_")
                tax = tmp[0]
                seqid[tax] = long_id

    return seqid


def changeSpeciesShortNamesWithGeneNames(seqid, tree):

    for sp in seqid.keys():
        tree = tree.replace(sp, seqid[sp])

    return(tree)

def preTree(outdir, seqid, tree):
    with open(tree, 'r') as tr:
        tree_str = tr.read()
    for hy in ['H0', 'H1']:
        o = open(outdir + "/" + hy + ".tree", "wt")
        o.write(changeSpeciesShortNamesWithGeneNames(seqid, tree_str) + "\n")
        o.close()


def setDefaultCtl(model, hy, o):
    if model == 'branch':
        if hy == 'H0':
            o.write("noisy = 3\nverbose = 1\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\nmodel = 0\nNSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 2.5\nfix_omega = 0\nomega = 0.2\n")
        elif hy == 'freeRatio':
            o.write("noisy = 3\nverbose = 1\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\nmodel = 1\nNSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 2.5\nfix_omega = 0\nomega = 0.2\n")
        else:
            o.write("noisy = 3\nverbose = 1\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\nmodel = 2\nNSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 2.5\nfix_omega = 0\nomega = 0.2\n")
    if model == 'branchSite':
        if hy == 'H0':
            o.write("noisy = 3\nverbose = 1\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\nmodel = 2\nNSsites = 2\nicode = 0\nfix_kappa = 0\nkappa = 2\nfix_omega = 1\nomega = 1\n")
        else:
            o.write("noisy = 3\nverbose = 1\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\nmodel = 2\nNSsites = 2\nicode = 0\nfix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 1.5\n")


def preCTL(outdir, seq):
    with open(outdir + "/H0.ctl", 'w') as o:
        o.write("seqfile = " + seq  + "\ntreefile = " +  "H0.tree\noutfile = " + "H0.mlc\n")
        setDefaultCtl('branchSite', 'H0', o)

    with open(outdir + "/H1.ctl", 'w') as o:
        o.write("seqfile = " + seq  + "\ntreefile = " +  "H1.tree\noutfile = " + "H1.mlc\n")
        setDefaultCtl('branchSite', 'H1', o)


def preShell(outdir):
    with open(outdir + "/paml.sh", 'w') as fh:
        fh.write("/hwfssz1/ST_DIVERSITY/PUB/USER/panhailin/bin/PAML20191217/PAML/bin/paml4.9i/bin/codeml H0.ctl\n")
        fh.write("/hwfssz1/ST_DIVERSITY/PUB/USER/panhailin/bin/PAML20191217/PAML/bin/paml4.9i/bin/codeml H1.ctl")

def main():
    tree = sys.argv[1]
    aligned_seq = sys.argv[2]
    outdir = sys.argv[3]
    outdir = os.path.abspath(outdir)
    aligned_seq = os.path.abspath(aligned_seq)

    tree = os.path.abspath(tree)
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    seqid = get_seqid(aligned_seq)
    preTree(outdir, seqid, tree)
    preCTL(outdir, aligned_seq)
    preShell(outdir)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("python3 {} <*.tree> <*.cds.align> <outdir>".format(sys.argv[0]))

    main()
