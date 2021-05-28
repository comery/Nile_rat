import sys

if len(sys.argv) < 2:
    sys.exit("python3 {} *.indel".format(sys.argv[0]))
def parse_pos(string):
    if "-" in string:
        pos1 = int(string.split("-")[0])
        pos2 = int(string.split("-")[1])
    else:
        pos1 = int(string)
        pos2 = pos1

    return (pos1, pos2)

with open(sys.argv[1], 'r') as fh:
    outi = open("ins.mat.pos.bed", 'w')
    outd = open("del.pat.pos.bed", 'w')
    for i in fh:
        if i.startswith("#"):
            continue
        tmp = i.strip().split("\t")
        if tmp[4] == 'ins':
            ref_pos1 = int(tmp[2].replace("_", ""))
            ref_pos2 = ref_pos1 + 1
            (que_pos1, que_pos2) = parse_pos(tmp[3])
            print(f"{tmp[1]}\t{que_pos1}\t{que_pos2}", file=outi)
        else:
            que_pos1 = int(tmp[3].replace("_", ""))
            que_pos2 = que_pos1 + 1
            (ref_pos1, ref_pos2) = parse_pos(tmp[2])
            print(f"{tmp[0]}\t{ref_pos1}\t{ref_pos2}", file=outd)

    outi.close()
    outd.close()
