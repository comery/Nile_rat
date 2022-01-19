#!/usr/bin/env python3
import sys
import os
if len(sys.argv) < 2:
    sys.exit("python3 {} <$a.1coords>".format(sys.argv[0]))


def main():
    link = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            pair = tmp[-2] + "-" + tmp[-1]
            if pair not in link.keys():
                link[pair] = [tmp[0:4],]
            else:
                link[pair].append(tmp[0:4])


    for x in link.keys():
        ids = x.split("-")
        target_start = link[x][0][0]
        target_end = link[x][-1][1]
        que_start = link[x][0][2]
        que_end = link[x][-1][-1]
        out = [ ids[0], target_start, target_end, ids[1], que_start, que_end ]
        print("\t".join(out))



if __name__ == '__main__':
    main()
