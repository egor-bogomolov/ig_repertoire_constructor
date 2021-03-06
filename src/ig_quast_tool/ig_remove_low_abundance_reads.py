#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys


import contextlib
@contextlib.contextmanager
def smart_open(filename, mode="r"):
    """
    From http://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    """
    import gzip
    import re
    from sys import stdout, stdin

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE= "a"
    else:
        MODE = "r"

    if filename != '-':
        if re.match(r"^.*\.gz$", filename):
            assert(MODE != "a")
            fh = gzip.open(filename, mode=MODE)
        else:
            fh = open(filename, mode=mode)
    else:
        assert(MODE != "a")
        fh = stdout if MODE == "w" else stdin
    try:
        yield fh
    finally:
        if fh is not stdout and fh is not stdin:
            fh.close()


def parse_abundance(s):
    import re

    m = re.match(r"^.*size___(\d+)$", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

assert(parse_abundance("size___10") == 10)


if __name__ == "__main__":
    parser = ArgumentParser(description="Remove low abundance reads")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA file")
    parser.add_argument("--limit", "-l",
                        type=int,
                        default=5,
                        help="abundance limit (default: %(default)s)")

    args = parser.parse_args()

    print "Command line: %s" % " ".join(sys.argv)

    result = []
    num_all_clusters = 0
    with smart_open(args.input, "r") as fin:
        for record in SeqIO.parse(fin, "fasta"):
            abundance = parse_abundance(str(record.id))
            num_all_clusters += 1
            if abundance >= args.limit:
                result.append(record)

    print str(len(result)) + " antibody clusters have abundance >= " + str(args.limit)
    print str(num_all_clusters - len(result)) + " lowly abundant antibody clusters will be discarded"

    with smart_open(args.output, "w") as fout:
        SeqIO.write(result, fout, "fasta")

    print "Highly abundant clusters were written to " + args.output
