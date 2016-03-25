#!/usr/bin/env python2

import igblast_utils
import os
import os.path
import csv
from Bio import SeqIO
import time
import sys
from argparse import ArgumentParser
import tempfile
import gzip


class FakeLog:
    def info(self, s):
        print s

def linear_search(obj, item, start=0):
    for i in range(start, len(obj)):
        if obj[i] == item:
            return i
    return -1

def idFormatByFileName(fname):
    import re
    if re.match(r"^.*\.fa(sta)?(\.gz)?$", fname):
        return "fasta"
    elif re.match(r"^.*\.((fq)|(fastq))(\.gz)?$", fname):
        return "fastq"
    else:
        raise "Unrecognized file type"


assert idFormatByFileName("fq.fa") == "fasta"
assert idFormatByFileName("fq.fa.gz") == "fasta"
assert idFormatByFileName("fq.fasta") == "fasta"
assert idFormatByFileName("fq.fasta.gz") == "fasta"
assert idFormatByFileName("fq.fq.gz") == "fastq"

def ilen(iterable):
    return sum(1 for i in iterable)


def fastX_len(fname):
    with smart_open(fname, "rU") as fh:
        return ilen(SeqIO.parse(fh, idFormatByFileName(fname)))

def fa2fq(input_file, output_file):
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, idFormatByFileName(input_file)):
            # record.letter_annotations["phred_quality"] = [40] * len(record)
            record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]
            SeqIO.write(record, fout, "fastq")

def fq2fa(input_file, output_file):
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        parser = SeqIO.parse(fh, idFormatByFileName(input_file))
        SeqIO.write(parser, fout, "fasta")


class HitTableRowVJF:
    def __init__(self, _type, query_id, subject_id, start, end):
        self.type = _type
        self.query_id = query_id
        self.subject_id = subject_id
        self.start = start
        self.end = end

def parse_vjf_output(filename, readfile):
    from collections import defaultdict

    with smart_open(readfile, "rU") as fh:
        parser = SeqIO.parse(fh, idFormatByFileName(readfile))
        descr_to_ind = { str(record.description).replace(" ", "_"): i for i, record in enumerate(parser) }

    result = defaultdict(dict)
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        headers = reader.next()

        id_col = linear_search(headers, "id")
        Vstart_col = linear_search(headers, "Vstart")
        Vend_col = linear_search(headers, "Vend")
        Vgene_col = linear_search(headers, "Vid")
        Jgene_col = linear_search(headers, "Jid")
        Jstart_col = linear_search(headers, "Jstart")
        Jend_col = linear_search(headers, "Jend")
        for line in reader:
            # print line[Vstart_col], line[Jend_col]
            # print Vstart_col
            # print line
            desc = line[id_col]

            Vstart = int(line[Vstart_col])
            Vend = int(line[Vend_col])
            Jstart = int(line[Jstart_col])
            Jend = int(line[Jend_col])

            Vgene = line[Vgene_col]
            Vgene = Vgene[:Vgene.find(" ")]
            Jgene = line[Jgene_col]
            Jgene = Jgene[:Jgene.find(" ")]

            ind = descr_to_ind[desc]
            result[desc]["V"] = HitTableRowVJF("V", desc, Vgene, Vstart, Vend)
            result[desc]["J"] = HitTableRowVJF("J", desc, Jgene, Jstart, Jend)
            result[ind] = result[desc]

        return result


class Empty:
    pass

def smart_open(filename, mode="r"):
    import gzip
    import re

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE= "a"
    else:
        MODE = "r"

    if re.match(r"^.*\.gz$", filename):
        assert(MODE != "a")
        fh = gzip.open(filename, mode=MODE)
    else:
        fh = open(filename, mode=mode)
    return fh

def md5_file(fname):
    import hashlib

    hash_md5 = hashlib.md5()
    with smart_open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def hash_file(fname, len=7):
    import os.path

    hash = md5_file(fname)
    base = os.path.basename(fname)

    return base.split(".")[0] + "_" + hash[:len]


if __name__ == "__main__":
    parser = ArgumentParser(description="Benchmark VJFinder vs IgBLAST")
    parser.add_argument("input",
                        type=str,
                        default="./SRA/SRR1383463.fq", #.gz
                        help="input FASTQ file (default: %(default)s)")
    parser.add_argument("--bad-reads", "-b",
                        type=str,
                        default="",
                        help="output FASTQ file for suspicious reads, <empty> for non-producing (default: <empty>)")
    # parser.add_argument("--path", "-p",
    #                     type=str,
    #                     default="/home/ashlemov/Git/ig_repertoire_constructor/",
    #                     help="path to installed IgReC (default: %(default)s)")
    parser.add_argument("--options", "-o",
                        type=str,
                        default="-t 16",
                        help="additional options for VJFinder (default: %(default)s)")
    # parser.add_argument("--igblast-output", "-I",
    #                     type=str,
    #                     default="out.blast",
    #                     help="IgBLAST output file (default: %(default)s)")
    # parser.add_argument("--vjfinder-output", "-V",
    #                     type=str,
    #                     default="out",
    #                     help="VJFinder output dir (default: %(default)s)")
    parser.add_argument("--rerun-igblast", "-G",
                        action="store_true",
                        help="perform IgBLAST")
    parser.add_argument("--do-not-run-vjfinder", "-F",
                        action="store_true",
                        help="perform VJFinder")
    parser.add_argument("--tmp-file", "-T",
                        type=str,
                        default="",
                        help="temporary FASTA file used as IgBLAST input, <empty> for temporary file usage (default: <empty>)")
    parser.add_argument("--workdir", "-w",
                        type=str,
                        default="",
                        help="working directory, <empty> for home directory of this script (default: <empty>)")
    parser.add_argument("--storage-dir", "-S",
                        type=str,
                        default=".",
                        help="storage directory for cached IgBlast output (default: %(default)s)")

    args = parser.parse_args()

    if not args.tmp_file:
        args.tmp_file = tempfile.mkstemp(suffix=".fa", prefix="vjf_benchmarking_")[1]

    if not args.workdir:
        args.workdir = os.path.dirname(os.path.realpath(__file__))

    args.path = args.workdir + "/../../../"

    germline_J_file = args.workdir +"/germline/human/IGHJ-allP.fa"
    with open(germline_J_file, "rU") as fh:
        germline_J_parser = SeqIO.parse(fh, "fasta")
        germline_J_map = { str(record.id): str(record.seq) for record in germline_J_parser }

    args.input_hash = hash_file(args.input)
    args.igblast_output = args.storage_dir + "/" + args.input_hash + ".blast"

    if args.rerun_igblast or not os.path.exists(args.igblast_output + ".gz"):
        print "IgBLAST output will be written to " + args.igblast_output + ".gz"
        fq2fa(args.input, args.tmp_file)
        igblast_time = time.time()
        rcode = os.system("bash %(workdir)s/blast.sh %(tmp_file)s %(igblast_output)s 2> /dev/null" % args.__dict__)
        if rcode != 0:
            print "IgBLAST failed with the code", rcode
            exit(rcode)
        igblast_time = time.time() - igblast_time
        os.unlink(args.tmp_file)
        rcode = os.system("gzip %s" % args.igblast_output)
        if rcode != 0:
            print "gzip failed with the code", rcode
            exit(rcode)

        print "IgBLAST time:", igblast_time
    else:
        igblast_time = 0.

    args.vjfinder_output = args.storage_dir + "/" + args.input_hash + "_vjf"
    if not args.do_not_run_vjfinder or not os.path.exists(args.vjfinder_output):
        vjf_time = time.time()
        rcode = os.system("%(path)s/build/release/bin/ig_kplus_vj_finder --db-directory=%(workdir)s/germline -i %(input)s -o %(vjfinder_output)s --separator=tab -Z1 --loci=IGH --organism=human %(options)s" % args.__dict__)
        if rcode != 0:
            print "VJFinder failed with the code", rcode
            exit(rcode)
        vjf_time = time.time() - vjf_time
    else:
        vjf_time = 0.

    vjf_hits = parse_vjf_output("%s/alignment_info.csv" % args.vjfinder_output,
                                args.input)

    log = FakeLog()
    blast = igblast_utils.ParseIgBlastOutput(args.igblast_output + ".gz", log, smart_open)
    # Normilize blast_blocks
    igblast_hits = [line.hit_table for line in blast.blocks]

    # Read all seqs
    with smart_open(args.input, "rU") as fh:
        parser = SeqIO.parse(fh, idFormatByFileName(args.input))
        reads = list(parser)

    assert len(reads) == len(igblast_hits)
    ids = [str(record.description) for record in reads]
    assert len(set(ids)) == len(reads)


    contaminations = 0
    is_cont = []
    evalues = []

    RESULT = []
    for _i, line in enumerate(igblast_hits):
        genes = {}
        for row in line:
            genes[row.type] = row
        line.genes = genes


        RES = Empty()

        # if "V" not in genes or "J" not in genes or genes["V"].evalue > 0.001 or genes["J"].evalue > 100005000:
        if ("V" not in genes) or (genes["V"].evalue > 0.001):
            RES.contamination = True
            RES.evalue = genes["V"].evalue if "V" in genes else 99999
        else:
            RES.contamination = False
            RES.evalue = genes["V"].evalue

        RES.vjf_identified = (_i in vjf_hits)
        if RES.vjf_identified:
            vhit = vjf_hits[_i]["V"]
            jhit = vjf_hits[_i]["J"]
            RES.vjf_start = vhit.start
            RES.vjf_end = jhit.end

        if "V" not in genes:
            # print "Not Ig read"
            RES.igb_start = -9999999999
        if "J" not in genes:
            # print "Not Ig read"
            RES.igb_end = -9999999999
        if "V" in genes:
            RES.igb_start = genes["V"].q_start - genes["V"].s_start + 1
        if "J" in genes:
            RES.igb_end = genes["J"].q_end + len(germline_J_map[genes["J"].subject_id]) - genes["J"].s_end

        RESULT.append(RES)

    contminations = 0
    idnf_contaminations = 0
    missed_cont = 0
    cropped = 0
    vok = 0
    jok = 0
    vjok = 0
    all = len(RESULT)
    aligned = 0
    missed_conts = []

    bad_reads = []


    for _i, RES in enumerate(RESULT):
        if RES.contamination:
            contaminations += 1
        if RES.contamination and not RES.vjf_identified:
            idnf_contaminations += 1
        if not RES.contamination and not RES.vjf_identified:
            cropped += 1
        if RES.contamination and RES.vjf_identified:
            missed_conts.append(_i)
            missed_cont += 1

        if RES.vjf_identified and not RES.contamination:
            aligned += 1
            if RES.vjf_start == RES.igb_start:
                vok += 1
            if RES.vjf_end == RES.igb_end:
                jok += 1
            if RES.vjf_start == RES.igb_start and RES.vjf_end == RES.igb_end:
                vjok += 1
            if RES.vjf_start != RES.igb_start or RES.vjf_end != RES.igb_end:
                bad_reads.append(_i)

    if args.bad_reads:
        with smart_open(args.bad_reads, "w") as f:
            SeqIO.write([reads[_] for _ in bad_reads], f, idFormatByFileName(args.bad_reads))

    print "Overall reads %d" % all
    print "Found contaminations %d from %d" % (idnf_contaminations, contaminations)
    print "Missed contaminations %d " % missed_cont
    print "Discarded (ill-cropped) reads %d " % cropped
    print "Aligned %d" % aligned
    print "V OK, J OK, VJ OK %d %d %d" % (vok, jok, vjok)
    print "OK rate %f" % (float(vjok) / float(aligned))
    print "Ill rate %f" % (1. - float(vjok) / float(aligned))
    print "Missed contaminations", missed_conts
