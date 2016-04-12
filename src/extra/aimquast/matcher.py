#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys


igrc_path = "/ssd/ig_repertoire_constructor/"

trie_comp = "%s/build/release/bin/ig_trie_compressor " % igrc_path
ig_matcher_bin = "%s/build/release/bin/ig_matcher " % igrc_path

def parse_multiplicity(s):
    import re

    m = re.match(r".*_multiplicity_(\d+)_", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

def parse_size(s):
    import re

    m = re.match(r".*___size___(\d+)", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

def mult2mult(clustering_fa, reference_fa, match_filename):
    with open(clustering_fa) as f:
        clustering_mults = [parse_size(rec.id) for rec in SeqIO.parse(f, "fasta")]

    with open(reference_fa) as f:
       reference_mults = [parse_multiplicity(rec.id) for rec in SeqIO.parse(f, "fasta")]

    unmatched_references = set(range(len(reference_mults)))
    result = []

    with open(match_filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            cl_mult = clustering_mults[i]
            ref_mult = 0
            if line:
                dist = -int(line.split()[0])
                if dist == 0:
                    neibs = map(int, line.split()[1:])
                    if len(neibs):
                        ref_mult = sum(reference_mults[j] for j in neibs)
                        unmatched_references.difference_update(neibs)

            result.append((cl_mult, ref_mult))

    for i in unmatched_references:
        result.append((0, reference_mults[i]))

    return result

def match_fas(clustering_fa, reference_fa, ig_matcher_bin=ig_matcher_bin):
    from os import system

    # TODO Use proper tmps
    system("%s -i %s -I %s -o m1.match -O m2.match -k 10 --tau 1" % (ig_matcher_bin,
                                                                     reference_fa,
                                                                     clustering_fa))

    m2m = mult2mult(clustering_fa, reference_fa, "m2.match")

    with open(clustering_fa) as f:
        clustering_mults = [parse_size(rec.id) for rec in SeqIO.parse(f, "fasta")]

    with open(reference_fa) as f:
       reference_mults = [parse_multiplicity(rec.id) for rec in SeqIO.parse(f, "fasta")]

    reference_mults_sorted = sorted(reference_mults)
    def ref_count(limit):
        return how_many_greater_or_equal(limit, reference_mults_sorted)

    clustering_mults_sorted = sorted(clustering_mults)
    def cluster_count(limit):
        return how_many_greater_or_equal(limit, clustering_mults_sorted)

    clustering_paired_mults, reference_paired_mults = zip(*m2m)
    minimum_paired_mults = [min(clust_mult, ref_mult) for clust_mult, ref_mult in m2m]

    minimum_paired_mults_sorted = sorted(minimum_paired_mults)
    def sensitivity(limit):
        assert limit > 0
        d = how_many_greater_or_equal(limit, minimum_paired_mults_sorted)
        dd = ref_count(limit)
        if dd == 0:
            return 0.
        else:
            return float(d) / float(dd)

    def specificity(limit):
        assert limit > 0
        d = how_many_greater_or_equal(limit, minimum_paired_mults_sorted)
        dd = cluster_count(limit)
        if dd == 0:
            return 0.
        else:
            return float(d) / float(dd)

    def fdr(limit):
        return 1.0 - specificity(limit)


    def median_rate(limit):
        # TODO Speedup it
        import numpy as np
        assert limit > 0
        rates = [float(clust_mult) / float(ref_mult) for clust_mult, ref_mult in m2m if min(ref_mult, clust_mult) >= limit]

        return np.median(rates)

    class Res:
        pass

    res = Res()
    res.sensitivity, res.specificity, res.fdr = sensitivity, specificity, fdr

    res.median_rate = median_rate
    res.m2m = m2m

    return res

def how_many_greater_or_equal(limit, x):
    import bisect
    return len(x) - bisect.bisect_left(x, limit)

assert how_many_greater_or_equal(4, [4, 4, 4, 4]) == 4
assert how_many_greater_or_equal(4, [4]) == 1
assert how_many_greater_or_equal(4, []) == 0
assert how_many_greater_or_equal(4, [1, 2, 3, 4, 5]) == 2


res = match_fas("igrc_res_lam2.fa", "reference.fa")

from rand import make_ideal_rcm, rcm_vs_rcm

make_ideal_rcm("final_repertoire.rcm", "ideal.rcm")

print rcm_vs_rcm("final_repertoire.rcm", "ideal.rcm")

def error_prof(dirname, **kwargs):
    from rand import error_profile
    return error_profile(dirname + "/final_repertoire.rcm",
                         dirname + "/cleaned_reads.fa",
                         dirname + "/final_repertoire.fa",
                         **kwargs)



ep = error_prof("out")


def plot_error_profile(ep):
    import seaborn as sns
    sns.distplot(ep.errors01)
    sns.plt.show()


plot_error_profile(ep)
