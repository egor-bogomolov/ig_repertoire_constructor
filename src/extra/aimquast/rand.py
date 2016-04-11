


def comb2(n):
    return n * (n - 1) // 2


def randstats(X, Y):
    from collections import defaultdict

    m = defaultdict(int)
    mX = defaultdict(int)
    mY = defaultdict(int)

    for x, y in zip(X, Y):
        m[(x, y)] += 1
        mX[x] += 1
        mY[y] += 1

    S00 = S01 = S10 = S11 = 0

    assert len(Y) == len(X)
    N = len(X)

    for x, y in zip(X, Y):
        _m = m[(x, y)]
        _mX = mX[x]
        _mY = mY[y]
        S00 += _m - 1
        S01 += _mX - _m
        S10 += _mY - _m
        S11 += N - _mX - _mY + _m


    S00 //= 2
    S01 //= 2
    S10 //= 2
    S11 //= 2

    assert S00 + S11 + S01 + S10 == N * (N - 1) // 2

    return S00, S11, S01, S10


def rand_adj(X, Y):
    from collections import defaultdict

    m = defaultdict(int)
    mX = defaultdict(int)
    mY = defaultdict(int)

    for x, y in zip(X, Y):
        m[(x, y)] += 1
        mX[x] += 1
        mY[y] += 1

    S00 = S01 = S10 = S11 = 0

    assert len(Y) == len(X)
    N = len(X)

    for x, y in zip(X, Y):
        _m = m[(x, y)]
        _mX = mX[x]
        _mY = mY[y]
        S00 += _m - 1
        S01 += _mX - _m
        S10 += _mY - _m
        S11 += N - _mX - _mY + _m


    # Compute the ARI using the contingency data
    sum_comb_c = sum(comb2(_) for _ in mY.itervalues())
    sum_comb_k = sum(comb2(_) for _ in mX.itervalues())
    sum_comb = sum(comb2(_) for _ in m.itervalues())

    prod_comb = (sum_comb_c * sum_comb_k) / float(comb2(N))
    mean_comb = (sum_comb_k + sum_comb_c) / 2.

    return (sum_comb - prod_comb) / (mean_comb - prod_comb)


def rand(X, Y):
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00 + S11) / float(S00 + S11 + S10 + S01)

def FM_index(X, Y):
    import math
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00) / math.sqrt((S00 + S10) * (S00 + S01))

def jaccard(X, Y):
    import math
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00) / (S00 + S10 + S01)

assert rand([0, 1], [0, 1]) == 1

import sklearn.metrics



assert sklearn.metrics.adjusted_rand_score([0, 1, 1], [0, 1, 2]) ==  rand_adj([0, 1, 1], [0, 1, 2])
assert sklearn.metrics.adjusted_rand_score([0, 1, 1], [0, 1, 1]) == rand_adj([0, 1, 1], [0, 1, 1])
assert sklearn.metrics.adjusted_rand_score([0, 0, 1, 2], [0, 0, 1, 1]) == rand_adj([0, 0, 1, 2], [0, 0, 1, 1])



def rcm2Rand(filename="./igrc_out/final_repertoire.rcm"):
    import sklearn
    import sklearn.metrics
    import re

    reference = []
    clustering = []

    with open(filename) as fh:
        for line in fh:
            id, cluster = line.strip().split("\t")
            cluster = int(cluster)
            m = re.match("^antibody_(\\d+)_", id)
            ant = int(m.groups()[0])
            clustering.append(ant)
            reference.append(cluster)

    # find large >= clusteres
    from collections import defaultdict

    clust2mult = defaultdict(int)
    for c in clustering:
        clust2mult[c] += 1

    reference_large = []
    clustering_large = []
    for ref, cluster in zip(reference, clustering):
        if clust2mult[cluster] >= 5:
            reference_large.append(ref)
            clustering_large.append(cluster)

    # print reference
    # print clustering
    assert rand(clustering, clustering) == FM_index(clustering, clustering) == jaccard(clustering, clustering) == 1

    return FM_index(clustering, reference), rand(clustering, reference), rand_adj(clustering, reference), jaccard(clustering, reference), FM_index(clustering_large, reference_large), rand(clustering_large, reference_large), rand_adj(clustering_large, reference_large),jaccard(clustering_large, reference_large)


if __name__ == "__main__":
    print rcm2Rand()
