#!/usr/bin/env python

import argparse
import functools
import multiprocessing
#import numpypy as np
import numpy as np
import numpy.random as npr
import sqlite3
import struct
import sys
import zlib
from itertools import izip
from collections import defaultdict

# number of js-divergence samples to generate
num_samples = 500


def unpack_samples_blob(blob):
    """ Unpack compressed samples into an numpy array. """
    data = zlib.decompress(blob)
    fmt = '{0}f'.format(len(data) / 4)
    return np.array(struct.unpack(fmt, data), dtype=np.float32)


def read_samples(cur, db_name):
    """ Return a dict mapping transcript_ids to numpy arrays of samples. """
    transcript_samples = {}
    q = 'select transcript_id, effective_length, samples from {0}.samples'.format(db_name)
    rows = cur.execute(q)
    for transcript_id, effective_length, samples in rows:
        samples = unpack_samples_blob(samples)
        samples /= effective_length
        transcript_samples[transcript_id] = samples
    return transcript_samples


def kldiv(xs, ys):
    """ Kullback-Leibler divergence. """
    return sum(x * np.log(x/y) for (x,y) in izip(xs, ys) if x > 0.0)


def jsdiv(xs, ys):
    """ Jenson-Shannon divergence between two multinomial distributions. """
    return (kldiv(xs, ys) + kldiv(ys, xs)) / 2.0


def kludge_test(data):
    (gene_id, sample_a, sample_b) = data

    n = len(sample_a)
    assert(len(sample_b) == n)

    # number of samples
    m_a = len(sample_a[0])
    m_b = len(sample_b[0])

    aexpr_samples = np.empty(num_samples, dtype=np.float32)
    bexpr_samples = np.empty(num_samples, dtype=np.float32)
    div_samples = np.empty(num_samples, dtype=np.float32)
    amix = np.empty(n, dtype=np.float32)
    bmix = np.empty(n, dtype=np.float32)

    eps = 1e-10

    for i in xrange(num_samples):
        k_a = npr.randint(0, m_a)
        k_b = npr.randint(0, m_b)
        for j in xrange(n):
            amix[j] = eps + sample_a[j][k_a]
            bmix[j] = eps + sample_b[j][k_b]
        aexpr_samples[i] = sum(amix)
        bexpr_samples[i] = sum(bmix)
        amix /= aexpr_samples[i]
        bmix /= bexpr_samples[i]
        div_samples[i] = jsdiv(amix, bmix)

    return (gene_id, np.median(div_samples),
            np.median(aexpr_samples),
            np.median(bexpr_samples),
            np.median(np.log2(aexpr_samples / bexpr_samples)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('afn', metavar='sample_a.db')
    ap.add_argument('bfn', metavar='sample_b.db')
    args = ap.parse_args()

    db = sqlite3.connect(':memory:')
    cur = db.cursor()

    cur.execute('attach "{0}" as sample_a'.format(args.afn))
    cur.execute('attach "{0}" as sample_b'.format(args.bfn))

    sample_a = read_samples(cur, 'sample_a')
    sample_b = read_samples(cur, 'sample_b')

    # gene_id -> trascript_id set mappings. This assumes they are the same in
    # sample_a and sample_b.
    ids = defaultdict(set)
    q = 'select gene_id, transcript_id from sample_a.ids'
    for gene_id, transcript_id in cur.execute(q):
        ids[gene_id].add(transcript_id)

    sys.stderr.write('Testing...\n')
    data = ((gene_id, [sample_a[tid] for tid in sorted(tids)],
                      [sample_b[tid] for tid in sorted(tids)])
            for (gene_id, tids) in ids.iteritems())
    pool = multiprocessing.Pool()
    print('gene_id,median_div,median_aexpr,median_bexpr,median_log2fc')
    for result in pool.imap_unordered(kludge_test, data):
        print('{0},{1},{2},{3},{4}'.format(*result))


if __name__ == '__main__': main()


