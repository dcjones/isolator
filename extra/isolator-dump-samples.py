#!/usr/bin/env python2

# Convert the sqlite3 files output by isolator to a simpler plain text format.

import numpy as np
import argparse
import sqlite3
import struct
import zlib

def unpack_samples_blob(blob):
    """ Unpack compressed samples into an numpy array. """
    data = zlib.decompress(blob)
    fmt = '{0}f'.format(len(data) / 4)
    return np.array(struct.unpack(fmt, data), dtype=np.float32)


def main():
    ap = argparse.ArgumentParser(
        description='Convert output from isolator to plaintext, tab-delimited data.')
    ap.add_argument('db_fn', metavar='samples.db')
    ap.add_argument('--gene_id', default=None)
    args = ap.parse_args()

    db = sqlite3.connect(args.db_fn)
    cur = db.cursor()

    if args.gene_id:
        q = '''
            select gene_id, ids.transcript_id, effective_length, samples
            from samples, ids
            where samples.transcript_id = ids.transcript_id and
                  ids.gene_id = "{0}"
            '''.format(args.gene_id)
    else:
        q = '''
            select gene_id, ids.transcript_id, effective_length, samples
            from samples, ids
            where samples.transcript_id = ids.transcript_id
            '''

    for gene_id, transcript_id, effective_length, samples in cur.execute(q):
        samples = unpack_samples_blob(samples)
        samples /= effective_length
        print('{gene_id}\t{transcript_id}\t{samples}'.format(
            gene_id=gene_id,
            transcript_id=transcript_id,
            samples=','.join(map(str, samples))))


if __name__ == '__main__':
    main()


