#!/usr/bin/env python

import sqlite3, sys

if __name__ == '__main__':

    import argparse

    description = '''A tool to extract sequences from a phlawd db'''

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--phlawd-db', type=sqlite3.connect, \
        required=True, help='The phlawd database from which to extract sequences.')

    parser.add_argument('-l', '--gilist', type=file, \
        required=False, help='A file containing a list with the gi numbers to be extracted. If omitted, stdin will be used.')

    args = parser.parse_args()

    phlawd_db = args.phlawd_db
    c = phlawd_db.cursor()

    src = None
    if args.gilist is not None:
        src = args.gilist
    else:
        src = sys.stdin

    gis = []
    for i in src.readlines():
        try:
            gis.append(int(i))
        except ValueError:
            continue

    for i in gis:
        r = c.execute('SELECT * FROM seq WHERE ncbi_id LIKE ?', (i,))
        print(r.fetchone())