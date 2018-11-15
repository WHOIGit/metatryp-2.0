"""
name: query_by_sequence.py

usage: query_by_sequence.py [--max-distance=0] sequence_file 

commissioned by : Dr. Makoto Saito, 2013-03

authorship: adorsk, 2013-05, updated by David Gaylord 2016-12

description: This script queries a peptides database for the given set of
peptide sequences.

Outputs: a CSV document to stdout whose rows contains:
    query_sequence | taxon_id | levenshtein_distance | match_sequence
"""

"""
Imports and setup.
"""
from proteomics import db
import argparse
import logging

"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Query database for peptide sequences'))
argparser.add_argument('--max-distance', type=int, nargs='?', default=0,
                       help=(
                           'List of FASTA files containing protein sequences.'
                       ))
argparser.add_argument('--sequence-file', help=(
    'File containing one amino acid sequence per line.'))

argparser.add_argument('--sequence', help='Amino acid sequence')

argparser.add_argument('--type', help='Type of search to perform')

"""
Main method.
"""
def main():
    args = argparser.parse_args()

    logger = logging.getLogger('query_by_sequence')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    # Read in sequences to query.
    sequences = []
    max_dist = args.max_distance
    if args.sequence_file:
        with open(args.sequence_file, 'rb') as f:
            sequences = [line.strip() for line in f.readlines()]
    elif args.sequence:
        sequences = [args.sequence]

    # Read in whether to query just genomes ('g'), just metagenomes ('g') or both ('b').
    # Genomes is the default search
    type = 'g'
    if args.type:
        type = args.type

    if not sequences:
        argparser.error("Provide a query sequence via the '--sequence' option, "
                        "or a set of sequences via the --sequence-file option")

    # Print headers.
    headers = ['query', 'taxon', 'lev_distance', 'match']
    print(','.join(headers))


    # Execute query for each sequence and print results.
    cur = db.get_psycopg2_cursor()
    for seq in sequences:
        if type == 'g' or type == 'b':
            print('GENOMIC RESULTS FOR')
            cur.execute("select * from genomic_query_by_peptide_sequence(%s)", (seq,))
            for row in cur.fetchall():
                print(','.join([str(s) for s in [seq] + list(row)]))
        if type == 'm' or type == 'b':
            print('\n\n');
            print('METAGENOMNIC RESULTS')
            cur.execute("select * from metagenomic_query_by_peptide_sequence(%s)", (seq,))
            for row in cur.fetchall():
                print(','.join([str(s) for s in [seq] + list(row)]))
    db.psycopg2_connection.commit()
if __name__ == '__main__':
    main()
