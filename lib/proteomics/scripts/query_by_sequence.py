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
from proteomics.services.query_by_sequence import QueryBySequenceTask
from proteomics.services.query_by_sequence_lca import QueryBySequenceLCATask
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

argparser.add_argument('--lca', help='Perform LCA analysis', action="store_true")

"""
Main method.
"""

def main():
    args = argparser.parse_args()
    logger = logging.getLogger('query_by_sequence')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    if args.lca:
        task = QueryBySequenceLCATask(
            logger=logger,
            args=args,
        
        )
        stats = task.run()
    else:
        task = QueryBySequenceTask(
            logger=logger,
            args=args,
        
        )
        stats = task.run()            

if __name__ == '__main__':

    main()

