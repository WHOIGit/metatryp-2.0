"""
name: query_by_sequence.py

usage: query_by_sequence.py [--max-distance=0] sequence_file 

commissioned by : Dr. Makoto Saito, 2013-03

authorship: adorsk, 2013-05 modified for api by dgaylord 2016-12

description: This script queries a peptides database for the given set of
peptide sequences.

Outputs: a JSON file containing
    query_sequence | taxon_id | levenshtein_distance | match_sequence
"""

"""
Imports and setup.
"""
from proteomics import db
from proteomics.models import (Peptide, TaxonDigestPeptide, TaxonDigest)
import logging
from sqlalchemy.sql import func
from flask import Flask, request
from flask_restful import Resource, Api


class QueryBySequence(Resource):
    def get(self):
        args = request.args
        logger = logging.getLogger('query_by_sequence')
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)

        # Define levenshtein function in SQLite.
        try:
            def levenshtein(s1,s2):
                l1 = len(s1)
                l2 = len(s2)
                matrix = [range(l1 + 1)] * (l2 + 1)
                for zz in range(l2 + 1):
                  matrix[zz] = range(zz,zz + l1 + 1)
                for zz in range(0,l2):
                  for sz in range(0,l1):
                    if s1[sz] == s2[zz]:
                      matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
                    else:
                      matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
                return matrix[l2][l1]

            connection = db.get_connection()
            #connection.connection.create_function("LEVENSHTEIN", 2, levenshtein)
        except Exception as e:
            logger.exception('Could not define Levenshtein distance function: %s' % e)
            raise e

        session = db.get_session(bind=connection)

        # Read in sequences to query.
        results = ''
        sequences = {}
        max_distance = 1;
        if args:

            if args.get('sequence_file') is not None:

                with open(args.get('sequence_file'), 'rb') as f:
                    sequences = [line.strip() for line in f.readlines()]
            elif args.get('sequence') is not None:
                sequences = {args.get('sequence')}

            if not sequences:
                results = "Provide a query sequence via the '--sequence' option, or a set of sequences via the --sequence-file option"
                return results

            #set a default max levenshtein distance
            if args.get('max_distance'):
                max_distance =  args.get('max_distance')
        # Print headers.
        headers = ['query', 'taxon', 'lev_distance', 'match']

        # Execute query for each sequence and print results.

        results = {
            "results": []
        }

        for seq in sequences:
            data = {"search_sequence":seq}

            lev_dist = func.levenshtein(Peptide.sequence, str(seq))
            q = (session.query(TaxonDigest.taxon_id, lev_dist,
                               Peptide.sequence)
                 .select_from(Peptide)
                 .join(TaxonDigestPeptide)
                 .join(TaxonDigest)
                 .filter(lev_dist <= max_distance)
                )
            for row in q:
                data["taxon"] = row[0]
                data["lev_distance"] = row[1]
                data["match"] = row[2]

        results["results"].append(data)
        return results


