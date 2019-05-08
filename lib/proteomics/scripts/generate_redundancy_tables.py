"""
name: generate_redundancy_tables.py

commissioned by : Dr. Makoto Saito, 2013-03

authorship: adorsk, 2013-03

description: This script generates redundancy tables for a given set of taxon
digests.

Assumptions:
    - The redundancy db has already been created and is readable.
"""

"""
Imports and setup.
"""
from proteomics import config
from proteomics import db
from proteomics.models import (Taxon, Protease, Digest, TaxonDigest)
from proteomics.services import redundancy
import argparse
import logging
import csv
import os


"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Generate redundancy tables for the given taxon digests.'))
argparser.add_argument('--taxon-id-file', help=(
    'A file containing a list of taxon IDs to include in the redundancy table,'
    ' one taxon id per line. The default digest will be used to search for'
    'TaxonDigests.'))
argparser.add_argument('--taxon-ids', nargs='*', help=(
    'List of taxon IDs to include in the redundancy table. The default digest'
    ' will be used to search for TaxonDigests. This option will override the'
    ' --taxons-file option.'))
argparser.add_argument('--output-dir', required=True, help=(
    'Output directory. CSV tables will be written to this directory.'))

"""
Main method.
"""
def main():
    args = argparser.parse_args()

    logger = logging.getLogger('redundancy_tables')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    # Check that taxon ids or taxon id file were provided.
    if not (args.taxon_ids or args.taxon_id_file):
        raise Exception("Must provide --taxon-ids or --taxon-id-file option")

    # Get taxons.
    if args.taxon_ids:
        taxon_ids = args.taxon_ids
    else:
        with open(args.taxon_id_file, 'rb') as f:
            taxon_ids = [row[0] for row in csv.reader(f)]
    logger.info("Taxon Ids: %s" % (taxon_ids))
    # Get the digest.
    digest = get_digest(logger, config.DEFAULT_DIGEST_DEFINITION)
    cur = db.get_psycopg2_cursor();
    cur.execute( "select * from taxon_digest td where td.taxon_id in (select t.id from taxon t where t.id = any(%s));", (taxon_ids,))

    taxon_digests = []

    for record in cur:
        td = TaxonDigest(id=record[0], taxon=record[1], digest=record[2])
        logger.info("Taxon Digest: %s" % (td.id) )

        taxon_digests.append(td)
    db.psycopg2_connection.commit()



    # Generate the redundancy tables.
    tables = redundancy.generate_redundancy_tables(taxon_digests, logger=logger)

    # Create output dir if it does not exist.
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Output tables.
    for table_id, table in list(tables.items()):
        table_file = os.path.join(args.output_dir, table_id + '.csv')
        logger.info("Writing '%s'..." % table_file)
        with open(table_file, 'w',  newline='') as f:
            w = csv.writer(f)
            for row in table:
                w.writerow(row)

    logger.info("Done.")

""" Helper methods. """
def get_digest(logger, digest_def):
    """ Fetch or create a digest from a digest definition."""
    #session = db.get_session()

    # Get or create protease.
    protease = Protease(**digest_def['protease'])
    protease_id = str(digest_def['protease']['id'])
    cur = db.get_psycopg2_cursor()
    cur.execute("select * from protease where protease.id=%s;", (protease_id,))
    results = cur.fetchone()

    if results is None:
        logger.info(
            "No protease exists for the given definition, creating...")
    else:
        protease = Protease(id=results[0], cleavage_rule=results[1])
    db.psycopg2_connection.commit()

    # Get or create digest object.
    cur = db.get_psycopg2_cursor()

    #not all possible digestion parameters will have a value so build the query to account for this
    query_params = [protease.id]
    digest_query = "select * from digest where digest.protease_id = %s";
    if digest_def.get('max_missed_cleavages') is not None:
        digest_query = digest_query + " and digest.max_missed_cleavages = %s "
        query_params.append(digest_def.get('max_missed_cleavages'))
    if digest_def.get('min_acids') is not None:
        digest_query = digest_query + " and digest.min_acids = %s "
        query_params.append(digest_def.get('min_acids'))
    if digest_def.get('max_acids') is not None:
        digest_query = digest_query + " and digest.max_acids = %s "
        query_params.append(digest_def.get('max_acids'))

    cur.execute(digest_query, (query_params))
    results = cur.fetchone()
    db.psycopg2_connection.commit
    if results is None:
    #if not digest:
        logger.info(
            "No digest exists for the given definition.")
    else:
        digest = Digest(id=results[0], protease = protease, max_missed_cleavages=results[2], min_acids = results[3], max_acids = results[4])
    db.psycopg2_connection.commit()
    return digest

if __name__ == '__main__':
    main()
