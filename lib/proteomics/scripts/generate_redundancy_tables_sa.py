"""
name: generate_redundancy_tables_sa.py

commissioned by : Dr. Makoto Saito


description: This script generates redundancy tables for a given set of specialized assemblies

Assumptions:
    - The redundancy db has already been created and is readable.
"""

"""
Imports and setup.
"""
from proteomics import config
from proteomics import db
from proteomics.models import (Specialized_Assembly)
from proteomics.services import redundancy
import argparse
import logging
import csv
import os


"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Generate redundancy tables for the given specialized assembly digests.'))
argparser.add_argument('--sa-id-file', help=(
    'A file containing a list of specialized assembly (i.e MAGs) IDs to include in the redundancy table,'
    ' one taxon id per line. The default digest will be used to search for'
    'Specialized Assembly Digests.'))
argparser.add_argument('--sa-ids', nargs='*', help=(
    'List of taxon IDs to include in the redundancy table. The default digest'
    ' will be used to search for Specialized Assembly Digests. This option will override the'
    ' --sa-file option.'))
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

    # Check that specialized assembly genome names or file were provided.
    if not (args.sa_ids or args.sa_id_file):
        raise Exception("Must provide --sa-ids or --sa-id-file option")

    # Get specialized assemblies.
    if args.sa_ids:
        sa_ids = args.sa_ids
    else:
        with open(args.sa_id_file, 'r') as f:
            sa_ids = [row[0] for row in csv.reader(f)]
    logger.info("Specialized Assembly Ids: %s" % (sa_ids))

    cur = db.get_psycopg2_cursor();

    cur.execute( "select sa.id, sa.genome_name from specialized_assembly sa where sa.genome_name = any(%s);", (sa_ids,))
    sa_digests = []

    for record in cur:
        sa = Specialized_Assembly(id=record[0],genome_name=record[1])

        logger.info("Specialized Assembly: %s  %s" % (sa.genome_name, sa.id) )

        sa_digests.append(sa)
    db.psycopg2_connection.commit()



    # Generate the redundancy tables.
    tables = redundancy.generate_redundancy_tables_sa(sa_digests, logger=logger)

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



if __name__ == '__main__':
    main()
