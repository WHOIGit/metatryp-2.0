"""
name: digest_and_ingest.py

usage: digest_and_ingest.py [--digest-config=digest_config_file] fasta1 fasta2 ...

commissioned by : Dr. Makoto Saito, 2013-03

authorship: adorsk, 2013-03

description: This script digests peptides from protein sequences stored in
FASTA files. 
Assumptions:
    - Each FASTA file is assumed to contain the proteome for one taxon.
    - FASTA file names contain the taxon id e.g. 'syn5802.fasta'.
    - The redundancy db has already been created and is writeable.
"""

"""
Imports and setup.
"""
from proteomics import db
from proteomics import config
from proteomics.models import (Protease, Digest)
from proteomics.services.digest_and_ingest_specialized_assembly import DigestAndIngestSpecializedAssemblyTask
import argparse
import logging
import json
import time
import os

"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Digest and ingest peptides from FASTA protein sequence files.'))
argparser.add_argument('--digest-def', help=(
    'JSON file containing a digest definition. If not provided, default digest'
    'will be used'))
argparser.add_argument('--digest_dir', help=(
    'Directory containing all the FASTA files to be processed.  Can be used for batch processing'))
argparser.add_argument('--annotation-file', help=(
    'File containing anotations.  Should be used in conjunction with annotation-source'  ))
argparser.add_argument('--annotation-source', help=(
    'Source of annotations. Valid values: "venter"'  ))
argparser.add_argument('--fasta_files', nargs='+',
                    help=('List of FASTA files containing protein sequences.'))
"""
Main method.
"""
def main():
    start_time = time.time()
    args = argparser.parse_args()

    logger = logging.getLogger('digest_and_ingest')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    fasta_files = []
    if args.fasta_files:
        fasta_files = args.fasta_files
    if args.digest_dir:
        try:
            directory = args.digest_dir

            for filename in os.listdir(directory):
                if filename.endswith(".faa") or filename.endswith(".fasta"):
                    fullpath=os.path.join(directory, filename)
                    logger.info("FASTA file: %s" % fullpath)
                    fasta_files.append(fullpath)
        except:
            logger.exception("Could not parse digest_dir file '%s'" % (
                args.digest_dir))
    # Parse digest definition if given.
    if args.digest_def:
        try:
            digest_def = json.loads(args.digest_def)
        except:
            logger.exception("Could not parse digest_def file '%s'" % (
                args.digest_def))
    # Otherwise use default digest.
    else:
        digest_def = config.DEFAULT_DIGEST_DEFINITION

    # if args.annotation-file:
    #     try:
    #         directory = args.digest_dir
    #
    #         for filename in os.listdir(directory):
    #             if filename.endswith(".faa") or filename.endswith(".fasta"):
    #                 fullpath=os.path.join(directory, filename)
    #                 logger.info("FASTA file: %s" % fullpath)
    #                 fasta_files.append(fullpath)
    #     except:
    #         logger.exception("Could not parse digest_dir file '%s'" % (
    #             args.digest_dir))
    # Get or create the digest.
    digest = get_digest(logger, digest_def)
    logger.info("Digest 2'%s'." % digest.id)
    # Run the digest/ingest task.
    task = DigestAndIngestSpecializedAssemblyTask(
        logger=logger,
        fasta_paths=fasta_files,
        digest=digest,
        #get_connection=db.get_connection,
    )
    stats = task.run()
    logger.info("Statistics on records created: %s" % stats)
    total_time = time.time() - start_time
    logger.info("total run time: %s" % (total_time))
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
        protease = Protease(**digest_def['protease'])
        cur.execute("insert into protease (id, cleavage_rule) values( %s, %s);", (str(digest_def['protease']['id']), str(digest_def['protease']['cleavage_rule']),))
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
            "No digest exists for the given definition, creating...")
        digest_kwargs = {}
        digest_kwargs.update(digest_def)
        digest_kwargs['protease'] = protease
        digest = Digest(**digest_kwargs)
        cur = db.get_psycopg2_cursor()


        cur.execute("select * from digest_insert( %s, %s, %s, %s);", (protease.id, digest.max_missed_cleavages, digest.min_acids, digest.max_acids,))
        digest_result = cur.fetchone()

        if digest_result:
            digest = Digest(id=digest_result[0], protease = protease, max_missed_cleavages=digest_result[2], min_acids = digest_result[3], max_acids = digest_result[4])
    else:
        digest = Digest(id=results[0], protease = protease, max_missed_cleavages=results[2], min_acids = results[3], max_acids = results[4])
    db.psycopg2_connection.commit()
    return digest

if __name__ == '__main__':
    main()
