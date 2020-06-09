from proteomics import config
from proteomics import db
import logging


def main():
    logger = logging.getLogger('metaomic_taxons')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    cur = db.get_psycopg2_cursor();
    cur.execute("select mt.tax_species from metagenome_taxon mt where mt.ncbi_id in \
        (select distinct ma.contig_tax_id from metagenome_annotations ma) or \
        mt.ncbi_id in (select distinct ma.orf_tax_id from metagenome_annotations ma)") 

    logger.info("Meta-omic Assembly Taxons")
    for record in cur:
        logger.info("%s" % (record[0]))
    db.psycopg2_connection.commit()

if __name__ == '__main__':
    main()

