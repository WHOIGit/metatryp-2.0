from proteomics import config
from proteomics import db
import logging


def main():
    logger = logging.getLogger('taxons')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    cur = db.get_psycopg2_cursor();
    cur.execute("select t.id from taxon t;")

    logger.info("Taxons");
    for record in cur:
        logger.info("%s" % (record[0]))
    db.psycopg2_connection.commit()

if __name__ == '__main__':
    main()

