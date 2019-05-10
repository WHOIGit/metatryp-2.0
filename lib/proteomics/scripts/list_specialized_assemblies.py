from proteomics import config
from proteomics import db
import logging


def main():
    logger = logging.getLogger('specialized_assemblies')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    cur = db.get_psycopg2_cursor();
    cur.execute("select sa.genome_name, sa.type_flag from specialized_assembly sa;")

    logger.info("Specialized Assemblies")
    logger.info("Genome Name            Type")
    for record in cur:
        logger.info("%s     %s" % (record[0], record[1]))
    db.psycopg2_connection.commit()

if __name__ == '__main__':
    main()

