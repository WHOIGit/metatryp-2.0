from proteomics import db
from proteomics.models import Taxon

import logging
import csv
import time

def main():
    start_time = time.time()

    logger = logging.getLogger('update_taxons')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    filename = 'C:/Users/David/Projects/metatryp-2.0/data_files/NAH_microbialgenomes/171211_NAH_updated_taxontable_withlabels_update.csv'
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        try:
            for row in reader:

                cur = db.get_psycopg2_cursor()
                #check to see if taxon file matching the fasta file key is in database

                cur.execute("select * from taxon where id = %s;", (row[16],))
                taxon_result = cur.fetchone()

                db.psycopg2_connection.commit()
                if taxon_result is not None:
                    print row
                    s = "."
                    seq = (row[6], row[7], row[8], row[9],row[10], row[11], row[12])
                    hierachy = s.join(seq)
                    hierachy = hierachy.replace(" ", "_")
                    print hierachy

                    cur = db.get_psycopg2_cursor()
                    cur.execute("update taxon set tax_oid = %s,  status = %s, study_name = %s,"
                                "genome_name= %s, sequencing_center= %s, img_genome_id= %s, tax_form = %s, tax_domain = %s, tax_phylum= %s, tax_class= %s, "
                                "tax_order= %s, tax_family= %s, tax_genus= %s, tax_species= %s, genome_size= %s, gene_count= %s, hierachy= %s, ncbi_id= %s, ncbi_url= %s, img_url= %s"
                                "where id= %s", (long(row[0]), row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],long(row[14]),long(row[15]), hierachy, long(row[17]), row[18], row[20], row[16]))
                    db.psycopg2_connection.commit()



        except csv.Error as e:
           logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))


if __name__ == '__main__':
    main()
