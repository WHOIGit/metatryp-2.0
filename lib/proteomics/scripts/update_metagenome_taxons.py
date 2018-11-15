from proteomics import db
from proteomics.models import Taxon

import logging
import csv
import time

def main():
    start_time = time.time()

    logger = logging.getLogger('update_metagenome_taxons')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    filename = 'C:/Users/David/Projects/metatryp-2.0/data_files/ProteOMZ/new_orf_taxa_keys.csv'
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        try:
            count=1;
            for row in reader:

                 #print count
                # count = count+1
                 s = "."
                 seq = (row[1], row[2], row[3], row[4],row[5], row[6], row[7], row[8])
                 hierachy = s.join(seq)
                 hierachy = hierachy.replace(" ", "_")
                 hierachy = hierachy.replace("-", "_")
                 hierachy = hierachy.replace("(", "_")
                 hierachy = hierachy.replace(")", "_")
                 hierachy = hierachy.replace("[", "")
                 hierachy = hierachy.replace("]", "")
                 hierachy = hierachy.replace("'", "")
                 hierachy = hierachy.replace("'", "")
                 #print hierachy

                 cur = db.get_psycopg2_cursor()
                 cur.execute("select mt.ncbi_id from metagenome_taxon mt where mt.ncbi_id =%s;", (row[0],))
                 taxon_result = cur.fetchone()

                 if taxon_result is None:
                     ncbi_url = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+row[0]
                     print(ncbi_url)
                     cur.execute("insert into metagenome_taxon (ncbi_id, tax_group, tax_kingdom, tax_phylum, tax_class, "
                                 "tax_order, tax_family, tax_genus, tax_species, ncbi_url, hierachy)"
                                 " values(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",  (int(row[0]), row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],ncbi_url, hierachy))
                     db.psycopg2_connection.commit()

        except csv.Error as e:
           logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))


if __name__ == '__main__':
    main()
