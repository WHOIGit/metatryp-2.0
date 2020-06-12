from proteomics import db

import logging
import csv
import time
import argparse
import os

def main():
    """
    Process arguments.
    """
    argparser = argparse.ArgumentParser(description=(
        'Update metagenome lineage from CSV file.'))
    argparser.add_argument('--filepath', help=(
        'CSV file containing metagenome lineage'))

    start_time = time.time()
    logger = logging.getLogger('update_mategenome_taxons')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    args = argparser.parse_args()

    if args.filepath:
        filename = args.filepath
        print(filename)

        with open(filename, 'rt') as f:
            reader = csv.reader(f)
            try:
                count=1;
                for row in reader:
                    if count > 1:
                        #fasta_file_name = row[0]
                        #genome_id = os.path.splitext(os.path.basename(fasta_file_name))[0]
                        taxon_name = row[0]
                        ncbi_id = row[1]
                        tax_group = row[2]
                        tax_kingdom = row[3]
                        tax_phylum = row[4]
                        tax_class = row[5]
                        tax_order = row[6]
                        tax_family = row[7]
                        tax_genus = row[8]
                        tax_species = row[9].replace(".", "")
                        ncbi_taxon_name = row[10]

                        print("Taxon Name: ", taxon_name)
                        print("ID: ", ncbi_id)
                        s = "."
                        seq = (tax_group, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family, tax_genus, tax_species)
                        hierachy = s.join(seq)
                        hierachy = hierachy.replace(" ", "_")
                        hierachy = hierachy.replace("-", "_")
                        hierachy = hierachy.replace("+", "_")
                        hierachy = hierachy.replace("(", "_")
                        hierachy = hierachy.replace(")", "_")
                        hierachy = hierachy.replace("/", "_")
                        hierachy = hierachy.replace("=", "")
                        hierachy = hierachy.replace("[", "")
                        hierachy = hierachy.replace("]", "")
                        hierachy = hierachy.replace("'", "")
                        hierachy = hierachy.replace("'", "")
                        #print hierachy

                        cur = db.get_psycopg2_cursor()

                        ncbi_url = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" + ncbi_id
                        print("URL: ", ncbi_url)

                        cur.execute("select mt.ncbi_id from metagenome_taxon mt where mt.ncbi_id =%s;", (ncbi_id,))
                        taxon_result = cur.fetchone()
                        if taxon_result is None:
                            cur.execute(
                                "insert into metagenome_taxon(ncbi_id, tax_group, tax_kingdom, tax_phylum, tax_class, "
                                "tax_order, tax_family, tax_genus, tax_species, hierachy, ncbi_url)"
                                " values(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);", (
                                ncbi_id, tax_group, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family,
                                tax_genus, tax_species, hierachy, ncbi_url))
                            db.psycopg2_connection.commit()
                        else:
                            cur.execute(
                                "update metagenome_taxon set ncbi_id = %s, tax_group = %s, tax_kingdom = %s, tax_phylum = %s, tax_class = %s, "
                                "tax_order = %s, tax_family = %s, tax_genus = %s, tax_species = %s, hierachy = %s, ncbi_url = %s where ncbi_id=%s",
                                (
                                    ncbi_id, tax_group, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family,
                                    tax_genus, tax_species, hierachy, ncbi_url, ncbi_id))
                            db.psycopg2_connection.commit()


                    count = count + 1
            except csv.Error as e:
                logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))


if __name__ == '__main__':
    main()

