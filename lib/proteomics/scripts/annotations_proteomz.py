from proteomics import db
import csv
import argparse
import logging
import time


argparser = argparse.ArgumentParser(description=('Import Annotation file'))
argparser.add_argument('--annotation_file', help=('Annotation file.'))

#filename = "data_files/ProteOMZ/2018_0515_ProteOMZ_for_ingest.csv"
def main():
    start_time = time.time()
    args = argparser.parse_args()

    logger = logging.getLogger('proteomz_annotations')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    filename = args.annotation_file

    # Parse digest definition if given.
    if filename:
        logger.info(
            "Starting annotation ingest")
        metagenome_annotations = {}
        metagenome_sequence_ids = []
        accession_numbers = []
        best_hit_annotations = []
        orf_tax_ids = []  #contig id?
        orf_taxonomies = []  #best hit species
        orf_tax_levels = []  # default to all being species

        #open the file and iterate through it line by line
        total_annotations = 0
        with open(filename, 'rb') as f:
            reader = csv.reader(f)
            try:
                batch_count = 0;  #handle ingestion in batches to reduce number of hits on the db
                for annot in reader:

                    accession_number = annot[0]
                    annotation = annot[1]
                    if annotation == '#N/A' or annotation == '':
                        annotation = None

                    taxon = annot[2]
                    if taxon == '#N/A' or taxon == '':
                        taxon = None

                    tax_id = annot[3]
                    if tax_id is not None:
                        if tax_id == '#N/A' or tax_id == '':
                            tax_id = None
                        else:
                            tax_id = int(tax_id)

                    tax_level = None
                    if tax_id is not None:
                        tax_level = "species"  # this is specific to the proteomz dataset as all annotations are categorized at the species level

                    #get all the annotations we need from the csv and save them into a dict for later use
                    ma = {"accession_number" : accession_number, "annotation" : annotation, "orf_tax_level" :  tax_level, "orf_tax_id" : tax_id, "orf_taxonomy" : taxon }
                    accession_numbers.append(accession_number)
                    metagenome_annotations[accession_number] = ma

                    batch_count = batch_count + 1

                    if batch_count == 1000:
                        metagenome_sequence_ids = []
                        a_numbers = []
                        best_hit_annotations = []
                        orf_tax_ids = []  # contig id?
                        orf_taxonomies = []  # best hit species
                        orf_tax_levels = []  # default to all being species
                        #add uniprot ids ???????????????/

                        #look up the already populated sequence id to match it to the accession number
                        cur = db.get_psycopg2_cursor()
                        cur.execute("select id, sequence_id from metagenome_sequence ms where ms.sequence_id = any(%s);",
                                             (accession_numbers,))

                        db.psycopg2_connection.commit()

                        for seq_id, a_number in cur:

                            if seq_id is not None:
                                batch_count = batch_count+1

                                #populate the arrays that will be passed to postgres for a bulk insert
                                metagenome_sequence_ids.append(seq_id)
                                meta_annon = metagenome_annotations[a_number];
                                #print meta_annon["accession_number"]

                                a_numbers.append(meta_annon["accession_number"])
                                best_hit_annotations.append(meta_annon["annotation"])
                                #print meta_annon["orf_tax_id"]
                                orf_tax_ids.append(meta_annon["orf_tax_id"])
                                orf_taxonomies.append(meta_annon["orf_taxonomy"])
                                orf_tax_levels.append(meta_annon["orf_tax_level"])

                        cur = db.get_psycopg2_cursor()
                        cur.execute("select * from metagenome_annotation_insert(%s, %s, %s, %s, %s, %s);",
                                    (a_numbers, metagenome_sequence_ids,best_hit_annotations,
                                     orf_tax_ids, orf_taxonomies, orf_tax_levels))

                        db.psycopg2_connection.commit()
                        logger.info("ingested: %s annotations" % (total_annotations))
                        # logger.info("annotations" % (the_annotations))
                        metagenome_annotations = {}
                        accession_numbers = []
                        total_annotations = total_annotations+batch_count
                        batch_count = 0


            except csv.Error as e:
                logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))
                print e

if __name__ == '__main__':
    main()