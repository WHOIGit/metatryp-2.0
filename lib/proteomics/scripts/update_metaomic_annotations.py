from proteomics import db
import csv
import argparse
import logging
import time


argparser = argparse.ArgumentParser(description=('Import Annotation file'))
argparser.add_argument('--annotation_file', help=('Annotation file.'))

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
        orf_ids = [] #i.e. node_123_orf format
        contig_ncbi_tax_ids = []  
        contig_tax_names = []
        orf_ncbi_tax_ids = [] 
        orf_tax_names = [] 
       

        #open the file and iterate through it line by line
        total_annotations = 0
        line_count=0
        with open(filename, 'rt') as f:
            reader = csv.reader(f)
            line_count = len(list(reader))
            
        with open(filename, 'rt') as f:
            reader = csv.reader(f)
         
            try:
                count = 1;
     
                batch_count = 0;  #handle ingestion in batches to reduce number of hits on the db
                for annot in reader:
                    if count > 1:                   
                        orf_id = annot[0]
                        contig_tax_name = annot[1]
                        if not contig_tax_name:
                            contig_tax_name = None

                        contig_ncbi_tax_id = annot[2]
                        #print(contig_ncbi_tax_id)
                        if  not contig_ncbi_tax_id or contig_ncbi_tax_id == '#N/A' or  contig_ncbi_tax_id == '':
                            contig_ncbi_tax_id = None
                        else:
                            contig_ncbi_tax_id = int(contig_ncbi_tax_id)
                        orf_tax_name = annot[3]
                        if not orf_tax_name:
                            orf_tax_name = None

                        orf_ncbi_tax_id = annot[3]
                        if  not orf_ncbi_tax_id or orf_ncbi_tax_id == '#N/A' or  orf_ncbi_tax_id == '':
                            orf_ncbi_tax_id = None
                        else:
                            orf_ncbi_tax_id = int(orf_ncbi_tax_id)

                        #get all the annotations we need from the csv and save them into a dict for later use
                        ma = {"orf_id" : orf_id, "contig_tax_name" : contig_tax_name, "contig_tax_id" :  contig_ncbi_tax_id, "orf_tax_name" : orf_tax_name,  "orf_tax_id" : orf_ncbi_tax_id }
                        orf_ids.append(orf_id)
                        metagenome_annotations[orf_id] = ma

                        batch_count = batch_count + 1
                        
                        if batch_count == 1000 or count == line_count:
                            a_numbers = []
                            metagenome_sequence_ids = []
                            contig_ncbi_tax_ids = []  
                            contig_tax_names = []
                            orf_ncbi_tax_ids = [] 
                            orf_tax_names = [] 
                            
                            #look up the already populated sequence id to match it to the accession number
                            cur = db.get_psycopg2_cursor()
                            cur.execute("select id, sequence_id from metagenome_sequence ms where ms.sequence_id = any(%s);",
                                                (orf_ids,))

                            db.psycopg2_connection.commit()

                            for seq_id, a_number in cur:
                                
                                if seq_id is not None:
                                    batch_count = batch_count+1

                                    #populate the arrays that will be passed to postgres for a bulk insert
                                    metagenome_sequence_ids.append(seq_id)
                                    meta_annon = metagenome_annotations[a_number];
                                    #print meta_annon["accession_number"]

                                    a_numbers.append(meta_annon["orf_id"])
                                    contig_ncbi_tax_ids.append(meta_annon["contig_tax_id"])
                                    contig_tax_names.append(meta_annon["contig_tax_name"])
                                    orf_ncbi_tax_ids.append(meta_annon["orf_tax_id"])
                                    orf_tax_names.append(meta_annon["orf_tax_name"])
                                

                            cur = db.get_psycopg2_cursor()
                            cur.execute("select * from metagenome_annotation_insert(%s::text[], %s::numeric[], %s::numeric[], %s::text[], %s::numeric[], %s::text[]);",
                                        (a_numbers, metagenome_sequence_ids, contig_ncbi_tax_ids, contig_tax_names,
                                        orf_ncbi_tax_ids, orf_tax_names))

                            db.psycopg2_connection.commit()
                            logger.info("ingested: %s annotations" % (total_annotations))
                            # logger.info("annotations" % (the_annotations))
                            metagenome_annotations = {}
                            orf_ids = []
                            total_annotations = total_annotations+batch_count
                            batch_count = 0
                    count = count + 1

            except csv.Error as e:
                logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))
                print(e)

if __name__ == '__main__':
    main()