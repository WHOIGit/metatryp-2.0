from proteomics import dbfrom proteomics.models import Taxonimport loggingimport csvimport timedef main():    start_time = time.time()    logger = logging.getLogger('update_ncbi_taxons')    logger.addHandler(logging.StreamHandler())    logger.setLevel(logging.INFO)    filename = 'C:/Users/David/Projects/metatryp2_python3/metatryp-2.0/data_files/Meren_MAGs_with_NCBI_taxon_IDs2.csv'    with open(filename, 'rt') as f:        reader = csv.reader(f)        try:            count=1;            for row in reader:                #get the ncbi_id and genome (MAG or SAG) name and update the specialized_assembly table with this info                if count > 1:                    ncbi_id = row[9]                    genome_name = row[1]                    print("Genome Name: ",genome_name)                    print("ID: ", ncbi_id)                    cur = db.get_psycopg2_cursor()                    cur.execute("update specialized_assembly set ncbi_id = %s where genome_name = %s", (ncbi_id, genome_name))                    db.psycopg2_connection.commit()                    s = "."                    seq = (row[2], row[3], row[4],row[5], row[6], row[7], row[8])                    hierachy = s.join(seq)                    hierachy = hierachy.replace(" ", "_")                    hierachy = hierachy.replace("-", "_")                    hierachy = hierachy.replace("(", "_")                    hierachy = hierachy.replace(")", "_")                    hierachy = hierachy.replace("[", "")                    hierachy = hierachy.replace("]", "")                    hierachy = hierachy.replace("'", "")                    hierachy = hierachy.replace("'", "")                    print(hierachy)                    cur = db.get_psycopg2_cursor()                    cur.execute("select nt.ncbi_id from ncbi_taxonomy nt where nt.ncbi_id = %s;", (ncbi_id,))                    taxon_result = cur.fetchone()                    if taxon_result is None:                        ncbi_url = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+ncbi_id                        print(ncbi_url)                        ##################out of order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        cur.execute("insert into ncbi_taxonomy(ncbi_id, tax_kingdom, tax_phylum, tax_class, "                                     "tax_order, tax_family, tax_genus, tax_species, hierachy, ncbi_url)"                                     " values(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",  (ncbi_id,row[2],row[3],row[4],row[5],row[6],row[7],row[8], hierachy, ncbi_url))                        db.psycopg2_connection.commit()                count = count + 1        except csv.Error as e:           logger.error('file %s, line %d: %s' % (filename, reader.line_num, e))if __name__ == '__main__':    main()