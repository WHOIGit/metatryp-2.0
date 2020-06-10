import logging
import argparse 
from proteomics import db



class QueryBySequenceLCATask(object):
    def __init__(self, logger=logging.getLogger(),
                 args=None, **kwargs):
        self.logger = logger
        self.args = args

   
    def run(self):
        
        # Read in sequences to query.
        sequences = []
        max_dist = self.args.max_distance

        if self.args.sequence_file:
            with open(self.args.sequence_file, 'rb') as f:
                sequences = [line.strip() for line in f.readlines()]

        elif self.args.sequence:
            sequences = [self.args.sequence]

        # Read in whether to query just genomes ('g'), just metagenomes ('g') or both ('b').
        # Genomes is the default search

        type = 'g'
        if self.args.type:
            type = self.args.type

        if not sequences:
            argparser.error("Provide a query sequence via the '--sequence' option, "
                            "or a set of sequences via the --sequence-file option")

        # Print headers.
        headers = ['query', 'taxon', 'lev_distance', 'match']
        print(','.join(headers))

        # Execute query for each sequence and print results.
        cur = db.get_psycopg2_cursor()

        for seq in sequences:
            if type == 'g' or type == 'all':
                hierachy = []
                results=[]
                taxon_lca=''
                print('GENOMIC RESULTS')
                print('LCA,search sequence,id,name')

                if max_dist == 0:
                    cur.execute("select id, genome_name, hierachy from genomic_query_taxon_by_peptide_sequence_new(%s)", (seq,))
                else:
                    #change this query when we add fuzzy matching
                    cur.execute("select id, genome_name, hierachy from genomic_query_taxon_by_peptide_sequence_new(%s)", (seq,))

                for row in cur.fetchall():
                    hierachy.append(row[0])
                    results.append(','.join([seq, str(row[0]), row[1]]))
                #print(','.join([str(s) for s in [seq] + list(row)]))

                cur.execute("select * from genomic_lca(%s);", [hierachy])
                if hierachy:
                    lca = cur.fetchone()
                    if lca is not None:
                        taxonHierachy = lca[0].split(".")
                        # iterate through taxon hierachy backwards until we find the first record that is not "unclassified"
                        for l in taxonHierachy[::-1]:

                            if l.lower() != "unclassified":
                                taxon_lca = l
                                break
                    else:
                        taxon_lca = "Unknown"
                for r in results:
                    print(','.join([taxon_lca, r]))
            if type == 'sa' or type == 'all':
                hierachy=[]
                results=[]
                print('\n');
                print('SPECIALIZED ASSEMBLY RESULTS')
                print('LCA,search sequence,genome name, sequence id, NCBI ID')
                cur.execute(
                    "select specialized_assembly_name, specialized_assembly_sequence, ncbi_id from specialized_assembly_taxon_query_by_peptide_sequence(%s)",
                    (seq,))
                for row in cur.fetchall():
                    ncbi_id = row[2]
                    hierachy.append(ncbi_id)
                    results.append(','.join([str(s) for s in [seq] + list(row)]))
                if hierachy:
                    assembly_lca = ""
                    # get the least common ancester for hierachy list
                    cur.execute("select * from specialized_assembly_lca(%s);", [hierachy])
                    lca = cur.fetchone()

                    if lca is not None:
                        # print lca[0]
                        taxonHierachy = lca[0].split(".")
                        # iterate through taxon hierachy backwards until we find the first record that is not "unclassified"
                        for l in taxonHierachy[::-1]:

                            if l.lower() != "unclassified":
                                assembly_lca = l
                                break
                    else:
                        assembly_lca = "Unknown"

                for r in results:
                    print(','.join([assembly_lca, r]))
            if type == 'm' or type == 'all':
                hierachy = []
                results = []
                print('\n');
                print('METAGENOMNIC RESULTS')
                print('LCA, search sequence, metagenome name, NCBI ID')
                cur.execute("select metagenome_name, contig_tax_id, orf_tax_id from metagenomic_query_by_peptide_sequence(%s)", (seq,))
                for row in cur.fetchall():
                    name = row[0]
                    contig_tax_id=row[1]
                    orf_tax_id=row[2]
                    if contig_tax_id:
                        prefered_tax_id = contig_tax_id
                    elif orf_tax_id:
                        prefered_tax_id = orf_tax_id
                    hierachy.append(prefered_tax_id)
                    results.append(','.join([seq,name,str(prefered_tax_id)]))
                    # get the least common ancester for hierachy list
                if hierachy:
                    cur.execute("select * from metagenomic_lca2(%s);", [hierachy])
                    lca = cur.fetchone()
                    metagenome_lca = ""
                    if lca is not None:
                        taxonHierachy = lca[0].split(".")
                        # iterate through taxon hierachy backwards until we find the first record that is not "unclassified"
                        for l in taxonHierachy[::-1]:
                            # print l
                            if l.lower() != "unclassified":
                                metagenome_lca = l
                                break
                    else:
                        metagenome_lca = "Unknown"
                for r in results:
                    print(','.join([metagenome_lca, r]))

        db.psycopg2_connection.commit()



    def get_child_logger(self, name=None, base_msg=None, parent_logger=None):
        if not parent_logger:
            parent_logger = self.logger
        logger = logging.getLogger("%s_%s" % (id(self), name))
        formatter = logging.Formatter(base_msg + ' %(message)s.')
        log_handler = LoggerLogHandler(parent_logger)
        log_handler.setFormatter(formatter)
        logger.addHandler(log_handler)
        logger.setLevel(parent_logger.level)
        return logger
