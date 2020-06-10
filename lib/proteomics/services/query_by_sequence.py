import logging
from proteomics import db


class QueryBySequenceTask(object):
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
                print('GENOMIC RESULTS')
                print('search sequence,id,name')
                #cur.execute("select taxon_digest_taxon_id from genomic_query_by_peptide_sequence(%s)", (seq,))
                cur.execute("select id, genome_name from genomic_query_taxon_by_peptide_sequence_new(%s)", (seq,))
                for row in cur.fetchall():
                    print(','.join([str(s) for s in [seq] + list(row)]))
            if type == 'sa' or type == 'all':
                print('\n');
                print('SPECIALIZED ASSEMBLY RESULTS')
                print('search sequence,genome name, sequence id')
                cur.execute(
                    "select specialized_assembly_name, specialized_assembly_sequence from specialized_assembly_taxon_query_by_peptide_sequence(%s)",
                    (seq,))
                for row in cur.fetchall():
                    print(','.join([str(s) for s in [seq] + list(row)]))
            if type == 'm' or type == 'all':
                print('\n');
                print('METAGENOMNIC RESULTS')
                print('search sequence, metagenome name')
                cur.execute("select metagenome_name from metagenomic_query_by_peptide_sequence(%s)", (seq,))
                for row in cur.fetchall():
                    print(','.join([str(s) for s in [seq] + list(row)]))



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
