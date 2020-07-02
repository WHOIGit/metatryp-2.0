from proteomics.models import (Metagenome, Protein, 
                               TaxonProtein, ProteinDigest, Peptide, 
                               ProteinDigestPeptide, TaxonDigestPeptide, 
                               TaxonDigest, Digest, Protease)
from proteomics import db
from proteomics.util.logging_util import LoggerLogHandler
import os
import logging
import traceback

from collections import defaultdict


class ClearSpecializedAssemblyDataTask(object):
    def __init__(self, logger=logging.getLogger(), specialized_assembly_ids=[], 
                 get_connection=None, **kwargs):
        self.logger = logger
        self.specialized_assembly_ids = specialized_assembly_ids

    def run(self):

        try:
            # Get session.
            cur = db.get_psycopg2_cursor()
            self.logger.info("Clearing data for specialized assemblies '%s'" % self.specialized_assembly_ids)

            cur.execute("select sa.id from specialized_assembly sa where sa.genome_name in %s", (tuple(self.specialized_assembly_ids),))
            specialized_assembly_results = cur.fetchall()

            if specialized_assembly_results is None:
                self.logger.info("No matching taxons found.  Nothing was changed")
                exit()

            for specialized_assembly in specialized_assembly_results:
                self.logger.info("Clearing data for specialized assembly '%s'" % specialized_assembly)
                self.clear_data_for_specialized_assembly(specialized_assembly)
        except Exception as e:
            self.logger.error("Problem removing specialized_assemblies: '%s'" % e)
            traceback.print_exc()
            db.psycopg2_connection.commit()
        else:
            db.psycopg2_connection.commit()

    def clear_data_for_specialized_assembly(self, specialized_assembly):
        self.logger.info("Clearing data for specialized assembly '%s'" % specialized_assembly)
        cur = db.get_psycopg2_cursor()
        try:
            # Get Specialized Assembly Sequences.
            cur.execute("select sas.id from specialized_assembly_sequence sas where sas.specialized_assembly_id = %s", (specialized_assembly,))
            specialized_assembly_sequences = cur.fetchall()
            if specialized_assembly_sequences is not None:
                # Delete Sequences, Annotations and Digests
                self.logger.info("Deleting Specialized Assembly Sequences and Digests")
                for sas in specialized_assembly_sequences:
                   # self.logger.info("SAS '%s'" % sas)
                    cur.execute("delete from specialized_assembly_digest_peptide where specialized_assembly_sequence_id = %s", (sas[0],))
                    
                db.psycopg2_connection.commit()
                cur = db.get_psycopg2_cursor()
                for sas in specialized_assembly_sequences:
                    cur.execute("delete from specialized_assembly_sequence where id = %s", (sas[0],))
           
                db.psycopg2_connection.commit()
                cur = db.get_psycopg2_cursor()
            # Delete Specialized Assembly
            cur.execute("delete from specialized_assembly where id = %s;", (specialized_assembly,))
        except Exception as e:
            self.logger.error("Problem removing '%s'" % specialized_assembly)
            traceback.print_exc()
        else:
            # Commit the deletes.
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
