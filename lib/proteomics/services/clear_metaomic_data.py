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


class ClearMetaomicDataTask(object):
    def __init__(self, logger=logging.getLogger(), metaomic_ids=[], 
                 get_connection=None, **kwargs):
        self.logger = logger
        self.metaomic_ids = metaomic_ids

    def run(self):

        try:
            # Get session.
            cur = db.get_psycopg2_cursor();
            self.logger.info("Clearing data for meta-omic assemblies '%s'" % self.metaomic_ids)

            cur.execute("select m.id from metagenome m where m.id in %s", (tuple(self.metaomic_ids),))
            metaomic_results = cur.fetchall()

            if metaomic_results is None:
                self.logger.info("No matching taxons found.  Nothing was changed")
                exit()

            for metaomic_assembly in metaomic_results:
                self.logger.info("Clearing data for meta-omic assembly '%s'" % metaomic_assembly)
                self.clear_data_for_metaomic_assembly(metaomic_assembly)
        except Exception as e:
            self.logger.error("Problem removing metaomic_assemblies: '%s'" % e)
            traceback.print_exc()
            db.psycopg2_connection.commit()
        else:
            db.psycopg2_connection.commit()

    def clear_data_for_metaomic_assembly(self, metaomic_assembly):
        self.logger.info("Clearing data for meta-omic assembly '%s'" % metaomic_assembly)
        cur = db.get_psycopg2_cursor()
        try:
            # Get MetagenomeSequences.
            cur.execute("select ms.id from metagenome_sequence ms where ms.metagenome_id = %s", (metaomic_assembly,))
            metagenome_sequences = cur.fetchall()
            if metagenome_sequences is not None:
                # Delete Sequences, Annotations and Digests
                self.logger.info("Deleting TaxonDigestPeptides and TaxonDigests")
                for ms in metagenome_sequences:
        
                    cur.execute("delete from metagenome_sequence_digest_peptide where metagenome_sequence_id = %s", (ms[0],))
                    cur.execute("delete from metagenome_annotations where metagenome_sequence_id = %s", (ms[0],))
                db.psycopg2_connection.commit()
                cur = db.get_psycopg2_cursor()
                for ms in metagenome_sequences:    
                    cur.execute("delete from metagenome_sequence where metagenome_id = %s", (ms[0],))
           
                db.psycopg2_connection.commit()
                cur = db.get_psycopg2_cursor()
            # Delete Metagenome
            cur.execute("delete from metagenome where id = %s;", (metaomic_assembly,))
        except Exception as e:
            self.logger.error("Problem removing '%s'" % metaomic_assembly)
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
