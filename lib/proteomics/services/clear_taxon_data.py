from proteomics.models import (Taxon, Protein, 
                               TaxonProtein, ProteinDigest, Peptide, 
                               ProteinDigestPeptide, TaxonDigestPeptide, 
                               TaxonDigest, Digest, Protease)
from proteomics import db
from proteomics.util.logging_util import LoggerLogHandler
import os
import logging
import traceback

from collections import defaultdict


class ClearTaxonDataTask(object):
    def __init__(self, logger=logging.getLogger(), taxon_ids=[], 
                 get_connection=None, **kwargs):
        self.logger = logger
        self.taxon_ids = taxon_ids

    def run(self):

        try:
            # Get session.
            cur = db.get_psycopg2_cursor();
            self.logger.info("Clearing data for taxon '%s'" % self.taxon_ids)

            cur.execute("select t.id from taxon t where t.id in %s", (tuple(self.taxon_ids),))
            taxon_results = cur.fetchall()

            if taxon_results is None:
                self.logger.info("No matching taxons found.  Nothing was changed")
                exit()

            for taxon in taxon_results:
                self.logger.info("Clearing data for taxon '%s'" % taxon)
                self.clear_data_for_taxon(taxon)
        except Exception as e:
            self.logger.error("Problem removing taxons: '%s'" % e)
            traceback.print_exc()
            db.psycopg2_connection.commit()
        else:
            db.psycopg2_connection.commit()

    def clear_data_for_taxon(self, taxon):
        self.logger.info("Clearing data for taxon '%s'" % taxon)
        cur = db.get_psycopg2_cursor();
        try:
            # Get TaxonDigests.
            cur.execute("select td.id from taxon_digest td where td.taxon_id = %s", (taxon,))
            taxon_digests = cur.fetchall()
            if taxon_digests is not None:
                # Delete TaxonDigestPeptides and TaxonDigests
                self.logger.info("Deleting TaxonDigestPeptides and TaxonDigests")
                for td in taxon_digests:
                    cur.execute("delete from taxon_digest_peptide where taxon_digest_id = %s", (td[0],))
                    cur.execute("delete from taxon_digest where id = %s", (td[0],))

            # Delete TaxonProteins.
            self.logger.info("Deleting TaxonProteins")
            cur.execute("delete from taxon_protein where taxon_id = %s", (taxon,))

            # Delete Taxon
            cur.execute("delete from taxon where id = %s;", (taxon,))
        except Exception as e:
            self.logger.error("Problem removing '%s'" % taxon)
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
