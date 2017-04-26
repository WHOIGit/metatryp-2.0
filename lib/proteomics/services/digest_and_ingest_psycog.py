from proteomics.models import (Taxon, Protein,
                               TaxonProtein, ProteinDigest, Peptide,
                               ProteinDigestPeptide, TaxonDigestPeptide,
                               TaxonDigest, Digest, Protease)
from proteomics import db
from proteomics.util.digest import cleave
from proteomics.util.logging_util import LoggerLogHandler
from proteomics.util.mass import get_aa_sequence_mass
from proteomics.util import fasta
import os
import hashlib
import logging

from collections import defaultdict
from datetime import datetime
import time


class DigestAndIngestTask(object):
    def __init__(self, logger=logging.getLogger(), fasta_paths=[],
                 digest=None, get_connection=None, **kwargs):
        self.logger = logger
        self.fasta_paths = fasta_paths
        self.digest = digest
        self.get_connection = get_connection

    def run(self):

        # Initialize stats dict.
        self.stats = defaultdict(int)

        # Process FASTA files.
        for path in self.fasta_paths:
            self.process_fasta_file(path)

        self.logger.info("Digest and ingest task complete.")
        return self.stats

    def process_fasta_file(self, path):
        base_msg = "Processing file '%s'..." % path
        file_logger = self.get_child_logger(id(path), base_msg,
                                            self.logger)

        # Get taxon from filename.
        taxon_id = os.path.splitext(os.path.basename(path))[0]

        # Get taxon object from db or create a new one.
        #taxon = self.session.query(Taxon).get(taxon_id)
        taxon = Taxon(id=taxon_id)
        cur = db.get_psycopg2_cursor();
        cur.execute("select t.id from taxon t where t.id = %s;", (taxon_id,))
        taxon_result = cur.fetchone()
        db.psycopg2_connection.commit()
        if taxon_result is None:
            #add a taxon to the DB

            #  file_logger.info("Taxon Results '%s'" % taxon_results)
            cur.execute("insert into taxon (id) values(%s);", (taxon_id,))
            db.psycopg2_connection.commit()
            self.stats['Taxon'] += 1
            file_logger.info("Created taxon '%s'" % taxon_id)

        # Check if TaxonDigest record exists in db.
        cur.execute("select t.id from taxon_digest t where t.taxon_id = %s and t.digest_id = %s;", (taxon_id,self.digest.id,))
        db.psycopg2_connection.commit()
        taxon_digest_result = cur.fetchone()
        taxon_digest = TaxonDigest(taxon=taxon, digest=self.digest)
        if taxon_digest_result:
        # If digest has been run on this taxon, don't do anything.
        #if taxon_digest:
            file_logger.info((
                                 "Taxon '%s' has already been digested with"
                                 " digest '%s', skipping."
                             ) % (taxon_id, self.digest))
            return
        else:
            # Otherwise create a new TaxonDigest.
            cur.execute("insert into taxon_digest (taxon_id, digest_id) values(%s,%s);", (taxon_digest.taxon.id,taxon_digest.digest.id,))
            db.psycopg2_connection.commit()
            self.stats['TaxonDigest'] += 1

        # Process protein sequences in batches.
        file_logger.info("Counting # of protein sequences...")
        num_proteins = 0
        for metadata, sequence in fasta.read(path):
            num_proteins += 1
        file_logger.info("%s total protein sequences." % num_proteins)
        batch_size = 500
        batch_counter = 0
        batch = []
        protein_logger = self.get_child_logger(
            "%s_proteins" % id(file_logger), "Processing proteins...",
            file_logger
        )
        protein_logger.info("")
        for metadata, sequence in fasta.read(path):
            batch.append((metadata, sequence,))
            batch_counter += 1
            if (batch_counter % batch_size) == 0:
                self.process_protein_batch(
                    batch, taxon, logger=protein_logger)
                protein_logger.info(
                    ("%s of %s (%.1f%%)") % (
                        batch_counter, num_proteins,
                        100.0 * batch_counter / num_proteins
                    )
                )
                batch = []
        self.process_protein_batch(
            batch, taxon, logger=protein_logger)

        batch_size = 1e4
        cur = db.get_psycopg2_cursor();
        cur.execute("select * from get_peptide_count(%s, %s);", (self.digest.id, taxon.id,))


        tdp_batch = []
        tdp_counter = 0
        for row in cur.fetchall():
            tdp_counter += 1
            tdp_batch.append(row)
            if (tdp_counter % batch_size) == 0:
                self.process_taxon_digest_peptide_batch(
                    taxon_digest, tdp_batch, logger=file_logger)
                tdp_batch = []
        self.process_taxon_digest_peptide_batch(
            taxon_digest, tdp_batch, logger=file_logger)
        self.stats['TaxonDigestPeptide'] += tdp_counter

        self.logger.info("Done processing file '%s'" % path)

    def get_checksum(self, path):
        sha1 = hashlib.sha1()
        with open(path, 'rb') as f:
            while True:
                data = f.read(8192)
                if not data: break
                sha1.update(data)
        return sha1.hexdigest()

    def process_protein_batch(self, batch, taxon, logger=None):
        """ Process a batch of proteins with the given digest. """
        if not batch:
            return
        if not logger:
            logger = self.logger
        # Get existing proteins by searching for sequences.
        existing_proteins = {}
        existing_protein_ids = []
        cur = db.get_psycopg2_cursor()
        sequences = []
        for metadata, sequence in batch:
            if sequence != "No sequence found":
                sequences.append(sequence)

        cur.execute("select * from protein where protein.sequence in %s", (tuple(sequences),))

        for record in cur.fetchall():
            protein = Protein(id=record[0], sequence=record[1], mass=record[2])
            existing_proteins[protein.sequence] = protein
            existing_protein_ids.append(record[0]);
        db.psycopg2_connection.commit()
        # Initialize collection of undigested proteins.
        undigested_proteins = {}
        digested_proteins = {}
        protein_sequences = []
        protein_masses = []
        #testing now, convert to stored procedure
        if existing_proteins:
            cur = db.get_psycopg2_cursor()
            cur.execute("select * from protein join protein_digest on protein.id = protein_digest.protein_id where protein.id in %s and protein_digest.digest_id = %s", ( tuple(existing_protein_ids), self.digest.id,))

            for record in cur.fetchall():
                protein = Protein(id=record[0], sequence=record[1], mass=record[2])
                digested_proteins[protein.sequence] = protein
            db.psycopg2_connection.commit()
        for protein in existing_proteins.values():
            if protein.sequence not in digested_proteins:
                undigested_proteins[protein.sequence] = protein

        # Create proteins which do not exist in the db and add to undigested
        # collection.

        start_time = time.time()
        num_new_proteins = 0
        for metadata, sequence in batch:
            if sequence not in existing_proteins and sequence != "No sequence found":
                try:
                    mass = get_aa_sequence_mass(sequence)
                except Exception as e:
                    logger.exception("Error processing protein, skipping")
                    continue
                num_new_proteins += 1
                # add sequence and mass to their respective lists to be passed to postgres stored procedure
                if (sequence not in protein_sequences):
                    protein_sequences.append(sequence)
                    protein_masses.append(mass)

        logger.info("creating %s new proteins..." % (
            num_new_proteins))
        cur = db.get_psycopg2_cursor()
        cur.execute("select * from protein_insert(%s, %s);", (protein_sequences, protein_masses))
        # iterate through the protein records returned from the insert and build a protein object
        for record in cur:
            try:
                protein = Protein(id=record[0], sequence=record[1], mass=record[2])
            except Exception as e:
                logger.exception("Error processing protein, skipping")
                continue
            undigested_proteins[record[1]] = protein
            existing_proteins[record[1]] = protein

        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("time elapsed: %s" % (total_time))
        self.stats['Protein'] += num_new_proteins

        # Digest undigested proteins.
        if undigested_proteins:
            num_undigested = len(undigested_proteins)
            logger.info("digesting %s proteins" % num_new_proteins)
            undigested_batch = {}
            peptide_counter = 0
            protein_digests = []

            for protein in undigested_proteins.values():
                protein_digest = ProteinDigest(protein=protein,
                                               digest=self.digest)
                protein_digests.append(protein_digest)
                peptide_sequences = cleave(
                    protein.sequence,
                    self.digest.protease.cleavage_rule,
                    self.digest.max_missed_cleavages,
                    min_acids=self.digest.min_acids,
                    max_acids=self.digest.max_acids,
                )
                peptide_counter += len(peptide_sequences)
                undigested_batch[protein.id] = {
                    'peptide_sequences': peptide_sequences,
                    'protein_digest': protein_digest,
                }
                if (peptide_counter > 1e4):
                    self.process_peptide_batch(undigested_batch, logger)
                    peptide_counter = 0
            self.process_peptide_batch(undigested_batch, logger)

        # Create taxon protein instances in bulk.
        taxon_protein_dicts = []

        taxon_protein_ids = []
        taxon_ids = []
        metadatas = []

        for metadata, sequence in batch:
            if sequence != "No sequence found":
                try:
                    protein = existing_proteins[sequence]
                except Exception as e:
                    logger.exception("Error processing protein, sequence does not"
                                 " exist in db, skipping")
                    continue
                taxon_protein_dicts.append({
                    'protein_id': protein.id,
                    'taxon_id': taxon.id,
                    'metadata': metadata,
                })
                taxon_protein_ids.append(protein.id)
                taxon_ids.append(taxon.id)
                metadatas.append(metadata)
        logger.info("Creating %s new taxon proteins..." % (
            len(taxon_protein_dicts)))

        cur = db.get_psycopg2_cursor()
        cur.execute("select * from taxon_protein_insert(%s, %s, %s);", (taxon_protein_ids, taxon_ids, metadatas))
        db.psycopg2_connection.commit()

        self.stats['TaxonProtein'] += len(taxon_protein_dicts)

    def process_peptide_batch(self, batch, logger=None):
        if not logger:
            logger = self.logger

        # Assemble combined peptide sequences and protein digests.
        combined_peptide_sequences = set()

        protein_ids = []
        digest_ids = []
        protein_digests = []
        protein_digests_dict = {}
        for proteinId, data in batch.items():
            for sequence in data['peptide_sequences']:
                combined_peptide_sequences.add(sequence)

            pd = data['protein_digest']
            protein_ids.append(pd.protein.id)
            digest_ids.append(pd.digest.id)

        cur = db.get_psycopg2_cursor()
        cur.execute("select * from protein_digest_insert(%s, %s);", (protein_ids, digest_ids))
        # iterate through the protein_digest records returned from the insert and build a protein_digest object
        for record in cur:

                try:
                    protein_digest = ProteinDigest(id=record[0], protein=record[1], digest=record[2])

                    protein_digests.append(protein_digest)
                    batch_record = batch.get(record[1])
                    protein_digests_dict[record[1]] = {
                        'peptide_sequences': batch_record['peptide_sequences'],
                        'protein_digest': protein_digest,
                    }
                except Exception as e:
                    logger.exception("Error processing protein digest, skipping")
                    continue

        db.psycopg2_connection.commit()
        self.stats['ProteinDigest'] += len(protein_digests)

        # Get existing peptides.
        existing_peptides = {}
        existing_peptides_batch = []
        existing_peptides_counter = 0
        for sequence in combined_peptide_sequences:
            existing_peptides_counter += 1
            existing_peptides_batch.append(sequence)
            if (existing_peptides_counter % 500) == 0:
                self.update_existing_peptides_(
                    existing_peptides_batch, existing_peptides)
                existing_peptides_batch = []
        self.update_existing_peptides_(
            existing_peptides_batch, existing_peptides)

        # Create non-existent peptides in bulk.
        start_time = time.time()
        num_new_peptides = 0
        peptide_dicts = []
        peptide_sequences = []
        peptide_masses = []
        for sequence in combined_peptide_sequences:
            if sequence not in existing_peptides:
                num_new_peptides += 1
                mass = get_aa_sequence_mass(sequence)
                peptide_dicts.append({
                    'sequence': sequence,
                    'mass': mass,
                })
                peptide_sequences.append(sequence)
                peptide_masses.append(mass)
        logger.info("Creating %s new peptides..." % num_new_peptides)
        cur = db.get_psycopg2_cursor()
        cur.execute("select peptide_insert(%s, %s);", (peptide_sequences, peptide_masses))

        db.psycopg2_connection.commit()


        self.stats['Peptide'] += num_new_peptides

        # Get newly created peptide objects and add to existing peptides.
        created_peptides_batch = []
        created_peptides_counter = 0
        for peptide_dict in peptide_dicts:
            created_peptides_counter += 1
            created_peptides_batch.append(peptide_dict['sequence'])
            if (created_peptides_counter % 500) == 0:
                self.update_existing_peptides_(created_peptides_batch,
                                               existing_peptides)
                created_peptides_batch = []
        self.update_existing_peptides_(
            created_peptides_batch, existing_peptides)

        # Create histogram of peptide sequence occurences for each protein.
        num_peptide_instances = 0

        for proteinId, data in protein_digests_dict.items():
            peptides_histogram = defaultdict(int)
            for sequence in data['peptide_sequences']:
                peptides_histogram[sequence] += 1
            data['peptide_histogram'] = peptides_histogram
            # Update number of peptide instances.
            num_peptide_instances += len(peptides_histogram)
        total_time = time.time() - start_time
        logger.info("peptide time elapsed: %s" % (total_time))
        # Create protein digest peptide instances in bulk.
        logger.info("Creating %s new protein digest peptides..." % (
            num_peptide_instances))

        start_time = time.time()
        pdp_batch = []
        pdp_peptide_ids = []
        pdp_protein_digest_ids = []
        pdp_peptide_count = []
        pdp_counter = 0
        for proteinId, data in protein_digests_dict.items():
            for sequence, count in data['peptide_histogram'].items():
                pdp_counter += 1
                peptide = existing_peptides[sequence]
                pdp_peptide_ids.append(peptide.id)
                pdp_protein_digest_ids.append(data['protein_digest'].id)
                pdp_peptide_count.append(count)
        total_time = time.time() - start_time
        logger.info("protein digest loop time elapsed: %s" % (total_time))
        cur = db.get_psycopg2_cursor()
        cur.execute("select protein_digest_peptide_insert(%s, %s, %s);",
                    (pdp_peptide_ids, pdp_protein_digest_ids, pdp_peptide_count))
        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("protein digest time elapsed: %s" % (total_time))
        self.stats['ProteinDigestPeptide'] += num_peptide_instances

    def update_existing_peptides_(self, sequences, existing_peptides):
        if not sequences:
            return

        cur = db.get_psycopg2_cursor()
        cur.execute("select * from peptide where peptide.sequence in %s", (tuple(sequences),))
        for record in cur.fetchall():
            peptide = Peptide(id=record[0], sequence=record[1], mass=record[2])
            existing_peptides[peptide.sequence] = peptide
        db.psycopg2_connection.commit()

    def process_taxon_digest_peptide_batch(self, taxon_digest, batch,
                                           logger=None):
        if not logger:
            logger = self.logger
        start_time = time.time()
        dicts = []
        taxon_digest_ids = []
        pepdide_ids = []
        peptide_count = []
        logger.info("Creating %s new taxon digest peptides..." % (
            len(batch)))
        for row in batch:
            dicts.append({
                'taxon_digest_id': row[2],
                'peptide_id': row[0],
                'count': row[1],
            })
            taxon_digest_ids.append(row[2])
            pepdide_ids.append(row[0])
            peptide_count.append(row[1])

        cur = db.get_psycopg2_cursor()
        cur.execute("select taxon_digest_peptide_insert(%s, %s, %s);",
                    (pepdide_ids, taxon_digest_ids, peptide_count))
        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("taxon digest time elapsed: %s" % (total_time))

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
