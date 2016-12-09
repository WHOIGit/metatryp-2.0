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
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
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
        # Get session.
        self.session = db.get_session(bind=self.get_connection())
        self.digest = self.session.merge(self.digest)


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
        taxon = self.session.query(Taxon).get(taxon_id)

      #  cur = db.get_psycopg2_cursor();
      #  cur.execute("select t.id from Taxon t where t.id = %s;", (taxon_id,))
      #  taxon_results = cur.fetchone()
      #  db.psycopg2_connection.commit()
      #  file_logger.info("Taxon Results '%s'" % taxon_results)
        if not taxon:
            taxon = Taxon(id=taxon_id)
            self.session.add(taxon)
            self.session.commit()
            self.stats['Taxon'] += 1
            file_logger.info("Created taxon '%s'" % taxon_id)

        # Check if TaxonDigest record exists in db.
        taxon_digest = (
            self.session.query(TaxonDigest)
                .filter(TaxonDigest.taxon == taxon)
                .filter(TaxonDigest.digest == self.digest)
        ).first()

        # If digest has been run on this taxon, don't do anything.
        if taxon_digest:
            file_logger.info((
                                 "Taxon '%s' has already been digested with"
                                 " digest '%s', skipping."
                             ) % (taxon_id, self.digest))
            return

        # Otherwise create a new TaxonDigest.
        if not taxon_digest:
            taxon_digest = TaxonDigest(taxon=taxon, digest=self.digest)
            self.session.add(taxon_digest)
            self.session.commit()
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

        # Generate TaxonDigestPeptides
        q = (
            self.session.query(
                Peptide.id,
                func.sum(ProteinDigestPeptide.count)
            )
                .select_from(Peptide)
                .join(ProteinDigestPeptide)
                .join(ProteinDigest)
                .join(Digest)
                .join(Protein)
                .join(TaxonProtein)
                .join(Taxon)
                .filter(Taxon.id == taxon.id)
                .filter(Digest.id == self.digest.id)
                .group_by(Peptide.id)
        )

        batch_size = 1e4
        tdp_batch = []
        tdp_counter = 0
        for row in db.get_batched_results(q, batch_size):
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
        for protein in (
                self.session.query(Protein)
                        .filter(Protein.sequence.in_(
                    [sequence for metadata, sequence in batch])
                )
        ):
            existing_proteins[protein.sequence] = protein

        # Initialize collection of undigested proteins.
        undigested_proteins = {}
        digested_proteins = {}
        protein_sequences = []
        protein_masses = []
        if existing_proteins:
            for protein in (
                    self.session.query(Protein)
                            .filter(Protein.id.in_(
                        [protein.id for protein in existing_proteins.values()]))
                            .join(ProteinDigest)
                            .filter(ProteinDigest.digest == self.digest)
            ):
                digested_proteins[protein.sequence] = protein
        for protein in existing_proteins.values():
            if protein.sequence not in digested_proteins:
                undigested_proteins[protein.sequence] = protein

        # Create proteins which do not exist in the db and add to undigested
        # collection.

        start_time = time.time()
        num_new_proteins = 0
        for metadata, sequence in batch:
            if sequence not in existing_proteins:
                try:
                    mass = get_aa_sequence_mass(sequence)
                    #protein = Protein(sequence=sequence, mass=mass)

                except Exception as e:
                    logger.exception("Error processing protein, skipping")
                    continue
                #self.session.add(protein)
                num_new_proteins += 1
                if (sequence not in protein_sequences):
                    protein_sequences.append(sequence)
                    protein_masses.append(mass)
                #undigested_proteins[sequence] = protein
                #existing_proteins[sequence] = protein
        logger.info("creating %s new proteins..." % (
            num_new_proteins))
        cur = db.get_psycopg2_cursor();
        cur.execute("select * from protein_insert(%s, %s);", (protein_sequences, protein_masses))
     #   all_records = cur.fetchall()
        for record in cur:
            try:
                protein = Protein(id=record[0],sequence=record[1], mass=record[2])

            except Exception as e:
                logger.exception("Error processing protein, skipping")
                continue
            # self.session.add(protein)
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
                undigested_batch[protein] = {
                    'peptide_sequences': peptide_sequences,
                    'protein_digest': protein_digest,
                }
                if (peptide_counter > 1e4):
                    self.process_peptide_batch(undigested_batch, logger)
                    peptide_counter = 0
            self.process_peptide_batch(undigested_batch, logger)

        # Create taxon protein instances in bulk.
        taxon_protein_dicts = []
        for metadata, sequence in batch:
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
        logger.info("Creating %s new taxon proteins..." % (
            len(taxon_protein_dicts)))
        self.session.execute(
            db.tables['TaxonProtein'].insert(), taxon_protein_dicts)
        self.session.commit()
        self.stats['TaxonProtein'] += len(taxon_protein_dicts)

    def process_peptide_batch(self, batch, logger=None):
        if not logger:
            logger = self.logger

        # Assemble combined peptide sequences and protein digests.
        combined_peptide_sequences = set()
        combined_protein_digests = []
        protein_ids = []
        protein_digests = []
        for protein, data in batch.items():
            for sequence in data['peptide_sequences']:
                combined_peptide_sequences.add(sequence)
            combined_protein_digests.append(data['protein_digest'])
            pd = data['protein_digest']
            protein_ids.append(pd.protein.id)
            protein_digests.append(pd.digest.id)

        # Add protein digests to db.
        logger.info("Creating %s new protein digests..." % (
            len(combined_protein_digests)))
        #self.session.add_all(combined_protein_digests)
        #self.session.commit()
        cur = db.get_psycopg2_cursor();
        cur.execute("select * from protein_digest_insert(%s, %s);", (protein_ids, protein_digests))

        db.psycopg2_connection.commit()
        self.stats['ProteinDigest'] += len(combined_protein_digests)

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
        cur = db.get_psycopg2_cursor();
        cur.execute("select peptide_insert(%s, %s);", (peptide_sequences, peptide_masses))

        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("peptide time elapsed: %s" % (total_time))
        #self.session.execute(db.tables['Peptide'].insert(), peptide_dicts)
        #self.session.commit()
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
        for protein, data in batch.items():
            peptides_histogram = defaultdict(int)
            for sequence in data['peptide_sequences']:
                peptides_histogram[sequence] += 1
            data['peptide_histogram'] = peptides_histogram
            # Update number of peptide instances.
            num_peptide_instances += len(peptides_histogram)

        # Create protein digest peptide instances in bulk.
        logger.info("Creating %s new protein digest peptides..." % (
            num_peptide_instances))
        start_time = time.time()
        pdp_batch = []
        pdp_peptide_ids = []
        pdp_protein_digest_ids = []
        pdp_peptide_count = []
        pdp_counter = 0
        for protein, data in batch.items():
            for sequence, count in data['peptide_histogram'].items():
                pdp_counter += 1
                peptide = existing_peptides[sequence]
                pdp_batch.append({
                    'peptide_id': peptide.id,
                    'protein_digest_id': data['protein_digest'].id,
                    'count': count,
                })
                pdp_peptide_ids.append(peptide.id)
                pdp_protein_digest_ids.append(data['protein_digest'].id)
                pdp_peptide_count.append(count)
                if (pdp_counter % 1e4) == 0:
                    #self.session.execute(
                    #    db.tables['ProteinDigestPeptide'].insert(),
                    #    pdp_batch)
                    #self.session.commit()
                    cur = db.get_psycopg2_cursor();
                    cur.execute("select protein_peptide_digest_insert(%s, %s, %s);", (pdp_peptide_ids, pdp_protein_digest_ids, pdp_peptide_count))
                    db.psycopg2_connection.commit()
        #self.session.execute(
        #    db.tables['ProteinDigestPeptide'].insert(), pdp_batch)
        #self.session.commit()
        cur = db.get_psycopg2_cursor();
        cur.execute("select protein_digest_peptide_insert(%s, %s, %s);",
                    (pdp_peptide_ids, pdp_protein_digest_ids, pdp_peptide_count))
        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("protein digest time elapsed: %s" % (total_time))
        self.stats['ProteinDigestPeptide'] += num_peptide_instances

    def update_existing_peptides_(self, sequences, existing_peptides):
        if not sequences:
            return
        for peptide in (
                self.session.query(Peptide).filter(Peptide.sequence.in_(sequences))
        ):
            existing_peptides[peptide.sequence] = peptide

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
                'taxon_digest_id': taxon_digest.id,
                'peptide_id': row[0],
                'count': row[1],
            })
            taxon_digest_ids.append(taxon_digest.id)
            pepdide_ids.append(row[0])
            peptide_count.append(row[1])
       # self.session.execute(db.tables['TaxonDigestPeptide'].insert(),
        #                     dicts)
        #self.session.commit()
        cur = db.get_psycopg2_cursor();
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