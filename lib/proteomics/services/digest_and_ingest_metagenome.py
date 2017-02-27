from proteomics.models import (Taxon, Protein,
                               TaxonProtein, ProteinDigest, Peptide,
                               ProteinDigestPeptide, TaxonDigestPeptide,
                               TaxonDigest, Digest, Protease, Metagenome, Metagenome_Annotation,
                               Metagenome_Sequence_Peptide, Metagenome_Sequence_Digest, Metagenome_Sequence, Metagenome_Taxon)
from proteomics import db
from proteomics.util.digest import cleave
from proteomics.util.logging_util import LoggerLogHandler
from proteomics.util.mass import get_aa_sequence_mass
from proteomics.util import fasta
from proteomics.services.annotations import extract_venter_annotations

import os
import hashlib
import logging
import numpy

from collections import defaultdict
from datetime import datetime
import time


class DigestAndIngestMetagenomeTask(object):
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

        # Get metagenome name from filename.
        # This may not be the best way to do this, might be better to allow user to input a name?
        metagenome_name = os.path.splitext(os.path.basename(path))[0]

        # Get taxon object from db or create a new one.
        metagenome_id = None

        cur = db.get_psycopg2_cursor();
        cur.execute("select m.id from metagenome m where m.name = %s;", (metagenome_name,))
        metagenome_result = cur.fetchone()
        db.psycopg2_connection.commit()
        if metagenome_result is None:
            #add a metagenome to the DB
            #cur.execute("insert into metagenome(name) values(%s);", (metagenome_name,))
            cur.execute("select * from metagenome_insert(%s);", (metagenome_name,))
            db.psycopg2_connection.commit()
            metagenome_result = cur.fetchone()
            self.stats['Metagenome'] += 1
            file_logger.info("Created metagenome '%s'" % metagenome_name)
        metagenome = Metagenome(id=metagenome_result[0], name=metagenome_name)
        # Check if metagenome has already been digested with given digestion agent.
        cur.execute("select md.digest_id from metagenome_sequence_digest md where md.metagenome_sequence_id in (select ms.id from metagenome_sequence ms where ms.metagenome_id = %s) and md.digest_id = %s;", (metagenome.id,self.digest.id,))
        db.psycopg2_connection.commit()
        metagenome_digest_result = cur.fetchone()

        if  metagenome_digest_result:
        # If digest has been run on this metagenome, don't do anything.
        #if taxon_digest:
            file_logger.info((
                                 "Metagenome '%s' has already been digested with"
                                 " digest '%s', skipping."
                             ) % (metagenome_name, self.digest))
            return

        # Process metagenome sequences in batches.
        file_logger.info("Counting # of metagenome sequences...")
        num_proteins = 0
        for metadata, sequence in fasta.read(path):
            num_proteins += 1
        file_logger.info("%s total metagenome sequences." % num_proteins)
        batch_size = 500
        batch_counter = 0
        batch = []
        protein_logger = self.get_child_logger(
            "%s_proteins" % id(file_logger), "Processing metagenome sequences...",
            file_logger
        )
        protein_logger.info("")
        for metadata, sequence in fasta.read(path):
            batch.append((metadata, sequence,))
            batch_counter += 1
            if (batch_counter % batch_size) == 0:
                self.process_metagenome_sequence_batch(
                    batch,metagenome, logger=protein_logger)
                protein_logger.info(
                    ("%s of %s (%.1f%%)") % (
                        batch_counter, num_proteins,
                        100.0 * batch_counter / num_proteins
                    )
                )
                batch = []
        self.process_metagenome_sequence_batch(
            batch, metagenome, logger=protein_logger)



    def get_checksum(self, path):
        sha1 = hashlib.sha1()
        with open(path, 'rb') as f:
            while True:
                data = f.read(8192)
                if not data: break
                sha1.update(data)
        return sha1.hexdigest()

    def process_metagenome_sequence_batch(self, batch, metagenome, logger=None):
        """ Process a batch of metagenome sequences with the given digest. """
        if not batch:
            return
        if not logger:
            logger = self.logger
        # Get existing metagenome sequences (proteins) by searching for sequences.
        existing_sequences = {}
        existing_sequence_ids = []
        cur = db.get_psycopg2_cursor()
        sequences = []
        metadataList = []
        for metadata, sequence in batch:
            sequences.append(sequence)


        cur.execute("select * from metagenome_sequence where metagenome_sequence.sequence in %s", (tuple(sequences),))

        for record in cur.fetchall():
            meta_seq = Metagenome_Sequence(id=record[0], sequence=record[1], metagenome_id=record[2])
            existing_sequences[meta_seq.sequence] = meta_seq
            existing_sequence_ids.append(record[0]);
        db.psycopg2_connection.commit()
        # Initialize collection of undigested proteins.
        undigested_sequences = {}
        digested_sequences = {}
        metagenome_sequences = []
        metagenome_digest_ids = []
        metagenome_ids = []
        metagenome_accesion_ids = {}

        #testing now, convert to stored procedure
        if existing_sequences:
            cur = db.get_psycopg2_cursor()
            #check to see if sequence has previously been digested
            cur.execute(
                "select ms.* from metagenome_sequence ms join metagenome_sequence_digest msd on ms.id = msd.metagenome_sequence_id where ms.id in %s and msd.digest_id = %s",
                (tuple(existing_sequence_ids), self.digest.id,))

            for record in cur.fetchall():
                meta_seq = Metagenome_Sequence(id=record[0], sequence=record[1], metagenome_id=record[2], sequence_id=record[3])
                digested_sequences[meta_seq.sequence] = meta_seq
            db.psycopg2_connection.commit()
        for meta_seq in existing_sequences.values():
       #     if meta_seq.sequence not in digested_sequences:
                undigested_sequences[meta_seq.sequence] = meta_seq

        # Create proteins which do not exist in the db and add to undigested
        # collection.

        start_time = time.time()
        num_new_sequences = 0
        for metadata, sequence in batch:

            if sequence not in existing_sequences:
                 num_new_sequences += 1
                 # add sequence and mass to their respective lists to be passed to postgres stored procedure
                 #if (sequence not in metagenome_sequences):
                 metagenome_sequences.append(sequence)
                 metadataList.append(metadata)
                 metagenome_digest_ids.append(self.digest.id)
                 metagenome_ids.append(metagenome.id)

        logger.info("creating %s new metagenome sequences..." % (
            num_new_sequences))
        cur = db.get_psycopg2_cursor()
        cur.execute("select * from metagenome_sequence_insert(%s, %s, %s);", (metagenome_sequences, metagenome_ids, metadataList))
        # iterate through the protein records returned from the insert and build a protein object
        for record in cur:
            try:
                meta_seq = Metagenome_Sequence(id=record[0], sequence=record[1], metagenome_id=record[2], sequence_id=record[3])
            except Exception as e:
                logger.exception("Error processing metagenome sequence, skipping")
                continue
            undigested_sequences[record[0]] = meta_seq
            existing_sequences[record[0]] = meta_seq
            metagenome_accesion_ids[meta_seq.sequence_id]= meta_seq.id;

        db.psycopg2_connection.commit()
        total_time = time.time() - start_time
        logger.info("time elapsed: %s" % (total_time))
        self.stats['Protein'] += num_new_sequences

        # Digest undigested proteins.
        if undigested_sequences:
            num_undigested = len(undigested_sequences)
            logger.info("digesting %s proteins" % num_new_sequences)
            undigested_batch = {}
            peptide_counter = 0
            meta_seq_digests = []

            for metagenome_sequence in undigested_sequences.values():
                metagenome_sequence_digest = Metagenome_Sequence_Digest(metagenome_sequence=metagenome_sequence,
                                               digest=self.digest)
                meta_seq_digests.append(metagenome_sequence_digest)
                peptide_sequences = cleave(
                    metagenome_sequence.sequence,
                    self.digest.protease.cleavage_rule,
                    self.digest.max_missed_cleavages,
                    min_acids=self.digest.min_acids,
                    max_acids=self.digest.max_acids,
                )
                peptide_counter += len(peptide_sequences)
                undigested_batch[metagenome_sequence.id] = {
                    'peptide_sequences': peptide_sequences,
                    'metagenome_sequence_digest': metagenome_sequence_digest,
                }
                if (peptide_counter > 1e4):
                    self.process_peptide_batch(undigested_batch, logger)
                    peptide_counter = 0
            self.process_peptide_batch(undigested_batch, logger)
            annotations = self.get_venter_annotations(metagenome_accesion_ids, logger)
            logger.info("annotations: %s" % (annotations))

    def process_peptide_batch(self, batch, logger=None):
        if not logger:
            logger = self.logger

        # Assemble combined peptide sequences and metagenome digests.
        combined_peptide_sequences = set()

        metagenome_sequence_ids = []
        digest_ids = []
        metagenome_sequence_digests = []
        metagenome_sequence_digests_dict = {}
        for proteinId, data in batch.items():
            for sequence in data['peptide_sequences']:
                combined_peptide_sequences.add(sequence)

            msd = data['metagenome_sequence_digest']
            metagenome_sequence_ids.append(msd.metagenome_sequence.id)
            digest_ids.append(msd.digest.id)

        cur = db.get_psycopg2_cursor()
        cur.execute("select * from metagenome_sequence_digest_insert(%s, %s);", (metagenome_sequence_ids, digest_ids))
        # iterate through the metagenome_sequence_digest records returned from the insert and build a metagenome_sequence_digest object
        for record in cur:

                try:
                    metagenome_sequence_digest = Metagenome_Sequence_Digest(id=record[0], metagenome_sequence=record[1], digest=record[2])
                    metagenome_sequence_digests.append(metagenome_sequence_digest)
                    batch_record = batch.get(record[1])
                    metagenome_sequence_digests_dict[record[1]] = {
                        'peptide_sequences': batch_record['peptide_sequences'],
                        'metagenome_sequence_digest': metagenome_sequence_digest,
                    }
                except Exception as e:
                    logger.exception("Error processing metagenome sequence digest, skipping")
                    continue

        db.psycopg2_connection.commit()
        self.stats['MetagenomeSequenceDigest'] += len(metagenome_sequence_digests)

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

        for sequenceId, data in metagenome_sequence_digests_dict.items():
            peptides_histogram = defaultdict(int)
            for sequence in data['peptide_sequences']:
                peptides_histogram[sequence] += 1
            data['peptide_histogram'] = peptides_histogram
            # Update number of peptide instances.
            num_peptide_instances += len(peptides_histogram)
        total_time = time.time() - start_time
        logger.info("peptide time elapsed: %s" % (total_time))
        # Create protein digest peptide instances in bulk.
        logger.info("Creating %s new metagenome sequence digest peptides..." % (
            num_peptide_instances))

        start_time = time.time()
        pdp_peptide_ids = []
        pdp_protein_digest_ids = []
        pdp_peptide_count = []
        pdp_counter = 0
        for sequenceId, data in metagenome_sequence_digests_dict.items():
            for sequence, count in data['peptide_histogram'].items():
                pdp_counter += 1
                peptide = existing_peptides[sequence]
                pdp_peptide_ids.append(peptide.id)
                pdp_protein_digest_ids.append(data['metagenome_sequence_digest'].id)
                pdp_peptide_count.append(count)
        total_time = time.time() - start_time
        logger.info("sequence digest loop time elapsed: %s" % (total_time))
        cur = db.get_psycopg2_cursor()
        cur.execute("select metagenome_sequence_digest_peptide_insert(%s, %s, %s);",
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

    def get_venter_annotations(self, seq_ids, logger=None):
        if not logger:
            logger = self.logger
        metagenome_sequence_ids = []
        the_annotations = extract_venter_annotations(self, seq_ids.keys(), logger)
        logger.info("annotations" % (the_annotations))
        accession_numbers =[]
        scaffold_ids =[]
        orf_ids =[]
        orf_nums =[]
        annotations =[]
        gene_names =[]
        orf_tax_levels =[]
        orf_taxonomies =[]
        orf_tax_ids =[]
        contig_tax_ids =[]
        contig_taxonomies =[]
        contig_tax_levels =[]
        metagenome_annotation_ids = []
        for annot in the_annotations:
            accession_numbers.append(annot[0])
            metagenome_sequence_ids.append(seq_ids[annot[0]])
            scaffold_ids.append(annot[1])
            orf_ids.append(annot[2])
            if annot[3] is not None:
                orf_nums.append(int(annot[3]))
            else:
                orf_nums.append(annot[3])
            annotations.append(annot[4])
            gene_names.append(annot[5])
            orf_tax_levels.append(annot[6])
            orf_taxonomies.append(annot[7])
            if annot[8] is not None:
                orf_tax_ids.append(int(annot[8]))
            else:
                orf_tax_ids.append(annot[8])
            contig_tax_levels.append(annot[9])
            contig_taxonomies.append(annot[10])
            if annot[11] is not None:
                contig_tax_ids.append(int(annot[11]))
            else:
                contig_tax_ids.append(annot[11])
        cur = db.get_psycopg2_cursor()
        cur.execute("select * from metagenome_annotation_insert(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",
                    (accession_numbers, metagenome_sequence_ids, scaffold_ids, orf_ids, orf_nums, annotations, gene_names,
                     orf_tax_levels,orf_taxonomies,orf_tax_ids))

        for record in cur.fetchall():
            metagenome_annotation_ids.append(record[0])
        cur.execute("select metagenome_contig_taxonomy_insert(%s, %s, %s, %s);",
                    (metagenome_annotation_ids, contig_tax_ids, contig_taxonomies,contig_tax_levels))
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
