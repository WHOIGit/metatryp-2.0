"""
Domain models.
"""

class File(object):
    def __init__(self, id=None, basename=None):
        self.id = id
        self.basename = basename

class FileDigest(object):
    def __init__(self, file_=None, digest=None):
        self.file = file_
        self.digest = digest

class Taxon(object):
    def __init__(self, id=None, metadata=None):
        self.id = id
        self.metadata = metadata

class TaxonDigest(object):
    def __init__(self, id=None, taxon=None, digest=None):
        self.id = id
        self.taxon = taxon
        self.digest = digest

class Protein(object):
    def __init__(self, id=None, sequence=None, mass=None, metadata=None):
        self.id = id
        self.sequence = sequence
        self.mass = mass
        self.metadata = metadata

class ProteinDigest(object):
    def __init__(self, id=None, protein=None, digest=None):
        self.id = id
        self.protein = protein
        self.digest = digest

class TaxonProtein(object):
    """
    A taxon protein is a single occurence of a protein
    that occurs within a taxon's proteome.
    We use protein records because the same protein can appear multiple
    times w/in a proteome.
    """
    def __init__(self, id=None, protein=None, taxon=None, 
                 metadata=None):
        self.id = id
        self.protein = protein
        self.taxon = taxon
        self.metadata = metadata

class Peptide(object):
    def __init__(self, id=None, sequence=None, mass=None, metadata=None):
        self.id = id
        self.sequence = sequence
        self.mass = mass
        self.metadata = metadata

class ProteinDigestPeptide(object):
    """
    A protein digest peptide is a count of how many times a peptide 
    sequence appears in the the digestion of a protein.
    """
    def __init__(self, id=None, peptide=None, protein_digest=None, count=None): 
        self.id = id
        self.peptide = peptide
        self.protein_digest = protein_digest
        self.count = count

class TaxonDigestPeptide(object):
    """
    A taxon digest peptide is a count of how many times a peptide
    sequence appears in the the digestion of a taxon.
    """
    def __init__(self, id=None, peptide=None, taxon_digest=None, count=None): 
        self.id = id
        self.peptide = peptide
        self.taxon_digest = taxon_digest
        self.count = count

class Protease(object):
    def __init__(self, id=None, cleavage_rule=None):
        self.id = id
        self.cleavage_rule = cleavage_rule

class Digest(object):
    def __init__(self, id=None, protease=None, max_missed_cleavages=0,
                 min_acids=0, max_acids=None):
        self.id = id
        self.protease = protease 
        self.max_missed_cleavages = max_missed_cleavages
        self.min_acids = min_acids
        self.max_acids = max_acids

class Metagenome(object):
    def __init__(self, id=None, name=None, sampling_date=None,
                 depth=None, filter_size=None, filter_cutoff=None, expedition_number=None, sample_type=None):
        self.id = id
        self.name = name
        self.sampling_date = sampling_date
        self.depth = depth
        self.filter_size = filter_size
        self.filter_cutoff = filter_cutoff
        self.expedition_number = expedition_number
        self.sample_type = sample_type

class Metagenome_Sequence(object):
    def __init__(self, id=None, sequence=None, metagenome_id=None, sequence_id=None):
        self.id = id
        self.sequence =sequence
        self.metagenome_id = metagenome_id
        self.sequence_id = sequence_id

class Metagenome_Sequence_Digest(object):
    def __init__(self, id=None, metagenome_sequence=None, digest=None):
        self.id = id
        self.metagenome_sequence = metagenome_sequence
        self.digest = digest

class Metagenome_Annotation(object):
    def __init__(self, id=None, metagenome_sequence_id=None, contig_taxonomy_id=None):
        self.id = id
        self.metagenome_sequence_id = metagenome_sequence_id
        self.contig_taxonomy_id = contig_taxonomy_id
        #will add additional annotation fields once list is finalized


class Metagenome_Sequence_Peptide(object):
    def __init__(self, id=None, metagenome_sequence_digest_id=None, peptide_id=None, count=None):
        self.id = id
        self.metagenome_sequence_digest_id =metagenome_sequence_digest_id
        self.peptide_id = peptide_id
        self.count = count

class Metagenome_Taxon(object):
    def __init__(self, id=None, contig_tax_id=None, contig_taxon=None, contig_tax_level=None):
        self.id = id
        self.contig_tax_id =contig_tax_id
        self.contig_taxon = contig_taxon
        self.contig_tax_level = contig_tax_level

class Specialized_Assembly(object):
    def __init__(self, id=None, type_flag=None, status=None,
                 study_name=None, genome_name=None, fasta_file_key=None, sequencing_center=None):
        self.id = id
        self.type_flag = type_flag
        self.status = status
        self.study_name = study_name
        self.genome_name = genome_name
        self.ffasta_file_key = fasta_file_key
        self.sequencing_center = sequencing_center

class Specialized_Assembly_Sequence(object):
    def __init__(self, id=None, sequence=None, specialized_assembly_id=None, sequence_id=None):
        self.id = id
        self.sequence =sequence
        self.specialized_assembly_id = specialized_assembly_id
        self.sequence_id = sequence_id

class Specialized_Assembly_Digest_Peptide(object):
    def __init__(self, id=None, peptide_id=None, specialized_assembly_sequence_id=None, count=None, digest_id_id=None):
        self.id = id
        self.peptide_id = peptide_id
        self.specialized_assembly_sequence_id = specialized_assembly_sequence_id
        self.count = count
        self.digest_id = digest_id_id

class Redundancy_Helper(object):
    def __init__(self, id=None, genome_name=None, type = None):
        self.id = id
        self.genome_name = genome_name
        self.type = type

