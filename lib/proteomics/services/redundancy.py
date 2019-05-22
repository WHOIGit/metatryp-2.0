import itertools
import logging
from proteomics import db
from proteomics.models import Redundancy_Helper,Specialized_Assembly

def count_common_peptides(taxon_digests=[], logger=None):

    taxon_digest_ids = [taxon_digest.id for taxon_digest in taxon_digests]
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from taxon_count_common_peptides(%s, %s)", (taxon_digest_ids, len(taxon_digests)))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_common_peptides_ids(taxon_digest_ids=[], logger=None):
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from taxon_count_common_peptides(%s, %s)", (taxon_digest_ids, len(taxon_digest_ids)))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_peptide_union(taxon_digests=[], logger=None):
    taxon_digest_ids = [taxon_digest.id for taxon_digest in taxon_digests]
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from taxon_count_peptide_union(%s)", (taxon_digest_ids,))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_peptide_union_ids(taxon_digest_ids=[], logger=None):
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from taxon_count_peptide_union(%s)", (taxon_digest_ids,))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_common_peptides_sa(specialized_assemblies=[], logger=None):
    sa_ids = [sa.id for sa in specialized_assemblies]
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sadp_count_common_peptides(%s, %s)", (sa_ids, len(sa_ids)))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_common_peptides_sa_ids(sa_ids=[], logger=None):
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sadp_count_common_peptides(%s, %s)", (sa_ids, len(sa_ids)))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_peptide_union_sa(specialized_assemblies=[], logger=None):
    sa_ids = [sa.id for sa in specialized_assemblies]
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sadp_count_peptide_union(%s)", (sa_ids,))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count
def count_peptide_union_sa_ids(sa_ids=[], logger=None):
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sadp_count_peptide_union(%s)", (sa_ids,))
    count = len(cur.fetchall());
    db.psycopg2_connection.commit()
    return count

def count_common_peptides_combined(sa_ids=[], td_ids=[], logger=None):
    #logger.info("In Count Common Pepdides Combined for SA %s and TD %s" % (sa_ids, td_ids))
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sa_taxon_count_common_peptides(%s, %s, %s)", (sa_ids, td_ids, 2))
    # cur.execute("SELECT count(peptide_id), peptide_id from "\
    #     "((SELECT distinct sadp.peptide_id AS peptide_id "\
    #     "FROM specialized_assembly_digest_peptide sadp "\
    #     "JOIN specialized_assembly_sequence ON specialized_assembly_sequence.id = sadp.specialized_assembly_sequence_id "\
    #     "join specialized_assembly on specialized_assembly.id = specialized_assembly_sequence.specialized_assembly_id "\
    #     "WHERE specialized_assembly_sequence.specialized_assembly_id = any(array[%s])) "\
    #     "union all "\
    #     "(SELECT distinct taxon_digest_peptide.peptide_id AS peptide_id "\
    #     "FROM taxon_digest_peptide JOIN taxon_digest ON taxon_digest.id = taxon_digest_peptide.taxon_digest_id "\
    #     "WHERE taxon_digest.id = any(array[%s]))) as peptide_unions "\
    #     "group by peptide_id "\
    #     "having count(peptide_id) = %s;", (sa_ids, td_ids, 2))
    count = len(cur.fetchall());
    #logger.info("Combined Count %s " % (count))
    db.psycopg2_connection.commit()
    return count

def count_peptide_union_combined(sa_ids=[], td_ids=[], logger=None):
    #logger.info("In Union Common Pepdides Combined for SA %s and TD %s" % (sa_ids, td_ids))
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from sa_taxon_count_peptide_union(%s, %s)", (sa_ids, td_ids))
    # cur.execute("SELECT peptide_id from "\
    #     "((SELECT distinct sadp.peptide_id AS peptide_id "\
	#     "FROM specialized_assembly_digest_peptide sadp "\
	#     "JOIN specialized_assembly_sequence ON specialized_assembly_sequence.id = sadp.specialized_assembly_sequence_id "\
	#     "join specialized_assembly on specialized_assembly.id = specialized_assembly_sequence.specialized_assembly_id "\
	#     "WHERE specialized_assembly_sequence.specialized_assembly_id = any(array[%s])) "\
	#     "union all "\
    #     "(SELECT distinct taxon_digest_peptide.peptide_id AS peptide_id "\
	#     "FROM taxon_digest_peptide JOIN taxon_digest ON taxon_digest.id = taxon_digest_peptide.taxon_digest_id "\
	#     "WHERE taxon_digest.id = any(array[%s]))) as peptide_unions "\
    #     "group by peptide_id;", (sa_ids, td_ids))
    count = len(cur.fetchall());
    #logger.info("Union Count %s " % (count))
    db.psycopg2_connection.commit()
    return count

def generate_redundancy_tables(taxon_digests=[], logger=None):
    """
    Generates tables of:
        - counts: counts of peptides in common between pairs of taxon digests
        - union percents: |td1 ^ td2|/|td1 + td2|
        - pairwise percents: |td1 ^ td2|/|td1|
    """

    if not logger:
        logger = logging.getLogger()

    # Generate pairs.
    combinations = [c for c in itertools.combinations(taxon_digests, 2)]

    # Get redundancies and sums by querying db.
    # Calculate percents.
    values = {
        'intersection_counts': {},
        'union_percents': {},
        'individual_percents': {},
        'individual_counts': {},
    }
    redundancies = {}
    percents = {}
    for combo in combinations:
        logger.info("Counting peptides in common for %s" % (str([
            "(taxon: %s, digest: %s)" % (
                taxon_digest.taxon, taxon_digest.digest
            ) for taxon_digest in combo])))
        # Sorted combo for keying.
        combo_key = get_td_combo_key(combo)
        # Get intersection counts and union percentages.
        num_in_intersection = count_common_peptides( combo, logger)
        num_in_union = count_peptide_union(combo, logger)
        if num_in_union:
            percent_in_common = 100.0 * num_in_intersection/num_in_union
        else:
            percent_in_common = 0
        values['intersection_counts'][combo_key] = num_in_intersection
        values['union_percents'][combo_key] = percent_in_common
        # Get individual counts and percentages.
        for td in combo:
            if td not in values['individual_counts']:
                values['individual_counts'][td] = count_common_peptides([td], logger)
            num_in_td = values['individual_counts'][td]
            if num_in_td:
                values['individual_percents'][(td,combo_key,)] = \
                        100.0 * num_in_intersection/num_in_td

    # Sort taxon digests.
    sorted_taxon_digests = sorted(taxon_digests, key=lambda td: td.taxon)

    # Assemble tables.
    tables = {
        'intersection_counts': [],
        'union_percents': [],
        'individual_percents': [],
        'individual_counts': []
    }

    # Assemble individual counts table.
    for i, td in enumerate(sorted_taxon_digests):
        label = "|%s|" % td.taxon
        value = values['individual_counts'][td]
        tables['individual_counts'].append([label, value])

    # Assemble intersection_count and union_pct tables.
    # These tables have one row per combination.
    for td1, td2 in combinations:
        combo_key = get_td_combo_key((td1,td2,))
        # Intersection count.
        intersection_label = '|%s ^ %s|' % (td1.taxon, td2.taxon)
        intersection_count = values['intersection_counts'][combo_key]
        tables['intersection_counts'].append(
            [intersection_label, intersection_count])
        # Union percents.
        union_label = '|%s U %s|' % (td1.taxon, td2.taxon)
        union_pct_label = "%s/%s" % (intersection_label, union_label)
        union_pct = values['union_percents'][combo_key]
        tables['union_percents'].append([union_pct_label, union_pct])

    # Assemble individual percents table.
    # This table has one row per permutation.
    for i, td1 in enumerate(sorted_taxon_digests):
        for j, td2 in enumerate(sorted_taxon_digests):
            combo_key = get_td_combo_key((td1,td2,))
            pct_key = (td1,combo_key,)
            if pct_key in values['individual_percents']:
                label = "|%s ^ %s|/|%s|" % (td1.taxon, td2.taxon,
                                            td1.taxon)
                value = values['individual_percents'][pct_key]
                tables['individual_percents'].append([label, value])

    return tables

def generate_redundancy_tables_sa(specialized_assemblies=[], logger=None):
    """
       Generates tables of:
           - counts: counts of peptides in common between pairs of specialized assemblies
           - union percents: |td1 ^ td2|/|td1 + td2|
           - pairwise percents: |td1 ^ td2|/|td1|
       Needs to be revised if more than one digest is ever used
       """

    if not logger:
        logger = logging.getLogger()

    # Generate pairs.
    combinations = [c for c in itertools.combinations(specialized_assemblies, 2)]

    # Get redundancies and sums by querying db.
    # Calculate percents.
    values = {
        'intersection_counts': {},
        'union_percents': {},
        'individual_percents': {},
        'individual_counts': {},
    }
    redundancies = {}
    percents = {}
    count = 0
    current_sa = Specialized_Assembly(id=None)
    for combo in combinations:
        if combo[0].id != current_sa.id or current_sa == None:
            current_sa = combo[0]
            count = count+1
            logger.info("MAG Count %s" % count)
        logger.info("Counting peptides in common for %s" % (str([
            "(specialized assembly: %s)" % (
                specialized_assembly.genome_name
            ) for specialized_assembly in combo])))
                # Sorted combo for keying.
        combo_key = get_sa_combo_key(combo)

        # Get intersection counts and union percentages.
        num_in_intersection = count_common_peptides_sa(combo, logger)
        num_in_union = count_peptide_union_sa(combo, logger)
        if num_in_union:
            percent_in_common = 100.0 * num_in_intersection / num_in_union
        else:
            percent_in_common = 0
        values['intersection_counts'][combo_key] = num_in_intersection
        values['union_percents'][combo_key] = percent_in_common
        # Get individual counts and percentages.
        for sa in combo:
            if sa not in values['individual_counts']:
                values['individual_counts'][sa] = count_common_peptides_sa([sa], logger)
            num_in_td = values['individual_counts'][sa]
            if num_in_td:
                values['individual_percents'][(sa, combo_key,)] = \
                    100.0 * num_in_intersection / num_in_td

    # Sort specialized assemblies.
    sorted_specialized_assemblies = sorted(specialized_assemblies, key=lambda sa: sa.genome_name)

    # Assemble tables.
    tables = {
        'intersection_counts': [],
        'union_percents': [],
        'individual_percents': [],
        'individual_counts': []
    }

    # Assemble individual counts table.
    for i, sa in enumerate(sorted_specialized_assemblies):
        label = "|%s|" % sa.genome_name
        value = values['individual_counts'][sa]
        tables['individual_counts'].append([label, value])

    # Assemble intersection_count and union_pct tables.
    # These tables have one row per combination.
    for sa1, sa2 in combinations:
        combo_key = get_sa_combo_key((sa1, sa2,))
        # Intersection count.
        intersection_label = '|%s ^ %s|' % (sa1.genome_name, sa2.genome_name)
        intersection_count = values['intersection_counts'][combo_key]
        tables['intersection_counts'].append(
            [intersection_label, intersection_count])
        # Union percents.
        union_label = '|%s U %s|' % (sa1.genome_name, sa2.genome_name)
        union_pct_label = "%s/%s" % (intersection_label, union_label)
        union_pct = values['union_percents'][combo_key]
        tables['union_percents'].append([union_pct_label, union_pct])

    # Assemble individual percents table.
    # This table has one row per permutation.
    for i, sa1 in enumerate(sorted_specialized_assemblies):
        for j, sa2 in enumerate(sorted_specialized_assemblies):
            combo_key = get_sa_combo_key((sa1, sa2,))
            pct_key = (sa1, combo_key,)
            if pct_key in values['individual_percents']:
                label = "|%s ^ %s|/|%s|" % (sa1.genome_name, sa2.genome_name,
                                            sa1.genome_name)
                value = values['individual_percents'][pct_key]
                tables['individual_percents'].append([label, value])

    return tables

def generate_redundancy_tables_combined(taxon_digests=[], specialized_assemblies=[], logger=None):
    """
       Generates tables of:
           - counts: counts of peptides in common between pairs of specialized assemblies
           - union percents: |td1 ^ td2|/|td1 + td2|
           - pairwise percents: |td1 ^ td2|/|td1|
       Needs to be revised if more than one digest is ever used
       """

    if not logger:
        logger = logging.getLogger()

    # Generate pairs.
   # combinations = [c for c in itertools.combinations(specialized_assemblies, 2)]
    combinations = []
    td_helpers = []
    sa_helpers =[]
    for td in taxon_digests:
        logger.info("TD: %s " % td.id)
        td_helpers.append(Redundancy_Helper(td.id,td.taxon,"Genome"))
    for sa in specialized_assemblies:
        logger.info("SA: %s " % sa.id)
        sa_helpers.append(Redundancy_Helper(sa.id,sa.genome_name,"Specialized Assembly"))

    redundancy_helpers = td_helpers + sa_helpers
    combinations_orig = [c for c in itertools.combinations(redundancy_helpers, 2)]
    # Get redundancies and sums by querying db.
    # Calculate percents.
    values = {
        'intersection_counts': {},
        'union_percents': {},
        'individual_percents': {},
        'individual_counts': {},
    }
    logger.info("Before len: %s" % len(combinations))
    for combo in combinations_orig:
        td_ids = []
        sa_ids = []
        for rh in combo:
            if rh.type == "Genome":
                td_ids.append(rh.id)
            elif rh.type == "Specialized Assembly":
                sa_ids.append(rh.id)
        if len(td_ids) == 1 and len(sa_ids) == 1:
            # check the taxaon digests
            combinations.append(combo)

    logger.info("After len: %s" % len(combinations))
    for combo in combinations:
        logger.info("Combo %s" % (str([
            "( %s, %s)" % (
                rh.genome_name, rh.type
            ) for rh in combo])))
    for combo in combinations:
        s_ids = []
        t_ids = []
        if combo[0].id != current_sa.id or current_sa == None:
            current_sa = combo[0]
            count = count+1
            logger.info("MAG Count %s" % count)
        logger.info("Counting peptides in common for %s" % (str([
            "(genome: %s)" % (
                rh.genome_name
            ) for rh in combo])))
        # Sorted combo for keying.
        combo_key = get_combined_combo_key(combo)

        # Get intersection counts and union percentages.
        for rh in combo:
            if rh.type == "Genome":
                t_ids.append(rh.id)
            elif rh.type == "Specialized Assembly":
                s_ids.append(rh.id)

        #get the combined results
        num_in_intersection = count_common_peptides_combined(s_ids, t_ids, logger)
        num_in_union = count_peptide_union_combined(s_ids, t_ids, logger)
        if num_in_union:
            percent_in_common = 100.0 * num_in_intersection / num_in_union
        else:
            percent_in_common = 0
        values['intersection_counts'][combo_key] = num_in_intersection
        values['union_percents'][combo_key] = percent_in_common
        #Get individual counts and percentages.
        for rh in combo:
            if rh not in values['individual_counts']:
                if rh.type == "Genome":
                    values['individual_counts'][rh] = count_common_peptides_ids([rh.id], logger)
                elif rh.type == "Specialized Assembly":
                    values['individual_counts'][rh] = count_common_peptides_sa_ids([rh.id], logger)
            num_in_td = values['individual_counts'][rh]
            if num_in_td:
                values['individual_percents'][(rh, combo_key,)] = \
                    100.0 * num_in_intersection / num_in_td

    # Sort specialized assemblies.
    sorted_redundancy_helpers = sorted(redundancy_helpers, key=lambda rh: rh.genome_name)

    # Assemble tables.
    tables = {
        'intersection_counts': [],
        'union_percents': [],
        'individual_percents': [],
        'individual_counts': []
    }

    # Assemble individual counts table.
    logger.info("Sorted Redundancy Helpers: %s" % sorted_redundancy_helpers)
    for i, rh in enumerate(sorted_redundancy_helpers):

        label = "|%s|" % rh.genome_name
        value = values['individual_counts'][rh]
        tables['individual_counts'].append([label, value])

    # Assemble intersection_count and union_pct tables.
    # These tables have one row per combination.
    for sa1, sa2 in combinations:
        combo_key = get_combined_combo_key((sa1, sa2,))
        # Intersection count.
        intersection_label = '|%s ^ %s|' % (sa1.genome_name, sa2.genome_name)
        intersection_count = values['intersection_counts'][combo_key]
        tables['intersection_counts'].append(
            [intersection_label, intersection_count])
        # Union percents.
        union_label = '|%s U %s|' % (sa1.genome_name, sa2.genome_name)
        union_pct_label = "%s/%s" % (intersection_label, union_label)
        union_pct = values['union_percents'][combo_key]
        tables['union_percents'].append([union_pct_label, union_pct])

    # Assemble individual percents table.
    # This table has one row per permutation.
    for i, sa1 in enumerate(sorted_redundancy_helpers):
        for j, sa2 in enumerate(sorted_redundancy_helpers):
            combo_key = get_sa_combo_key((sa1, sa2,))
            pct_key = (sa1, combo_key,)
            if pct_key in values['individual_percents']:
                label = "|%s ^ %s|/|%s|" % (sa1.genome_name, sa2.genome_name,
                                            sa1.genome_name)
                value = values['individual_percents'][pct_key]
                tables['individual_percents'].append([label, value])

    return tables


def get_td_combo_key(td_combo):
    """ Get a key for a taxon digest combo."""
    return tuple(sorted(td_combo, key=lambda td: td.taxon))

def get_sa_combo_key(sa_combo):
    """ Get a key for a specialized assemblies combo."""
    return tuple(sorted(sa_combo, key=lambda sa: sa.genome_name))

def get_combined_combo_key(redundancy_combo):
    """ Get a key for a redundancy helper combo."""
    return tuple(sorted(redundancy_combo, key=lambda redundancy: redundancy.genome_name))