from proteomics.models import (TaxonDigestPeptide, TaxonDigest, Peptide)

import itertools
import logging
from proteomics import db

def count_common_peptides(taxon_digests=[], logger=None):
    taxon_digest_ids = [taxon_digest.id for taxon_digest in taxon_digests]
    cur = db.get_psycopg2_cursor();
    cur.execute("select * from taxon_count_common_peptides(%s, %s)", (taxon_digest_ids, len(taxon_digests)))
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

def get_td_combo_key(td_combo):
    """ Get a key for a taxon digest combo."""
    return tuple(sorted(td_combo, key=lambda td: td.taxon))
