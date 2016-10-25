"""
Proteomics DB.
The goal of the db schema is to make proteomics data ingesta and analysis
modular, flexible, and fast.
"""
from conf import db_config
from sqlalchemy import (Table, Column, Integer, String, ForeignKey,
                        Float)
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import mapper, relationship
from sqlalchemy.orm import object_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.util import has_identity

from proteomics import models

db_connection=db_config.DB_CONN_PREFIX+'+psycopg2://+'+db_config.DB_USER+':'+db_config.DB_PASS+'@'+db_config.DB_HOST+':'+db_config.DB_PORT+'/'+db_config.DB_NAME
engine = create_engine(db_connection)

def get_connection():
    return engine.connect()

def get_session(bind=None):
    if not bind:
        bind = get_connection()
    return sessionmaker(bind=bind)()

def init_db(bind=engine):
    metadata.create_all(bind=bind, checkfirst=True)

def clear_db(bind=engine):
    metadata.drop_all(bind=bind)

def get_session_w_external_trans(orig_session):
    con = orig_session.bind.connect()
    trans = con.begin()
    new_session = sessionmaker()(bind=con)
    return con, trans, new_session

def get_obj_state(obj):
    if object_session(obj) is None and not has_identity(obj):
        return 'transient'
    elif object_session(obj) is not None and not has_identity(obj):
        return 'pending'
    elif object_session(obj) is None and has_identity(obj):
        return 'detached'
    elif object_session(obj) is not None and has_identity(obj):
        return 'persistent'

def get_batched_results(q, batch_size):
    """ Return an iterator that batches query results.
    This is usually used in order to avoid overloading memory.
    Ideally we would use window functions, but not all dbs support this.
    """
    total = q.count()
    returned = 0
    batch = None
    while returned < total:
        if not batch:
            batch = q.limit(batch_size).offset(returned)
        for row in batch:
            returned += 1
            if (returned % batch_size) == 0:
                batch = None
            yield row


# Define tables.
metadata = MetaData()

tables = {}

tables['Taxon'] = Table(
    'taxon', metadata,
    Column('id', String, primary_key=True),
)
mapper(models.Taxon, tables['Taxon'])

tables['TaxonDigest'] = Table(
    'taxon_digest', metadata,
    Column('id', Integer, primary_key=True),
    Column('taxon_id', String, ForeignKey('taxon.id')),
    Column('digest_id', Integer, ForeignKey('digest.id')),
)
mapper(models.TaxonDigest, tables['TaxonDigest'], properties={
    'taxon': relationship(models.Taxon),
    'digest': relationship(models.Digest)
})

tables['Protein'] = Table(
    'protein', metadata,
    Column('id', Integer, primary_key=True),
    Column('sequence', String, index=True),
    Column('mass', Float),
)
mapper(models.Protein, tables['Protein'])

tables['ProteinDigest'] = Table(
    'protein_digest', metadata,
    Column('id', Integer, primary_key=True),
    Column('protein_id', Integer, ForeignKey('protein.id'), index=True),
    Column('digest_id', Integer, ForeignKey('digest.id'), index=True)
)
mapper(models.ProteinDigest, tables['ProteinDigest'], properties={
    'protein': relationship(models.Protein),
    'digest': relationship(models.Digest)
})

tables['TaxonProtein'] = Table(
    'taxon_protein', metadata,
    Column('id', Integer, primary_key=True),
    Column('taxon_id', String, ForeignKey('taxon.id'), index=True),
    Column('protein_id', Integer, ForeignKey('protein.id'), index=True),
    Column('metadata', String),
)
mapper(models.TaxonProtein, tables['TaxonProtein'], properties={
    'protein': relationship(models.Protein),
    'taxon': relationship(models.Taxon),
})

tables['Peptide'] = Table(
    'peptide', metadata,
    Column('id', Integer, primary_key=True),
    Column('sequence', String, index=True),
    Column('mass', Float),
)
mapper(models.Peptide, tables['Peptide'])

tables['ProteinDigestPeptide'] = Table(
    'protein_digest_peptide', metadata,
    Column('id', Integer, primary_key=True),
    Column('peptide_id', Integer, ForeignKey('peptide.id'), index=True),
    Column('protein_digest_id', Integer, ForeignKey('protein_digest.id'),
           index=True),
    Column('count', Integer),
)
mapper(
    models.ProteinDigestPeptide, 
    tables['ProteinDigestPeptide'], 
    properties={
        'peptide': relationship(models.Peptide),
        'protein_digest': relationship(models.ProteinDigest),
    }
)

tables['TaxonDigestPeptide'] = Table(
    'taxon_digest_peptide', metadata,
    Column('id', Integer, primary_key=True),
    Column('peptide_id', Integer, ForeignKey('peptide.id'), index=True),
    Column('taxon_digest_id', Integer, ForeignKey('taxon_digest.id'), 
           index=True),
    Column('count', Integer),
)
mapper(
    models.TaxonDigestPeptide, 
    tables['TaxonDigestPeptide'],
    properties={
        'peptide': relationship(models.Peptide),
        'taxon_digest': relationship(models.TaxonDigest),
    }
)

tables['Protease'] = Table(
    'protease', metadata,
    Column('id', String, primary_key=True),
    Column('cleavage_rule', String),
)
mapper(models.Protease, tables['Protease'])

tables['Digest'] = Table(
    'digest', metadata,
    Column('id', Integer, primary_key=True),
    Column('protease_id', String, ForeignKey('protease.id')),
    Column('max_missed_cleavages', Integer),
    Column('min_acids', Integer),
    Column('max_acids', Integer),
)
mapper(models.Digest, tables['Digest'], properties={
    'protease': relationship(models.Protease),
})
