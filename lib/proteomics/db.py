"""
Proteomics DB.
The goal of the db schema is to make proteomics data ingesta and analysis
modular, flexible, and fast.
"""
import psycopg2
from . import db_config




def get_psycopg2_cursor():
    return psycopg2_connection.cursor()

#def get_psycopg2_read_cursor():
 #   return psycopg2_read_connection.cursor()

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



#psycopg2_read_connection =  psycopg2.connect("user="+db_config.DB_READ_USER + " password=" + db_config.DB_READ_PASS + " host=" + db_config.DB_HOST + " dbname="+ db_config.DB_NAME + " sslmode=prefer")
psycopg2_connection =  psycopg2.connect("user="+db_config.DB_USER + " password=" + db_config.DB_PASS+ " host=" + db_config.DB_HOST + " dbname="+ db_config.DB_NAME + " sslmode=prefer")