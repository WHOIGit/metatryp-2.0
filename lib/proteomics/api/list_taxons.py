"""
name: list_taxons.py


commissioned by : Dr. Makoto Saito, 2013-03

authorship: dgaylord, 2016-02

description: This script queries the database for all available taxons

Outputs: a JSON file containing a list of taxons:

"""

"""
Imports and setup.
"""
from proteomics import db
from proteomics.models import (Peptide, TaxonDigestPeptide, TaxonDigest)
import logging
from flask import Flask, request
from flask_restful import Resource, Api


class QueryTaxons(Resource):
    def get(self):
        args = request.args
        logger = logging.getLogger('list_taxons')
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)

        connection = db.get_connection()
        session = db.get_session(bind=connection)

        results = {
            "taxons": []
        }

        # Execute query and rerurn results.
        query = "select * from taxon"
        queryResults = session.execute(query)

        for row in queryResults:
            results["taxons"].append(row[0])
        return results


