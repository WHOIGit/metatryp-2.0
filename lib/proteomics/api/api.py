from flask import Flask, request
from flask_restful import Resource, Api
from query_by_sequence import QueryBySequence

app = Flask(__name__)
api = Api(app)

#Query sequences for taxon matches:
#You can query the database for exact and fuzzy matches. e.g.:  /sequences/query?sequence=MGFPCNR&max_distance=1
#default max leveinstein distance is 1
api.add_resource(QueryBySequence, '/sequence/query')

if __name__ == '__main__':
    app.run(debug=True)