from flask import Flask, Response, request, jsonify
from flask_cors import CORS

from dotenv import load_dotenv
import os

load_dotenv()

from handlers.praproses import praproses
from handlers.proses import proses, get_taxonomy_from_string_handler
from handlers.visualisasi import get_pos_and_nx_data,get_embeddings_entities_and_dict_insect
from handlers.praproses_dari_proses import pra_proses_dari_proses
from handlers.enhancement import get_musuh_alami_handler, get_fact_handler, get_picture_handler, get_abstract_handler, get_wd_id_handler,get_relatives_handler
import time,json
import pandas as pd

app = Flask(__name__)
CORS(app)

@app.route("/")
def index():
    return "Hai, ini percobaan untuk GET request"

@app.route("/submit", methods=["POST"])
def submit():
    return "Hai, ini percobaan untuk POST request"

@app.route('/praproses/<virus>')
def praproses_endpoint(virus):
    return Response(praproses(virus), mimetype='text/event-stream')


@app.route('/visualisasi/pos', methods=["POST"])
def visualisasi_get_pos_from_df():
    data = request.get_json()
    split_node=data['node']
    split_edge=data['edge']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    return get_pos_and_nx_data(df_node,df_edge)

@app.route('/visualisasi/embedding', methods=["POST"])
def visualisasi_get_embeddings_entities_and_dict_insect():
    data = request.get_json()
    split_node=data['node']
    split_edge=data['edge']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    
    df_node, df_edge = pra_proses_dari_proses(df_node, df_edge)
    return get_embeddings_entities_and_dict_insect(df_node)

@app.route('/proses/data_to_count', methods=["POST"])
def get_data_to_count():
    data = request.get_json()
    split_node=data['node']
    split_edge=data['edge']
    acuan_=data['acuan_']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    return proses(df_node,df_edge,acuan_)

@app.route('/proses/get_taxonomy_from_string', methods=["POST"])
def get_taxonomy_from_string():
    data = request.get_json()
    _virus = data['virus']
    return get_taxonomy_from_string_handler(_virus)

@app.route('/enhancement/musuh-alami/<serangga>', methods=["GET"])
def enhancement_musuh_alami(serangga):
    return get_musuh_alami_handler(serangga)

@app.route('/enhancement/wd-id/<serangga>', methods=["GET"])
def enhancement_wd_id(serangga):
    return get_wd_id_handler(serangga)

@app.route('/enhancement/fact/<wd_id>', methods=["GET"])
def enhancement_fact(wd_id):
    return get_fact_handler(wd_id)

@app.route('/enhancement/picture/<wd_id>', methods=["GET"])
def enhancement_picture(wd_id):
    return get_picture_handler(wd_id)

@app.route('/enhancement/abstract/<wd_id>', methods=["GET"])
def enhancement_abstract(wd_id):
    return get_abstract_handler(wd_id)

@app.route('/enhancement/relatives/<ncbi_taxon_id>', methods=["GET"])
def enhancement_relatives(ncbi_taxon_id):
    return get_relatives_handler(ncbi_taxon_id)





# ini tidak terpakai tapi klao mau cari fungsinya di modul.visualisasiHelper.py
# tidak terpakai stw ini
# @app.route('/visualisasi/pie/', methods=["POST"])
# def visualisasi_pie_proses():
#     data = request.get_json()
#     df_node = pd.DataFrame.from_dict(json.loads(data['node']), orient='columns')
#     df_edge = pd.DataFrame.from_dict(json.loads(data['edge']), orient='columns')
#     virus_txt = data['search']
#     return pie_proporsi(df_node)

# tidak terpakai stw ini
# @app.route('/visualisasi/plotly/graph', methods=["POST"])
# def visualisasi_plotly_graph():
#     data = request.get_json()
#     split_node=data['node']
#     split_edge=data['edge']
#     print(split_node.keys())

#     df_node = pd.read_json(json.dumps(split_node), orient='split')
#     df_edge = pd.read_json(json.dumps(split_edge), orient='split')

#     return plotly_graph(df_node,df_edge)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port="8009", debug=True)
    

# pakai SSL
# from flask_sslify import SSLify
# # setelah deklarasi app
# if not app.debug and not app.testing:
#     sslify = SSLify(app)
# # dibagian akhir
# if __name__ == "__main__":
#     app.run(host="0.0.0.0", port="8009", ssl_context=('certificate_ssl/cert.pem', 'certificate_ssl/key.pem'),debug=True)