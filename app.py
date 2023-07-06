from flask import Flask, Response, request, jsonify
from flask_cors import CORS
from handlers.praproses import praproses
from handlers.visualisasi import pie_proporsi, plotly_graph, get_pos_and_nx_data
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
def stream(virus):
    return Response(praproses(virus), mimetype='text/event-stream')

# tidak terpakai stw ini
@app.route('/visualisasi/pie/', methods=["POST"])
def visualisasi_pie_proses():
    data = request.get_json()
    df_node = pd.DataFrame.from_dict(json.loads(data['node']), orient='columns')
    df_edge = pd.DataFrame.from_dict(json.loads(data['edge']), orient='columns')
    virus_txt = data['search']
    return pie_proporsi(df_node)


@app.route('/visualisasi/plotly/graph', methods=["POST"])
def visualisasi_plotly_graph():
    data = request.get_json()
    split_node=data['node']
    split_edge=data['edge']
    print(split_node.keys())

    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')

    return plotly_graph(df_node,df_edge)

@app.route('/visualisasi/pos', methods=["POST"])
def visualisasi_get_pos_from_df():
    data = request.get_json()
    split_node=data['node']
    split_edge=data['edge']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    return get_pos_and_nx_data(df_node,df_edge)



if __name__ == "__main__":
    app.run(host="0.0.0.0", port="8009")