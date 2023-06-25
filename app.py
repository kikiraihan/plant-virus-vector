from flask import Flask, Response, request, jsonify
from flask_cors import CORS
from handlers.praproses import praproses
from handlers.visualisasi import pie_proporsi
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

if __name__ == "__main__":
    app.run(host="0.0.0.0", port="8009")