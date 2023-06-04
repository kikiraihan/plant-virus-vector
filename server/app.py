from flask import Flask, Response
from flask_cors import CORS
from handlers.praproses import praproses
import time,json

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

if __name__ == "__main__":
    app.run(host="0.0.0.0", port="8009")