from fastapi import FastAPI, Response, HTTPException, Depends
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from dotenv import load_dotenv
from pymongo import MongoClient
import pandas as pd
import json
import os
import asyncio
import uvicorn

from handlers.praproses import praproses, praproses_from_cache
from handlers.proses import proses, get_taxonomy_from_string_handler
from handlers.visualisasi import get_pos_and_nx_data, get_embeddings_entities_and_dict_insect
from handlers.praproses_dari_proses import pra_proses_dari_proses
from handlers.enhancement import get_musuh_alami_handler, get_fact_handler, get_picture_handler, get_abstract_handler, get_wd_id_handler, get_relatives_handler

app = FastAPI()
# app.mount("/static", StaticFiles(directory="static"), name="static")

# Load environment variables
load_dotenv()
MONGODB_URL = os.environ.get("MONGODB_URL")
JENA_URL_MAINDB = os.environ.get("JENA_URL_MAINDB")

# Setup MongoDB
client = MongoClient(MONGODB_URL)
mongodb = client["thesisdbmongo"]

# Setup CORS
origins = ["*"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Hello, this is a test for a GET request"}

@app.post("/submit")
def read_item():
    return {"message": "Hello, this is a test for a POST request"}

@app.get('/praproses/{virus}')
async def praproses_endpoint(virus: str):
    ncbi_server_url = f'{JENA_URL_MAINDB}/query'
    virus = virus.lower()
    hasil = mongodb.praproses.find_one({"key": virus})
    if hasil is None:
        return StreamingResponse(praproses(virus, mongodb, ncbi_server_url), media_type='text/event-stream')
    print('ini cache')
    return StreamingResponse(praproses_from_cache(hasil), media_type='text/event-stream')

@app.post('/visualisasi/pos')
async def visualisasi_get_pos_from_df(data: dict):
    split_node = data['node']
    split_edge = data['edge']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    return get_pos_and_nx_data(df_node, df_edge)

@app.post('/visualisasi/embedding')
async def visualisasi_get_embeddings_entities_and_dict_insect(data: dict):
    split_node = data['node']
    split_edge = data['edge']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    
    df_node, df_edge = pra_proses_dari_proses(df_node, df_edge)
    return get_embeddings_entities_and_dict_insect(df_node)

@app.post('/proses/data_to_count')
async def get_data_to_count(data: dict):
    ncbi_ontology_url = f'{JENA_URL_MAINDB}/query'
    split_node = data['node']
    split_edge = data['edge']
    acuan_ = data['acuan_']
    df_node = pd.read_json(json.dumps(split_node), orient='split')
    df_edge = pd.read_json(json.dumps(split_edge), orient='split')
    return proses(df_node, df_edge, acuan_, ncbi_ontology_url)

@app.post('/proses/get_taxonomy_from_string')
async def get_taxonomy_from_string(data: dict):
    ncbi_ontology_url = f'{JENA_URL_MAINDB}/query'
    _virus = data['virus']
    return get_taxonomy_from_string_handler(_virus, ncbi_ontology_url)

@app.get('/enhancement/musuh-alami/{serangga}')
async def enhancement_musuh_alami(serangga: str):
    url_ncbi_endpoint = f'{JENA_URL_MAINDB}/query'
    return get_musuh_alami_handler(serangga, url_ncbi_endpoint)

@app.get('/enhancement/wd-id/{serangga}')
def enhancement_wd_id(serangga: str):
    return get_wd_id_handler(serangga)

@app.get('/enhancement/fact/{wd_id}')
def enhancement_fact(wd_id: str):
    return get_fact_handler(wd_id)

@app.get('/enhancement/picture/{wd_id}')
def enhancement_picture(wd_id: str):
    return get_picture_handler(wd_id)

@app.get('/enhancement/abstract/{wd_id}')
def enhancement_abstract(wd_id: str):
    return get_abstract_handler(wd_id)

@app.get('/enhancement/relatives/{ncbi_taxon_id}')
def enhancement_relatives(ncbi_taxon_id: str):
    url_ncbi_endpoint = f'{JENA_URL_MAINDB}/query'
    return get_relatives_handler(ncbi_taxon_id, url_ncbi_endpoint)

# Add other routes for enhancements...

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8009, log_level="info")
