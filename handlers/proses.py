import os

import pandas as pd
import numpy as np
import networkx as nx

from modul.vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id_via_string
from modul.helper_umum import minmax, std_scale
from modul.grafHelper import _set_networkx_graph
from modul.embeddingHelper import df_serangga_to_rdf, rdf_KG_to_embeddings, df_to_dictionary_taxon
from modul.custom_degree_centrality import degree_centrality_custom
from handlers.praproses_dari_proses import pra_proses_dari_proses

JENA_URL = os.environ.get("JENA_URL")
JENA_URL_MAINDB = os.environ.get("JENA_URL_MAINDB")

print("JENA_URL_MAINDB",JENA_URL_MAINDB)

def proses(df_node,df_edge,acuan_):
    ncbi_ontology_url = f'{JENA_URL_MAINDB}/query'
    df_node, df_edge = pra_proses_dari_proses(df_node, df_edge, acuan_, True)

    if(len(df_node[df_node['group']=="serangga"])<=2):
        return {
            'status': '403',
            'message': 'Hanya dua serangga ditemukan'
        }
    
    print('praproses dari proses selesai')
    
    #konversi graph 
    gnx = _set_networkx_graph(df_node, df_edge)

    # # Analisis Interaksi
    # Degree Centrality Custom
    virus_utama_ids=list(df_node[df_node['virus_utama']==True].taxon_id)
    serangga_ids=list(df_node[df_node['group']=="serangga"].taxon_id)
    results_dc = degree_centrality_custom(gnx,virus_utama_ids,serangga_ids,print_relasi=False)
    allnodes = gnx.nodes

    # # Analisis Taksonomi
    # Ambil  data acuan
    # sudah dinamis kalo di input famili yang di proses familinya
    data_acuan=get_taxon_vector(acuan_,ncbi_ontology_url,False)
    #6 #konversi node networkx ke RDF
    # input : df_node, URL, data_acuan
    URL = "http://pyRDF2Vec"
    df_serangga = df_node[df_node['group']=="serangga"]
    CUSTOM_KG = df_serangga_to_rdf(df_serangga, URL, data_acuan)

    #7 #embedding
    list_serangga=df_node[df_node['group']=="serangga"].taxon_id.to_list()
    # list entity yang akan diembedd. serangga acuan urutan terakhir
    list_of_entities = [ URL+"#"+taxon_id for taxon_id in list_serangga ]
    list_of_entities.append(f"{URL}#SERANGGA_ACUAN")
    transformer, embeddings, _ = rdf_KG_to_embeddings(CUSTOM_KG, list_of_entities)
    # list_of_entities == transformer._entities # True # artinya sama dua2nya

    # dictionary serangga
    dictionary_serangga = df_to_dictionary_taxon(df_serangga)
    # output
    # embeddings, list_of_entities, dictionary_serangga


    # data to count making
    # buat dataframe
    data_to_count=pd.DataFrame(embeddings, columns=list(range(0,100)))
    # buat kolom label
    ent=[data['label'] for index,data in gnx.nodes(data=True) if(data['group']=='serangga')] #jika serangga
    ent.append("#SERANGGA_ACUAN")
    data_to_count['label']=ent
    #buat kolom entity
    data_to_count['entity']=[i.replace("http://pyRDF2Vec#","") for i in transformer._entities]
    # buat kolom hasil dc
    for idx,row in data_to_count.iterrows(): #jika serangga acuan maka DC di isi nilai 1
        data_to_count.loc[idx,['dc_result']] = results_dc[row['entity']] if(row['entity']!="SERANGGA_ACUAN") else 1
    
    #8 euclidean distance
    #ambil koordinat acuan
    acuan=next(data_to_count[data_to_count['label']=='#SERANGGA_ACUAN'].iterrows())[1]
    acuan=np.array(tuple(acuan[i] for i in range(0,100)))
    #hitung ED
    for idx, row in data_to_count.iterrows():
        temp = np.array(tuple(row[i] for i in range(0,100)))
        data_to_count.loc[idx,['ed_result']] = np.linalg.norm(temp - acuan)
    #drop data acuan
    data_to_count.drop(data_to_count[data_to_count.label=="#SERANGGA_ACUAN"].index,inplace=True)
    #drop kolom embedding
    data_to_count.drop(columns=list(range(0,100)), inplace=True)

    # scaling
    # scaling ed_result
    data_to_count['ed_result_scaled'] = minmax(data_to_count['ed_result'])
    # scaling dc 
    data_to_count['dc_result_scaled'] = minmax(data_to_count['dc_result'])

    #9 #hitung kombinasi
    data_to_count['result'] = (1+data_to_count['dc_result_scaled']) / (1+data_to_count['ed_result_scaled'])
    # simple scaling result (final/kombinasi)
    data_to_count['result'] = data_to_count['result'] / data_to_count['result'].max()

    # urutkan berdasarkan final score kombinasi
    data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)

    return {
        'status':200,
        'message': 'Kirim json', 
        'data_to_count':data_to_count.to_json(orient="split"),
    }


def get_taxonomy_from_string_handler(virus_name):
    ncbi_ontology_url = f'{JENA_URL_MAINDB}/query'
    print("ncbi_ontology_url : ",ncbi_ontology_url)
    data = get_taxon_vector(virus_name,ncbi_ontology_url,False)
    print('smppe disini aman')
    if data == False:
        return {
            'status': '404',
            'message': 'Tidak ditemukan taksonomi',
            'taxonomy': None,
        }
    data = dict(data)
    return {
        'status':200,
        'message': 'Kirim json',
        'taxonomy':data,
    }
