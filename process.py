from tqdm import tqdm
from owlready2 import *
#from sklearn.manifold import TSNE
from pyrdf2vec import RDF2VecTransformer
from pyrdf2vec.graphs import KG, Vertex
from pyrdf2vec.embedders import FastText,Word2Vec
from pyrdf2vec.walkers import RandomWalker
from pyvis.network import Network

import pandas as pd
import numpy as np
import requests
import os
import networkx as nx
import matplotlib.pyplot as plt

from modul.vectorReferenced import get_taxon_vector


def dataframeToNetworkxGraph(df_node,df_edge):
    gnx = nx.MultiDiGraph()
    #node
    for i,a in df_node.iterrows():
        #filter partisi untuk multipartite
        if(a['superkingdom']=='NCBI:10239_Viruses'):
            grup='virus'
            warna='#671f92' #ungu
        elif(a['kingdom']=='NCBI:33090_Viridiplantae'):
            grup='tanaman'
            warna='#1f922b' #hijau
        #ini jadi dilema, ada kalau ingin memasukan artropoda lain misal, laba2 maka pake or. untuk saat ini fokus ke insect
        elif(a['kingdom']=='NCBI:33208_Metazoa' and (a['phylum']=='NCBI:6656_Arthropoda' and a['class']=='NCBI:50557_Insecta') ):
            grup='serangga'
            warna='#b22222' #merah
        else:
            grup='nogroup'
            warna='#EADDCA' #abu-abu
        #mulai disini akan digunakan taksonomi bahasa indonesia pada data.
        gnx.add_node(
            a['taxon_id'],
            label=a['taxon_name'],
            superkingdom=a['superkingdom'],
            kingdom=a['kingdom'],
            filum=a['phylum'],
            kelas=a['class'],
            ordo=a['order'],
            famili=a['family'],
            genus=a['genus'],
            spesies=a['species'],
            group=grup,
            color=warna,
        )
    #edge
    for i,a in df_edge.iterrows():
        gnx.add_edge(
            a['source_taxon_id'],
            a['target_taxon_id'],
            label=a['interaction_type'],
        )
    return gnx

def makeDegreeCentralityList(gnx):
    results = nx.degree_centrality(gnx)
    allnodes = gnx.nodes
    dc_serangga=[]
    for node_id, rank in sorted(results.items(), key=lambda item: item[1], reverse=True):
        if allnodes[node_id]['group'] in ['serangga']:
            label = gnx.nodes[node_id]['label']
            dc_serangga.append((rank, label, node_id))
    return dc_serangga, results, allnodes


def networkxGraphToRdf(gnx, data_acuan, URL):
    CUSTOM_KG = KG()

    takson=[i[0] for i in data_acuan]
    for i in ["superkingdom","kingdom","filum","kelas"]:
        takson.remove(i) 

    # memasukan RDF serangga acuan
    subj = Vertex(f"{URL}#SERANGGA_ACUAN")
    for i,j in data_acuan:
        if(i not in ["superkingdom","kingdom","filum","kelas"]):
            j = j.replace(' ','-')
            obj = Vertex((URL+"#"+j))
            pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
            #pred = Vertex((URL+"#taxon_path_ids"), predicate=True, vprev=subj, vnext=obj)
            CUSTOM_KG.add_walk(subj, pred, obj)

    # proses konversi 
    for index,data in gnx.nodes(data=True):
        if(data['group']=='serangga'): #jika serangga
            subj = Vertex(URL+"#"+index)
            for i in takson:
                #if(isinstance(data[i], str)): #jika dia string atau tidak nan/kosong.
                #if(i not in ["superkingdom","kingdom","filum","kelas"]):
                    id_takson=data[i].replace(' ','-')#.split('_')[0]
                    obj = Vertex((URL+"#"+id_takson))
                    pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
                    #pred = Vertex((URL+"#taxon_path_ids"), predicate=True, vprev=subj, vnext=obj)
                    CUSTOM_KG.add_walk(subj, pred, obj)
    # CUSTOM_KG.literals=[
    #         [f"{URL}#taxon_path_ids"],
    #     ]
    CUSTOM_KG.literals = [[URL+"#"+i] for i in takson]
    return CUSTOM_KG, takson

def embedding(CUSTOM_KG,gnx,URL):
    # Ensure the determinism of this script by initializing a pseudo-random number.
    RANDOM_STATE = 22
    transformer = RDF2VecTransformer(
        # Use one worker threads for Word2Vec to ensure random determinism.
        # Must be used with PYTHONHASHSEED.
        Word2Vec(epochs=1000),
        # Extract a maximum of 10 walks of a maximum depth of 4 for each entity
        # using two processes and use a random state to ensure that the same walks
        # are generated for the entities.
        walkers=[RandomWalker(2, 5, n_jobs=2, with_reverse=False, random_state=RANDOM_STATE)],
        #verbose=1,
    )
    # transformer = RDF2VecTransformer(verbose=1)
    # list entity yang akan diembedd. serangga acuan urutan terakhir
    ent = [ URL+"#"+index for index,data in gnx.nodes(data=True) if(data['group']=='serangga') ] #jika serangga
    ent.append(f"{URL}#SERANGGA_ACUAN")
    # Fit the transformer to the knowledge graph and the entities.
    embeddings, _ = transformer.fit_transform(
        CUSTOM_KG, #the KG
        ent, #entity
    )
    return embeddings, _, transformer


def hitungED(embeddings, entities, dc_dict, gnx):
    # buat dataframe
    data_to_count=pd.DataFrame(embeddings, columns=list(range(0,100)))
    # buat kolom label
    ent=[data['label'] for index,data in gnx.nodes(data=True) if(data['group']=='serangga')] #jika serangga
    ent.append("#SERANGGA_ACUAN")
    data_to_count['label']=ent
    #buat kolom entity
    data_to_count['entity']=[i.replace("http://pyRDF2Vec#","") for i in entities]
    # buat kolom hasil dc
    for idx,row in data_to_count.iterrows(): #jika serangga acuan maka DC di isi nilai 1
        data_to_count.loc[idx,['dc_result']] = dc_dict[row['entity']] if(row['entity']!="SERANGGA_ACUAN") else 1
    #ambil koordinat acuan
    acuan=next(data_to_count[data_to_count['label']=='#SERANGGA_ACUAN'].iterrows())[1]
    acuan=np.array(tuple(acuan[i] for i in range(0,100)))
    acuan
    #hitung ED
    for idx, row in data_to_count.iterrows():
        temp = np.array(tuple(row[i] for i in range(0,100)))
        data_to_count.loc[idx,['ed_result']] = np.linalg.norm(temp - acuan)
    #drop data acuan
    data_to_count.drop(data_to_count[data_to_count.label=="#SERANGGA_ACUAN"].index,inplace=True)
    # minmax scaling dc dan ed
    from sklearn.preprocessing import MinMaxScaler
    for i in ['dc_result', 'ed_result']:
        scaler = MinMaxScaler()
        scaler.fit(data_to_count[i].to_numpy().reshape(1, -1))
        scaler.transform(data_to_count[i].to_numpy().reshape(1, -1))
    
    return data_to_count


def cek_bfs(gnx, to_search):
    for edge in nx.bfs_edges(gnx.to_undirected(), source=to_search, depth_limit=2):
        s_id, o_id = edge

        s_label = gnx.nodes[s_id]['label'] +' '+s_id
        o_label = gnx.nodes[o_id]['label'] +' '+o_id

        # skip subjek to_search
        if s_id == to_search:
            print(s_label,'-->', o_label)

def nx_to_pyviz(gnx):
    warna={
        'virus':'#671f92', #ungu
        'tanaman':'#1f922b', #hijau
        'serangga':'#b22222',#merah
        'nogroup':'#EADDCA', #abu-abu
    }   
    nt = Network('500px', '900px',directed=True,notebook=True)
    # nt.show_buttons(filter_=['physics'])
    nt.toggle_physics(True)
    
    # node
    for i,a in gnx.nodes(data=True):
        nt.add_node(
            i,
            label=a['label'],
            #taxon_path=a['taxon_path'],
            #taxon_path_ids=a['taxon_path_ids'],
            #group=a['group'],
            size='12',
            color=warna[a['group']]
        )

    #edge
    for a in gnx.edges(data=True):
        nt.add_edge(
            a[0],
            a[1],
            label=a[2]['label'],
        )
        
    return nt




def allProcess(data_, acuan_, endpoint_url, bobot_dc=1, bobot_ed=1):
    #1
    #baca data
    print('#1 Baca data interaksi')
    df_node=pd.read_csv('dari_praproses/'+data_+'_node.csv',index_col=0) 
    df_edge=pd.read_csv('dari_praproses/'+data_+'_edge.csv',index_col=0)

    #2
    #isi data kosong
    print('#2 isi data kosong')
    takson=[
        'superkingdom','kingdom','phylum','class','order','family','genus','species'
    ]
    for x,i in enumerate(takson):
        if (i!='superkingdom'): #selain superkingdom update dengan data sebelumnya
            for idx, row in df_node[pd.isnull(df_node[i])].iterrows():
                df_node.loc[idx,[i]] = row[takson[x-1]]+'^'+i
        else: 
            for idx, row in df_node[pd.isnull(df_node[i])].iterrows():
                df_node.loc[idx,[i]] = row[takson[x+1]]+'^'+i

    #3
    #konversi graph 
    print('#3 konversi graph ')
    gnx = dataframeToNetworkxGraph(df_node,df_edge)

    #4 
    # Degree Centrality
    print('#4 Degree Centrality')
    dc_serangga, dc_dict,allnodes = makeDegreeCentralityList(gnx)

    #5
    # data acuan
    print('#5 ambil data acuan')
    data_acuan=get_taxon_vector(acuan_,endpoint_url)

    #6
    #konversi node networkx ke RDF
    print('#6 konversi node networkx ke RDF')
    URL = "http://pyRDF2Vec"
    CUSTOM_KG, takson = networkxGraphToRdf(gnx, data_acuan, URL)

    #7
    #embedding
    print('#7 embedding')
    embeddings, _, transformer = embedding(CUSTOM_KG,gnx,URL)

    #8
    #euclidean distance
    print('#8 hitung euclidean distance')
    data_to_count = hitungED(embeddings,transformer._entities, dc_dict, gnx)
    #drop kolom embedding
    data_to_count.drop(columns=list(range(0,100)), inplace=True)
    
    #9
    #hitung kombinasi
    print('#9 hitung kombinasi')
    for idx, row in data_to_count.iterrows():    
        data_to_count.loc[idx,['result']] = (bobot_dc*row['dc_result'])/( (bobot_ed*row['ed_result']) if row['ed_result']!=0 else 1)
    # urutkan
    data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)

    return data_to_count.to_json(), gnx, data_acuan, takson, df_node, df_edge









