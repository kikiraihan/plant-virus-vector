from tqdm import tqdm
from pyrdf2vec import RDF2VecTransformer
from pyrdf2vec.graphs import KG, Vertex
from pyrdf2vec.embedders import FastText,Word2Vec
from pyrdf2vec.walkers import RandomWalker
from pyvis.network import Network
from sklearn.manifold import TSNE
from umap import UMAP
from SPARQLWrapper import SPARQLWrapper
from jcopml.plot import plot_missing_value

import plotly.express as px
import pandas as pd
import numpy as np
import requests
import os
import networkx as nx
import matplotlib.pyplot as plt

from vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id


data=[
    ('Myzus','1cucu','Cucumovirus'),
    ('Bemisia','2cri','Crinivirus'),
    ('Graminella','3wai','Waikavirus'),
    ('Aleyrodidae','begomovirus_22_mei','Begomovirus'),
    ('Schizaphis','5pol','Polerovirus'),
    ('Acyrthosiphon','6pea-nama','Enamovirus'),
    ('Frankliniella','7ort','Orthotospovirus'),
    ('Thrips','8capchlo','Orthotospovirus'),
    ('Laodelphax','9ten','Tenuivirus'),
    ('Sogatella','10fiji','Fijivirus'),
    ('Nilaparvata','+11tung','Tungrovirus'),
    ('Myzus','+12pol','Polerovirus'),
    ('Myzus','+13Poty','Potyvirus'),
]

for each in data:


    acuan_,data_,search_=each # vektor acuan  #data virus
    bobot_ed=1;
    bobot_dc=1;
    # link enpoint sparql ncbi_ontology
    endpoint_url = 'http://localhost:3030/mydataset/query'

    #1
    #baca data
    df_node=pd.read_csv('dari_praproses/'+data_+'_node.csv',index_col=0) 
    df_edge=pd.read_csv('dari_praproses/'+data_+'_edge.csv',index_col=0)
    
    #2
    #isi data kosong. mengisi takson kosong, dengan takson sebelumnya, untuk tambalan
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

    # 3
    # pengelompokan
    # Binning of the data based on serangga, virus, tanaman, nogroup
    filter_tanaman = df_node['kingdom']=='NCBI:33090_Viridiplantae' 
    filter_virus = (
        (df_node['superkingdom']=='NCBI:10239_Viruses')
        | (df_node.taxon_name.str.contains(r'\b(virus\w*|\w*virus)\b',case=False))
        | (df_node.taxon_path.str.contains(r'\b(virus\w*|\w*virus)\b', case=False)) 
        #jika berawalan atau berakhiran kata virus
    )
    filter_serangga = ((df_node['class']=='NCBI:50557_Insecta') )

    df_node.loc[filter_tanaman, ['group','color']] = ["tanaman",'#1f922b'] #hijau
    df_node.loc[filter_virus, ['group','color']] = ['virus','#671f92'] #ungu
    df_node.loc[filter_serangga, ['group','color']] = ['serangga','#b22222'] #merah
    df_node.loc[(
        (filter_tanaman==False) & 
        (filter_virus==False) &
        (filter_serangga==False) 
        ),['group','color']] = ['nogroup','#EADDCA'] #abu-abu

    #3
    #konversi graph 
    gnx = nx.MultiDiGraph()
    #node
    for i,a in df_node.iterrows():
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
            group=a['group'],
            color=a['color'],
        )
    #edge
    for i,a in df_edge.iterrows():
        gnx.add_edge(
            a['source_taxon_id'],
            a['target_taxon_id'],
            label=a['interaction_type'],
        )

    #4 
    # Degree Centrality
    results = nx.degree_centrality(gnx)
    allnodes = gnx.nodes

    # TAMBAHAN UNTUK CEK KUALITAS RELASI SERANGGA
    def bfs_relasi_ke_virus_utama(gnx, to_search, virus_utama_ids):
        counter=0
        for edge in nx.bfs_edges(gnx.to_undirected(), source=to_search, depth_limit=2):
            s_id, o_id = edge

            s_label = gnx.nodes[s_id]['label'] +' '+s_id
            o_label = gnx.nodes[o_id]['label'] +' '+o_id

            # skip subjek to_search
            if (s_id == to_search) & (o_id in virus_utama_ids):
                print(s_label,'-->', o_label)
                counter+=1
        
        return counter

    # tandai virus utama
    search_virus,taxon_,ncbi_id_=cek_ncbi_id_by_wiki_id(search_)
    print("keyword virus utama: ",(search_virus,taxon_,ncbi_id_))
    df_node.loc[df_node[taxon_].str.contains(search_virus), ['virus_utama']] = True

    # hitung relasi ke virus utama setiap serangga
    virus_utama_ids=list(df_node[df_node['virus_utama']==True].taxon_id)
    for idx,data in df_node[(df_node['group']=='serangga')].iterrows():
        print(idx,data.taxon_name,data.taxon_id)
        _relasi = bfs_relasi_ke_virus_utama(gnx,data.taxon_id,virus_utama_ids)
        print(_relasi)
        df_node.loc[idx,'relasi_ke_virus_utama'] = _relasi
        print("=================")
        # update DC pake bobot
        # reset_n=(len(gnx.nodes)-1)/(len([node for node, data in gnx.nodes(data=True) if data.get('group') == "serangga"])-1)
        # results[data.taxon_id] = 1+(results[data.taxon_id]*_relasi*reset_n) #1+(CM*w) #kalo pake jumlah serangga sebagai pembagi
        results[data.taxon_id] = 1+(results[data.taxon_id]*_relasi) #1+(CM*w)

    dc_serangga=[]
    for node_id, rank in sorted(results.items(), key=lambda item: item[1], reverse=True):
        if allnodes[node_id]['group'] in ['serangga']:
            label = gnx.nodes[node_id]['label']
            dc_serangga.append((rank, label, node_id))

    #5
    # Ambil data NCBI
    # data acuan
    data_acuan=get_taxon_vector(acuan_,endpoint_url)
    data_acuan

    #6
    #konversi node networkx ke RDF
    URL = "http://pyRDF2Vec"
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

    #7
    #embedding
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

    #8
    #euclidean distance
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
        data_to_count.loc[idx,['dc_result']] = results[row['entity']] if(row['entity']!="SERANGGA_ACUAN") else 1
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
    #drop kolom embedding
    data_to_count.drop(columns=list(range(0,100)), inplace=True)


    #9
    #hitung kombinasi
    for idx, row in data_to_count.iterrows():
        _dc = row['dc_result']
        _ed=( (row['ed_result']) if row['ed_result']!=0 else 1)
        data_to_count.loc[idx,['result']] = _dc/_ed
    # urutkan
    data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)
    data_to_count[['label','dc_result','ed_result','result']]
    data_to_count[['label','dc_result']].sort_values('dc_result',ascending=False).reset_index(drop=True)
    data_to_count[['label','entity','ed_result']].sort_values('ed_result',ascending=True).reset_index(drop=True)


    # # Pengujian
    id_hasil=data_to_count.iloc[0].entity
    cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
    cek_acuan= { k:v for k,v in data_acuan if k in takson }
    print(acuan_,'->', data_)
    print('acuan ',cek_acuan)
    print('hasil ',cek_hasil)
    cek=0
    for i in reversed(takson):
        cekk=cek_hasil[i]==cek_acuan[i]
        cek+=cekk
        print(i, cekk)


    print(cek/len(takson))

