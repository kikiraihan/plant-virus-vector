
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

from modul.vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id_via_string
from modul.filterNodeEdge import removeNodeAndEdgeByFilter,removeEdgesNotInNodes
from modul.handler_umum import contains_string_entire_column,contains_string_entire_column_boolean

datas=[
    ('1cucu','Cucumber mosaic virus','Aphididae','Myzus'),
    ('2cri','Tomato chlorosis virus','Aleyrodidae','Bemisia'),
    ('3wai','Maize chlorotic dwarf virus','Cicadellidae','Graminella'),
    ('4beg','Tomato yellow leaf curl China virus','Aleyrodidae','Bemisia'),
    ('5pol','Cereal yellow dwarf virus','Aphididae','Schizaphis'),
    ('6pea','Pea enation mosaic virus 1','Aphididae','Acyrthosiphon pisum'),
    ('7cucur','Cucurbit yellow stunting disorder virus','Aleyrodidae','Bemisia'),
    ('8ten','Rice stripe tenuivirus','Delphacidae','Laodelphax'),
    ('9fiji','Southern rice black-streaked dwarf virus','Delphacidae','Sogatella'),
    ('10capchlo','Capsicum chlorosis orthotospovirus','Thripidae','Thrips Palmi'),
    ('11barley','Barley yellow dwarf virus GAV','Aphididae','Sitobion avenae'),
    ('12tospot','Tomato spotted wilt orthotospovirus','Thripidae','Frankliniella'),
]

for data in datas:
    data_,virus_utama,acuan_,ujian_=data # vektor acuan  #data virus
    bobot_ed=1;
    bobot_dc=1;

    print('================================================================================================')
    print('data : ',data_)
    print('virus : ',virus_utama)


    # link enpoint sparql ncbi_ontology
    endpoint_url = 'http://localhost:3030/mydataset/query'
    #1
    #baca data
    df_node=pd.read_csv('dari_praproses/'+data_+'_node.csv',index_col=0) 
    df_edge=pd.read_csv('dari_praproses/'+data_+'_edge.csv',index_col=0)

    # # eksperimen tambahan. bikin fakta tambahan, yaitu relasi virus utama dengan acuan 
    virus_utama=df_node[df_node.virus_utama==True].taxon_id.to_list()
    serangga_acuan=contains_string_entire_column(df_node,acuan_).taxon_id.to_list()
    for i in virus_utama:
        for j in serangga_acuan:
            df_edge.loc[len(df_edge),['source_taxon_id','target_taxon_id','interaction_type']] = [i,j,'pathogenOf']


    if(len(df_node[df_node['group']=="serangga"])<=2):
        print("cuma dua serangga")


    # #### Konversi graf
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
                # print(s_label,'-->', o_label)
                counter+=1
        
        return counter

    # hitung relasi ke virus utama setiap serangga
    virus_utama_ids=list(df_node[df_node['virus_utama']==True].taxon_id)
    for idx,data in df_node[(df_node['group']=='serangga')].iterrows():
        _relasi = bfs_relasi_ke_virus_utama(gnx,data.taxon_id,virus_utama_ids)
        # if(_relasi>0):
        #     print(idx,data.taxon_name,data.taxon_id)
        #     print(_relasi)
        df_node.loc[idx,'relasi_ke_virus_utama'] = _relasi
        # print("=================")
        # update DC pake bobot
        # reset_n=(len(gnx.nodes)-1)/(len([node for node, data in gnx.nodes(data=True) if data.get('group') == "serangga"])-1)
        # results[data.taxon_id] = 1+(results[data.taxon_id]*_relasi*reset_n) #1+(CM*w) #kalo pake jumlah serangga sebagai pembagi
        results[data.taxon_id] = 1+(results[data.taxon_id]*_relasi) #1+(CM*w)

    dc_serangga=[]
    for node_id, rank in sorted(results.items(), key=lambda item: item[1], reverse=True):
        if allnodes[node_id]['group'] in ['serangga']:
            label = gnx.nodes[node_id]['label']
            dc_serangga.append((rank, label, node_id))

    #visualisasi data
    print("mengecek BFS Degree tertinggi")
    for to_search  in [dc_serangga[0][2]]:#["NCBI:7038"]:#"NCBI:33377","NCBI:7036",
        for edge in nx.bfs_edges(gnx.to_undirected(), source=to_search, depth_limit=1):
                s_id, o_id = edge

                s_label = gnx.nodes[s_id]['label'] +' '+s_id
                o_label = gnx.nodes[o_id]['label'] +' '+o_id
                o_grup = gnx.nodes[o_id]['group']

                # skip subjek to_search
                # if (o_grup == 'tanaman'):
                print(s_label,'-->', o_label)
        print('===============')


    print('detail takson serangga DC tertinggi :')
    for i in [dc_serangga[0][2]]:#["NCBI:33377","NCBI:7036","NCBI:7038","NCBI:65032"]:
        print(gnx.nodes[i])



    #5
    # Ambil data NCBI
    # data acuan
    print('ambil data acuan dan ujian')
    data_acuan=get_taxon_vector(acuan_,endpoint_url)
    print("data_acuan : ",data_acuan)
    data_ujian=get_taxon_vector(ujian_,endpoint_url)
    print("data_ujian : ",data_ujian)

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


    # visualisasi
    # Reduce the dimensions of entity embeddings to represent them in a 2D plane.
    # X= UMAP().fit_transform(embeddings)
    # df_umap=pd.DataFrame(X,columns=['feature-vector-1','feature-vector-2'])


    # text=[]
    # labels=[]
    # for x in transformer._entities:
    #     if(x!="http://pyRDF2Vec#SERANGGA_ACUAN"):
    #         text.append(gnx.nodes[x.split("#")[-1]]['famili'].split('_')[-1])
    #         labels.append(gnx.nodes[x.split("#")[-1]]['label'])
    #     else:
    #         text.append("#TITIK_VEKTOR_ACUAN")
    #         labels.append("#TITIK_VEKTOR_ACUAN")
    # df_umap['text']=text
    # df_umap['labels']=labels

    # # # gnx.nodes[x.split("#")[-1]]['label']
    # # df_umap['text']=list(map(lambda x: x.split("#")[-1],transformer._entities))
    # fig = px.scatter(df_umap, x='feature-vector-1',y='feature-vector-2',text='text',hover_name='labels')
    # fig.update_traces(textposition='top center')
    # fig.update_layout(
    #     height=650,
    #     title_text='reduced word2vec visualization'
    # )
    # fig.show()



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
        

    # DC
    print('hasil DC')
    data_to_count=data_to_count.sort_values('dc_result',ascending=False).reset_index(drop=True)
    print(data_to_count[['label','entity','dc_result']])

    # Pengujian dc
    print('pengujian DC')
    for urutan in [0]:
        takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
        id_hasil=data_to_count.iloc[urutan].entity
        cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
        cek_ujian= { k:v for k,v in data_ujian if k in takson }
        print(acuan_,'->', data_)
        print('ujian ',cek_ujian)
        print('hasil ',cek_hasil)
        cek=0
        for i in reversed(takson):
            cekk=cek_hasil[i]==cek_ujian[i]
            cek+=cekk
            print(i, cekk)
        print(cek/len(takson))


    # ED
    print('hasil ED')
    data_to_count=data_to_count.sort_values('ed_result',ascending=True).reset_index(drop=True)
    print(data_to_count[['label','entity','ed_result']])
    # Pengujian ed
    print('pengujian ED')
    for urutan in [0]:
        takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
        id_hasil=data_to_count.iloc[urutan].entity
        cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
        cek_ujian= { k:v for k,v in data_ujian if k in takson }
        print(acuan_,'->', data_)
        print('ujian ',cek_ujian)
        print('hasil ',cek_hasil)
        cek=0
        for i in reversed(takson):
            cekk=cek_hasil[i]==cek_ujian[i]
            cek+=cekk
            print(i, cekk)
        print(cek/len(takson))


    # final score
    print('pengujian final score')
    data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)
    print(data_to_count[['label','dc_result','ed_result','result']])
    # Pengujian kombinasi
    print('pengujian final score')
    for urutan in [0]:
        takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
        id_hasil=data_to_count.iloc[urutan].entity
        cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
        cek_ujian= { k:v for k,v in data_ujian if k in takson }
        print(acuan_,'->', data_)
        print('ujian ',cek_ujian)
        print('hasil ',cek_hasil)
        cek=0
        for i in reversed(takson):
            cekk=cek_hasil[i]==cek_ujian[i]
            cek+=cekk
            print(i, cekk)
        print(cek/len(takson))
    
    print('================================================================================================')

