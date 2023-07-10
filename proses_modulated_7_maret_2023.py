#!/usr/bin/env python
# coding: utf-8

# #### catatan: 
# beberapa data taksonomi memang kosong dari NCBInya, misal data interaksi virus pycv malah terhubung ke genus capsicum bukan spesies capsicum. dan itu memang dari data GloBInya langsung

# In[1]:


# # reset import package
# def reloadPackageOwn():
#     from importlib import reload  
#     import os # we use os.path.join, os.path.basename
#     import sys # we use sys.path
#     import glob # we use glob.glob
#     import importlib # we use importlib.import_module

#     import_folder = os.getcwd()
#     sys.path.append(import_folder) # this tells python to look in `import_folder` for imports
#     for src_file in glob.glob(os.path.join(import_folder, '*.py')):
#         name = os.path.basename(src_file)[:-3]
#         importlib.import_module(name)
#         reload(sys.modules[name])
#         importlib.import_module(name)
        
# reloadPackageOwn()


# In[2]:


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
from modul.helper_umum import contains_string_entire_column,contains_string_entire_column_boolean
#from process import cek_bfs, nx_to_pyviz
from modul.grafHelper import _set_networkx_graph, _plot_nx_by_matplotlib
from modul.visualisasiHelper import embeddingPlot,plotly_graph
from modul.embeddingHelper import df_serangga_to_rdf, rdf_KG_to_embeddings, df_to_dictionary_taxon
from modul.custom_degree_centrality import degree_centrality_custom


# #### Parameter

# In[3]:


data=[
    ('begomovirus_contoh_hasil','Pepper yellow leaf curl virus','Aleyrodidae','Bemisia Tabaci'),
    ('1cucu','Cucumber mosaic virus','Aphididae','Myzus persicae'),
    ('2cri','Tomato chlorosis virus','Aleyrodidae','Bemisia Tabaci'),
    ('3wai','Maize chlorotic dwarf virus','Cicadellidae','Graminella nigrifrons'),
    ('4beg','Tomato yellow leaf curl China virus','Aleyrodidae','Bemisia Tabaci'),
    ('5pol','Cereal yellow dwarf virus','Aphididae','Schizaphis graminum'),
    ('6pea','Pea enation mosaic virus 1','Aphididae','Acyrthosiphon pisum'),
    ('7cucur','Cucurbit yellow stunting disorder virus','Aleyrodidae','Bemisia Tabaci'),
    ('8ten','Rice stripe tenuivirus','Delphacidae','Laodelphax striatellus'),
    ('9fiji','Southern rice black-streaked dwarf virus','Delphacidae','Sogatella furcifera'),
    ('10capchlo','Capsicum chlorosis orthotospovirus','Thripidae','Thrips Palmi'),
    ('11barley','Barley yellow dwarf virus GAV','Aphididae','Sitobion avenae'),
    ('12tospot','Tomato spotted wilt orthotospovirus','Thripidae','Frankliniella occidentalis'),
    ('13svyv','squash vein yellowing virus','Aleyrodidae','Bemisia Tabaci'),
    ('14sbmv','soybean mosaic virus','Aphididae','Aphis glycines'),
    ('15blv','bean leafroll virus','Aphididae','Acyrthosiphon pisum'),
    ('16rgdv','rice gall dwarf virus','Cicadellidae','Recilia dorsalis'), #sedikit
    ('17srbsdv','southern rice black-streaked dwarf virus','Delphacidae','Sogatella furcifera'),
    ('18tsrv','tomato severe rugose virus','Aleyrodidae','Bemisia tabaci'),
    ('19gbnv','groundnut bud necrosis virus','Thripidae','Thrips palmi'),
    ('20wbnv','Watermelon bud necrosis virus','Thripidae','Thrips palmi'),
    # error dibawah ini
    # ('+13Poty','Potyvirus','Aphididae','Myzus'),
    # ('+11tung','Tungrovirus','Nilaparvata','Nilaparvata'),
]

data_,nama_virus,acuan_,ujian_=data[0] # vektor acuan  #data virus
# link enpoint sparql ncbi_ontology
ncbi_ontology_url = 'http://localhost:3030/mydataset/query'


# #### input data

# In[4]:


#1
#baca data
df_node=pd.read_csv('dari_praproses/'+data_+'_node.csv',index_col=0) 
df_edge=pd.read_csv('dari_praproses/'+data_+'_edge.csv',index_col=0)


# In[5]:


df_node[df_node['group']=='virus']
# df_node['group'].unique()

# df_node


# In[6]:


acuan_


# In[7]:


# pra-proses khusus proses
# hapus serangga yg cuma famili (mengikuti acuan). soalnya klo cuma tampil famili apa gunanya?
filter_genus_sampai_species_null=(
    (df_node.genus.isnull()) &
    (df_node.species.isnull()) &
    (df_node['class']=='NCBI:50557_Insecta')
)
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_genus_sampai_species_null], df_node,df_edge)


# In[8]:


print(len(df_edge))
df_edge.drop_duplicates(inplace=True)
print(len(df_edge))


# In[9]:


# pra-proses khusus proses
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


# In[10]:


# pra-proses khusus proses
# eksperimen tambahan. bikin fakta tambahan, yaitu relasi virus utama dengan acuan 
# virus_utama=df_node[df_node.virus_utama==True].taxon_id.to_list()
# serangga_acuan=contains_string_entire_column(df_node,acuan_).taxon_id.to_list()
# print(len(df_edge))
# for i in virus_utama:
#     for j in serangga_acuan:
#         dict = {'source_taxon_id':i,'target_taxon_id':j,'interaction_type':'pathogenOf'}
#         df_edge = pd.concat([pd.DataFrame(dict,index=[0]), df_edge], ignore_index = True)
#         # df_edge.loc[len(df_edge.index),['source_taxon_id','target_taxon_id','interaction_type']] = [i,j,'pathogenOf']
# print(len(df_edge))

#satu saja
virus_utama=df_node[df_node.virus_utama==True].taxon_id.to_list()
serangga_acuan=contains_string_entire_column(df_node,acuan_).taxon_id.to_list()
print(len(df_edge))
for j in serangga_acuan:
    dict = {'source_taxon_id':virus_utama[0],'target_taxon_id':j,'interaction_type':'pathogenOf'}
    df_edge = pd.concat([pd.DataFrame(dict,index=[0]), df_edge], ignore_index = True)
    # df_edge.loc[len(df_edge.index),['source_taxon_id','target_taxon_id','interaction_type']] = [i,j,'pathogenOf']
print(len(df_edge))


# In[11]:


# Ini harusnya di praproses. tapi belum fix baiknya hapus atau tidak
# hapus yang bukan virus utama, terakhir akurasi 0.85, kalau berkurang hapus saja ini
bukan_virus_utama=(df_node['group']=="virus") & (df_node.virus_utama!=True)
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[bukan_virus_utama], df_node,df_edge)


# In[12]:


if(len(df_node[df_node['group']=="serangga"])<=2):
    print("cuma dua serangga")


# # Konversi graf

# In[13]:


#3
#konversi graph 
gnx = _set_networkx_graph(df_node, df_edge)

# # cuma tampilan, visualisasi graf
# _plot_nx_by_matplotlib(gnx)


# # visualisasi Graf

# In[14]:


plotly_graph(gnx)

# # # cuma tampilan, visualisasi pyviz
# from pyvis.network import Network
# nt = Network('500px', '900px',directed=True,notebook=True)
# # nt.show_buttons(filter_=['physics'])
# nt.toggle_physics(True)

# for i,data in df_node.iterrows():
#     nt.add_node(
#         data['taxon_id'],
#         label=data['taxon_name'],
#         superkingdom=data['superkingdom'],
#         kingdom=data['kingdom'],
#         filum=data['phylum'],
#         kelas=data['class'],
#         ordo=data['order'],
#         famili=data['family'],
#         genus=data['genus'],
#         spesies=data['species'],
#         # group=data['group'],
#         color=data['color'],
#         )
    
# for i,data in df_edge.iterrows():
#     nt.add_edge(
#         data['source_taxon_id'],
#         data['target_taxon_id'],
#         label=data['interaction_type'])

# # nt.show_buttons(filter_=['physics'])
# nt.show("tmp.fig02.html")


# In[15]:


# visualisasi embedding serangga

# embedding
# input : df_node, URL
URL = "http://pyRDF2Vec"
df_serangga = df_node[df_node['group']=="serangga"]
CUSTOM_KG = df_serangga_to_rdf(df_serangga, URL)
# proses
list_serangga=df_node[df_node['group']=="serangga"].taxon_id.to_list()
list_of_entities = [ URL+"#"+taxon_id for taxon_id in list_serangga ]
transformer, embeddings, _ = rdf_KG_to_embeddings(CUSTOM_KG, list_of_entities)
# sama dua2nya # list_of_entities == transformer._entities 

# dictionary serangga
dictionary_serangga = df_to_dictionary_taxon(df_serangga)
# output
# embeddings, list_of_entities, dictionary_serangga

# plotting
embeddingPlot(embeddings, list_of_entities, dictionary_serangga)


# # Analisis Interaksi

# In[16]:


# from modul.custom_degree_centrality import degree_centrality_custom
# import importlib, sys
# importlib.reload(sys.modules['modul.custom_degree_centrality'])

#4 
# Degree Centrality Custom
virus_utama_ids=list(df_node[df_node['virus_utama']==True].taxon_id)
serangga_ids=list(df_node[df_node['group']=="serangga"].taxon_id)
results_dc = degree_centrality_custom(gnx,virus_utama_ids,serangga_ids)
allnodes = gnx.nodes


# In[17]:


urutan=sorted(results_dc.items(), key=lambda item: item[1], reverse=True)
urutan


# In[18]:


#visualisasi data
# mengecek BFS Degree tertinggi
for to_search  in [urutan[0][0]]:#["NCBI:7038"]:#"NCBI:33377","NCBI:7036",
    for edge in nx.bfs_edges(gnx.to_undirected(), source=to_search, depth_limit=1):
            s_id, o_id = edge

            s_label = gnx.nodes[s_id]['label'] +' '+s_id
            o_label = gnx.nodes[o_id]['label'] +' '+o_id
            o_grup = gnx.nodes[o_id]['group']

            # skip subjek to_search
            if (o_grup == 'virus'):
                print(s_label,'-->', o_label)
    print('===============')


# # Analisis Taksonomi

# In[19]:


#5
# Ambil data NCBI
# data acuan
data_acuan=get_taxon_vector(acuan_,ncbi_ontology_url,False)
print(data_acuan)
data_ujian=get_taxon_vector(ujian_,ncbi_ontology_url)
print(data_ujian)


# In[20]:


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



# visualisasi
# input
# embeddings, list_of_entities, dictionary_serangga
# embeddingPlot(embeddings, list_of_entities, dictionary_serangga)


# In[21]:


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
    data_to_count.loc[idx,['dc_result']] = results_dc[row['entity']] if(row['entity']!="SERANGGA_ACUAN") else 1

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


# In[22]:


#drop kolom embedding
data_to_count.drop(columns=list(range(0,100)), inplace=True)


# In[23]:





# simple scaling ed_result
# data_to_count["ed_result_scaled"] = data_to_count["ed_result"] / data_to_count["ed_result"].max()
data_to_count['ed_result_scaled'] = minmax(data_to_count['ed_result'])


# scaling dc 
# # dengan std
# data_to_count['dc_result_scaled'] = data_to_count['dc_result'] - data_to_count['dc_result'].mode()[0]
# # data_to_count['dc_result_scaled'] = data_to_count['dc_result']
# # lagi dengan minmax
# data_to_count["dc_result_scaled"] = data_to_count["dc_result_scaled"] / data_to_count["dc_result_scaled"].max()

# Perform scaling using standard deviation
# data_to_count['dc_result_scaled'] = std_scale(data_to_count['dc_result'])
# data_to_count['dc_result_scaled'] = minmax(data_to_count['dc_result_scaled'])

data_to_count['dc_result_scaled'] = minmax(data_to_count['dc_result'])


# In[24]:


data_to_count


# In[25]:


#9
#hitung kombinasi
# for idx, row in data_to_count.iterrows():
#     _dc = row['dc_result_scaled']
#     _ed=( (row['ed_result_scaled']) if row['ed_result_scaled']!=0 else 1)
#     data_to_count.loc[idx,['result']] = _dc/_ed

data_to_count['result'] = (1+data_to_count['dc_result_scaled']) / (1+data_to_count['ed_result_scaled'])

# simple scaling result (final/kombinasi)
data_to_count['result'] = data_to_count['result'] / data_to_count['result'].max()


# # Pengujian

# In[26]:


# DC
data_to_count=data_to_count.sort_values('dc_result',ascending=False).reset_index(drop=True)
data_to_count[['label','entity','dc_result']]


# In[27]:


# Pengujian dc
for urutan in range(0,3):
    takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
    id_hasil=data_to_count.iloc[urutan].entity
    cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
    cek_ujian= { k:v for k,v in data_ujian if k in takson }
    # print(acuan_,'->', data_)
    # print('ujian ',cek_ujian)
    print('hasil ',cek_hasil)
    cek=0
    for i in reversed(takson):
        cekk=cek_hasil[i]==cek_ujian[i]
        cek+=cekk
        # print(i, cekk)
    print(cek/len(takson))


# In[28]:


# ED
data_to_count=data_to_count.sort_values('ed_result',ascending=True).reset_index(drop=True)
data_to_count[['label','entity','ed_result']]


# In[29]:


# Pengujian ed
for urutan in range(0,3):
    takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
    id_hasil=data_to_count.iloc[urutan].entity
    cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
    cek_ujian= { k:v for k,v in data_ujian if k in takson }
    # print(acuan_,'->', data_)
    # print('ujian ',cek_ujian)
    print('hasil ',cek_hasil)
    cek=0
    for i in reversed(takson):
        cekk=cek_hasil[i]==cek_ujian[i]
        cek+=cekk
        # print(i, cekk)
    print(cek/len(takson))


# In[30]:


# final score kombinasi
data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)
data_to_count[['label','dc_result','ed_result','result']]


# In[31]:


# Pengujian kombinasi
for urutan in range(0,3):
    takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
    id_hasil=data_to_count.iloc[urutan].entity
    cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
    cek_ujian= { k:v for k,v in data_ujian if k in takson }
    # print(acuan_,'->', data_)
    # print('ujian ',cek_ujian)
    print('hasil ',cek_hasil)
    cek=0
    for i in reversed(takson):
        cekk=cek_hasil[i]==cek_ujian[i]
        cek+=cekk
        # print(i, cekk)
    print(cek/len(takson))


# In[32]:


# Pengujian kombinasi
# cek=1
# takson=[i[0] for i in data_ujian if i[0] not in ["superkingdom","kingdom","filum","kelas"]]
# urutan=0
# while cek/len(takson)>0:
#     id_hasil=data_to_count.iloc[urutan].entity
#     cek_hasil= { k:v for k,v in reversed(allnodes[id_hasil].items()) if k in takson }
#     cek_ujian= { k:v for k,v in data_ujian if k in takson }
#     # print(acuan_,'->', data_)
#     # print('ujian ',cek_ujian)
#     print('hasil ', urutan ,cek_hasil)
#     cek=0
#     for i in reversed(takson):
#         cekk=cek_hasil[i]==cek_ujian[i]
#         cek+=cekk
#         # print(i, cekk)
#     print(cek/len(takson))
#     urutan+=1


# # Visualisasi Hasil Analisis

# In[33]:


import plotly.graph_objects as go

to_itter=data_to_count.sort_values('dc_result',ascending=False).reset_index(drop=True)[:10]

label=[]
degree=[]
for index,data in to_itter.sort_values('dc_result',ascending=True).reset_index(drop=True).iterrows():
    label.append(data.label)
    degree.append(data.dc_result)

fig = go.Figure(data=go.Bar(
    x=degree,
    y=label,
    orientation='h'
))

fig.update_layout(
    title='Degree centrality',
    xaxis_title='Degree',
    yaxis_title='Insect'
)

fig.show()


# In[34]:


to_itter=data_to_count.sort_values('ed_result',ascending=True).reset_index(drop=True)[:10]

label=[]
degree=[]
for index,data in to_itter.sort_values('ed_result',ascending=False).reset_index(drop=True).iterrows():
    label.append(data.label)
    degree.append(data.ed_result)

fig = go.Figure(data=go.Bar(
    x=degree,
    y=label,
    orientation='h'
))

fig.update_layout(
    title='Euclidean distance',
    xaxis_title='Distance',
    yaxis_title='Insect'
)

fig.show()


# In[35]:


import plotly.graph_objects as go

to_itter=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)[:10]

label=[]
degree=[]
for index,data in to_itter.sort_values('result',ascending=True).reset_index(drop=True).iterrows():
    label.append(data.label)
    degree.append(data.result)

fig = go.Figure(data=go.Bar(
    x=degree,
    y=label,
    orientation='h'
))

fig.update_layout(
    title='final score',
    xaxis_title='Score',
    yaxis_title='Insect'
)

fig.show()


# # Dibawah ini tidak termasuk

# In[36]:


# #visualisasi data
# # mengecek BFS Degree tertinggi
# for edge in nx.bfs_edges(gnx.to_undirected(), source="NCBI:7038", depth_limit=1):
#         s_id, o_id = edge

#         s_label = gnx.nodes[s_id]['label'] +' '+s_id
#         o_label = gnx.nodes[o_id]['label'] +' '+o_id
#         o_grup = gnx.nodes[o_id]['group']

#         # skip subjek to_search

#         print(s_label,'-->', o_label)
# print('===============')


# In[37]:


df_node[df_node['taxon_name'].str.contains('Thrips')]


# In[38]:


df_edge.interaction_type.unique()


# In[39]:


ini=df_edge[df_edge['interaction_type'].isin(['visitFlowersOf'])]
ini


# In[40]:


df_node[df_node['taxon_id'].isin(
    ini.source_taxon_id.to_list()+
    ini.target_taxon_id.to_list()
)]


# In[41]:


contains_string_entire_column(df_edge,'NCBI:1341303')


# In[42]:


contains_string_entire_column(df_node,'NCBI:12315')


# In[ ]:




