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

# In[299]:


#1
#baca data
df_node=pd.read_csv('dari_praproses/'+data_+'_node.csv',index_col=0) 
df_edge=pd.read_csv('dari_praproses/'+data_+'_edge.csv',index_col=0)


# In[300]:


df_node[df_node['group']=='virus']
# df_node['group'].unique()

# df_node


# In[301]:


acuan_


# In[302]:


# pra-proses khusus proses
# hapus serangga yg cuma famili (mengikuti acuan). soalnya klo cuma tampil famili apa gunanya?
filter_genus_sampai_species_null=(
    (df_node.genus.isnull()) &
    (df_node.species.isnull()) &
    (df_node['class']=='NCBI:50557_Insecta')
)
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_genus_sampai_species_null], df_node,df_edge)


# In[303]:


print(len(df_edge))
df_edge.drop_duplicates(inplace=True)
print(len(df_edge))


# In[304]:


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


# In[305]:


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


# In[306]:


# Ini harusnya di praproses. tapi belum fix baiknya hapus atau tidak
# hapus yang bukan virus utama, terakhir akurasi 0.85, kalau berkurang hapus saja ini
bukan_virus_utama=(df_node['group']=="virus") & (df_node.virus_utama!=True)
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[bukan_virus_utama], df_node,df_edge)


# In[307]:


if(len(df_node[df_node['group']=="serangga"])<=2):
    print("cuma dua serangga")


# # Konversi graf

# In[308]:


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


# In[309]:


# # cuma tampilan, visualisasi graf
# G=gnx

# fig, ax = plt.subplots(figsize=(20, 20))

# # Generate layout for visualization
# # pos = nx.kamada_kawai_layout(G)
# # pos = nx.spring_layout(G)
# pos = nx.nx_agraph.graphviz_layout(G, prog="neato", args="")

# # Visualize graph components
# nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='g')
# nx.draw_networkx_nodes(G, pos, node_color=list(nx.get_node_attributes(G, "color").values()), alpha=0.9)

# #node label
# # for i in ['#b22222','#671f92','#1f922b','#EADDCA']: # filtering dengan bedakan warna node
# #     label_options = {"ec": i, "fc": 'white', "alpha": 0.7}
# #     nx.draw_networkx_labels(
# #         nx.subgraph_view(G, filter_node=lambda n1: G.nodes(data=True)[n1].get("color", True) == i),
# #         pos, 
# #         font_size=10, 
# #         bbox=label_options
# #     )

# #edge labels
# edge_labels={x:i for i,x in zip(nx.get_edge_attributes(G, "label").values(),G.edges())}
# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)


# # Title/legend
# font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
# ax.set_title("Graf Interaksi "+nama_virus, font)
# # Change font color for legend
# font["color"] = "r"

# ax.text(
#     0.80,
#     0.10,
#     "hijau = Tanaman",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )
# ax.text(
#     0.80,
#     0.08,
#     "merah = Serangga",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )

# ax.text(
#     0.80,
#     0.06,
#     "ungu = Virus",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )

# ax.text(
#     0.80,
#     0.04,
#     "abu-abu = Nogroup",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )

# # Resize figure for label readibility
# ax.margins(0.1, 0.05)
# fig.tight_layout()
# plt.axis("off")
# plt.show()


# In[285]:


import plotly.graph_objects as go
G=gnx
pos = nx.nx_agraph.graphviz_layout(G, prog="neato", args="")

edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line= {"width":0.5, "color":'#888'},
    hoverinfo='none',
    mode='lines')

node_x = []
node_y = []
node_colors = []
node_text = []
for node,data in G.nodes(data=True):
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_colors.append(data['color'])
    node_text.append(data['label'])

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker={
        # 'showscale':True,
        # 'colorscale':'Reds',
        'reversescale':True,
        'color':[],
        'size':10,
        # 'colorbar':{
        #     # 'thickness':15,
        #     # 'title':'Node Connections',
        #     # 'xanchor':'left',
        #     # 'titleside':'right'
        # },
        'line_width':2   
    }
)
node_trace.marker.color = node_colors
node_trace.text = node_text


# In[293]:


fig = go.Figure(
    data=[edge_trace, node_trace],
    layout=go.Layout(
        title='Network graph made with Python',
        titlefont_size=16,
        showlegend=False,
        hovermode='closest',
        margin={
            'b':20,'l':5,'r':5,'t':40
        },
        annotations=[{
            "text":"Insect-virus-plant",
            'showarrow':False,
            'xref':"paper", 
            'yref':"paper",
            'x':0.005, 
            'y':-0.002 
        }],
        xaxis={'showgrid':False, 'zeroline':False, 'showticklabels':False},
        yaxis={'showgrid':False, 'zeroline':False, 'showticklabels':False}
    )
)
fig.show()


# In[261]:


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


# # Analisis Interaksi

# In[138]:


from modul.custom_degree_centrality import degree_centrality_custom
import importlib, sys
importlib.reload(sys.modules['modul.custom_degree_centrality'])

#4 
# Degree Centrality Custom
virus_utama_ids=list(df_node[df_node['virus_utama']==True].taxon_id)
serangga_ids=list(df_node[df_node['group']=="serangga"].taxon_id)
results_dc = degree_centrality_custom(gnx,virus_utama_ids,serangga_ids)
allnodes = gnx.nodes


# In[139]:


urutan=sorted(results_dc.items(), key=lambda item: item[1], reverse=True)
urutan


# In[140]:


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

# In[141]:


#5
# Ambil data NCBI
# data acuan
data_acuan=get_taxon_vector(acuan_,ncbi_ontology_url)
print(data_acuan)
data_ujian=get_taxon_vector(ujian_,ncbi_ontology_url)
print(data_ujian)


# In[142]:


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


# In[143]:


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


# In[144]:


# visualisasi

# Reduce the dimensions of entity embeddings to represent them in a 2D plane.
X= UMAP().fit_transform(embeddings)
df_umap=pd.DataFrame(X,columns=['feature-vector-1','feature-vector-2'])


text=[]
labels=[]
for x in transformer._entities:
    if(x!="http://pyRDF2Vec#SERANGGA_ACUAN"):
        text.append(gnx.nodes[x.split("#")[-1]]['famili'].split('_')[-1])
        labels.append(gnx.nodes[x.split("#")[-1]]['label'])
    else:
        text.append("#TITIK_VEKTOR_ACUAN")
        labels.append("#TITIK_VEKTOR_ACUAN")
df_umap['text']=text
df_umap['labels']=labels

# # gnx.nodes[x.split("#")[-1]]['label']
# df_umap['text']=list(map(lambda x: x.split("#")[-1],transformer._entities))
fig = px.scatter(df_umap, x='feature-vector-1',y='feature-vector-2',text='text',hover_name='labels')
fig.update_traces(textposition='top center')
fig.update_layout(
    height=650,
    title_text='reduced word2vec visualization'
)
fig.show()


# In[145]:


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


# In[146]:


#drop kolom embedding
data_to_count.drop(columns=list(range(0,100)), inplace=True)


# In[147]:


def minmax(data):
    return (data - data.min())/ (data.max() - data.min())

def std_scale(data):
    return (data - data.mean()) / data.std()


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


# In[148]:


data_to_count


# In[149]:


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

# In[150]:


# DC
data_to_count=data_to_count.sort_values('dc_result',ascending=False).reset_index(drop=True)
data_to_count[['label','entity','dc_result']]


# In[151]:


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


# In[152]:


# ED
data_to_count=data_to_count.sort_values('ed_result',ascending=True).reset_index(drop=True)
data_to_count[['label','entity','ed_result']]


# In[153]:


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


# In[154]:


# final score kombinasi
data_to_count=data_to_count.sort_values('result',ascending=False).reset_index(drop=True)
data_to_count[['label','dc_result','ed_result','result']]


# In[155]:


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


# In[156]:


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

# In[157]:


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


# In[158]:


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


# In[159]:


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

# In[160]:


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


# In[161]:


df_node[df_node['taxon_name'].str.contains('Thrips')]


# In[162]:


df_edge.interaction_type.unique()


# In[163]:


ini=df_edge[df_edge['interaction_type'].isin(['visitFlowersOf'])]
ini


# In[164]:


df_node[df_node['taxon_id'].isin(
    ini.source_taxon_id.to_list()+
    ini.target_taxon_id.to_list()
)]


# In[165]:


contains_string_entire_column(df_edge,'NCBI:1341303')


# In[166]:


contains_string_entire_column(df_node,'NCBI:12315')


# In[ ]:




