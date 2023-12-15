#!/usr/bin/env python
# coding: utf-8

# In[44]:


# import importlib, sys
# importlib.reload(sys.modules['modul.standardization_usingsparql'])
import requests
from tqdm import tqdm
import pandas as pd
# from owlready2 import get_ontology
import os


from modul.standardization_usingsparql import addTaxonColumn, buat_kolom_taxon_awal
from modul.disambiguation_optimized import buat_kamus_kosong, update_kamus_pake_wikidata, update_df_pake_kamus, update_df_pake_path_ujung, removeOtherThanNCBI, __chunk_list
from modul.preprocess import cleaning, splitInteractionToNodeEdge, pagination_search_globi, test_pagination_link
from modul.filterNodeEdge import removeNodeAndEdgeByFilter,takeNodeAndEdgeByFilter,removeEdgesNotInNodes
from modul.helper_umum import contains_string_entire_column,contains_string_entire_column_boolean, pemecah_generator
from modul.vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id_via_string


# In[45]:


data_init=[
    ('begomovirus_contoh_hasil_baru','Pepper yellow leaf curl virus','Aleyrodidae','Bemisia Tabaci'),
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


# In[46]:


#0 inisiasi parameter
ini_data=data_init[0]
# parameter
nama_file = ini_data[0]
virus_txt = ini_data[1].replace(' ','%20')
tipe_interaksi_virus = 'hasHost' #pathogenOf, pake relasi hasHost lebih dapat banyak relasi dari pada pathogenOf
tipe_interaksi_tanaman = 'hostOf' #hasPathogen, pake relasi hostOf lebih dapat banyak relasi dari pada hasPathogen
tipe_interaksi_serangga_ke_tanaman = 'hasHost' 
tipe_interaksi_serangga_ke_virus = 'hostOf' 
ncbi_server_url = 'http://localhost:3030/mydataset/query'
offset_limit_pertama=30
offset_limit_kedua=3


# In[47]:


virus_search = get_taxon_vector(virus_txt,ncbi_server_url)
if (virus_search==False):
    print('virus tidak ditemukan')
    # virus_search=[('unknown',virus_txt)]
virus_search


# In[48]:


#1 BFS Data 
# interaksi virus --pathogenOf-> tanaman dan serangga
kolom=[
    'source_taxon_external_id',
    'source_taxon_name',
    'source_taxon_path',
    'source_taxon_path_ids',
    'source_taxon_path_ranks',
    
    'interaction_type',
    
    'target_taxon_external_id',
    'target_taxon_name',
    'target_taxon_path',
    'target_taxon_path_ids',
    'target_taxon_path_ranks',
]
interactionType=tipe_interaksi_virus

# inisiasi dataframe
df_init=pd.DataFrame(columns = kolom)

# list pencarian
if (virus_search!=False):
    list_source_taxon_virus = []
    for i in range(len(virus_search)):
        if (
            virus_search[i][0] in ['famili','genus','spesies'] and 
            len(list_source_taxon_virus) < 2 # maksimal 2 pencarian
        ):
            search = virus_search[i][1].split('_')[0]
            list_source_taxon_virus.append(search)
    text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(list_source_taxon_virus)
else:
    # jika virus_search tidak ditemukan pencarian menggunakan text saja
    text_source_taxon = "sourceTaxon=" + virus_txt
print(text_source_taxon)

# pencarian data
link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+"&targetTaxon=Viridiplantae&targetTaxon=Insecta"+"&fields="+(','.join(kolom))
print(link)
df = pemecah_generator(pagination_search_globi(link, df_init, offset_limit_pertama))

if len(df)==0:
    print('tidak ada data')
    # hentikan disini.
    # exit()


# In[49]:


#2 splitting layer 1 interaksi virus
df_node, df_edge = splitInteractionToNodeEdge(df)


# In[50]:


# cleaning_after_split layer 1 interaksi virus
# drop duplikat
df_node.drop_duplicates(inplace=True)
# no_ncbi dan path_null
no_ncbi_and_path_null=(df_node.taxon_id.str.contains('NCBI')==False) & (df_node.taxon_path_ids.isnull())
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[no_ncbi_and_path_null], df_node,df_edge) 
# drop duplikat
df_edge.drop_duplicates(inplace=True)


# In[51]:


# tandai virus utama
filter_virus_utama=(
    (df_node.taxon_name.str.contains(r'\b(virus\w*|\w*virus)\b',case=False))
    | (df_node.taxon_path.str.contains(r'\b(virus\w*|\w*virus)\b', case=False))  
    #jika berawalan atau berakhiran kata virus
)
# df_node.loc[filter_virus_utama, ['virus_utama']] = True
# virus_utama=[data.taxon_id for idx,data in df_node[filter_virus_utama].iterrows()]
virus_utama = df_node[filter_virus_utama].taxon_id.unique().tolist()


# In[9]:


#3 disambiguasi layer 1 interaksi virus
kamus_ncbi = buat_kamus_kosong(df_node)
kamus_ncbi = pemecah_generator(update_kamus_pake_wikidata(kamus_ncbi))
#update dataframe pake kamus
df_node,df_edge = pemecah_generator(update_df_pake_kamus(kamus_ncbi,df_node,df_edge))
df_node,df_edge = pemecah_generator(update_df_pake_path_ujung(df_node,df_edge))
#tambah kolom takson pake data NCBI
df_node = buat_kolom_taxon_awal(df_node) #buat kolom taxon, default none
df_node = addTaxonColumn(df_node,'http://localhost:3030/mydataset/query') # isi pake ncbi


# In[10]:


# cleaning_after_disambiguasi layer 1
df_node, df_edge = removeOtherThanNCBI(df_node,df_edge)# Hapus kalo masih ada selain NCBI
df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus edge yang tidak ada nodenya
filter_kingdom_atau_class_null=( (df_node.kingdom.isnull()) | (df_node['class'].isnull()) )
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_kingdom_atau_class_null], df_node,df_edge)


# In[11]:


# inisiasi
df_to_add=pd.DataFrame(columns = kolom)

#4.1 BFS interaksi tanaman --hostOf-> serangga dan virus
df_plant=df_node[df_node.kingdom=='NCBI:33090_Viridiplantae']
if not df_plant.empty:
    interactionType = tipe_interaksi_tanaman
    list_source_taxon = df_plant.taxon_id.unique()
    list_source_taxon = list(set(list_source_taxon)) #unique
    # list_target_taxon = list_source_taxon_virus + ['Insecta']
    list_target_taxon = ['Insecta','Viruses']

    text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(list_source_taxon)
    text_target_taxon = "&targetTaxon=" + "&targetTaxon=".join(list_target_taxon)
    # pencarian data
    link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI"
    print(link)

    # test dan eksekusi
    test_link = pemecah_generator(test_pagination_link(link))
    print(test_link)
    if test_link=='aman':
        df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='pecah':
        print('perlu dipecah')
        for i in __chunk_list(list_source_taxon, 1):
            text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(i)
            link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI"
            df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='skip':
        print('skip, link error')


# In[12]:


#4.2 BFS interaksi serangga --pathogenOf-> tanaman
df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
if not df_insect.empty:
    interactionType = tipe_interaksi_serangga_ke_tanaman
    list_source_taxon = df_insect.taxon_id.unique()
    list_source_taxon = list(set(list_source_taxon)) #unique
    list_target_taxon = ['Viridiplantae']

    text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(list_source_taxon)
    text_target_taxon = "&targetTaxon=" + "&targetTaxon=".join(list_target_taxon)
    # pencarian data
    link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI" # klo tda bagus hapus NCBI prefix itu
    print(link)

    # test dan eksekusi
    test_link = pemecah_generator(test_pagination_link(link))
    print(test_link)
    if test_link=='aman':
        df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='pecah':
        print('perlu dipecah')
        for i in __chunk_list(list_source_taxon, 1):
            text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(i)
            link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI"
            df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='skip':
        print('skip, link error')


# In[13]:


#4.3 BFS interaksi serangga --hostOf-> virus
df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
if not df_insect.empty:
    interactionType = tipe_interaksi_serangga_ke_virus
    list_source_taxon = df_insect.taxon_id.unique()
    list_source_taxon = list(set(list_source_taxon)) #unique
    list_target_taxon = ['Viruses']

    text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(list_source_taxon)
    text_target_taxon = "&targetTaxon=" + "&targetTaxon=".join(list_target_taxon)
    # pencarian data
    link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI" # klo tda bagus hapus NCBI prefix itu
    print(link)

    # test dan eksekusi
    test_link = pemecah_generator(test_pagination_link(link))
    print(test_link)
    if test_link=='aman':
        df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='pecah':
        print('perlu dipecah')
        for i in __chunk_list(list_source_taxon, 1):
            text_source_taxon = "sourceTaxon=" + "&sourceTaxon=".join(i)
            link="https://api.globalbioticinteractions.org/interaction?"+text_source_taxon+"&interactionType="+interactionType+text_target_taxon+"&fields="+(','.join(kolom))+"taxonIdPrefix=NCBI"
            df_to_add = pemecah_generator(pagination_search_globi(link, df_to_add, offset_limit_kedua))
    elif test_link=='skip':
        print('skip, link error')


# In[14]:


#5 splitting depth 2 interaksi serangga dan tanaman
node_to_add, edge_to_add = splitInteractionToNodeEdge(df_to_add)


# In[15]:


df_to_add.shape


# In[16]:


# cleaning_after_split depth 2 interaksi serangga dan tanaman

# hapus edge inverse
kebalikan={
    'hostOf':'hasHost',
    'hasPathogen':'pathogenOf', 
    'pollinatedBy':'pollinates', 
    'flowersVisitedBy':'visitFlowersOf',
    'visitedBy':'visit'
}
for i,data in edge_to_add[edge_to_add['interaction_type'].isin([
    'hostOf', 'hasPathogen', 'pollinatedBy', 'flowersVisitedBy','visitedBy'
    ])].iterrows():
    edge_to_add.iloc[i]['interaction_type']=kebalikan[data['interaction_type']]
    edge_to_add.iloc[i]['source_taxon_id']=data['target_taxon_id']
    edge_to_add.iloc[i]['target_taxon_id']=data['source_taxon_id']

# drop duplikat
print('node_to_add.drop_duplicates')
node_to_add.drop_duplicates(inplace=True)
print(len(node_to_add),len(edge_to_add))

# hapus no_ncbi_and_path_null
no_ncbi_and_path_null=(node_to_add.taxon_id.str.contains('NCBI')==False) & (node_to_add.taxon_path_ids.isnull())
node_to_add,edge_to_add = removeNodeAndEdgeByFilter(node_to_add[no_ncbi_and_path_null], node_to_add,edge_to_add) 

# hapus edge duplikat
print('edge_to_add.drop_duplicates')
edge_to_add.drop_duplicates(inplace=True)
print(len(node_to_add),len(edge_to_add))


# In[17]:


# 6 disambiguasi layer 2
kamus_ncbi = buat_kamus_kosong(node_to_add)
kamus_ncbi = pemecah_generator(update_kamus_pake_wikidata(kamus_ncbi))
#update dataframe pake kamus
node_to_add,edge_to_add = pemecah_generator(update_df_pake_kamus(kamus_ncbi,node_to_add,edge_to_add))
node_to_add,edge_to_add = pemecah_generator(update_df_pake_path_ujung(node_to_add, edge_to_add))
# tambah kolom takson pake data NCBI
node_to_add = buat_kolom_taxon_awal(node_to_add) #buat kolom taxon, isi none dan isi dari path
node_to_add = addTaxonColumn(node_to_add,'http://localhost:3030/mydataset/query') #isi kolom taxon, pake NCBI


# In[18]:


#untuk laporan
#kode database
print(kamus_ncbi.keys())

key=[]
val=[]
for i in kamus_ncbi:
    key.extend(list(kamus_ncbi[i].keys()))
    val.extend(list(kamus_ncbi[i].values()))
    
df_kamus=pd.DataFrame({'key':key,'val':val})

# csv
# df_kamus.to_csv('output.xlsx', index=False)

# semua 
print(df_kamus.count())
# yang kosong
# df_kamus[(df_kamus.val != '') & (df_kamus.key.str.contains("NBN"))]
print(df_kamus[(df_kamus.val != '')].count())


# In[19]:


# cleaning_after_disambiguasi depth 2
node_to_add,edge_to_add = removeOtherThanNCBI(node_to_add,edge_to_add) #hapus selain NCBI


# In[20]:


#7 konkatenasi tabel
df_node=pd.concat([df_node,node_to_add], axis=0, ignore_index=True)
df_edge=pd.concat([df_edge,edge_to_add], axis=0, ignore_index=True)


# In[21]:


# praproses tambahan

# hapus duplikat
df_node.drop_duplicates(inplace=True)
df_edge.drop_duplicates(inplace=True)

# pengelompokan
# Binning of the data based on serangga, virus, tanaman, nogroup
filter_tanaman = df_node['kingdom']=='NCBI:33090_Viridiplantae' 
filter_virus = (
    (df_node['superkingdom']=='NCBI:10239_Viruses')
    | (df_node.taxon_name.str.contains(r'\b(virus\w*|\w*virus)\b',case=False))
    | (df_node.taxon_path.str.contains(r'\b(virus\w*|\w*virus)\b', case=False)) 
    #jika berawalan atau berakhiran kata virus
)
filter_serangga = ((df_node['class']=='NCBI:50557_Insecta') )#& (df_node['order']!='NCBI:7399_Hymenoptera')) #dan bukan lebah hymenoptera

df_node.loc[filter_tanaman, ['group','color']] = ["tanaman",'#1f922b'] #hijau
df_node.loc[filter_virus, ['group','color']] = ['virus','#671f92'] #ungu
df_node.loc[filter_serangga, ['group','color']] = ['serangga','#b22222'] #merah
df_node.loc[(
    (filter_tanaman==False) & 
    (filter_virus==False) &
    (filter_serangga==False) 
    ),['group','color']] = ['nogroup','#EADDCA'] #abu-abu


# In[22]:


# cleaning setelah pengelompokan, sebelum konversi graf
# hapus node yang ordo sampai specie isi null
filter_ordo_sampai_species_null=(
    (df_node.order.isnull()) & 
    (df_node.family.isnull()) & 
    (df_node.genus.isnull()) &
    (df_node.species.isnull())
)
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_ordo_sampai_species_null], df_node,df_edge)

# hapus duplikasi
df_node.drop_duplicates(subset=["taxon_id"], keep='last',inplace=True)
df_edge = removeEdgesNotInNodes(df_node, df_edge)#edge menyesuaikan

# hapus kingdom isi null
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[(df_node.kingdom.isnull()) & (df_node.group!='virus')], df_node,df_edge) 
df_edge = removeEdgesNotInNodes(df_node, df_edge) #edge menyesuaikan #cuma memastikan saja

# 13 # hapus yang no group
df_node,df_edge = removeNodeAndEdgeByFilter(df_node[df_node.group=="nogroup"], df_node,df_edge) 
df_edge = removeEdgesNotInNodes(df_node, df_edge) #edge menyesuaikan #cuma memastikan saja

#14 # hapus node yang tidak punya edge
print('hapus node yang tidak ada di edge (tidak punya edge)')
print('sebelum',len(df_node))
df_node = df_node[
    (df_node.taxon_id.isin(df_edge.source_taxon_id.unique())) 
    | (df_node.taxon_id.isin(df_edge.target_taxon_id.unique()))
]
print('sesudah',len(df_node))

#reset index
df_node.reset_index(drop=True)

# masukan tanda virus utama
df_node.loc[df_node.taxon_id.isin(virus_utama), ['virus_utama']] = True

print(df_node.shape, df_edge.shape)


# In[23]:


# akhir pra proses
# save file
df_edge.reset_index(drop=True,inplace=True)
df_node.reset_index(drop=True,inplace=True)
df_node.to_csv(os.getcwd()+'/dari_praproses/'+nama_file+'_node.csv')
df_edge.to_csv(os.getcwd()+'/dari_praproses/'+nama_file+'_edge.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# # Dibawah ini tidak masuk pra proses

# In[24]:


import pandas as pd


# In[25]:


df_node=pd.read_csv('dari_praproses/'+nama_file+'_node.csv',index_col=0) 
df_edge=pd.read_csv('dari_praproses/'+nama_file+'_edge.csv',index_col=0)


# ## proporsi

# In[26]:


# cuma tampilan
import plotly.graph_objects as go

data = df_node.groupby(['group','color']).agg({'group': ['count'], }).reset_index().sort_values(
    ('group', 'count'),ascending=False
).reset_index(drop=True).values
labels = [i[0] for i in data]
colors = [i[1] for i in data]
slices = [i[2] for i in data]

fig = go.Figure(data=[go.Pie(labels=labels,values=slices)])
fig.update_traces(hoverinfo='label+percent', textinfo='value+percent', textfont_size=20, marker=dict(colors=colors, line=dict(color='#000000', width=0.1)))
fig.show()


# In[27]:


slices,labels


# In[28]:


#cek


# In[29]:


import networkx as nx
import matplotlib.pyplot as plt

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


# In[30]:


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
# ax.set_title("Interaksi Tanaman-Serangga-Virus", font)
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


# In[31]:


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


# In[32]:


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


# In[33]:


def sub_generator():
    for i in range(1, 4):
        yield i
    return "Data dari sub-generator"

def main_generator():
    yield "Awal"
    result = yield from sub_generator()
    # yield result  # Menggunakan hasil yang dikembalikan dari sub-generator
    print('ini kembalian', result) 
    yield "Akhir"

# Menggunakan main_generator
for item in main_generator():
    print(item)


# In[34]:


def sub_generator():
    for i in range(1, 4):
        yield i
    return "Data dari sub-generator"

def middle_generator():
    result = yield from sub_generator()
    yield "Data dari middle_generator"
    return "kembalian middle_generator dan " + result

def main_generator():
    yield "Awal"
    kembalian = yield from middle_generator()
    print('ini kembalian', kembalian)
    yield "Akhir"

# Menggunakan main_generator
for item in main_generator():
    print(item)


# In[ ]:





# In[ ]:





# In[35]:


link= "https://api.globalbioticinteractions.org/interaction?sourceTaxon=Alternanthera%20sessilis&sourceTaxon=Gossypium%20raimondii&sourceTaxon=Ruellia%20blechum&sourceTaxon=Scrophulariaceae&sourceTaxon=Brassica%20rapa&sourceTaxon=Spinacia%20oleracea&sourceTaxon=Jatropha%20integerrima&sourceTaxon=Solanum%20lycopersicum&sourceTaxon=Capsicum%20annuum&sourceTaxon=Lablab%20purpureus&sourceTaxon=Physalis%20pubescens&sourceTaxon=Physalis%20ixocarpa&sourceTaxon=Jacquemontia&sourceTaxon=Phaseolus%20vulgaris&sourceTaxon=Solanum%20aculeatissimum&sourceTaxon=Triumfetta&sourceTaxon=Rosa&sourceTaxon=Sauropus%20androgynus&sourceTaxon=Sonchus%20oleraceus&sourceTaxon=Clerodendrum%20philippinum&sourceTaxon=Acalypha%20indica&sourceTaxon=Passiflora%20edulis%20f.%20flavicarpa&sourceTaxon=Corchorus%20siliquosus&sourceTaxon=Sida&sourceTaxon=Clerodendrum%20cyrtophyllum&sourceTaxon=Capsicum%20frutescens&sourceTaxon=Solanum%20tuberosum&sourceTaxon=Macroptilium&sourceTaxon=Gossypium%20darwinii&sourceTaxon=Ludwigia%20octovalvis&sourceTaxon=Solanum%20nigrum&sourceTaxon=Manihot%20glaziovii&sourceTaxon=Pueraria%20montana&sourceTaxon=Nicotiana&sourceTaxon=Jatropha%20curcas&sourceTaxon=Ipomoea%20tiliacea&sourceTaxon=Croton%20glandulosus&sourceTaxon=Crassocephalum%20crepidioides&sourceTaxon=Luffa&sourceTaxon=Ipomoea%20purpurea&sourceTaxon=Cucurbita%20moschata&sourceTaxon=Phaseolus%20acutifolius&sourceTaxon=Desmodium%20glabrum&sourceTaxon=Stachytarpheta&sourceTaxon=Momordica%20charantia&sourceTaxon=Nicotiana%20debneyi&sourceTaxon=Ocimum%20gratissimum&sourceTaxon=Carica%20papaya&sourceTaxon=Ipomoea%20trifida&sourceTaxon=Ipomoea%20indica&sourceTaxon=Dicliptera&sourceTaxon=Zinnia%20elegans&sourceTaxon=Malva%20parviflora&sourceTaxon=Sechium%20edule&sourceTaxon=Senecio%20scandens&sourceTaxon=Pueraria%20montana%20var.%20lobata&sourceTaxon=Nicandra%20physalodes&sourceTaxon=Macrotyloma%20uniflorum&sourceTaxon=Phaseolus&sourceTaxon=Dysphania%20ambrosioides&sourceTaxon=Euphorbia%20tithymaloides&sourceTaxon=Ipomoea%20cordatotriloba&sourceTaxon=Clerodendrum&sourceTaxon=Gossypium%20stocksii&sourceTaxon=Cynanchum%20acutum&sourceTaxon=Asystasia%20gangetica&sourceTaxon=Dalechampia&sourceTaxon=Nicotiana%20glutinosa&sourceTaxon=Gossypium%20somalense&sourceTaxon=Vigna%20radiata%20var.%20radiata&sourceTaxon=Piperaceae&sourceTaxon=Malvastrum&sourceTaxon=Alternanthera%20philoxeroides&sourceTaxon=Mirabilis%20jalapa&sourceTaxon=Premna%20serratifolia&sourceTaxon=Abutilon%20permolle&sourceTaxon=Deinbollia%20borbonica&sourceTaxon=Nicotiana%20rustica&sourceTaxon=Vigna%20unguiculata&sourceTaxon=Ipomoea%20lacunosa&sourceTaxon=Rhynchosia%20minima&sourceTaxon=Ageratum%20conyzoides&sourceTaxon=Nicotiana%20benthamiana&sourceTaxon=Pedilanthus&sourceTaxon=Sinapis%20arvensis&sourceTaxon=Pouzolzia%20zeylanica&sourceTaxon=Boehmeria%20nivea&sourceTaxon=Dicliptera%20vahliana&sourceTaxon=Brassica%20oleracea%20var.%20capitata&sourceTaxon=Ipomoea%20batatas&sourceTaxon=Cleome&sourceTaxon=Jacquemontia%20tamnifolia&sourceTaxon=Pisum%20sativum&sourceTaxon=Boehmeria&sourceTaxon=Vigna%20mungo&sourceTaxon=Oxalis%20debilis&sourceTaxon=Gossypium%20hirsutum&sourceTaxon=Gossypium%20lobatum&sourceTaxon=Cajanus%20cajan&sourceTaxon=Eupatorium&sourceTaxon=Vigna%20unguiculata%20subsp.%20unguiculata&sourceTaxon=Gomphrena%20globosa&sourceTaxon=Jatropha%20multifida&sourceTaxon=Gossypium%20gossypioides&sourceTaxon=Lamiaceae&sourceTaxon=Ageratum&sourceTaxon=Cucurbita&sourceTaxon=Abelmoschus%20esculentus&sourceTaxon=Melochia&sourceTaxon=Sida%20ciliaris&sourceTaxon=Corchorus%20capsularis&sourceTaxon=Daucus%20carota&sourceTaxon=Desmodium&sourceTaxon=Asystasia&sourceTaxon=Hibiscus%20sabdariffa&sourceTaxon=Trichosanthes%20dioica&sourceTaxon=Eustoma%20grandiflorum&sourceTaxon=Capsicum%20baccatum&sourceTaxon=Solanum%20melongena&sourceTaxon=Solanum&sourceTaxon=Croton&sourceTaxon=Lonicera%20japonica&sourceTaxon=Malvastrum%20coromandelianum&sourceTaxon=Amaranthus%20cruentus&sourceTaxon=Lycianthes%20biflora&sourceTaxon=Amaranthus%20hypochondriacus&sourceTaxon=Petunia&sourceTaxon=Citrullus%20lanatus%20subsp.%20vulgaris&sourceTaxon=Cyanthillium%20cinereum&sourceTaxon=Senna%20occidentalis&sourceTaxon=Leucaena%20leucocephala&sourceTaxon=Erechtites%20valerianifolius&sourceTaxon=Sida%20acuta&sourceTaxon=Sida%20urens&sourceTaxon=Artemisia%20carvifolia&sourceTaxon=Allamanda%20cathartica&sourceTaxon=Onagraceae&sourceTaxon=Sigesbeckia%20orientalis&sourceTaxon=Ipomoea%20setosa&sourceTaxon=Hedyotis%20uncinella&sourceTaxon=Glycine%20max&sourceTaxon=Jatropha%20gossypiifolia&sourceTaxon=Dolichos&sourceTaxon=Passiflora%20edulis&sourceTaxon=Jatropha&sourceTaxon=Blainvillea%20rhomboidea&sourceTaxon=Wissadula&sourceTaxon=Mucuna%20pruriens%20var.%20utilis&sourceTaxon=Hibiscus%20rosa-sinensis&sourceTaxon=Datura%20inoxia&sourceTaxon=Stachytarpheta%20jamaicensis&sourceTaxon=Nicotiana%20clevelandii&sourceTaxon=Fabaceae&sourceTaxon=Lonicera&sourceTaxon=Callianthe%20sellowiana&sourceTaxon=Oxalis%20corniculata&sourceTaxon=Ipomoea%20alba&sourceTaxon=Luffa%20acutangula&sourceTaxon=Mucuna&sourceTaxon=Capsicum&sourceTaxon=Ipomoea%20nil&sourceTaxon=Crotalaria%20juncea&sourceTaxon=Jacquemontia%20pentanthos&sourceTaxon=Vigna%20radiata&sourceTaxon=Telfairia%20occidentalis&sourceTaxon=Capsicum%20baccatum%20var.%20pendulum&sourceTaxon=Pombalia%20attenuata&sourceTaxon=Lindernia%20procumbens&sourceTaxon=Synedrella%20nodiflora&sourceTaxon=Solanum%20aethiopicum&sourceTaxon=Raphanus%20sativus&sourceTaxon=Datura%20stramonium&sourceTaxon=Brassica%20oleracea&sourceTaxon=Alcea%20rosea&sourceTaxon=Eclipta%20prostrata&sourceTaxon=Sonchus%20arvensis&sourceTaxon=Gymnanthemum%20amygdalinum&sourceTaxon=Cucurbita%20pepo&sourceTaxon=Capsicum%20annuum%20var.%20annuum&sourceTaxon=Acmella%20paniculata&sourceTaxon=Distimake%20quinquefolius&sourceTaxon=Euphorbia%20heterophylla&sourceTaxon=Solanum%20lycopersicum%20var.%20cerasiforme&sourceTaxon=Mimosa&sourceTaxon=Physalis&sourceTaxon=Centrosema%20brasilianum&sourceTaxon=Boerhavia%20coccinea&sourceTaxon=Solanum%20americanum&sourceTaxon=Cucumis%20melo&sourceTaxon=Benincasa%20hispida&sourceTaxon=Malvastrum%20americanum&sourceTaxon=Mucuna%20pruriens&sourceTaxon=Plantaginaceae&sourceTaxon=Solanum%20pennellii&sourceTaxon=Phaseolus%20lunatus&sourceTaxon=Cucurbita%20maxima&sourceTaxon=Coccinia%20grandis&sourceTaxon=Catharanthus%20roseus&sourceTaxon=Andrographis%20paniculata&sourceTaxon=Croton%20bonplandianus&sourceTaxon=Cnidoscolus%20urens&sourceTaxon=Urena%20lobata&sourceTaxon=Pavonia&sourceTaxon=Verbena&sourceTaxon=Abutilon&sourceTaxon=Ipomoea%20lobata&sourceTaxon=Hibiscus%20cannabinus&sourceTaxon=Alcea&sourceTaxon=Capraria%20biflora&sourceTaxon=Sida%20cordifolia&sourceTaxon=Solanum%20pimpinellifolium&sourceTaxon=Eupatorium%20makinoi&sourceTaxon=Euphorbia%20pulcherrima&sourceTaxon=Ipomoea%20aquatica&sourceTaxon=Brassica%20oleracea%20var.%20botrytis&sourceTaxon=Corchorus%20olitorius&sourceTaxon=Ageratum%20houstonianum&sourceTaxon=Citrullus%20lanatus&sourceTaxon=Cleome%20affinis&sourceTaxon=Sidastrum%20micranthum&sourceTaxon=Sida%20rhombifolia&sourceTaxon=Manihot%20esculenta&sourceTaxon=Cucumis%20sativus&sourceTaxon=Nicotiana%20tabacum&sourceTaxon=Macroptilium%20lathyroides&sourceTaxon=Jatropha%20podagrica&sourceTaxon=Ludwigia%20hyssopifolia&sourceTaxon=Leonurus%20sibiricus&sourceTaxon=Malachra%20capitata&sourceTaxon=Wissadula%20amplissima&sourceTaxon=Orobanchaceae&sourceTaxon=Rhynchosia&sourceTaxon=Hemidesmus%20indicus&sourceTaxon=Sigmoidotropis%20elegans&sourceTaxon=Sidastrum&sourceTaxon=Emilia%20sonchifolia&sourceTaxon=Duranta%20erecta&sourceTaxon=Gossypium%20davidsonii&sourceTaxon=Dimorphotheca%20sinuata&sourceTaxon=Ipomoea%20carnea%20subsp.%20fistulosa&sourceTaxon=Solanum%20subgen.%20Lycopersicon&interactionType=hostOf&targetTaxon=NCBI:881944&targetTaxon=NCBI:10814&targetTaxon=Insecta&fields=source_taxon_external_id,source_taxon_name,source_taxon_path,source_taxon_path_ids,source_taxon_path_ranks,interaction_type,target_taxon_external_id,target_taxon_name,target_taxon_path,target_taxon_path_ids,target_taxon_path_rankstaxonIdPrefix=NCBI"
response = requests.get(link)
res=response.json()


# In[ ]:


response.status_code


# In[ ]:




