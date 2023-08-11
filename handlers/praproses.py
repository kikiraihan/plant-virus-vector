import requests
from tqdm import tqdm
import pandas as pd
# from owlready2 import get_ontology
import os, json

from modul.standardization_usingsparql import addTaxonColumn, buat_kolom_taxon_awal
from modul.disambiguation_optimized import buat_kamus_kosong, update_kamus_pake_wikidata, update_df_pake_kamus, update_df_pake_path_ujung, removeOtherThanNCBI
from modul.preprocess import cleaning, splitInteractionToNodeEdge
from modul.filterNodeEdge import removeNodeAndEdgeByFilter,takeNodeAndEdgeByFilter,removeEdgesNotInNodes
from modul.helper_umum import contains_string_entire_column,contains_string_entire_column_boolean
from modul.vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id_via_string
from handlers.praproses_dari_proses import pra_proses_dari_proses

def report_back(progress,message):
    data=json.dumps({
        "progress":progress,
        'message':message,
        })
    return f'data: {data}\n\n'

def praproses(virus_txt):
    # 0 parameter
    yield report_back(5,'inisiasi parameter')
    nama_file = "tes"
    virus_txt = virus_txt.replace(' ','%20')
    tipe_interaksi_virus = 'hasHost' #pathogenOf, pake relasi hasHost lebih dapat banyak relasi dari pada pathogenOf
    tipe_interaksi_tanaman = 'hostOf' #hasPathogen, pake relasi hostOf lebih dapat banyak relasi dari pada hasPathogen
    tipe_interaksi_serangga_ke_tanaman = 'hasHost' 
    tipe_interaksi_serangga_ke_virus = 'hostOf' 
    ncbi_server_url = 'http://localhost:3030/mydataset/query'

    virus_search = get_taxon_vector(virus_txt,ncbi_server_url)
    if (virus_search==False):
        virus_search=[('unknown',virus_txt)]

    yield report_back(10,'BFS Data virus')
    #1 BFS Data
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
    print('layer 1 pertama',virus_search[0][1].split('_'))
    link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+virus_search[0][1].split('_')[0]+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
    response = requests.get(link)
    if response.status_code==500:
        print('ADA ERROR DARI SERVER GLOBI')
    res=response.json()
    if not res['data']:
        print(virus_search, ': kosong dari GloBI')
    # JSON To Pandas Dataframe
    df = pd.json_normalize(res, record_path =['data'])

    if(len(virus_search) == 1):
        pass
    elif(virus_search[1][0] not in ['famili','genus']):
        pass
    else:
        print('layer 1 kedua',virus_search[1][1].split('_')) # yg kedua dari famili atau genusnya.
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+virus_search[1][1].split('_')[0]+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        if not res['data']:
            print('kosong')
        # JSON To Pandas Dataframe
        df_ = pd.json_normalize(res, record_path =['data'])
        df=pd.concat([df,df_],ignore_index=True)

    df.columns = kolom

    yield report_back(20,'Splitting layer 1 interaksi virus')
    #2 splitting layer 1 interaksi virus
    df_node, df_edge = splitInteractionToNodeEdge(df)

    # cleaning_after_split layer 1 interaksi virus
    # drop duplikat
    df_node.drop_duplicates(inplace=True)
    # no_ncbi dan path_null
    no_ncbi_and_path_null=(df_node.taxon_id.str.contains('NCBI')==False) & (df_node.taxon_path_ids.isnull())
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[no_ncbi_and_path_null], df_node,df_edge) 
    # drop duplikat
    df_edge.drop_duplicates(inplace=True)

    # tandai virus utama
    filter_virus_utama=(
        (df_node.taxon_name.str.contains(r'\b(virus\w*|\w*virus)\b',case=False))
        | (df_node.taxon_path.str.contains(r'\b(virus\w*|\w*virus)\b', case=False))  
        #jika berawalan atau berakhiran kata virus
    )
    # df_node.loc[filter_virus_utama, ['virus_utama']] = True
    virus_utama=[data.taxon_id for idx,data in df_node[filter_virus_utama].iterrows()]


    yield report_back(30,'disambiguasi layer 1 interaksi virus, ini akan memakan waktu sedikit lama')
    #3 disambiguasi layer 1 interaksi virus
    kamus_ncbi = buat_kamus_kosong(df_node)
    kamus_ncbi = update_kamus_pake_wikidata(kamus_ncbi)
    #update dataframe pake kamus
    df_node,df_edge = update_df_pake_kamus(kamus_ncbi,df_node,df_edge)
    df_node,df_edge = update_df_pake_path_ujung(df_node,df_edge)
    #tambah kolom takson pake data NCBI
    df_node = buat_kolom_taxon_awal(df_node) #buat kolom taxon, default none
    df_node = addTaxonColumn(df_node,'http://localhost:3030/mydataset/query') # isi pake ncbi




    # cleaning_after_disambiguasi layer 1
    df_node, df_edge = removeOtherThanNCBI(df_node,df_edge)# Hapus kalo masih ada selain NCBI
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus edge yang tidak ada nodenya
    filter_kingdom_atau_class_null=( (df_node.kingdom.isnull()) | (df_node['class'].isnull()) )
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_kingdom_atau_class_null], df_node,df_edge)



    yield report_back(40,'BFS interaksi tanaman')
    #4.1 BFS interaksi tanaman
    df_to_add=pd.DataFrame(columns = kolom)
    df_plant=df_node[df_node.kingdom=='NCBI:33090_Viridiplantae']
    interactionType=tipe_interaksi_tanaman
    for idx,i in tqdm(df_plant.iterrows(), total=df_plant.shape[0]):
        search=i.taxon_name.replace(' ','%20')
        # print('depth 2 tanaman :', search)
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+search+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        if not res['data']:
            print(i.taxon_name, ': kosong dari GloBI')
            continue
        # JSON To Pandas Dataframe
        temp_to_add=pd.json_normalize(res, record_path =['data'])
        temp_to_add.columns = kolom
        # add to sebelumnya
        df_to_add = pd.concat([
            df_to_add,
            temp_to_add
        ], ignore_index = True)



    yield report_back(45,'BFS interaksi serangga -> tanaman')
    #4.2 BFS interaksi serangga -> tanaman
    df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
    interactionType = tipe_interaksi_serangga_ke_tanaman
    for idx,i in tqdm(df_insect.iterrows(), total=df_insect.shape[0]):
        search=i.taxon_name.replace(' ','%20')
        # print('depth 2 tanaman :', search)
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+search+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        if not res['data']:
            print(i.taxon_name, ': kosong dari GloBI')
            continue
        # JSON To Pandas Dataframe
        temp_to_add=pd.json_normalize(res, record_path =['data'])
        temp_to_add.columns = kolom
        # add to sebelumnya
        df_to_add = pd.concat([
            df_to_add,
            temp_to_add
        ], ignore_index = True)



    yield report_back(48,'BFS interaksi serangga -> virus')
    #4.3 BFS interaksi serangga -> virus
    df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
    interactionType = tipe_interaksi_serangga_ke_virus
    for idx,i in tqdm(df_insect.iterrows(), total=df_insect.shape[0]):
        search=i.taxon_name.replace(' ','%20')
        # print('depth 2 tanaman :', search)
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+search+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        if not res['data']:
            print(i.taxon_name, ': kosong dari GloBI')
            continue
        # JSON To Pandas Dataframe
        temp_to_add=pd.json_normalize(res, record_path =['data'])
        temp_to_add.columns = kolom
        # add to sebelumnya
        df_to_add = pd.concat([
            df_to_add,
            temp_to_add
        ], ignore_index = True)



    yield report_back(50,'splitting depth 2 interaksi serangga dan tanaman')
    #5 splitting depth 2 interaksi serangga dan tanaman
    node_to_add, edge_to_add = splitInteractionToNodeEdge(df_to_add)




    # cleaning_after_split depth 2 interaksi serangga dan tanaman
    # hapus edge inverse
    kebalikan={
        'hostOf':'hasHost',
        'hasPathogen':'pathogenOf', 
        'pollinatedBy':'pollinates', 
        'flowersVisitedBy':'visitFlowersOf',
        'visitedBy':'visit'
    }
    for i,data in  edge_to_add[edge_to_add['interaction_type'].isin([
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



    yield report_back(60,'disambiguasi layer 2, ini akan memakan waktu sedikit lama')
    # 6 disambiguasi layer 2
    kamus_ncbi = buat_kamus_kosong(node_to_add)
    kamus_ncbi = update_kamus_pake_wikidata(kamus_ncbi)
    #update dataframe pake kamus
    node_to_add,edge_to_add = update_df_pake_kamus(kamus_ncbi,node_to_add,edge_to_add)
    node_to_add,edge_to_add = update_df_pake_path_ujung(node_to_add, edge_to_add)
    # tambah kolom takson pake data NCBI
    node_to_add = buat_kolom_taxon_awal(node_to_add) #buat kolom taxon, isi none dan isi dari path
    node_to_add = addTaxonColumn(node_to_add,'http://localhost:3030/mydataset/query') #isi kolom taxon, pake NCBI



    # cleaning_after_disambiguasi depth 2
    node_to_add,edge_to_add = removeOtherThanNCBI(node_to_add,edge_to_add) #hapus selain NCBI



    yield report_back(70,'konkatenasi tabel')
    #7 konkatenasi tabel
    df_node=pd.concat([df_node,node_to_add], axis=0, ignore_index=True)
    df_edge=pd.concat([df_edge,edge_to_add], axis=0, ignore_index=True)



    yield report_back(80,'praproses tambahan')
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

    # hapus yang no group
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[df_node.group=="nogroup"], df_node,df_edge) 
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #edge menyesuaikan #cuma memastikan saja

    # hapus node yang tidak punya edge
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

    df_node, df_edge = pra_proses_dari_proses(df_node, df_edge)

    #12
    yield report_back(98,'save file')
    # akhir pra proses
    # save file
    # Mengirimkan data hasil proses sebagai respons JSON
    # format orientasi # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_json.html
    # orient="split" cocok dengan dataframe js.
    result = json.dumps({
        'status':200,
        'message': 'Kirim json', 
        'df_node':df_node.to_json(orient="split"),
        'df_edge':df_edge.to_json(orient="split"),
        })
    yield f'data: {(result)}\n\n'

    # selesai
    yield report_back(100,'selesai')