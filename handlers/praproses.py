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

def report_back(progress,message):
    data=json.dumps({
        "progress":progress,
        'message':message,
        })
    return f'data: {data}\n\n'

def praproses(virus_txt):
    # 1 parameter
    yield report_back(1,'inisiasi parameter')
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

    #2 ambil data dari globi
    yield report_back(5,'ambil data interaksi virus dari globi')
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
    print('depth 1 pertama',virus_search[0][1].split('_'))
    link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+virus_search[0][1].split('_')[0]+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
    response = requests.get(link)
    res=response.json()
    if not res['data']:
        print('stop kosong')
    # JSON To Pandas Dataframe
    df = pd.json_normalize(res, record_path =['data'])

    if(len(virus_search) == 1):
        pass
    elif(virus_search[1][0] not in ['famili','genus']):
        pass
    else:
        print('depth 1 kedua',virus_search[1][1].split('_'))
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+virus_search[1][1].split('_')[0]+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        if not res['data']:
            print('kosong')
        # JSON To Pandas Dataframe
        df_ = pd.json_normalize(res, record_path =['data'])
        df=pd.concat([df,df_])

    df.columns = kolom

    yield report_back(5,'split data interaksi virus')
    #3 splitting depth 1 interaksi virus
    df_node, df_edge = splitInteractionToNodeEdge(df)

    yield report_back(10,'cleaning setelah ambil data depth 1')
    #4 cleaning_after_get depth 1 interaksi virus
    df_node.drop_duplicates(inplace=True)
    no_ncbi_and_path_null=(df_node.taxon_id.str.contains('NCBI')==False) & (df_node.taxon_path_ids.isnull())
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[no_ncbi_and_path_null], df_node,df_edge) 
    df_edge.drop_duplicates(inplace=True)

    #4a tandai virus utama
    filter_virus_utama=(
        (df_node.taxon_name.str.contains('virus',case=False))
        | (df_node.taxon_path.str.contains('virus', case=False)) 
        #jika berawalan atau berakhiran kata virus
    )
    # df_node.loc[filter_virus_utama, ['virus_utama']] = True
    virus_utama=[data.taxon_id for idx,data in df_node[filter_virus_utama].iterrows()]

    yield report_back(15,'disambiguasi data depth 1 (interaksi virus)')
    #5 disambiguasi depth 1 interaksi virus
    kamus_ncbi = buat_kamus_kosong(df_node)
    kamus_ncbi = update_kamus_pake_wikidata(kamus_ncbi)
    #update dataframe pake kamus
    df_node,df_edge = update_df_pake_kamus(kamus_ncbi,df_node,df_edge)
    df_node,df_edge = update_df_pake_path_ujung(df_node,df_edge)

    #6 standarisasi depth 1 interaksi virus
    df_node = buat_kolom_taxon_awal(df_node) #buat kolom taxon, default none
    df_node = addTaxonColumn(df_node,'http://localhost:3030/mydataset/query') # isi pake ncbi
    df_node, df_edge = removeOtherThanNCBI(df_node,df_edge)# Hapus kalo masih ada selain NCBI
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus edge yang tidak ada nodenya

    yield report_back(20,'cleaning setelah disambiguasi depth 1 (interaksi virus)')
    # cleaning_after_disambiguasi depth 1
    # 7 hapus node yang ordo_sampai_species_null
    filter_ordo_sampai_species_null=(
        (df_node.order.isnull())
        & (df_node.family.isnull())
        & (df_node.genus.isnull())
        & (df_node.species.isnull())
    )
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_ordo_sampai_species_null], df_node,df_edge)

    yield report_back(25,'ambil data tanaman')
    #7 pengambilan interaksi tanaman
    df_plant=df_node[df_node.kingdom=='NCBI:33090_Viridiplantae']
    interactionType=tipe_interaksi_tanaman
    node_to_add=pd.DataFrame(columns = [
        'taxon_id',
        'taxon_name',
        'taxon_path',
        'taxon_path_ids',
        'taxon_path_rank',
    ])
    edge_to_add=pd.DataFrame(columns = [
        'source_taxon_id',
        'target_taxon_id',
        'interaction_type',
    ])
    kebalikan={
        'hostOf':'hasHost',
        'hasPathogen':'pathogenOf', 
        'pollinatedBy':'pollinates', 
        'flowersVisitedBy':'visitFlowersOf',
        'visitedBy':'visit'
    }
    #tanaman
    for idx,i in tqdm(df_plant.iterrows(), total=df_plant.shape[0]):
        plant=i.taxon_name.replace(' ','%20')
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+plant+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        #list data
        taxon_id=[]
        taxon_name=[]
        taxon_path=[]
        taxon_path_ids=[]
        taxon_path_rank=[]
        interaction_type=[]
        for x in res['data']:
            taxon_id.append(x[6])
            taxon_name.append(x[7])
            taxon_path.append(x[8])
            taxon_path_ids.append(x[9])
            taxon_path_rank.append(x[10])
            if x[5] in ['hostOf', 'hasPathogen', 'pollinatedBy', 'flowersVisitedBy','visitedBy']:
                interaction_type.append(kebalikan[x[5]])
            else:
                interaction_type.append('patogennya')
        ## NODE
        #concat ke dataframe lama
        node_to_add = pd.concat([
            node_to_add,
            pd.DataFrame({
                'taxon_id':taxon_id,
                'taxon_name':taxon_name,
                'taxon_path':taxon_path,
                'taxon_path_ids':taxon_path_ids,
                'taxon_path_rank':taxon_path_rank,
            })
        ], ignore_index = True)
        node_to_add.reset_index()
        ## EDGE
        #concat ke dataframe lama
        edge_to_add = pd.concat([
            edge_to_add,
            pd.DataFrame({
                'source_taxon_id': taxon_id, # patogen
                'target_taxon_id': i.taxon_id, # tanaman
                'interaction_type': interaction_type})
        ], ignore_index = True)
        edge_to_add.reset_index()

    yield report_back(30,'ambil data serangga->tanaman')
    # 7 pengambilan interaksi serangga -> tanaman
    df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
    interactionType = tipe_interaksi_serangga_ke_tanaman
    #serangga
    for idx,i in tqdm(df_insect.iterrows(), total=df_insect.shape[0]):
        insect=i.taxon_name.replace(' ','%20')
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+insect+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        #list data
        taxon_id=[]
        taxon_name=[]
        taxon_path=[]
        taxon_path_ids=[]
        taxon_path_rank=[]
        interaction_type=[]
        for x in res['data']:
            taxon_id.append(x[6])
            taxon_name.append(x[7])
            taxon_path.append(x[8])
            taxon_path_ids.append(x[9])
            taxon_path_rank.append(x[10])
            interaction_type.append(x[5])
        ## NODE
        #concat ke dataframe lama
        node_to_add = pd.concat([
            node_to_add,
            pd.DataFrame({
                'taxon_id':taxon_id,
                'taxon_name':taxon_name,
                'taxon_path':taxon_path,
                'taxon_path_ids':taxon_path_ids,
                'taxon_path_rank':taxon_path_rank,
            })
        ], ignore_index = True)
        node_to_add.reset_index()
        ## EDGE
        #concat ke dataframe lama
        edge_to_add = pd.concat([
            edge_to_add,
            pd.DataFrame({
                'source_taxon_id': i.taxon_id, # serangga
                'target_taxon_id': taxon_id, # inagnya
                'interaction_type': interaction_type
                })
        ], ignore_index = True)
        edge_to_add.reset_index()

    yield report_back(35,'ambil data serangga->virus')
    # 7 pengambilan interaksi serangga -> virus
    df_insect = df_node[df_node['class']=='NCBI:50557_Insecta']
    interactionType = tipe_interaksi_serangga_ke_virus
    #serangga
    for idx,i in tqdm(df_insect.iterrows(), total=df_insect.shape[0]):
        insect=i.taxon_name.replace(' ','%20')
        link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+insect+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
        response = requests.get(link)
        res=response.json()
        #list data
        taxon_id=[]
        taxon_name=[]
        taxon_path=[]
        taxon_path_ids=[]
        taxon_path_rank=[]
        interaction_type=[]
        for x in res['data']:
            taxon_id.append(x[6])
            taxon_name.append(x[7])
            taxon_path.append(x[8])
            taxon_path_ids.append(x[9])
            taxon_path_rank.append(x[10])
            if x[5] in ['hostOf', 'hasPathogen', 'pollinatedBy', 'flowersVisitedBy','visitedBy']:
                interaction_type.append(kebalikan[x[5]])
            else:
                interaction_type.append('patogennya')
        ## NODE
        #concat ke dataframe lama
        node_to_add = pd.concat([
            node_to_add,
            pd.DataFrame({
                'taxon_id':taxon_id,
                'taxon_name':taxon_name,
                'taxon_path':taxon_path,
                'taxon_path_ids':taxon_path_ids,
                'taxon_path_rank':taxon_path_rank,
            })
        ], ignore_index = True)
        node_to_add.reset_index()
        ## EDGE
        #concat ke dataframe lama
        edge_to_add = pd.concat([
            edge_to_add,
            pd.DataFrame({
                'source_taxon_id': taxon_id, # virus
                'target_taxon_id': i.taxon_id, # serangga
                'interaction_type': interaction_type
                })
        ], ignore_index = True)
        edge_to_add.reset_index()

    print(len(node_to_add),len(edge_to_add))

    yield report_back(40,'cleaning setelah ambil data depth 2')
    #8 cleaning_after_get depth 2 interaksi tanaman dan serangga
    node_to_add.drop_duplicates(inplace=True)
    hapus=[i for i,d in node_to_add[
        (node_to_add.taxon_id.str.contains('NCBI')==False) & (node_to_add.taxon_path_ids.isnull())
    ].iterrows()]
    node_to_add.drop(hapus, inplace=True)
    edge_to_add.drop_duplicates(inplace=True)

    yield report_back(45,'disambiguasi data depth 2 (interaksi serangga dan tanaman)')
    # 9 disambiguasi depth 2
    kamus_ncbi = buat_kamus_kosong(node_to_add)
    kamus_ncbi = update_kamus_pake_wikidata(kamus_ncbi)
    #update dataframe pake kamus
    node_to_add,edge_to_add = update_df_pake_kamus(kamus_ncbi,node_to_add,edge_to_add)
    node_to_add,edge_to_add = update_df_pake_path_ujung(node_to_add, edge_to_add)

    #10 standarisasi depth 2
    node_to_add = buat_kolom_taxon_awal(node_to_add) #buat kolom taxon, isi none dan isi dari path
    node_to_add = addTaxonColumn(node_to_add,'http://localhost:3030/mydataset/query') #isi kolom taxon, pake NCBI
    node_to_add,edge_to_add = removeOtherThanNCBI(node_to_add,edge_to_add) #hapus selain NCBI

    yield report_back(50,'concat data virus <-> serangga, tanaman')
    #10a # concat dengan data lama
    df_node=pd.concat([df_node,node_to_add], axis=0)
    df_edge=pd.concat([df_edge,edge_to_add], axis=0)

    yield report_back(60,'cleaning tambahan')
    #12 praproses tambahan
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

    yield report_back(70,'cleaning_after_disambiguasi dilakukan setelah pengelompokan')
    # cleaning_after_disambiguasi layer 2
    # 10b hapus node yang ordo sampai specie isi null
    filter_ordo_sampai_species_null=(
        (df_node.order.isnull()) & 
        (df_node.family.isnull()) & 
        (df_node.genus.isnull()) &
        (df_node.species.isnull())
    )
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_ordo_sampai_species_null], df_node,df_edge)

    # 10c pra proses akhir
    df_node.drop_duplicates(subset=["taxon_id"], keep='last',inplace=True)#hapus duplikasi node
    df_edge = removeEdgesNotInNodes(df_node, df_edge)#edge menyesuaikan
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[(df_node.kingdom.isnull()) & (df_node.group!='virus')], df_node,df_edge) # hapus kingdom isi null
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus lagi edge kalo masih ada edge tidak ada di node #cuma memastikan saja
    # reset index
    df_node.reset_index(drop=True,inplace=True)
    df_edge.reset_index(drop=True,inplace=True)

    # 13 # hapus yang no group
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[df_node.group=="nogroup"], df_node,df_edge) 
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus lagi edge kalo masih ada edge tidak ada di node

    #14 # tambahan, hapus node yang tidak punya edge
    yield report_back(80,'hapus node yang tidak ada di edge (tidak punya edge)')
    print('hapus node yang tidak ada di edge (tidak punya edge)')
    sebelum=len(df_node)
    df_node = df_node[
        (df_node.taxon_id.isin(df_edge.source_taxon_id.unique())) 
        | (df_node.taxon_id.isin(df_edge.target_taxon_id.unique()))
    ]
    print(sebelum, '->', len(df_node))

    yield report_back(95,'cleaning tambahan')
    # 15 # masukan tanda virus utama
    df_node.loc[df_node.taxon_id.isin(virus_utama), ['virus_utama']] = True

    #12
    yield report_back(98,'save file')
    # akhir pra proses
    # save file
    # Mengirimkan data hasil proses sebagai respons JSON
    # format orientasi # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_json.html
    # orient="split" cocok dengan dataframe js.
    result = json.dumps({
        'message': 'Kirim json', 
        'df_node':df_node.to_json(orient="split"),
        'df_edge':df_edge.to_json(orient="split"),
        })
    yield f'data: {(result)}\n\n'

    # selesai
    yield report_back(100,'selesai')