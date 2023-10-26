from SPARQLWrapper import SPARQLWrapper
from tqdm import tqdm
import pandas as pd

def __query_text(join_ids_text, pred):
    if(pred=='iniwd'):
        join_ids_text = join_ids_text.replace("('", "(wd:").replace("')", ")")
        query="""
        SELECT ?id ?ncbiLabel ?patogenLabel ?rankLabel
        WHERE {
            VALUES (?entity) { """+join_ids_text+""" } 
            ?entity wdt:P105 ?rank.
            ?entity wdt:P685 ?ncbi.
            ?entity wdt:P225 ?patogen.
            BIND(REPLACE(STR(?entity), ".*Q", "Q") AS ?id) 

            SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
        }
        """
    else:  
        query="""
        SELECT ?id ?ncbiLabel ?patogenLabel ?rankLabel
        WHERE {
            VALUES (?id) { """+join_ids_text+""" }
            ?patogen """+pred+""" ?id.
            ?patogen wdt:P105 ?rank.
            OPTIONAL{?patogen wdt:P685 ?ncbi}.

            SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
        }
        """
    return query.replace("\n"," ").strip().replace("  ","")

def __querying(endpoint_url,query,format_):
    spw=SPARQLWrapper(endpoint_url)
    spw.setQuery(query)
    spw.setReturnFormat(format_)
    hasil=spw.query().convert()
    return hasil

def __updating(hasil,group,kamus_ncbi):
    #update kamus
    for i in tqdm(hasil['results']['bindings']):
        key = group + ':' + i['id']['value']
        if('ncbiLabel' in i):
            kamus_ncbi[group][key] = "NCBI:"+i['ncbiLabel']['value']
    
    return kamus_ncbi

def __chunk_list(lst, chunk_size):
    # my_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # print(__chunk_list(my_list, 3))
    # [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10]]
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

def buat_kamus_kosong(df_node):
    # masking bukan NCBI dan null 
    masking = ((df_node["taxon_id"].str.contains("NCBI") == False) & (df_node['taxon_id'].isnull()==False))

    #kamus berdasarkan database biodiversity pada seluruh dataframe
    # untuk taxon_id
    kamus_ncbi={ a.taxon_id.split(':')[0]:{} for i,a in df_node[masking].iterrows() if a.taxon_id != ""}
    # untuk taxon_path_ids
    for idx,data in df_node[masking].iterrows():
        for a in data.taxon_path_ids.replace(" ","").split('|'):
            if a!="":
                kamus_ncbi[a.split(':')[0]]={}
        # versi comprehension tapi malah menimpa
        # kamus_ncbi={ a.split(':')[0]:{} for idx,data in df_node[masking].iterrows() for a in data.taxon_path_ids.replace(" ","").split('|') if a!="" }

    #{ db:{ kode:arti } }
    for idx,data in df_node[masking].iterrows():
        # untuk taxon_id
        i = data.taxon_id
        kamus_ncbi[i.split(':')[0]][i]=''
        # untuk taxon_path_ids
        for i in data.taxon_path_ids.replace(" ","").split('|'):
            if i != '':
                kamus_ncbi[i.split(':')[0]][i]=''
    
    return kamus_ncbi

def update_kamus_pake_wikidata(kamus_ncbi):
    pred={
        'GBIF':'wdt:P846',
        'EOL':'wdt:P830',
        'EOL_V2':'wdt:P830',
        'ITIS':'wdt:P815',
        'COL':'wdt:P10585',
        'INAT_TAXON':'wdt:P3151',
        'IRMNG':'wdt:P5055',
        'NBN':'wdt:P3240',
        'WD':'iniwd',
        }
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    
    # looping
    print(list(kamus_ncbi), len(kamus_ncbi),' database, ', len(kamus_ncbi),' kali perulangan akses NCBI')
    for group in list(kamus_ncbi):

        if group not in pred:
            print(group,': tidak diketahui predikatnya')
            continue

        ids = ' '.join(["('"+i.split(":")[1]+"')" for i in list(kamus_ncbi[group])])
        query = __query_text(ids, pred[group])
        
        # kalau len <=130
        print(group,': jumlah id',len(kamus_ncbi[group])) # print(len(query))
        if len(kamus_ncbi[group]) <= 130:
            hasil = __querying(endpoint_url,query,format_)
            # skip kalo kosong
            if(not hasil['results']['bindings']):
                continue
            kamus_ncbi = __updating(hasil,group,kamus_ncbi)
        
        # kalau len > 130
        else:
            print(group, ': query terlalu panjang, dilakukan chunk')
            for list_chunk in __chunk_list(list(kamus_ncbi[group]), 130):
                ids = ' '.join(["('"+i.split(":")[1]+"')" for i in list_chunk])
                query = __query_text(ids, pred[group])

                #querying
                hasil = __querying(endpoint_url,query,format_)
                # skip kalo kosong
                if(not hasil['results']['bindings']):
                    continue
                kamus_ncbi = __updating(hasil,group,kamus_ncbi)

    return kamus_ncbi




def update_df_pake_kamus(kamus_ncbi,df_node,df_edge,printOutput=False):
    # masking bukan NCBI dan null 
    masking = (
        (df_node["taxon_id"].str.contains("NCBI") == False) & 
        (df_node['taxon_id'].isnull()==False)
    )
    masking_path_ids = (
        (df_node["taxon_path_ids"].str.contains("NCBI") == False) & 
        (df_node['taxon_path_ids'].isnull() == False)
    )
    
    # update data
    for idx, data in df_node[(masking_path_ids) & (masking)].iterrows():
        taxon_id = data.taxon_id

        # taxon_id
        if taxon_id.split(":")[0] not in ("NCBI",""):
            isi = kamus_ncbi[taxon_id.split(":")[0]][taxon_id]
            if isi != "":
                if printOutput:
                    print(idx)
                    print('sebelumnya', data.taxon_id)
                    print('menjadi', isi)
                #node
                df_node.loc[idx,'taxon_id'] = isi
                #edge
                df_edge.replace(
                    taxon_id,
                    isi,
                    inplace=True
                )

        # taxon_path_ids
        list_isi=[]
        for i in data.taxon_path_ids.replace(" ","").split('|'):
            if i.split(":")[0] not in ("NCBI",""):
                isi = kamus_ncbi[i.split(":")[0]][i]
                if isi != "":
                    list_isi.append(isi)
                else:
                    list_isi.append(i)
        if printOutput:
            print(idx)
            print('sebelumnya', data.taxon_path_ids)
            print('menjadi', (" | ".join(list_isi)))
        df_node.loc[idx,'taxon_path_ids'] = " | ".join(list_isi)
    
    df_node.reset_index(drop=True,inplace=True)
    df_edge.reset_index(drop=True,inplace=True)
    return df_node,df_edge

def update_df_pake_path_ujung(df_node, df_edge,printOutput=False):

    #masking lama tida berlaku, harus update masking baru, karena dataframe baru diupdate
    #masking = ((df_node["taxon_id"].str.contains("NCBI") == False) & (df_node['taxon_id'].isnull()==False))
    #masking_path_ids = (
    #    (df_node["taxon_path_ids"].str.contains("NCBI") == False) & 
    #    (df_node['taxon_path_ids'].isnull() == False)
    #)

    # update taxon_id, yang taxon_idnya bukan NCBI tapi pathnya punya NCBI
    for idx,data in df_node[(
        (df_node["taxon_id"].str.contains("NCBI")==False) & 
        (df_node["taxon_path_ids"].str.contains("NCBI"))
    )].iterrows():
        taxon_path_ids = data.taxon_path_ids.replace(" ","").split("|")
        if pd.isna(data.taxon_path_rank):
            taxon_path_rank = [""]*len(taxon_path_ids)
        else:
            taxon_path_rank = data.taxon_path_rank.replace(" ","").split("|") 
        zipkan = zip( 
            taxon_path_ids, 
            taxon_path_rank
        )
        iterkan = list(zipkan)
        #update dengan nilai NCBI pertama dari belakang
        counter=-1
        cek = iterkan[counter]
        while "NCBI" not in cek[0]:
            counter = counter - 1
            try:
                cek = iterkan[counter]
            except IndexError:
                break
        if printOutput:
            print(idx)
            print('sebelumnya', (data.taxon_id, data.taxon_path_rank.replace(" ","").split("|")[-1]) )
            print('menjadi', cek)
        df_node.loc[idx,'taxon_id'] = cek[0]
        df_node.loc[idx,'taxon_rank'] = cek[1]
        df_edge.replace(
                    data.taxon_id,
                    cek[0],
                    inplace=True
                )
    
    df_node.reset_index(drop=True,inplace=True)
    df_edge.reset_index(drop=True,inplace=True)
    return df_node,df_edge


def removeOtherThanNCBI(node,edge):
    print("removeOtherThanNCBI")
    print("sebelum :",len(node),len(edge))
    #drop selain NCBI
    node = node[node["taxon_id"].str.contains("NCBI")]
    edge = edge[edge["target_taxon_id"].str.contains("NCBI")]
    edge = edge[edge["source_taxon_id"].str.contains("NCBI")]

    print("sesudah :",len(node),len(edge))
    node.reset_index(drop=True,inplace=True)
    edge.reset_index(drop=True,inplace=True)
    return (node,edge)