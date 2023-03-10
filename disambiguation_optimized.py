from SPARQLWrapper import SPARQLWrapper
from tqdm import tqdm

def buat_kamus_kosong(df_node):
    # masking bukan NCBI dan null 
    masking = ((df_node["taxon_id"].str.contains("NCBI") == False) & (df_node['taxon_id'].isnull()==False))

    #kamus berdasarkan database biodiversity pada seluruh dataframe
    kamus_ncbi={ a.taxon_id.split(':')[0]:{} for i,a in df_node[masking].iterrows() if a.taxon_id != ""}
    kamus_ncbi={ a.split(':')[0]:{} for idx,data in df_node[masking].iterrows() for a in data.taxon_path_ids.replace(" ","").split('|') if a!="" }

    #{ db:{ kode:arti } }
    for idx,data in df_node[masking].iterrows():
        for i in data.taxon_path_ids.replace(" ","").split('|'):
            if i != '':
                kamus_ncbi[i.split(':')[0]][i]=''
    
    return kamus_ncbi

def query_text(join_ids_text, pred):
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
        'WD':'',
        }
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    
    # looping
    print(list(kamus_ncbi), len(kamus_ncbi),' database, ', len(kamus_ncbi),' kali perulangan akses NCBI')
    for group in list(kamus_ncbi):
        ids = ' '.join(["('"+i.split(":")[1]+"')" for i in list(kamus_ncbi[group])])
        query = query_text(ids, pred[group])
        
        # kalau query lebih kecil dari 2000 lengthnya
        print(len(query))
        if len(query) < 1000:
            #querying
            spw=SPARQLWrapper(endpoint_url)
            spw.setQuery(query)
            spw.setReturnFormat(format_)
            hasil=spw.query().convert()
            # skip kalo kosong
            if(not hasil['results']['bindings']):
                continue
            #update kamus
            for i in tqdm(hasil['results']['bindings']):
                key = group + ':' + i['id']['value']
                if('ncbiLabel' in i):
                    kamus_ncbi[group][key] = "NCBI:"+i['ncbiLabel']['value']
        
        # kalau query lebih kecil dari 2000 lengthnya
        else:
            cek = handle_pake_chunk(kamus_ncbi,pred,group)
            if not cek :
                continue
            else :
                kamus_ncbi = cek

    return kamus_ncbi

def handle_pake_chunk(kamus_ncbi,pred,group):
    print(group, ' query terlalu panjang, skip')
    #querying 
    # skip kalo kosong
    # update kamus
    return False
    # return kamus_ncbi

def update_df_pake_kamus(kamus_ncbi,df_node,df_edge,printOutput=False):
    # masking bukan NCBI dan null 
    masking = ((df_node["taxon_id"].str.contains("NCBI") == False) & (df_node['taxon_id'].isnull()==False))
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
        df_node["taxon_path_ids"].str.contains("NCBI")
    )].iterrows():
        iterkan = list(zip(data.taxon_path_ids.replace(" ","").split("|"),data.taxon_path_rank.replace(" ","").split("|")))
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
        
    return df_node,df_edge


