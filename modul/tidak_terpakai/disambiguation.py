#kalau pake SPARQL get clasification ke arah parent
from SPARQLWrapper import SPARQLWrapper, JSON, N3, TURTLE

def removeOtherThanNCBI(node,edge):
    #Cek apakah masih ada selain NCBI?
    print({a.taxon_id.split(':')[0] for i,a in node[node["taxon_id"].str.contains("NCBI") == False].iterrows()})
    print('------------------------------------')

    #preview yang dihapus selain NCBI
    print('node')
    print('sebelum ', len(node))
    print('sesudah ', len(node) - len(node[node["taxon_id"].str.contains("NCBI") == False]))
    print('edge')
    print('target: sebelum', len(edge))
    print('target: sesudah', len(edge)- len(edge[edge["target_taxon_id"].str.contains("NCBI") == False]))
    print('source: sebelum', len(edge))
    print('source: sesudah', len(edge)- len(edge[edge["source_taxon_id"].str.contains("NCBI") == False]))
    #drop selain NCBI
    node = node[node["taxon_id"].str.contains("NCBI")]
    edge = edge[edge["target_taxon_id"].str.contains("NCBI")]
    edge = edge[edge["source_taxon_id"].str.contains("NCBI")]

    return (node,edge)

# TIDAK TEROPTIMASI, BANYAK ERROR REQUEST DI SPARQL
def getConvertNCBI(id_,type_):
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    
    pred={
        'GBIF':'wdt:P846',
        'EOL':'wdt:P830',
        'ITIS':'wdt:P815',
        'COL':'wdt:P10585',
        'INAT_TAXON':'wdt:P3151',
        'IRMNG':'wdt:P5055',
        'NBN':'wdt:P3240',
        'WD':'',
    }
    
    if(type_ not in pred):
        print("tidak ditemukan tipe: ",type_)
        return -1
    
    
    if(type_=='WD'):
        query="""
        SELECT ?patogenLabel ?rankLabel ?ncbiLabel 
        WHERE {
            wd:"""+id_+""" wdt:P105 ?rank.
            wd:"""+id_+""" wdt:P685 ?ncbi.
            wd:"""+id_+""" wdt:P225 ?patogen.

            SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
        }
        """
    else:    
        query="""
        SELECT ?patogenLabel ?rankLabel ?ncbiLabel
        WHERE {
            ?patogen """+pred[type_]+""" '"""+id_+"""'.
            ?patogen wdt:P105 ?rank.
            OPTIONAL{?patogen wdt:P685 ?ncbi}.

            SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
        }
        """
    
    #querying
    spw=SPARQLWrapper(endpoint_url)
    spw.setQuery(query)
    spw.setReturnFormat(format_)
    hasil = spw.query().convert()
    
    if(not hasil['results']['bindings']):
        print("tidak ditemukan id: ",id_)
        return -2
    
    hasil = hasil['results']['bindings'][0]
    
    if 'ncbiLabel' in hasil:
        return (
            hasil['patogenLabel']['value'],
            hasil['rankLabel']['value'],
            "NCBI:"+hasil['ncbiLabel']['value']
            )
    else: 
        print("tidak ada id NCBI")
        return -3



def disambiguation_lama(node,edge):
    
    iterkan=node[node["taxon_id"].str.contains("NCBI") == False]
    print('total : ',iterkan.shape[0])
    
    for idx,dat in iterkan.iterrows(): #tqdm(iterkan.iterrows(),total=iterkan.shape[0])
        print(idx,"-======-")
        cek=dat['taxon_id']
        print(dat['taxon_name'], cek.split(':'))
        hasil=getConvertNCBI(cek.split(':')[1],cek.split(':')[0])
        
        # jika hasil berupa int dan kode error kurang dari -2 
        # (tidak ada id NCBI atau tidak ditemukan id), cek NCBI pada tiap parent.
        if(type(hasil)==int):
            a=-2
            tpti=dat['taxon_path_ids'].split('|')
            while type(hasil)==int and a>=len(tpti)*-1:
                if(hasil>-2): #kalau -3 brrti masalah type, tidak masuk.
                    break
                cek=tpti[a].strip()
                if(cek.split(':')==['']): # kalau di taxon path ada nilai kosong, pass
                    a=a-1
                    continue
                print(dat['taxon_name'], cek.split(':'))
                hasil=getConvertNCBI(cek.split(':')[1],cek.split(':')[0])                
                a=a-1

        #update kalo dapat (tuple)
        if(type(hasil)==tuple):
            node.loc[idx,['taxon_name','taxon_rank','taxon_id']] = list(hasil)
            #node.loc[idx,['taxon_path','taxon_path_ids','taxon_path_rank']] = [None,None,None]
            edge.replace(
                dat['taxon_id'], 
                hasil[2],
                inplace=True)

        print('perubahan : ',hasil)
        
    return (node,edge)