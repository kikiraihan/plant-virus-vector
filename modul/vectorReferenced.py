from modul.standardization_usingsparql import getTaxonomy
from SPARQLWrapper import SPARQLWrapper
import requests

# pencarian ID berdasarkan nama entity
def get_wiki_id_by_name_api(param_string):
  try:
    url = f"https://www.wikidata.org/w/api.php?action=wbsearchentities&search={param_string}&language=en&format=json"
    data = requests.get(url).json()
    # Return the first id (Could upgrade this in the future)
    return data['search'][0]['id']
  except:
    return False

def cek_ncbi_id_by_wiki_id_via_string(param_string):    
    id_=get_wiki_id_by_name_api(param_string)
    if not id_:
        return False
    
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    query="""
        SELECT ?patogenLabel ?rankLabel ?ncbiLabel 
        WHERE {
            wd:"""+id_+""" wdt:P105 ?rank.
            wd:"""+id_+""" wdt:P685 ?ncbi.
            wd:"""+id_+""" wdt:P225 ?patogen.

            SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
        }
        """
    
    #querying
    spw=SPARQLWrapper(endpoint_url)
    spw.setQuery(query)
    spw.setReturnFormat(format_)
    hasil=spw.query().convert()
    
    if(not hasil['results']['bindings']):
        print("tidak ditemukan id: ",id_)
        return False
    
    hasil=hasil['results']['bindings'][0]
    
    if 'ncbiLabel' in hasil:        
        return (hasil['patogenLabel']['value'],
                hasil['rankLabel']['value'],
                "NCBI:"+hasil['ncbiLabel']['value'])
    else: 
        print("tidak ada id NCBI")
        return False



# untuk vector
def get_taxon_vector(param_string, endpoint_url, taxon_indonesia=True):
    # input 
    # acuan_ = 'Aphididae' (string)
    # endpoint_url = 'https://localhost:3030/mydatabase/query' (string) --> ncbi_ontology_url 
    # taxon_indonesia = True (boolean)
    # output
    # [('family', 'NCBI:27482_Aphididae'), ...]
    cek = cek_ncbi_id_by_wiki_id_via_string(param_string) #(Aphididae, family, NCBI:27482)
    if not cek:
        return False

    #langsung pake bahasa indonesia karena di graf nanti akan pake bahasa indonesia
    terjemahan={
        'superkingdom':'superkingdom',
        'kingdom':'kingdom',
        'phylum':'filum',
        'class':'kelas',
        'order':'ordo',
        'family':'famili',
        'genus':'genus',
        'species':'spesies'
    }
    data=[]
    for nama,rank,kode in getTaxonomy(cek[2],endpoint_url):
        #nama=nama.replace(" ","-")
        if rank not in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
            continue
        if (taxon_indonesia):
            data.append((terjemahan[rank],kode+'_'+nama)) #('family', 'NCBI:27482_Aphididae')
        else :
            data.append((rank,kode+'_'+nama))
    return data