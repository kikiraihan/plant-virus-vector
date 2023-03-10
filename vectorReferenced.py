from standardization import getTaxonomy
from SPARQLWrapper import SPARQLWrapper
import requests

# pencarian ID berdasarkan nama entity
def call_wiki_api(item):
  try:
    url = f"https://www.wikidata.org/w/api.php?action=wbsearchentities&search={item}&language=en&format=json"
    data = requests.get(url).json()
    # Return the first id (Could upgrade this in the future)
    return data['search'][0]['id']
  except:
    return False

def cek_vector_wiki_id(acuan_):    
    id_=call_wiki_api(acuan_)
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
def get_taxon_vector(acuan_, obo):
    cek = cek_vector_wiki_id(acuan_) #(Aphididae, family, NCBI:27482)
    if not cek:
        return False

    #langsung pake bahasa indonesia karena di graf nanti akan pake bahasa indonesia
    taxon_dict={
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
    for nama,rank,kode in getTaxonomy(cek[2],obo):
        #nama=nama.replace(" ","-")
        if rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
            data.append((taxon_dict[rank],kode+'_'+nama)) #('family', 'NCBI:27482_Aphididae')
    return data