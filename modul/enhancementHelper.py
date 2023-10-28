import os

from modul.vectorReferenced import get_taxon_vector,cek_ncbi_id_by_wiki_id_via_string
from modul.filterNodeEdge import removeNodeAndEdgeByFilter,removeEdgesNotInNodes
from modul.preprocess import cleaning, splitInteractionToNodeEdge
from modul.disambiguation_optimized import buat_kamus_kosong, update_kamus_pake_wikidata, update_df_pake_kamus, update_df_pake_path_ujung, removeOtherThanNCBI
from modul.standardization_usingsparql import addTaxonColumn, buat_kolom_taxon_awal, getTaxonomy, getDescendant
import requests
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON, N3
from modul.disambiguation_optimized import __querying
from modul.helper_umum import pemecah_generator

def getDFMusuhAlami(search, url_ncbi_endpoint):
    # parameter
    # search --> string nama serangga: bemisia, atau kode ncbi : NCBI:7038
    # url_ncbi_endpoint='http://localhost:3030/mydataset/query' --> string url endpoint sparql ncbi
    # output
    # df_node, df_edge
    # none,none --> jika tidak ada hasil

    #2 ambil data dari globi
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
    interactionType="hostOf"
    link="https://api.globalbioticinteractions.org/interaction?sourceTaxon="+search+"&interactionType="+interactionType+"&fields="+(','.join(kolom))
    response = requests.get(link)
    if response.status_code==500:
        print('ADA ERROR DARI SERVER GLOBI')
    res=response.json()
    if not res['data']:
        return None,None
    # JSON To Pandas Dataframe
    df_serangga_hasil = pd.json_normalize(res, record_path =['data'])
    df_serangga_hasil.columns=kolom

    # split interaksi
    df_node, df_edge = splitInteractionToNodeEdge(df_serangga_hasil)
    # # df_edge.interaction_type.unique() # tipe interaksi cuma dua hostOf dan hasPathogen sama2 objeknya patogen dari serangga bemisia
    # # df_edge.source_taxon_id.unique() # node bemisia semua yang ada di source

    #4 cleaning_after_get layer 1 interaksi virus
    df_node.drop_duplicates(inplace=True)
    df_node.reset_index(drop=True,inplace=True)
    no_ncbi_and_path_null=(df_node.taxon_id.str.contains('NCBI')==False) & (df_node.taxon_path_ids.isnull())
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[no_ncbi_and_path_null], df_node,df_edge) 
    df_edge.drop_duplicates(inplace=True)
    df_edge.reset_index(drop=True,inplace=True)

    #5 disambiguasi layer 1 interaksi virus
    kamus_ncbi = buat_kamus_kosong(df_node)
    kamus_ncbi = pemecah_generator(update_kamus_pake_wikidata(kamus_ncbi))
    #update dataframe pake kamus
    df_node,df_edge = pemecah_generator(update_df_pake_kamus(kamus_ncbi,df_node,df_edge))
    df_node,df_edge = pemecah_generator(update_df_pake_path_ujung(df_node,df_edge))
    #standarisasi layer 1 interaksi virus
    df_node = buat_kolom_taxon_awal(df_node) #buat kolom taxon, default none
    df_node = addTaxonColumn(df_node,url_ncbi_endpoint) # isi pake ncbi

    # cleaning after disambiguasi
    df_node, df_edge = removeOtherThanNCBI(df_node,df_edge)# Hapus kalo masih ada selain NCBI
    df_edge = removeEdgesNotInNodes(df_node, df_edge) #hapus edge yang tidak ada nodenya

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
    
    return df_node, df_edge

def getWdId(search):
    # input
    # search = 'Bemisia tabaci'/NCBI:7038
    # output
    # wd id
    item=search
    url = f"https://www.wikidata.org/w/api.php?action=wbsearchentities&search={item}&language=en&format=json"
    data = requests.get(url).json()
    if len(data['search'])==0:
        print("tidak ditemukan entitas")
        return None
    print("ditemukan entitas sebanyak :",len(data['search']))
    print("entitas pertama :",data['search'][0])
    return data['search'][0]['id']

def _getFactFromWd(wd_id):
    # input
    # wd='Q104882'
    # output
    # dataframe
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    query="""
    SELECT ?property ?value
    WHERE {
        #dibagi jika value nya IRI atau literal, terus digabungkan lagi
        {   # IRI values
            wd:"""+wd_id+""" ?prop ?v .
            FILTER(CONTAINS(STR(?prop),"http://www.wikidata.org/prop/direct/")) # direct/truthy properties
            ?fullprop wikibase:directClaim ?prop .
            ?fullprop rdfs:label ?property . # get the label of direct properties
            FILTER(LANG(?property) = "en") # EN label of property
            FILTER(!(CONTAINS(?property,"ID"))) # get rid of identifier properties
            FILTER NOT EXISTS { ?fullprop wikibase:propertyType wikibase:ExternalId }
            FILTER(ISIRI(?v)) # IRI value
            ?v rdfs:label ?value . # label of the IRI value
            FILTER(LANG(?value) = "en") # EN label of IRI value
        } UNION
        {  # literal values
            wd:"""+wd_id+""" ?prop ?v . 
            FILTER(CONTAINS(STR(?prop),"http://www.wikidata.org/prop/direct/")) # direct/truthy properties
            ?fullprop wikibase:directClaim ?prop .
            ?fullprop rdfs:label ?property . # get the label of direct properties
            FILTER(LANG(?property) = "en") # EN label of property
            FILTER(!(CONTAINS(?property,"ID"))) # get rid of identifier properties
            FILTER NOT EXISTS { ?fullprop wikibase:propertyType wikibase:ExternalId }
            FILTER(!(ISIRI(?v))) # non-IRI value aka literal
            BIND(?v AS ?value) # variable renaming
        }
    } ORDER BY ?property ?value # sort by property then value
    """
    #querying
    spw=SPARQLWrapper(endpoint_url)
    spw.setQuery(query)
    spw.setReturnFormat(format_)
    hasil=spw.query().convert()
    if hasil['results']['bindings']==[]:
        return None
    return hasil

# tidak terpakai
def getFactDf(wd_id):
    hasil=_getFactFromWd(wd_id)
    if hasil is None:
        return None
    #df convert
    ls = []
    for h in hasil['results']['bindings']:
        ls.append((h['property']['value'],h['value']['value']))
    df=pd.DataFrame(ls, columns=['property','value'])
    return df
    # from df to dict
    # dict_fact=df_fact.groupby('property').agg(lambda x: np.unique(x).tolist())['value'].to_dict()
    # dict_fact

def getFactDict(wd_id):
    hasil=_getFactFromWd(wd_id)
    if hasil is None:
        return None
    # dict convert
    kembalian={}
    for i in hasil['results']['bindings']:
        prop = i['property']['value']
        value = i['value']['value']
        if prop in kembalian.keys():
            kembalian[prop].append(value)
        else:
            kembalian[prop]=[value]
    return kembalian


# gambar
def getPicture(wd_id):
    format_='json'
    endpoint_url='https://query.wikidata.org/sparql'
    query = """
        SELECT ?pic
        WHERE
        {
            wd:"""+wd_id+""" wdt:P18 ?pic.
        }
    """
    hasil=__querying(endpoint_url,query,format_)
    if hasil['results']['bindings']==[]:
        return None;
    if(len(hasil['results']['bindings'])==0):
        return None;

    return hasil['results']['bindings'][0]['pic']['value']

# abstak
def getAbstract(wd_id):
    format_='json'
    endpoint_url='https://dbpedia.org/sparql'
    query = """
    PREFIX dbo: <http://dbpedia.org/ontology/>
    PREFIX dbr: <http://dbpedia.org/resource/>
    PREFIX dbp: <http://dbpedia.org/property/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX wd: <http://www.wikidata.org/entity/>

    SELECT ?insect ?label ?abstract
    WHERE {
    ?insect a dbo:Insect ;
            rdfs:label ?label ;
            dbo:abstract ?abstract ;
            owl:sameAs wd:"""+wd_id+""" .
    FILTER (LANG(?label) = 'en' && LANG(?abstract) = 'en')
    }
    LIMIT 10
    """
    hasil=__querying(endpoint_url,query,format_)
    if hasil['results']['bindings']==[]:
        return None
    return hasil['results']['bindings'][0]['abstract']['value']


def getSeranggaKerabatNCBI(taxon_id, url_ncbi_endpoint):
    
    try:
        one_level_parent=getTaxonomy(taxon_id,url_ncbi_endpoint)[1][2]
        data=getDescendant(one_level_parent, url_ncbi_endpoint)
        # df = pd.DataFrame(data, columns= ['Name','Taxon','Taxon ID'])
        return data
    except:
        return None
