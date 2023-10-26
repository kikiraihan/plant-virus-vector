import networkx as nx
from modul.embeddingHelper import df_serangga_to_rdf, rdf_KG_to_embeddings, df_to_dictionary_taxon

from modul.grafHelper import _set_networkx_graph
import json
from umap import UMAP

def get_pos_and_nx_data(df_node,df_edge):
    # ambil pos dan nx_data untuk visualisasi graf
    # input : df_node, df_edge
    # output : pos, nx_node, nx_edge
    G = _set_networkx_graph(df_node, df_edge)
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato", args="")
    # edge_trace, node_trace= _buat_tracer(G,pos)

    return {
        "status":"200",
        "pos":pos,
        "nx_node":{i:data for i,data in G.nodes(data=True)},
        "nx_edge":list(G.edges(data=True)),
    }

def get_embeddings_entities_and_dict_insect(df_node):
    # visualisasi embedding serangga
    # input : df_node

    # input : df_node, URL
    URL = "http://pyRDF2Vec"
    df_serangga = df_node[df_node['group']=="serangga"]
    CUSTOM_KG = df_serangga_to_rdf(df_serangga, URL)
    
    # embedding
    list_serangga=df_node[df_node['group']=="serangga"].taxon_id.to_list()
    list_of_entities = [ URL+"#"+taxon_id for taxon_id in list_serangga ]
    transformer, embeddings, _ = rdf_KG_to_embeddings(CUSTOM_KG, list_of_entities)
    # sama dua2nya # list_of_entities == transformer._entities 
    if (len(embeddings)<=2):
        return {
            "status":"403",
            "message":"serangga, hanya dua"
        }

    # dictionary serangga
    dictionary_serangga = df_to_dictionary_taxon(df_serangga)
    
    # UMAP
    embeddings = UMAP().fit_transform(embeddings)
    # ubah isinya ke list supaya bisa dijson
    embeddings  = [arr.tolist() for arr in embeddings]

    return {
        "status":"200",
        "embeddings":embeddings, 
        "list_of_entities":list_of_entities,
        "dictionary_serangga":dictionary_serangga
    }