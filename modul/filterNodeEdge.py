def removeNodeAndEdgeByFilter(data_terfilter,node,edge):   
    print("removeNodeAndEdgeByFilter")
    print("sebelum :",len(node),len(edge))
    # hapus edge 
    hapusEdge=[i.taxon_id for idx, i in data_terfilter.iterrows()]
    edge.drop(edge[edge.source_taxon_id.isin(hapusEdge)].index,inplace=True)
    edge.drop(edge[edge.target_taxon_id.isin(hapusEdge)].index,inplace=True)

    # hapus node 
    node.drop( data_terfilter.index,inplace=True)
    
    print("sesudah :",len(node),len(edge))
    node.reset_index(drop=True, inplace=True)
    edge.reset_index(drop=True, inplace=True)
    return node, edge

def takeNodeAndEdgeByFilter(data_terfilter,node,edge):    
    print("takeNodeAndEdgeByFilter")
    print("sebelum :",len(node),len(edge))
    # hapus edge lama
    # hapusEdge=[i.taxon_id for idx, i in data_terfilter.iterrows()]
    # edge.drop(edge[edge.source_taxon_id.isin(hapusEdge)==False].index,inplace=True)
    # edge.drop(edge[edge.target_taxon_id.isin(hapusEdge)==False].index,inplace=True)
    # hapus edge baru
    edgeAda=[i.taxon_id for idx, i in node.iterrows()]
    edge=edge[ ( edge.source_taxon_id.isin(edgeAda) & edge.target_taxon_id.isin(edgeAda) )  ]

    # hapus node fungi
    node = data_terfilter
    
    print("sesudah :",len(node),len(edge))
    node.reset_index(drop=True, inplace=True)
    edge.reset_index(drop=True, inplace=True)
    return node, edge


#delete edges that are not in nodes
def removeEdgesNotInNodes(node, edge):
    print("removeEdgesNotInNodes")
    print("sebelum : ",len(edge))
    edgeAda = [i.taxon_id for idx, i in node.iterrows()]
    edge = edge[ ( edge.source_taxon_id.isin(edgeAda) & edge.target_taxon_id.isin(edgeAda) )  ]
    
    print("sesudah : ",len(edge))
    edge.reset_index(drop=True, inplace=True)
    return edge;