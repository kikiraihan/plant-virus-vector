import networkx as nx

def bfs_relasi_virus_utama_langsung(gnx, serangga_to_search, virus_utama_ids):
    counter=0
    for edge in nx.bfs_edges(gnx.to_undirected(), source=serangga_to_search, depth_limit=2):
        s_id, o_id = edge

        s_label = gnx.nodes[s_id]['label'] +' '+s_id
        o_label = gnx.nodes[o_id]['label'] +' '+o_id

        # subjek serangga_to_search dan objek harus virus utama
        if (s_id == serangga_to_search) & (o_id in virus_utama_ids):
            # print(s_label,'-->', o_label)
            counter+=1
    
    return counter

def bfs_relasi_virus_utama_melalui_tanaman(gnx, serangga_to_search, virus_utama_ids):
    counter=0
    for edge in nx.bfs_edges(gnx.to_undirected(), source=serangga_to_search, depth_limit=2):
        s_id, o_id = edge

        s_label = gnx.nodes[s_id]['label'] +' '+s_id
        o_label = gnx.nodes[o_id]['label'] +' '+o_id
        s_group = gnx.nodes[s_id]['group']
        o_group = gnx.nodes[o_id]['group']

        # subjek serangga_to_search dan objek harus virus utama
        if (s_group == 'tanaman') & (o_id in virus_utama_ids):
            # print(s_label,'-->', o_label)
            counter+=1
    
    return counter

def bfs_relasi_virus_utama_langsung_dan_melalui_tanaman(gnx, serangga_to_search, virus_utama_ids):
    langsung=0
    melalui_tanaman=0
    for edge in nx.bfs_edges(gnx.to_undirected(), source=serangga_to_search, depth_limit=2):
        s_id, o_id = edge

        # s_label = gnx.nodes[s_id]['label'] +' '+s_id
        # o_label = gnx.nodes[o_id]['label'] +' '+o_id
        s_group = gnx.nodes[s_id]['group']
        # o_group = gnx.nodes[o_id]['group']

        # subjek serangga_to_search dan objek  virus utama
        if (s_id == serangga_to_search) & (o_id in virus_utama_ids):
            # print(s_label,'-->', o_label)
            langsung+=1

        # subjek tanaman dan objek virus utama
        if (s_group == 'tanaman') & (o_id in virus_utama_ids):
            # print(s_label,'-->', o_label)
            melalui_tanaman+=1

    return {'langsung':langsung, 'melalui_tanaman':melalui_tanaman, 'sum':langsung+melalui_tanaman}


# dari networkx dicustom
def degree_centrality_custom(G,virus_utama_ids,serangga_ids):
    """
    Parameters
    ----------
    G : graph
      A networkx graph
    Returns
    -------
    nodes : dictionary
       Dictionary of nodes with degree centrality as the value.
    """
    if len(G) <= 1:
        return {n: 1 for n in G}
    
    node_data = G.nodes(data=True)

    s = 1.0 / (len(G) - 1.0)
    # reset_n=(len(G.nodes)-1)/(len([node for node, data in G.nodes(data=True) if data.get('group') == "serangga"])-1)
    # kalo pake serangga sebagai pembagi

    centrality = {}
    for n, d in G.degree(serangga_ids):
        _relasi_virus = bfs_relasi_virus_utama_langsung_dan_melalui_tanaman(G,n,virus_utama_ids)
        if(_relasi_virus['sum']>0):
            print(node_data[n]['label'],n,'| degree:',d,'| langsung:',_relasi_virus['langsung'],'| melalui tanaman:',_relasi_virus['melalui_tanaman'])
        centrality[n] = s * d *  _relasi_virus['langsung'] * (1+_relasi_virus['melalui_tanaman']) #(CM * w1 * w2) # 

    return centrality
