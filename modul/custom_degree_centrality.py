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


def _relasi_langsung_dan_tidak_langsung(gnx, virus_utama_ids, serangga_ids):
    # buat dict serangga yang berelasi langsung dan tidak langsung dengan virus
    serangga_langsung=dict(zip(serangga_ids, [0]*len(serangga_ids)))
    serangga_tidak_langsung=dict(zip(serangga_ids, [0]*len(serangga_ids)))
    # print(serangga_tidak_langsung)

    for virus_source in virus_utama_ids:
        for s_id, o_id in nx.bfs_edges(gnx.to_undirected(), source=virus_source, depth_limit=2):
            s_group = gnx.nodes[s_id]['group']
            o_group = gnx.nodes[o_id]['group']

            if (s_group=='virus' and o_group=='serangga'):
                serangga_langsung[o_id] += 1
            elif (s_group=='tanaman' and o_group=='serangga'):
                serangga_tidak_langsung[o_id] += 1

    return serangga_langsung, serangga_tidak_langsung


# dari networkx dicustom
def degree_centrality_custom(G,virus_utama_ids,serangga_ids, print_relasi=True):
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

    # asli
    # s = 1.0 / (len(G) - 1.0)
    # ketika pakai serangga sebagai pembagi
    # s =(len(G.nodes)-1)/(len([node for node, data in G.nodes(data=True) if data.get('group') == "serangga"])-1)
    # pembaruan karena ada bobot juga
    s = 1.0 / ( ((len(G) - 1.0)**2) * len(G) )
    # kalo pake serangga sebagai pembagi
    
    # cek semua serangga berelasi langsung dan tidak langsung, dari sudut pandang virus
    _langsung, _tidak_langsung = _relasi_langsung_dan_tidak_langsung(G, virus_utama_ids, serangga_ids)

    centrality = {}
    for n, d in G.degree(serangga_ids):
        # _relasi_virus = bfs_relasi_virus_utama_langsung_dan_melalui_tanaman(G,n,virus_utama_ids)
        # if(_relasi_virus['sum']>0 and print_relasi):
        #     print(node_data[n]['label'],n,'| degree:',d,'| langsung:',_relasi_virus['langsung'],'| melalui tanaman:',_relasi_virus['melalui_tanaman'])
        # centrality[n] = s * d *  _relasi_virus['langsung'] * (1+_relasi_virus['melalui_tanaman']) #(CM * w1 * w2) # 

        if( (_langsung[n]>0 or _tidak_langsung[n]>0) and print_relasi):
            print(node_data[n]['label'],n,'| degree:',d,'| langsung:',_langsung[n],'| melalui tanaman:',_tidak_langsung[n])
        centrality[n] = s * d *  _langsung[n] * (1+_tidak_langsung[n]) #(CM * w1 * w2) # 

    return centrality
