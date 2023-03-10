import pandas as pd

# kalo ba update disini, bken juga di praproses tanaman
def cleaning(data, dropNonNCBI=False):
    
    df_tan=data.copy()
    
    #drop duplicate
    print('drop duplicate')
    print('sebelum : ',len(df_tan))
    df_tan.drop_duplicates(inplace=True)
    print('setelah: ', len(df_tan))
    
    #jika taxon_external_id tidak mengandung NCBI dan taxon_path_id-nya null
    print('drop tidak mengandung NCBI dan taxon pathnya null')
    print('sebelum : ',len(df_tan))
    hapus=[i for i,d in df_tan[
        ((df_tan.source_taxon_external_id.str.contains('NCBI')==False) & (df_tan.source_taxon_path_ids.isnull())) |
        ((df_tan.target_taxon_external_id.str.contains('NCBI')==False) & (df_tan.target_taxon_path_ids.isnull()))
    ].iterrows()]
    df_tan.drop(hapus, inplace=True)
    print('setelah: ', len(df_tan))

    #dropNonNCBI
    if(dropNonNCBI):
        print('drop non NCBI')
        print('sebelum : ',len(df_tan))
        #drop target yang punya data selain NCBI
        df_tan=df_tan.query('target_taxon_external_id.str.contains("NCBI")', engine='python')
        #drop source yang punya data selain NCBI
        df_tan=df_tan.query('source_taxon_external_id.str.contains("NCBI")', engine='python')
        print('setelah: ', len(df_tan))


    # MODIF
    # BERSIHKAN DUPLIKASI PER NAMA TAPI BEDA ID.
    #semua nama duplikat, di kumpulkan datanya per nama.
    to_node={}
    for i,a in df_tan[[
        'source_taxon_external_id',
        'source_taxon_name',
        'source_taxon_path',
        'source_taxon_path_ids',
        'source_taxon_path_ranks',
    ]].drop_duplicates().iterrows():
        if a.source_taxon_name in to_node:
            #utamakan NCBI
            if 'NCBI' in a['source_taxon_external_id']:
                to_node[a.source_taxon_name][0]=a
            else:
                to_node[a.source_taxon_name].append(a)
        else:
            to_node[a.source_taxon_name]=[a]

    to_node_pat={}
    for i,a in df_tan[[
        'target_taxon_external_id',
        'target_taxon_name',
        'target_taxon_path',
        'target_taxon_path_ids',
        'target_taxon_path_ranks',
    ]].drop_duplicates().iterrows():
        if a.target_taxon_name in to_node_pat:
            #utamakan NCBI
            if 'NCBI' in a['target_taxon_external_id']:
                to_node_pat[a.target_taxon_name][0]=a
            else:
                to_node_pat[a.target_taxon_name].append(a)
        else:
            to_node_pat[a.target_taxon_name]=[a]


    #source : update duplikasi dengan nilai pertama
    df_tan[[
        'source_taxon_external_id',
        'source_taxon_name',
        'source_taxon_path',
        'source_taxon_path_ids',
        'source_taxon_path_ranks',
       ]]=df_tan[[
        'source_taxon_external_id',
        'source_taxon_name',
        'source_taxon_path',
        'source_taxon_path_ids',
        'source_taxon_path_ranks',
       ]].apply(lambda x: to_node[x.source_taxon_name][0], axis=1)

    #target : update duplikasi dengan nilai pertama
    df_tan[[
        'target_taxon_external_id',
        'target_taxon_name',
        'target_taxon_path',
        'target_taxon_path_ids',
        'target_taxon_path_ranks',
       ]]=df_tan[[
        'target_taxon_external_id',
        'target_taxon_name',
        'target_taxon_path',
        'target_taxon_path_ids',
        'target_taxon_path_ranks',
       ]].apply(lambda x: to_node_pat[x.target_taxon_name][0], axis=1)

    #masih ada duplikasi, jadi hapus
    print('drop duplicate nama sama beda id')
    print('sebelum : ',len(df_tan))
    df_tan.drop_duplicates(inplace=True)
    print('setelah: ', len(df_tan))

    #reset index
    df_tan = df_tan.reset_index(drop=True)
    

    return df_tan







def removeNodeRootEntity(node,edge):
    #hapus kalo cuma root
    node=node[node['taxon_id']!='NCBI:1']
    edge=edge[edge['source_taxon_id']!='NCBI:1']
    edge=edge[edge['target_taxon_id']!='NCBI:1']
    
    return node,edge






def splitInteractionToNodeEdge(df):
    #pemisahan
    df_source=df[[
    'source_taxon_external_id',
    'source_taxon_name',
    'source_taxon_path',
    'source_taxon_path_ids',
    'source_taxon_path_ranks',
    ]].drop_duplicates().reset_index(drop=True).copy()
    df_target=df[[
        'target_taxon_external_id',
        'target_taxon_name',
        'target_taxon_path',
        'target_taxon_path_ids',
        'target_taxon_path_ranks',
    ]].drop_duplicates().reset_index(drop=True).copy()
    df_edge=df[[
        'source_taxon_external_id',
        'target_taxon_external_id',
        'interaction_type'
    ]].drop_duplicates().reset_index(drop=True).copy()

    #rename kolom
    df_source.rename(columns={
        'source_taxon_external_id': 'taxon_id',
        'source_taxon_name': 'taxon_name',
        'source_taxon_path': 'taxon_path',
        'source_taxon_path_ids' :'taxon_path_ids',
        'source_taxon_path_ranks' :'taxon_path_rank',
    },inplace=True)
    df_target.rename(columns={
        'target_taxon_external_id': 'taxon_id',
        'target_taxon_name': 'taxon_name',
        'target_taxon_path': 'taxon_path',
        'target_taxon_path_ids' :'taxon_path_ids',
        'target_taxon_path_ranks' :'taxon_path_rank',
    },inplace=True)
    df_edge.rename(columns={
        'source_taxon_external_id': 'source_taxon_id',
        'target_taxon_external_id': 'target_taxon_id',
        'interaction_type': 'interaction_type',    
    },inplace=True)

    #gabung jadi node
    df_node=pd.concat([df_target,df_source]).drop_duplicates()

    #hapus kalo cuma root
    df_node,df_edge = removeNodeRootEntity(df_node,df_edge);

    #reset index
    df_node=df_node.reset_index(drop=True)
    
    return df_node,df_edge

