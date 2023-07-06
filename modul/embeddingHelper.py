from pyrdf2vec import RDF2VecTransformer
from pyrdf2vec.graphs import KG, Vertex
from pyrdf2vec.embedders import FastText,Word2Vec
from pyrdf2vec.walkers import RandomWalker


def nxnode_to_rdf(gnx, URL, data_acuan=None):
    #konversi node networkx ke RDF
    # parameter ------------------------------------------------------------------------------------------------
    # df_serangga adalah dataframe yang berisi data serangga
    # URL adalah URL dari KG
    # data acuan
    # output ---------------------------------------------------------------------------------------------------
    # KG() . rdf kg
    #------------------------------------------------------------------------------------------------------------

    CUSTOM_KG = KG()
    keindo={
        'superkingdom':'superkingdom',
        'kingdom':'kingdom',
        'phylum':'filum',
        'class':'kelas',
        'order':'ordo',
        'family':'famili',
        'genus':'genus',
        'species':'spesies'
    }
    keinggris={
        'superkingdom':'superkingdom',
        'kingdom':'kingdom',
        'filum':'phylum',
        'kelas':'class',
        'ordo':'order',
        'famili':'family',
        'genus':'genus',
        'spesies':'species',
    }
    takson_to_embedd=[
        #  'superkingdom','kingdom','filum','kelas',
        'ordo','famili','genus','spesies'
    ];

    # kalo ada data acuan
    if(data_acuan is not None):
        # kalau takson_to_embedd not in data_acuan
        listnya=[key for key, value in data_acuan]
        takson_to_embedd=[x for x in takson_to_embedd if x in listnya]
        
        # memasukan RDF serangga acuan
        subj = Vertex(f"{URL}#SERANGGA_ACUAN")
        for i,j in data_acuan:
            if(i in takson_to_embedd):
                j = j.replace(' ','-')
                obj = Vertex((URL+"#"+j))
                pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
                #pred = Vertex((URL+"#taxon_path_ids"), predicate=True, vprev=subj, vnext=obj)
                CUSTOM_KG.add_walk(subj, pred, obj)

    # proses konversi 
    for index,data in gnx.nodes(data=True):
        if(data['group']=='serangga'): #jika serangga
            subj = Vertex(URL+"#"+index)
            for i in takson_to_embedd:
                    id_takson=data[keindo[i]].replace(' ','-')#.split('_')[0]
                    obj = Vertex((URL+"#"+id_takson))
                    pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
                    CUSTOM_KG.add_walk(subj, pred, obj)
    # CUSTOM_KG.literals=[
    #         [f"{URL}#taxon_path_ids"],
    #     ]
    CUSTOM_KG.literals = [[URL+"#"+i] for i in takson_to_embedd]


def df_serangga_to_rdf(df_serangga, URL, data_acuan=None):
    #konversi dataframe serangga ke RDF KG
    # parameter ------------------------------------------------------------------------------------------------
    # df_serangga adalah dataframe yang berisi data serangga
    # URL adalah URL dari KG
    # data acuan
    # output ---------------------------------------------------------------------------------------------------
    # KG() . rdf kg
    #------------------------------------------------------------------------------------------------------------
    
    CUSTOM_KG = KG()
    takson_to_embedd=[
        #  'superkingdom', 'kingdom', 'phylum','class', 
        'order', 'family', 'genus', 'species'
    ];

    # kalo ada data acuan
    if(data_acuan is not None):
        # kalau takson_to_embedd not in data_acuan
        listnya=[key for key, value in data_acuan]
        takson_to_embedd=[x for x in takson_to_embedd if x in listnya]
        
        # memasukan RDF serangga acuan
        subj = Vertex(f"{URL}#SERANGGA_ACUAN")
        for i,j in data_acuan:
            if(i in takson_to_embedd):
                j = j.replace(' ','-')
                obj = Vertex((URL+"#"+j))
                pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
                #pred = Vertex((URL+"#taxon_path_ids"), predicate=True, vprev=subj, vnext=obj)
                CUSTOM_KG.add_walk(subj, pred, obj)

    # proses konversi 
    for index,data in df_serangga.iterrows():
        # if(data['group']=='serangga'): #jika serangga
        subj = Vertex(URL+"#"+data.taxon_id)
        for i in takson_to_embedd:
            #if(isinstance(data[i], str)): #jika dia string atau tidak nan/kosong.
            #if(i not in ["superkingdom","kingdom","filum","kelas"]):
                id_takson=data[i].replace(' ','-')#.split('_')[0]
                obj = Vertex((URL+"#"+id_takson))
                pred = Vertex((URL+"#"+i), predicate=True, vprev=subj, vnext=obj)
                #pred = Vertex((URL+"#taxon_path_ids"), predicate=True, vprev=subj, vnext=obj)
                CUSTOM_KG.add_walk(subj, pred, obj)
    # CUSTOM_KG.literals=[
        #         [f"{URL}#taxon_path_ids"],
        #     ]
    CUSTOM_KG.literals = [[URL+"#"+i] for i in takson_to_embedd]
    return CUSTOM_KG




def rdf_KG_to_embeddings(CUSTOM_KG, list_entity):
    # embedding
    # parameter ------------------------------------------------------------------------------------------------
    #  custom_kg : KG() dari konversi node networkx ke RDF
    #  URL : URL dari KG
    #  list_id_serangga : list id serangga yang akan diembedd # ['NCBI:355109',...,'NCBI:262664','NCBI:355026',]
    # output ---------------------------------------------------------------------------------------------------
    # transformer. # untuk mengambil transformer.entities # transformer.entities sama dengan ent
    # embeddings, 
    # _
    #------------------------------------------------------------------------------------------------------------

    # Ensure the determinism of this script by initializing a pseudo-random number.
    RANDOM_STATE = 22
    transformer = RDF2VecTransformer(
        # Use one worker threads for Word2Vec to ensure random determinism.
        # Must be used with PYTHONHASHSEED.
        Word2Vec(epochs=1000),
        # Extract a maximum of 10 walks of a maximum depth of 4 for each entity
        # using two processes and use a random state to ensure that the same walks
        # are generated for the entities.
        walkers=[RandomWalker(2, 5, n_jobs=2, with_reverse=False, random_state=RANDOM_STATE)],
        #verbose=1,
    )
    # transformer = RDF2VecTransformer(verbose=1)
    # Fit the transformer to the knowledge graph and the entities.
    embeddings, _ = transformer.fit_transform(
        CUSTOM_KG, #the KG
        list_entity, #entity
    )

    return transformer, embeddings, _