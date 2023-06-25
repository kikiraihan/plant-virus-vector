import pandas as pd

def extractNameClassId(obo_class):
    nama=obo_class.label.first()
    kode=obo_class.name.replace('NCBITaxon_','NCBI:')
    
    if obo_class.has_rank.first() == None:
        kelas="X"
    elif(type(obo_class.has_rank.first())==str):
        kelas=obo_class.has_rank.first().replace("http://purl.obolibrary.org/obo/NCBITaxon_","")
    else:
        kelas=obo_class.has_rank.first().label.first()
    
    return (nama, kelas, kode)



def checkClass(obo_class, obo):
    #obo.NCBITaxon_1 adalah think,
    #jika obo class tidak none, dan dia juga bukan suatu kelas yang hanya think saja.
    if obo_class is not None and obo_class.is_a.first() != obo.NCBITaxon_1:
        return True
    else:
        return False



def getTaxonomy(id_globi, obo):
    
    id_='NCBITaxon_'+id_globi.split(':')[-1]
    obo_class=getattr(obo, id_)

    if(checkClass(obo_class,obo)): 
        anc=obo_class.is_a.first()
        rank=None
        temp=[extractNameClassId(obo_class)]

        while rank!=obo.NCBITaxon_superkingdom and checkClass(anc,obo):
            temp.append(extractNameClassId(anc))
            anc=anc.is_a.first()
            rank=anc.has_rank.first()

        #print terakhir
        if checkClass(anc,obo):
            temp.append((anc.label.first(),anc.has_rank.first().label.first(),anc.name.replace('NCBITaxon_','NCBI:')))
        return temp
    else:
        return [(None,None,None)]


def addTaxonColumn(node, obo): # node adalah dataframe, obo adalah kelas ontology NCBI
    # masih per baris data
    data={}
    for idx,row in node[node["taxon_id"].str.contains("NCBI")].iterrows():
        for nama,rank,kode in getTaxonomy(row.taxon_id,obo):
            if rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
                node.loc[idx,rank]=kode+'_'+nama
    
    return node



def buat_kolom_taxon_awal(df_node):
    #buat kolom kosong
    for i in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
        df_node[i]=None
    
    #isi kolom taxon, pake data yang ada
    for idx, data in df_node.iterrows():
        if pd.isnull(data.taxon_path) | pd.isnull(data.taxon_path_ids) | pd.isnull(data.taxon_path_rank) :
            continue
        taxon_path = data.taxon_path.replace(' ','').split('|')
        taxon_path_ids = data.taxon_path_ids.replace(' ','').split('|')
        taxon_path_rank = data.taxon_path_rank.replace(' ','').split('|')
        for taxon, taxon_id, taxon_rank in zip(taxon_path,taxon_path_ids,taxon_path_rank):
            if taxon_rank in ['superkingdom','kingdom','phylum','class','order','family','genus','species']:
                df_node.loc[idx,taxon_rank] = taxon_id+"_"+taxon
    
    return df_node