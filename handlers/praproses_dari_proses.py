
import pandas as pd
from modul.filterNodeEdge import removeNodeAndEdgeByFilter
from modul.helper_umum import contains_string_entire_column


def pra_proses_dari_proses(df_node,df_edge, acuan_=None,tambah_relasi_dengan_acuan=False):
    # pra-proses khusus proses
    # hapus serangga yg cuma famili (mengikuti acuan). soalnya klo cuma tampil famili apa gunanya?
    filter_genus_sampai_species_null=(
        (df_node.genus.isnull()) &
        (df_node.species.isnull()) &
        (df_node['class']=='NCBI:50557_Insecta')
    )
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[filter_genus_sampai_species_null], df_node,df_edge)

    df_edge.drop_duplicates(inplace=True)

    # pra-proses khusus proses
    #isi data kosong. mengisi takson kosong, dengan takson sebelumnya, untuk tambalan
    takson=[
        'superkingdom','kingdom','phylum','class','order','family','genus','species'
    ]
    for x,i in enumerate(takson):
        if (i!='superkingdom'): #selain superkingdom update dengan data sebelumnya
            for idx, row in df_node[pd.isnull(df_node[i])].iterrows():
                df_node.loc[idx,[i]] = row[takson[x-1]]+'^'+i
        else: 
            for idx, row in df_node[pd.isnull(df_node[i])].iterrows():
                df_node.loc[idx,[i]] = row[takson[x+1]]+'^'+i

    if (tambah_relasi_dengan_acuan):
        # pra-proses khusus proses
        # eksperimen tambahan. bikin fakta tambahan, yaitu relasi virus utama dengan acuan 
        # virus_utama=df_node[df_node.virus_utama==True].taxon_id.to_list()
        # serangga_acuan=contains_string_entire_column(df_node,acuan_).taxon_id.to_list()
        # print(len(df_edge))
        # for i in virus_utama:
        #     for j in serangga_acuan:
        #         dict = {'source_taxon_id':i,'target_taxon_id':j,'interaction_type':'pathogenOf'}
        #         df_edge = pd.concat([pd.DataFrame(dict,index=[0]), df_edge], ignore_index = True)
        #         # df_edge.loc[len(df_edge.index),['source_taxon_id','target_taxon_id','interaction_type']] = [i,j,'pathogenOf']
        # print(len(df_edge))

        #satu saja
        virus_utama=df_node[df_node.virus_utama==True].taxon_id.to_list()
        serangga_acuan=contains_string_entire_column(df_node,acuan_).taxon_id.to_list()
        print(len(df_edge))
        for j in serangga_acuan:
            dict = {'source_taxon_id':virus_utama[0],'target_taxon_id':j,'interaction_type':'pathogenOf'}
            df_edge = pd.concat([pd.DataFrame(dict,index=[0]), df_edge], ignore_index = True)
            # df_edge.loc[len(df_edge.index),['source_taxon_id','target_taxon_id','interaction_type']] = [i,j,'pathogenOf']
        print(len(df_edge))

    # Ini harusnya di praproses. tapi belum fix baiknya hapus atau tidak
    # hapus yang bukan virus utama, terakhir akurasi 0.85, kalau berkurang hapus saja ini
    bukan_virus_utama=(df_node['group']=="virus") & (df_node.virus_utama!=True)
    df_node,df_edge = removeNodeAndEdgeByFilter(df_node[bukan_virus_utama], df_node,df_edge)

    print('praproses dari proses selesai')
    return df_node,df_edge
