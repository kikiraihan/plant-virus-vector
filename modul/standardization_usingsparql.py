import os

from SPARQLWrapper import SPARQLWrapper
import pandas as pd

JENA_URL = os.environ.get("JENA_URL")
JENA_URL_MAINDB = os.environ.get("JENA_URL_MAINDB")

def getDescendant(id_globi, endpoint_url):
	id_ncbi=id_globi.split(':')[-1]
	query="""
	PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX owl: <http://www.w3.org/2002/07/owl#>
	PREFIX obo: <http://purl.obolibrary.org/obo/>
	PREFIX obo2: <http://purl.obolibrary.org/obo/ncbitaxon#>

	SELECT ?Objlabel ?Ranklabel (REPLACE(STR(?obj), "http://purl.obolibrary.org/obo/NCBITaxon_", "NCBI:") AS ?objWithoutPrefix) WHERE {
	?obj rdfs:subClassOf* obo:NCBITaxon_"""+id_ncbi+""".
	?obj rdfs:label ?Objlabel.
	?obj obo2:has_rank ?rank.
	?rank rdfs:label ?Ranklabel.
	}
	"""
	spw=SPARQLWrapper(endpoint_url)
	spw.setQuery(query)
	spw.setReturnFormat('json')
	hasil=spw.query().convert()

	kembalian=[
		(
			i['Objlabel']['value'],
			i['Ranklabel']['value'],
			i['objWithoutPrefix']['value']
		) 
		for i in hasil['results']['bindings']
	]
	
	return kembalian

def getTaxonomy(id_globi, endpoint_url):
	id_ncbi=id_globi.split(':')[-1].strip()
	query="""
	PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX owl: <http://www.w3.org/2002/07/owl#>
	PREFIX obo: <http://purl.obolibrary.org/obo/>
	PREFIX obo2: <http://purl.obolibrary.org/obo/ncbitaxon#>

	SELECT ?Objlabel ?Ranklabel (REPLACE(STR(?obj), "http://purl.obolibrary.org/obo/NCBITaxon_", "NCBI:") AS ?objWithoutPrefix) WHERE {
	obo:NCBITaxon_"""+id_ncbi+""" rdfs:subClassOf* ?obj.
	?obj rdfs:label ?Objlabel.
	?obj obo2:has_rank ?rank.
	?rank rdfs:label ?Ranklabel.
	}
	"""
	spw=SPARQLWrapper(endpoint_url)
	spw.setQuery(query)
	spw.setReturnFormat('json')
	hasil=spw.query().convert()

	kembalian=[
		(
			i['Objlabel']['value'],
			i['Ranklabel']['value'],
			i['objWithoutPrefix']['value']
		) 
		for i in hasil['results']['bindings']
	]
	
	return kembalian

def addTaxonColumn(node, endpoint_url = f"{JENA_URL_MAINDB}/query"): # node adalah dataframe, obo adalah kelas ontology NCBI
	# masih per baris data
	data={}
	for idx,row in node[node["taxon_id"].str.contains("NCBI")].iterrows():
		for nama,rank,kode in getTaxonomy(row.taxon_id, endpoint_url):
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