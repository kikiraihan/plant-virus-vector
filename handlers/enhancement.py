from modul.enhancementHelper import getDFMusuhAlami, getWdId, getFactDict, getPicture, getAbstract, getSeranggaKerabatNCBI

def get_musuh_alami_handler(search):
    df_node, df_edge = getDFMusuhAlami(search)
    if(df_node is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan musuh alami',
        }
    return {
        'status':200,
        'message': 'Kirim json', 
        'node': df_node.to_json(orient="split"), 
        'edge': df_edge.to_json(orient="split"), 
    }

def get_fact_handler(wd_id):# NCBI:7032 atau #bemisia tabaci
    fact = getFactDict(wd_id)
    if(fact is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan fact box',
        }
    return {
        'status': 200,
        'message': 'Kirim json',
        'fact': fact,
    }

def get_picture_handler(wd_id):
    picture = getPicture(wd_id)
    if(picture is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan gambar',
        }
    
    return {
        'status': 200,
        'message': 'Kirim json',
        'picture': picture,
    }

def get_abstract_handler(wd_id):
    abstract = getAbstract(wd_id)
    if(abstract is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan abstrak',
        }
    
    return {
        'status': 200,
        'message': 'Kirim json',
        'abstract': abstract,
    }


def get_wd_id_handler(search):
    wd_id = getWdId(search)
    if(wd_id is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan entitas',
        }
    
    return {
        'status': 200,
        'message': 'Kirim json',
        'wd_id': wd_id,
    }
    

def get_relatives_handler(taxon_id):
    relatives = getSeranggaKerabatNCBI(taxon_id)
    if(relatives is None):
        return {
            'status': 404,
            'message': 'Tidak ditemukan entitas',
        }
    
    return {
        'status': 200,
        'message': 'Kirim json',
        'relatives': relatives,
    }