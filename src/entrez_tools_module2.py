'''
 NAME
      entrez_tools_module2.py
VERSION
        8t1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Funcion 1: extrae los ids asociados a db de articulos
        Funcion 2: extrae los titulos y abstracts de un id y
        una db especificados.
        Funcion 3: extrae ids de articulos que han citado 
        al articulo dado
USAGE
       Diseniado para usarse como modulo
         
ARGUMENTS
     
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/entrez_tools_module_2.py
'''

# Importar modulos necesarios
from Bio import Entrez, Medline
from pprint import pprint as pt
from io import StringIO

# Indicar e.mail
Entrez.email = 'ulisesp@lcg.unam.mx'

# Funcion para obtner los ids asociados a articulos
def get_ids(dic_org):
    '''
    Returns a dictionary with ids and titles of articles
        Parameters:    
            dict_org (dic): dictionary with {Orgn: {db:[ids]}} 
        Returns:
            tuple of dicts
    '''
    # Inicializar diccionarios de ids y titulos
    id_article = {}
    title_articles = {}
    # Recorrer el diccionario de entrada por organismo
    # y posteriormente los db para acceder a los ids
    for organism, dbs in dic_org.items():
        dic_db = dbs
        dbID = {}
        db_articles = {}
        for db,ids in dic_db.items():
            # Tomar los ids que corresponden a las db con articulos
            if ((db == 'pubmed') or (db == 'pmc')):
                dbID[db] = ids
                # Obtner los titulos
                handle = Entrez.efetch(db = db, id= ids, rettype='medline', retmode='text')
                # Leer la info
                data = handle.read()
                handle.close()
                # Generar una lista manejable
                info_list = data.split(sep='\n\n')
                # Lista donde se almacenan los titles
                title_list = []
                count = 0
                # Generar los titulos por id
                for id in ids:
                    title = info_list[count].split(sep='\nTI  - ')[1].split('.')[0].replace('\n      ', '')
                    count = count + 1
                    title_list.append(title)
        # Acceder la informacion a los diccionarios 
                db_articles[db] = title_list      
        id_article[organism] = dbID
        title_articles[organism] = db_articles
        
    # Retornar la tupla        
    return(id_article,title_articles)
    
    
    ## Funcion que extrae los titulos y abstracts
    ## Funcion que extrae los titulos
def get_titles_abstracts(db, id):
    '''
    Returns a list containing titles and abstract
        Parameters:    
            db (str): db to look for in
            id(str): id of article 
        Returns:
            list of title and abstract
    '''
    # Buscar en las db y leer los datos
    handle = Entrez.efetch(db=db, id = id,  retmode = 'text' ,rettype = 'medline') 
    record = handle.read()
    rec_file = StringIO(record)
    medline_rec = Medline.read(rec_file)
    # Preguntar por la informacion y accesarla a la lista
    if 'AB' in medline_rec:
        TIAB_list = [medline_rec['TI'], medline_rec['AB']]
    else:
        TIAB_list = [medline_rec['TI']]
    # Retornar el valor
    return(TIAB_list)
    
    # Funcion que extrae los ids que han citado el articulo dado  
def get_citedby_ids(id_artdb_dic):
    '''
    Returns a dictionary with ids of citing articles
        Parameters:    
            id_artdb_dic (dic): dictionary with {Orgn: {db_articles:[ids]}} 
        Returns:
            dict of citing articles
    '''
    art_cited_dict = {}
    for organism in id_artdb_dic.keys():
        # Extraer los articulos citados en pubmed
        dbs = ['pubmed','pmc']
        
        for db in dbs:
            # Asignar el linkname
            if (db == 'pubmed'):
                linkname = f'{db}_{db}'
            else:
                linkname = f'{db}_{db}_cites'
            # Buscar los ids de los articulos citantes
            ids = id_artdb_dic[organism][db]
            results = Entrez.read(Entrez.elink(dbfrom=db, db=db, from_uid=ids, linkname = linkname))
                
            # Recorrer la lista de los resultados por ID
            indv_art_dict = {}
            
            for i,result in enumerate(results):
                #Verificar que el id tenga articulos que lo hayan citado
                # Accesar los ids en tal caso
                if ((result["LinkSetDb"])):
                    data_list = []
                    for link in result["LinkSetDb"][0]["Link"]:
                        data_list.append(link['Id'])

                    # Llamar a la funcion 
                else:
                    data_list = []
                # Asociar los articulos citantes al citado
                indv_art_dict[ids[i]] = data_list
            # Reasigamos el linkname de acuerdo con la db
        # Asociar el organismo a sus ids    
        art_cited_dict[organism] = indv_art_dict
    # Regresar el diccionario
    return(art_cited_dict)