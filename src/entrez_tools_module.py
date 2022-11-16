'''
 NAME
      entrez_tools.py
VERSION
        6t1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Funcion 1: genera un termino a usar en busquedas del 
        modulo Entrez
        Funcion 2. Genera un diccionario con las bases de datos
        y los ids de los documentos obtenidos a partir de un 
        un termino de busqueda generado automaticamente
USAGE
       Diseniado para usarse como modulo
         
ARGUMENTS
     
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/entrez_tools_module.py
'''

# Importar modulos necesarios
from Bio import Entrez
import pprint

# Indicar a Entrez el e.mail.
Entrez.email = 'ulisesp@lcg.unam.mx'


def create_term(info_organism):
    '''
    Returns a term to make a search with Entrez module
        Parameters
            info_organism (str): string with organism 
        Returns:
            term(str or list of strs): search term
    '''
    # Verificar la presencia del caracter de separacion de organismos
    if ';' in info_organism:
        # Crear una lista a partir del caracter de separacion
        # recorrer la lista e ir generando terminos de busqueda
        organisms = info_organism.split(sep=';')
        term_list = []
        for organism in organisms:
            # Agregar el nombre del organismo al termino
            term = organism.split(sep=':')[0] + '[orgn]' + ' AND ' + '('
            gene_list = organism.split(sep=':')[1].split(',')
            # Agregar los nombres de los genes al termino
            # agrupandolos y usando el OL OR
            for i in range(len(gene_list)):
                if(i < (len(gene_list) - 1)):
                    gene = gene_list[i] + ' OR '
                else:
                    gene = gene_list[i] + ')'
                term = term + gene
            term_list.append(term)
        return(term_list)
    # Si solo se indico un organismo, no generar la lista
    else:
        organism = info_organism
        # Agregar el nombre del organismo al termino
        term = organism.split(sep=':')[0] + '[orgn]' + ' AND ' + '('
        gene_list = organism.split(sep=':')[1].split(',')
        # Agregar los nombres de los genes al termino
        # agrupandolos y usando el OL OR
        for i in range(len(gene_list)):
            if(i < (len(gene_list) - 1)):
                gene = gene_list[i] + ' OR '
            else:
                gene = gene_list[i] + ')'
            term = term + gene
        # Regresar el termino generado
        return(term)


def get_ID_perDB(terminos):
    '''
    Return a dictionary containing the DB and the
    IDs associated to documents related to the search
    given a search term
        Parameters:
            terminos(list): search term
        Returns:
            organism_dic(dict): ids associated to a search term
    '''
    # Crear diccionario vacio
    organism_dic = {}
    # Recorrer los terminos de busqueda
    for termino in terminos:
        dictDbId = {}
        # Buscar db que contienen informacion
        # del termino dado
        handle = Entrez.egquery(term=termino)
        record = Entrez.read(handle)
        handle.close()

        # Recorrer la lista de resultados
        # de egquery y extraer las bases de datos y los conteos
        # de los match encontrados
        # Excluir 0s y Errores
        for db in record['eGQueryResult']:
            dbase = (db['DbName'])
            count = db['Count']
            if((count == 'Error') or (count == '0')):
                continue

            # Por cada base de datos extraer la lista de IDs
            # Agregar la db y los ids a un diccionario
            handle = Entrez.esearch(db=dbase, term=termino, retmax=count)
            db_record = Entrez.read(handle)
            handle.close()
            dictDbId[dbase] = db_record['IdList']
        organism_dic['\n' + termino.split('[orgn]')[0] + '\n'] = dictDbId
    return(organism_dic)