'''
 NAME
      entrez_tools.py
VERSION
        6t1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Genera un diccionario con las bases de datos
        y los ids de los documentos correspondientes a la base 
        de datos dato un termino de busqueda generado automaticamente
USAGE
       python entrez_tools.py 'organism1:gene1,gene2;organism2:gene1,gene2;...'
         
ARGUMENTS
     info     information separating organism from genes by ':' genes by ','
              and organisms by ';' Example:
              organism1:gene1,gene2;organism2:gene1,gene2
        
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/entrez_tools.py

'''

# Importar los modulos necesarios
import argparse
from Bio import Entrez
import pprint

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Return IDs associated with DB \
    in a search using the entrez module")
parser.add_argument('info', metavar='info', help='information separating organism from genes by \':\'\
    genes by \',\' and organisms by \';\' Example: organism1:gene1,gene2;organism2:gene1,gene2')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
info = args.info


# Indicar e.mail
Entrez.email = 'ulisesp@lcg.unam.mx'

# Definir excepciones del programa


class FormatError(Exception):
    pass


try:
    # Funcion que genera el termino de busqueda
    def create_term(info_organism):
        '''
            Return the term to be used
            in entrez module functions
                    Parameters:
                            tuple: organism(string), genes(list)
                    Returns:
                            term ready to use in entrez module (str)
            '''
        # Agregar el nombre del organismo al termino
        term = info_organism[0] + ' AND ' + '('
        # Agregar los nombres de los genes al termino
        # agrupandolos y usando el OL OR
        for i in range(len(info_organism[1])):
            if(i < (len(info_organism[1]) - 1)):
                gene = info_organism[1][i] + ' OR '
            else:
                gene = info_organism[1][i] + ')'
            term = term + gene
        # Regresar el termino generado
        return(term)

    # Funcion get_ID_perDB

    def get_ID_perDB(termino):
        '''
        Return a dictionary containing the DB and the
        IDs associated to documents related to the search
        given a search term
            Parameters:
                termino: string
            Returns:
                dictDbID: dictionary
        '''
        # Buscar db que contienen informacion
        # del termino dado
        handle = Entrez.egquery(term=termino)
        record = Entrez.read(handle)
        handle.close()

        # Crear diccionario vacio
        dictDbId = {}
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
        # Regresar el diccionario
        return(dictDbId)

    # Crear diccionario con la informacion proveida por el usuario
    organisms = info.split(sep=';')
    dictOrganism = {}
    for organism in organisms:
        dictOrganism[organism.split(sep=':')[0]] = organism.split(sep=':')[
            1].split(sep=',')

    # Llamar a la funciones para generar los IDs
    # asociados a cada DB recorriendo el diccionario
    # De los organismos
    for item in dictOrganism.items():
        termino = create_term(item)
        dictIDs = get_ID_perDB(termino)
        if dictIDs == {}:
            raise FormatError('You probably entered the info\
 with an incorrect format')
        print(f'\nDict of {item[0]}')
        pprint.pprint(dictIDs)

# Verificar si existen excepciones
except FormatError as ex:
    print('FormatError: ' + ex.args[0])
