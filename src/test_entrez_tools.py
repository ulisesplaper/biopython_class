'''
 NAME
      test_entrez_tools.py
VERSION
        7.2t1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Genera un diccionario con las bases de datos
        y los ids de los documentos correspondientes a la base 
        de datos dato un termino de busqueda generado automaticamente
USAGE
       python test_entrez_tools.py 'organism1:gene1,gene2;organism2:gene1,gene2;...'
         
ARGUMENTS
     info     information separating organism from genes by ':' genes by ','
              and organisms by ';' Example:
              organism1:gene1,gene2;organism2:gene1,gene2
        
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/test_entrez_tools.py
'''
# Importar los modulos necesarios
import argparse
import entrez_tools_module as et
from pprint import pprint

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Return IDs associated with DB \
    in a search using the entrez module")
parser.add_argument('info', metavar='info', help='information separating organism from genes by \':\'\
    genes by \',\' and organisms by \';\' Example: organism1:gene1,gene2;organism2:gene1,gene2')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
info = args.info

# Llamar a las funciones del modulo entrez_tools_module.py
pprint(et.get_ID_perDB(et.create_term(info)))



