'''
 NAME
      test_entrez_tools2.py
VERSION
        8.2t1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
     Prueba las funciones del modulo entrez_tools_module2.py   
USAGE
    python test_entrez_tools.py 'organism1:gene1,gene2;organism2:gene1,gene2;...'
         
ARGUMENTS
     info     information separating organism from genes by ':' genes by ','
              and organisms by ';' Example:
              organism1:gene1,gene2;organism2:gene1,gene2
        
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/test_entrez_tools2.py
'''
# Importar los modulos necesarios
import argparse
import entrez_tools_module as et
import entrez_tools_module2 as et2
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
ids = (et.get_ID_perDB(et.create_term(info)))

# Obtener los ids y titulos de los ids asociados a articulos
print("Ids asociados a articulos\n")
print(et2.get_ids(ids))

# Obtener los ids de los articulos que los citan
print('\narticulos citantes:\n')
citingArt = et2.get_citedby_ids(ids)
pprint(citingArt)

# Imprimir abstracts y titulos de articulos citantes
print("\nTitulos y abstracts\n")
for keys, values in citingArt.items():
    for secID in values:
        pprint(et2.get_titles_abstracts('pubmed', secID))