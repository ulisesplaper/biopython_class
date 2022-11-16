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
        https://github.com/ulisesplaper/python_class/blob/master/src/entrez_tools.py
'''

# Importar modulos necesarios

# Indicar a Entrez el e.mail.

def create_term(info_organism):
    '''
    Returns a term to make a search with Entrez module
        Parameters
            info_organism (str): string with organism 
        Returns:
            term(str or list of strs): search term
    '''
    # Verificar la presencia del caracter de separacion de organismos
        # Crear una lista a partir del caracter de separacion
        # recorrer la lista e ir generando terminos de busqueda
        
            # Agregar el nombre del organismo al termino
            
            # Agregar los nombres de los genes al termino
            # agrupandolos y usando el OL OR
    # Si solo se indico un organismo, no generar la lista
    
        # Agregar el nombre del organismo al termino
    
        # Agregar los nombres de los genes al termino
        # agrupandolos y usando el OL OR
        
        # Regresar el termino generado
        


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
    
        # Buscar db que contienen informacion
        # del termino dado
        

        # Recorrer la lista de resultados
        # de egquery y extraer las bases de datos y los conteos
        # de los match encontrados
        # Excluir 0s y Errores
            # Por cada base de datos extraer la lista de IDs
            # Agregar la db y los ids a un diccionario
            