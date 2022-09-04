''' NAME
        fasta_tools.py
VERSION
        1.1
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
    Este modulo contiene funciones que:
    1. Extraen la secuencia de aminoacidos o nucleotidos de un archivo fasta
    2. Extraen la secuencia de nucleotidos de un archivo fasta.
    
USAGE
    import fasta_tools
    
ARGUMENTS
    none
SEE ALSO
    none
GitHub link
    https://github.com/ulisesplaper/python_class/blob/master/src/fasta_tools.py
        
'''
# Importar librerias necesarias
import re

# Funcion get sequence


def get_sequence(path):
    '''
    Returns sequence extracted from a fasta file 
        Parameters
            path (str): path of the fasta file
        Returns:
            sequence(str): nucleotide or aminoacid sequence extracted from
            a fasta file
    '''
    # Abrir el archivo, guardar cada linea en una lista y cerrar el archivo
    file = open(path)
    fasta = file.readlines()
    file.close()
    # Recorrer cada elemento de la lista y, mientras no sea
    # el encabezado, aniadir cada elemento a un string
    # eliminar caracteres que no son secuencia
    # Devolver la secuencia
    sequence = ""
    for line in fasta:
        if re.search("^>", line):
            continue
        sequence += line
    sequence = sequence.replace("\n", "")
    return sequence
