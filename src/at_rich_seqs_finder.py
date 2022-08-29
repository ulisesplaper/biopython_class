'''
 NAME
      at_rich_seqs_finder.py
VERSION
        1.0.1
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
       Encuentra regiones ricas en AT de tamanio especificado en un archivo fasta
USAGE
        usage: at_rich_seqs_finder.py [-h] -s size [--version] path
         
ARGUMENTS
        positional arguments:
         path                  path of the nucleotide sequence fasta file

        options:
         -h, --help            show this help message and exit
         -s size, --size size  size of minimal region rich in AT
         --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/at_rich_seqs_finder.py

'''
# Importar librerias necesarias para el programa
import re
import argparse
import fasta_tools as ft

# Definir error para cuando no se encuentre una secuencia rica en at.
class NotATRegionError(Exception):
        pass

# Determinar argumentos opcionales y posicionales del programa
parser = argparse.ArgumentParser(description="Determines the sequences rich in adenines/timines")
parser.add_argument('path', metavar='path',help='path of the nucleotide sequence fasta file')
parser.add_argument('-s','--size',metavar = 'size',help='size of minimal region rich in AT', required=True)
parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')

# leer argumentos desde la terminal y asignarlos a variables
try:
        args = parser.parse_args()
        sequence = args.path
        at_region_size = args.size

        # Extraer la secuencia de nts del archivo fasta
        sequence = ft.get_nucleotide_sequence(sequence)

        # Crear funcion para validar la secuencia del archivo fasta
        def validate_sequence(sequence):
                non_base_char = re.search(r"[^ATGC]", sequence)
                if non_base_char:
                        non_base_chars = re.finditer(r"[^ATGC]", sequence)
                        return non_base_chars
                else:
                        return None

        # Validar la secuencia con la funcion creada
        # Si la secuencia tiene bases ambiguas, informar su posicion y
        # su identidad
        ambiguos_bases = validate_sequence(sequence)
        if ambiguos_bases:
                for base in ambiguos_bases:
                        ambiguos_base = base.group()
                        pos = base.start()
                        print(f"\nSe encontro la base ambigua {ambiguos_base} en la posicion {pos+1}\n")


        # Si no contiene bases ambiguas, encontrar las secuencias que 
        # contengan AT en la minima cantidad indicada por el usuario y 
        # asignarlas a un iterable
        else:
                # Previamente determinar si hay match para secuencias ricas
                # en at del tamanio indicado por el usuario.
                # En caso de que no, levantar un error.
                if (not(re.search("(A|T){"+at_region_size+",}",sequence))):
                        raise NotATRegionError

                sequences = re.finditer("(A|T){"+at_region_size+",}",sequence)
                sequences = re.finditer("(A|T){"+at_region_size+",}",sequence)

        ## Recorrer el iterable y determinar el contenido de cada 
        ## uno de sus elementos (regiones ricas en AT) y su posicion en 
        # la secuencia original.
                for sequence in sequences:
                        at_seq = sequence.group()
                        position = sequence.start()
                        print(f"\nSe encontro la secuencia {at_seq} en la posicion {position+1}\n")
except FileNotFoundError:
        print('\n FileNotFoundError: El nombre o ruta del archivo son incorrectos\n')
except NotATRegionError:
        print("\nNo se encontraron regiones ricas en AT del tamanio especificado\n")

