'''
 NAME
      extract_seqs.py
VERSION
        4t1.0.0
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Extrae la secuencias de un archivo fastq 
        cuyas bases esten por arriba de un score
        dado
USAGE
        usage: py extract_seqs.py [-h] -th threshold -o output [--version] path
         
ARGUMENTS
        positional arguments:
        path                  path of the fastq file

        optional arguments:
        -h, --help            show this help message and exit
        -th threshold, --threshold threshold
                                Minimal threshold for each BP
        -o output, --output output
                                name of the output file
        --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/get_largest_protein.py

'''
# Importar los modulos necesarios
from Bio import SeqIO
import numpy as np
import argparse

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Extract sequences whose BP are above \
a given threshold")
parser.add_argument('file', metavar='path', help='path of the fastq file')
parser.add_argument('-th', '--threshold', metavar='threshold', help='Minimal threshold for\
 each BP', required=True)
parser.add_argument('-o', '--output', metavar='output',
                    help='name of the output file', required=True)
parser.add_argument('--version', action='version', version='%(prog)s 4t.1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
file = args.file
umbral = int(args.threshold)
output = args.output
try:

    # Funcion get_seq_umbral
    def get_seq_umbral(file, umbral, output):
        ''' Print the sequence number whose
                BP are above a given threshold
                and the name of the output file
                with such sequences
                    Parameters:
                            file(string): fastq file path
                            umbral(int): threshold of BP
                    Returns:
                            None
        '''
        # Inicializar un contador
        # Abrir el archivo y recorrer cada uno de los SeqObj
        # Extraer los scores de cada SeqObj y determinar el mino
        # Si el minimo esta por arriba del umbral, escribir la
        # secuencia en un archivo e incrementar el contador
        n = 0
        with open(output, "w") as out_handle:
            for record in SeqIO.parse(file, 'fastq'):
                score = record.letter_annotations["phred_quality"]
                small_score = min(score)
                if umbral < small_score:
                    n += 1
                    out_handle.write(record.id + '\n' + str(record.seq) + '\n')
        # Imprimir el numero de secuencias que pasan el umbral
        # y un mensaje para el usuario
        print(f"Secuencias con puntajes mayor al umbral: {n}")
        print('Archivo de salida generado exitosamente')

    # Llamar a la funcion
    get_seq_umbral(file, umbral, output)

except FileNotFoundError:
    print('\n FileNotFoundError: El nombre o ruta del archivo son incorrectos\n')
