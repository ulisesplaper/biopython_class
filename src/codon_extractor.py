'''
 NAME
      codon_extractor.py
VERSION
        v3_1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Extrae los codones de una secuencia de DNA iniciando con 
        el primer codon de inicio encontrado considerando 6 marcos
        de lectura
USAGE
        usage: py codon_extractor.py [-h] [--version] path
         
ARGUMENTS
        positional arguments:
         path                  path of the DNA sequence fasta file

        options:
        optional arguments:
        -h, --help            show this help message and exit
        --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/codon_extractor.py
'''
# Importar modulos necesarios
import argparse
from Bio.Seq import Seq, MutableSeq
from Bio import SeqUtils
from Bio import SeqIO

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Extract the codons from a \
DNA sequence starting with the first start codon found in six ORFs")
parser.add_argument('path', metavar='path', help='path of the DNA\
sequence fasta file')
parser.add_argument('--version', action='version',
                    version='%(prog)s v3_2.1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
path = args.path

try:
    # Ejercicio 1. Parte basica
    # Convertir el archivo a un diccionario
    id_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
    id_list = list(id_dict.keys())
    # Recorre el diccionario y por cada id extraer la secuencia de DNA
    # Determinar el primer codon de inicio (ya que puede haver varios)
    # Recorrer la secuenca e ir imprimiendo los codones desde el codon de inico
    # Terminar el recorrido si se alcanza un codon de paro
    for id in id_list:
        print(f'>{id_dict[id].id}')
        start_codon = SeqUtils.nt_search(str(id_dict[id].seq), 'ATG')[1]
        for j in range(start_codon, len(id_dict[id].seq), 3):
            if id_dict[id].seq[j: j + 3] == ('TGA') or id_dict[id].seq[j: j + 3] == 'TAA' or id_dict[id].seq[j: j + 3] == 'TAG':
                break
            print(id_dict[id].seq[j: j + 3], end=' ')
        print('\n')

    # Ejercicio 2. Parte avanzada

    # Convertir el archivo a un diccionario
    id_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
    id_list = list(id_dict.keys())
    # Recorrer el diccionario por cada indice
    # Determinar los codones de inicio, verificando que efectivamente
    # haya al menos un codon
    for id in id_list:
        start_codons = SeqUtils.nt_search(str(id_dict[id].seq), 'ATG')
        # Aniadimos una cadena vacia a start_codons como verificador
        start_codons.append('')
        if start_codons[1] == '':
            print(f"Not start codon found for sequence {id_dict[id].id}")
            continue
        # Por cada ID, extraer los codones en tres marcos de lectura
        # Por cada marco de lectura, extraer los codones e imprimirlos hasta
        # alcanzar un codon de paro
        for i in range(3):
            print(f'>{id_dict[id].id} ORF {i+1}')
            start_codon = start_codons[1] + i
            for j in range(start_codon, len(id_dict[id].seq), 3):
                if id_dict[id].seq[j: j + 3] == ('TGA') or id_dict[id].seq[j: j + 3] == 'TAA' or id_dict[id].seq[j: j + 3] == 'TAG':
                    break
                print(id_dict[id].seq[j: j + 3], end=' ')
            print('\n')

        # Para los marcos de lectura de la cadena reverse,
        # obtener la cadena complementaria y proceder de la misma forma
        # que para la secuencia anterior
        seq_complement = (id_dict[id].seq).complement()
        start_codons = SeqUtils.nt_search(str(seq_complement), 'ATG')
        # Aniadimos una cadena vacia a start_codons como verificador
        start_codons.append('')
        if start_codons[1] == '':
            print(
                f"Not start codon found for complementary sequence {id_dict[id].id}")
            continue
        for i in range(3):
            print(f'>{id_dict[id].id} ORF {i+4}')
            start_codon = start_codons[1] + i
            for j in range(start_codon, len(id_dict[id].seq), 3):
                if id_dict[id].seq[j: j + 3] == ('TGA') or id_dict[id].seq[j: j + 3] == 'TAA' or id_dict[id].seq[j: j + 3] == 'TAG':
                    break
                print(id_dict[id].seq[j: j + 3], end=' ')
            print('\n')
except FileNotFoundError:
    print('\n FileNotFoundError: El nombre o ruta del archivo son\
 incorrectos\n')
