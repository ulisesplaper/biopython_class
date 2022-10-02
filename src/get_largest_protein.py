'''
 NAME
      get_largest_protein.py
VERSION
        t3.1.0.0.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Determina la secuencia proteica de mayor longitud en
        cualquiera de los posibles ORFs de una secuencia dada
USAGE
        usage: py get_largest_protein.py [-h] [--version] path
         
ARGUMENTS
        positional arguments:
         path                  path of the mRNA sequence fasta file

        options:
        optional arguments:
        -h, --help            show this help message and exit
        --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/get_largest_protein.py

'''  # Importar los modulos necesarios
import argparse
from Bio.Seq import Seq, MutableSeq
from Bio import SeqUtils
from Bio import SeqIO


# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Determine the largest sequence\
protein given a DNA sequence")
parser.add_argument('path', metavar='path', help='path of the DNA\
sequence fasta file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
DNA_file = args.path

try:
    # Ejercicio 1. Parte Avanzada
    # Funcion get_largest_protein
    def get_largest_protein(dna_seq):
        '''
        Return the largest sequence\
        protein given a DNA sequence
                Parameters:
                        dna_seq(Seq_obj): DNA sequence
                Returns:
                        largest_prot(str): largest protein sequence
        '''
        # Determinar las posiciones de los codones de inicio
        # de las cadenas complementarias
        codon_posic = SeqUtils.nt_search(str(dna_seq), 'ATG')
        codon_posic_com = SeqUtils.nt_search(str(dna_seq.complement), 'ATG')

        # Recorrer la lista con los codenes de inicio
        # Por cada codon de inicio encontrado, traducir
        # la secuencia en sus 3 diferentes orfs
        # Agregar las secuencias proteicas a una lista
        prot_list = []
        for i in range(1, len(codon_posic)):
            for j in range(3):
                protein = dna_seq[(codon_posic[i]) +
                                  j:].translate(to_stop=True)
                prot_list.append(protein)
        # Proceder con la cadena complementaria de igual forma
        for i in range(1, len(codon_posic_com)):
            for j in range(3):
                protein = dna_seq[(codon_posic_com[i]) +
                                  j:].translate(to_stop=True)
                prot_list.append(protein)
        # Determinar la secuencia proteica mas larga
        # y regresarla
        largest_prot = max(prot_list, key=len)
        return(largest_prot)
    # Parte 3.
    # Crear diccionario con datos de cada secuencia
    id_dict = SeqIO.to_dict(SeqIO.parse(DNA_file, "fasta"))
    id_list = list(id_dict.keys())
    # Recorrer cada elemento del diccionario, extraer la secuencia
    # de DNA y determinar la secuencia proteica mas grande
    for id in id_list:
        biggest_prot = get_largest_protein(id_dict[id].seq)
        print(
            f"The biggest protein that can be obtained from the sequence {id_dict[id].id} is " + get_largest_protein(id_dict[id].seq))

except FileNotFoundError:
    print('\n FileNotFoundError: El nombre o ruta del archivo son\
 incorrectos\n')
