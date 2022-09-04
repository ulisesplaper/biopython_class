'''
 NAME
      mrna_translator.py
VERSION
        0.1.0
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
       Traduce una secuencia de mRNA en una secuencia proteica
USAGE
        usage: py mrna_translator.py [-h] [--version] path
         
ARGUMENTS
        positional arguments:
         path                  path of the mRNA sequence fasta file

        options:
         -h, --help            show this help message and exit
         --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/mrna_translator.py

'''
# Importar los modulos necesarios
import argparse
import fasta_tools as ft

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Translate a sequences of mRNA into a protein sequence")
parser.add_argument('path', metavar='path',help='path of the mRNA sequence fasta file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
mRNA_file = args.path

# Abrir el archivo que contiene la secuencia de mRNA, asignarlo a una variable y cerrarlo
# (usar modulo fasta tools).
mRNA_seq = ft.get_sequence(mRNA_file)

# Crear un diccionario con el codigo genetico {codon:aminoacido}
# Se deja el codigo con Ts dado que las seqs de mRNA bajadas del NCBI
# usan Ts y no Us
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 
    'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 
    'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 
    'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 
    'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 
    'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 
    'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 
    'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 
    'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 
    'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

# Convertir la variable tipo string que contiene la seq. de mRNA en una lista de codones
def mRNA_string_converter(mRNA_chain):
    length =len(mRNA_chain)
    mRNA_list=[]
    j = 0
    for i in range(3,length+3,3):
        mRNA_list.append(mRNA_chain[j:i])
        j = i
    return(mRNA_list)
mRNA_list = mRNA_string_converter(mRNA_seq)

# Recorre la lista de codones y, con base en el diccionario, llenar una lista nueva
# con la secuencia proteica.
protein = []
for codon in mRNA_list:
    protein.append(gencode.get(codon))

# Imprimir la secuencia proteica en pantalla
Strprotein = "".join(protein)
print(Strprotein)