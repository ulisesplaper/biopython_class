'''
 NAME
      mrna_translator.py
VERSION
        v1.0.1.
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
       Traduce una secuencia de mRNA en una secuencia proteica
       e imprime la secuencia en pantalla o la guarda en un
       archivo fasta
USAGE
        usage: py mrna_translator.py [-h] [--version] path
         
ARGUMENTS
        positional arguments:
         path                  path of the mRNA sequence fasta file

        options:
        optional arguments:
        -h, --help            show this help message and exit
        -he heading           Heading of the output fasta file
        -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path of the output fasta file with the protein
                        sequence
  --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/mrna_translator.py

'''
# Importar los modulos necesarios
import argparse
import fasta_tools as ft

# Definir error para cuando se indique un output path,
#  pero no encabezado.
class NotHeadingError(Exception):
        pass

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Translate a sequences\
of mRNA into a protein sequence")
parser.add_argument('path', metavar='path',help='path of the mRNA\
sequence fasta file')
parser.add_argument('-he', metavar='heading', default=0, help="Heading\
 of the output fasta file")
parser.add_argument('-o','--output_path', default=0, help="Path of the\
 output fasta file with the protein sequence")
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
mRNA_file = args.path
output_path = args.output_path
heading = args.he
try:
        # Abrir el archivo que contiene la secuencia de mRNA,
        #  asignarlo a una variable y cerrarlo
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

        # Convertir la variable tipo string que contiene la seq.
        #  de mRNA en una lista de codones
        def mRNA_string_converter(mRNA_chain):
                length =len(mRNA_chain)
                mRNA_list=[]
                j = 0
                for i in range(3,length+3,3):
                        mRNA_list.append(mRNA_chain[j:i])
                        j = i
                return(mRNA_list)
        mRNA_list = mRNA_string_converter(mRNA_seq)

        # Recorre la lista de codones y, con base en el diccionario,
        #  llenar una lista nueva con la secuencia proteica.
        protein = []
        for codon in mRNA_list:
                protein.append(gencode.get(codon))

        # Convertir la lista de aminoacidos a una string de aminoacidos
        Strprotein = "".join(protein)

        # Comprobar si el usuario desea generar un archivo fasta 
        # generarlo con el modulo create_fasta de ser asi
        if output_path:
                if not(heading):
                        raise NotHeadingError
                ft.create_fasta(f'{output_path}', Strprotein, f'{heading}' )

        # Imprimir la secuencia proteica en pantalla si el usuario
        # no indica si desea generar un archivo fasta
        else:
                print(f"Secuencia porteica:\n{Strprotein}")

except FileNotFoundError:
        print('\n FileNotFoundError: El nombre o ruta del archivo son\
 incorrectos\n')
except NotHeadingError:
        print("\nNotHeadingError: Indique un encabezado para el archivo fasta\
con la secuencia proteica\n")
