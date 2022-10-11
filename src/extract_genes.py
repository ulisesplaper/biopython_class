'''
 NAME
      extract_gen_info.py
VERSION
        5t1.0.0
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Extrae informacion de un archivo en formato gene bank, 
        asi como de genes especificados por el usuario
        (nombre, producto, primeras bases del gen, transcrito y proteina)
USAGE
        extract_genes.py [-h] -gl GENE_LIST [GENE_LIST ...] [-o output file] path
         
ARGUMENTS
     -h, --help            show this help message and exit
    -gl GENE_LIST [GENE_LIST ...], --gene_list GENE_LIST [GENE_LIST ...]
                        list of genes to extract info
    -o output file, --output output file
                        path of the output file
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/extract_genes.py

'''
# Importar los modulos necesarios
from Bio import SeqIO
import argparse

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(
    description="Extract info of a genbank file and given genes")
parser.add_argument('file', metavar='path', help='path of the genbank file')
parser.add_argument('-gl', '--gene_list', nargs='+',
                    help='list of genes to extract info', required=True)
parser.add_argument('-o', '--output', metavar='output file',
                    help='path of the output file')


# Leer los argumentos desde la terminal
args = parser.parse_args()
file = args.file
gene_list = args.gene_list
output_file = args.output

try:
    # Definir la funcion resumen
    def resumen(file, genes):
        # parsear el archivo con SeqIO
        for gb_record in SeqIO.parse(file, "genbank"):
            continue
        # Agregar a un archivo output informacion general del archivo
        # obtenida de annotations o la feature 'source'
        with open(output_file, 'w') as output:
            output.write(
                'organism:' + gb_record.annotations['organism'] + "\n")
            output.write('date:' + gb_record.annotations['date'] + "\n")
            output.write(
                'country: ' + str(gb_record.features[0].qualifiers['country']) + "\n")
            if ('isolate' in gb_record.features[0].qualifiers):
                output.write('isolate number' +
                             gb_record.features[0].qualifiers['isolated'] + "\n")

        # Recorrer el las features del archivo y cuando se encuentre
        # con el gen buscado, aniadir la informacion del gen y agregarla
        # al archivo output
        for gene in gene_list:
            for genetic_element in gb_record.features:
                if (genetic_element.type == 'source'):
                    continue
                if (genetic_element.type == 'CDS'):
                    if (genetic_element.qualifiers['gene'] == [gene]):
                        dna_seq = gb_record.seq[genetic_element.location.
                                                nofuzzy_start:genetic_element.
                                                location.nofuzzy_start + 15]
                        with open(output_file, 'a') as output:
                            output.write('gene: ' + gene + '\n')
                            output.write(
                                'product: ' + str(genetic_element.qualifiers['product'][0]) + '\n')
                            output.write('ADN:' + str(dna_seq) + '\n')
                            output.write(
                                'ARN:' + str(dna_seq.transcribe()) + '\n')
                            output.write(
                                'Prot:' + str(dna_seq.translate()) + '\n')

    # Llamar a la funcion con los argumentos correspondientes
    resumen(file, gene_list)
except FileNotFoundError:
    print('\n FileNotFoundError: El nombre o ruta del archivo son\
 incorrectos\n')
