'''
 NAME
      genbank_inf_extract.py
VERSION
        4t_2.1.0.0
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Extrae informacion de archivos
        tipo gb
USAGE
        usage: py genbank_info_extract.py path
         
ARGUMENTS
        positional arguments:
        path                  path of the fastq file

        --version             show program's version number and exit
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/genbank_info_extract.py

'''
# Importar los modulos necesarios
from Bio import SeqIO
import argparse

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Extract info of a genbank file")
parser.add_argument('file', metavar='path', help='path of the genbank file')


# Leer los argumentos desde la terminal
args = parser.parse_args()
file = args.file

try:
    # Recorrer cada SeqObj en el archivo
    # Extraer la informacion (fecha, organismo y pais)
    # Imprimir
    for gb_record in SeqIO.parse(file, "genbank"):
        fecha = gb_record.annotations["date"]
        organismo = gb_record.annotations['organism']
        features = gb_record.features[0].qualifiers['country']
        print(
            f"information:\nFecha: {fecha}\n organismo:{organismo}\nPais:{features}")

except FileNotFoundError:
    print('\n FileNotFoundError: El nombre o ruta del archivo son incorrectos\n')
