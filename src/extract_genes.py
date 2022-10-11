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


# Definir argumentos opcionales y posicionales

# Leer los argumentos desde la terminal

# Definir la funcion resumen

# parsear el archivo con SeqIO
# Agregar a un archivo output informacion general del archivo
# obtenida de annotations o la feature 'source'

# Recorrer el las features del archivo y cuando se encuentre
# con el gen buscado, aniadir la informacion del gen y agregarla
# al archivo output

# Llamar a la funcion con los argumentos correspondientes
