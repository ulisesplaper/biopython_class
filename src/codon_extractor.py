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

# Definir argumentos opcionales y posicionales


# Leer los argumentos desde la terminal
# Ejercicio 1. Parte basica
# Convertir el archivo a un diccionario
# Recorre el diccionario y por cada id extraer la secuencia de DNA
# Determinar el primer codon de inicio (ya que puede haver varios)
# Recorrer la secuenca e ir imprimiendo los codones desde el codon de inico
# Terminar el recorrido si se alcanza un codon de paro

# Ejercicio 2. Parte avanzada
# Convertir el archivo a un diccionario
# Recorrer el diccionario por cada indice
# Determinar los codones de inicio, verificando que efectivamente
# haya al menos un codon

# Por cada ID, extraer los codones en tres marcos de lectura
# Por cada marco de lectura, extraer los codones e imprimirlos hasta
# alcanzar un codon de paro

# Para los marcos de lectura de la cadena reverse,
# obtener la cadena complementaria y proceder de la misma forma
# que para la secuencia anterior
