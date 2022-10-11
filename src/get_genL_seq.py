'''
 NAME
      get_genL_seq.py
VERSION
        5t_2_1.0.0
AUTHOR
        Victor Ulises Plascencia Perez
DESCRIPTION
        Extrae secuencia de DNA, RNA y proteina
        del gen L del archivo data/virus.gb
        
USAGE
        usage: py get_genL_seq.py
         
ARGUMENTS
        None
SEE ALSO
GitHub link
        https://github.com/ulisesplaper/python_class/blob/master/src/get_genL_seq.py
'''
# Importar los modulos necesarios
from Bio import SeqIO
# Extraer el objeto seq record
for gb_record in SeqIO.parse('data/virus.gb', "genbank"):
    continue
# Extrae la secuencia del gen L
# cuando el qualifiers gene coincida con su nombre
# Traducir y transcribir la secuencia
for genetic_element in gb_record.features:
    if (genetic_element.type == 'source'):
        continue
    if (genetic_element.type == 'gene'):
        if (genetic_element.qualifiers['gene'] == ['L']):
            gen_L_cds = gb_record.seq[genetic_element.location.nofuzzy_start:
                                      genetic_element.location.nofuzzy_end]
# Guardae los resultados en un archivo
with open("results/genL_seq.txt", 'w') as file:
    file.write('DNA: ' + str(gen_L_cds) + '\n\n')
    file.write('protein: ' +
               str(gen_L_cds.translate()) + '\n\n')
    file.write('RNA: ' + str(gen_L_cds.transcribe()) + '\n\n')
print(f"\narchivo results/genL_seq.txt generado exitosamente\n")
