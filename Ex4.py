from Bio import SeqIO, Entrez
from constants import EX4_OUTPUT_DIR, ALIGNMENT_HEADER, EMAIL
from file_helper import file_exists
from error_helper import exit_with_error

Entrez.email = EMAIL  # Always tell NCBI who you are

def run_exercise_4(blast_file, pattern):
    """
    Escribir un script para analizar (parsear) un reporte de salida de blast que
    identifique los hits que en su descripción aparezca un Pattern determinado que 
    le damos como parámetro de entrada. El pattern puede ser una palabra. Punto extra: 
    pueden a su vez parsear cuál es el ACCESSION del hit identificado (donde hay una 
    coincidencia del Pattern) y con el módulo Bio::DB::GenBank obtener la secuencia 
    completa del hit en formato FASTA y escribirla a un archivo, es decir, levantar 
    las secuencias originales completas de los hits seleccionados.
    """
    if not file_exists(blast_file):
        exit_with_error("Blast report file with path %s does not exists." % blast_file)

    f = open(blast_file, "r")
    spaceless_pattern = pattern.replace(' ', '_')
    output = open(EX4_OUTPUT_DIR + spaceless_pattern + "_hits.report", "w")

    line = f.readline()
    while line:
        if line == ALIGNMENT_HEADER:
            # read next line (description)
            line = f.readline()
            # remove the "sequence: " portion of the line to get only the description
            descr = line.split(' ', 1)[1]
            # check if the pattern is present
            if pattern.lower() in descr.lower():
                # Pattern was found => save accession and hit
                handle_alignment_hit(output, f, line, spaceless_pattern)
    
        line = f.readline()
        
    # close the files
    f.close()
    output.close()


def handle_alignment_hit(output, input_file, line, pattern):
    """
    Will save the alignment hit to the output file and fetch the accession fasta file.
    """
    print("[INFO] Hit found. Processing...")

    output.write(ALIGNMENT_HEADER)
    output.write(line)
    
    # Get the accession and fetch the complete sequence of it
    line = input_file.readline()
    accession_id = line.split(' ')[1].rstrip() # Format = Accession: <accession_id>
    output.write(line)
    fetch_accession_sequence(accession_id, pattern)
    
    # While the end of the hit is not found, print the lines
    while line != '\n':
        line = input_file.readline()
        output.write(line)
    
    output.write("\n")


def fetch_accession_sequence(accession_id, pattern):
    """
    Given a accession id and pattern, will retrieve the fasta file with said id using
    Entrez and save the record to another fasta file with the pattern as part of the name.
    """
    handle = Entrez.efetch(db="protein", id=accession_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    SeqIO.write(record, EX4_OUTPUT_DIR + "/%s_%s.fasta" % (pattern, accession_id), "fasta")
    handle.close()