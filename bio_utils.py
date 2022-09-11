from Bio.Seq import Seq

def transcribe_dna(sequence):
    # Nucleotide sequences are normally read from the 5' to 3; direction
    coding_dna = sequence
    template_dna = coding_dna.reverse_complement()

    # Want to transcribe the coding strand into the corresponding mRNA
    messenger_rna = coding_dna.transcribe()
    
    # If you do want to do a true biological transcription starting with 
    # the template strand, then this becomes a two-step process:
    template_dna.reverse_complement().transcribe() 
    
    return messenger_rna

# Accepts mRNA sequence of coding DNA
def translate_sequence_to_protein(sequence, table = "Standard", stop_symbol = "*"):
    return sequence.translate(table=table, stop_symbol=stop_symbol)