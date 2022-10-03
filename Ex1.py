import math
from Bio import SeqIO
from file_helper import save_file
from bio_utils import translate_sequence_to_protein, check_if_valid_orf_with_cds, get_cds_protein_from_record

codon_size = 3
orf_count = 6

def run_exercise_1(genbank_file, fasta_nuc_output_file, fasta_prot_output_file, fasta_final_protein_output_file):

    nucleotide_handle = open(fasta_nuc_output_file, "w")
    proteins_handle = open(fasta_prot_output_file, "w")

    orf_status = [None] * (orf_count + 1)
    
    orf = 1

    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        print('Working with GenBank sequence record %s' % seq_record.id)

        sequence = seq_record.seq
        reverse_sequence = sequence.reverse_complement()

        for i in range(0, 3):
            orf_status = handle_new_open_reading_frame(orf, sequence, i, False, orf_status, nucleotide_handle, proteins_handle)
            orf += 1
            orf_status = handle_new_open_reading_frame(orf, reverse_sequence, i, True, orf_status, nucleotide_handle, proteins_handle)
            orf += 1

        for j in range(1, len(orf_status)):
            status = orf_status[j] 
            
            if check_if_valid_orf_with_cds(seq_record, status["protein"]):
                print("########################## START ##########################")
                print("Found correct Open Reading Frame based on CDS comparison:")
                print("ORF: %i" % j)
                print("Start: %i" % (status["start"]))
                print("End: %i" % (status["end"]))
                print("Is reverse complement sequence: ", (status["is_reverse"]))
                print("Protein sequence: %s" % (status["protein"]))
                print("########################### END ###########################")
                
                # Create valid protein fasta sequence
                save_file(fasta_final_protein_output_file, ">Seq1\n" + str(get_cds_protein_from_record(seq_record)))
                break

    nucleotide_handle.close()
    proteins_handle.close()


def handle_new_open_reading_frame(orf_id, nucleotides, start_offset, is_reverse, orf_status, nucleotide_handle, protein_handle):
    # Calculate where the reading frame starts and ends       
    start = start_offset
    end = math.floor((len(nucleotides) - start_offset)/codon_size) * codon_size + start_offset
    # translate sequence to protein for a given ORF
    protein_sequence = translate_sequence_to_protein(nucleotides[start:end], cds=False)
    # Save the translation to compare with the CDS
    orf_status[orf_id] = {
        "start": start,
        "end": end,
        "is_reverse": is_reverse,
        "protein": protein_sequence
    }
    # Write the ORF to the file
    protein_handle.write('>lcl|ORF%i\n%s\n' % (orf_id, protein_sequence))
    nucleotide_handle.write('>lcl|ORF%i\n%s\n' % (orf_id, nucleotides[start:end]))

    return orf_status