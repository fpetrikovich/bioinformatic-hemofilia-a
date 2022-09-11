from Bio import SeqIO

from bio_utils import translate_sequence_to_protein, get_cds_feature_from_sequence_record

def run_exercise_1(genbank_file, fasta_nuc_output_file, fasta_prot_output_file):

    nucleotide_handle = open(fasta_nuc_output_file, "w")
    protein_handle = open(fasta_prot_output_file, "w")

    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        print('Working with GenBank sequence record %s' % seq_record.id)

        cds_feature = get_cds_feature_from_sequence_record(seq_record)
        
        gene_sequence = cds_feature.extract(seq_record.seq)
        protein_sequence = translate_sequence_to_protein(seq=gene_sequence, cds=True)

        # Halt if the translation does not match
        assert protein_sequence == cds_feature.qualifiers["translation"][0]

        nucleotide_handle.write('>lcl\n%s' % (gene_sequence))
        protein_handle.write('>lcl\n%s' % (protein_sequence))

 
    nucleotide_handle.close()
    protein_handle.close()
