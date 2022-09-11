from Bio import SeqIO

def run_exercise_1(genbank_file, fasta_output_file):

    output_handle = open(fasta_output_file, "w")

    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        print('Working with GenBank sequence record %s' % seq_record.id)
        


    # orf = 1
    # for seq_record in SeqIO.parse(input_handle, format='genbank'):
    #     print('Dealing with GenBank record %s' % seq_record.id)
    #     for seq_feature in seq_record.features:
    #         # A CoDing Sequence (CDS) is a region of DNA or RNA whose sequence determines
    #         # the sequence of amino acids in a protein
    #         if seq_feature.type == 'CDS':
    #             output_handle.write('>lcl|ORF%s\n%s' % (orf, seq_feature.qualifiers['translation'][0]))
    #             orf += 1

    output_handle.write(seq_record.id)
    output_handle.close()
