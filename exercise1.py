from Bio import SeqIO

def run_exercise_1(genbank_file):
    for seq_record in SeqIO.parse("Inputs/GenBank_F8_mRNA/NM_000132.4_sequence.gb", "genbank"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))