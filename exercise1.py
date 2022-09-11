from Bio import SeqIO

def run_exercise_1(genbank_file):
    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))