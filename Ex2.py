import argparse

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

E_VALUE_THRESHOLD = 0.04

def run_exercise_2(fasta_file, online_file_report, local_file_report):
	# Read file
	fasta_string = open(fasta_file).read()

	# Sanity string by removing empty lines 
	fasta_string = '\n'.join([i for i in fasta_string.split('\n') if len(i) > 0])

	sequences = separate_sequences(fasta_string)

	index = 1

	for sequence in sequences:
		blast_online_record = get_blast_online_record(sequence)
		# blast_offline_record = get_blast_offline_record(sequence)

		online_report = analyze_blast_record(blast_online_record)
		# offline_report = analyze_blast_record(blast_offline_record)

		save_file(online_file_report + "_ORF" + str(index), online_report)
		# save_file(local_file_report + "_ORF" + str(index), offline_report)

		index += 1

	return

def separate_sequences(fasta_string): 
	fasta_array = [i for i in fasta_string.split('\n')]
	return [fasta_array[i] + '\n' + fasta_array[i+1] for i in range(0, len(fasta_array), 2)]

def get_blast_online_record(sequence):
	result_handle = NCBIWWW.qblast('blastp', 'nr', sequence)

	return NCBIXML.read(result_handle)

def get_blast_offline_record(sequence):
	blastx_cline = NcbiblastxCommandline(query=sequence, db="nr", evalue=0.001, outfmt=5, out="tmp.xml")
	blastx_cline()

	result_handle = open("tmp.xml")
	return NCBIXML.read(result_handle)

def analyze_blast_record(blast_record):
	output = ""

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESHOLD:
				output += "------Alignment------\n"
				output += "Sequence: %s\n" % alignment.hit_def.split(' >')[0]
				output += "Accession: %s\n" % alignment.hit_id.split('|')[1]
				output += "Length: %d" % alignment.length
				output += "Score: %s\n" % str(hsp.score)
				output += "Identiy: %d/%d(%.2f%%)\n" % (hsp.identities, hsp.align_length, (100 * hsp.identities / hsp.align_length))
				output += "E Value: %f\n" % hsp.expect
				output += "Query: %s\n" % hsp.query
				output += "Match: %s\n" % hsp.match
				output += "Subject: %s\n\n" % hsp.sbjct

	return output

def save_file(file_name, content): 
	f = open(file_name, "w")
	f.write(content)
	f.close()
