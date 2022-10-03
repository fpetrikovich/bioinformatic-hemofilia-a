import argparse

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from file_helper import run_bash_file, save_file, create_bash_file

E_VALUE_THRESHOLD = 0.01

def run_exercise_2(fasta_file, online_file_report, local_file_report):
	# Read file
	fasta_string = open(fasta_file).read()

	# Sanity string by removing empty lines 
	fasta_string = '\n'.join([i for i in fasta_string.split('\n') if len(i) > 0])

	sequences = separate_sequences(fasta_string)

	index = 1

	for sequence in sequences:
		print("----------\nStarting BLAST Query for ORF " + str(index))
		
		print("ONLINE QUERY IN PROGRESS...")
		blast_online_record = get_blast_online_record(sequence)

		online_report = analyze_blast_record(blast_online_record)

		save_file(online_file_report + "_ORF" + str(index) + ".report", online_report)
		
		print("OFFLINE QUERY IN PROGRESS...")
		run_blast_offline_query(local_file_report + "_ORF" + str(index) + ".report", sequence, str(index))

		index += 1

	return

def run_blast_offline_query(report_name, sequence, orf_index):
	# create file with the ORF sequence
	orf_file_name = "outputs/orfs/f8_protein_orf_" + orf_index + ".faa"
	save_file(orf_file_name, sequence)

	# name of the script to run local blast
	script_name = "outputs/bash/script_orf" + orf_index

	# Locate where local blastp command is
	command = "./ncbi-blast-2.13.0+/bin/blastp "
	# Specify the database
	command += "-db ncbi-blast-2.13.0+/data/swissprot "
	# file with the query
	command += "-query " + orf_file_name + " "
	# E value limit
	command += "-evalue " + str(E_VALUE_THRESHOLD) + " "
	# output file
	command += "-out " + report_name

	create_bash_file(script_name, command)
	run_bash_file("./" + script_name)


def separate_sequences(fasta_string): 
	fasta_array = [i for i in fasta_string.split('\n')]
	# Want to remove the >... line to only keep the sequence itself
	return [fasta_array[i+1] for i in range(0, len(fasta_array), 2)]

def get_blast_online_record(sequence):
	result_handle = NCBIWWW.qblast('blastp', 'nr', sequence)

	return NCBIXML.read(result_handle)

def analyze_blast_record(blast_record):
	output = ""

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESHOLD:
				output += "------Alignment------\n"
				output += "Sequence: %s\n" % alignment.hit_def.split(' >')[0]
				output += "Accession: %s\n" % alignment.hit_id.split('|')[1]
				output += "Length: %d\n" % alignment.length
				output += "Bits: %s\n" % str(hsp.bits)
				output += "Score: %s\n" % str(hsp.score)
				output += "Identiy: %d/%d(%.2f%%)\n" % (hsp.identities, hsp.align_length, (100 * hsp.identities / hsp.align_length))
				output += "E Value: %f\n" % hsp.expect
				output += "Query: %s\n" % hsp.query
				output += "Match: %s\n" % hsp.match
				output += "Subject: %s\n\n" % hsp.sbjct

	return output

