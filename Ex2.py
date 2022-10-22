import time

from Bio.Blast import NCBIWWW, NCBIXML
from file_helper import run_bash_file, save_file, create_bash_file

E_VALUE_THRESHOLD = 0.04

def run_exercise_2(fasta_file, file_report, is_local, online_blast_db):
	# Read file
	fasta_string = open(fasta_file).read()
	# Sanity string by removing empty lines 
	fasta_string = '\n'.join([i for i in fasta_string.split('\n') if len(i) > 0])
	# Remove all > lines and return an array of the sequences
	sequences = separate_sequences(fasta_string)
	# Run all blast consults and generate the reports
	run_blast_for_sequences(sequences, file_report, is_local, online_blast_db)


def run_blast_for_sequences(sequences, file_report, is_local, online_blast_db):
	"""
    Saves each sequence in a separate file and runs a local and remote blast 
	consult. Creates a report for each type of consult.
    Arguments:
        sequences: Array of aminoacid sequences
        online_file_report: main name the remote reports will have
		local_file_report: main name the local reports will have.
    """
	index = 1

	for sequence in sequences:
		# Create a separate file for the ORF to avoid multiple sequences in one file
		file_name = create_sequence_fasta_file(str(index), sequence)

		print("----------\nStarting " + ("Local" if is_local else "Online") + " BLAST Query for ORF " + str(index))

		if not is_local:			
			blast_online_record = get_blast_online_record(sequence, online_blast_db)

			online_report = analyze_blast_record(blast_online_record)

			save_file(file_report + "_ORF" + str(index) + ".report", online_report)
		
		else:			
			run_blast_offline_query(file_report + "_ORF" + str(index) + ".report", file_name, str(index))

		index += 1

	print("----------\n")


def create_sequence_fasta_file(orf_index, sequence):
	"""
    Creates a new fasta file for the ORF sequence.
    Arguments:
        orf_index: integer representing the ORF
        sequence: sequence to add to the fasta file
	Returns: Name of the file created
    """
	orf_file_name = "outputs/orfs/f8_protein_orf_" + orf_index + ".faa"
	save_file(orf_file_name, sequence)
	return orf_file_name


def run_blast_offline_query(report_name, file_name, orf_index):
	"""
    Creates a bash script with the command to run the local Blast 
	consult and executes the script.
    Arguments:
        report_name: name of the file where the consult response will be stored
        file_name: fasta file to use as blast query
		orf_index: index of the ORF being queried
    """

	# name of the script to run local blast
	script_name = "outputs/bash/script_orf" + orf_index

	# Locate where local blastp command is
	command = "./ncbi-blast-2.13.0+/bin/blastp "
	# Specify the database
	command += "-db ncbi-blast-2.13.0+/data/swissprot "
	# file with the query
	command += "-query " + file_name + " "
	# E value limit
	command += "-evalue " + str(E_VALUE_THRESHOLD) + " "
	# output file
	command += "-out " + report_name

	create_bash_file(script_name, command)
	run_bash_file("./" + script_name)


def separate_sequences(fasta_string): 
	"""
    Retrieves only the sequences from the contents of a fasta file.
    Arguments:
        fasta_string: contents of the fasta file
	Returns: array of sequences
    """
	fasta_array = [i for i in fasta_string.split('\n')]
	# Want to remove the >... line to only keep the sequence itself
	return [fasta_array[i+1] for i in range(0, len(fasta_array), 2)]


def get_blast_online_record(sequence, online_blast_db):
	"""
    Uses the swissprot database to query an aminoacid sequence.
    Arguments:
        sequence: sequence to query through a remote consult
	Returns: results of query
    """
	result_handle = NCBIWWW.qblast('blastp', online_blast_db, sequence)
	return NCBIXML.read(result_handle)


def analyze_blast_record(blast_record):
	"""
    Parses the blast consult response to facilitate its reading.
    Arguments:
        blast_record: response of the online blast consult.
	Returns: formatted string ready to be read
    """
	output = ""
	records_found = 0

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESHOLD:
				records_found += 1
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

	print("Found %i records.\n" % (records_found))

	return output

