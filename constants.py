from enum import Enum

ALIGNMENT_HEADER = "------Alignment------\n"

## CONFIGURATION + ARGUMENTS CONSTANTS ##

CONFIG_FILE = "configuration.ini"
CONFIG_DICT = {
    1: "ex1",
    2: "ex2",
    3: "ex3",
    4: "ex4",
    5: "ex5"
}

EX1_GB = "genbank"
EX2_DB = "database"
EX2_QUERY = "query"
EX2_LOCAL = "local"
EX2_REPORT ="report"
EX3_SEQS = "sequences"
EX3_OUT = "output"
EX4_BLAST = "blast"
EX4_PATTERN = "pattern"
EX5_SEQ = "sequence"
EX5_OUT = "outputseq"

ARGUMENTS = {
    EX1_GB: None,
    EX2_DB: None,
    EX2_QUERY: None,
    EX2_LOCAL: None,
    EX2_REPORT: None,
    EX3_SEQS: None,
    EX3_OUT: None,
    EX4_BLAST: None,
    EX4_PATTERN: None,
    EX5_SEQ: None,
    EX5_OUT: None,
}

## DIRECTORY CONSTANTS ##

EX1_OUTPUT_DIR = "fasta_outputs/"
EX4_OUTPUT_DIR = "pattern_hits/"
BLAST_REPORTS_DIR = "reports/"
EMBOSS_DIR = "emboss/"
GENBANK_DIR = "genbank/"
MSA_OUTPUT_DIR = "msa/"

## FILE HANDLING CONSTANTS ##

FASTA_EXTENSION = ".fasta"
FASTA_TYPES = ['fas', 'fasta', 'faa']
OUTPUT_TYPES = ['txt', 'output', 'out']
ORFS_FILE_SUFFIX = "_proteins"
CORRECT_ORF_FILE_SUFFIX = "_final_protein"

# mRNA coagulation factor VIII GenBank Files
class mRNA_GENBANK_FILES(Enum):
    ISOFORM_A_PREPROTEIN = "NM_000132.4_sequence.gb"
    ISOFORM_B = "NM_019863.3_sequence.gb"

class FASTA_FILES(Enum):
    ORFS_ISOFORM_A_PREPROTEIN = "1A" + ORFS_FILE_SUFFIX + FASTA_EXTENSION
    ORFS_ISOFORM_B = "1B" + ORFS_FILE_SUFFIX + FASTA_EXTENSION
    ISOFORM_A_PREPROTEIN = "1A" + CORRECT_ORF_FILE_SUFFIX + FASTA_EXTENSION
    ISOFORM_B = "1B" + CORRECT_ORF_FILE_SUFFIX + FASTA_EXTENSION

FILE_IDENTIFIERS = {
    "1A": mRNA_GENBANK_FILES.ISOFORM_A_PREPROTEIN.value,
    "1B": mRNA_GENBANK_FILES.ISOFORM_B.value,
    "2A": FASTA_FILES.ORFS_ISOFORM_A_PREPROTEIN.value,
    "2B": FASTA_FILES.ORFS_ISOFORM_B.value,
    "2AProtein": FASTA_FILES.ISOFORM_A_PREPROTEIN.value,
    "2BProtein": FASTA_FILES.ISOFORM_B.value
}

FOLDER_IDENTIFIERS = {
    "1A": GENBANK_DIR,
    "1B": GENBANK_DIR,
    "2A": EX1_OUTPUT_DIR,
    "2B": EX1_OUTPUT_DIR,
    "2AProtein": EX1_OUTPUT_DIR,
    "2BProtein": EX1_OUTPUT_DIR,
}

## OTHER CONSTANTS ##

START_CODON = "AUG"
END_CODONS = ["UAA", "UAG", "UGA"]

EMAIL = "fpetrikovich@itba.edu.ar"

