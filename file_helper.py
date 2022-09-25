from constants import mRNA_GENBANK_FILES


def file_param_to_file_name(param):
    file_names = {
        "1A": mRNA_GENBANK_FILES.ISOFORM_A_PREPROTEIN.value,
        "1B": mRNA_GENBANK_FILES.ISOFORM_B.value,
        "2": mRNA_GENBANK_FILES.FASTA_FILE.value
    }

    folder_names = {
        "1A": "inputs/",
        "1B": "inputs/",
        "2": "outputs/",
    }
 
    return folder_names.get(param, "") + file_names.get(param, "")

def generate_output_path(filename):
    return "outputs/" + filename
