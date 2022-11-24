import shutil
from constants import EMBOSS_DIR
from file_helper import delete_file, run_bash_file_with_arguments

EX5_SCRIPT = EMBOSS_DIR + "emboss_script.sh"

def run_exercise_5(input_path, orf_file, prosite_file, orf_minsize):
    """
    Escribir un script que llame a algún programa EMBOSS para que realice un análisis 
    sobre una secuencia de nucleótidos fasta (del Ej. 1). Por ejemplo, que calcule los
    ORFs y obtenga las secuencias de proteínas posibles. Luego bájense la bases de datos
    PROSITE (archivo prosite.dat) de dominios/motivos funcionales conocidos, por medio 
    del llamado a otro programa EMBOSS realizar el análisis de dominios de las secuencias 
    de aminoácidos obtenidas y escribir los resultados en un archivo de salida.
    """
    # Copy the file to the emboss dir to have it in the same folder as the script
    try:
        file_name = copy_input_to_emboss_dir(input_path)
    except Exception as e:
        print("[ERROR] Error while copying file to the emboss directory.\n[MESSAGE] " + str(e))
        exit(1)
        
    # Run the emboss script
    run_bash_file_with_arguments(EX5_SCRIPT, [EMBOSS_DIR, file_name, orf_file, prosite_file, str(orf_minsize)])
    # Delete the file we copied to the emboss dir to avoid duplicate files
    delete_file(EMBOSS_DIR + file_name)


def copy_input_to_emboss_dir(input_path):
    """
    Given the path to a file, will copy the file to the emboss directory.
    """
    file_path_array = input_path.split('/')
    file_name = file_path_array[len(file_path_array)-1]
    shutil.copyfile(input_path, EMBOSS_DIR + file_name)
    return file_name