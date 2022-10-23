import shutil
from constants import EMBOSS_DIR
from file_helper import delete_file, run_bash_file_with_arguments

EX5_SCRIPT = EMBOSS_DIR + "emboss_script.sh"

def run_exercise_5(input_path, output_file):
    file_name = copy_input_to_emboss_dir(input_path)
    run_bash_file_with_arguments(EX5_SCRIPT, [file_name, output_file, "400"])
    delete_file(EMBOSS_DIR + file_name)

def copy_input_to_emboss_dir(input_path):
    file_path_array = input_path.split('/')
    file_name = file_path_array[len(file_path_array)-1]
    shutil.copyfile(input_path, EMBOSS_DIR + file_name)
    return file_name