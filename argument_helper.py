import configparser

from file_helper import file_exists
from constants import ARGUMENTS, CONFIG_DICT, EX1_GB, EX2_DB, EX2_QUERY, EX2_LOCAL, EX2_REPORT, EX3_SEQS, EX3_OUT, EX4_BLAST, EX4_PATTERN, EX5_OUT, EX5_SEQ, EX5_SIZE

def create_config(config_file, exercise):
    if not file_exists(config_file):
        print("[ERROR] File %s does not exists in the root directory." % (config_file))
        exit(1)
    config = configparser.ConfigParser()
    config.read(config_file)
    exercise_config = config[CONFIG_DICT[exercise]]
    print(exercise_config)
    return exercise_config

def handle_arguments(exercise, args):
    use_config = bool(args.use_config)
    arguments = ARGUMENTS
    if (use_config):
        config_file = str(args.config_file)
        # Get the arguments from the config file
        config = create_config(config_file, exercise)
        for key in config: 
            value = config.get(key)
            # Convert the value too boolean if necessary
            if value == "yes" or value == "no": value = config.getboolean(key) 
            arguments[key] = value
    else:
        # Get the arguments from the command line
        arguments[EX1_GB] = args.genbank
        arguments[EX2_DB] = args.database
        arguments[EX2_LOCAL] = bool(args.local)
        arguments[EX2_QUERY] = args.query
        arguments[EX2_REPORT] = args.report
        arguments[EX3_SEQS] = args.sequences
        arguments[EX3_OUT] = args.output
        arguments[EX4_BLAST] = args.blast
        arguments[EX4_PATTERN] = args.pattern
        arguments[EX5_SEQ] = args.sequence
        arguments[EX5_OUT] = args.outputseq
        arguments[EX5_SIZE] = int(args.minsize)
    
    return arguments

def print_curr_arguments(exercise, arguments):
    print("[INFO] Using the following arguments:")

    if exercise == 1:
        print('\t- %s : %s' % (EX1_GB, arguments[EX1_GB]))
    if exercise == 2:
        print('\t - %s : %s' % (EX2_DB, arguments[EX2_DB]))
        print('\t - %s : %s' % (EX2_QUERY, arguments[EX2_QUERY]))
        print('\t - %s : %s' % (EX2_REPORT, arguments[EX2_REPORT]))
        print('\t - %s : %s' % (EX2_LOCAL, arguments[EX2_LOCAL]))
    if exercise == 3:
        print('\t - %s : %s' % (EX3_SEQS, arguments[EX3_SEQS]))
        print('\t - %s : %s' % (EX3_OUT, arguments[EX3_OUT]))
    if exercise == 4:
        print('\t - %s : %s' % (EX4_BLAST, arguments[EX4_BLAST]))
        print('\t - %s : %s' % (EX4_PATTERN, arguments[EX4_PATTERN]))
    if exercise == 5:
        print('\t - %s : %s' % (EX5_SEQ, arguments[EX5_SEQ]))
        print('\t - %s : %s' % (EX5_OUT, arguments[EX5_OUT]))
        print('\t - %s : %s' % (EX5_SIZE, arguments[EX5_SIZE]))    
    
    print() # new line
        
    