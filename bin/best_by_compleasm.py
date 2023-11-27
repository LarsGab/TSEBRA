#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
from inspect import currentframe, getframeinfo
import shutil
from bbb_functions import parse_busco, find_input_files, find_genemark_gtf, run_busco, run_getanno

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023. All rights reserved."
__credits__ = ""
__license__ = "Artistic License 2.0"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

argparser = argparse.ArgumentParser(description = 'Find or build the best gene set generated with BRAKER and ' + 
                                    'TSEBRA minimizing missing BUSCOs.')
argparser.add_argument('-m', '--tmp_dir', type=str, required = True,
                          help = 'Temporary directory where intermediate files will be written')
argparser.add_argument('-d', '--braker_working_dir', type=str, required = True,
                            help = 'Output directory of BRAKER run')
argparser.add_argument('-y', '--tsebra', type=str, required = False,
                            help = 'Location of tesbra.py on system')
argparser.add_argument('-f', '--getanno', type=str, required = False,
                            help = 'Location of getAnnoFastaFromJoingenes.py on system')
argparser.add_argument('-g', '--genome', type=str, required = True,
                            help = 'Genome FASTA file')
argparser.add_argument('-t', '--threads', type=str, required = False, default = 1,
                       help = 'Number of threads to use for running BUSCO')
argparser.add_argument('-p', '--busco_db', type=str, required = True,
                            help = 'BUSCO lineage for running BUSCO')

args = argparser.parse_args()

def parse_busco(file):
    """
    Parse the BUSCO statistics from the specified file.

    Args:
        file (str): Path to the file containing BUSCO statistics.

    Returns:
        float: percentage of missing BUSCOs.

    Raises:
        SystemExit: If the file cannot be opened.

    """
    missing = 0
    try:
        with open(file, "r") as f:
            for line in f:
                stat_pattern = r'\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)\%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)'
                if re.search(stat_pattern, line):
                    print("BUSCO statistics found in file: " + file)
                    print(line)
                    missing = float(re.search(stat_pattern, line).group(5))
    except IOError:
        print("ERROR: Could not open file: " + file)
        sys.exit(1)
    return missing


def find_input_files(args):
    """
    Find all the required input files for the process.

    Args:
        args (argparse.Namespace): Command-line arguments.

    Returns:>=
        dict: Dictionary containing the file paths.

    """
    file_paths = {}
    file_paths["braker_aa"] = check_file(os.path.join(args.braker_working_dir, "braker.aa"))
    file_paths["augustus_aa"] = check_file(os.path.join(args.braker_working_dir, "Augustus", "augustus.hints.aa"))
    file_paths["genome"] = check_file(args.genome)  # required to generate genemark protein file
    file_paths["hints"] = check_file(os.path.join(args.braker_working_dir, "hintsfile.gff"))  # required to possibly run TSEBRA
    file_paths["braker_gtf"] = check_file(os.path.join(args.braker_working_dir, "braker.gtf"))
    file_paths["augustus_gtf"] = check_file(os.path.join(args.braker_working_dir, "Augustus", "augustus.hints.gtf"))

    return file_paths

def find_genemark_gtf(braker_working_dir):
    """
    Find the genemark.gtf file and corresponding training.gtf file in the specified directory.

    Args:
        braker_working_dir (str): Path to the BRAKER working directory.

    Returns:
        tuple: Tuple containing the paths to the genemark.gtf and training.gtf files.

    Raises:
        SystemExit: If no genemark.gtf file is found in the directory.

    """
    genemark_directories = ["GeneMark-ETP", "GeneMark-EP", "GeneMark-ET", "GeneMark-ES"]

    for genemark_dir in genemark_directories:
        genemark_gtf_file = os.path.join(braker_working_dir, genemark_dir, "genemark.gtf")
        training_gtf_file = os.path.join(braker_working_dir, genemark_dir, "training.gtf")
        if os.path.isfile(genemark_gtf_file):
            return check_file(genemark_gtf_file), check_file(training_gtf_file)

    error_msg = f"ERROR: No genemark.gtf file found in {braker_working_dir}"
    print(error_msg)
    sys.exit(1)


def run_busco(protein_files, threads, busco_db, tmp_dir):
    """
    Run BUSCO on the specified protein files using a bash script.

    Args:
        protein_files (list): List of protein files to run BUSCO on.
        threads (int): Number of threads to use.
        busco_db (str): Path to the BUSCO database.
        tmp_dir (str): Temporary directory to store the output.

    Returns:
        dict: Dictionary containing the paths to the BUSCO result files.

    Raises:
        SystemExit: If there is an error in writing the bash script, executing it, or removing the script.

    """
    busco_script = os.path.join(tmp_dir, "busco.sh")
    try:
        with open(busco_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("source activate busco_env\n")
            for protein_file in protein_files:
                # identify the gene prediction program from the protein file name
                tool = re.search(r'^([^.]+)\.', os.path.basename(protein_file)).group(1)
                f.write(
                    f"busco -c {threads} -i {protein_file} -o {tmp_dir}/{tool}_busco "
                    f"-l {busco_db} -m protein\n"
                )
    except IOError:
        error_msg = f"ERROR: Could not write to file {busco_script}"
        print(error_msg)
        sys.exit(1)

    busco_cmd = ["bash", busco_script]
    run_simple_process(busco_cmd)

    result_dict = {}
    for protein_file in protein_files:
        # identify the gene prediction program from the protein file name
        tool = re.search(r'^([^.]+)\.', os.path.basename(protein_file)).group(1)
        result_dict[tool] = check_file(tmp_dir + "/" + tool + "_busco/short_summary.specific." + 
                                      busco_db + "." + tool + "_busco.txt")
    
    try:
        # remove the buso script
        os.remove(busco_script)
    except OSError:
        print("ERROR: Could not remove file: " + busco_script)
        sys.exit(1)
    
    return result_dict

def run_getanno(annobin, genome_file, gtf, output_dir):
    """
    Run the getAnnoFastaFromJoingenes.py tool to process the Genemark GTF file and generate protein sequences.

    Args:
        annobin: getAnnoFastaFromJoingenes.py script.
        genome_file (str): Path to the genome file.
        gtf (str): Path to the GTF file.
        output_dir (str): Output directory to store the generated protein sequences

    """
    tool = re.search(r'^([^.]+)\.', os.path.basename(gtf)).group(1)
    cmd = [annobin, '-g', genome_file, '-f', gtf, '-o', output_dir + "/" + tool]
    run_simple_process(cmd)
    return check_file(output_dir + "/" + tool + ".aa")

def run_simple_process(args_lst):
    """
    Execute a subprocess command with the provided arguments.

    Args:
        args_lst (list): List of command-line arguments.

    Returns:
        CompletedProcess: Result of the subprocess execution.

    Raises:
        SystemExit: If the subprocess returns a non-zero exit code.

    """
    try:
        print(str(datetime.now().date()) + " " +
          str(datetime.now().time()) + " Trying to execute the following command:")
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(str(datetime.now().date()) + " " +
          str(datetime.now().time()) + " Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " +
                  str(result.returncode) + str(result.args))
            sys.exit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexc.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        sys.exit(1)


def check_binary(binary, name):
    """
    Check if the provided binary is executable or available in the system's PATH.

    Args:
        binary (str): Path to the binary file or None.
        name (str): Name of the binary/tool.

    Returns:
        str: Path to the binary if found and executable.

    Raises:
        SystemExit: If the binary is not found, not executable, or not in the PATH.

    """        # check if binary is in PATH
    if binary is not None:
        if os.path.isfile(binary):
            # check if binary is executable
            if os.access(binary, os.X_OK):
                return binary
            else:
                print("ERROR: " + name + " binary is not executable: " + binary)
                sys.exit(1)
        else:
            print("ERROR: " + name + " binary not found at " + binary)
            sys.exit(1)
    else:
        # check if binary is in PATH
        if shutil.which(name) is not None:
            return shutil.which(name)
        else:
            print("ERROR: " + name + " binary not found in PATH")
            sys.exit(1)


def check_file(file):
    """
    Check if the specified file exists.

    Args:
        file (str): Path to the file.

    Returns:
        str: Path to the file if it exists.

    Raises:
        SystemExit: If the file is not found.

    """
    if os.path.isfile(file):
        return file
    else:
        print("ERROR: File not found: " + file)
        sys.exit(1)


def check_dir(dir):
    """
    Check if the specified directory exists.
# Step 1: Find all input files
file_paths = find_input_files(args)

# Step 2: find out which GeneMark was used: ETP or EP or ET or ES
file_paths["genemark_gtf"], file_paths["training_gtf"] = find_genemark_gtf(args.braker_working_dir)

# Step 3: Check if all dependencies are available
args.tsebra = check_binary(args.tsebra, "tsebra.py")
args.getanno = check_binary(args.getanno, "getAnnoFastaFromJoingenes.py")
# we do not check for BUSCO because it is hidden in an environment

# Step 4: Check if temporary directory exists
# if not, create it
if not os.path.exists(args.tmp_dir):
    try:
        os.makedirs(args.tmp_dir)
    except OSError:
        print("ERROR: Creation of the directory %s failed" % args.tmp_dir)
        sys.exit(1)

# Step 5: Create protein sequene file for GeneMark gtf file
file_paths["genemark_aa"] = run_getanno(args.getanno, args.genome, file_paths["genemark_gtf"], args.tmp_dir)

# Step 6: run BUSCO
protein_file_list = [file_paths["genemark_aa"], file_paths["braker_aa"], file_paths["augustus_aa"]]
busco_out_dict = run_busco(protein_file_list, args.threads, args.busco_db, args.tmp_dir)

# Step 7: parse and compare the number of missing BUSCOs
genemark_missing = parse_busco(busco_out_dict["genemark"])
augustus_missing = parse_busco(busco_out_dict["augustus"])
braker_missing = parse_busco(busco_out_dict["braker"])

# Step 8: Decide whether the provided BRAKER gene set is good enough
if braker_missing <= augustus_missing and braker_missing <= genemark_missing and braker_missing <= 5:
    print("The BRAKER gene set " + file_paths["braker_gtf"] + " is the best one. It lacks " + str(braker_missing) + " BUSCOs.")
    sys.exit(0)
elif augustus_missing >= genemark_missing:
    tsebra_force = file_paths["genemark_gtf"]
    not_tsebra_force = file_paths["augustus_gtf"]
else:
    tsebra_force = file_paths["augustus_gtf"]
    not_tsebra_force = file_paths["genemark_gtf"]

# Step 9: Run TSEBRA and enforce the best gene set
tsebra_cmd = [args.tsebra, "-k", file_paths["training_gtf"] + "," + tsebra_force,
              "-g", not_tsebra_force, "-e", file_paths["hints"], "-o", args.tmp_dir + "/better_braker.gtf"]
run_simple_process(tsebra_cmd)

# Step 10: generate protein sequence file for the new BRAKER gene set
better_braker_gtf = check_file(args.tmp_dir + "/better_braker.gtf")
bb_aa = run_getanno(args.getanno, args.genome, better_braker_gtf, args.tmp_dir)

# Step 11: Run BUSCO on the new BRAKER gene set
secondary_busco_out = run_busco([bb_aa], args.threads, args.busco_db, args.tmp_dir)

# Step 12: parse BUSCO output and report numbers
better_braker_missing = parse_busco(list(secondary_busco_out.values())[0])
if better_braker_missing < braker_missing:
    print("The new best BRAKER gene set is " + better_braker_gtf + ". BUSCO results:")
    parse_busco(better_braker_missing)
else:
    print("WARNING: The new BRAKER gene set is not better than the original one :-(")

    Returns:
        str: Path to the directory if it exists.

    Raises:
        SystemExit: If the directory is not found.

    """
    if os.path.isdir(dir):
        return dir
    else:
        print("ERROR: Directory not found: " + dir)
        sys.exit(1)

def main():
    """
    Execute workflow according to the following steps:

    Step 1: Find all input files
    Step 2: Determine the GeneMark variant used
    Step 3: Check if all dependencies are available
    Step 4: Create the temporary directory if it doesn't exist
    Step 5: Generate protein sequence file for the GeneMark GTF file
    Step 6: Run BUSCO on multiple protein files
    Step 7: Parse and compare the number of missing BUSCOs
    Step 8: Decide if the provided BRAKER gene set is sufficient. If not:
    Step 9: Run TSEBRA and enforce the best gene set
    Step 10: Generate protein sequence file for the new BRAKER gene set
    Step 11: Run BUSCO on the new BRAKER gene set
    Step 12: Parse the BUSCO output and report the numbers
    """
    
    # Step 1: Find all input files
    file_paths = find_input_files(args)

    # Step 2: find out which GeneMark was used: ETP or EP or ET or ES
    file_paths["genemark_gtf"], file_paths["training_gtf"] = find_genemark_gtf(args.braker_working_dir)

    # Step 3: Check if all dependencies are available
    args.tsebra = check_binary(args.tsebra, "tsebra.py")
    args.getanno = check_binary(args.getanno, "getAnnoFastaFromJoingenes.py")
    # we do not check for BUSCO because it is hidden in an environment

    # Step 4: Check if temporary directory exists
    # if not, create it
    if not os.path.exists(args.tmp_dir):
        try:
            os.makedirs(args.tmp_dir)
        except OSError:
            print("ERROR: Creation of the directory %s failed" % args.tmp_dir)
            sys.exit(1)

    # Step 5: Create protein sequene file for GeneMark gtf file
    file_paths["genemark_aa"] = run_getanno(args.getanno, args.genome, file_paths["genemark_gtf"], args.tmp_dir)

    # Step 6: run BUSCO
    protein_file_list = [file_paths["genemark_aa"], file_paths["braker_aa"], file_paths["augustus_aa"]]
    busco_out_dict = run_busco(protein_file_list, args.threads, args.busco_db, args.tmp_dir)

    # Step 7: parse and compare the number of missing BUSCOs
    genemark_missing = parse_busco(busco_out_dict["genemark"])
    augustus_missing = parse_busco(busco_out_dict["augustus"])
    braker_missing = parse_busco(busco_out_dict["braker"])
    
    # Step 8: Decide whether the provided BRAKER gene set is good enough
    print("BRAKER is missing " + str(braker_missing) + " BUSCOs.")
    print("GeneMark is missing " + str(genemark_missing) + " BUSCOs.")
    print("Augustus is missing " + str(augustus_missing) + " BUSCOs.")

    if braker_missing <= augustus_missing and braker_missing <= genemark_missing:
        print("The BRAKER gene set " + file_paths["braker_gtf"] + " is the best one. It lacks " + str(braker_missing) + "% BUSCOs.")
        sys.exit(0)
    elif augustus_missing >= genemark_missing:
        print("We end up in the second case.")
        tsebra_force = file_paths["genemark_gtf"]
        not_tsebra_force = file_paths["augustus_gtf"]
    else:
        print("We end up in the third case.")
        tsebra_force = file_paths["augustus_gtf"]
        not_tsebra_force = file_paths["genemark_gtf"]

    # Step 9: Run TSEBRA and enforce the best gene set
    tsebra_cmd = [args.tsebra, "-k", file_paths["training_gtf"] + "," + tsebra_force,
                  "-g", not_tsebra_force, "-e", file_paths["hints"], "-o", args.tmp_dir + "/better_braker.gtf"]
    run_simple_process(tsebra_cmd)

    # Step 10: generate protein sequence file for the new BRAKER gene set
    better_braker_gtf = check_file(args.tmp_dir + "/better_braker.gtf")
    bb_aa = run_getanno(args.getanno, args.genome, better_braker_gtf, args.tmp_dir)

    # Step 11: Run BUSCO on the new BRAKER gene set
    secondary_busco_out = run_busco([bb_aa], args.threads, args.busco_db, args.tmp_dir)

    # Step 12: parse BUSCO output and report numbers
    better_braker_missing = parse_busco(list(secondary_busco_out.values())[0])
    if better_braker_missing < braker_missing:
        print("The new best BRAKER gene set is " + better_braker_gtf + ".")
    else:
        # if the new gene set is not superior, produce an output that tells the user what
        # of the previously existing gene sets had the lowest percentage of missing BUSCOs
        print("WARNING: The new BRAKER gene set is not better than the original one :-(")
        print("The best gene set produced by the original BRAKER run is:")
        gene_sets = [file_paths["braker_gtf"], file_paths["genemark_gtf"], file_paths["augustus_gtf"]]
        # Create a dictionary to map gene set names to their respective missing values
        gene_set_values = {
            file_paths["braker_gtf"]: braker_missing,
            file_paths["genemark_gtf"]: genemark_missing,
            file_paths["augustus_gtf"]: augustus_missing
        }
        # Find the minimum value among the three
        min_value = min(gene_set_values.values())
        # Find the gene set name associated with the minimum value
        min_gene_set = None
        for gene_set, value in gene_set_values.items():
            if value == min_value:
                min_gene_set = gene_set
                break
        # Print the name of the gene set with the lowest number of missing BUSCOs
        print(min_gene_set)

if __name__ == "__main__":
    main()