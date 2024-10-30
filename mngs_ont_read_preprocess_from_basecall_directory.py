"""
You've run ONT sequencing, and you've basecalled the reads.
As a result you have a directory containing fastq_pass, fastq_fail, pod5, etc, ... as subdirectories.

Given that directory path as an input, the current script will do 

1) create the raw read fastq.gz file per barcode 
2) create the QC-passed fastq.gz file per barcode
3) create the host-depleted QC-passed fastq.gz file per barcode
"""


import os
import argparse
from datetime import datetime
import subprocess
import gzip
import tarfile


def get_arguments_from_command_line():
    # Setting
    parser = argparse.ArgumentParser(description = "Script for preprocessing ONT raw reads basecalled by MINKNOW; Doing fastq pooling per barcode, quality filtering, host depletion")

    # input directory
    parser.add_argument("--bd", dest="basecall_dir", required=True, type=str, help="basecalling output directory; the directory should be most likely /some/path/to/fastq_pass/, which can be found next to other directories such as fastq_fail/, pod5/, ... after running MINKNOW")

    # input barcode - sample name - sample type matching table
    parser.add_argument("--barcode_def", dest="barcode_definition_tsv", required=False, default = "_ns", type=str, help="[optional] [default = not provided, all barcodes will be regarded as 'clinical' samples, i.e., the patient samples] give a two-column tsv file, without header row, column 1 barcode number (e.g., 'barcode01'), column 2 sample name (e.g., 'patient_1_blood_240831'), column 3 sample type (possible values = 'clinical', 'negative', 'positive'); You can exclude some barcode numbers from the table in this file, then we will consider those barcode as not used, and will not create fastq file for those missing barcodes.")

    # threads and executable external programs
    parser.add_argument("--threads", dest="num_threads_str", required=False, default = "4", type=str, help="[optional] [default = 4] number of threads to use.")
    parser.add_argument("--minimap_path", dest="executable_minimap_path", required=False, type=str, default="minimap2", help="path to executable minimap2 [default = minimap2]")
    parser.add_argument("--nanoq_path", dest="executable_nanoq_path", required=False, type=str, default="nanoq", help="path to executable nanoq [default = nanoq]")
    parser.add_argument("--pigz_path", dest="executable_pigz_path", required=False, type=str, default="pigz", help="path to executable pigz [default = pigz]")

    # read quality filtering parameters
    parser.add_argument("--min_read_len", dest="min_read_length_threshold", required=False, default = 500, type=int, help="[default = 500] Minimum threshold read length, for raw read QC filtering")
    parser.add_argument("--min_read_qual", dest="min_read_quality_threshold", required=False, default=12, type=float, help="[default = 12] Minimum threshold mean quality score of a read, for raw read QC filtering")

    # host removal parameters
    parser.add_argument("--host_genome_mmi", dest="host_genome_mmi_file_path", required=True, type=str, help="host reference genome sequence minimap index file path")
    parser.add_argument("--minimap_x", dest="minimap_preset_to_use", required=False, type=str, default="map-ont", help="minimap preset (-x) to use during alignment; it is not effective if your reference DB minimap index is not built with this preset too. [default = map-ont]")
    parser.add_argument("--host_map_idpct", dest="min_aln_idpct_host_removal", required=False, default=80, type=float, help="[default = 80] Alignment percent identity cutoff, to remove the reads when it is aligned to host genome")
    parser.add_argument("--host_map_len", dest="min_aln_length_host_removal", required=False, default=100, type=float, help="[default = 100] Alignment length bp cutoff, to remove the reads then it is aligned to host genome")

    # output file path designation
    parser.add_argument("--outdir", dest="output_dir_path", required=True, type=str, help="output files directory (if directory does not exist, new directory will be created at this path)")

    # 2. Do parsing
    args = parser.parse_args()

    return args


def check_tool_availability(tool_name, command):
    is_normal = False
    what_happen = ""
    try:
        # Run the version command for the tool
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Check if the command was successful
        if result.returncode == 0:
            is_normal = True
            what_happen = tool_name + " runs fine"
        else:
            is_normal = False
            what_happen = tool_name + " is not working"
    except Exception as e:
        is_normal = False
        what_happen = tool_name + " checking failed"
    return is_normal, what_happen


def just_count_read_and_base_from_fastq(input_fastq_file):
    read_count = 0
    total_base_count = 0

    # Open the qc-passed reads fastq and just count the reads and bases
    open_func = gzip.open if input_fastq_file.endswith('.gz') else open
    with open_func(input_fastq_file, 'rt') as fastq_file:
        while True:
            # Read each set of four lines in a FASTQ format
            header = fastq_file.readline()
            if not header:
                break  # End of file
            sequence = fastq_file.readline().strip()
            plus = fastq_file.readline()  # skip '+'
            quality = fastq_file.readline()  # skip quality scores
            # Count sequence and bases
            read_count += 1
            total_base_count += len(sequence)

    return read_count, total_base_count


def find_tag_from_paf_columns(column_value_list, target_tag):
    tag_str = ''
    for colzbi in range(12, len(column_value_list)):
        if column_value_list[colzbi].startswith(target_tag):
            tag_str = column_value_list[colzbi]
            break
    return tag_str


def write_fastq_depleted_of_aligned_reads(input_fastq, output_fastq, alignment_paf, min_idpct, min_length):
    dict_target_readid = {}
    fr = open(alignment_paf, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        query_seq = ls[0]

        # for alignment length inspection, we use the reference / target sequence alignment length
        ref_aligned_length = int(ls[8]) + 1 - int(ls[7])
        if ref_aligned_length < min_length:
            continue
        # for identity inspection, we use dv:f tag as the divergence and convert it to identity percentage
        dvtag = find_tag_from_paf_columns(ls, 'dv:')
        dv_ident = 100*(1 - float(dvtag.split(':')[2]))
        if dv_ident < min_idpct:
            continue

        dict_target_readid[query_seq] = True
    fr.close()            
    num_target_read_id = len(dict_target_readid)

    n_kept = 0
    n_total = 0
    n_excl = 0

    fr = open(input_fastq, "r")
    fw = open(output_fastq, "w")
    line = fr.readline()
    while line != '':
        readid = line.strip()[1:].split(' ')[0]
        n_total += 1
        if readid.endswith('/1') or readid.endswith('/2'):
            readid = readid[:-2]
        if readid not in dict_target_readid:
            # read & write
            n_kept += 1
            fw.write(line.strip() + "\n")   # id
            line = fr.readline()
            fw.write(line.strip() + "\n")   # seq
            line = fr.readline()
            fw.write(line.strip() + "\n")   # separator
            line = fr.readline()
            fw.write(line.strip() + "\n")   # quality        
            
            line = fr.readline()
            # next iter (id)
        else:
            # just read
            n_excl += 1
            # id
            line = fr.readline()
            # seq
            line = fr.readline()
            # separator
            line = fr.readline()
            # quality
            
            line = fr.readline()
            # next iter (id)
    fr.close()
    fw.close()

    print("  ...  Total reads = " + str(n_total))
    print("  ...  Kept reads = " + str(n_kept))
    print("  ...  Removed reads = " + str(n_excl))


def pigz_compress_a_file(pigz_path, input_file, threads):
    pigz_command = []
    pigz_command.append(pigz_path)
    pigz_command.append("--force")
    pigz_command.append("--processes")
    pigz_command.append(threads)
    pigz_command.append(input_file)
    try:
        print("# Run: " + " ".join(pigz_command))
        subprocess.run(pigz_command, check = True)
        return True
    except subprocess.CalledProcessError as e:
        print("[mngs_ont_read_preprocess_from_basecall_directory.py] the subprocess was aborted, pigz decompression was not successful")
        return False




"""
Main script starts
"""
starting_time = datetime.now()
starting_time_fstr = starting_time.strftime("%Y-%m-%d %H:%M:%S")
print("[TIME LOG] start of the script mngs_ont_read_preprocess_from_basecall_directory.py = " + starting_time_fstr)


args = get_arguments_from_command_line()
# input fasta(.gz) file
basecall_dir = args.basecall_dir
barcode_definition_tsv = args.barcode_definition_tsv
num_threads_str = args.num_threads_str
executable_minimap_path = args.executable_minimap_path
executable_pigz_path = args.executable_pigz_path
executable_nanoq_path = args.executable_nanoq_path
min_read_length_threshold = args.min_read_length_threshold
min_read_quality_threshold = args.min_read_quality_threshold
host_genome_mmi_file_path = args.host_genome_mmi_file_path
minimap_preset_to_use = args.minimap_preset_to_use
min_aln_idpct_host_removal = args.min_aln_idpct_host_removal
min_aln_length_host_removal = args.min_aln_length_host_removal
output_dir_path = args.output_dir_path


# file paths and script parameters that can be derived from user inputs
basecall_fastq_pass_dir = basecall_dir # os.path.join(basecall_dir, "fastq_pass")
use_barcode_definition = False
if barcode_definition_tsv != "_ns":
    use_barcode_definition = True
output_barcode_information_tsv = os.path.join(output_dir_path, "barcode_information.tsv")
output_barcode_raw_fastq_dir = os.path.join(output_dir_path, "per_barcode_raw_fastq")
output_barcode_qcpass_fastq_dir = os.path.join(output_dir_path, "per_barcode_qcpass_fastq")
output_barcode_hostrm_fastq_dir = os.path.join(output_dir_path, "per_barcode_hostrm_fastq")
output_tarball = os.path.join(output_dir_path, "host_removed_fastq_and_sample_information.tar.gz")

# basic sanity check about user-provided inputs
if not os.path.isdir(basecall_dir):
    print("[Incorrect operation!] The input basecall result directory (--bd) is not found at the provided path: " + basecall_dir)
    quit()
if not os.path.isdir(basecall_fastq_pass_dir):
    print("[Incorrect operation!] The input basecall result directory (--bd) does not contain fastq_pass/ subdirectory")
    quit()
if not os.path.isfile(host_genome_mmi_file_path):
    print("[Incorrect operation!] The input host genome minimap index mmi file (--host_genome_mmi) is not found at the provided path: " + host_genome_mmi_file_path)
    quit()
minimap_working, what_happend_to_minimap = check_tool_availability("minimap2", f"{executable_minimap_path} --version")
pigz_working, what_happend_to_pigz = check_tool_availability("pigz", f"{executable_pigz_path} --version")
nanoq_working, what_happend_to_nanoq = check_tool_availability("pigz", f"{executable_nanoq_path} --version")
if not minimap_working:
    print("[Incorrect operation!] external command line tool minimap is not working with the given path: " + what_happend_to_minimap)
    quit()
if not pigz_working:
    print("[Incorrect operation!] external command line tool pigz is not working with the given path: " + what_happend_to_pigz)
    quit()
if not nanoq_working:
    print("[Incorrect operation!] external command line tool nanoq is not working with the given path: " + what_happend_to_nanoq)
    quit()

if use_barcode_definition:
    # 1. sample type assignment must make sense; possible / allowed values are strictly among these: 'clinical', 'negative', 'positive'
    # 2. at least one barcode name used in the definition table should be present as the subdirectory in basecall directory's fastq_pass/ directory.
    dict_allowed_sample_type = {'clinical': 0, 'negative': 1, 'positive': 2}
    list_barcode_dir_names_in_basecall = os.listdir(basecall_fastq_pass_dir)

    list_illegal_sample_type = []
    list_barcode_directory_present = []
    fr = open(barcode_definition_tsv, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        barcode = ls[0]
        sample_name = ls[1]
        sample_type = ls[2]
        if sample_type not in dict_allowed_sample_type:
            list_illegal_sample_type.append(sample_type)
        if barcode in list_barcode_dir_names_in_basecall:
            list_barcode_directory_present.append(barcode)
            print("# basecalled result direictory presence check, in fastq_pass/:  " + barcode + "\t" + sample_name + "\tyes")
        else:
            print("# basecalled result direictory presence check, in fastq_pass/:  " + barcode + "\t" + sample_name + "\tno")
    fr.close()
    if len(list_illegal_sample_type) > 0:
        print("[Incorrect operation!] Sample type value other than the allowed types (clinical, negative, positive) is used in your --barcode_def input table")
        quit()
    if len(list_barcode_directory_present) < 1:
        print("[Unexpected situation!] None of the user-specified barcodes is find in the --bd directory under fastq_pass/ subdirectory")
        quit()


# create directories for writing output files
if not os.path.isdir(output_dir_path):
    os.makedirs(output_dir_path)
if not os.path.isdir(output_barcode_raw_fastq_dir):
    os.makedirs(output_barcode_raw_fastq_dir)
if not os.path.isdir(output_barcode_qcpass_fastq_dir):
    os.makedirs(output_barcode_qcpass_fastq_dir)
if not os.path.isdir(output_barcode_hostrm_fastq_dir):
    os.makedirs(output_barcode_hostrm_fastq_dir)


# create data slots to fill with sample information and read statistics per barcode
list_barcode_dir_names_in_basecall = os.listdir(basecall_fastq_pass_dir)
list_barcodeno = []
dict_barcodeno_index = {}
buff_barcodeno_index = -1
list_barcodeno_samplename = []
list_barcodeno_sampletype = []

if use_barcode_definition:
    fr = open(barcode_definition_tsv, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        barcode = ls[0]
        sample_name = ls[1]
        sample_type = ls[2]
        if barcode not in list_barcode_dir_names_in_basecall:
            continue
        list_barcodeno.append(barcode)
        buff_barcodeno_index += 1
        dict_barcodeno_index[barcode] = buff_barcodeno_index
        list_barcodeno_samplename.append(sample_name)
        list_barcodeno_sampletype.append(sample_type)
    fr.close()
else:
    for dir_name in list_barcode_dir_names_in_basecall:
        if not dir_name.startswith('barcode'):
            continue
        list_barcodeno.append(dir_name)
        buff_barcodeno_index += 1
        dict_barcodeno_index[dir_name] = buff_barcodeno_index
        list_barcodeno_samplename.append(dir_name)
        list_barcodeno_sampletype.append("clinical")

n_barcodeno = len(list_barcodeno)
list_barcodeno_n_read_raw = [0]*n_barcodeno
list_barcodeno_total_bp_raw = [0]*n_barcodeno
list_barcodeno_n_read_qcpass = [0]*n_barcodeno
list_barcodeno_total_bp_qcpass = [0]*n_barcodeno
list_barcodeno_n_read_hostrm = [0]*n_barcodeno
list_barcodeno_total_bp_hostrm = [0]*n_barcodeno

list_barcodeno_raw_read_fastq = ['']*n_barcodeno
list_barcodeno_qcpass_read_fastq = ['']*n_barcodeno
list_barcodeno_hostrm_read_fastq = ['']*n_barcodeno
list_barcodeno_host_align_paf = ['']*n_barcodeno
list_barcodeno_sample_info_txt = ['']*n_barcodeno
for barcodeno_index in range(n_barcodeno):
    barcodeno = list_barcodeno[barcodeno_index]
    list_barcodeno_raw_read_fastq[barcodeno_index] = os.path.join(output_barcode_raw_fastq_dir, barcodeno + ".raw.fastq")
    list_barcodeno_qcpass_read_fastq[barcodeno_index] = os.path.join(output_barcode_qcpass_fastq_dir, barcodeno + ".qcpass.fastq")
    list_barcodeno_hostrm_read_fastq[barcodeno_index] = os.path.join(output_barcode_hostrm_fastq_dir, barcodeno + ".hostrm.fastq")
    list_barcodeno_host_align_paf[barcodeno_index] = os.path.join(output_barcode_hostrm_fastq_dir, barcodeno + ".host_genome_align.paf")
    list_barcodeno_sample_info_txt[barcodeno_index] = os.path.join(output_barcode_hostrm_fastq_dir, barcodeno + ".sample_info.txt")


"""
From each barcode, first create a pooled fastq file, not zipped, count the read number and total bases.
"""
midproc_time = datetime.now()
midproc_time_fstr = midproc_time.strftime("%Y-%m-%d %H:%M:%S")
time_flyen_sec = (midproc_time - starting_time).total_seconds()
time_flyen_sec_str = f"{time_flyen_sec:.3f}"
print("[TIME LOG] Now we start to create fastq per barcode: " + midproc_time_fstr + " | time passed from the script start (seconds) = " + time_flyen_sec_str)

for barcodeno_index in range(n_barcodeno):
    barcodeno = list_barcodeno[barcodeno_index]
    output_raw_read_fastq = list_barcodeno_raw_read_fastq[barcodeno_index]
    input_fastq_gz_filename_list = os.listdir(os.path.join(basecall_fastq_pass_dir, barcodeno))

    # Remove the output file if it already exists
    if os.path.isfile(output_raw_read_fastq):
        os.remove(output_raw_read_fastq)

    read_count = 0
    total_base_count = 0

    # Open the output file in append mode and process each fastq.gz file
    with open(output_raw_read_fastq, 'a') as outfile:
        for filename in input_fastq_gz_filename_list:
            input_fastq_gz_filepath = os.path.join(basecall_fastq_pass_dir, barcodeno, filename)

            # Determine if the file is gzipped
            open_func = gzip.open if input_fastq_gz_filepath.endswith('.gz') else open

            # Open the file and read four lines at a time
            with open_func(input_fastq_gz_filepath, 'rt') as fastq_file:
                while True:
                    # Read each set of four lines in a FASTQ format
                    header = fastq_file.readline()
                    if not header:
                        break  # End of file
                    sequence = fastq_file.readline().strip()
                    plus = fastq_file.readline()  # skip '+'
                    quality = fastq_file.readline()  # skip quality scores

                    outfile.write(header.strip() + "\n" + sequence + "\n" + "+" + "\n" + quality.strip() + "\n")

                    # Count sequence and bases
                    read_count += 1
                    total_base_count += len(sequence)

    list_barcodeno_n_read_raw[barcodeno_index] = read_count
    list_barcodeno_total_bp_raw[barcodeno_index] = total_base_count
    print("  ...  " + str(barcodeno_index + 1) + "/" + str(n_barcodeno) + " " + barcodeno + ": collected " + str(read_count) + " reads, " + str(total_base_count) + " bases")


"""
From each barcode, now filter the reads by length and mean quality.
"""
midproc_time = datetime.now()
midproc_time_fstr = midproc_time.strftime("%Y-%m-%d %H:%M:%S")
time_flyen_sec = (midproc_time - starting_time).total_seconds()
time_flyen_sec_str = f"{time_flyen_sec:.3f}"
print("[TIME LOG] Now we start to filter reads by length and quality per barcode: " + midproc_time_fstr + " | time passed from the script start (seconds) = " + time_flyen_sec_str)

for barcodeno_index in range(n_barcodeno):
    barcodeno = list_barcodeno[barcodeno_index]
    raw_read_fastq = list_barcodeno_raw_read_fastq[barcodeno_index]
    qcpass_read_fastq = list_barcodeno_qcpass_read_fastq[barcodeno_index]

    # nanoq -i ${input_fastq} --output-type u -l ${min_length} -q ${min_mean_q} -o ${output_fastq}
    nanoq_command = []
    nanoq_command.append(executable_nanoq_path)
    nanoq_command.append("-i")
    nanoq_command.append(raw_read_fastq)
    nanoq_command.append("-o")
    nanoq_command.append(qcpass_read_fastq)
    nanoq_command.append("--output-type")
    nanoq_command.append("u")
    nanoq_command.append("-l")
    nanoq_command.append(str(min_read_length_threshold))
    nanoq_command.append("-q")
    nanoq_command.append(str(min_read_quality_threshold))
    try:
        print("# Run: " + " ".join(nanoq_command))
        subprocess.run(nanoq_command, check = True)
    except subprocess.CalledProcessError as e:
        print("[mngs_ont_read_preprocess_from_basecall_directory.py] the subprocess was aborted, nanoq was not successful")
        exit(1)

    read_count, total_base_count = just_count_read_and_base_from_fastq(qcpass_read_fastq)
    list_barcodeno_n_read_qcpass[barcodeno_index] = read_count
    list_barcodeno_total_bp_qcpass[barcodeno_index] = total_base_count
    print("  ...  " + str(barcodeno_index + 1) + "/" + str(n_barcodeno) + " " + barcodeno + ": QC filtered the reads - to keep " + str(read_count) + " reads, " + str(total_base_count) + " bases")


"""
From each barcode, now remove the reads that map to host genome.
"""
midproc_time = datetime.now()
midproc_time_fstr = midproc_time.strftime("%Y-%m-%d %H:%M:%S")
time_flyen_sec = (midproc_time - starting_time).total_seconds()
time_flyen_sec_str = f"{time_flyen_sec:.3f}"
print("[TIME LOG] Now we start to remove host-aligned reads, per barcode: " + midproc_time_fstr + " | time passed from the script start (seconds) = " + time_flyen_sec_str)

for barcodeno_index in range(n_barcodeno):
    barcodeno = list_barcodeno[barcodeno_index]
    qcpass_read_fastq = list_barcodeno_qcpass_read_fastq[barcodeno_index]
    hostrm_read_fastq = list_barcodeno_hostrm_read_fastq[barcodeno_index]
    host_mapping_paf = list_barcodeno_host_align_paf[barcodeno_index]

    # minimap2 -x map-ont -t 4 -o ${output_paf} --secondary=no ${host_genome_mmi} ${input_fastq}
    minimap_command = []
    minimap_command.append(executable_minimap_path)
    minimap_command.append("-x")
    minimap_command.append(minimap_preset_to_use)
    minimap_command.append("-t")
    minimap_command.append(num_threads_str)
    minimap_command.append("-o")
    minimap_command.append(host_mapping_paf)
    minimap_command.append("--secondary=no")
    minimap_command.append(host_genome_mmi_file_path)
    minimap_command.append(qcpass_read_fastq)
    try:
        print("# Run: " + " ".join(minimap_command))
        subprocess.run(minimap_command, check = True)
    except subprocess.CalledProcessError as e:
        print("[mngs_ont_read_preprocess_from_basecall_directory.py] the subprocess was aborted, minimap against host genome was not successful")
        exit(1)

    # remove using paf the aligned reads from the fastq file
    write_fastq_depleted_of_aligned_reads(qcpass_read_fastq, hostrm_read_fastq, host_mapping_paf, min_aln_idpct_host_removal, min_aln_length_host_removal)

    read_count, total_base_count = just_count_read_and_base_from_fastq(hostrm_read_fastq)
    list_barcodeno_n_read_hostrm[barcodeno_index] = read_count
    list_barcodeno_total_bp_hostrm[barcodeno_index] = total_base_count
    print("  ...  " + str(barcodeno_index + 1) + "/" + str(n_barcodeno) + " " + barcodeno + ": Filtered host-aligned reads - " + str(read_count) + " reads, " + str(total_base_count) + " bases")


"""
Finally, gzip compress all fastq files
"""
midproc_time = datetime.now()
midproc_time_fstr = midproc_time.strftime("%Y-%m-%d %H:%M:%S")
time_flyen_sec = (midproc_time - starting_time).total_seconds()
time_flyen_sec_str = f"{time_flyen_sec:.3f}"
print("[TIME LOG] Now we start to compress all fastq output files: " + midproc_time_fstr + " | time passed from the script start (seconds) = " + time_flyen_sec_str)


for barcodeno_index in range(n_barcodeno):
    decomp_raw = pigz_compress_a_file(executable_pigz_path, list_barcodeno_raw_read_fastq[barcodeno_index], num_threads_str)
    decomp_qcpass = pigz_compress_a_file(executable_pigz_path, list_barcodeno_qcpass_read_fastq[barcodeno_index], num_threads_str)
    decomp_hostrm = pigz_compress_a_file(executable_pigz_path, list_barcodeno_hostrm_read_fastq[barcodeno_index], num_threads_str)
        

"""
Write barcode sample information
"""
fw = open(output_barcode_information_tsv, 'w')
fw.write("barcode\tsample_name\tsample_type\tn_reads_raw\tn_reads_qcpass\tn_reads_hostrm\ttotal_bp_raw\ttotal_bp_qcpass\ttotal_bp_hostrm\n")
for barcodeno_index in range(n_barcodeno):
    barcodeno = list_barcodeno[barcodeno_index]
    samplename = list_barcodeno_samplename[barcodeno_index]
    sampletype = list_barcodeno_sampletype[barcodeno_index]
    n_read_raw = list_barcodeno_n_read_raw[barcodeno_index]
    n_read_qcpass = list_barcodeno_n_read_qcpass[barcodeno_index]
    n_read_hostrm = list_barcodeno_n_read_hostrm[barcodeno_index]
    total_bp_raw = list_barcodeno_total_bp_raw[barcodeno_index]
    total_bp_qcpass = list_barcodeno_total_bp_qcpass[barcodeno_index]
    total_bp_hostrm = list_barcodeno_total_bp_hostrm[barcodeno_index]
    fw.write(barcodeno)
    fw.write("\t" + samplename)
    fw.write("\t" + sampletype)
    fw.write("\t" + str(n_read_raw))
    fw.write("\t" + str(n_read_qcpass))
    fw.write("\t" + str(n_read_hostrm))
    fw.write("\t" + str(total_bp_raw))
    fw.write("\t" + str(total_bp_qcpass))
    fw.write("\t" + str(total_bp_hostrm))
    fw.write("\n")
    #
    output_indiv_sample_info = list_barcodeno_sample_info_txt[barcodeno_index]
    fw_indiv = open(output_indiv_sample_info, 'w')
    fw_indiv.write(">" + barcodeno + "\n")
    fw_indiv.write(":\tsample_name\t" + samplename + "\n")
    fw_indiv.write(":\tsample_type\t" + sampletype + "\n")
    fw_indiv.write(":\tn_read_hostrm\t" + str(n_read_hostrm) + "\n")
    fw_indiv.write(":\ttotal_bp_hostrm\t" + str(total_bp_hostrm) + "\n")
    fw_indiv.write(":\tn_read_qcpass\t" + str(n_read_qcpass) + "\n")
    fw_indiv.write(":\ttotal_bp_qcpass\t" + str(total_bp_qcpass) + "\n")
    fw_indiv.write(":\tn_read_raw\t" + str(n_read_raw) + "\n")
    fw_indiv.write(":\ttotal_bp_raw\t" + str(total_bp_raw) + "\n")
    fw_indiv.close()
fw.close()


# Package a tarball that contains host-removed reads fastq.gz files and sample information txt files.
files_to_package = []
files_to_package.append(output_barcode_information_tsv)
for barcodeno_index in range(n_barcodeno):
    files_to_package.append(list_barcodeno_hostrm_read_fastq[barcodeno_index] + ".gz")
    files_to_package.append(list_barcodeno_sample_info_txt[barcodeno_index])

# Create a tarball and add each specified file
with tarfile.open(output_tarball, "w:gz") as tar:
    for file_path in files_to_package:
        # Use os.path.basename to get only the file name
        arcname = os.path.basename(file_path)
        tar.add(file_path, arcname=arcname)

print(f"Tarball created at {output_tarball}")


midproc_time = datetime.now()
midproc_time_fstr = midproc_time.strftime("%Y-%m-%d %H:%M:%S")
time_flyen_sec = (midproc_time - starting_time).total_seconds()
time_flyen_sec_str = f"{time_flyen_sec:.3f}"
print("[TIME LOG] Everything is completed: " + midproc_time_fstr + " | time passed from the script start (seconds) = " + time_flyen_sec_str)

