import os
import subprocess
from pathlib import Path
import numpy as np # type: ignore
import re

HUMAS_BASH = "/g/korbel/hain/centromere/humas/HumAS-HMMER_for_AnVIL/hmmer-run.sh"
HUMAS_HMM = "/g/korbel/hain/centromere/humas/HumAS-HMMER_for_AnVIL/AS-HORs-hmmer3.4-071024.hmm"
STV_BASH = "/g/korbel/hain/centromere/stv/stv/stv.sh"

length_to_color = {
        1:(109,50,110),
        2:(211,105,80),
        3:(195,95,90),
        4:(240,220,120),
        5:(125,190,175),
        6:(165,36,90),    
        7:(144,196,155),
        8:(137,63,137),
        9:(100,100,169),    
        10:(150,195,135),
        11:(169,61,95),
        12:(55,150,195),
        13:(40,160,205),
        14:(92,175,165),
        15:(55,160,155),
        16:(66,178,205),
        17:(163,60,95),
        18:(167,25,95),
        19:(60,100,160),
        20:(190,220,150),
        21:(190,215,115),
        22:(225,225,120),
        23:(125,178,93),
        24:(180,80,95),
        25:(64,224,208),
        26:(247,225,149),
        30:(226,209,164)
    }

def _run_command(command):
    # Execute an external command and fail immediately if it returns a non-zero exit code.
    subprocess.run(command, check=True)
    
def remove_file_if_exists(path):
	# Delete temporary files only when present to avoid cleanup errors.
	if os.path.exists(path):
		os.remove(path)

def expand_repetitions(s):
    # replace repetition patterns like (AB){3} with ABABAB
    pattern = re.compile(r'\(([^)]+)\)\{(\d+)\}')

    def replacer(match):
        content = match.group(1)
        times = int(match.group(2))
        return content * times

    return pattern.sub(replacer, s)

def replace_mixed_monomers (monomer):
    # replace certain mixed monomer notations with a single monomer
    return monomer.replace("6/4", "6").replace("2/6", "6").replace("2/4", "6").replace("6/2", "6")

def count_monomers (hor, chromosome):

    working_hor = hor.split(".")[1]
    num_monomers = 0

    for block in expand_repetitions(working_hor).split("_"):
        # no minus
        if len(block.split("-")) == 1:
            num_monomers += 1
        # multiple ordered monomers
        elif len(block.split("-")) == 2:
            start = np.amin([int(replace_mixed_monomers(block.split("-")[0])), int(replace_mixed_monomers(block.split("-")[1]))])
            end = np.amax([int(replace_mixed_monomers(block.split("-")[0])), int(replace_mixed_monomers(block.split("-")[1]))])
            num_monomers += end - start + 1
        else:
            print(f"Cannot count HOR: {hor} - {block}")
            
    ### fix for chr1, chr19
    if chromosome in ["chr1", "chr19"] and ("{" in working_hor or num_monomers > 30):

        working_hor = expand_repetitions(working_hor)
        working_hor = working_hor.replace("6/4", "6").replace("2/6", "6").replace("2/4", "6").replace("(", "").replace(")", "")
        working_hor = working_hor.replace("5-6", "5_6")
        working_hor = working_hor.replace("4-5", "4_5")

        # test for 4_5 repeat
        tmp = working_hor.replace("4_5", "A")
        tmp = tmp.replace("_", "")
        tmp_4_5_frac = tmp.count("A") / len(tmp)

        # test for 5_6 repeat
        tmp = working_hor.replace("5_6", "A")
        tmp = tmp.replace("_", "")
        tmp_5_6_frac = tmp.count("A") / len(tmp)

        if tmp_4_5_frac > 0.8:
            num_monomers = 2
        elif tmp_5_6_frac > 0.8:
            num_monomers = 2
        else:
            num_monomers = len(working_hor.split("_"))

    return num_monomers

def humas_annotate (fasta_path:str, output_dir:str, sample_name:str, chromosome:str):
    
    ### move fasta file to working directory
    workdir = f"{output_dir}/{sample_name}/" 
    _run_command(["mkdir", "-p", workdir])
    _run_command(["cp", fasta_path, workdir])
    
    # run humas
    _run_command([HUMAS_BASH, workdir, HUMAS_HMM, "16"])

    # build output file names
    humas_output_file = f"AS-HOR-vs-{Path(fasta_path).stem}.bed"
    humas_plus_sf_output_file = f"AS-HOR+SF-vs-{Path(fasta_path).stem}.bed"
    
    # move output files into workdir
    _run_command(["mv", humas_output_file, workdir])
    
    # run stv
    _run_command([STV_BASH, f"{workdir}/{humas_output_file}"])
    
    # build name of output file
    outname = f"{output_dir}/{sample_name}_{chromosome}.contig_harmonized.bed"
    
    # rename contig into chromosome AND recolor HOR calls by number of monomers
    with open(outname, "w") as fo:
        for line in open("stv.bed"):
            Z = line.rstrip().split("\t")
            Z[0] = chromosome
            hor_length = count_monomers (Z[3], chromosome)
            hor_color = (0,0,0)
            if hor_length in length_to_color:
                hor_color = length_to_color[hor_length]
            Z[-1] = ",".join([str(x) for x in hor_color])
            Z[4] = str(hor_length)
            fo.write("\t".join(Z) + "\n")
            
    # remove workdir
    _run_command(["rm", "-r", workdir])
            
    # remove unused files
    remove_file_if_exists(humas_plus_sf_output_file)
    remove_file_if_exists("stv.bed")
    remove_file_if_exists("stv_raw.bed")
    remove_file_if_exists("stv_stats.table")
    remove_file_if_exists("stv_stats.tsv")

    
"""
conda create -n XYZ -c conda-forge -c bioconda python=3.11 hmmer bedtools bedops
conda activate centromeres
"""
    

if __name__ == "__main__":
    humas_annotate (
        fasta_path="/scratch/hain/humas_test/test_chr3.fa", 
        output_dir="/scratch/hain/humas_test/output", 
        sample_name="CHM13", chromosome="chr3")