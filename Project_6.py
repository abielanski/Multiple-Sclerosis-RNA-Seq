import os 
import Bio
from Bio import Entrez

# 1
# Create a new directory for the project
os.system("mkdir Project_6")
os.chdir("Project_6")

# Generate the log file
log_file = open("Project_6.log", "w")
log_file.write("Log file generated.")
log_file.close()

# Retrieving the reference transcriptome of transgenic mice
os.system("wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz")

# Running Kallisto and Retrieving the Dr. Chen's Data
Chen_data = "/home/data/Chen_data"

# Reference transcriptome path
Ref_T_path = "/home/amishra1/Project_6"

# Indexing command 
os.system(f"time kallisto index -i MS.idx Mus_musculus.GRCm39.cdna.all.fa.gz")

# Set the index path to be used later in the Kallisto quantification step
index_path = "MS.idx"

# Define the sample groups and their corresponding FASTQ files
# Each tuple represents a sample with the following structure:
# (sample_name, fastq_file)
samples = [
    ("CFA_1", "YC-C1_S5_L004_R1_001.fastq.gz"),
    ("CFA_2", "YC-C2_S6_L004_R1_001.fastq.gz"),
    ("CFA_3", "YC-C3_S7_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_1", "YC-V1_S1_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_2", "YC-V2_S2_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_3", "YC-V3_S3_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_1", "YC-S3_S4_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_2", "YC-4_S8_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_3", "YC-5_S9_L004_R1_001.fastq.gz"),
]

# Set the Chen data path
Chen_data = "/home/data/Chen_data"

# For each sample, create a directory for the sample (if it doesn't exist),
# and run Kallisto quantification using the provided index and FASTQ files
output_folder = "kallisto_results"
for sample, fastq_file in samples:
    os.system(f"mkdir -p {output_folder}/{sample}")
    os.system(f"kallisto quant -i {index_path} -o {output_folder}/{sample} --single -l 200 -s 20 -b 100 -t 4 {Chen_data}/{fastq_file}")

# Next step is to run sleuth on the kallisto output files in an R script.



