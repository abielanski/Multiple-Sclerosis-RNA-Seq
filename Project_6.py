import os 
import Bio
from Bio import Entrez

# 1
# Create a new directory for the project
#os.system("mkdir Project_6")
#os.chdir("Project_6")

# Generate the log file
log_file = open("Project_6.log", "w")
log_file.write("Log file generated.")
log_file.close()

# Retrieving the reference transcriptome of transgenic mice
#os.system("wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz")

# Running Kallisto and Retrieving the Dr. Chen's Data
Chen_data = "/home/data/Chen_data"

# Reference transcriptome path
Ref_T_path = "/home/amishra1/Project_6"

# Indexing command 
#os.system(f"time kallisto index -i MS.idx Mus_musculus.GRCm39.cdna.all.fa.gz")

# Set the index path to be used later in the Kallisto quantification step
index_path = "MS.idx"

# Define the sample groups and their corresponding FASTQ files
#Each tuple represents a sample with the following structure:
#(group_name, sample_name, reference_transcriptome_file, fastq_file)
samples = [
    ("CFA_group", "CFA_1", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-C1_S5_L004_R1_001.fastq.gz"),
    ("CFA_group", "CFA_2", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-C2_S6_L004_R1_001.fastq.gz"),
    ("CFA_group", "CFA_3", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-C3_S7_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_group", "EAE_Vehicle_1", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-V1_S1_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_group", "EAE_Vehicle_2", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-V2_S2_L004_R1_001.fastq.gz"),
    ("EAE_Vehicle_group", "EAE_Vehicle_3", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-V3_S3_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_group", "EAE_Sephin1_1", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-S3_S4_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_group", "EAE_Sephin1_2", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-4_S8_L004_R1_001.fastq.gz"),
    ("EAE_Sephin1_group", "EAE_Sephin1_3", "Mus_musculus.GRCm39.cdna.all.fa.gz", "YC-5_S9_L004_R1_001.fastq.gz"),
]

# For each sample, create a directory for the group (if it doesn't exist),
# and run Kallisto quantification using the provided index and FASTQ files
for group, sample, file_name_R1, file_name_R2 in samples:
    os.system(f"mkdir -p {group}")
    os.system(f"kallisto quant -i {index_path} -o {group}/{sample} -b 100 -t 4 {Ref_T_path}/{file_name_R1} {Chen_data}/{file_name_R2}")


