# Import the following packages 
import os 
import Bio
from Bio import Entrez

# Create a new directory for the project
os.system("mkdir Multiple_Sclerosis_Kallisto")
os.chdir("Multiple_Sclerosis_Kallisto")

# Generate a log file
log_file = open("Multiple_Sclerosis_Kallisto.log", "w")
log_file.write("Log file generated.")
log_file.close()

# Retrive the reference transcriptome, Transgenic Mice, from Ensembl
os.system("wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz")

# Retrive Dr. Chen's Data from the CompBio class server 
Chen_data = "/home/data/Chen_data"

# Create the reference transcriptome path
# NOTE: users will have to enter their Loyola UNIVID used to access the class server 
Ref_T_path = "/home/<yourUNIVID>/Multiple_Sclerosis_Kallisto"

# Create the Index 
os.system(f"time kallisto index -i MS.idx Mus_musculus.GRCm39.cdna.all.fa.gz")

# Set the index path to be used later in the Kallisto quantification step
index_path = "MS.idx"

# Define the sample groups (from Dr. Chen's data) and their corresponding FASTQ files
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
# Store the reuslts in a folder called "Kallisto_results"

output_folder = "kallisto_results"

for sample, fastq_file in samples:
  os.system(f"mkdir -p {output_folder}/{sample}")
os.system(f"kallisto quant -i {index_path} -o {output_folder}/{sample} --single -l 200 -s 20 -b 100 -t 4 {Chen_data}/{fastq_file}")

# The next step is to run Sleuth on the kallisto output files in an R script.


