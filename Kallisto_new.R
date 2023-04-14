# Generate the log file
log_file = open("Project_6.log", "w")
log_file.write("Log file generated.")
log_file.close()

# to create an index of the reference transciptome 
# the index gets stored in MS.idx
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
time kallisto index -i MS.idx Mus_musculus.GRCm39.cdna.all.fa.gz
index_path = "MS.idx"

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

# Quantification 
mkdir -p {output_folder}/{sample}
kallisto quant -i {index_path} -o {output_folder}/{sample} --single -l 200 -s 20 -b 100 -t 4 {Chen_data}/{fastq_file}