#!/bin/bash
echo "script start: download and initial sequencing read quality control"
#Import: It is assumed that the script is run from the analyses subdirectory of the gitrepo.
date

#Start a screen sessions:
screen

#Make a subdirectory for RAW FASTQ files:
mkdir ../data/sra_fastq

#export list of sequencing run identifiers in a file which is used to download the corresponding FASTQ files:
sqlite3 -batch -noheader -csv  /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "SELECT run_accession FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b USING(patient_code) WHERE username=='mtronstad';" > "mtronstad_run_accession.txt"

#Load sra-tools:
module load sra-tools

#download corresponding FASTQ-files:
cat mtronstad_run_accession.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip --outdir ../data/sra_fastq/ --disable-multithreading --split-e
date

#Check the number of files:
ls ../data/sra_fastq | wc -l

#Print the number of reads for each file in sra_fastq corresponding to "my" samples:
for f in $(ls ../data/sra_fastq/); do srun --cpus-per-task=1 --time=00:1:00 zgrep -c  "@" "../data/sra_fastq/$f" >> "number_of_reads.txt"; done

#Print meta_data (run_accession, patient_code, total_reads, read_count, base_count) for "my" samples
sqlite3 -batch /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "SELECT run_accession, patient_code, total_reads, read_count, base_count FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b USING(patient_code) WHERE username=='mtronstad';"

#Load seqkit module:
module load seqkit

#Check for duplicates, and store cleaned files in a new directory:
for f in $(ls ../data/sra_fastq/); do srun --cpus-per-task=1 zcat ../data/sra_fastq/$f | seqkit --threads 1 rmdup -s -o /cleaned_fastq_files/clean.$f; done

#Check for matches on AGATCGGAAGAGCACACGTCTGAACTCCAGTCA:
for f in $(ls ../data/sra_fastq/); do srun --cpus-per-task=1 zcat ../data/sra_fastq/$f | seqkit --threads 1 grep -s -i -p "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" | seqkit --threads 1 stats; done

#Check for matches on AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT:
for f in $(ls ../data/sra_fastq/); do srun --cpus-per-task=1 zcat ../data/sra_fastq/$f | seqkit --threads 1 grep -s -i -p "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" | seqkit --threads 1 stats; done

#Load FastQC module:
module load fastqc

#Create QC reports with FastQC:
srun --cpus-per-task=2 --time=00:10:00 xargs -I{} -a mtronstad_run_accession.txt fastqc --outdir ./fastqc/ --threads 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz

#Copy files onto local drive from remote machine:
scp mtronstad@core.cluster.france-bioinformatique.fr:/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/fastqc/*.html  /Users/magnus/Library/CloudStorage/OneDrive-KarolinskaInstitutet/project/documents/courses/med_bio_info/applied_bioinformatics/pre_course_assignments/fastqc_html_files/

#Merge forward and reverse reads:
srun --cpus-per-task=2 --time=00:30:00 xargs -a mtronstad_run_accession.txt -n 1 -I{} flash2 --threads=2 --max-overlap=150 -z --output-directory=../data/merged_pairs/ --output-prefix={}.flash ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a mtronstad_flash2.log
#61.63% combined and number of seqs went from 884,598 to 545,214 on the ERR6913179 files

#Load the bowtie2 module:
module load bowtie2

#Download GENBANK:
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

#Look at a particular sequence:
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

#Create a bowtie2 indexed database from the reference squence
mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

#Align the FASTQ merged reads against the above bowtie2 index DB:
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/mtronstad_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/mtronstad_bowtie_merged2PhiX.log
#This resulted in 0% alignments.

#Look at a particular sequence:
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna

#Create bowti2 indexed database from the SC2 sequence:
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB

#Align the FASTQ merged reads against the above bowtie2 index DB:
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz -S bowtie/mtronstad_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/mtronstad_bowtie_merged2SC2.log
#This resulted in:
#7326629 reads; of these:
#  7326629 (100.00%) were unpaired; of these:
#    7318620 (99.89%) aligned 0 times
#    8009 (0.11%) aligned exactly 1 time
#    0 (0.00%) aligned >1 times
#0.11% overall alignment rate
#I.e., hits on SARS-Cov2!

#Load multiqc module
module load multiqc

#Create combined report for all samples analyzed, including a .HTML file:
srun multiqc --force --title "mtronstad sample sub-set" ../data/merged_pairs/ ./fastqc/ ./mtronstad_flash2.log ./bowtie/

#Copy files onto local drive from remote machine:
scp mtronstad@core.cluster.france-bioinformatique.fr:/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/mtronstad-sample-sub-set_multiqc_report.html  /Users/magnus/Library/CloudStorage/OneDrive-KarolinskaInstitutet/project/documents/courses/med_bio_info/applied_bioinformatics/pre_course_assignments/multiqc_report_data/

date