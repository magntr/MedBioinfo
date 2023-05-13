#!/bin/bash
echo "script start: download and initial sequencing read quality control"
date

#export list of sequencing run identifiers in a file which is used to download the corresponding FASTQ files:
sqlite3 -batch -noheader -csv  /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "SELECT run_accession FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b USING(patient_code) WHERE username=='mtronstad';" > "mtronstad_run_accession.txt"

#download corresponding FASTQ-files:
cat mtronstad_run_accession.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip --outdir ../data/sra_fastq/ --disable-multithreading --split-e
date

#Print the number of reads for each file in sra_fastq:
for f in $(ls ../data/sra_fastq/); do srun --cpus-per-task=1 --time=00:1:00 zgrep -c  "@" "../data/sra_fastq/$f" >> "number_of_reads.txt"; done

#Print meta_data (run_accession, patient_code, total_reads, read_count, base_count) for "my" samples
sqlite3 -batch /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "SELECT run_accession, patient_code, total_reads, read_count, base_count FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b USING(patient_code) WHERE username=='mtronstad';"

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