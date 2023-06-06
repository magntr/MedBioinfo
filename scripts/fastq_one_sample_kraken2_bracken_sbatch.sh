#!/bin/bash
#
#SBATCH --partition=fast             # long, fast, etc.
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=1            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=90GB                    # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=0-01:00               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o ./outputs/slurm.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names)
#SBATCH -e ./outputs/slurm.%A.%a.err   # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
# /!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
## SBATCH --array=1-10                 # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
#SBATCH --job-name=kraken2_bracken_one_sample      # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
## SBATCH --mail-type END            # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user magnus.tronstad@ki.se

#################################################################
# Preparing work (cd to working dir, get hold of input data, convert/un-compress input data when needed etc.)
workdir_kraken2="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/kraken2"
workdir_bracken="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/bracken"
datadir_kraken2="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/sra_fastq"
datadir_bracken=${workdir_kraken2}

input_file1_kraken2=${datadir_kraken2}/ERR6913276_1.fastq.gz
input_file2_kraken2=${datadir_kraken2}/ERR6913276_2.fastq.gz

output_file_kraken2=${workdir_kraken2}/fastq_one_sample_kraken2

output_file_bracken=${workdir_bracken}/fastq_one_sample_bracken

db="/shared/projects/2314_medbioinfo/kraken2/arch_bact_vir_hum_protoz_fung/"


echo START: `date`

module load kraken2 #as required
module load bracken

mkdir -p ${workdir_kraken2}      # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)
mkdir -p ${workdir_bracken}      # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)

####################################################################################################
# Start work

srun --job-name=kraken2_one_sample kraken2 -db ${db} --threads ${SLURM_CPUS_PER_TASK} --paired ${input_file1_kraken2} ${input_file2_kraken2} --output ${output_file_kraken2}_out --report ${output_file_kraken2}_report

srun --job-name=bracken_one_sample bracken -d ${db} -i ${output_file_kraken2}_report -o ${output_file_bracken}_out  -w ${output_file_bracken}_report

####################################################################################################


echo END: `date`