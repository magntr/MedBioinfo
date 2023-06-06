#!/bin/bash
#
#SBATCH --partition=fast             # long, fast, etc.
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=4            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=80GB                    # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=0-01:00               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o ./outputs/slurm.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names)
#SBATCH -e ./outputs/slurm.%A.%a.err   # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
# /!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
#SBATCH --array=1-10                 # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
#SBATCH --job-name=global_kraken2_bracken_krona  # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
## SBATCH --mail-type END            # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user magnus.tronstad@ki.se

#################################################################
# Preparing work (cd to working dir, get hold of input data, convert/un-compress input data when needed etc.)
accnum_file="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/mtronstad_run_accession.txt"
workdir_kraken2="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/kraken2"
workdir_bracken="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/bracken"
workdir_krona="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/krona"
datadir_kraken2="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/sra_fastq"
datadir_bracken=${workdir_kraken2}

accnum=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${accnum_file})
echo Accnum is ${accnum}

input_file1_kraken2="${datadir_kraken2}/${accnum}_1.fastq.gz"
input_file2_kraken2="${datadir_kraken2}/${accnum}_2.fastq.gz"

output_file_kraken2=${workdir_kraken2}/${accnum}_fastq_full_kraken2
output_file_bracken=${workdir_bracken}/${accnum}_fastq_full_bracken

db="/shared/projects/2314_medbioinfo/kraken2/arch_bact_vir_hum_protoz_fung/"


echo START: `date`

module load kraken2 #as required
module load bracken

mkdir -p ${workdir_kraken2}      # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)
mkdir -p ${workdir_bracken}      # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)

# This extracts the item number SLURM_ARRAY_TASK_ID from the file of accnums


####################################################################################################
# Start work

## Kraken2
srun --job-name=${accnum} kraken2 -db ${db} --threads ${SLURM_CPUS_PER_TASK} --paired ${input_file1_kraken2} ${input_file2_kraken2} --output ${output_file_kraken2}_out --report ${output_file_kraken2}_report

## Bracken
srun --job-name=bracken_full bracken -d ${db} -i ${output_file_kraken2}_report -o ${output_file_bracken}_out  -w ${output_file_bracken}_report -r 50 -l S -t 5

## Convert bracken report to krona
srun --job-name=bracken2krona /shared/projects/2314_medbioinfo/kraken2/KrakenTools/kreport2krona.py -r ${output_file_bracken}_report -o ${workdir_krona}/${accnum}_full_bracken_to_krona_output

## Remove prefixes:
srun --job-name=ksed sed 's/[a-z]__//g' ${workdir_krona}/${accnum}_full_bracken_to_krona_output > ${workdir_krona}/${accnum}_full_bracken_to_krona_output_without_prefix

## Create Krona html
srun --job-name=krona2html /shared/projects/2314_medbioinfo/kraken2/bin/ktImportText -o ${workdir_krona}/text.${accnum}_without_prefix.html ${workdir_krona}/${accnum}_full_bracken_to_krona_output_without_prefix

####################################################################################################

echo END: `date`