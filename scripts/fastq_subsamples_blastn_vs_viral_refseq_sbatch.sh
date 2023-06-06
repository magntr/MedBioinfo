#!/bin/bash
#
#SBATCH --partition=fast             # long, fast, etc.
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=1            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=2GB                    # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=0-00:15               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o ./outputs/slurm.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names)
#SBATCH -e ./outputs/slurm.%A.%a.err   # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
# /!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
## SBATCH --array=1-8                # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
#SBATCH --job-name=blastn_viral       # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
## SBATCH --mail-type END              # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user magnus.tronstad@ki.se

#################################################################
# Preparing work (cd to working dir, get hold of input data, convert/un-compress input data when needed etc.)
workdir="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/analyses/blastn"
datadir="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/merged_pairs"
input_file1="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/merged_pairs/ERR6913277.flash.extendedFrags_100.fna"
input_file2="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/merged_pairs/ERR6913277.flash.extendedFrags_1000.fna"
input_file3="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/merged_pairs/ERR6913277.flash.extendedFrags_10000.fna"
output_file=${workdir}/fastq_subsample_blastn_vs_refseq_viral
db="/shared/projects/2314_medbioinfo/magnus/MedBioinfo/data/blast_db/refseq_viral_genomic"

echo START: `date`

module load blast #as required

mkdir -p ${workdir}    # -p because it creates all required dir levels **and** doesn't throw an error if the dir exists :)
cd ${workdir}

###################################################################################################
# Start work

srun --job-name=blastn_100sq blastn -num_threads ${SLURM_CPUS_PER_TASK} -query ${input_file1} -db ${db} -evalue 10 -max_target_seqs 5 -outfmt 6 -perc_identity 70 -out ${output_file}_100.tsv
srun --job-name=blastn_1000sq blastn -num_threads ${SLURM_CPUS_PER_TASK} -query ${input_file2} -db ${db} -evalue 10 -max_target_seqs 5 -outfmt 6 -perc_identity 70 -out ${output_file}_1000.tsv
srun --job-name=blastn_1000sq blastn -num_threads ${SLURM_CPUS_PER_TASK} -query ${input_file3} -db ${db} -evalue 10 -max_target_seqs 5 -outfmt 6 -perc_identity 70 -out ${output_file}_10000.tsv

####################################################################################################
echo END: `date`