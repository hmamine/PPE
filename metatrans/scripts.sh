#!/bin/bash

###script to analyze RNAseq data of E.amylovora, Pseudomonas1079 & Pantoea1B4 in liquid cultures 
###written by M Amine Hassani - Feb 2022 - CAES New Haven.
 
#SBATCH --job-name=Sal
#SBATCH --partition=serial
#SBATCH --ntasks=12
#SBATCH --nodes=2
#SBATCH --time=5-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mohamed-amine.hassani@ct.gov
#SBATCH --output=Sal.out
#SBATCH --error=Sal.err

#loading modules & setting paths
module load anaconda 
source activate salmon
path_data=/scratch/quz19005/Amine/EPPZ/PPE
path_out=/home/quz19005/Amine/EPPZ

#salmon index -t /home/quz19005/Amine/EPPZ/Ea.ffn -i Ea_index
#salmon index -t /home/quz19005/Amine/EPPZ/Pa.ffn -i Pa_index
#salmon index -t /home/quz19005/Amine/EPPZ/Ps.ffn -i Ps_index

DB_Ea=/home/quz19005/Amine/EPPZ/Ea_index
DB_Ps=/home/quz19005/Amine/EPPZ/Ps_index
DB_Pa=/home/quz19005/Amine/EPPZ/Pa_index


ln -s /home/quz19005/Amine/Truseq3-all.fa Adapt
for file in *_R1.fastq.gz ; do
	base=$(basename ${file} _R1.fastq.gz) ; 
	echo " 1 ==> Trimming adaptor sequences from sample ${base} ...";
	trimmomatic PE $path_data/${base}_R1.fastq.gz $path_data/${base}_R2.fastq.gz\
	$path_data/${base}_R1_trim.fastq.gz $path_data/${base}_U1.fastq.gz\
	$path_data/${base}_R2_trim.fastq.gz $path_data/${base}_U2.fastq.gz\
	ILLUMINACLIP:Adapt:2:30:10\
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30;
	
	rm -f $path_data/${base}_U1.fastq.gz;
	rm -f $path_data/${base}_U2.fastq.gz

	echo " 2 ==> Mapping reads from ${base} to Ea genome ... "
	salmon quant -i $DB_Ea -l A -p 12\
	-1 $path_data/${base}_R1_trim.fastq.gz\
	-2 $path_data/${base}_R2_trim.fastq.gz\
	--validateMappings -o $path_out/${base};
	
	done


