#!/bin/bash

#Add options to define rnaspades or metaspades upon user input

for srr in "$@"
do
	echo "Downloading $srr now."
	mkdir $srr
	echo "fastq-dump -I --split-files --outdir ./$srr $srr"
	#######
	echo "Quality control of $srr."
	echo "trim_galore -o ./${srr} --paired ./${srr}/${srr}_1.fastq ./${srr}/${srr}_2.fastq"
	#######
	echo "Assembly of $srr"
	#mkdir ./${srr}/metaspades_out
	mkdir ./${srr}/rnaspades_out
	#echo "metaspades.py -o ./${srr}/metaspades-out/ -1 ./${srr}/${srr}_1_val_1.fq -2 ./${srr}/${srr}_2_val_2.fq"
	echo "rnaspades.py -o ./${srr}/rnaspades-out/ -1 ./${srr}/${srr}_1_val_1.fq -2 ./${srr}/${srr}_2_val_2.fq"
	#######
	echo "Assembly evaluation of $srr"
	mkdir ./${srr}/quast_check
	mkdir ./${srr}/metaquast_check
	echo "metaquast.py --rna-finding -o ./${srr}/metaquast_check ./${srr}/rnaspades_out/transcripts.fasta"
	echo "quast.py --rna-finding -o ./${srr}/quast_check ./${srr}/rnaspades_out/transcripts.fasta"
done


#qsub -I -q bioforce-6 -l nodes=8:ppn=2,walltime=15:00:00,pmem=2gb