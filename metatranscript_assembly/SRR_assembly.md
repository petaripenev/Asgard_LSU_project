# SRR assembly tutorial

## Download a SRR into forward and reverse paired reads:

	fastq-dump -I --split-files SRR{number}

This generates two files - SRR{number}_1.fastq and SRR{number}_2.fastq

## Quality Control
You can check initial fastq files with **FastQC** or directly run **Trim Galore** to remove adapters and low quality regions.

	trim_galore --paired SRR{number}_1.fastq SRR{number}_2.fastq

This will generate two validated .fq files SRR{number}_1_val_1.fq and SRR{number}_2_val_2.fq;
You can run FastQC again to check what has improved or add the option --fastqc to trim_galore command to automatically run on validated reads.

## Assembly with **metaSPAdes** or **rnaSPAdes**
Requires a lot of space, ideally run on server/cluster. **rnaSPAdes** seems to be better.

	mkdir metaspades-out
	mkdir rnaspades-out
	metaspades.py -o metaspades-out/ -1 SRR{number}_1_val_1.fq -2 SRR{number}_2_val_2.fq
	rnaspades.py -o rnaspades-out/ -1 SRR{number}_1_val_1.fq -2 SRR{number}_2_val_2.fq


## Assembly evaluation
**MetaQUAST** or **QUAST**, **rnaQUAST** is not as good

	mkdir metaquast_check
	mkdir quast_check
	metaquast.py --rna-finding -o ./metaquast_check contigs.fasta
	quast.py --rna-finding -o ./quast_check transcripts.fasta

Might need to add option --max-ref-number 0 for metaquast, since often all of the references it will find do not exist in NCBI's assembly database. But always run without this option first to do BLAST detection.
It would be good to add our own reference.