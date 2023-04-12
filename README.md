# ACP_bulkRNAseq

A pipeline to perform a simple PCA on RNA-seq data. 
Developped for Todorov et al. 2018 published in RIP journal Genomics and Computational Biology. 

	README.txt - 07.02.2018

	The zip file NextGenBender.zip contains the following files.
	Except installation.sh and TimeTable.ods, they have to be put 
	in the same folder:

	- Installation.sh
		The step-by-step installation of all tools used by 
		NextGenBender. Do not make it executable, but rather
		install step by step the tools that you are missing on
		your system. Some commands are not from sh but have to
		executed or typed in other environement like nano or R.
	- NextGenBender_Slurm.py
		The pipeline to use with the Slurm job scheduler. 
	- NextGenBender_LSF.py
		The pipeline to use in a LSF job scheduler. 
	- NextGenBender.sh
		The pipeline to use without scheduler. 
	- automated-fastqc2.sh
		The quality control step for the fastq files. 
	- NormEnrichments_parallel.R
		The final step of the pipeline that generates the information
		table.
	- wig.R
		R script to generate wiggle tracks for each bam file. 
	- TimeTable.ods
		File containing a comparison of duration of pipeline steps
		depending on the number of threads used for computation. 


