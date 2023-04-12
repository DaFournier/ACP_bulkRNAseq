# Installation manual for tools used by NexGenBender.
# David Fournier, Robert Deelen 2018
 
# Instruction: 
# Attention! Not to be launched as an executable script. 
# Some commands have to be executed inside of specific tools (R, nano). 

# some built essentials for python
sudo apt-get update
sudo apt-get install build-essential checkinstall
sudo apt-get install libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev

# python 2.7 (used by MACS2 only)
cd /usr/src
sudo wget https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz
sudo tar xzf Python-2.7.14.tgz
cd Python-2.7.14
sudo ./configure
sudo make altinstall

# python 3 
sudo apt-get install python3.6

# "pip" command for python-based tools installs
sudo apt-get install python-pip python-dev build-essential 

# installation of MACS2 -- python 2.7 required
sudo apt-get install python-numpy # numpy library required
pip install MACS2

# sratoolkit
cd /home/big_hd/software
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz

# samtools
#  method 1
sudo apt-get install samtools
#  method 2
wget https://sourceforge.net/projects/samtools/files/samtools/1.6/samtools-1.6.tar.bz2 # version 1.6
tar -xjvf samtools-1.6.tar.bz2
cd samtools-1.6
sudo make
sudo make prefix=/usr/local/bin install

# bowtie2
sudo apt-get install bowtie2

# R
sudo apt-get install r-base

# QuasR library
sudo apt-get install libssl-dev
sudo apt-get install r-cran-rmysql
sudo apt-get install libxml2-dev
sudo apt-get install r-bioc-shortread
R
  source("https://bioconductor.org/biocLite.R") # two R commands
  biocLite("QuasR")



## software for quality control

#fastqc
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip
sudo nano /home/username/.bashrc # may require this language settings
 export LC_CTYPE=en_US.UTF-8 # to add at the end of .bashrc
 export LC_ALL=en_US.UTF-8

# cutadapt
sudo pip install cutadapt

# trim Galore!
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip
sudo unzip 0.4.5.zip




