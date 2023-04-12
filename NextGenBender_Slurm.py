import csv
import time
import subprocess
from datetime import datetime
from os import path
import os
import shlex
import locale
import threading

'''
The input file has to be in the CSV format and named "input.csv". The file has to be made of 3 columns: column 1 is the name of the feature, column 2 the filename of the corresponding sequencing file (SRA or Fastq format), and column 3 contains an optional URL to download the file 
feature     filename    URL
example input:
EZH2,SRR568340.sra,https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR568/SRR568340/SRR568340.sra

To start the pipeline from any point, change the pipeline counter value (0: full pipeline with download,
 2: without download, 4: without fastqDump (if you have sam files), 5: samtools View, 6: merging for files of the same feature, 7: samtoolsSort, 8: samtoolsIndex,
 9: peak calling with MACS2, 10: creates the masterFile, 11: only read conting)

Personalized settings
after each value in " " hast to be a space. npt, cpt, download_delay, and pipeline_counter have to be integers.
'''
pipeline_counter = 0                                                                #startpoint
download_delay = 1                                                                  #seconds for each download
npt = "1 "                                                                          #nodes per task
cpt = "8 "                                                                          #cpu per task
tasks = 8                                                                           #number of parallel tasks
pipe = True                                                                         #pipes bowtie directly to samtools
delFastQ = True                                                                     #deletion of fastq files 
delSam = True                                                                       #deletion of SAM files
quality = True                                                                      #quality check
enrichment = True                                                                   #enrichment step 
reference = "/dir/to/ref/genome/hg19 "                                              #Path to the reference genome for bowtie
folder = ""                                                                         # working directory, leave empty if the same as script folder
fastqDump = "/dir/to/sratoolkit/bin/fastq-dump "                                     #path to fastqDump
bowtie = "bowtie2 -p " + cpt + " -x " + reference                                   #bowtie2 settings
samtoolsView = "samtools view -b -S "                                               #samtools settings
samtoolsSort = "samtools sort "
samtoolsMerge = "samtools merge "
samtoolsIndex = "samtools index "
macs2 = "macs2 callpeak -t "                                                        #macs parameters
macsModule = "module load bio/MACS2/2.1.0.20150731-foss-2017a-Python-2.7.13\n"      #module required for MACS2
bedtools = "/dir/to/bedtools2/bin/"                                                 #path to bedtools


def submit(caller, ignore_errors = False):
    """function for calling tools through slurm"""

    n_threads = os.environ['SLURM_CPUS_PER_TASK']
    os.environ['OMP_NUM_THREADS'] = n_threads
    if ("-o " in caller):
        if int(n_threads) > 1:
            caller = 'srun -n 1 -c %s --hint=multithread --cpu_bind=q %s' % (n_threads, caller)
        else:
            caller = 'srun -n 1 %s' % caller

        caller = shlex.split(caller)
        process = subprocess.Popen(caller, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode(locale.getdefaultlocale()[1])
        err = err.decode(locale.getdefaultlocale()[1])
        if (not ignore_errors) and (process.returncode):
            print("call failed, call was: %s" % ' '.join(caller))
            print("Message was: %s" % str(out))
            print("Error code was %s, stderr: %s" % (process.returncode, err))
        return process.returncode, out, err

    else:
        process = subprocess.call(caller, shell=True)
        return


cmdlist = []
filtered_cmdlist = []
filenames = []
new_filenames = []
filtered_new_filenames = []
download_list = []
merge_strings = []
names = []
files = []
name_files_dict = {}

sh = False

# Assigns each column a key, corresponding to the structure of "input.csv".
header_row = ["name", "File1", "URL1"]
with open('input.csv', 'rU') as csvfile:
    list = csv.DictReader(csvfile, header_row)
    for row in list:
        new_filenames.append(str(row["name"]) + "_" + str(row["File1"]))
        download_list.append([row["URL1"]])
        # Name and available files are stored in lists at same positions.
        files.append([row["File1"]])
        names.append(row["name"])
# Using the two lists "names" and "files" creates a dictionary with the name of e.g. Histone modification as a key
# and corresponding files from the table as values.
for r in range(len(names)):
    name_files_dict[names[r]] = files[r]
# Filter for empty elements and create a new filtered list
for x in new_filenames:
    if len(x) < 8:
        continue
    else:
        filtered_new_filenames.append(x)
tmp_list = []
while pipeline_counter <= 11:
    while(threading.active_count() > 1):
        continue
    if pipeline_counter in [1, 2, 3, 4, 5, 7]:
        for j in range(len(filtered_new_filenames)):
            #renames downloaded files
            if pipeline_counter == 1:
                for fil in files:
                    for x in fil:
                        for filt in filtered_new_filenames:
                            if (str(x) not in str(filt)):
                                continue
                            else:
                                cmdlist.append("mv " + folder + str(x) + " " + folder + str(filt))

            elif pipeline_counter == 2:
                file = path.exists(folder + filtered_new_filenames[j][:-4] + ".fastq")
                if (file):
                    continue
                #sratoolkit
                cmdlist.append("-o out.log --threads=" + cpt + fastqDump + folder + filtered_new_filenames[j])

            elif pipeline_counter == 3 and quality == True:
            #quality control for the fastq files
                cmdlist.append("-o out.log --threads=" + cpt + "./automated_fastqc2.sh")
                quality = False

            # If a file with new ending is needed, string manipulation is used. filtered_new_filenames[j][:-4] would be
            # "H3K27me3_SRR933992", filtered_new_filenames[j][:-14] would be "H3K27me3". Command in pipeline_counter = 4
            # merges data from different experiments but same modification, therefore from pipeline_counter = 5 onwards
            # specific experiment IDs aren't needed anymore.
            elif pipeline_counter == 4:
                file = path.exists(folder + filtered_new_filenames[j][:-4] + ".sam")
                if (file):
                    continue
                #bowtie2
                if(pipe):
                    cmdlist.append("-o out.log --threads=" + cpt + bowtie +
                                   "-U " + folder + filtered_new_filenames[j][:-4]
                                   + ".fastq -S | " + samtoolsView
                                   + ".sam -o " + folder + filtered_new_filenames[j][:-4] +
                                   ".bam ")
                else:
                    cmdlist.append("-o out.log --threads=" + cpt + bowtie +
                                   "-U " + folder + filtered_new_filenames[j][:-4]
                                   + ".fastq -S " + folder + filtered_new_filenames[j][:-4] + ".sam")

            elif pipeline_counter == 5 and not pipe:
                file = path.exists(folder + filtered_new_filenames[j][:-4] + ".bam")
                if (file):
                    continue
                if(delFastQ == True):
                    subprocess.call(["rm", folder + "*.fastq"])
                    delFastQ = False
                #samtools
                cmdlist.append("-o out.log --threads=" + cpt + samtoolsView
                               + folder + filtered_new_filenames[j][:-4] + ".sam -o " + folder + filtered_new_filenames[j][:-4] + ".bam ")

            elif pipeline_counter == 7:
            #samtoolsSort
                if(delSam == True):
                    subprocess.call(["rm", folder + "*.sam"])
                    delSam = False
                file = path.exists(folder
                                   + filtered_new_filenames[j][:-14] + "_sort.bam")
                if not file:
                    #samtools
                    cmdlist.append("-o out.log --threads=" + cpt + samtoolsSort + folder + filtered_new_filenames[j][:-14] + ".bam -o "
                                   + folder + filtered_new_filenames[j][:-14] + "_sort.bam " )
                else:
                    continue

    # Downloads all files from ftp URL's specified in input.csv
    elif pipeline_counter == 0:
        for file_URLs in download_list:
            for URL in file_URLs:
                cmdlist.append("wget " + str(URL))
        print(cmdlist)

    # Merges data from experiments with same genome modification
    elif pipeline_counter == 6:
        for genome_modification_name, files_to_merge in name_files_dict.items():
            # Filters all files for a specific genome modification. name_files_dict has empty values as of now.
            # If a column in input.csv is not empty, add it to tmp_list.
            for x in files_to_merge:
                if x != "":
                    tmp_list.append(x)
            # Reassigns name_files_dict[genome_modification_name] with updated list of files without empty elements.
            name_files_dict[genome_modification_name] = tmp_list
            # Empties tmp_list for next iteration
            tmp_list = []
        # Create executeable merge command
        for genome_modification in name_files_dict:
            # for each file corresponging to a genome modification
            for experiments in range(len(name_files_dict[genome_modification])):
                # create list of all files that need to be merged
                tmp_list.append(folder + genome_modification + "_" + name_files_dict[genome_modification][experiments][:-4]
                                + ".bam ")
                # create initial commandstring with merged output filename e.g. H3K27me3.bam
                cmd_string = "-o out.log --threads=" + cpt + samtoolsMerge + folder + genome_modification \
                             + ".bam "
                for x in range(len(name_files_dict[genome_modification])):
                    # Append all input files, who will be merged, to initial commandstring
                    cmd_string += folder + genome_modification + "_" + name_files_dict[genome_modification][x][:-4] \
                                  + ".bam "
                    # if all input files are appended, create an executeable commandstring to merge and add to cmdlist
                    if x == len(name_files_dict[genome_modification]) - 1:
                        cmdlist.append(cmd_string)
            del tmp_list[:]

    elif pipeline_counter == 8:
        for a in name_files_dict:
            # Checks if merged and indexed bam file is present. If not tries to index the merged and sorted
            # file e.g. H3K27me3.bam
            file = path.exists(folder + a + "_sort.bam")
            if not file:
                cmdlist.append("-o out.log --threads=" + cpt + samtoolsIndex + folder + a + "_sort.bam")
                # Adds a new row e.g. "H3K27" to indexed_filelist.txt (not a clever way via shell)
                # not needed, just a control
                subprocess.call(["echo " + folder + a + " >> indexed_filelist.txt"], shell=True)
            else:
                continue

    #peak calling, needs to be changed to calling bash script to work around version conflicts
    elif pipeline_counter == 9:
        if sh:
            continue
        else:
            sh = True
            bash = open("benderBash.sh", "w+")
            bash.write("#!/bin/bash\n"
                       "\n" +
                       macsModule)
            for a in name_files_dict:
                file = path.exists(folder + a + "_summits.bed")
                if file:
                    continue
                bash.write(macs2 + a + "_sort.bam -n " + a + "\n")
            subprocess.call(["chmod", "+x", "benderBash.sh"])
            cmdlist.append("-o out.log --threads=" + cpt + "./benderBash.sh")
            bash.close()

    elif pipeline_counter ==11 and enrichment == True:
    #calls the R script for the enrichment step
        inp = True
        bamlist = open("bamFiles.txt", "w+")
        bamlist.write("FileName\tSampleName\n")
        for x in names:
            if inp:
                bamlist.write(str(x) + "_sort.bam\tInput\n")
                inp = False
            bamlist.write(str(x) + "_sort.bam\t" + str(x) + "\n")
        bamlist.close()
        cmdlist.append("-o out.log --threads=" + cpt + "R CMD BATCH --args bamFiles.txt master.bed " + cpt + "outfile NormEnrichments_parallel.R")
        enrichment = False

    elif pipeline_counter == 10:
    # Creates a txt file with names of experiments already added to master.bed
        with open ("done_files.txt", 'r+') as donefiles:
            # stores names that are done
            done_files = [x for x in donefiles if x not in names]
        done_names = []
        for h in range(len(done_files)):
            done_names.append(done_files[h][:-1])
        unique_names = set(names)
        # removes names of done files from names to execute
        for s in done_names:
            unique_names.remove(s)
        print(done_names)
        print(len(unique_names))
        for x in unique_names:
            cmdlist.append([    # take first 5 colums with relevant data from narrowpeak file, save tmp_peaksfile
                                "cut -f 1-5 " + x + "_peaks.narrowPeak > tmp_peaksfile.bed",
                                # finds peaks that are already in master.bed with an overlap > 50 %
                                # saves in overlapping_peaks.bed
                                "-o overlapping_peaks.bed --threads=" + cpt + bedtools +
                                "intersectBed "
                                "-a tmp_peaksfile.bed -b master.bed -f 0.5 -r -wa",
                                # finds peaks from tmp_peaksfile.bed which are not in overlapping_peaks.bed so only new peaks are present
                                # in new_peaks.bed
                                "-o new_peaks.bed --threads=" + cpt + bedtools +
                                "subtractBed "
                                "-a tmp_peaksfile.bed -b overlapping_peaks.bed -f 1 -r",
                                # merges new peaks with overlap < 50 % with existing peaks to master.bed
                                "cat new_peaks.bed >> master.bed",
                                #sorting of the masterFile
                                "-o out.log --threads=" + cpt + bedtools +
                                "sortBed "
                                "-i master.bed",
                                # appends a new line in done_files.txt with modification name to see which files are done
                                "echo " + x + " >> done_files.txt"])
        # removes temporary files to reset pipeline
        remove_list = [folder + "/new_peaks.bed",
                       folder + "/overlapping_peaks.bed"]

    # filter cmdlist for empty false "wget " e.g. when no URL 5 present
    for x in cmdlist:
        if len(x) < 15:
            continue
        else:
            filtered_cmdlist.append(x)
    print(filtered_cmdlist)

    if filtered_cmdlist == []:
        filtered_cmdlist.append("echo your commandlist is empty")

    with open('cmdlist.csv', 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL, delimiter="\n")
        wr.writerow(filtered_cmdlist)

    print(len(filtered_cmdlist))

    # output_list creates a log file to see which commands have been executed.
    output_list = []

    # counts how many commands have been executed
    count = 0

    while (threading.activeCount() > 1):
        continue

    while True:

        # extra condition for downloading with input user delay
        if pipeline_counter == 0:

            for r in range(len(filtered_cmdlist)):
                caller = filtered_cmdlist[r]
                submit(caller)
                time.sleep(download_delay)



        else:

            print(count)

            # parallel command execution

            for l in range(0, len(filtered_cmdlist)):
                print(count)

                while (threading.active_count() == (tasks + 1)):
                    continue

                try:

                    if (count <= (len(filtered_cmdlist) -1)):
                        output = filtered_cmdlist[count]
                        t = threading.Thread(target=submit, name="Thread-" + str(count), args=(output,))
                        print("Thread-" + str(count))
                        t.start()

                        output_list.append(output)
                        output_list.append("cmd executed at: " + str(datetime.now()))
                        print("cmd: " + filtered_cmdlist[count] + " executed at: " + str(datetime.now()))
                        count += 1
                        print(count)

                finally:
                    print("inside finished ")
                continue

        # resets commandlists for next pipeline_counter
        cmdlist = []
        filtered_cmdlist = []

        pipeline_counter += 1

        break
