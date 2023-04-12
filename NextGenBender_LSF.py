import csv
import time
import subprocess
from datetime import datetime
from os import path
import os

'''
The input file has to be in the CSV format and named "input.csv". The file has to be made of 3 columns: column 1 is the name of the feature, column 2 the filename of the corresponding sequencing file (SRA or Fastq format), and column 3 contains an optional URL to download the file 
feature     filename    URL
example input:
EZH2,SRR568340.sra,https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR568/SRR568340/SRR568340.sra

To start the pipeline from any point, change the pipeline counter value (0: full pipeline with download,
 2: without download, 4: without fastqDump (if you have sam files), 5: samtools View, 6: merging for files of the same feature, 7: samtoolsSort, 8: samtoolsIndex,
 9: peak calling with MACS2, 10: creates the masterFile, 11: only read conting)

Personalized settings
after each value in " " hast to be a space. npt, cpt, download_delay, and pipeline_counter have to be integer numbers.
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
reference = "/dir/to/ref/genome/hg19 "                                              #path to the reference genome for bowtie
folder = ""                                                                         #working directory, leave empty if the same as script folder
fastqDump = "/dir/tosratoolkit/bin/fastq-dump "                                     #path to fastqDump
bowtie = "bowtie2 -p " + cpt + " -x " + reference                                   #bowtie2 settings
samtoolsView = "samtools view -b "                                                  #samtools settings
samtoolsSort = "samtools sort "
samtoolsMerge = "samtools merge "
samtoolsIndex = "samtools index "
macs2 = "macs2 callpeak -t "                                                        #macs parameters
macsModule = "module load bio/MACS2/2.1.0.20150731-foss-2017a-Python-2.7.13\n"      #odule required for MACS2
bedtools = "/dir/to/bedtools2/bin/"                                                 #path to bedtools


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

tmp_master_bed = []


# Assigns each column a key, corresponding to the structure of "ChipSeqFilelist.csv".
header_row = ["name", "File1", "URL1"]



with open('ChipSeq_Filelist_3.csv', 'rU') as csvfile:
    list = csv.DictReader(csvfile, header_row)
    for row in list:
        new_filenames.append(row["name"] + "_" + row["File1"])

        download_list.append(row["URL1"])

        # Name and available files are stored in lists at same positions.
        files.append(row["File1"])
        names.append(row["name"])

f = open('done_files.txt')
f.close()

# Using the two lists "names" and "files" creates a dictionary with the name of e.g. Histone modification as a key
# and corresponding files from the table as values.
for r in range(len(names)):
    name_files_dict[names[r]] = files[r]


# Filter for empty elements and create a new filtered list
for x in new_filenames:
    if len(x) < 15:
        continue
    else:
        filtered_new_filenames.append(x)


pipeline_counter = 0
tmp_list = []

# specify amound of parallel commands
parallel_cmd_amount = int(input("# of simultaneous commands: "))

while pipeline_counter <= 11:
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

            # Creates an executeable shell command for each .sra file using filtered_new_filenames.
            # filtered_new_filenames example element: H3K27me3_SRR933992.sra

            if pipeline_counter == 2:
                cmdlist.append("bsub -G jgu-adprp  -q short -app Reserve1800M "
                               "-o out.log" + fastqDump
                               + filtered_new_filenames[j])

            elif pipeline_counter == 3 and quality == True:
            #quality control for the fastq files
                cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M -n 20 "
                               "-o out.log ./automated_fastqc2.sh")
                quality = False
            # If a file with new ending is needed, string manipulation is used. filtered_new_filenames[j][:-4] would be
            # "H3K27me3_SRR933992", filtered_new_filenames[j][:-14] would be "H3K27me3". Command in pipeline_counter = 4
            # merges data from different experiments but same modification, therefore from pipeline_counter = 5 onwards
            # specific experiment IDs aren't needed anymore.
            elif pipeline_counter == 4:
                cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M -n 20 "
                               "-o out.log ")
                if(pipe):
                    cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M -n 20 "
                                   "-o out.log "+ bowtie +
                                   "-U " + folder + filtered_new_filenames[j][:-4]
                                   + ".fastq -S | " + samtoolsView
                                   + ".sam -o " + folder + filtered_new_filenames[j][:-4] +
                                   ".bam ")
                else:
                    cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M -n 20 "
                                   "-o out.log " + bowtie +
                                   "-U " + folder + filtered_new_filenames[j][:-4]
                                   + ".fastq -S " + folder + filtered_new_filenames[j][:-4] + ".sam")

            elif pipeline_counter == 5 and not pipe:
                if(delFastQ == True):
                    subprocess.call(["rm", folder + "*.fastq"])
                    delFastQ = False
                cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M "
                               "-o out.log " + samtoolsView
                               + filtered_new_filenames[j][:-4] + ".bam " + filtered_new_filenames[j][:-4] + ".sam")

            elif pipeline_counter == 7:
                if(delSam == True):
                    subprocess.call(["rm", folder + "*.sam"])
                    delSam = False
                file = path.exists(folder
                                   + filtered_new_filenames[j][:-14] + "_sort.bam")

                if not file:
                    cmdlist.append("bsub -G jgu-adprp -W 120 -q short "
                                   "-app Reserve1800M " + samtoolsSort +
                                   "-m 1000000000 " + filtered_new_filenames[j][:-14] + ".bam -o "
                                   + filtered_new_filenames[j][:-14] + "_sort.bam")

                else:
                    continue




    # Downloads all files from ftp URL's specified in ChipSeqFilelist.csv
    elif pipeline_counter == 0:

        for file_URLs in download_list:
            for URL in file_URLs:
                cmdlist.append("wget " + URL)

    # Merges data from experiments with same genome modification
    elif pipeline_counter == 6:

        for genome_modification_name, files_to_merge in name_files_dict.items():

            # Filters all files for a specific genome modification. name_files_dict has empty values as of now.
            # If a column in ChipSeqFilelist.csv is not empty, add it to tmp_list.
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
                tmp_list.append(genome_modification + "_" + name_files_dict[genome_modification][experiments][:-4]
                                + "_sort.bam ")

                # create initial commandstring with merged output filename e.g. H3K27me3.bam
                cmd_string = "bsub -G jgu-adprp -q short -W 120 -app Reserve1800M " \
                             "-o out.log " + samtoolsMerge + genome_modification \
                             + ".bam "

                for x in range(len(name_files_dict[genome_modification])):

                    # Append all input files, who will be merged, to initial commandstring
                    cmd_string += genome_modification + "_" + name_files_dict[genome_modification][x][:-4] \
                                  + "_sort.bam "

                    # if all input files are appended, create an executeable commandstring to merge and add to cmdlist
                    if x == len(name_files_dict[genome_modification]) - 1:
                        cmdlist.append(cmd_string)
            del tmp_list[:]



    elif pipeline_counter == 8:
        for a in name_files_dict:

            # Checks if merged and indexed bam file is present. If not tries to index the merged and sorted
            # file e.g. H3K27me3.bam

            file = path.exists("~/" + a + "_sort.bam")
            if not file:
                cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M "
                               "-o out.log " + samtoolsIndex + a + "_sort.bam")

                # Adds a new row e.g. "H3K27" to indexed_filelist.txt (not a clever way via shell)
                # not needed, just a control
                subprocess.run(["echo " + a + " >> indexed_filelist.txt"], shell=True)

            else:
                continue



    #using MACS2 for the call peaks
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
            subprocess.call(["chmod", "+x", "benderBash"])
            cmdlist.append("bsub -G jgu-adprp -W 120 -q short "
                           "-app Reserve1800M ./benderBash.sh")
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
        cmdlist.append("bsub -G jgu-adprp -W 120 -q short -app Reserve1800M "
                       "-o out.log R CMD BATCH --args bamFiles.txt master.bed " + cpt + "outfile NormEnrichments_parallel.R")
        enrichment = False

    #creating the master.bed
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
                                "bsub -G jgu-adprp -q short -W 120 -o subtract_from_query_tmp.bed "
                                "-app Reserve10G " + bedtools +
                                "intersectBed "
                                "-a tmp_peaksfile.bed -b master.bed -f 0.5 -r -wa",
                                # finds peaks from tmp_peaksfile.bed which are not in overlapping_peaks.bed so only new peaks are present
                                # in new_peaks.bed
                                "bsub -G jgu-adprp -q short -W 120 -o add_to_master_tmp.bed "
                                "-app Reserve10G " + bedtools +
                                "subtractBed "
                                "-a tmp_peaksfile.bed -b overlapping_peaks.bed -f 1 -r",
                                # merges new peaks with overlap < 50 % with existing peaks to master.bed
                                "cat new_peaks.bed >> master.bed",
                                #sorting of the masterFile
                                "bsub -G jgu-adprp -q short -W 120 -o out6.log "
                                "-app Reserve10G " + bedtools +
                                "sortBed "
                                "-i master.bed",
                                # appends a new line in done_files.txt with modification name to see which files are done
                                "echo " + x + " >> done_files.txt"])
        # removes temporary files to reset pipeline
        remove_list = [folder + "/new_peaks.bed",
                       folder + "/overlapping_peaks.bed"]


        # cmdlist is a list of lists, iterates through cmdlist and executes commands in cmdlist[j], when intersecting or
        # subtracting delays by longer period than if a simpy append command is executed. delays increase each iteration
        # because filesizes increase too. prints last 50 lines of master.bed after each iteration to check for errors.
        # you should see only bed-type lines, no error message. lastly prints number of lines in master.bed
        for j in range(len(cmdlist)):
            for k in range(len(cmdlist[j])):
                subprocess.call([cmdlist[j][k]], shell=True)
                print(cmdlist[j][k])

                if "bsub" in cmdlist[j][k]:
                    print("long wait")
                    time.sleep(j + 90)
                else:
                    print("short wait")
                    time.sleep(30)
            print("Durchlauf Nr: " + str(j))

            subprocess.call(["tail -50 master.bed"], shell=True)
            time.sleep(20)
            subprocess.call(["wc -l master.bed"], shell=True)

            try:
                for e in remove_list:
                    remove(e)
                    print(e)
                    time.sleep(10)
            except FileNotFoundError:
                pass


    # filter cmdlist for empty false "wget " e.g. when no URL 5 present
    for x in set(cmdlist):
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


    while True:

        # extra condition for downloading with input user delay
        if pipeline_counter == 0:

            for r in range(len(filtered_cmdlist)):
                subprocess.run([filtered_cmdlist[r]], shell=True)
                time.sleep(download_delay)


        else:
            print(count)

            # parallel command execution
            for l in range(0, len(filtered_cmdlist)//parallel_cmd_amount):
                # bjobs lists all jobs currently active
                check_jobs = subprocess.run(["bjobs"], stdout=subprocess.PIPE, shell=True)
                cmdoutput = check_jobs.stdout.decode('utf-8')
                print(len(cmdoutput))
                print(cmdoutput)
                time.sleep(0.5)
                print(count)

                # If bjobs returns "no jobs currently active" begin executing parallel_cmd_amount commands from
                # filterd_cmdlist with a delay of 2 seconds inbetween. Also generates output to inform user.
                # Increments count by 1 each time a command is executed
                if len(cmdoutput) < 30:
                    try:
                        for k in range(0, parallel_cmd_amount):
                            print("inside parallel")
                            output = subprocess.run([filtered_cmdlist[count]], shell=True)
                            output_list.append(output)
                            output_list.append("cmd executed at: " + str(datetime.now()))
                            print("cmd: " + filtered_cmdlist[count] + " executed at: " + str(datetime.now()))
                            count += 1
                            print(count)
                            time.sleep(0.5)
                    finally:
                        print("inside finished ")
                    continue
                else:
                    time.sleep(120)
                    print("wait here")

            # If count reaches a point where count + parallel_cmd_amount > current command list length,
            # executes the remaining cmds
            if count == len(filtered_cmdlist) - (len(filtered_cmdlist) % parallel_cmd_amount):
                for m in range(len(filtered_cmdlist) % parallel_cmd_amount):
                    try:
                        output = subprocess.run(filtered_cmdlist[count], shell=True)
                        output_list.append(output)
                        output_list.append("cmd executed at: " + str(datetime.now()))
                        print("cmd: " + filtered_cmdlist[count] + " executed at: " + str(datetime.now()))
                        count += 1
                    except FileNotFoundError:
                        print(filtered_cmdlist[count] + " could not be executed")

            # Checks "bjobs" for running jobs. If bjobs return a smaller string than 30 ("no jobs running"), breaks
            while True:
                check_jobs = subprocess.run(["bjobs"], stdout=subprocess.PIPE, shell=True)
                cmdoutput = check_jobs.stdout.decode('utf-8')
                print(len(cmdoutput))
                print(cmdoutput)

                if len(cmdoutput) < 30:
                    break
                else:
                    time.sleep(60)


        # resets commandlists for next pipeline_counter
        cmdlist = []
        filtered_cmdlist = []

        pipeline_counter += 1

        break
