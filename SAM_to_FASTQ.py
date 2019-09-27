#!/usr/bin/env python
# coding: utf-8



'''
put all of your reference seqs in ONE fasta, we will handle the rest
STEPS:
    1. read in and parse reference fasta 
    2. simulate reads from fasta using ART ILLUMINA 
    3. read in the SAM files created and simulate AR
    4. output a FASTQ

'''
# pysam docs
# https://pysam.readthedocs.io/en/latest/index.html

import os
import argparse
from Bio import SeqIO
import numpy as np
import pysam
from collections import Counter
import csv


# ## Need to set some variables here:
#     -- FASTA Input path 
#     -- SAM files
#     -- freqs
#     -- ref names
#     -- % AR to simulate
#####################################


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="path to your reference fasta files", type=str)
parser.add_argument("-fr", "--freqs", help="list of frequencies for each virus in sample", type=str)
parser.add_argument("-ar", "--artrecomb", help="the PERCENT of artificial recombination you want in your fastQ file", type=int)

args = parser.parse_args()

ref_fasta = args.file
data = SeqIO.parse(ref_fasta, "fasta")

fasta_files = []
num_refs = 0
for record in data:
    fasta_files.append(record.id)
    num_refs += 1
    out_fas = record.id + ".fasta"
    with open(out_fas, "w") as out:
        out.write(">{}\n{}\n".format(record.id, str(record.seq)))

freqs = [float(args.freqs.split(",")[i]) for i in range(num_refs)]
per_ar = float((args.artrecomb)/100)

print("Your reference fasta: ",ref_fasta, "\nThe number of references in the file: ", num_refs, "\nThe frequencies of each of your viruses: ", freqs, "\nThe % of AR that you want to simulate: ",per_ar)


## make this dicitonary here where we have {virus_names: freq}
reads_to_sim = 10000
for i in range(num_refs):
     
    virus_name = fasta_files[i]
     
    my_cmd_1 = "echo simulating from %s" % fasta_files[i]
    my_cmd_2 = "../art_src_MountRainier_MacOS/art_illumina -ss HS25 --samout -i %s -l 120 -s 50 -c %d -o %s" % (fasta_files[i] +".fasta", reads_to_sim, virus_name)
    
    #os.system(my_cmd_2)
    
    ## change this so that if it is made you can remake it, but leave it for now ##
    if os.path.exists(fasta_files[i] + ".sam") == True:
        print("%s already made" % (fasta_files[i] + ".sam"))
    else:
        os.system(my_cmd_1)
        os.system(my_cmd_2)

## gathering number of reads to find ##
over_estimate = int(reads_to_sim * .50)
desired_AR_reads = int(reads_to_sim * per_ar)
print("these are desired reads --> ", desired_AR_reads)

## this will put all of our sam file reads into an N (number of sam files) list of list ##
data = {}
for i in range(num_refs):
    d = {}
    temp = pysam.AlignmentFile(fasta_files[i]+".sam", "r")
    reads = list(temp.fetch())
    d[fasta_files[i]] = (freqs[i], reads)
    data.update(d)

## try to grab it all now... want {(read_a, read_b, number): index} #
sampled_AR_data = []
counter = 0
for number in range(desired_AR_reads): 
    #print(number, counter)
    viruses = list(data.keys())

    seq_grabber = np.random.randint(reads_to_sim)
    i = np.random.choice(num_refs, 1, True, np.array(freqs))[0]
    #print("first read comes from --> ", fasta_files[i])
    
    ## this is getting read A to be the start of the AR ##
    read_A = data[fasta_files[i]][1][seq_grabber]
    
    read_A_coords = (read_A.reference_start, read_A.reference_end)
    #print(read_A_coords)

    #print(viruses[i])
    viruses.pop(i)
    #print("popped ",viruses)
    
    freqs = np.array(freqs)
    freq_reweight_list = freqs[np.arange(len(freqs)) !=i]
    freq_reweight = freq_reweight_list / np.sum(freq_reweight_list)
    j = np.random.choice(len(freq_reweight), 1, True, freq_reweight)[0]
    #print(j)
    #print(viruses[j]) 
    new_key = viruses[j]

    read_B = data[new_key][1][seq_grabber]
    #print(read_B)
    read_B_coords = (read_B.reference_start, read_B.reference_end)
    
    #print(read_B_coords)
    
    if read_A_coords == read_B_coords:
        sampled_AR_data.append((read_A, read_B))
        counter += 1
    
    else:
        continue

    if counter == per_ar:
        print("Gathered all our ARs!!")
        break
    else:
        continue


with open("tester.fastq", "w") as out:
    for i in range(reads_to_sim):
        if i < len(sampled_AR_data):
            # write AR reads to file
            #count += 1
        
            ## randomly grab a seq ## 
            seq_getter = np.random.randint(len(sampled_AR_data))
            r_a, r_b = sampled_AR_data[i][0], sampled_AR_data[i][1]

            '''
            here I need to:
                1. index one seq from the other
                2. combine the seqs
                    a. label
                    b. joined seq
                    c. joined q-score
            '''
            seq_splitter = np.random.randint(len(r_a.seq))
            #print(seq_splitter)
            r_start = r_a.seq[:seq_splitter]
            r_end = r_b.seq[seq_splitter:]

            joined_label = r_a.qname + "+" + r_b.qname
            joined_seq = r_start+r_end
            joined_q = r_a.qual[:seq_splitter] + r_b.qual[seq_splitter:]

            out.write("{}\n{}\n+\n{}\n".format(joined_label, joined_seq, joined_q)) 
            
        else:
            ## jsut write everything else here
            viruses = list(data.keys())
            seq_getter = np.random.randint(reads_to_sim)
            key_get = np.random.randint(num_refs)
            
            r_a = data[viruses[key_get]][1][seq_getter]
            key_get = np.random.randint(num_refs)
            r_b = data[viruses[key_get]][1][seq_getter]
            
            out.write("{}\n{}\n+\n{}\n".format(r_a.qname, r_a.seq, r_a.qual)) 
            out.write("{}\n{}\n+\n{}\n".format(r_b.qname, r_b.seq, r_b.qual)) 

        


