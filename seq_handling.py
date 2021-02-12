#######################
#
# author: Robert Schwarz
# mail: rschwarz@leibniz-fli.de
# initial date:  Fri Sep 27 14:40:18 CEST 2019
#
# description: ....
#

##########
# import #
##########

import sys
import random
import string
import re
import math
import numpy as np
import simulation as sim
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter

###########
# methods #
###########

###################
# loading methods #
###################

# this parts handles with sequencing files

# extracts only the quality strings from a fastq file and returns them in a
# list
def read_fastq_qual(fastq_file):

    print("read quality strings...")

    quality_strings = []
    
    for name, seq, qual in FastqGeneralIterator(open(fastq_file)):
            
            quality_strings.append(qual)

    print("read quality strings done...")

    return quality_strings

   
#transform the phred_score into p_error value
def p_error(phred_score):

    return 10 **(-(ord(phred_score)-33)/10)

# a univeral method to load either a fasta or a fastq file
# Biopython is necessary to load the respective file. The method stores the
# records objects into a list and returns the list.
def load_seq_file(seq_file):

    print("load %s" %seq_file) 
    records = [] 


    if seq_file.split(".")[-1] in ["fa", "fasta"]:
        
        for rec in SeqIO.parse(seq_file, "fasta"):
            records.append(rec)

    elif seq_file.split(".")[-1] in ["fq", "fastq"]:

        for rec in SeqIO.parse(seq_file, "fastq"):
            records.append(rec)

    else:
        print("no seq file")

    return records

# This method reads a sequencing file and works with the yield function. In
# that way the sequence is gone after processing. This method safes memory.
def read_seq_file(seq_file):

    if seq_file.split(".")[-1] in ["fa", "fasta"]:
        
        for rec in SeqIO.parse(seq_file, "fasta"):
            yield rec

    elif seq_file.split(".")[-1] in ["fq", "fastq"]:

        for rec in SeqIO.parse(seq_file, "fastq"):
            yield rec

    else:
        print("no seq file")

#######
# load the reference sequences for the simulation, and filter for
# a certain length and determine how much windows of a certain size
# fits to the sequence length. Random assignment of differentially expression
# factors for simulation of a data set with differentially expressed instances.
#
###

def load_reference_seq(refFile, windowSize = 100, minLength = 100, 
        diff = False, diffFraction = 0.05, diffFactorCap = 10, fracUpDown = 0.5):
   
    records = {}
    windowsAll = 0

    for rec in read_seq_file(refFile):
        
        if rec.__len__() >= minLength:
            
            rec.windows = math.ceil(rec.__len__()/windowSize)
            
            if diff:

                rec = sim.is_diff_expressed(rec, 
                                            diffFraction, 
                                            diffFactorCap,
                                            fracUpDown)
                
            else:

                rec.diff = False
                rec.diffFactor = 1
                rec.up = False
            
            rec.reads = {}
            records[rec.id] = rec
            windowsAll += rec.windows


    return records, windowsAll


###################
# writing methods #
###################

def write_fastq(list_of_records, output="fastq_out", file_name ="sim.fq"):

    with open("{}/{}.fastq".format(output, file_name), "w") as output_handle:
       
        for rec in list_of_records:
            
            output_handle.write("@%s\n" %rec.id)
            output_handle.write("%s\n" %rec.fastqSeq)
            output_handle.write("+\n")
            output_handle.write("%s\n" %rec.qualStr)

def write_paired_fastqs(list_of_records, output="fastq_out", file_prefix = "sample"):

    with open("{}/{}_R1.fastq".format(output, file_prefix), "w") as R1_handle,\
         open("{}/{}_R2.fastq".format(output, file_prefix), "w") as R2_handle:

            for rec in list_of_records:

                if rec == 'R1':
                    for read in list_of_records[rec]:
                        R1_handle.write("@%s\n" %read.id)
                        R1_handle.write("%s\n" %read.fastqSeq)
                        R1_handle.write("+\n")
                        R1_handle.write("%s\n" %read.qualStr)

                else:
                    for read in list_of_records[rec]:
                        R2_handle.write("@%s\n" %read.id)
                        R2_handle.write("%s\n" %read.fastqSeq)
                        R2_handle.write("+\n")
                        R2_handle.write("%s\n" %read.qualStr)

#########################
# sequence manipulation #
#########################

# The base N can mutate to ATGC or remains N
def mutate_read(read, random_qual_string):

    mutated_read = ''
    mutated_bases = 0

    for idx, base in enumerate(read):

        if base in "ACGT" and random.random() < p_error(random_qual_string[idx]):

            mutated_read = mutated_read + random.choice("ATGC".replace(base, ""))
            mutated_bases += 1
        
        elif base == "N":

            mutated_read = mutated_read + random.choice("ATGCN")
            mutated_bases += 1

        else:

            mutated_read = mutated_read + base

    assert mutated_bases == sum(1 for a, b in zip(read, 
            mutated_read) if a != b), 'different number of mutations than expected'

    return mutated_read, mutated_bases


def sub_sequence(seq, start, length):
    
    if start + length > len(seq):
        err = ValueError("sub sequence exceeds end of the sequence")
        raise err
    else:
        return seq[start:start+length]

# This method creates random reads of an input sequence, 
# a fashion normal distribution.
def draw_sub_sequences(seq, read_count, read_length):
   
    sub_seq_list = []
    possible_start_pos = len(seq)-read_length
    
    for read in range(0, read_count):
        # check To-Do for random choice of a number from a list
        read = sub_sequence(seq, \
                random.choice(list(range(0,possible_start_pos))),\
                read_length)
        
        assert len(read) == read_length, "the read hast not the required \
        length"
        
        sub_seq_list.append(read)

    return sub_seq_list


def draw_fragments(seq, count, readLength,  mean = 200, sd = 100):
    "Draw given number of sub sequences from a reference sequence \
     the fragment length has to be at least as long as 25 % of the \
     read length"

    fragments = []
    
    while len(fragments) < count:

        fragmentLength = int(np.random.normal(mean, sd))

        if fragmentLength >= len(seq):
            fragment = seq
        elif fragmentLength < 0.25 * readLength:
            fragment = None
        else:
            fragmentStart = random.choice(list(range(0,len(seq) -  
                fragmentLength)))

            fragment = sub_sequence(seq, fragmentStart, fragmentLength)

        if fragment != None:
            fragments.append(fragment)

    return fragments

# To-Do: avoid the crash of the script by the raised error
def rev_comp(seq):
    
    if is_it_DNA(seq):
        return Seq(seq).reverse_complement()
    else:
        err = ValueError("seq is not a dna string")
        raise err

def rev_comp_chance(read_rec, probability):

    if random.random() >= probability:
        read_rec.seq = rev_comp(str(read_rec.seq))

######################
# checking sequences #
######################

# check if the seq is a DNA-String
def is_it_DNA(seq):
    
    nucleotides = ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']

    anti_nucleotides = list(string.ascii_uppercase)
    
    for an in list(string.ascii_lowercase):
       anti_nucleotides.append(an)
    for n in nucleotides:
       anti_nucleotides.remove(n)

    res = any(ele in seq for ele in anti_nucleotides)
    if res:
        return False
    else:
        return True

def n_check(seq):

    if re.search(r'N', str(seq)) or re.search(r'n', str(seq)):
        return True
    else:
        return False

def base_composition(seq):

    base_compo = Counter(seq)

    return base_compo
