#!/gsc/biosw/bin/python3.7

##############
#
# This script simulates reads from reference sequences which are stored in a
# fasta file. It stores the reads in a fasta or fastq file. To simulate a
# fastq file a fastq file of a real experiment is needed. It extracts the
# quality strings. A quality string is randomly selected for each sequence.
# Each base of the read-sequence can be changed on chance which depends on the
# qualtiy value at the same position in the quality string. After the
# simulation process the quality string and the mutated sequence are merged and
# stored in a fastq-file.
# I used DNemulator as template <http://cbrc3.cbrc.jp/~martin/dnemulator/>. 
# 
# It is possible to set a seed to get the same random file.
#
# Usage: simulator.py [flag][value]
#
# options:
#   -f fasta-file
#   -fq fastq-file
#   -l read length, the reference quality strings need to be as long the reads
#   you want to simulate
#   -d additional set with differentially expressed elements
#   -o fasta
#   -sz number replicates 
#   -r library size, the final files will not contain the exactly number of
#   reads
#
# author: Robert Schwarz
# date: 23.04.19
#
#######


########
# To-DO
#  - make the read length variable in the single end simulation
#  - add a checkpoint where the length of the reads and the qual strings is
#    compared
#  - write a log file with information like: processed sequences, time,
#  mutated bases and so on
#  - add the possibility to simulate a fasta file --> this does not need the mutation process for instance
#  - make the ref seq length restriction variable and set the default to 100 bp
#  - think about threads
#########

##########
# import #
##########

import time
import seq_handling 
import simulation as sim
import argparse
import random
import sys
import re
import dataSet
import yaml
import numpy as np

def parse_options():

    parser=argparse.ArgumentParser(prog="Fastq simulator", description="This\
            script simulates a fastq file. Therfore it needs a fasta file \
                    which contains the reference sequences and a fastq file\
                    from a real experiment.")
   
    parser.add_argument("-f", type=argparse.FileType('r'), dest="fastaFile")
    parser.add_argument("-t", type=int, default = 100, dest="min_length",
            help="set the min length of the instances that are simulated \
            (default = 100)")
    parser.add_argument("-fq", type=argparse.FileType('r'), dest="fastqFile")
    parser.add_argument("-l", type=int, default = 100, dest="readLength")
    parser.add_argument("-o", type=str, default="simulated_fastqs", dest="output")
    parser.add_argument("-w", action='store_true', dest="warn",
            help="set on warnings")
    parser.add_argument("-s", type=int, default = -1, dest="seed")
    parser.add_argument("-r", type=int, default=1000000, dest="librarySize")
    parser.add_argument("-sz", type=int, default = 1, dest="sampleSize")
    parser.add_argument("-d", action="store_true", dest="diff")
    parser.add_argument("-p", action="store_true", dest="paired")
    args = parser.parse_args()
   
    return parser, args
          

def create_time_stamp():

    stamp = time.strftime("%Y%m%d_%H%M%S", time.gmtime())

    return stamp

def fastq_simulator(args, seed):

    print("Initialize the data set...")
    
    references, windows = seq_handling.load_reference_seq(args.fastaFile.name,
            diff = args.diff)
    
    data = dataSet.dataSet(seed,
            args.fastqFile.name,
            references,
            windows,
            args.sampleSize,
            args.readLength,
            args.min_length,
            args.librarySize,
            args.diff,
            args.paired,
            args.output)
            
def set_seed(seed):

    if seed == -1:

        seed = random.choice(list(range(1,1000)))

    random.seed(seed) 
    np.random.seed(seed)

    
def main():

    parser, args = parse_options()
    
    time_stamp = create_time_stamp()
    
    seed = set_seed(args.seed)

    fastq_simulator(args, seed)
   

if __name__ == "__main__":
    
    main()
