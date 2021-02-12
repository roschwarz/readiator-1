import sys
import pandas as pd
import numpy as np
import re
import seq_handling
from Bio.SeqRecord import SeqRecord
import random
from collections import Counter
from time import sleep

class simRecord(SeqRecord):

    def __init__(self, readSeq, readID, qualStr, R = 'R1'):

        SeqRecord.__init__(self, readSeq, readID)
        self.qualStr = qualStr
        self.R = R
         
        if self.R == 'R2':
           
           self.seq = self.seq.reverse_complement()

        
        self.sim_fastq_entry()

    def sim_fastq_entry(self):
        
        self.fastqSeq, self.mutatedBases = seq_handling.mutate_read(self.seq,
                self.qualStr)


def single_end_simulation(references, qualityStrings, counts, readLength):
    

    sample = counts.columns[0]
    print('Simulate {}'.format(sample))
    reads = []

    for element in counts.itertuples(index = True, name = 'Pandas'):

        record = references[element[0]]
        count = 0 
        readCounter = 0

        for read in seq_handling.draw_sub_sequences(record.seq,
                record.reads[sample],
                readLength):

            count += 1
            readCounter += 1

            if random.random() < 0.5:

                read = seq_handling.rev_comp(str(read))

                readID = '{}_{}_rev_{}'.format(str(record.id[0:50]),
                        str(count), sample)
            else: 

                readID = '{}_{}_fwd_{}'.format(str(record.id[0:50]),
                        str(count), sample)
            
            r = simRecord(read, str(readID), random.choice(qualityStrings)) 
            
        
            reads.append(r)

        assert readCounter == record.reads[sample], 'simulated number unequal \
        requested number'

    return reads

def extract_paired_reads(fragment, readID, qualityStrings, readLength):
  
    if len(fragment) < readLength:

        R1_qual = random.choice(qualityStrings)[0:len(fragment)]
        R2_qual = random.choice(qualityStrings)[0:len(fragment)]

        R1_seq = fragment
        R2_seq = fragment
    
        assert len(R1_seq) == len(fragment), "fragment too short, read length \
should have as long as the fragment"
    
        assert len(R1_qual) == len(fragment), "fragment too short, qual str \
length should have as long as the fragment"
        
        assert len(R1_seq) != 1, "read length is 1"

    else:

        R1_qual = random.choice(qualityStrings)
        R2_qual = random.choice(qualityStrings)

        R1_seq = fragment[0:readLength]
        R2_seq = fragment[len(fragment)-readLength:len(fragment)]

    if len(R1_qual) == len(R1_seq) and len(R2_qual) == len(R2_seq):
        
        R1 = simRecord(R1_seq, readID, R1_qual, "R1")
        R2 = simRecord(R2_seq, readID, R2_qual, "R2")
    
        return (R1, R2)

    else:
        print("Fragment length: {}".format(len(fragment)))
        sys.stderr.write("Read {} has not the same length as quality \
string\n".format(readID))
        return None

def paired_end_simulation(references, qualityStrings, counts, readLength):

    sample = counts.columns[0]
    print('Simulate {}'.format(sample))
    reads = {'R1':[], 'R2':[]}
    fragmentLength = []

    for element in counts.itertuples(index = True, name = 'Pandas'):
        
        record = references[element[0]]
        count = 0 
        fragCounter = 0    

        if len(record.seq) >= readLength: 
       
            for fragment in seq_handling.draw_fragments(record.seq, 
                record.reads[sample], readLength):
                
                fragCounter +=1         
                
                fragmentLength.append(len(fragment))

                count += 1
                readID = '{}_{}_{}'.format(str(record.id[0:50]), 
                        str(count),
                        sample)

                read_pair = extract_paired_reads(fragment, 
                        readID, 
                        qualityStrings,
                        readLength)
                
                assert len(read_pair[0]) == len(read_pair[1]), 'read pairs does\
                not have the same length'

                if read_pair != None:
                
                    reads['R1'].append(read_pair[0])
                    reads['R2'].append(read_pair[1])

        assert fragCounter == record.reads[sample], "The true fragment count \
                and given fragment count is different"

    return reads, Counter(fragmentLength)

def set_mean_expression(record, diff, librarySize):

    if diff and record.diff:

        if record.up:
                record.meanExpression = int(record.cpkm * librarySize) 
                record.meanExpressionDiff = round(record.meanExpression * \
                record.diffFactor)
        else:

                record.meanExpression = int(record.cpkm * librarySize) 
                record.meanExpressionDiff = round(record.meanExpression / \
                        record.diffFactor)
    else:

        record.meanExpression = int(record.cpkm * librarySize)
        record.meanExpressionDiff = record.meanExpression


def set_sample_counts(record, sample):

    if record.diff and re.search(r'diff', sample):
        
        record.reads[sample] = get_count_NB(record.meanExpressionDiff)[0]

    else:

        record.reads[sample] = get_count_NB(record.meanExpression)[0]

def get_count_NB(meanExpression, requestedValues = 1, rsize = 10):
    
    if meanExpression != 0:
    
        probability = rsize/(rsize+meanExpression)
        mu = (meanExpression*probability)/(1-probability)
    
        return np.random.negative_binomial(mu, probability, requestedValues)
    
    else:

        return np.zeros(requestedValues, dtype=int)


def is_diff_expressed(record, fraction = 0.05, diffFactorCap = 10, fracUpDown = 0.5):

    if random.random() <= fraction:

        record.diff = True
        record.diffFactor = random.randrange(2,diffFactorCap)
        
        if random.random() > fracUpDown:
            record.up = True
        else:
            record.up = False

    else:

        record.diff = False
        record.diffFactor = 1
        record.up = False


    return record
