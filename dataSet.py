
import seq_handling as seq_handling
import simulation as sim
import numpy as np
import pandas as pd
import math
import subprocess
import shlex
from collections import Counter

class dataSet():

    def __init__(self,
            seed,
            refFastq,
            references,
            components,
            sampleSize = 1,
            readLength = 50,
            minRefLength = 100,
            librarySize = 1000000,
            diff = False,
            paired = False,
            outDir = "simulated_data_set"):

        self.seed = seed
        self.refFastq = refFastq
        self.references = references
        self.sampleSize = sampleSize
        self.readLength = readLength
        self.librarySize = librarySize
        self.diff = diff
        self.paired = paired
        self.outDir = create_folder(outDir)
        self.minRefLength = minRefLength
        self.components = components
        self.fragmentLength = {}

        self.gen_sample_list()
        self.set_cpkms()
        self.set_counts()
        self.gen_count_table()
        self.load_quality_strings()
        self.simulate_fastq()
        self.write_count_table()
        
    def gen_sample_list(self):
        
        self.samples = []

        if self.diff:
            
            counter = 0

            for i in range(1,2*self.sampleSize+1):
                

                if i <= self.sampleSize:
                    self.samples.append("sample_{}".format(i))
                else:
                    counter += 1
                    self.samples.append("sample_diff_{}".format(counter))

        else:

            for i in range(1,self.sampleSize+1):
                
                self.samples.append("sample_{}".format(i))

    def set_cpkms(self):

        diriList = dirichletList(self.components) 
        array_start = 0
        
        for record in self.references:
            
            record = self.references[record]
            sumComp= sum(diriList[0][array_start:array_start + \
                record.windows])

            array_start = array_start + record.windows
           
            
            record.cpkm = sumComp  
            sim.set_mean_expression(record, self.diff, self.librarySize) 
   
    def set_counts(self):

        for sample in self.samples:
            
            for rec in self.references:

                record = self.references[rec]

                sim.set_sample_counts(record, sample)
                
    def gen_count_table(self):

        tabIndex = [ ]
        meanExpression = [ ]
        
        if self.diff:

            ct = {'meanExpression':[],
                  'meanExpressionDiff' : [],
                  'diff' : [],
                  'factor': [],
                  'direction': []}
        else:

            ct = {'meanExpression':[]}


        for sample in self.samples:
            
            ct[sample] = [ ]

        for rec in self.references:
            
            record = self.references[rec]
            
            tabIndex.append(record.id)

            ct['meanExpression'].append(record.meanExpression)
             
            if self.diff:

                ct['meanExpressionDiff'].append(record.meanExpressionDiff)
                ct['diff'].append(record.diff)
                ct['factor'].append(record.diffFactor)
                ct['direction'].append(record.up)
           
            for sample in self.samples:

                ct[sample].append(record.reads[sample])
        

        self.cntTab = pd.DataFrame(ct,index = tabIndex)  
            

    def write_count_table(self):

        self.cntTab.to_csv("{}/countTable.csv".format(self.outDir), sep = '\t')

    def head_count_table(self):

        print(self.cntTab.head())

    def load_quality_strings(self):  
        
        print("Load the quality strings.")        
        self.qualStringPool = {}
        self.qualStrings = seq_handling.read_fastq_qual(self.refFastq)

        for qualStr in self.qualStrings:

            lenQualStr = len(qualStr)

            if lenQualStr in self.qualStringPool:
                self.qualStringPool[lenQualStr].append(qualStr)
            else:
                self.qualStringPool[lenQualStr] = [qualStr]
        
    
    def simulate_fastq(self):

        self.fragHisto = {}

        if self.paired:
            print('paired end simulation is started.')
            for sample in self.samples:
                
                reads, self.fragHisto[sample] = sim.paired_end_simulation(self.references, 
                        self.qualStringPool[self.readLength], 
                        self.cntTab[[sample]],
                        self.readLength)
                

                seq_handling.write_paired_fastqs(reads, self.outDir, sample)

                reads = None
        else:
            print('single end simulation is started.')

            for sample in self.samples:
                
                reads = sim.single_end_simulation(self.references,
                        self.qualStringPool[self.readLength],
                        self.cntTab[[sample]],
                        self.readLength)

                seq_handling.write_fastq(reads, self.outDir, sample)
                reads = None
        writeHistoNumbers(self.fragHisto,
                '{}/histonumbers'.format(self.outDir))

def writeHistoNumbers(fragHisto, name):

    fragHisto = pd.DataFrame.from_dict(fragHisto)
    fragHisto = fragHisto.sort_index()
    fragHisto.to_csv('{}.csv'.format(name), index_label = 'fragment_length')

def dirichletList(components):

    return np.random.dirichlet(np.ones(components), size=1)

def create_folder(aim_dict):
    
    mkdir_command = "mkdir -p "
    mkdir_command = mkdir_command + aim_dict
    
    print("create the storage directory\n%s" %mkdir_command)
    
    subprocess.run(shlex.split(mkdir_command))
    
    return aim_dict
