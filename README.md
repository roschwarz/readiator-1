# Readiator

The readiator is a script to simulate RNA-Seq data. The basic idea is to take
reference sequences and extract a number sub-sequence, with a requested read 
length, out of it. The tool needs a fastq-file from a real experiment to 
simulate sequencing errors. The script can simulate single end as well as 
paired end data sets with replicates. Additionally, it
is possible to generate data sets that contain differentially expressed genes.

# Arguments

-f references.fa
-fq fastq-file
-l read length
-r library size (default 1.000.000)
-sz number of replicates
-s seed (default random number between 1 and 1000)
-t filter for min length of reference sequences that are considered (default 100)
-o output directory (default simulated_fastq)
-p to run a paired-end simulation
-d simulate a data set that contains differentially expressed genes

Example:

single end data:

python3.7 readiator.py -f test_data_set/referenceSequences.fa 
-fq test_data_set/example.fastq
-l 50
-r 10000
-sz 3
-s 5
-t 100

paired end data:

python3.7 readiator.py -f test_data_set/referenceSequences.fa 
-fq test_data_set/example.fastq
-l 50
-r 10000
-sz 3
-s 5
-t 100
-p

Important is that the fastq file contains quality strings with a length at least
as long as the requested read length

## Generation of random read counts

Let **n_g** be the number of kilo bases per gene **g**. such that

$$
\sum_i^{g} n_{g} = N
$$

Thus, *N* is the number of kilo bases of all genes to be simulated, denoted as 
components in the following. We assume, that the first **k** components belong 
to the first **q** genes, i.e. the components are sorted with respect to the 
genes.  Now, get a random dirichlet decomposition $D$ of all $N$ components 
with size 1. In the following, we sum up over the decomposition for each gene. 
Such that,

$$
d_g = \sum_{j=0}^{n_g} D[o+j]
$$

and $d_g$ is the sum of the elements of the decompositions for all components 
of gene $g$. Since the components are ordered, the offset is defined as 

$$
o = \sum_{k=0}^{g-1} n_g
$$

such that $\sum D = 1$. Now, we have a $d_g$ for each gene and only need to 
multiply the value $d_g$ with the desired number of reads in the formulate, 
i.e. $10^6$ (-r).



