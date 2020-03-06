# 1. DiffReg

DiffReg is a tool to predict disease genes based on SNP-splicing association altering between case and control groups. DiffReg takes 4 input files : (1) genotype data(single nucleotide polymorphism,SNP), providing SNP type for every samples. (2) transcript isoform expression data, providing all isoforms RPKM values for the same samples. (3) gene location data, providing gene id, chromosome number, gene start and end position of isoforms corresponding gene . (4)samples information, providing grouping information of samples in genotype and isoform expression data.


# 2. Dependencies

- ***R*** environment and R package ***sQTLseekeR***. Installation of R package ***sQTLseekeR*** is provided in  [https://github.com/guigolab/sQTLseekeR ](https://github.com/guigolab/sQTLseekeR).
- ***Python3*** environment and Python package ***click***, ***pandas*** ,***numpy***. 

# 3. usage

refer to [runExample.sh](https://github.com/genemine/DiffReg/blob/master/runExample.sh)

## 3.1 details of input

To run this tool, the following input data have to prepared in specific format.

- Transcript isoform expression data. It should be a DataFrame-like text file, first two columns should be transcript isoform id and corresponding gene id, what’s more, the columns name must be ‘trId’ and ‘geneId’. For each column of the other columns is transcript isoform expression values of a person(or a sample). You’d better input the RPKM values of this expression matrix.(Reference: “./demoData_diffreg/transcript.tsv”)
- Single Nucleotide Polymorphism(SNP) data. It should be a DataFrame-like text file, first four columns should be chromosome number, SNP start position, SNP end position and SNP id. The corresponding columns name must be ‘chr’, ‘start’, ‘end’ and ‘snpId’. The following columns is genotype of each sample, and genotype is coded as 0(ref/ref),1(ref/alt),2(alt/alt),-1(missing value). What’s more, this file needs to be ordered per **chr** and **start** position. (Reference: “./demoData_diffreg/snp.tsv”)
- Gene location data. BED-like format file. It contains four columns: chromosome number, gene start position, gene end position and gene id. No column names in this file. (Reference: “./demoData_diffreg/bed.tsv”)
- Samples group information. A python dict type format file. Samples in transcript expression file and SNP data file is identical. This key-value file describes each sample belonging to disease group or control disease. 0 means control and 1 means disease. (Reference: “./demoData_diffreg/samples.txt”)

Since you have prepared the four files, put them together in a folder. Because this folder will be a parameter in command line interface. And output file will be stored in this folder.

## 3.2 command line interface parameter 

- w: Work directory that stored the four input files.
- b: File name of the gene location information with BED-like format data.
- i: File name of transcript isoform expression data
- p: File name of SNP data.
- s: File name of samples information diction.
- t: Number of threads that running DiffReg kernel code. (It may takes long time to run the permutation test,multiprocess is used in the pipeline. parameter “t” can be set as N_CPUs -1)
- -n: Times of permutation test in DiffReg. (1x10^6 is recommend for a reliable permutation test)



## 3.3 output

According to the information show in console when running this pipeline, two result file can be accessed.

### 3.3.1 diffReg_result.tsv

Result of DiffReg for all the significant SNP-Gene pairs. Here is the explanation of fields in the file.

- delta: the score of association altering between disease group and control group. And result ordered by delta.
- pv: p values of permutation test of ‘delta’
- p_adj:  adjusted p values

### 3.3.2 significant_predict_gene.tsv

Subset of “diffReg_result.tsv”. Retain lines that delta > 0 and p_adj < 0.05, meanwhile, keep the maximum delta as the score of the gene if there has two or more transcript isoforms of a Gene show in the result table. Field “gene” is considered as the significant disease genes by pipeline DiffReg. On the other hand, you can customize the thresholds in “diffReg_result.tsv” to obtain significant predicted genes in different extent.
