# eGap: BWT and LCP computation for sequence collections in external memory

This software is an implementation of the eGap algorithm described in 
*External memory BWT and LCP computation for sequence collections with applications* by
L. Egidi, F. A. Louza, G. Manzini, and G. P. Telles. Copyright 2017-2018 by the authors. 


## Prerequisites

* A relatively recent version of *gcc*
* Python 3.X


## Install

```sh
git clone https://github.com/felipelouza/egap.git
cd egap
```

## Compile

```sh
make 
```

## Quick test

```sh
eGap -l dataset/reads.fastq
```

This will produce the file `dataset/reads.fastq.bwt` and `dataset/reads.fastq.2.lcp` containing the BWT and LCP array (the latter using 2 byte per entry)


## Description


Tool to build the BWT and optionally the LCP and DA array for a collection  of sequences in external memory. There are two different usages depending on whether you already have the BWT of the input files:

* If you do have the BWTs use option -b: you must specify the file names on the command line  and use the option -o to specify an output basename. 
For example:
 `  eGap  -bl  -o merge  file1.bwt file2.bwt`
will produce the output files: *merge.bwt*, *merge.2.lcp*, *merge.da*. Globbing is accepted: multiple file names can be denoted for example as *file?.bwt*
 
* If you don't have the BWTs then your input must consists of a single file with extension 
 `   .fasta`  (one input document per sequence)
  `  .fastq`  (one input document per sequence)
  `  .txt`    (one input document for line)
and it is not mandatory to specify the output basename. For example:
  `   eGap -l  file.fasta` 
will produce the output files: *file.fasta.bwt*, *files.fasta.2.lcp*

All input and output files are uncompressed. The value 0 is used as the eof symbol in the output BWT.


## Command line options

*-h, --help*      
  show usage

*-o, --out*        
  specify basename for output and temporary files

*-l, --lcp*          
  compute LCP Array
  
*-b, --bwt*          
  inputs are bwt files (requires -o)

*-m, --mem*     
  specify available memory in MB (def. 4096)

*--lbytes*      
  number of bytes for each LCP entry (def. 2)

*-v*       
  verbose output in the log file
 
