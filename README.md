# eGap: BWT and LCP computation for sequence collections in external memory

This software is an implementation of the eGap algorithm described in 
[*External memory BWT and LCP computation for sequence collections with applications*](https://doi.org/10.4230/LIPIcs.WABI.2018.10) by
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
./eGap --lcp dataset/reads.fastq
```

This will produce the file `dataset/reads.fastq.bwt` and `dataset/reads.fastq.2.lcp` containing the BWT and LCP array (the latter using 2 bytes per entry)


## Description

Tool to build the BWT and optionally the LCP and DA array for a collection  of sequences in external memory. There are two different usages depending on whether you already have the BWT of the input files:

* If you do have the BWTs use option -b: you must specify the file names on the command line  and use the option -o to specify an output basename. 
For example:
 `  eGap  -b --lcp -o merge  file1.bwt file2.bwt`
will produce the output files: *merge.bwt*, *merge.2.lcp*, *merge.da*. Globbing is accepted: multiple file names can be denoted for example as *file?.bwt*
 
* If you don't have the BWTs then your input must consists of a single file with extension 
  `  .fasta/.fa`  (one input document per sequence)
  `  .fastq`      (one input document per sequence)
  `  .txt`        (one input document per line)
and it is not mandatory to specify the output basename. For example:
  `   eGap --lcp  file.fasta` 
will produce the output files: *file.fasta.bwt*, *files.fasta.2.lcp*

All input and output files are uncompressed. The value 0 is used as the eof symbol in the output BWT.


## Main command line options

*-h, --help*      
  show usage

*-o, --out*        
  specify basename for output and temporary files

*-l, --lcp*          
  compute LCP Array
  
*-d, --da*          
  compute Document Array
  
*-b, --bwt*          
  inputs are bwt files (requires -o)

*-m, --mem*     
  specify available memory in MB (def. 4096). **Note:** do not assign all the available RAM to the algorithm: leave *at least* 5% to the operating system.
  
*--lbytes*      
  number of bytes for each LCP entry (def. 2)

*--dbytes*      
  number of bytes for each DA entry (def. 4)

*--rev*      
  compute data structures for the reversed string  

*-v*       
  verbose output in the log file



## Datasets


We have compared eGap with the available BWT/LCP construction tools on the following collections


Name         |SizeGBs|Num Docs    |Max DocLen|Ave DocLen|Max LCP| Ave LCP | Download Link
-------------|:-----:|-----------:|---------:|---------:|------:|--------:|-----------
Shortreads   | 8.0   | 85,899,345 | 100      | 100      | 99    | 27.90   | [.tar.gz](https://drive.google.com/open?id=199dUcf-NgCV4WaWTs96siJtibd0GsDM2)
Longreads    | 8.0   | 28,633,115 | 300      | 300      | 299   | 90.28   | [.tar.gz](https://drive.google.com/open?id=1uck1L79ERqkX4G26_-3LYYlGkw0r2Qxe)
Pacbio.1000  | 8.0   | 8,589,934  | 1,000    | 1,000    | 876   | 18.05   | [.tar.gz](https://drive.google.com/open?id=1ehqbYJmRedwiR2iLMYEP1TerkvxhhXZV)
Pacbio       | 8.0   | 942,248    | 71,561   | 9,116    | 3,084 | 18.32   | [.tar.gz](https://drive.google.com/open?id=1JER4Ci1DyZtQERqILNbrebWBXQdVQrW4)


We have also used versions of the above collections shortened to 1GB. The shortened versions can be obtained by the above files using [simple command line instructions](https://drive.google.com/open?id=1rjObN6fzXU_LrOLadCgxQ0bXh5mlTwNq). Check all [md5sums](https://drive.google.com/open?id=1CgoVBpElte6iQ6I1XkvYvi56lHqvq1kK) after dowloading and extraction.

The results of our experiments are reported on the above WABI paper and on an extended journal version (in preparation).

