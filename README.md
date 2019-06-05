# eGap: BWT and LCP computation for sequence collections in external memory

This software is an implementation of the eGap algorithm described in 
[*External memory BWT and LCP computation for sequence collections with 
applications*](https://doi.org/10.1186/s13015-019-0140-0) by
L. Egidi, F. A. Louza, G. Manzini, and G. P. Telles, Algorithms for 
Molecular Biology (2019).

Copyright 2017-2019 by the authors. 


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
./eGap --lcp -m 4096 dataset/reads.fastq
```

This will produce the file `dataset/reads.fastq.bwt` and `dataset/reads.fastq.2.lcp` containing the BWT and LCP array (the latter using 2 bytes per entry). The computation will use 4GB (4096 MB) of RAM


## Description

Tool to build the **BWT** and optionally the **LCP, DA and SA array** for a collection  of sequences in external memory. There are two different usages depending on whether you already have the BWT of the input files:

* If you do have the BWTs use option -b: you must specify the file names on the command line  and use the option -o to specify an output basename. 
For example:
 `  eGap  -b --lcp -o merge -m 4096 file1.bwt file2.bwt`
will produce the output files: *merge.bwt*, *merge.2.lcp*, *merge.da*. Globbing is accepted: multiple file names can be denoted for example as *file?.bwt*
 
* If you don't have the BWTs then your input must consists of a single file with extension 
  `  .fasta/.fa`  (one input document per sequence)
  `  .fastq`      (one input document per sequence)
  `  .txt`        (one input document per line)
and it is not mandatory to specify the output basename. For example:
  `   eGap --lcp -m 4096 file.fasta` 
will produce the output files: *file.fasta.bwt*, *files.fasta.2.lcp*

All input and output files are uncompressed. The value 0 is used as the eof symbol in the output BWT.


## Main command line options

*-m, --mem*     
  specify memory assigned to the algorithm in MB. This is a *mandatory* option. **Note:** do not assign all the available RAM to the algorithm: leave at least 5% to the operating system.

*-o, --out*        
  specify basename for output and temporary files

*-b, --bwt*          
 inputs are bwt files (requires -o)

*-l, --lcp*          
  compute LCP Array
 
*--rev*      
  compute data structures for the reversed string  
    
*--lbytes*      
  number of bytes for each LCP entry (def. 2)

*-v*       
  verbose output in the log file

*-h, --help*      
  show usage


## Suffix array and document array computation 

Use the options: 

*-d, --da*          
  compute Document Array
  
*-s, --sa*          
  compute Suffix Array

*--dbytes*      
  number of bytes for each DA entry (def. 4)

*--sbytes*      
  number of bytes for each SA entry (def. 4)

### Merging BWT files

In the case you want to **merge BWT files** and later compute the **Document Array**, you must provide DA and `file.docs` with the following option:

*--da --docs*      
  compute Document Array and output the number of documents into `file.docs` (required to use --bwt --da)

**Example**

```sh
./eGap -m 4096 dataset/file1.fastq -o file1 --da --docs
./eGap -m 4096 dataset/file2.fastq -o file2 --da --docs

./eGap -m 4096 --bwt -o merge file1.bwt file2.bwt --da
```


## Truncated LCP values and de Bruijn graph info 

The running time of eGap can be decreased if, instead of the true 
LCP values, the user settles for computing the LCP values up to a certain 
threshold *k*. Using the option *--trlcp k*, as an altenative to *--lcp*, 
the algorithm computes an LCP array in which all values greater than *k* are
replaced by the value *k*.

*--trlcp*
  compute LCP values only up to TRLCP (truncated LCP)


Another option offered by eGap, alternative to (truncated) LCP, 
is to compute the info required for the construction of the succinct (BOSS)
representation of the de Bruijn graph associated to the input sequences. 
Using the option *--deB k* eGap compute two bitfiles with extension 
.lcpbit0 and .lcpbit1 which, together with the BWT, can be used to compute 
the BOSS representation of the de Bruijn graph as described in the 
[Application](https://almob.biomedcentral.com/articles/10.1186/s13015-019-0140-0#Sec14)
section of the [AMB paper](https://doi.org/10.1186/s13015-019-0140-0). 

*--deB*      
  compute info for order-DEB deBruijn graph


## Datasets


We have compared eGap with the available BWT/LCP construction tools on the following collections


Name         |SizeGBs|Num Docs    |Max DocLen|Ave DocLen|Max LCP| Ave LCP | Download Link
-------------|:-----:|-----------:|---------:|---------:|------:|--------:|-----------
Shortreads   | 8.0   | 85,899,345 | 100      | 100      | 99    | 27.90   | [.tar.gz](https://drive.google.com/open?id=199dUcf-NgCV4WaWTs96siJtibd0GsDM2)
Longreads    | 8.0   | 28,633,115 | 300      | 300      | 299   | 90.28   | [.tar.gz](https://drive.google.com/open?id=1uck1L79ERqkX4G26_-3LYYlGkw0r2Qxe)
Pacbio.1000  | 8.0   | 8,589,934  | 1,000    | 1,000    | 876   | 18.05   | [.tar.gz](https://drive.google.com/open?id=1ehqbYJmRedwiR2iLMYEP1TerkvxhhXZV)
Pacbio       | 8.0   | 942,248    | 71,561   | 9,116    | 3,084 | 18.32   | [.tar.gz](https://drive.google.com/open?id=1JER4Ci1DyZtQERqILNbrebWBXQdVQrW4)


We have also used versions of the above collections shortened to 1GB. The shortened versions can be obtained by the above files using [simple command line instructions](https://drive.google.com/open?id=1rjObN6fzXU_LrOLadCgxQ0bXh5mlTwNq). Check all [md5sums](https://drive.google.com/open?id=1CgoVBpElte6iQ6I1XkvYvi56lHqvq1kK) after dowloading and extraction.

The results of our experiments are reported on the above AMB paper.

