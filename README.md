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
  specify memory assigned to the algorithm in MB. Default is 95% of the available RAM

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

### Document array requirements

If you want to simultaneaouly **merge BWT files** and compute the **Document Array** for each input BWT you must provide, in addition to the DA, also a `.docs` file for containing the number of documents in the file in 64 bits little endian format. The `.docs` file is automatically computed when the option *-d, --da* is used.

**Example**

```sh
./eGap -m 4096 dataset/file1.fastq -o file1 --da
./eGap -m 4096 dataset/file2.fastq -o file2 --da

./eGap -m 4096 --bwt -o merge file1.bwt file2.bwt --da
```

The first two commands compute `file1.bwt`, `file1.da`, `file1.docs` and `file2.bwt`, `file2.da`, `file2.docs` which are used by the third command to compute `merge.bwt`, `merge.da`, and `merge.docs`


## Applications

### Truncated LCP values

The running time of eGap can be decreased if, instead of the true 
LCP values, the user settles for computing the LCP values up to a certain 
threshold *k*. Using the option *--trlcp k*, as an altenative to *--lcp*, 
the algorithm computes an LCP array in which all values greater than *k* are
replaced by the value *k*.

*--trlcp k*
  truncate LCP values to the value *k*

### Bruijn graph info (BOSS)

Another option offered by eGap, is to compute the info required for the
construction of the succinct (BOSS) representation of the de Bruijn graph
associated to the input sequences. 
Using the option *--deB k* eGap compute two bitfiles with extension 
.lcpbit0 and .lcpbit1 which, together with the BWT, can be used to compute 
the BOSS representation of the de Bruijn graph as described in the 
[Application](https://almob.biomedcentral.com/articles/10.1186/s13015-019-0140-0#Sec14)
section of the [AMB paper](https://doi.org/10.1186/s13015-019-0140-0). 

*--deB K*
  compute info for building the order-K deBruijn graph

**Notice:** if the options *--trlcp* or *--deB* are used, suffixes are sorted only up the first *k* symbols so the resulting BWT *will not* be the standard one.

### Quality score (QS) sequences

eGap can output the Quality Score (QS) sequences of a FASTQ file permuted according to the BWT symbols (allowed only for `.fastq` files).

*--qs*
  QS permuted according to the BWT

**Example**

```sh
./eGap -m 4096 dataset/reads.fastq -o file1 --qs
```

The output files are `file1.bwt`, `file1.bwt.qs`, and `file1.eGap.log`

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

## Thanks

Thanks to [Pierre Peterlongo](https://github.com/pierrepeterlongo) by helpful suggestions and debugging.
