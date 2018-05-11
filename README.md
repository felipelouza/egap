# egap

This code is an implementation of eGAP, which computes the Burrows-Wheeler transform (BWT) and the LCP-array for a string collection in external memory.

## Install

```sh
git clone https://github.com/felipelouza/egap.git
cd egap
```

## Compile

```sh
make 
```

## Run:

To run a test with K=100 strings from INPUT=dataset/input.txt using RAM=1024MB, type:

```sh
make 
make run INPUT=dataset/input.txt K=100 RAM=1024
```

Obs. K=0 gives all strings in INPUT.

## Output

The multi-string BWT and LCP-array are stored in the same directory of the input file with the names:

```sh
INPUT.2.lcp
INPUT.bwt
```
