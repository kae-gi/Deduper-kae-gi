
# Strategy for a Reference Based PCR Duplicate Removal Tool
## Defining the problem

PCR duplicates are produced by a polymerase chain reaction during the amplication step during library prep. These PCR duplicates need to be removed, as they can cause bias/overrepresentation in the gene expression of our sequences. Thus is is important to normalize our data. In the RNA-seq experimental workflow, we would want to remove these duplicates from our SAM file post alignment to the genome (hence "reference based"). 

To find a PCR duplicate there are a few key tells: 
  - Same alignment position
    - same chromosome
    - same (start) position 
    - same strandedness
  - Soft clipping
  - Same Unique Molecular Index (UMI)

This tool will be designed for single-end data with 96 known UMIs, a list of which is in the repo under the ```STL96.txt``` file.

Within the SAM format, there are some fields that are important due to their correspondence with how we know we have found a PCR duplicate:
- RNAME (col 3) -> chromosome
- POS (col 4) -> position
- FLAG (col 2) -> strandedness
- CIGAR (col 6) -> soft clipping
- QNAME (col 1) -> UMI


## Examples

## Pseudocode
initialize dict for counting dupes {(chrom, stpos, str, umi): count}
open files (SAM, STL96.txt)


## High level functions

def softclip

def umi



