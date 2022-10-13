
# Strategy for a Reference Based PCR Duplicate Removal Tool
## Defining the problem

PCR duplicates are produced by a polymerase chain reaction during the amplication step during library prep. These PCR duplicates need to be removed to only retain a single copy of each read, as they can cause bias/overrepresentation in the gene expression of our sequences. Thus is is important to normalize our data. In the RNA-seq experimental workflow, we would want to remove these duplicates from our SAM file post alignment to the genome (hence "reference based"). 

To find a PCR duplicate there are a few key things to look for: 
  - Same alignment position
    - same chromosome
    - same (start) position 
    - same strandedness
  - Soft clipping (based on unclipped position)
  - Same Unique Molecular Index (UMI)

This tool will be designed for single-end data with 96 known UMIs, a list of which is in the repo under the ```STL96.txt``` file.

Within the SAM format, there are some fields that are important due to their correspondence with how we know we have found a PCR duplicate:
- RNAME (col 3) -> chromosome
- POS (col 4) -> position
- FLAG (col 2) -> strandedness
- CIGAR (col 6) -> soft clipping
- QNAME (col 1) -> UMI

## Pseudocode
```
possibly sort the input SAM file first based on alignment position elements?
output SAM maybe 3 files?: a "deduped" SAM, a dupe SAM, and UMI error SAM

open files (input SAM, output SAMs, STL96.txt)
    create set of all UMI from STL96.txt
    initialize empty dict for counting dupes {(chrom, stpos, strand, umi): count}
    from input SAM, lines not including ^@:
        for each line:
            strip \n and split (on \t?)
            umi = line[0], then split based on ':' and take [-1]
            strand = strandednessCheck(line[1])
            chrom = line[2]
            pos = line[3]
            cigar = line[5]
            stpos = stposCheck(strand, pos, cigar)
            check if umi not in umiset:
                add line to UMI error SAM
            else:
                check if (chrom, stpos, strand, umi) not in dupe dict
                    add line to deduped SAM
                    add {(chrom, stpos, strand, umi):0} to dict for counting dupes
                else:
                    add line to dupe SAM
                    {(chrom, stpos, strand, umi):+1}
close files
```
## High level functions
```
def strandednessCheck(int: FLAG) -> int:
    """ 
    Given the SAM bitwise FLAG column input, checks to see if fwd or rev strand.
    16 = 16, - strand (rev comp) -> 1
    0 = 16, + strand (fwd) -> 0
    
    EX 1:
        NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M	*	0	0	\
        TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	\
        6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	\
        MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
        -> strandednessCheck(0)
        -> 0
    EX 2:
        NS500451:154:HWKTMBGXX:1:11101:25533:1187:GTTCACCT	16	2	76743835	36	71M	*	0	0	\
        CTTGGTAACTTTCAGAGAATTAGTCACAACTTCTGAAGCAACCACAGTCCATGCAAGTCGACTGGTTTCTC	\
        6AEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEE<EEEEEEEEEEEEEEEEEEEE	\
        MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
        -> strandednessCheck(16)
        -> 1
    """
    return 1 or 0
```
```
def stposCheck(strand, pos, cigar) ->:
    """ 
    Given the strandedness, position, and CIGAR string as input, returns the starting position.
    If strand is 0 (fwd strand), 
        stpos = pos (no S in CIGAR, M)
        stpos = pos - x(S in CIGAR, M)
    Else if strand is 1 (rev comp strand),
        stpos = pos (no start S in CIGAR) + y(M, N, etc.) - 1 (no end S in CIGAR)
        stpos = pos (no start S in CIGAR) + y(M, N, etc.) + z(end S in CIGAR) - 1
        stpos = pos + x(start S in CIGAR) + y(M, N, etc.) - z(end S in CIGAR)
        
    EX 1:
        NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	3	76814284	36	71M	*	0	0	\
        TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	\
        6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	\
        MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
        -> stposCheck(0, 76814284, '71M')
        -> 76814284
    EX 2:
        NS500451:154:HWKTMBGXX:1:11101:25533:1187:GTTCACCT	16	2	76743835	36	71M3S	*	0	0	\
        CTTGGTAACTTTCAGAGAATTAGTCACAACTTCTGAAGCAACCACAGTCCATGCAAGTCGACTGGTTTCTC	\
        6AEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEE<EEEEEEEEEEEEEEEEEEEE	\
        MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
        -> stposCheck(1, 76743835, '71M3S')
        -> 76743835 + 71 + 3 - 1 = 76743908
    """
  return stpos
```


