#!/usr/bin/env python
"""
Author: Kaetlyn Gibson
Deduper
Bi624
"""
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="A program for removing all PCR duplicates (retains a single copy of each read) \
                                                  given a sorted SAM file of uniquely mapped reads.")
    parser.add_argument('-f', '--infile', help="Designates absolute file path to sorted sam file.", type=str, required=True)
    parser.add_argument('-o', '--outfile', help="Designates absolute file path to deduplicated sorted sam file.", type=str, required=True)
    parser.add_argument('-u', '--umi', help="Designates file containing the list of UMIs.", type=str, required=True)
    return parser.parse_args()	
args = get_args()

def getKnownUMIs(umi_filename:str) -> set:
    """
    Given a file of known UMIs as input, returns a set
    of the known UMIs.
    """
    with open(umi_filename, 'r') as f:
        return set(line.strip() for line in f)

def recordParser(record: str) -> tuple:
    """
    Given a record from the sorted SAM file, parses the record.
    Returns a tuple consisting of the record's:
    (chromosome, adjusted starting position, strandedness, umi)
    """
    umi = record[0].split(':')[-1]
    chrom = record[2]
    strand = strandednessCheck(int(record[1]))
    pos = int(record[3])
    cigar = record[5]
    stpos = stposCheck(strand, pos, cigar)
    return chrom, stpos, strand, umi


def strandednessCheck(flag:int) -> bool:
    """
    Given an input of the SAM bitwise FLAG field, returns
    True if reverse complement or False if forward strand.
        16 = 16, - strand (rev comp) -> True
        0 = 16, + strand (fwd) -> False
    """
    return (flag & 16) == 16

def stposCheck(strand:bool, pos:int, cigar:str) -> int:
    """
    Given the strandedness, position, and CIGAR string as input,
    returns the adjusted start position.
        Forward strand (strand = False):
            In the CIGAR string, only subtract the value of the leftmost 'S' to the position.
            EX: 20S10M -> pos - 20 
        Reverse strand (strand = True):
            In the CIGAR string, only add the values of rightmost 'S', 
            any values attached to an 'D', any values attached to a 'M',
            and any values attached to a 'N'.
            EX: 3S47M6I111D1000N30M6S -> pos + 47 + 111 + 1000 + 30 + 6 
    """
    dmns_search = re.findall('\d+[DMNS]', cigar)
    if not strand:
        if 'S' in dmns_search[0]:
            pos -= int(re.findall('\d+', dmns_search[0])[0])
    else:
        if 'S' in dmns_search[-1]:
            pos += int(re.findall('\d+', dmns_search[-1])[0])
        pos += sum([int(re.findall('\d+', i)[0]) for i in dmns_search if 'S' not in i])
    return pos

# initialize
known_UMIs = getKnownUMIs(args.umi)
unknown_umis, dups_removed, header_lines, unique_records, total_records = 0, 0, 0, 0, 0
# deduplicate process
with open(args.outfile, 'w') as outf, open(args.infile, 'r') as inf:
    curr_chrom = ''
    curr_chrom_records = set()
    for line in inf:
        temp_line = line.strip()
        # sort & write header lines out
        if temp_line[0] == '@':
            header_lines += 1
            outf.write(line)
        # sort & write record lines out
        elif temp_line[0] != '@':
            identifier = recordParser(temp_line.split('\t'))
            total_records += 1
            # umi is unknown, ignore the record
            if identifier[3] not in known_UMIs:
                unknown_umis += 1
            # umi is known, deduplicate
            else:
                if curr_chrom != identifier[0]:
                    curr_chrom = identifier[0]
                    curr_chrom_records.clear()
                if identifier not in curr_chrom_records:
                    curr_chrom_records.add(identifier)
                    unique_records += 1
                    outf.write(line)
                else:
                    dups_removed += 1 

print(f'Number of unknown UMIs found: {unknown_umis}')
print(f'Number of duplicates removed: {dups_removed}')
print(f'Number of unique records: {unique_records}')
print(f'Number of header lines: {header_lines}')
print(f'Number of total records: {total_records}')
