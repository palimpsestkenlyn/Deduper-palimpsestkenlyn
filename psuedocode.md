# Deduper

## The problem:
During library prep, PCR amplification makes many copies of the same original molecule. While intended for sequencing, in some types of analysis (when quantification is desired, like differential expression)these copies inflate read counts and bias downstream analyses. The problem is how to identify and filter out PCR duplicates. The goal of this tool is to take a sorted SAM file of uniquely mapped, single-end reads and write out only one copy of each original molecule (the first one encountered).

Defining/Identifying: Two reads are PCR duplicates if they match on all four of:

 * Chromosome (RNAME, column 3)
 * Strand (FLAG, column 2, bit 16)
 * 5' start position 
 * UMI (last field of QNAME, column 1)

Reads whose UMI is not in the list of 96 known UMIs (STL96.txt) are discarded. 

### Input: Sorted SAM file with PCR duplicates
### Output: Sorted SAM file excluding PCR duplicates

## Functions

```
def get_umi(qname: str) -> str:
    '''take the QNAME (column 1) of a SAM line and pull out the UMI, which is the last field after the final ":"'''
    split qname on ":" and keep the last piece
    returns the UMI as a string
Input: NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG
Expected output: TAGCAAGG
```
```
def get_strand(flag: int) -> str:
    '''take the FLAG (column 2) of a SAM line and decide which strand the read is on using bit 16'''
    if (flag & 16) == 16:
        strand = "-" #bit 16 set means the read was reverse complemented, so minus strand
    else:
        strand = "+"
    returns "+" or "-"
Input: 0
Expected output: +
Input: 16
Expected output: -
```

logic for calculating 5' position. where on the genome does this read really start?

That comes down to two things:

 * Which end is the 5' end? 
    * Plus strand: the left end
    * Minus strand: the right end
 * How far is that end from POS?
    * Plus: POS is the left => 5' end on + strand. When bases are soft clipped on the left POS != 5' position. If left soft clipped, subtract or back up by the genome positions by the number of soft clipped bases to find true 5'
    * Minus: POS is the left => 3' end on - strand. have to go from POS all the way to the right end (up genome positions), counting everything the read covers on the genome (including right soft clips), to find the 5' position which is the right most end base. Add up all the occupied/marked positions from POS going right. Left soft clips don't count/matter. Insertions don't count(?). 5' position = POS + that total - 1 (minus 1 because POS itself is the first position counted)

```
def get_five_prime_start(pos: int, cigar: str, strand: str) -> int:
    '''take POS (column 4), CIGAR (column 6), and strand, and return the true 5' start of the read, adjusted for soft clipping and, on the minus strand, for everything the read covers on the genome'''
    break the CIGAR into number/letter pairs #ex: 3S68M becomes (3, S) and (68, M)
    if strand == "+":
        if the first letter is S:
            five_prime_start = pos - that number #undo the left soft clip
        else:
            five_prime_start = pos
    else strand == "-":
        total = add up the numbers for every M, D, and N #the letters that take up space on the genome
        if the last letter is S:
            total = total + that number/of soft clipped bases #add back the right soft clip
        five_prime_start = pos + soft clip total - 1 #-1 because pos itself is the first base
        #I is not counted, it adds bases to the read but not to the genome coverage/positions
        #a left S is not counted, the minus strand 5' start is on the right end
    returns five_prime_start as an int
Input: 76718927, "3S68M", "+"
Expected output: 76718924
Input: 130171653, "40M1I30M", "+"
Expected output: 130171653
Input: 76718924, "69M2S", "-"
Expected output: 76718994
Input: 76875969, "2S13M470N56M", "-"
Expected output: 76876507
```

## Pseudocode 

Loop through the known UMI file (STL96.txt) to create a set: umi_set

dup_info = (umi, strand, five_prime_start) for a read to be held for comparison
seen_set = set() #temporary bucket, holds dup_info already seen, ONLY for the current chromosome
current_chrom = None #chromosome of the last read, used to know when to clear seen_set and track each section/comparison

Open input sam to read and output sam to write

Loop through input sam one line at a time:
    if the line starts with "@":
        write the line to output sam #header lines get copied as-is
    else:
        split the line on tabs, then extract and store UMI and chromosome number
        umi = get_umi(QNAME, column 1)
        chrom = RNAME (column 3)

        if umi not in umi_set:
            do not write the read #unknown UMI, discard
        else:
            if chrom != current_chrom:
                seen_set = set() #new chromosome, nothing on the old one can match, so clear it
                current_chrom = chrom

            #extract strand and factor true 5' start for this read, then bundle them with umi into dup_info for the duplicate check=> dup_info = (umi, strand, five_prime_start)
            strand = get_strand(FLAG, column 2)
            five_prime_start = get_five_prime_start(POS column 4, CIGAR column 6, strand)
            dup_info = (umi, strand, five_prime_start) #chromosome not needed, seen_set only holds the current chromosome

            if dup_info in seen_set:
                do not write the read #PCR duplicate, the first copy was already written, skip read go to next
            else:
                add dup_info to seen_set #to store for looking at duplicates
                write the line to output sam #write out the first occurence 

Close both files