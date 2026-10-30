#!/usr/bin/env python

'''Reference based PCR duplicate removal for single-end reads with known UMIs.

Two reads are PCR duplicates if they match on all four of:
    1. chromosome (RNAME, column 3)
    2. strand (FLAG bit 16, column 2)
    3. true 5' start (POS column 4, adjusted using the CIGAR string in column 6)
    4. UMI (last "field" of QNAME, column 1)

The input SAM must be sorted (samtools sort, outside this script) so every read
from one chromosome sits together in one block. That's what lets seen_set get
emptied each time the chromosome changes, so memory never holds more than one
chromosome's worth of dup_info.

The first copy of a read is written out, every later copy is skipped.
Reads with a UMI that isn't in the known UMI file are thrown out.
'''

import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description="Removes PCR duplicates from a sorted, single-end SAM file of uniquely mapped reads. "
                                     "Two reads count as duplicates if they share chromosome, strand, 5' start (adjusted for soft clipping, and on the minus strand, for everything the read covers on the genome), and UMI. "
                                     "The first copy found is kept. Reads with a UMI not in the UMI file are thrown out. "
                                     "Example: ./hodapp_deduper.py -u STL96.txt -f in.sam -o out.sam")
    parser.add_argument("-f", "--file", help="absolute path to the SORTED input SAM file", required=True)
    parser.add_argument("-o", "--outfile", help="absolute path for the deduplicated output SAM file", required=True)
    parser.add_argument("-u", "--umi", help="text file of known UMIs, one per line (ex: STL96.txt)", required=True)
    return parser.parse_args()
    #-h/--help doesn't need its own line, argparse builds it automatically from the description and help= text above

args = get_args()
f = args.file #sorted input SAM
o = args.outfile #deduplicated output SAM
u = args.umi #known UMI list


def get_umi(qname: str) -> str:
    '''Takes the QNAME (column 1) of a SAM line and pulls out the UMI, which is the last field after the final ":"
    Input: NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCAAGG
    Expected output: TAGCAAGG'''
    umi = qname.split(":")[-1] #split makes a list of the pieces between colons, [-1] grabs the last piece no matter how many pieces there are
    return umi


def get_strand(flag: int) -> str:
    '''Takes the FLAG (column 2) of a SAM line and decides which strand the read is on using bit 16
    Input: 0      Expected output: +
    Input: 16     Expected output: -
    Input: 272    Expected output: -   (272 = 256 + 16, so bit 16 is set even though the number isn't 16)'''
    if (flag & 16) == 16: #& keeps only the bits that both numbers share, so this asks "is the 16 bit turned on inside flag?"
        strand = "-" #bit 16 set means the read was reverse complemented, so minus strand
    else:
        strand = "+"
    return strand


def get_five_prime_start(pos: int, cigar: str, strand: str) -> int:
    '''Takes POS (column 4), CIGAR (column 6), and strand, and returns the true 5' start of the read,
    adjusted for soft clipping and, on the minus strand, for everything the read covers on the genome
    Input: 76718927, "3S68M", "+"           Expected output: 76718924
    Input: 130171653, "40M1I30M", "+"       Expected output: 130171653
    Input: 76718924, "69M2S", "-"           Expected output: 76718994
    Input: 76875969, "2S13M470N56M", "-"    Expected output: 76876507
    Input: 76718927, "5H3S68M", "+"         Expected output: 76718924   (hard clip on the outside of the left soft clip)
    Input: 76718924, "68M3S5H", "-"         Expected output: 76718994   (hard clip on the outside of the right soft clip)'''

    #HARD CLIP (H) vs SOFT CLIP (S), since both can show up in a CIGAR:
    #both mean "these bases on the end of the read didn't line up with the genome"
    #S = the clipped bases are still in the SEQ column, just not aligned
    #H = the clipped bases were deleted from the SEQ column, the CIGAR only notes they used to be there
    #if a read has both, H is ALWAYS on the very outside, ex: 5H3S68M or 68M3S5H
    #that matters below: the code looks at the end pieces of the CIGAR for an S,
    #so if an H is sitting on the outside, it has to look one piece further in to find the S
    #H itself is never added or subtracted, those bases aren't in the read anymore,
    #and as long as every read gets the same rule, duplicates still match up


    cigar_pieces = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    #break the CIGAR into (number, letter) pairs
    #ex: "3S68M" becomes [("3", "S"), ("68", "M")]
    #ex: "5H3S68M" becomes [("5", "H"), ("3", "S"), ("68", "M")]
    #findall looks through the CIGAR string and pulls out every "number then letter" chunk, in order, as a list
    #the r"..." part is just the description of what one chunk looks like, some digits, followed by one CIGAR letter.
    #the numbers come out as STRINGS, so they need int() before any math

    if strand == "+":
        #plus strand 5' end is the LEFT end, which is where POS already is
        #only thing that can move it is a soft clip on the left, since POS skips over clipped bases
        first_number, first_letter = cigar_pieces[0] #unpacks the first pair, ex: ("3", "S") gives first_number = "3", first_letter = "S"
        if first_letter == "H":
            first_number, first_letter = cigar_pieces[1] #H is on the very outside, so any left soft clip would be the NEXT piece in
            #ex: 5H3S68M, cigar_pieces[0] is ("5", "H"), so this grabs cigar_pieces[1], which is ("3", "S")
            #without this, the 3S would get missed and the 5' start would be off by 3
        if first_letter == "S":
            five_prime_start = pos - int(first_number) #undo the left soft clip
        else:
            five_prime_start = pos

    else:
        #minus strand 5' end is the RIGHT end, so start at POS and move to the right across every reference position the read covers
        total = 0
        for number, letter in cigar_pieces:
            if letter in "MDN=X": #the letters that take up space on the genome
                total = total + int(number)
            #I is not counted, it adds bases to the read but not to the genome
            #S and H are not counted here, they didn't line up with the genome (the right S gets added back below)
            #a left S is never counted, the minus strand 5' start is on the right end
            #= and X are just M split into "matches" and "mismatches", so they count the same as M
        last_number, last_letter = cigar_pieces[-1] #[-1] = last pair
        if last_letter == "H":
            last_number, last_letter = cigar_pieces[-2] #H is on the very outside, so any right soft clip would be the piece just before it
            #ex: 68M3S5H, cigar_pieces[-1] is ("5", "H"), so this grabs cigar_pieces[-2], which is ("3", "S")
            #without this, the 3S would get missed and the 5' start would be off by 3
        if last_letter == "S":
            total = total + int(last_number) #add back the right soft clip
        five_prime_start = pos + total - 1 #move down the reference positions to the right end, -1 because pos itself is the first base

    return five_prime_start


#loop through the known UMI file to create a set: umi_set
#a set because the only question ever asked of it is "is this UMI in here?", and sets answer that fast
umi_set = set()
umi_fh = open(u, "r")
for line in umi_fh:
    umi_set.add(line.strip()) #strip takes off the \n at the end, otherwise "TAGCAAGG\n" would never match "TAGCAAGG"
umi_fh.close()


seen_set = set() #temporary bucket, holds dup_info already seen, ONLY for the current chromosome
current_chrom = None #chromosome of the last read, starts as None so the very first read always counts as a "new" chromosome
#counters for tracking final outputs:
wrong_umi_count = 0 #to count umi's not found in umi_set
duplicate_count = 0 # to count PCR duplicates
header_count = 0 #counts lines starting with @
read_count = 0 #counts every non-header line, before any filtering
unique_read_count = 0 #counts reads actually written to the output


sam_in = open(f, "r")
sam_out = open(o, "w")

#loop through input sam one line at a time, so only one line is ever held in memory
while True:
    line = sam_in.readline()
    if line == "":
        break #readline only gives back "" when the file is out of lines (a blank line in the middle would still be "\n")

    if line.startswith("@"):
        header_count += 1 #one more header line
        sam_out.write(line) #header lines get copied as-is, line still ends in \n so nothing gets added
    
    else:
        read_count += 1 #every read gets counted here, before the UMI or duplicate checks throw any out
        #split the line on tabs, then extract and store UMI and chromosome
        columns = line.split("\t") #columns[0] = QNAME, columns[1] = FLAG, columns[2] = RNAME, columns[3] = POS, columns[5] = CIGAR (python counts from 0)
        umi = get_umi(columns[0])
        chrom = columns[2]

        if umi not in umi_set:
            wrong_umi_count += 1 #unknown UMI, discard. but count for final, the line just never gets written and the loop moves on
        else:
            if chrom != current_chrom:
                seen_set = set() #new chromosome, nothing on the old one can match, so clear it
                current_chrom = chrom

            #get strand and true 5' start for this read, then bundle them with the UMI into dup_info for the duplicate check
            strand = get_strand(int(columns[1])) #every column comes out of split as text, int() turns FLAG into a number for the & math
            five_prime_start = get_five_prime_start(int(columns[3]), columns[5], strand) #same for POS
            dup_info = (umi, strand, five_prime_start) #chromosome not needed, seen_set only holds the current chromosome

            if dup_info in seen_set:
                duplicate_count += 1 #PCR duplicate, the first copy was already written, count for final output
            else:
                seen_set.add(dup_info)
                sam_out.write(line) #the original line, untouched, so the output stays a properly formatted SAM
                unique_read_count += 1 #one more read kept

sam_in.close()
sam_out.close()


print("Header lines:", header_count)
print("Total reads (non-header lines):", read_count)
print("Unique reads (written out):", unique_read_count)
print("Wrong UMIs (discarded):", wrong_umi_count) #reads whose UMI wasn't in umi_set
print("PCR duplicates (removed):", duplicate_count) #reads whose dup_info was already in seen_set