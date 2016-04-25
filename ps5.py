#! /usr/bin/env python

from collections import Counter, defaultdict
import pybedtools

# Problem 1 (lamina.bed)
loc = "../data-sets/bed/lamina.bed"
bedfile = pybedtools.BedTool(loc)

# 1.1: what chromosome has the region w/ the largest start position (col2)?

lg_start = 0

for record in bedfile:
    if lg_start < record.start:
        lg_start = record.start
        chrom_start = record.chrom
print 'answer-1.1', chrom_start

# 1.2: Report the region w/ largest end position on ChrY as chrom start
# end value region_length

lg_end = 0

for record in bedfile:
    if record.chrom == 'chrY':
        if lg_end < record.end:
            lg_end = record.end
            chrom_stop = record.chrom
            region_start = record.start
            region_value = record.fields[3]
            region_length = (lg_end - region_start)
print 'answer-1.2', chrom_stop, region_start, lg_end, region_value,
region_length

# Problem 2 (SP1.fq)
sp1 = "../data-sets/fastq/SP1.fq"

# 2.1: Which of the 1st 10 sequence records has the most Cs? Report record
# name.

# define a funciton to parse fastq and get name, seq & q score per record:
# M made this nice little modification to pass # of records as 2nd
# argument to parse_fq, totally borrowing that because it's clever.

def parse_fq(sp1, num_of_records):
    line_num = 0
    num_records = 0
    for line in open(sp1):
        if line_num > (num_of_records *4): break
        line_type = line_num % 4
        if line_type == 0:
            name = line.strip()
        elif line_type == 1:
            seq =line.strip()
        elif line_type == 3:
            quals = line.strip()

            yield name, seq, quals

        line_num += 1

# then a new function to convert ASCII qualscores into numbers

def sum_quals(qual):
    sum = 0
    for char in quals:
        sum += ord(char)
    return sum


# and lastly a reverse complement function to reversec sequences

def reversec(seq):
    comps = []
    for char in seq:
        if char == 'A':
            comps.append('T')
        elif char == 'T':
            comps.append('A')
        elif char == 'C':
            comps.append('G')
        elif char == 'G':
            comps.append('C')
    return ''.join(reversed(comps))

# now to the questions:

max_C = 0
for name, seq, quals in parse_fq(sp1, 10):
    if max_C < Counter(seq)['C']:
        max_C = Counter(seq)['C']
        name_max_C = name
print 'answer-2.1:', name_max_C

# 2.2: For each record, convert each quality score character to a number
# (using ord) and report the largest total quality score.

Quality = 0
for name, seq, quals in parse_fq(sp1, 10):
    if Quality < sum_quals(quals):
        Quality = sum_quals(quals)
print 'answer-2.2:', Quality

# 2.3: Report the reverse complement of each of the 1st 10 sequences.

rev_comp = []
for name, seq, quals in parse_fq(sp1, 10):
    rev_comp.append(reversec(seq))

print 'answer-2.3:', rev_comp
