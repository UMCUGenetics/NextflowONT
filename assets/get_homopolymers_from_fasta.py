#! /usr/bin/env python3
import sys

fasta_file = open(sys.argv[1], 'r')
homopolymer_length = int(sys.argv[2])
chromosome = base = start = end = -1

for line in fasta_file:
    if line[0] == '>':
        if end - start > homopolymer_length:
            print('\t'.join(map(str, [chromosome, start, end, end - start])))
        chromosome = line.strip().split()[0][1:]
        base = -1
        start = end = 0
    else:
        for b in line.strip():
            if b != base:
                if end - start > homopolymer_length:
                    print('\t'.join(map(str, [chromosome, start, end, end - start])))
                start = end
                base = b
            end += 1

if end - start > homopolymer_length:
    print('\t'.join(map(str, [chromosome, start, end, end - start])))
