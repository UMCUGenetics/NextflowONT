#! /usr/bin/env python3
import sys
import argparse

def get_homopolymers(fasta_file, homopolymer_length):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=argparse.FileType('r'), help='path to reference genome (fasta). Note: fasta will be loaded into memory')
    parser.add_argument('length', type=int, help='minimum length of homolopymer to consider [int]')
    args = parser.parse_args()
    get_homopolymers(args.reference, args.length)
