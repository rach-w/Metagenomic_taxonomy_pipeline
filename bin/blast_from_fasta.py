#!/usr/bin/env python3

import sys
import argparse

def parse_fasta(file):
    ids = set()
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                id = line[1:].split()[0]
                ids.add(id)
    return ids

def process_blast(blast_file, fasta_ids, output_misses):
    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip()
            fields = line.split('\t')
            if len(fields) < 12:
                print(f"error: unexpected format for blast file. Line: {line}", file=sys.stderr)
                continue
            id = fields[0]

            if (id in fasta_ids and not output_misses) or (id not in fasta_ids and output_misses):
                print(line)

def main():
    parser = argparse.ArgumentParser(description='Filter BLAST results based on FASTA file entries.')
    parser.add_argument('-f', '--fasta_file', type=str, help='FASTA file that was used as BLAST query.')
    parser.add_argument('-r', '--reverse', action='store_true', help='Output all those BLAST records without a FASTA counterpart.')
    parser.add_argument('blast_file', type=str, help='BLAST results file.')

    args = parser.parse_args()

    if not args.fasta_file:
        # Try to guess the fasta filename based on the blast results filename
        if args.blast_file.endswith('.bx_nr'):
            fasta_file = args.blast_file.rsplit('.', 2)[0] 
        elif args.blast_file.endswith('.blast'):
            fasta_file = args.blast_file.rsplit('.', 2)[0] + '.fa'
        else:
            print(f"error: unable to guess FASTA file name based on BLAST file name: {args.blast_file}", file=sys.stderr)
            sys.exit(1)
    else:
        fasta_file = args.fasta_file

    try:
        fasta_ids = parse_fasta(fasta_file)
    except FileNotFoundError:
        print(f"error: couldn't open FASTA file: {fasta_file}", file=sys.stderr)
        sys.exit(1)

    try:
        process_blast(args.blast_file, fasta_ids, args.reverse)
    except FileNotFoundError:
        print(f"error: couldn't open BLAST results file: {args.blast_file}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
