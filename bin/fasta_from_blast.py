#!/usr/bin/env python3

import argparse
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a BLAST file, output the records from the FASTA or FASTQ file used as input to the BLAST."
    )
    parser.add_argument('-e', '--max_evalue', type=float, help='Records with e-value more than this will be ignored')
    parser.add_argument('-m', '--max_mismatches', type=int, help='Records with more than this many mismatches in alignment will be ignored.')
    parser.add_argument('-l', '--min_len', type=int, help='Records with an alignment length less than this will be ignored.')
    parser.add_argument('-r', '--reverse', action='store_true', help='Reverse of normal behavior - output all those records without a hit at or below the specified e-value')
    parser.add_argument('-s', '--retrieve_subjects', action='store_true', help='Retrieve subjects of BLAST hits from FASTA file instead of queries of BLAST hits.')
    parser.add_argument('-f', '--fasta_file', help='The FASTA/Q file that was used as BLAST query. If not specified, the name will be inferred from the BLAST results file.')
    parser.add_argument('blast_file', type=str, help='The BLAST results file.')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    # Determine FASTA/Q filename
    if args.fasta_file is None:
        match = re.match(r'(.*\.(fasta|fastq|fq|fa))(?!.*\.(fasta|fastq|fq|fa))', args.blast_file)
        if match:
            fasta_file = match.group(1)
        else:
            sys.stderr.write(f"Error: Couldn't infer FASTA/Q file from BLAST results file: {args.blast_file}\n")
            sys.exit(1)
    else:
        fasta_file = args.fasta_file

    # Parse BLAST file
    ids = set()
    with open(args.blast_file, 'r') as blast_fh:
        header_skipped = False
        for line in blast_fh:
            fields = line.strip().split('\t')
            if not header_skipped:
                # Skip the header line
                if any(field.lower() == 'evalue' for field in fields):
                    header_skipped = True
                    continue

            if len(fields) < 12:
                sys.stderr.write(f"Warning: Unexpected format for BLAST file. Line: {line.strip()}\n")
                continue

            id = fields[1] if args.retrieve_subjects else fields[0]
            evalue = float(fields[10])
            mismatches = int(fields[4])
            aln_length = int(fields[3])

            if args.max_evalue and evalue > args.max_evalue:
                continue
            if args.max_mismatches and mismatches > args.max_mismatches:
                continue
            if args.min_len and aln_length < args.min_len:
                continue

            if args.reverse:
                ids.discard(id)
            else:
                ids.add(id)


    # Parse FASTA/Q file and output records
    with open(fasta_file, 'r') as fasta_fh:
        is_fastq = False
        header = ""
        seq = ""

        for line in fasta_fh:
            line = line.strip()
            if line.startswith('@'):
                is_fastq = True
                if header:
                    # Print the last record
                    if args.reverse:
                        if header not in ids:
                            if is_fastq:
                                print(f"{header}\n{seq}\n+")
                            else:
                                print(f"{header}\n{seq}")
                    else:
                        if header in ids:
                            if is_fastq:
                                print(f"{header}\n{seq}\n+")
                            else:
                                print(f"{header}\n{seq}")
                # Start new record
                header = ">" + line[1:].split(' ')[0]  # Extract the ID part
                seq = ""
            elif line.startswith('>'):
                if header:
                    # Print the last record
                    if args.reverse:
                        if header not in ids:
                            if is_fastq:
                                print(f"{header}\n{seq}\n+")
                            else:
                                print(f"{header}{seq}")
                    else:
                        if header in ids:
                            if is_fastq:
                                print(f"{header}\n{seq}\n+")
                            else:
                                print(f"{header}{seq}")
                # Start new record
                header = ">"+ line[1:].split(' ')[0]  # Extract the ID part
                seq = ""
            else:
                seq += f"\n{line}"

            if is_fastq and len(line) == 0:
                # If FASTQ and a blank line, assume end of the record
                continue

        # Handle the last record
        if header:
            if args.reverse:
                if header not in ids:
                    if is_fastq:
                        print(f"{header}\n{seq}\n+")
                    else:
                        print(f"{header}{seq}")
            else:
                if header in ids:
                    if is_fastq:
                        print(f"{header}\n{seq}\n+")
                    else:
                        print(f"{header}{seq}")

if __name__ == "__main__":
    main()
