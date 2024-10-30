import argparse
import os
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Process .tally files and output a tab-delimited matrix.')
    parser.add_argument('-t', '--taxids_file', type=str, help='Limit output to the taxids listed in this file.')
    parser.add_argument('-c', '--tally_cutoff', type=float, default=0, help='Only output a read count value if value is > this cutoff.')
    parser.add_argument('-v', '--virus_only', action='store_true', help='Only include virus taxids in matrix.')
    parser.add_argument('-p', '--no_phage', action='store_true', help='Exclude phage taxids.')
    parser.add_argument('tally_files', nargs='+', help='Tally files to process.')
    
    return parser.parse_args()

def read_taxids(taxids_file):
    taxid_name_map = {}
    with open(taxids_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            taxid = parts[0]
            name = parts[1]
            taxid_name_map[taxid] = name
    return taxid_name_map

def main():
    args = parse_args()
    
    taxid_name_map = {}
    taxid_common_name_map = {}
    taxids = []
    barcodes = []
    matrix = defaultdict(lambda: defaultdict(float))
    
    if args.taxids_file:
        taxid_name_map = read_taxids(args.taxids_file)
        taxids = list(taxid_name_map.keys())
    
    observed_barcodes = set()
    
    for tally_file in args.tally_files:
        barcode_match = re.search(r'_([ACGT]{7})_', tally_file)
        barcode = barcode_match.group(1) if barcode_match else tally_file
        
        if barcode not in observed_barcodes:
            observed_barcodes.add(barcode)
            barcodes.append(barcode)

        with open(tally_file, 'r') as f:
            header_line = True
            for line in f:
                if (header_line):
                    header_line = False  
                    continue   
                parts = line.strip().split("\t")
                
                if len(parts) < 5:
                    continue
                
                if len(parts) == 7:  # rank_column is present
                    taxid, name, common_name, kingdom, tally = parts[1], parts[3], parts[4], parts[5], float(parts[6])
                else:  # no rank_column
                    taxid, name, common_name, kingdom, tally = parts[0], parts[1], parts[2], parts[3], float(parts[4])

                if not taxid:
                    taxid = "X"
                    name = "Unknown Taxid"
                    kingdom = "X"

                if args.virus_only and kingdom != "Viruses":
                    continue

                if args.no_phage and "phage" in name.lower():
                    continue

                if taxid not in taxid_name_map:
                    taxid_name_map[taxid] = name
                    taxid_common_name_map[taxid] = common_name
                    if taxid not in taxids:
                        taxids.append(taxid)

                matrix[barcode][taxid] += tally
    
    # Filter taxids based on cutoff
    passing_taxids = [taxid for taxid in taxids if any(matrix[barcode][taxid] > args.tally_cutoff for barcode in barcodes)]
    
    # Output the matrix
    # First line (taxids)
    print("\t" + "\t".join(passing_taxids))
    
    # Second line (taxid names)
    print("\t" + "\t".join(taxid_name_map[taxid] for taxid in passing_taxids))
    
    # One line for each barcode
    for barcode in barcodes:
        line = [barcode] + [
            f"{matrix[barcode][taxid]:.1f}" if matrix[barcode][taxid] > args.tally_cutoff else ""
            for taxid in passing_taxids
        ]
        print("\t".join(line))

if __name__ == '__main__':
    main()
