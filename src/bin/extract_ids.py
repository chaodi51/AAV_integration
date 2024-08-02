import sys
import re

def find_genes_in_file(gene_file, data_file, output_file):
    # Read the gene names from file A into a set
    with open(gene_file, 'r') as file_a:
        ids = set(line.strip().split('\t')[1].split('.')[0] for line in file_a)

    # Open the output file for writing
    with open(output_file, 'w') as output:
        with open(data_file, 'r') as file_b:
            for line in file_b:
                # Split the line into columns
                tr_id = line.strip().split('\t')[1]
                if tr_id in ids:
                    # Write the line to the output file
                    output.write(line)

if len(sys.argv) != 4:
    print("Usage: python extract_genes.py gene_file data_file output_file")
    sys.exit(1)

gene_file = sys.argv[1]
data_file = sys.argv[2]
output_file = sys.argv[3]

# Call the function to find and write the matching rows
find_genes_in_file(gene_file, data_file, output_file)
