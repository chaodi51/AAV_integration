import sys

def extract_lines_not_in_linker(all_read_file, linker_read_file, nolinker_read_file):
    # Create a set to store the lines from linker_read_file
    linker_lines = set()
    with open(linker_read_file, 'r') as linker_read:
        for line in linker_read:
            linker_lines.add(line.strip())

    # Write the lines from all_read_file that are not in linker_read_file to the output file
    with open(all_read_file, 'r') as all_read, open(nolinker_read_file, 'w') as nolinker_read:
        for line in all_read:
            if line.strip() not in linker_lines:
                nolinker_read.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_all_read.txt input_linker_read.txt output_nolinker_read.txt")
    else:
        all_read_file = sys.argv[1]
        linker_read_file = sys.argv[2]
        nolinker_read_file = sys.argv[3]

        extract_lines_not_in_linker(all_read_file, linker_read_file, nolinker_read_file)

