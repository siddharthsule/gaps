#!/usr/bin/env python3

# fix any zero values in the first two columns of a YODA file
import argparse
import os


def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip lines that are not data lines
            if line.startswith('#') or line.startswith('BEGIN') or \
                line.startswith('Path') or line.startswith('ScaledBy') or \
                line.startswith('Title') or line.startswith('Type') or \
                    line.startswith('---') or line.startswith('END') or \
                    line.startswith('Total') or line.startswith('Error'):
                outfile.write(line)
                continue

            # Split the line into columns
            columns = line.split()

            # Check and modify the first and second columns if necessary
            if len(columns) >= 2:
                try:
                    first_col = float(columns[0])
                    second_col = float(columns[1])
                    if abs(first_col) < 1e-15:
                        columns[0] = '0.000000e+00'
                    if abs(second_col) < 1e-15:
                        columns[1] = '0.000000e+00'
                except ValueError:
                    # If conversion to float fails, just write the line as is
                    pass

            # Write the modified line to the output file
            outfile.write('\t'.join(columns) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Fix any zero values in the first two columns of a YODA file.')
    parser.add_argument('input_file', help='The input YODA file to process')
    args = parser.parse_args()

    input_file = args.input_file
    temp_output_file = input_file + '.tmp'

    # Process the file and write to a temporary file
    process_file(input_file, temp_output_file)

    # Replace the original file with the modified file
    os.replace(temp_output_file, input_file)


if __name__ == '__main__':
    main()
