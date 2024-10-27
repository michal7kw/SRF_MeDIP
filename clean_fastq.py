import gzip
import sys
from pathlib import Path

def clean_fastq(input_file, output_file):
    """
    Clean a FASTQ file by ensuring all sequence identifier lines start with '@'
    and maintaining the proper 4-line format.
    """
    with gzip.open(input_file, 'rt') as f_in, gzip.open(output_file, 'wt') as f_out:
        line_count = 0
        current_group = []
        
        for line in f_in:
            line_count += 1
            line = line.strip()
            
            # Start of a new sequence
            if line_count % 4 == 1:
                if not line.startswith('@'):
                    print(f"Warning: Found malformed header at line {line_count}: {line}")
                    line = '@' + line
            
            current_group.append(line)
            
            # Write complete sequence record
            if len(current_group) == 4:
                if (current_group[0].startswith('@') and 
                    current_group[2].startswith('+') and 
                    len(current_group[1]) == len(current_group[3])):
                    f_out.write('\n'.join(current_group) + '\n')
                else:
                    print(f"Skipping malformed record at line {line_count-3}")
                current_group = []

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = str(Path(input_file).with_suffix('')) + '_cleaned.fastq.gz'
    clean_fastq(input_file, output_file)
    print(f"Cleaned file saved as: {output_file}")

# python clean_fastq.py medip/Medip_N_input_r1.fastq.gz
# mv medip/Medip_N_input_r1_cleaned.fastq.gz medip/Medip_N_input_r1.fastq.gz