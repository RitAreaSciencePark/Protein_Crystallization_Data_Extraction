#!/usr/bin/env python3

def read_and_clean_fasta_file(fasta_file):
    """
    Reads a FASTA file, removes all '-' from sequences,
    and returns a dictionary {query_id: sequence}.
    Works with any type of header.
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    # Save previous sequence
                    fasta_dict[header] = ''.join(seq_lines).replace('-', '')
                header = line[1:].split()[0]  # Take entire header or first word as ID
                seq_lines = []
            else:
                seq_lines.append(line)
        # Save last sequence
        if header:
            fasta_dict[header] = ''.join(seq_lines).replace('-', '')
    return fasta_dict


def write_fasta_dict(fasta_dict, output_file):
    """Writes a dictionary {header: seq} to a FASTA file"""
    with open(output_file, 'w') as f:
        for header, seq in fasta_dict.items():
            f.write(f">{header}\n")
            # wrap sequence at 80 chars
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


