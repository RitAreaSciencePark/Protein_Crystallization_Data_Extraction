def clean_fasta_dots(input_fasta, output_fasta):
    """
    Remove leading/trailing '...' from FASTA sequences.
    """
    with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
        current_header = None
        sequence_parts = []

        def write_record():
            if current_header:
                clean_seq = "".join(sequence_parts)
                clean_seq = clean_seq.lstrip("-").rstrip("-")
                fout.write(current_header)
                fout.write(clean_seq + "\n")

        for line in fin:
            line = line.strip()

            if line.startswith(">"):
                write_record()
                current_header = line + "\n"
                sequence_parts = []
            elif line:
                sequence_parts.append(line)

        write_record()


# --------- usage ----------
input_fasta = "KaiB_MSA.fasta"
output_fasta = "cleaned_KaiB_MSA.fasta"

clean_fasta_dots(input_fasta, output_fasta)
