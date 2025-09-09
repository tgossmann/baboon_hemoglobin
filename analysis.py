import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def read_fasta_files(folder_path):
    """Reads all fasta files in a folder and renames headers with filenames."""
    combined_records = []

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fasta") or file_name.endswith(".fa") or file_name.endswith(".fna"):
            file_path = os.path.join(folder_path, file_name)
            
            for record in SeqIO.parse(file_path, "fasta"):
                # Rename the header with the filename (without extension)
                new_header = os.path.splitext(file_name)[0]
                record.id = new_header
                record.description = ""
                combined_records.append(record)

    return combined_records

def translate_sequences(records):
    """Translates CDS sequences into protein sequences using the standard genetic code."""
    translated_records = []
    for record in records:
        protein_seq = record.seq.translate(to_stop=True)  # Translate until stop codon
        translated_record = SeqRecord(protein_seq, id=record.id, description="Translated")
        translated_records.append(translated_record)
    return translated_records

def write_fasta_file(records, output_file):
    """Writes a list of SeqRecord objects to a fasta file."""
    with open(output_file, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")

def align_proteins_with_prank(input_file, output_file):
    """Aligns protein sequences using PRANK."""
    prank_command = [
        "prank", f"-d={input_file}", f"-o={output_file}", "-f=fa"
    ]
    print(prank_command)
    subprocess.run(prank_command, check=True)

def normalize_header(header):
    """Normalizes headers for comparison by replacing spaces with underscores and removing '_Translated'."""
    return header.replace(" ", "_").replace("_Translated", "")

def back_translate_alignment(aligned_protein_file, original_cds_records, output_file):
    """
    Back-translates the aligned protein sequences to nucleotide sequences, keeping gaps.
    Outputs the alignment in long PHYLIP format suitable for PAML analysis.
    """
    protein_alignment = list(SeqIO.parse(aligned_protein_file, "fasta"))
    back_translated_records = []

    for prot_record in protein_alignment:
        # Normalize protein header and find the matching CDS record
        normalized_id = normalize_header(prot_record.id)
        matching_cds = next((cds for cds in original_cds_records if normalize_header(cds.id) == normalized_id), None)
        if not matching_cds:
            raise ValueError(f"No matching CDS found for protein ID {prot_record.id}")

        back_translated_seq = []
        cds_index = 0  # Track position in the CDS sequence

        # Iterate through protein alignment and construct back-translated nucleotide sequence
        for aa in prot_record.seq:
            if aa == "-":
                # Add gap for alignment
                back_translated_seq.append("---")
            else:
                # Add corresponding codon from CDS
                back_translated_seq.append(str(matching_cds.seq[cds_index:cds_index + 3]))
                cds_index += 3

        # Join the codons to form the full nucleotide sequence
        back_translated_seq = "".join(back_translated_seq)

        # Replace spaces with underscores in sequence IDs
        back_translated_record = SeqRecord(
            Seq(back_translated_seq),
            id=matching_cds.id.replace(" ", "_"),
            description=""
        )
        back_translated_records.append(back_translated_record)

    # Write the back-translated sequences to a PHYLIP file in long format
    write_phylip_file_long(back_translated_records, output_file)

def write_phylip_file_long(records, output_file):
    """
    Writes a list of SeqRecord objects to a PHYLIP file in long format.
    """
    num_sequences = len(records)
    alignment_length = len(records[0].seq)

    with open(output_file, "w") as out_f:
        # Write the header line with the number of sequences and alignment length
        out_f.write(f"{num_sequences} {alignment_length}\n")

        # Write each sequence in long PHYLIP format
        for record in records:
            # Format: Full ID followed by the sequence
            out_f.write(f"{record.id} {record.seq}\n")



if __name__ == "__main__":
    # Set paths
    input_folder = "Seqs2"  # Replace with your folder path
    output_cds_file = "combined_cds.fasta"
    output_protein_file = "translated_proteins.fasta"
    aligned_protein_file = "aligned_proteins.fasta"
    back_translated_file = "back_translated_cds.fasta"

    # Step 1: Read and rename fasta files
    cds_records = read_fasta_files(input_folder)
    write_fasta_file(cds_records, output_cds_file)

    # Step 2: Translate sequences to proteins
    protein_records = translate_sequences(cds_records)
    write_fasta_file(protein_records, output_protein_file)

    # Step 3: Align protein sequences with PRANK
    align_proteins_with_prank(output_protein_file, aligned_protein_file)

    # Step 4: Back-translate the alignment
    back_translate_alignment(aligned_protein_file+'.best.fas', cds_records, back_translated_file)

    print("Pipeline complete. Outputs:")
    print(f" - Combined CDS: {output_cds_file}")
    print(f" - Translated Proteins: {output_protein_file}")
    print(f" - Aligned Proteins: {aligned_protein_file}")
    print(f" - Back-translated CDS: {back_translated_file}")
