import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_species_mapping(species_file, species_mod_file):
    """
    Parses species mapping from two files.
    Returns a dictionary mapping original species names to modified names.
    """
    with open(species_file, "r") as f1, open(species_mod_file, "r") as f2:
        original_names = [line.strip() for line in f1.readlines()]
        modified_names = [line.strip() for line in f2.readlines()]

    if len(original_names) != len(modified_names):
        raise ValueError("Mismatch in number of species between the two files.")

    return {
        orig.replace(" ", "_"): mod.replace(" ", "_")  # Add underscores to names
        for orig, mod in zip(original_names, modified_names)
    }

def modify_phylip_species(phylip_file, species_mapping, output_file):
    """
    Modifies the PHYLIP file to replace species names using the provided mapping.
    Ensures that each species name and sequence are separated by at least two spaces.
    """
    with open(phylip_file, "r") as infile:
        lines = infile.readlines()

    # First line: number of sequences and alignment length
    header = lines[0]

    # Update sequence names
    updated_lines = [header]
    for line in lines[1:]:
        parts = line.split(maxsplit=1)  # Split only into name and sequence
        if len(parts) != 2:
            raise ValueError("Invalid PHYLIP line format.")
        species_name, sequence = parts
        updated_name = species_mapping.get(species_name, species_name)  # Use mapped name or keep original
        # Note the two (or more) spaces after updated_name
        updated_lines.append(f"{updated_name}  {sequence}")

    # Write updated PHYLIP file
    with open(output_file, "w") as outfile:
        outfile.writelines(updated_lines)

if __name__ == "__main__":
    # Set paths
    aligned_phylip = "back_translated_cds.fasta"
    species_file = "species.txt"
    species_mod_file = "speciesMOD.txt"
    modified_phylip = "modified_back_translated_alignment.phy"

    # Step 1: Parse species mapping
    species_mapping = parse_species_mapping(species_file, species_mod_file)

    # Step 2: Modify PHYLIP file
    modify_phylip_species(aligned_phylip, species_mapping, modified_phylip)

    print("Pipeline complete. The modified PHYLIP file (suitable for PAML) is:")
    print(f" - {modified_phylip}")

