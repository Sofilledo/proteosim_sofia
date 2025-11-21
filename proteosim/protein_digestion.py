import re
enzyme_cleavage_patterns = {
        'LysC': r'(?<=K)',
        'LysN': r'(?=K)',
        'ArgC': r'(?<=R)',
        'Trypsin': r'(?<=[KR])(?!P)', 
    }




def digest_protein_sequence(protein_seq, cleave_pattern): # defines the new function
    """
    Function description
    Digest a protein sequence into peptides using a regex cleavage pattern.
    Parameters
    ----------
    protein_seq : str
        Amino-acid sequence of the protein (one-letter code).
    cleave_pattern : str
        Regular expression describing the cleavage sites, e.g.
        enzyme_cleavage_patterns['LysC'].
    Returns
    -------
    list of str
        List of peptide sequences produced by the in-silico digestion.
    """
    peptides = re.split(protein_seq, cleave_pattern)
    
    enzyme_cleavage_patterns = {
        'LysC': r'(?<=K)',
        'LysN': r'(?=K)',
        'ArgC': r'(?<=R)',
        'Trypsin': r'(?<=[KR])(?!P)', 
    }
    #sequence = protein_seq # aqui he cambiado algo

    cleave_pattern = enzyme_cleavage_patterns[cleave_pattern]


    print(peptides)
    print(f'Nr. of digested peptides: {len(peptides)}')



def digest_protein_collection(protein_map, enzyme_name, min_pep_len=5, max_pep_len=30):
    """
    Digest a collection of protein sequences into peptides using a regex cleavage pattern. 

    Parameters
    ----------
    protein_map : dict
        Dictionary mapping protein IDs to amino-acid sequences (one-letter code),
        e.g. {"P11802": "MKK....", ...}
    enzyme_name : str
        Name of the enzyme, one of: 'LysC', 'LysN', 'ArgC', 'Trypsin'.
    min_pep_len : int
        Minimum allowed peptide length (in amino acids).
    max_pep_len : int
        Maximum allowed peptide length (in amino acids).

    Returns
    -------
    dict
        Dictionary mapping protein IDs to a list of digested peptides that pass
        the length filters.
    """

    enzyme_cleavage_patterns = {
        'LysC':   r'(?<=K)',
        'LysN':   r'(?=K)',
        'ArgC':   r'(?<=R)',
        'Trypsin': r'(?<=[KR])(?!P)',   # cut after K/R, but not if followed by P
    }

    # Look up the regex pattern for the chosen enzyme
    pattern = enzyme_cleavage_patterns[enzyme_name]

    # This will store the final result
    digested_dict = {}

    # Loop over all proteins in the collection
    for protein_id, sequence in protein_map.items():
        # Split the sequence at the cleavage pattern
        peptides = re.split(pattern, sequence)

        # Filter peptides by length (keep only those in the desired range)
        filtered_peptides = [
            pep for pep in peptides
            if min_pep_len <= len(pep) <= max_pep_len
        ]

        if not filtered_peptides:
            print(f"{protein_id}: does not fit criteria")
            # go on to the next protein instead of stopping the whole function
            continue

        print(f"{protein_id}: {filtered_peptides}")
        print(f'Nr. of digested peptides: {len(filtered_peptides)}')

        # Save peptides for this protein in the output dictionary
        digested_dict[protein_id] = filtered_peptides

    # Only return once, after the loop is done
    return digested_dict



def compute_sequence_coverage(protein_seq, peptides):
    """
    Add a short description here.

    Parameters
    ----------
    protein_seq = based on the orginal sequence from FASTA
    peptides = digested using function 
    Returns
    -------
    """
    if not protein_seq:
        return 0.0

    covered = [False] * len(protein_seq)

    for pep in peptides:
        if not pep:
            continue

        # find first occurrence of the peptide in the protein
        idx = protein_seq.find(pep)
        if idx == -1:
            continue

        # mark covered positions
        for i in range(idx, idx + len(pep)):
            if i < len(protein_seq):
                covered[i] = True

    return sum(covered) / len(protein_seq) * 100.0
    