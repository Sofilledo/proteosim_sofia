def read_fasta(filepath):
    """
    Function converts text into a dictionary with seq and ids.

    Parameters
    ----------
     fasta_path : str or path-like
        Path to a FASTA file containing protein sequences. Header lines are
        expected to start with '>' and contain UniProt-style IDs separated
        by '|' (e.g. '>sp|P12345|Description').
    Returns
    -------
    dict of str to str
        Dictionary mapping protein IDs to their full amino-acid sequences
        (one continuous string per protein).
    """


    protein_map = {} # create an empty dictionary 
    current_id = None 
    current_sequence = [] #create an empty list 

    with open(filepath, 'r', encoding='utf-8') as fasta_handle: # opens the file and delets whitespaces
        for line in fasta_handle:
            stripped = line.strip()  
            if not stripped:
                continue
            if stripped.startswith('>'): # header line in FASTA
                if current_id is not None: #so we skip the first line as the first id is = 0 => None
                    protein_map[current_id] = ''.join(current_sequence) # 
                    current_sequence = []
                current_id = stripped.split('|')[1] # takes the second charachter and delets the bars
            else:
                current_sequence.append(stripped) #chooses the stripped part 

    if current_id is not None: # saves the previous protein
        protein_map[current_id] = ''.join(current_sequence)

    return protein_map