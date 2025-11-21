amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}



def calculate_mol_mass(peptide_seq, amino_acid_mass_dict=None):
    """
    Calculate the molecular mass of a single peptide sequence and
    return a dictionary {peptide_seq: mass}.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence (one-letter amino-acid codes).
    amino_acid_mass_dict : dict, optional
        Mapping from amino-acid letter to mass in Dalton. If None,
        amino_acid_mass_dalton is used.

    Returns
    -------
    dict
        {peptide_seq: total_mass_in_dalton}
    """
    if amino_acid_mass_dict is None:
        amino_acid_mass_dict = amino_acid_mass_dalton

    total_mass = 0.0
    for aa in peptide_seq:
        total_mass += amino_acid_mass_dict[aa]

    return {peptide_seq: total_mass}


def calculate_mol_mass_collection(peptides, amino_acid_mass_dict=None):
    """
    Calculate the molecular mass for a list of peptide sequences.

    Parameters
    ----------
    peptide_seqs : list of str
        Peptide sequences, e.g. ["MATSR", "PEPTIDE", ...]
    amino_acid_mass_dict : dict, optional
        Mapping {amino_acid_letter: mass_in_dalton}.
        If None, amino_acid_mass_dalton is used.

    Returns
    -------
    dict
        {peptide_sequence: total_mass_in_dalton}
    """
    if amino_acid_mass_dict is None:
        amino_acid_mass_dict = amino_acid_mass_dalton

    result = {}

    for peptide in peptides:
        total_mass = 0.0
        for aa in peptide:
    
            total_mass += amino_acid_mass_dict[aa]

        result[peptide] = total_mass

    return result


def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    Convert a peptide→mass dictionary to a peptide→m/z dictionary.

    Parameters
    ----------
    peptide_mass_dict : dict
        Mapping {peptide_sequence: neutral_mass_in_dalton}.
    charge : int, optional
        Charge state (z). Default is 2.
    proton_mass : float, optional
        Mass of a proton in Dalton. Default is 1.007.

    Returns
    -------
    dict
        {peptide_sequence: result}
    """
    mz_dict = {}

    for peptide, mass in peptide_mass_map.items():
        mz = (mass + charge * proton_mass) / charge
        mz_dict[peptide] = mz

    return mz_dict



import numpy as np
import matplotlib.pyplot as plt

def plot_ms1(mz_dict, random_count_range=(10, 100), seed=42):
    """
    Plot a simple mass spectrum from m/z values with random intensities
    drawn using NumPy.

    Parameters
    ----------
    mz_dict : dict
        {peptide_sequence: m_over_z_value}
    random_count_range : tuple (min_intensity, max_intensity)
    seed : int or None
        Random seed for reproducible intensities. If None, no seed is set.
    """
    if not isinstance(mz_dict, dict):
        raise TypeError(f"mz_dict must be a dict, got {type(mz_dict)}")

    # Optionally set the seed for reproducibility
    rng = np.random.default_rng(seed)

    # m/z values as a NumPy array (used on the x-axis)
    mz_values = np.array(list(mz_dict.values()), dtype=float)

    # Create random intensities
    min_count, max_count = random_count_range
    intensities = rng.integers(
        low=min_count,
        high=max_count + 1,   # upper bound is exclusive
        size=len(mz_values)
    )

    # Plot
    
    plt.bar(mz_values, intensities)
    plt.xlabel("m/z")
    plt.ylabel("Intensity (a.u.)")
    plt.title("Simulated MS1 Spectrum")
    
    plt.show()



def fragment_peptide(peptide):
    """
    Generate b- and y-ion fragment sequences for a peptide.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence, e.g. "PEPTIDE".

    Returns
    -------
    list of str
        Combined list of b- and y-ion sequences.
        First all b-ions (N-terminal), then all y-ions (C-terminal).
    """
    b_ions = []  # from N-terminus: P, PE, PEP, ...
    y_ions = []  # from C-terminus: E, DE, IDE, ...


    # i goes from 1 up to n (inclusive)
    for i in range(1, len(peptide)):
        # b-ion: take the first i amino acids
        b_ions.append(peptide[:i])

        # y-ion: take the last i amino acids
        y_ions.append(peptide[len(peptide) - i:]) # start to end

    # return one combined list: all b-ions followed by all y-ions
    return b_ions + y_ions

fragment_peptide("PEPT")