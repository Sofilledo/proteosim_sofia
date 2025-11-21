from pyteomics import achrom

def predict_lc_retention_times(peptides): 
    """
    peptides: list of peptide sequences, e.g. ["MATSR", "PEPTIDE", ...]
    Returns: dict {peptide_sequence: retention_time}
    """
    relative_RT = {
        p: achrom.calculate_RT(p, achrom.RCs_guo_ph7_0)
        for p in peptides  # p is e.g. "MATSR"
    }
    print(relative_RT)
    return relative_RT

import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Add a short description here.

    Parameters
    ----------

    Returns
    -------

    """
    if retention_times == dict : 
        retention_times = list()
    
    else : 

        plt.hist(retention_times, resolution)   # you can change bins
        plt.xlabel("")
        plt.ylabel("a.u(retention time)")
        plt.title("Histogram of retention times")
        plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    peptide_rt_map : dict
        {peptide_sequence: retention_time}
    lower_ret_time, upper_ret_time : float
        Allowed RT range.
    """
    filtered = {}  # will store only peptides within the RT range

    for pep, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            filtered[pep] = rt
        
    return filtered
    