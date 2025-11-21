from proteosim.mass_spectra_simulation import calculate_mol_mass_collection
from proteosim.mass_spectra_simulation import calculate_mz_collection
from proteosim.mass_spectra_simulation import fragment_peptide



def test_calculate_mol_mass_collection():
    peptides = ['MATSR',
                'TYLDK']
    expected = {'MATSR': 546.65,
                'TYLDK': 620.71}
    actual = calculate_mol_mass_collection(peptides)
    
    print(actual)
    print(expected)
    assert actual == expected

test_calculate_mol_mass_collection()

def test_calculate_mz_collection():
    # input: peptide -> neutral mass
    peptide_mass_map = {'MATSR': 546.65}

    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {'MATSR': 274.332}

    assert actual == expected


test_calculate_mz_collection()

def test_fragment_peptide():
    peptide = "PEPT"
    expected = {'P', 'PE', 'PEP', 'T', 'PT', 'EPT'}
    actual = set(fragment_peptide(peptide))

    assert actual == expected

test_fragment_peptide()