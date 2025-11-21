
from proteosim.protein_digestion import digest_protein_collection
from proteosim.protein_digestion import compute_sequence_coverage
import re

def test_digest_protein_collection():
    dummy_proteins = {
        "P1": "AKRPQ"
    }

    result = digest_protein_collection(
        protein_map=dummy_proteins,
        enzyme_name="Trypsin",
        min_pep_len=1,
        max_pep_len=10,
    )

    assert result["P1"] == ["AK", "RPQ"]


def test_compute_sequence_coverage():
    dummy_prot_seq = "AAAAA"
    dummy_peps = ["AAA"]

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)

    # "AAA" covers 3 out of 5 positions -> 60%
    assert coverage == 60.0