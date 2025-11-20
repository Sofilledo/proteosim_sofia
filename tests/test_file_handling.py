from proteosim.file_handling import read_fasta

def test_read_fasta():
    # path is relative to the tests/ folder
    tmp_fasta_path = "data/dummy_fasta.fasta"
    
    # call the function you want to test
    protein_map = read_fasta(tmp_fasta_path)

    # now check the content
    assert protein_map["P11802"] == "DUMMY"
    assert protein_map["A0A087WTH1"] == "DUMMY"
