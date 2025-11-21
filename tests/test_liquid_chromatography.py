from proteosim.liquid_chromatography import predict_lc_retention_times
from proteosim.liquid_chromatography import select_retention_time_window
from pyteomics import achrom



def test_predict_lc_retention_times():
    peptides = ["KK"]
    expected = achrom.calculate_RT("KK", achrom.RCs_guo_ph7_0)

    actual = predict_lc_retention_times(peptides)

    # actual is a dict {"KK": some_value}
    assert actual == {"KK": expected}

test_predict_lc_retention_times()




def test_select_retention_time_window():
    peptide_rt_map = {
        "PEP1": 5.0,
        "PEP2": 10.0,
        "PEP3": 15.0,
    }

    selected = select_retention_time_window(
        peptide_rt_map,
        lower_ret_time=6.0,
        upper_ret_time=14.0,
    )

    # Only PEP2 (RT = 10.0) is within [6.0, 14.0]
    assert selected == {"PEP2": 10.0}


test_select_retention_time_window()
