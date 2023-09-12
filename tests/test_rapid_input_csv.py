"""
Test module for rapid_utils/rapid_input_csv.py.
"""
import os
import numpy as np
from numpy.testing import assert_array_equal
from rapid_utils import rapid_input_csv

TEST_DIR = os.path.relpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(TEST_DIR, 'data')
OUTPUT_DIR = os.path.join(TEST_DIR, 'output')


def test_write_connectivity_file():
    """
    Verify that write_connectivity_file produces a connectivity CSV consistent
    with benchmark.
    """

    flowline_file = os.path.join(DATA_DIR, 'nhdplusv2_sample',
                                 'nhdplusv2_sample.shp')
    flowline_id_field_name = 'COMID'
    downstream_id_field_name = 'DWNCOMID'
    out_csv_file = os.path.join(OUTPUT_DIR,
                                'nhdplusv2_sample_rapid_connect_test.csv')

    rapid_input_csv.write_connectivity_file(flowline_file,
                                            flowline_id_field_name,
                                            downstream_id_field_name,
                                            out_csv_file)

    connectivity = np.genfromtxt(out_csv_file, delimiter=',', dtype=int)

    expected = [[3820329, 0, 2, 3820227, 3820229],
                [3820227, 3820329, 2, 3820161, 3820163],
                [3820161, 3820227, 0, 0, 0],
                [3820163, 3820227, 0, 0, 0],
                [3820229, 3820329, 0, 0, 0]]

    assert_array_equal(connectivity, expected)


def test_write_riv_bas_id_file():
    """
    Verify that write_riv_bas_id_file produces a basin CSV consistent with
    benchmark.
    """

    flowline_file = os.path.join(DATA_DIR, 'nhdplusv2_sample',
                                 'nhdplusv2_sample.shp')
    flowline_id_field_name = 'COMID'
    out_csv_file = os.path.join(OUTPUT_DIR,
                                'nhdplusv2_sample_riv_bas_id_test.csv')

    rapid_input_csv.write_riv_bas_id_file(flowline_file,
                                          flowline_id_field_name,
                                          out_csv_file)

    riv_bas_id = np.genfromtxt(out_csv_file, delimiter=',', dtype=int)

    expected = [3820329, 3820227, 3820161, 3820163, 3820229]

    assert_array_equal(riv_bas_id, expected)


def test_write_kfac_file():
    """
    Verify that write_kfac_file produces a kfac CSV consistent with benchmark.
    """

    flowline_file = os.path.join(DATA_DIR, 'nhdplusv2_sample',
                                 'nhdplusv2_sample.shp')
    connectivity_file = os.path.join(DATA_DIR, 'nhdplusv2_sample',
                                     'nhdplusv2_sample_rapid_connect.csv')
    out_csv_file = os.path.join(OUTPUT_DIR,
                                'nhdplusv2_sample_kfac_test.csv')
    flowline_id_field_name = 'COMID'
    length_field_name = 'LENGTHKM'
    slope_field_name = 'SLOPE'
    formula_type = 3
    input_length_units = 'km'

    rapid_input_csv.write_kfac_file(flowline_file, connectivity_file,
                                    out_csv_file, flowline_id_field_name,
                                    length_field_name, slope_field_name,
                                    formula_type=formula_type,
                                    input_length_units=input_length_units)

    kfac = np.genfromtxt(out_csv_file, delimiter=',')

    expected = [7934.5, 4338.8, 6032.0, 8023.5, 4534.1]

    assert_array_equal(kfac, expected)


def test_write_constant_x_file():
    """
    Verify that write_constant_x_file produces a constant Muskingum x CSV
    consistent with benchmark.
    """

    connectivity_file = os.path.join(DATA_DIR, 'nhdplusv2_sample',
                                     'nhdplusv2_sample_rapid_connect.csv')
    out_csv_file = os.path.join(OUTPUT_DIR,
                                'nhdplusv2_sample_x_test.csv')

    rapid_input_csv.write_constant_x_file(connectivity_file,
                                          out_csv_file,
                                          0.3)

    x_param = np.genfromtxt(out_csv_file, delimiter=',')

    expected = [0.3, 0.3, 0.3, 0.3, 0.3]

    assert_array_equal(x_param, expected)
