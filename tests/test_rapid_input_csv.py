"""
Test module for rapid_utils/rapid_input_csv.py.
"""
import os
from rapid_utils import rapid_input_csv
import numpy as np
from numpy.testing import assert_almost_equal

TEST_DIR = os.path.relpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(TEST_DIR, 'data')
OUTPUT_DIR = os.path.join(TEST_DIR, 'output')

def read_csv_file(fname, skip_header=0):
    """
    Read a CSV file into a NumPy array.

    Parameters
    ----------
    fname : str
        Name of file.
    skip_header : int, optional
        Number of lines to skip at the beginning of the file.

    Returns
    -------
    data : ndarray
        Contents of CSV as an array.
    """

    data = np.genfromtxt(fname, delimiter=',', skip_header=skip_header)

    return data

def compare_csv_files(output_file, benchmark_file, skip_header=0):
    """
    Check if the contents of two weight-table files match.

    Parameters
    ----------
    output_file : str
        Name of file to compare.
    benchmark_file : str
        Name of file to compare.
    skip_header : int, optional
        Number of lines to skip at the beginning of each file.
    """

    output = read_csv_file(output_file, skip_header=skip_header)
    benchmark = read_csv_file(benchmark_file, skip_header=skip_header)
    assert_almost_equal(output, benchmark)

def test_write_connectivity_file():
    """
    Verify that write_connectivity_file produces a connectivity CSV consistent
    with benchmark.
    """

    flowline_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                               'NHDPlusV2_sample.shp')
    flowline_id_field_name='COMID'
    downstream_id_field_name='DWNCOMID'
    out_csv_file=os.path.join(OUTPUT_DIR,
                              'NHDPlusV2_sample_rapid_connect_test.csv')

    rapid_input_csv.write_connectivity_file(flowline_file,
                                            flowline_id_field_name,
                                            downstream_id_field_name,
                                            out_csv_file)

    benchmark_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                'NHDPlusV2_sample_rapid_connect.csv')

    compare_csv_files(out_csv_file, benchmark_file, skip_header=0)

def test_write_riv_bas_id_file():
    """
    Verify that write_riv_bas_id_file produces a basin CSV consistent with
    benchmark.
    """

    flowline_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                               'NHDPlusV2_sample.shp')
    flowline_id_field_name='COMID'
    out_csv_file=os.path.join(OUTPUT_DIR,
                              'NHDPlusV2_sample_riv_bas_id_test.csv')

    rapid_input_csv.write_riv_bas_id_file(flowline_file,
                                          flowline_id_field_name,
                                          out_csv_file)

    benchmark_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                'NHDPlusV2_sample_riv_bas_id.csv')

    compare_csv_files(out_csv_file, benchmark_file, skip_header=0)

def test_write_kfac_file():
    """
    Verify that write_kfac_file produces a kfac CSV consistent with benchmark.
    """

    flowline_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                               'NHDPlusV2_sample.shp')
    connectivity_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                   'NHDPlusV2_sample_rapid_connect.csv')
    out_csv_file=os.path.join(OUTPUT_DIR,
                              'NHDPlusV2_sample_kfac_test.csv')
    flowline_id_field_name='COMID'
    length_field_name='LENGTHKM'
    slope_field_name='SLOPE'
    formula_type=3
    input_length_units='km'

    rapid_input_csv.write_kfac_file(flowline_file, connectivity_file,
                                    out_csv_file, flowline_id_field_name,
                                    length_field_name, slope_field_name,
                                    formula_type, input_length_units)

    benchmark_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                'NHDPlusV2_sample_kfac.csv')

    compare_csv_files(out_csv_file, benchmark_file, skip_header=0)

def test_write_constant_x_file():
    """
    Verify that write_constant_x_file produces a constant Muskingum x CSV
    consistent with benchmark.
    """

    connectivity_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                   'NHDPlusV2_sample_rapid_connect.csv')
    out_csv_file=os.path.join(OUTPUT_DIR,
                              'NHDPlusV2_sample_x_test.csv')

    rapid_input_csv.write_constant_x_file(connectivity_file,
                                          out_csv_file,
                                          0.3)

    benchmark_file=os.path.join(DATA_DIR,'NHDPlusV2_sample',
                                'NHDPlusV2_sample_x.csv')

    compare_csv_files(out_csv_file, benchmark_file, skip_header=0)
