"""Test module for rapid_utils/rapid_input_discharge.py."""

import os
import numpy as np
from netCDF4 import Dataset
from rapid_utils.rapid_input_discharge import RAPIDInputDischarge

TEST_DIR = os.path.relpath(os.path.dirname(__file__))
DATADIR = os.path.join(TEST_DIR, 'data')
OUTPUTDIR = os.path.join(TEST_DIR, 'output')

def test_write_nc():
    output_filename = os.path.join(
        OUTPUTDIR, 'rapid_input_discharge_test_write.nc')

    a = RAPIDInputDischarge(output_filename=output_filename)
    a.write_nc()

    d = Dataset(output_filename)

    keys = d.variables.keys()
    expected = ('Qout', 'rivid', 'time', 'lon', 'lat', 'crs')

    assert set(keys) == set(expected)

def test_parse_input_discharge_file():
    input_discharge_file = os.path.join(DATADIR, 'qout.nc')

    a = RAPIDInputDischarge(input_discharge_file=input_discharge_file)
    a.parse_input_discharge_file()

    discharge = a.input_discharge[:5]
    expected_discharge = [
        0.07708612, 0.02149609, 0.04948052, 0.91938674, 0.02969897]

    np.testing.assert_allclose(discharge, expected_discharge, atol=1e-8, rtol=0)

    rivid = a.input_discharge_rivid[:5]
    expected_rivid = [70563, 70564, 70618, 70625, 70626]

    assert np.array_equal(rivid, expected_rivid)

def test_parse_rivid_file():
    rivid_file = os.path.join(DATADIR, 'rapid_connect.csv')
    a = RAPIDInputDischarge(rivid_file=rivid_file)
    a.parse_rivid_file()

    rivid = a.rivid[:5]
    expected_rivid = [70563, 70564, 70618, 70625, 70626]

    assert np.array_equal(rivid, expected_rivid)

def test_sort_discharge_by_rivid():
    input_discharge_file = os.path.join(DATADIR, 'qout.nc')

    rivid = [70625, 70652]
    
    a = RAPIDInputDischarge(input_discharge_file=input_discharge_file,
                            rivid=rivid)
    a.parse_input_discharge_file()
    a.sort_discharge_by_rivid()

    discharge = a.discharge
    expected = [0.91938674, 0.5366134]

    np.testing.assert_allclose(discharge, expected, atol=1e-8, rtol=0)

def test_integrate_over_files_mean():
    input_file_list = [os.path.join(DATADIR, 'qout_1.nc'),
                       os.path.join(DATADIR, 'qout_2.nc')]

    a = RAPIDInputDischarge(input_file_list=input_file_list)
    a.integrate_over_files(integration_type='mean')

    input_discharge = a.input_discharge
    expected = [564.20994494, 703.41273011, 659.68159059]
    np.testing.assert_allclose(input_discharge, expected, atol=1e-8, rtol=0)

def test_main():
    rivid = [70625, 70652]
    input_discharge_file = os.path.join(DATADIR, 'qout.nc')

    output_filename = os.path.join(
        OUTPUTDIR, 'rapid_input_discharge_test_main.nc')

    a = RAPIDInputDischarge(input_discharge_file=input_discharge_file,
                            rivid=rivid, output_filename=output_filename)
    a.main()

    discharge = a.discharge
    expected = [0.91938674,  0.5366134]

    np.testing.assert_allclose(discharge, expected, atol=1e-8, rtol=0)

if __name__ == '__main__':
    # test_write_nc()
    # test_parse_input_discharge_file()
    # test_parse_rivid_file()
    # test_sort_discharge_by_rivid()
    # test_integrate_over_files_mean()
    test_main()
    
