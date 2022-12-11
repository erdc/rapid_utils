"""Test module for rapid_utils/rapid_namelist.py."""

import os
import rapid_utils.rapid_namelist as rn


TEST_DIR = os.path.relpath(os.path.dirname(__file__))
DATADIR = os.path.join(TEST_DIR, 'data')
OUTPUTDIR = os.path.join(TEST_DIR, 'output')

def test_formatted_warning():
    """Verify that a warning message is formatted correctly."""
    message = 'message'
    category = UserWarning
    filename = 'filename'
    lineno = 'lineno'
    f = None
    line = None

    out = rn.formatted_warning(message, category, filename, lineno,
                               file=f, line=line)

    expected = 'filename:lineno: UserWarning: message\n'

    assert out == expected

def test_read_namelist():
    """Verify that read_namelist correctly parses parameters from an
    existing namelist file.
    """
    input_namelist_file = os.path.join(DATADIR, 'rapid_namelist')

    a = rn.RAPIDNamelist()
    parsed = a.read_namelist(filename=input_namelist_file)

    expected = {'BS_opt_Qfinal': '.FALSE.',
                'BS_opt_Qinit': '.FALSE.',
                'BS_opt_dam': '.FALSE.',
                'BS_opt_for': '.FALSE.',
                'BS_opt_influence': '.FALSE.',
                'BS_opt_transpose_qout': '.FALSE.',
                'IS_dam_tot': '0',
                'IS_dam_use': '0',
                'IS_for_tot': '0',
                'IS_for_use': '0',
                'IS_max_up': '0',
                'IS_obs_tot': '0',
                'IS_obs_use': '0',
                'IS_opt_phi': '1',
                'IS_opt_routing': '1',
                'IS_opt_run': '1',
                'IS_riv_bas': '0',
                'IS_riv_tot': '0',
                'IS_strt_opt': '0',
                'Qfinal_file': "''",
                'Qfor_file': "''",
                'Qinit_file': "''",
                'Qobs_file': "''",
                'Qobsbarrec_file': "''",
                'QoutRabsmax_file': "''",
                'QoutRabsmin_file': "''",
                'Qout_file': "''",
                'Vlat_file': "''",
                'ZS_TauM': '172800',
                'ZS_TauO': '0',
                'ZS_TauR': '10800',
                'ZS_dtF': '0',
                'ZS_dtM': '10800',
                'ZS_dtO': '0',
                'ZS_dtR': '900',
                'ZS_knorm_init': '0',
                'ZS_phifac': '0',
                'ZS_xnorm_init': '0',
                'babsmax_file': "''",
                'dam_tot_id_file': "''",
                'dam_use_id_file': "''",
                'for_tot_id_file': "''",
                'for_use_id_file': "''",
                'k_file': 'input/k.csv',
                'kfac_file': "''",
                'obs_tot_id_file': "''",
                'obs_use_id_file': "''",
                'rapid_connect_file':
                'input/rapid_connect.csv',
                'riv_bas_id_file': 'input/riv_bas_id.csv',
                'x_file': 'input/x.csv',
                'xfac_file': "''"}

    assert parsed == expected

def test_write_default_namelist():
    """Verify that RAPIDNamelist can write a namelist with the correct
    default parameters.
    """
    a = rn.RAPIDNamelist()
    output_filename = os.path.join(OUTPUTDIR, 'rapid_namelist_default.txt')
    a.write_namelist(filename=output_filename)

    parsed = a.read_namelist(filename=output_filename)

    expected = {'BS_opt_Qfinal': '.FALSE.',
                'BS_opt_Qinit': '.FALSE.',
                'BS_opt_dam': '.FALSE.',
                'BS_opt_for': '.FALSE.',
                'BS_opt_influence': '.FALSE.',
                'BS_opt_transpose_qout': '.FALSE.',
                'IS_dam_tot': '0',
                'IS_dam_use': '0',
                'IS_for_tot': '0',
                'IS_for_use': '0',
                'IS_max_up': '2',
                'IS_obs_tot': '0',
                'IS_obs_use': '0',
                'IS_opt_phi': '1',
                'IS_opt_routing': '1',
                'IS_opt_run': '1',
                'IS_riv_bas': '0',
                'IS_riv_tot': '0',
                'IS_strt_opt': '0',
                'Qfinal_file': "''",
                'Qfor_file': "''",
                'Qinit_file': "''",
                'Qobs_file': "''",
                'Qobsbarrec_file': "''",
                'QoutRabsmax_file': "''",
                'QoutRabsmin_file': "''",
                'Qout_file': "''",
                'Vlat_file': "''",
                'ZS_TauM': '0',
                'ZS_TauO': '0',
                'ZS_TauR': '0',
                'ZS_dtF': '0',
                'ZS_dtM': '86400',
                'ZS_dtO': '0',
                'ZS_dtR': '900',
                'ZS_knorm_init': '0',
                'ZS_phifac': '0',
                'ZS_xnorm_init': '0',
                'babsmax_file': "''",
                'dam_tot_id_file': "''",
                'dam_use_id_file': "''",
                'for_tot_id_file': "''",
                'for_use_id_file': "''",
                'k_file': 'input/k.csv',
                'kfac_file': "''",
                'obs_tot_id_file': "''",
                'obs_use_id_file': "''",
                'rapid_connect_file': 'input/rapid_connect.csv',
                'riv_bas_id_file': 'input/riv_bas_id.csv',
                'x_file': 'input/x.csv',
                'xfac_file': "''"}

    assert parsed == expected

def test_update_params():
    """Verify that update_params correctly updates the `params` attribute
    from a specified dictionary.
    """
    new_params = {'ZS_TauM': 172800, 'ZS_TauR': 10800}

    a = rn.RAPIDNamelist()
    a.update_params(new_params)

    for new_key, new_value in new_params.items():
        key = new_key.lower()
        value = a.params[key]['value']
        assert value == new_value

def test_parse_riv_bas_id_file():
    """Verify that parse_riv_bas_id_file returns a dictionary with the
    correct value for `IS_riv_bas`.
    """
    filename = os.path.join(DATADIR, 'riv_bas_id.csv')

    a = rn.RAPIDNamelist()
    a.update_params({'riv_bas_id_file': filename})

    parsed = a.parse_riv_bas_id_file()

    expected = {'IS_riv_bas': 50}

    assert parsed == expected

def test_parse_connectivity_file():
    """Verify that parse_connectivity_file returns a dictionary with the
    correct values for `IS_riv_tot` and `IS_max_up`.
    """
    filename = os.path.join(DATADIR, 'rapid_connect.csv')

    a = rn.RAPIDNamelist()
    a.update_params({'rapid_connect_file': filename})

    parsed = a.parse_connectivity_file()

    expected = {'IS_riv_tot': 50,
                'IS_max_up': 2}

    assert parsed == expected

def test_parse_forcing_tot_file():
    """Verify that parse_forcing_tot_file returns a dictionary with the
    correct value for `IS_for_tot`.
    """
    filename = os.path.join(DATADIR, 'for_tot_id.csv')

    a = rn.RAPIDNamelist()
    a.update_params({'for_tot_id_file': filename})

    parsed = a.parse_forcing_tot_file()

    expected = {'IS_for_tot': 4}

    assert parsed == expected

def test_parse_forcing_use_file():
    """Verify that parse_forcing_use_file returns a dictionary with the
    correct value for `IS_for_use`.
    """
    filename = os.path.join(DATADIR, 'for_use_id.csv')

    a = rn.RAPIDNamelist()
    a.update_params({'for_use_id_file': filename})

    parsed = a.parse_forcing_use_file()

    expected = {'IS_for_use': 3}

    assert parsed == expected

def test_parse_vlat_file():
    """Verify that parse_vlat_file returns the correct values for
    `ZS_TauM` and `ZS_TauR`.
    """
    filename = os.path.join(DATADIR, 'inflow_lis.nc')

    a = rn.RAPIDNamelist()
    a.update_params({'Vlat_file': filename})

    parsed = a.parse_vlat_file()
    print(parsed)

    expected = {'ZS_TauM': 86400,
                'ZS_TauR': 10800}

    assert parsed == expected

def test_main():
    """Verify that the main method of RAPIDNamelist writes a valid
    namelist file using default parameters and parameters parsed from
    specified input files.
    """
    output_filename = os.path.join(OUTPUTDIR, 'rapid_namelist')
    input_file_dict = {'riv_bas_id_file':
                       os.path.join(DATADIR, 'riv_bas_id.csv'),
                       'for_tot_id_file':
                       os.path.join(DATADIR, 'for_tot_id.csv'),
                       'for_use_id_file':
                       os.path.join(DATADIR, 'for_use_id.csv'),
                       'rapid_connect_file':
                       os.path.join(DATADIR, 'rapid_connect.csv'),
                       'Vlat_file':
                       os.path.join(DATADIR, 'inflow_lis.nc')}

    a = rn.RAPIDNamelist(output_filename=output_filename, **input_file_dict)
    a.main()

    b = rn.RAPIDNamelist()
    parsed = b.read_namelist(output_filename)
    b.update_params(parsed, purge=True)

    params = {param_dict['rapid_name']:param_dict['value'] for param_dict
              in b.params.values()}

    expected = {'BS_opt_Qfinal': False,
                'BS_opt_Qinit': False,
                'BS_opt_dam': False,
                'BS_opt_for': False,
                'BS_opt_influence': False,
                'BS_opt_transpose_qout': False,
                'IS_dam_tot': 0,
                'IS_dam_use': 0,
                'IS_for_tot': 4,
                'IS_for_use': 3,
                'IS_max_up': 2,
                'IS_obs_tot': 0,
                'IS_obs_use': 0,
                'IS_opt_phi': 1,
                'IS_opt_routing': 1,
                'IS_opt_run': 1,
                'IS_riv_bas': 50,
                'IS_riv_tot': 50,
                'IS_strt_opt': 0,
                'Qfinal_file': None,
                'Qfor_file': None,
                'Qinit_file': None,
                'Qobs_file': None,
                'Qobsbarrec_file': None,
                'QoutRabsmax_file': None,
                'QoutRabsmin_file': None,
                'Qout_file': None,
                'Vlat_file': 'tests/data/inflow_lis.nc',
                'ZS_TauM': 86400,
                'ZS_TauO': 0,
                'ZS_TauR': 10800,
                'ZS_dtF': 0,
                'ZS_dtM': 86400,
                'ZS_dtO': 0,
                'ZS_dtR': 900,
                'ZS_knorm_init': 0,
                'ZS_phifac': 0,
                'ZS_xnorm_init': 0,
                'babsmax_file': None,
                'dam_tot_id_file': None,
                'dam_use_id_file': None,
                'for_tot_id_file': 'tests/data/for_tot_id.csv',
                'for_use_id_file': 'tests/data/for_use_id.csv',
                'k_file': 'input/k.csv',
                'kfac_file': None,
                'obs_tot_id_file': None,
                'obs_use_id_file': None,
                'rapid_connect_file': 'tests/data/rapid_connect.csv',
                'riv_bas_id_file': 'tests/data/riv_bas_id.csv',
                'x_file': 'input/x.csv',
                'xfac_file': None}

    assert params == expected
