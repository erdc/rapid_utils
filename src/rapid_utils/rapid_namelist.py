"""Tools for reading, writing, and modifying a RAPID namelist."""

import os
import warnings
import numpy as np
from netCDF4 import Dataset

def formatted_warning(message, category, filename, lineno, file=None,
                      line=None):
    """Format a warning message to give useful information on one line.

    Parameters
    ----------
    message : str
        Message.
    category : s

    Returns
    -------
    out : str
        Formatted message.
    """
    out = f'{filename}:{lineno}: {category.__name__}: {message}\n'

    return out

warnings.formatwarning = formatted_warning

def read_table_from_file(filename, ftype, delimiter=',', dtype=float,
                         usecols=None):
    """Convenience function for extracting data from a text-based file.

    Parameters
    ----------
    filename : str
        Name of file from which to extract data.

    Returns
    -------
    table : ndarray
        Data extracted from file.
    """
    if filename is None:
        warnings.warn(f'No {ftype} file specified.')
        table = None
    elif os.path.exists(filename):
        try:
            table = np.genfromtxt(
                filename, delimiter=delimiter, dtype=dtype, usecols=usecols)
        except IOError:
            warnings.warn(f'Unable to read {ftype} file {filename}.')
    else:
        warnings.warn(f'{ftype} file {filename} not found.')
        table = None

    return table

def boolean_to_fortran_string(python_bool):
    """Convert a python boolean variable to its Fortran equivalent.

    Parameters
    ----------
    python_bool : bool
        Python boolean variable.

    Returns
    -------
    fortran_str : str
        String representation of Fortran boolean variable.
    """

    if python_bool:
        fortran_str = '.TRUE.'
    else:
        fortran_str = '.FALSE.'

    return fortran_str

def fortran_string_to_boolean(fortran_str):
    """Map ".TRUE." to True and any other string to False.

    Parameters
    ----------
    fortran_str : str
        String representation of Fortran boolean variable.

    Returns
    -------
    python_bool : bool
        Python boolean variable.
    """

    python_bool = (fortran_str == '.TRUE.')

    return python_bool

class RAPIDNamelist:

    """Tools for reading, writing and updating a RAPID namelist file."""

    # TODO: Add methods for updating forcing and optimization parameters.
    # TODO: Add methods to parse "Qfor" file.

    def __init__(self,
                 output_filename='rapid_namelist',
                 simulation_time_step_s=None,
                 input_filename=None,
                 runoff_variable_name='m3_riv',
                 runoff_time_axis=0,
                 **kwargs):
        """Constructor for the RAPIDNamelist class.

        Parameters
        ----------
        output_filename : str
            Name of file to which the namelist is to be written.
        simulation_time_step_s : int (optional)
            Simulation time_step in seconds. This only needs to be specified
            in cases where the "Vlat" file does not contain a valid "time"
            variable.
        input_filename : str (optional)
            Name of existing namelist file from which to read parameters.
        runoff_variable_name : str (optional)
            Name of the runoff variable in "Vlat" netCDF file.
        runoff_time_axis : int (optional)
            Axis of the runoff variable corresponding to the time dimension.
            This is only used if the "Vlat" file does not contain a valid
            "time" variable.
        """
        self.output_filename = output_filename
        self.simulation_time_step_s = simulation_time_step_s or 0
        self.runoff_variable_name = runoff_variable_name
        self.runoff_time_axis = runoff_time_axis
        self.input_filename = input_filename
        self.input_params = kwargs

        self.default_params = {
            'bs_opt_qfinal':
                {'rapid_name': 'BS_opt_Qfinal', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'bs_opt_qinit':
                {'rapid_name': 'BS_opt_Qinit', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'bs_opt_dam':
                {'rapid_name': 'BS_opt_dam', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'bs_opt_for':
                {'rapid_name': 'BS_opt_for', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'bs_opt_influence':
                {'rapid_name': 'BS_opt_influence', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'bs_opt_transpose_qout':
                {'rapid_name': 'BS_opt_transpose_qout', 'value': False,
                 'cast':  fortran_string_to_boolean},
            'is_dam_tot':
                {'rapid_name': 'IS_dam_tot', 'value': 0, 'cast':  int},
            'is_dam_use':
                {'rapid_name': 'IS_dam_use', 'value': 0, 'cast':  int},
            'is_for_tot':
                {'rapid_name': 'IS_for_tot', 'value': 0, 'cast':  int},
            'is_for_use':
                {'rapid_name': 'IS_for_use', 'value': 0, 'cast':  int},
            'is_max_up':
                {'rapid_name': 'IS_max_up', 'value': 2, 'cast':  int},
            'is_obs_tot':
                {'rapid_name': 'IS_obs_tot', 'value': 0, 'cast':  int},
            'is_obs_use':
                {'rapid_name': 'IS_obs_use', 'value': 0, 'cast':  int},
            'is_opt_phi':
                {'rapid_name': 'IS_opt_phi', 'value': 1, 'cast':  int},
            'is_opt_routing':
                {'rapid_name': 'IS_opt_routing', 'value': 1, 'cast':  int},
            'is_opt_run':
                {'rapid_name': 'IS_opt_run', 'value': 1, 'cast':  int},
            'is_riv_bas':
                {'rapid_name': 'IS_riv_bas', 'value': 0, 'cast':  int},
            'is_riv_tot':
                {'rapid_name': 'IS_riv_tot', 'value': 0, 'cast':  int},
            'is_strt_opt':
                {'rapid_name': 'IS_strt_opt', 'value': 0, 'cast':  int},
            'qfinal_file':
                {'rapid_name': 'Qfinal_file', 'value': None, 'cast':  str},
            'qfor_file':
                {'rapid_name': 'Qfor_file', 'value': None, 'cast':  str},
            'qinit_file':
                {'rapid_name': 'Qinit_file', 'value': None, 'cast':  str},
            'qobs_file':
                {'rapid_name': 'Qobs_file', 'value': None, 'cast':  str},
            'qobsbarrec_file':
                {'rapid_name': 'Qobsbarrec_file', 'value': None, 'cast':  str},
            'qoutrabsmax_file':
                {'rapid_name': 'QoutRabsmax_file', 'value': None, 'cast':  str},
            'qoutrabsmin_file':
                {'rapid_name': 'QoutRabsmin_file', 'value': None, 'cast':  str},
            'qout_file':
                {'rapid_name': 'Qout_file', 'value': None, 'cast':  str},
            'vlat_file':
                {'rapid_name': 'Vlat_file', 'value': None, 'cast':  str},
            'zs_taum':
                {'rapid_name': 'ZS_TauM', 'value': 0, 'cast':  int},
            'zs_tauo':
                {'rapid_name': 'ZS_TauO', 'value': 0, 'cast':  int},
            'zs_taur':
                {'rapid_name': 'ZS_TauR', 'value': 0, 'cast':  int},
            'zs_dtf':
                {'rapid_name': 'ZS_dtF', 'value': 0, 'cast':  int},
            'zs_dtm':
                {'rapid_name': 'ZS_dtM', 'value': 86400, 'cast':  int},
            'zs_dto':
                {'rapid_name': 'ZS_dtO', 'value': 0, 'cast':  int},
            'zs_dtr':
                {'rapid_name': 'ZS_dtR', 'value': 900, 'cast':  int},
            'zs_knorm_init':
                {'rapid_name': 'ZS_knorm_init', 'value': 0, 'cast':  int},
            'zs_phifac':
                {'rapid_name': 'ZS_phifac', 'value': 0, 'cast':  int},
            'zs_xnorm_init':
                {'rapid_name': 'ZS_xnorm_init', 'value': 0, 'cast':  int},
            'babsmax_file':
                {'rapid_name': 'babsmax_file', 'value': None, 'cast':  str},
            'dam_tot_id_file':
                {'rapid_name': 'dam_tot_id_file', 'value': None, 'cast':  str},
            'dam_use_id_file':
                {'rapid_name': 'dam_use_id_file', 'value': None, 'cast':  str},
            'for_tot_id_file':
                {'rapid_name': 'for_tot_id_file', 'value': None, 'cast':  str},
            'for_use_id_file':
                {'rapid_name': 'for_use_id_file', 'value': None, 'cast':  str},
            'k_file':
                {'rapid_name': 'k_file', 'value': 'input/k.csv', 'cast':  str},
            'kfac_file':
                {'rapid_name': 'kfac_file', 'value': None, 'cast':  str},
            'obs_tot_id_file':
                {'rapid_name': 'obs_tot_id_file', 'value': None, 'cast':  str},
            'obs_use_id_file':
                {'rapid_name': 'obs_use_id_file', 'value': None, 'cast':  str},
            'rapid_connect_file':
                {'rapid_name': 'rapid_connect_file',
                 'value': 'input/rapid_connect.csv', 'cast':  str},
            'riv_bas_id_file':
                {'rapid_name': 'riv_bas_id_file',
                 'value': 'input/riv_bas_id.csv', 'cast':  str},
            'x_file':
                {'rapid_name': 'x_file', 'value': 'input/x.csv', 'cast':  str},
            'xfac_file':
                {'rapid_name': 'xfac_file', 'value': None, 'cast':  str}}

        self.params = self.default_params.copy()

        self.update_params(self.input_params)

        # Parameters to be determined.
        self.total_simulation_time_s = None
        self.input_time_array = np.array([])
        self.input_time_units = ''
        self.n_time_steps = None

    def read_namelist(self, filename=None, comment_char='!'):
        """Read an existing namelist file.

        Parameters
        ----------
        filename : str
            Name of existing namelist file to be read.
        comment_char : str (optional)
            Character indicating that a line should be ignored.

        Returns
        -------
        parsed : dict
             Dictionary with parameters parsed from the namelist file.
        """
        parsed = {}

        if filename is None:
            warnings.warn('No input filename specified. Returning empty '
                          + 'dictionary.')
        elif not os.path.exists(filename):
            warnings.warn(f'File {filename} not found. Returning empty '
                          + 'dictionary.')
        else:
            with open(filename, 'r', encoding='utf8') as f:
                for line in f:
                    if line.startswith(comment_char):
                        continue

                    try:
                        comps = line.split('=')
                        key = comps[0].strip()
                        value = comps[1].strip()
                        parsed[key] = value
                    except (IndexError, IOError):
                        continue

        return parsed

    def write_namelist(self, filename=None, header='&NL_namelist',
                       end_char='/'):
        """Write `params` key, value pairs to file."""

        if filename is not None:
            self.output_filename = filename

        with open(self.output_filename, 'w', encoding='utf8') as f:
            f.write(f'{header}\n')
            for key, item in self.params.items():
                name = item['rapid_name']
                value = item['value']
                cast = item['cast']

                if key not in self.default_params.keys():
                    warnings.warn(
                        f'Key {name} not recognized.\n' +
                        f'Writing {name} = {value} to the namelist.')

                if item['value'] is None:
                    value = "''"
                elif cast == str:
                    value = f"'{value}'"
                elif isinstance(value, bool):
                    value = boolean_to_fortran_string(value)

                f.write(f'{name} = {value}\n')

            f.write(end_char)

    def update_params(self, new_param_dict, purge=False):
        """Add or modify `params` dictionary entries.

        Parameters
        ----------
        new_param_dict : dict
            Dictionary with key, value pairs to update params dictionary.
        purge : bool, optional
            If True, purge the existing params dictionary.
        """

        if purge:
            self.params = {}

        for input_key, input_value in new_param_dict.items():
            key = input_key.lower()
            if key in self.default_params.keys():
                self.params[key] = self.default_params[key].copy()
            else:
                self.params[key] = {}
                self.params[key]['rapid_name'] = input_key
                self.params[key]['cast'] = str

            cast = self.params[key]['cast']

            if input_value in ["''", '']:
                value = None
            else:
                value = cast(input_value)

            if isinstance(value, str):
                value = value.strip("'")
                value = value.strip('"')

            self.params[key]['value'] = value

    def parse_riv_bas_id_file(self):
        """Read a "riv_bas_id" file to determine the value for the
        `IS_riv_bas` parameter.

        Returns
        -------
        parsed : dict
            Dictionary with `IS_riv_bas` parameter.
        """
        ftype = 'riv_bas_id'
        filename = self.params['riv_bas_id_file']['value']

        riv_bas_id_table = read_table_from_file(
            filename, ftype, delimiter=',', usecols=0, dtype=int)

        if riv_bas_id_table is None:
            is_riv_bas = 0
        else:
            is_riv_bas = int(riv_bas_id_table.size)

        parsed = {'IS_riv_bas': is_riv_bas}

        return parsed

    def parse_connectivity_file(self):
        """Read a "rapid_connect" file to determine the value for the
        `IS_riv_tot` and `IS_max_up` parameters.

        Returns
        -------
        parsed : dict
            Dictionary with `IS_riv_tot` and `IS_max_upstream` parameters.
        """
        ftype = 'rapid_connect'
        filename = self.params['rapid_connect_file']['value']

        rapid_connect_table = read_table_from_file(
            filename, ftype, delimiter=',', dtype=int)

        if rapid_connect_table is None:
            is_riv_tot = 0
            is_max_up = 0
        else:
            is_riv_tot = int(rapid_connect_table.shape[0])
            is_max_up = int(rapid_connect_table[:, 2].max())

        parsed = {'IS_riv_tot': is_riv_tot,
                  'IS_max_up': is_max_up}

        return parsed

    def parse_forcing_tot_file(self):
        """Read a "for_tot_id_file" file to determine the value for the
        `IS_for_tot` parameter.

        Returns
        -------
        parsed : dict
            Dictionary with `IS_for_tot` parameter.
        """
        ftype = 'for_tot_id'
        filename = self.params['for_tot_id_file']['value']

        for_tot_table = read_table_from_file(
            filename, ftype, delimiter=',', usecols=0, dtype=int)

        if for_tot_table is None:
            is_for_tot = 0
        else:
            is_for_tot = int(for_tot_table.size)

        parsed = {'IS_for_tot': is_for_tot}

        return parsed

    def parse_forcing_use_file(self):
        """Read a "for_tot_use_file" file to determine the value for the
        `IS_for_use` parameter.

        Returns
        -------
        parsed : dict
            Dictionary with `IS_for_tot` parameter.
        """
        ftype = 'for_use_id'
        filename = self.params['for_use_id_file']['value']

        for_use_table = read_table_from_file(
            filename, ftype, delimiter=',', usecols=0, dtype=int)

        if for_use_table is None:
            is_for_use = 0
        else:
            is_for_use = int(for_use_table.size)

        parsed = {'IS_for_use': is_for_use}

        return parsed

    def parse_vlat_file(self):
        """Read a "vlat" file to determine input time variables.

        Returns
        -------
        time : ndarry
            Input time variable.
        units : str
            Input time units.
        n_time_steps : int
            Number of time steps in the input runoff as determined by the
            "time" variable in the "Vlat" file, or, if no valid time variable
            is found, the length of the time dimension of the runoff variable.
        """
        ftype = 'vlat'
        filename = self.params['vlat_file']['value']
        runoff_dataset = None
        time_var = None
        runoff_var = None

        if filename is None:
            warnings.warn(f'No {ftype} file specified.')
        else:
            try:
                runoff_dataset = Dataset(filename)
            except (IOError, FileNotFoundError):
                warnings.warn(f'Unable to read {fytpe} file {filename}.')

        if runoff_dataset is not None:
            try:
                time_var = runoff_dataset['time']
            except:
                warnings.warn(
                    f'Variable "time" not found in {ftype} file {filename}.')
            try:
                runoff_var = runoff_dataset[self.runoff_variable_name]
            except:
                warnings.warn(f'Variable {self.runoff_variable_name} not ' +
                              f'found in {ftype} file {filename}.')

        if time_var is not None:
            try:
                self.n_time_steps = time_var.shape[0]
            except:
                warnings.warn('Unable to determine shape of "time" variable ' +
                              f'in {ftype} file {filename}.')
            try:
                self.simulation_time_step_s = time_var[1] - time_var[0]
            except:
                warnings.warn('Unable to determine time step size in ' +
                              f'{ftype} file {filename}.')
            try:
                self.input_time_units = time_var.units
            except AttributeError:
                warnings.warn('Unable to access "units" attribute of ' +
                              '"time" variable in {ftype} file {filename}.')
        elif runoff_var is not None:
            self.n_time_steps = runoff_var.shape[self.runoff_time_axis]
        else:
            self.n_time_steps = 0

        self.total_simulation_time_s = (
            self.n_time_steps * self.simulation_time_step_s)
        
        if self.total_simulation_time_s <= 0:
            warnings.warn('Unable to determine length of simulation time.')

        parsed = {'ZS_TauM': self.total_simulation_time_s,
                  'ZS_TauR': self.simulation_time_step_s}

        return parsed

    def main(self):
        """Update parameters from RAPID input files and write namelist."""

        new_params = {}

        if self.input_filename is not None:
            params = self.read_namelist(self.input_filename)
            new_params.update(params)

        parsers = [self.parse_riv_bas_id_file,
                   self.parse_connectivity_file,
                   self.parse_forcing_tot_file,
                   self.parse_forcing_use_file,
                   self.parse_vlat_file]

        for fn in parsers:
            parsed = fn()
            new_params.update(parsed)

        self.update_params(new_params)

        # User specified parameters will overwrite any parameters parsed from
        # file.
        self.update_params(self.input_params)

        self.write_namelist()

if __name__ == '__main__':
    n = RAPIDNamelist()
    n.main()
