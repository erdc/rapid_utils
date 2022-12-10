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
        self.simulation_time_step_s = simulation_time_step_s
        self.runoff_variable_name = runoff_variable_name
        self.runoff_time_axis = runoff_time_axis
        self.input_filename = input_filename
        self.input_params = kwargs

        self.default_params = {
            'BS_opt_Qfinal': False,
            'BS_opt_Qinit': False,
            'BS_opt_dam': False,
            'BS_opt_for': False,
            'BS_opt_influence': False,
            'BS_opt_transpose_qout': False,
            'IS_dam_tot': 0,
            'IS_dam_use': 0,
            'IS_for_tot': 0,
            'IS_for_use': 0,
            'IS_max_up': 2,
            'IS_obs_tot': 0,
            'IS_obs_use': 0,
            'IS_opt_phi': 1,
            'IS_opt_routing': 1,
            'IS_opt_run': 1,
            'IS_riv_bas': 0,
            'IS_riv_tot': 0,
            'IS_strt_opt': 0,
            'Qfinal_file': None,
            'Qfor_file': None,
            'Qinit_file': None,
            'Qobs_file': None,
            'Qobsbarrec_file': None,
            'QoutRabsmax_file': None,
            'QoutRabsmin_file': None,
            'Qout_file': None,
            'Vlat_file': None,
            'ZS_TauM': 0,
            'ZS_TauO': 0,
            'ZS_TauR': 0,
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
            'for_tot_id_file': None,
            'for_use_id_file': None,
            'k_file': 'input/k.csv',
            'kfac_file': None,
            'obs_tot_id_file': None,
            'obs_use_id_file': None,
            'rapid_connect_file': 'input/rapid_connect.csv',
            'riv_bas_id_file': 'input/riv_bas_id.csv',
            'x_file': 'input/x.csv',
            'xfac_file': None,}

        self.params = self.default_params.copy()

        self.rapid_str_format_dict = {
            k.lower(): k for k in self.default_params.keys()}

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

    def write_namelist(self, filename=None, header='&NL_namelist'):
        """Write `params` key, value pairs to file."""

        if filename is not None:
            self.output_filename = filename

        end_char = '/'

        with open(self.output_filename, 'w', encoding='utf8') as f:
            f.write(f'{header}\n')
            for key, value in self.params.items():
                if isinstance(value, bool):
                    value = boolean_to_fortran_string(value)
                elif value is None:
                    value = "''"
                if key not in self.default_params.keys():
                    warnings.warn(
                        f'Key {key} not recognized.\n' +
                        f'Writing {key} = {value} ' +
                        'to the namelist.')

                f.write(f'{key} = {value}\n')

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
            lkey = input_key.lower()
            if lkey in self.rapid_str_format_dict.keys():
                rapid_key = self.rapid_str_format_dict[lkey]
                self.params[rapid_key] = input_value
            else:
                self.params[input_key] = input_value

    def parse_riv_bas_id_file(self):
        """Read a "riv_bas_id" file to determine the value for the
        `IS_riv_bas` parameter.

        Returns
        -------
        parsed : dict
            Dictionary with `IS_riv_bas` parameter.
        """
        ftype = 'riv_bas_id'
        filename = self.params['riv_bas_id_file']

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
        filename = self.params['rapid_connect_file']

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
        filename = self.params['for_tot_id_file']

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
        filename = self.params['for_use_id_file']

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
        filename = self.params['Vlat_file']
        runoff_dataset = None
        time_var = None
        runoff_var = None

        print(f'***{filename}***')
        if filename is None:
            warnings.warn(f'No {ftype} file specified.')
        else:
            runoff_dataset = Dataset(filename)

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
            self.n_time_steps = None

        if (self.n_time_steps is not None) and (
                self.simulation_time_step_s is not None):
            self.total_simulation_time_s = (
                self.n_time_steps * self.simulation_time_step_s)
        else:
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
